#include "ScarfBlend.hpp"

#include "ClipperUtils.hpp"
#include "ExtrusionEntity.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Print.hpp"

#include <boost/log/trivial.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <set>

namespace Slic3r {

namespace {

static double points_length(const Points &pts)
{
    double L = 0.;
    for (size_t i = 1; i < pts.size(); ++i)
        L += (pts[i] - pts[i - 1]).cast<double>().norm();
    return L;
}

// Sub-path of a closed ring by arc length, wrapping around. `ring` must NOT
// repeat its first point at the end. t is in scaled units, any real value.
static Points ring_sub(const Points &ring, const std::vector<double> &cum, double total, double t0, double t1)
{
    Points out;
    if (ring.size() < 2 || total <= 0.)
        return out;

    auto wrap = [total](double t) {
        double r = std::fmod(t, total);
        return r < 0. ? r + total : r;
    };
    auto point_at = [&](double t) -> Point {
        t = wrap(t);
        size_t i = std::upper_bound(cum.begin(), cum.end(), t) - cum.begin();
        if (i == 0) i = 1;
        if (i >= ring.size() + 1) i = ring.size();
        const Point &a = ring[i - 1];
        const Point &b = ring[i % ring.size()];
        double seg = cum[i] - cum[i - 1];
        double f = seg > 0. ? (t - cum[i - 1]) / seg : 0.;
        return Point(coord_t(a.x() + (b.x() - a.x()) * f), coord_t(a.y() + (b.y() - a.y()) * f));
    };

    double span = t1 - t0;
    if (span <= 0.)
        return out;
    span = std::min(span, total);

    out.push_back(point_at(t0));
    // walk every ring vertex strictly inside (t0, t0 + span)
    double t = wrap(t0);
    size_t start_i = std::upper_bound(cum.begin(), cum.end(), t) - cum.begin();
    double walked = cum[start_i] - t;
    size_t i = start_i;
    while (walked < span - EPSILON) {
        // cum has ring.size() + 1 entries, so the segment leaving vertex i must be
        // looked up modulo ring.size() - i keeps counting past the wrap point.
        const size_t j = i % ring.size();
        out.push_back(ring[j]);
        walked += cum[j + 1] - cum[j];
        ++i;
        if (i - start_i > ring.size())   // never walk more than one full turn
            break;
    }
    out.push_back(point_at(t0 + span));

    // drop duplicates introduced at the joins
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

static ExtrusionPath make_path(const ExtrusionPath &tmpl, Points &&pts)
{
    ExtrusionPath p(tmpl.role(), tmpl.mm3_per_mm, tmpl.width, tmpl.height);
    // ExtrusionPath::polyline is a Polyline3 here and Polyline3(const Polyline&)
    // is explicit, so go through it the same way PerimeterGenerator does.
    p.polyline = Polyline3(Polyline(std::move(pts)));
    return p;
}

// TEMPORARY DIAGNOSTICS. Counts why loops are or are not split, so one slice
// shows where the colour is being lost. Remove once the behaviour is settled.
struct BlendStats
{
    size_t loops               = 0;
    size_t split               = 0;
    // bail reasons, in the order they are tested
    size_t bail_no_paths       = 0;
    size_t bail_thin_poly      = 0;
    size_t bail_one_arc        = 0;   // fewer than two fragments found on the loop
    size_t bail_one_region     = 0;   // fragments found, but all the same region
    size_t bail_zero_length    = 0;
    size_t bail_degenerate_arc = 0;
    size_t bail_empty_emit     = 0;
    // shape of what we did manage to classify
    size_t arcs_found          = 0;
    size_t loops_one_region    = 0;   // loops whose arcs all resolved to one region
    size_t loops_multi_region  = 0;
    size_t junction_taper      = 0;
    size_t junction_butt       = 0;
    // health of the tiling that replaced endpoint chaining
    size_t gaps_filled         = 0;
    size_t overlaps_trimmed    = 0;
};

// Where a run of the loop sits, as arc length along the loop.
struct Span
{
    size_t region_idx = 0;
    double t0 = 0.;
    double t1 = 0.;
};

// Arc-length position of the point on the ring nearest to p.
static double project_t(const Points &ring, const std::vector<double> &cum, const Point &p)
{
    const size_t n        = ring.size();
    double       best_d2  = std::numeric_limits<double>::max();
    double       best_t   = 0.;
    const Vec2d  pd       = p.cast<double>();
    for (size_t i = 0; i < n; ++i) {
        const Vec2d  a    = ring[i].cast<double>();
        const Vec2d  ab   = ring[(i + 1) % n].cast<double>() - a;
        const double len2 = ab.squaredNorm();
        const Vec2d  ap   = pd - a;
        const double u    = len2 > 0. ? std::clamp(ap.dot(ab) / len2, 0., 1.) : 0.;
        const double d2   = (ap - ab * u).squaredNorm();
        if (d2 < best_d2) {
            best_d2 = d2;
            best_t  = cum[i] + std::sqrt(len2) * u;
        }
    }
    return best_t;
}

// Cut one loop among the regions. Returns false if the loop should be left alone.
//
// The loop is parametrised by arc length once, and every fragment returned by
// the clipper is placed onto that parametrisation. An earlier version instead
// chained fragments by matching endpoints within a tolerance, which failed on
// 99% of the loops it rejected: the clipper fragments a loop far more finely
// than there are colour boundaries, and any gap, overlap or reversed fragment
// broke the chain and cost the whole loop its colour. Arc length has none of
// those failure modes - gaps and overlaps are resolved arithmetically, and a
// reversed fragment is detected by comparing its length against the two ways
// round the ring.
static bool split_loop(const ExtrusionLoop      &loop,
                       const LayerRegionPtrs    &layerms,
                       const std::vector<ExPolygons> &region_slices,
                       double                    scarf_scaled,
                       std::vector<ExtrusionEntitiesPtr> &out_by_region,
                       BlendStats               &st)
{
    ++st.loops;

    if (loop.paths.empty()) {
        ++st.bail_no_paths;
        return false;
    }

    Polygon poly = loop.polygon();
    if (poly.points.size() < 3) {
        ++st.bail_thin_poly;
        return false;
    }

    // ---- parametrise the loop ------------------------------------------------
    Points ring = poly.points;
    if (ring.size() > 1 && ring.front() == ring.back())
        ring.pop_back();
    if (ring.size() < 3) {
        ++st.bail_thin_poly;
        return false;
    }

    std::vector<double> cum(ring.size() + 1, 0.);
    for (size_t i = 0; i < ring.size(); ++i)
        cum[i + 1] = cum[i] + (ring[(i + 1) % ring.size()] - ring[i]).cast<double>().norm();
    const double total = cum.back();
    if (total <= 0.) {
        ++st.bail_zero_length;
        return false;
    }

    // ---- place every fragment onto that parametrisation ----------------------
    Polyline          loop_pl = poly.split_at_first_point();
    std::vector<Span> frags;
    for (size_t r = 0; r < layerms.size(); ++r) {
        if (region_slices[r].empty())
            continue;
        // no intersection_pl(Polyline, ExPolygons) overload exists - wrap the subject
        for (const Polyline &pl : intersection_pl(Polylines{loop_pl}, region_slices[r])) {
            if (pl.points.size() < 2)
                continue;
            const double a   = project_t(ring, cum, pl.points.front());
            const double b   = project_t(ring, cum, pl.points.back());
            const double len = points_length(pl.points);
            double fwd = b - a; if (fwd < 0.) fwd += total;
            double bwd = a - b; if (bwd < 0.) bwd += total;
            // The tiling below assumes every fragment starts inside [0, total).
            // project_t can return exactly `total` when a point lands on the
            // seam, which would leave the closing run zero-length.
            auto push = [&](double t0, double span_len) {
                if (t0 >= total) t0 -= total;
                if (t0 < 0.)     t0 += total;
                frags.push_back(Span{r, t0, t0 + span_len});
            };
            if (std::abs(fwd - len) <= std::abs(bwd - len))
                push(a, fwd);
            else
                push(b, bwd);
        }
    }

    st.arcs_found += frags.size();
    {
        std::set<size_t> seen;
        for (const Span &f : frags)
            seen.insert(f.region_idx);
        if (seen.size() < 2) {
            ++st.loops_one_region;
            ++st.bail_one_region;
            return false;                       // single colour, nothing to do
        }
        ++st.loops_multi_region;
    }
    if (frags.size() < 2) {
        ++st.bail_one_arc;
        return false;
    }

    // ---- resolve the fragments into a clean tiling of the ring ---------------
    std::sort(frags.begin(), frags.end(),
              [](const Span &x, const Span &y) { return x.t0 < y.t0; });

    std::vector<Span> spans;
    for (const Span &f : frags) {
        if (spans.empty()) {
            spans.push_back(f);
            continue;
        }
        Span &prev = spans.back();
        if (f.t0 < prev.t1 - EPSILON) {          // overlapping coverage
            ++st.overlaps_trimmed;
            if (f.t1 <= prev.t1)
                continue;                        // wholly inside the previous run
            spans.push_back(Span{f.region_idx, prev.t1, f.t1});
        } else if (f.t0 > prev.t1 + EPSILON) {   // nothing claimed this stretch
            ++st.gaps_filled;
            prev.t1 = f.t0;                      // let the previous run carry it
            spans.push_back(f);
        } else {
            spans.push_back(f);
        }
    }

    // close the ring - the last run absorbs whatever is left at the seam
    if (spans.size() < 2) {
        ++st.bail_one_region;
        return false;
    }
    spans.back().t1 = spans.front().t0 + total;

    // merge neighbours of the same region, including across the seam
    for (size_t i = 0; i + 1 < spans.size();) {
        if (spans[i].region_idx == spans[i + 1].region_idx) {
            spans[i].t1 = spans[i + 1].t1;
            spans.erase(spans.begin() + i + 1);
        } else
            ++i;
    }
    if (spans.size() > 1 && spans.front().region_idx == spans.back().region_idx) {
        spans.front().t0 = spans.back().t0 - total;
        spans.pop_back();
    }
    if (spans.size() < 2) {
        ++st.bail_one_region;
        return false;
    }

    // ---- per-junction scarf --------------------------------------------------
    const size_t        n_spans = spans.size();
    std::vector<double> half(n_spans, 0.);
    const double        min_arc = scale_(0.01);

    for (const Span &s : spans)
        if (s.t1 - s.t0 < min_arc) {
            ++st.bail_degenerate_arc;
            return false;
        }

    for (size_t k = 0; k < n_spans; ++k) {
        const Span  &prev   = spans[(k + n_spans - 1) % n_spans];
        const Span  &cur    = spans[k];
        // A junction may eat at most 40% of either run it joins, so the two
        // junctions of a run always leave a body behind.
        const double budget = 0.4 * std::min(prev.t1 - prev.t0, cur.t1 - cur.t0);
        half[k] = std::min(scarf_scaled * 0.5, budget);
        if (half[k] < scarf_scaled * 0.5)
            ++st.junction_butt;
        else
            ++st.junction_taper;
    }

    // ---- emit ----------------------------------------------------------------
    const ExtrusionPath &tmpl = loop.paths.front();
    using Slope = ExtrusionPathSloped::Slope;

    for (size_t k = 0; k < n_spans; ++k) {
        const Span  &s          = spans[k];
        const double half_start = half[k];
        const double half_end   = half[(k + 1) % n_spans];

        // Emit the three segments as separate heap-allocated paths, in order.
        //
        // They must NOT be wrapped in an ExtrusionEntityCollection: an island's
        // children are spliced straight into ObjectByExtruder::Island::Region::
        // perimeters (GCode.cpp:8312) and each one is handed to extrude_entity(),
        // which only accepts ExtrusionPath, ExtrusionMultiPath and ExtrusionLoop -
        // a collection there throws InvalidArgument (GCode.cpp:6168).
        //
        // ExtrusionMultiPath is no good either: its `paths` is a
        // std::vector<ExtrusionPath> by value, so an ExtrusionPathSloped stored in
        // it is sliced and GCode::_extrude's dynamic_cast finds no slope. Bare
        // paths satisfy extrude_entity() and keep the slope, since
        // ExtrusionPathSloped derives from ExtrusionPath.
        const size_t before = out_by_region[s.region_idx].size();

        Points lead = ring_sub(ring, cum, total, s.t0 - half_start, s.t0 + half_start);
        if (lead.size() >= 2)
            out_by_region[s.region_idx].emplace_back(
                new ExtrusionPathSloped(make_path(tmpl, std::move(lead)), Slope{1., 0.}, Slope{1., 1.}));

        Points body = ring_sub(ring, cum, total, s.t0 + half_start, s.t1 - half_end);
        if (body.size() >= 2)
            out_by_region[s.region_idx].emplace_back(
                new ExtrusionPath(make_path(tmpl, std::move(body))));

        Points tail = ring_sub(ring, cum, total, s.t1 - half_end, s.t1 + half_end);
        if (tail.size() >= 2)
            out_by_region[s.region_idx].emplace_back(
                new ExtrusionPathSloped(make_path(tmpl, std::move(tail)), Slope{1., 1.}, Slope{1., 0.}));

        if (out_by_region[s.region_idx].size() == before) {
            ++st.bail_empty_emit;
            return false;
        }
    }
    return true;
}

static void collect_loops(ExtrusionEntityCollection &coll, std::vector<ExtrusionLoop *> &out)
{
    for (ExtrusionEntity *e : coll.entities) {
        if (auto *c = dynamic_cast<ExtrusionEntityCollection *>(e))
            collect_loops(*c, out);
        else if (auto *l = dynamic_cast<ExtrusionLoop *>(e))
            out.push_back(l);
    }
}

// ExtrusionEntityCollection has no remove_entities(), and collect_loops() descends
// into nested collections, so erase the consumed loops wherever they actually live.
static void erase_entities(ExtrusionEntityCollection &coll, const std::set<const ExtrusionEntity *> &drop)
{
    auto &ents = coll.entities;
    for (auto it = ents.begin(); it != ents.end();) {
        if (auto *c = dynamic_cast<ExtrusionEntityCollection *>(*it)) {
            erase_entities(*c, drop);
            ++it;
        } else if (drop.count(*it)) {
            delete *it;
            it = ents.erase(it);
        } else
            ++it;
    }
}

// Every top-level entity of LayerRegion::perimeters is one island, and GCode.cpp
// static_casts it to ExtrusionEntityCollection without checking in release builds
// (see the assert at GCode.cpp:5086). So arcs are handed back inside an island,
// never appended to perimeters as bare top-level paths.
//
// The island is left sortable on purpose. Region::append() splices a sortable
// island's children into region.perimeters one by one, which is what we want;
// a no_sort island is instead pushed whole and would reach extrude_entity() as
// a collection, which throws (GCode.cpp:6168).
static ExtrusionEntityCollection *island_for(LayerRegion *lr)
{
    auto *coll = new ExtrusionEntityCollection();
    coll->no_sort = false;
    lr->perimeters.entities.emplace_back(coll);
    return coll;
}

} // namespace

bool apply_scarf_blend(LayerRegion *source, const LayerRegionPtrs &layerms, double scarf_width_mm)
{
    if (source == nullptr || layerms.size() < 2 || scarf_width_mm <= 0.) {
        BOOST_LOG_TRIVIAL(info) << "ScarfBlend: skipped, regions=" << layerms.size()
                                << " width=" << scarf_width_mm;
        return false;
    }

    std::vector<ExPolygons> region_slices(layerms.size());
    for (size_t i = 0; i < layerms.size(); ++i)
        region_slices[i] = to_expolygons(layerms[i]->slices.surfaces);

    const double scarf_scaled = scale_(scarf_width_mm);
    size_t       split_count  = 0;
    BlendStats   st;

    // Work island by island so the arcs stay grouped the way the rest of the
    // pipeline expects. Appending them to perimeters as separate top-level
    // entities would make every arc look like its own island.
    const ExtrusionEntitiesPtr islands = source->perimeters.entities;
    for (ExtrusionEntity *island_ee : islands) {
        auto *island = dynamic_cast<ExtrusionEntityCollection *>(island_ee);
        if (island == nullptr)
            continue;

        std::vector<ExtrusionLoop *> loops;
        collect_loops(*island, loops);
        if (loops.empty())
            continue;

        std::vector<ExtrusionEntitiesPtr> out_by_region(layerms.size());
        std::vector<ExtrusionLoop *>      consumed;

        for (ExtrusionLoop *loop : loops) {
            std::vector<ExtrusionEntitiesPtr> staged(layerms.size());
            if (split_loop(*loop, layerms, region_slices, scarf_scaled, staged, st)) {
                for (size_t i = 0; i < staged.size(); ++i)
                    for (ExtrusionEntity *e : staged[i])
                        out_by_region[i].emplace_back(e);
                consumed.push_back(loop);
            } else {
                for (auto &v : staged)
                    for (ExtrusionEntity *e : v)
                        delete e;
            }
        }

        if (consumed.empty())
            continue;

        // Drop the loops we replaced, keep the ones we could not handle.
        erase_entities(*island, std::set<const ExtrusionEntity *>(consumed.begin(), consumed.end()));
        split_count += consumed.size();
        st.split    += consumed.size();

        for (size_t i = 0; i < layerms.size(); ++i) {
            if (out_by_region[i].empty())
                continue;
            if (layerms[i] == source)
                // same region: the arcs belong to the island we just edited
                island->append(std::move(out_by_region[i]));
            else
                // other region: give them an island of their own over there
                island_for(layerms[i])->append(std::move(out_by_region[i]));
        }
    }

    // TEMPORARY DIAGNOSTICS - one line per layer per merged region group.
    // `1region` is the key figure: loops whose arcs all resolved to a single
    // region are loops where classification found no colour boundary at all.
    BOOST_LOG_TRIVIAL(info)
        << "ScarfBlend z=" << (source->layer() ? source->layer()->print_z : -1.)
        << " regions=" << layerms.size()
        << " loops=" << st.loops
        << " split=" << st.split
        << " arcs=" << st.arcs_found
        << " 1region=" << st.loops_one_region
        << " Nregion=" << st.loops_multi_region
        << " | bail nopaths=" << st.bail_no_paths
        << " thin=" << st.bail_thin_poly
        << " onearc=" << st.bail_one_arc
        << " oneregion=" << st.bail_one_region
        << " zerolen=" << st.bail_zero_length
        << " degen=" << st.bail_degenerate_arc
        << " emptyemit=" << st.bail_empty_emit
        << " | junctions taper=" << st.junction_taper
        << " butt=" << st.junction_butt
        << " | tiling gaps=" << st.gaps_filled
        << " overlaps=" << st.overlaps_trimmed;

    if (split_count == 0)
        return false;

    return true;
}

} // namespace Slic3r
