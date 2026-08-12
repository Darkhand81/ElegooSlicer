#include "ScarfBlend.hpp"

#include "ClipperUtils.hpp"
#include "ExtrusionEntity.hpp"
#include "ExtrusionEntityCollection.hpp"
#include "Print.hpp"

#include <boost/log/trivial.hpp>
#include <algorithm>
#include <cmath>
#include <memory>
#include <set>

namespace Slic3r {

namespace {

// One contiguous run of a loop that belongs to a single LayerRegion.
struct Arc
{
    size_t region_idx = 0;
    Points pts;            // ordered, lying on the loop
    double t_start = 0.;   // arc length along the assembled ring (scaled units)
    double t_end   = 0.;
};

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

// Cut one loop among the regions. Returns false if the loop should be left alone.
static bool split_loop(const ExtrusionLoop      &loop,
                       const LayerRegionPtrs    &layerms,
                       const std::vector<ExPolygons> &region_slices,
                       double                    scarf_scaled,
                       std::vector<ExtrusionEntitiesPtr> &out_by_region)
{
    if (loop.paths.empty())
        return false;

    Polygon poly = loop.polygon();
    if (poly.points.size() < 3)
        return false;

    Polyline loop_pl = poly.split_at_first_point();   // closed: first == last

    // ---- classify the loop into per-region arcs -----------------------------
    std::vector<Arc> arcs;
    for (size_t r = 0; r < layerms.size(); ++r) {
        if (region_slices[r].empty())
            continue;
        // no intersection_pl(Polyline, ExPolygons) overload exists - wrap the subject
        for (Polyline &pl : intersection_pl(Polylines{loop_pl}, region_slices[r])) {
            if (pl.points.size() >= 2)
                arcs.push_back(Arc{r, std::move(pl.points), 0., 0.});
        }
    }
    if (arcs.size() < 2)
        return false;                                   // single color, nothing to do

    // ---- chain the arcs into ring order -------------------------------------
    // The arcs tile the loop exactly, so each one starts where another ends.
    const double join_eps = scale_(0.02);
    std::vector<bool> used(arcs.size(), false);
    std::vector<Arc>  ordered;
    ordered.push_back(arcs[0]);
    used[0] = true;
    for (size_t guard = 0; guard + 1 < arcs.size(); ++guard) {
        const Point &tail = ordered.back().pts.back();
        size_t best = arcs.size();
        double best_d = join_eps;
        for (size_t i = 0; i < arcs.size(); ++i) {
            if (used[i]) continue;
            double d = (arcs[i].pts.front() - tail).cast<double>().norm();
            if (d < best_d) { best_d = d; best = i; }
        }
        if (best == arcs.size())
            return false;                               // chain broke - bail out
        used[best] = true;
        ordered.push_back(arcs[best]);
    }
    // the last arc must close back onto the first
    if ((ordered.front().pts.front() - ordered.back().pts.back()).cast<double>().norm() > join_eps)
        return false;

    // merge neighbours that belong to the same region (can happen at the loop seam)
    for (size_t i = 0; i + 1 < ordered.size();) {
        if (ordered[i].region_idx == ordered[i + 1].region_idx) {
            ordered[i].pts.insert(ordered[i].pts.end(),
                                  ordered[i + 1].pts.begin() + 1, ordered[i + 1].pts.end());
            ordered.erase(ordered.begin() + i + 1);
        } else ++i;
    }
    if (ordered.size() > 1 && ordered.front().region_idx == ordered.back().region_idx) {
        Arc tail = ordered.back();
        ordered.pop_back();
        tail.pts.insert(tail.pts.end(), ordered.front().pts.begin() + 1, ordered.front().pts.end());
        ordered.front() = tail;
    }
    if (ordered.size() < 2)
        return false;

    // ---- assemble the ring and parametrise it -------------------------------
    Points ring;
    for (Arc &a : ordered) {
        a.t_start = ring.empty() ? 0. : points_length(ring);
        // skip the duplicated joint vertex
        ring.insert(ring.end(), a.pts.begin() + (ring.empty() ? 0 : 1), a.pts.end());
        a.t_end = points_length(ring);
    }
    if (ring.size() > 1 && ring.front() == ring.back())
        ring.pop_back();

    std::vector<double> cum(ring.size() + 1, 0.);
    for (size_t i = 0; i < ring.size(); ++i)
        cum[i + 1] = cum[i] + (ring[(i + 1) % ring.size()] - ring[i]).cast<double>().norm();
    const double total = cum.back();
    if (total <= 0.)
        return false;

    const double half = scarf_scaled * 0.5;
    // every arc has to be long enough to carry two half-scarves plus a body
    for (const Arc &a : ordered)
        if (a.t_end - a.t_start < scarf_scaled * 1.5)
            return false;

    // ---- emit ---------------------------------------------------------------
    const ExtrusionPath &tmpl = loop.paths.front();
    using Slope = ExtrusionPathSloped::Slope;

    for (const Arc &a : ordered) {
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
        const size_t before = out_by_region[a.region_idx].size();

        Points lead = ring_sub(ring, cum, total, a.t_start - half, a.t_start + half);
        if (lead.size() >= 2)
            out_by_region[a.region_idx].emplace_back(
                new ExtrusionPathSloped(make_path(tmpl, std::move(lead)), Slope{1., 0.}, Slope{1., 1.}));

        Points body = ring_sub(ring, cum, total, a.t_start + half, a.t_end - half);
        if (body.size() >= 2)
            out_by_region[a.region_idx].emplace_back(
                new ExtrusionPath(make_path(tmpl, std::move(body))));

        Points tail = ring_sub(ring, cum, total, a.t_end - half, a.t_end + half);
        if (tail.size() >= 2)
            out_by_region[a.region_idx].emplace_back(
                new ExtrusionPathSloped(make_path(tmpl, std::move(tail)), Slope{1., 1.}, Slope{1., 0.}));

        if (out_by_region[a.region_idx].size() == before)
            return false;
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
    if (source == nullptr || layerms.size() < 2 || scarf_width_mm <= 0.)
        return false;

    std::vector<ExPolygons> region_slices(layerms.size());
    for (size_t i = 0; i < layerms.size(); ++i)
        region_slices[i] = to_expolygons(layerms[i]->slices.surfaces);

    const double scarf_scaled = scale_(scarf_width_mm);
    size_t       split_count  = 0;

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
            if (split_loop(*loop, layerms, region_slices, scarf_scaled, staged)) {
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

    if (split_count == 0)
        return false;

    BOOST_LOG_TRIVIAL(trace) << "ScarfBlend: split " << split_count
                             << " loops across " << layerms.size() << " regions";
    return true;
}

} // namespace Slic3r
