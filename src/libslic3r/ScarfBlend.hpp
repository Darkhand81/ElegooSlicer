#ifndef slic3r_ScarfBlend_hpp_
#define slic3r_ScarfBlend_hpp_

#include "Layer.hpp"

namespace Slic3r {

// ScarfBlend
//
// When several painted (per-filament) LayerRegions are perimeter-compatible,
// Layer::make_perimeters() generates ONE set of loops around their merged
// slices and parks them all on `source`. Normally that would print the whole
// wall in a single color.
//
// This pass cuts each of those loops at the color boundaries and hands the
// resulting open arcs to the LayerRegion that owns them, so the wall runs
// continuously around the part and simply changes filament partway round.
// It never turns inward, so no wall is generated along the color boundary.
//
// Each arc is extended half of `scarf_width_mm` past both of its junctions and
// those overlaps are emitted as ExtrusionPathSloped with e_ratio ramping 0->1
// and 1->0 (z_ratio stays 1.0, so this is a flow ramp in the plane, not a Z
// slope). The two neighbouring arcs traverse the same overlap in the same
// direction with complementary flow, so deposited material sums to one bead.
//
// Returns true if any loop was redistributed. On any geometric surprise it
// bails out for that loop and leaves it untouched, so a failure degrades to
// current behaviour rather than to a broken toolpath.
bool apply_scarf_blend(LayerRegion *source, const LayerRegionPtrs &layerms, double scarf_width_mm);

} // namespace Slic3r

#endif // slic3r_ScarfBlend_hpp_
