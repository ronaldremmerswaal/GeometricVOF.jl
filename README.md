# GeometricVOF

[![Build Status](https://github.com/ronaldremmerswaal/GeometricVOF.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ronaldremmerswaal/GeometricVOF.jl/actions/workflows/CI.yml?query=branch%3Amain)

Fast, two-dimensional geometric volume-of-fluid primitives for polygonal cells
from [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl). Coordinates and
areas retain Meshes.jl's Unitful units.

```julia
using GeometricVOF, Meshes, StaticArrays, Unitful

cell = Quadrangle((0, 0), (1, 0), (1, 1), (0, 1))
interface = PlanarHS(SVector(1.0, 0.0), 0.4u"m")

liquid_area = smeasure(interface, cell)  # 0.4 m²
liquid = intersect(cell, interface)      # clipped polygon, or `nothing`

# Construct the plane that contains a prescribed area in a cell.
offset = shift(cell, SVector(1.0, 0.0), 0.4u"m^2")
```

## Stable 2D API

`PlanarHS(normal, shift)` represents `normal ⋅ x ≤ shift`; `normal(p)` and
`distance(p, point)` expose its geometry. `smeasure` computes signed areas for
polygons, clipped half-spaces, and level sets. `shift` inverts clipped area to
a plane shift. `reconstruct` performs LVIRA reconstruction from a central
volume fraction and neighbouring cells.

The ordinary operations allocate only for values they return. In cell loops,
use `StaticNgon`, `intersect!`, and `reconstruct!` with reusable workspaces:

```julia
out = StaticNgon(cell)       # default capacity is safe for normal clipping
scratch = StaticNgon(cell)
intersect!(out, cell, interface)
area = smeasure(out)
```

`capacity(polygon)` reports a static workspace's vertex limit. Supply a larger
capacity when a workflow can create more intermediate vertices; operations
throw an `ArgumentError` rather than writing past it.

## Moment and parabolic reconstruction

`moments(region)` returns its signed area and global first moment.  `mof`
uses a volume fraction and first moment to reconstruct a planar interface:

```julia
liquid_moment = moments(interface, cell)[2]
plane = mof(initial_plane, fraction, liquid_moment, cell)
```

`Parabola(normal, shift, curvature, origin)` represents the liquid side of
`normal ⋅ (x - origin) - shift + curvature/2 * (tangent ⋅ (x - origin))^2 ≤ 0`.
Parabolic clipping returns a `StaticParabolicNgon`; its marked arc faces make
`smeasure` and `moments` exact without tessellating the curve. `pmof` and
`plvira` take a supplied curvature, while `prost` searches a bounded curvature
range along with the normal:

```julia
curve = pmof(initial_plane, fraction, liquid_moment, curvature, cell)
curve = plvira(initial_plane, fraction, curvature, cell, neighbor_fractions, neighbor_cells)
curve = prost(initial_plane, fraction, cell, neighbor_fractions, neighbor_cells;
              curvature_bounds=(-8u"m^-1", 8u"m^-1"))
```

Reusable `StaticParabolicNgon(cell, curve)` workspaces remove allocations from
repeated parabolic clipping and reconstruction calls. Its capacity defaults to
twice the cell's vertex count plus two, which accommodates a quadratic edge
intersection on every face.

## Performance checks

The reproducible suite covers planar and parabolic clipping, moments,
area-to-shift inversion, and all reconstruction APIs:

```sh
julia --project=. benchmark/benchmarks.jl
```

See [benchmark/README.md](benchmark/README.md) for comparison guidance.
