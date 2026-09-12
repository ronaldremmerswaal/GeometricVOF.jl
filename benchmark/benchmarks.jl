# Reproducible microbenchmarks for the public 2D VOF operations. Run from the
# repository root with `julia --project=. benchmark/benchmarks.jl`. The fixtures
# deliberately use a non-axis-aligned interface and include both the allocation-
# free workspace API and its allocating convenience counterpart. This makes
# regressions in the clipping/reconstruction hot path visible while keeping setup
# out of the measurements.
using BenchmarkTools
using GeometricVOF
using Meshes
using Printf
using StaticArrays
using Unitful

function reconstruction_fixture()
    mesh = CartesianGrid((3, 3), (-0.5, -0.5), (1 / 3, 1 / 3))
    central_index = 5
    normal = SVector(cos(0.73), sin(0.73))
    reference = PlanarHS(normal, 0.12u"m")
    fractions = [smeasure(reference, cell) / smeasure(cell) for cell in mesh]
    central = mesh[central_index]
    initial = PlanarHS(SVector(cos(1.2), sin(1.2)), 0u"m")
    cells = view(mesh, 1:9)
    static_central = StaticNgon(central)
    GeometricVOF.copy!(static_central, central)
    static_cells = [let polygon = StaticNgon(cell)
        GeometricVOF.copy!(polygon, cell)
        polygon
    end for cell in cells]
    cell_areas = smeasure.(cells)
    output = StaticNgon(central)
    workspace = StaticNgon(central)
    shifts = MVector{32, Float64}(undef)

    return (; mesh, central, cells, fractions, initial, static_central,
        static_cells, cell_areas, output, workspace, shifts)
end

const QUAD = Quadrangle((-0.1, -0.1), (0.9, -0.1), (0.9, 0.9), (-0.1, 0.9))
const HALFSPACE = PlanarHS(SVector(cos(0.61), sin(0.61)), 0.42u"m")
const CLIP_OUTPUT = StaticNgon(QUAD)
const CLIP_WORKSPACE = StaticNgon(QUAD)
const SHIFT_WORKSPACE = StaticNgon(QUAD)
const SHIFT_SCRATCH = MVector{32, Float64}(undef)
const SHIFT_NORMAL = SVector(cos(0.61), sin(0.61))
const SHIFT_AREA = 0.37u"m^2"
const CURVED_LEVELSET = (x, y) -> y - (0.35u"m" + 0.15 * (x / u"m")^2 * u"m")
const RECONSTRUCTION = reconstruction_fixture()
const PARABOLA_ORIGIN = Point(0u"m", 0u"m")
const PARABOLA = Parabola(SVector(0.0, 1.0), 0.4u"m", -0.4u"m^-1", PARABOLA_ORIGIN)
const PARABOLA_OUTPUT = StaticParabolicNgon(QUAD, PARABOLA)
const PARABOLA_AREA = smeasure(PARABOLA, QUAD)
const PARABOLA_MOMENT = moments(PARABOLA, QUAD)[2]
const MOF_INITIAL = PlanarHS(SVector(cos(1.2), sin(1.2)), 0u"m")
const MOF_REFERENCE = PlanarHS(SVector(cos(0.73), sin(0.73)), 0.12u"m")
const MOF_FRACTION = smeasure(MOF_REFERENCE, QUAD) / smeasure(QUAD)
const MOF_MOMENT = moments(MOF_REFERENCE, QUAD)[2]
const PMOF_INITIAL = PlanarHS(SVector(1.0, 0.0), 0u"m")
const PMOF_FRACTION = PARABOLA_AREA / smeasure(QUAD)
const STENCIL_PARABOLA = Parabola(SVector(0.0, 1.0), 0.1u"m", -0.5u"m^-1", PARABOLA_ORIGIN)
const STENCIL_FRACTIONS = [smeasure(STENCIL_PARABOLA, cell) / smeasure(cell)
    for cell in RECONSTRUCTION.cells]
const STENCIL_WORKSPACE = StaticParabolicNgon(RECONSTRUCTION.central, STENCIL_PARABOLA)

const SUITE = BenchmarkGroup()
SUITE["clip"]["convenience"] = @benchmarkable intersect($QUAD, $HALFSPACE)
SUITE["clip"]["workspace"] = @benchmarkable intersect!($CLIP_OUTPUT, $QUAD, $HALFSPACE)
SUITE["clip"]["parabola_workspace"] = @benchmarkable intersect!($PARABOLA_OUTPUT, $QUAD, $PARABOLA)
SUITE["area"]["halfspace"] = @benchmarkable smeasure($HALFSPACE, $QUAD)
SUITE["area"]["levelset"] = @benchmarkable smeasure($CURVED_LEVELSET, $QUAD)
SUITE["area"]["parabola"] = @benchmarkable smeasure($PARABOLA, $QUAD; workspace=$PARABOLA_OUTPUT)
SUITE["moments"]["parabola"] = @benchmarkable moments($PARABOLA, $QUAD; workspace=$PARABOLA_OUTPUT)
SUITE["shift"]["workspace"] = @benchmarkable shift(
    $QUAD, $SHIFT_NORMAL, $SHIFT_AREA;
    workspace=$SHIFT_WORKSPACE, shift_workspace=$SHIFT_SCRATCH,
)
SUITE["reconstruct"]["convenience"] = @benchmarkable reconstruct(
    $RECONSTRUCTION.initial,
    $RECONSTRUCTION.fractions[5],
    $RECONSTRUCTION.central,
    $RECONSTRUCTION.fractions,
    $RECONSTRUCTION.cells;
    workspace=$RECONSTRUCTION.workspace,
    shift_workspace=$RECONSTRUCTION.shifts,
)
SUITE["reconstruct"]["workspace"] = @benchmarkable reconstruct!(
    $RECONSTRUCTION.output,
    $RECONSTRUCTION.initial,
    $RECONSTRUCTION.fractions[5],
    $RECONSTRUCTION.static_central,
    $RECONSTRUCTION.fractions,
    $RECONSTRUCTION.static_cells,
    $RECONSTRUCTION.cell_areas;
    workspace=$RECONSTRUCTION.workspace,
    shift_workspace=$RECONSTRUCTION.shifts,
)
SUITE["reconstruct"]["mof"] = @benchmarkable mof(
    $MOF_INITIAL, $MOF_FRACTION, $MOF_MOMENT, $QUAD;
    workspace=$CLIP_WORKSPACE, shift_workspace=$SHIFT_SCRATCH,
)
SUITE["reconstruct"]["pmof"] = @benchmarkable pmof(
    $PMOF_INITIAL, $PMOF_FRACTION, $PARABOLA_MOMENT, $PARABOLA.curvature, $QUAD;
    origin=$PARABOLA_ORIGIN, workspace=$PARABOLA_OUTPUT,
) samples=1_000 evals=1
SUITE["reconstruct"]["plvira"] = @benchmarkable plvira(
    $PMOF_INITIAL, $STENCIL_FRACTIONS[5], $STENCIL_PARABOLA.curvature,
    $RECONSTRUCTION.central, $STENCIL_FRACTIONS, $RECONSTRUCTION.cells;
    cmeasures=$RECONSTRUCTION.cell_areas, origin=$PARABOLA_ORIGIN,
    workspace=$STENCIL_WORKSPACE,
) samples=1_000 evals=1
SUITE["reconstruct"]["prost"] = @benchmarkable prost(
    $PMOF_INITIAL, $STENCIL_FRACTIONS[5], $RECONSTRUCTION.central,
    $STENCIL_FRACTIONS, $RECONSTRUCTION.cells;
    cmeasures=$RECONSTRUCTION.cell_areas, origin=$PARABOLA_ORIGIN,
    curvature_bounds=(-1u"m^-1", 1u"m^-1"), maxiters=100,
    workspace=$STENCIL_WORKSPACE,
) samples=1_000 evals=1

if abspath(PROGRAM_FILE) == @__FILE__
    println("GeometricVOF 2D benchmark suite")
    println("Julia ", VERSION, " on ", Sys.MACHINE)
    results = run(SUITE; verbose=true)
    println("\nSummary (minimum time, allocation estimate)")
    for (group_name, group) in results
        for (name, trial) in group
            println(rpad("$(group_name)/$(name)", 28),
                rpad(BenchmarkTools.prettytime(minimum(trial).time), 14),
                @sprintf("%d bytes", minimum(trial).memory))
        end
    end
end
