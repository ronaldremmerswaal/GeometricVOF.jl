"""
Performance analysis for GeometricVOF.jl.

Run with:
    julia --project=. benchmark/performance_analysis.jl

Covers:
  1. High-level @benchmark: reconstruct (full LVIRA pipeline, preallocated workspaces)
  2. Mid-level @benchmark: smeasure(Φ,c), shift, intersect!
  3. Low-level @benchmark: lvira_costfun, donating_region!, brent_min, smeasure(StaticNgon)
  4. Allocation analysis: @allocated (preallocated vs default-arg) + Profile.Allocs on reconstruct
  5. Type stability: JET @report_call / @report_opt
"""

cd(joinpath(@__DIR__, ".."))

using GeometricVOF
using Meshes, Unitful, StaticArrays, LinearAlgebra
using FastGaussQuadrature: gausslegendre
using BenchmarkTools
using Profile
using JET

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

section(title)    = println("\n", "="^72, "\n  ", title, "\n", "="^72)
subsection(title) = println("\n── ", title, " ──")

# ──────────────────────────────────────────────────────────────────────────────
# ❶  SETUP
# ──────────────────────────────────────────────────────────────────────────────
section("SETUP — flower interface on N×N CartesianGrid, N=32")

const N = 32
R = 0.25u"m"
Φ(x, y) = (x - 0.5u"m")^2 + (y - 0.5u"m")^2 -
           R^2 * (1 + 0.1cos(0.1 + 5atan(y - 0.5u"m", x - 0.5u"m")))^2

mesh = CartesianGrid((N, N), (0.0, 0.0), (1.0/N, 1.0/N))
inds = LinearIndices(size(mesh))

println("  Building αs ($(N)×$(N) grid)…")
# ustrip: smeasure(Φ,c)/smeasure(c) is a dimensionless Unitful Quantity which is not <: Real;
# reconstruct requires T <: Real, so we store plain Float64 from the start.
αs = reshape([Float64(ustrip(smeasure(Φ, c) / smeasure(c))) for c ∈ mesh], (N, N))
println("  Interface cells: $(count(α -> 1e-6 < α < 1 - 1e-6, αs))")

# First interior interface cell — use let to avoid soft-scope assignment warnings
i0, j0 = let result = (0, 0)
    for i in 2:N-1, j in 2:N-1
        if 1e-4 < αs[i, j] < 1 - 1e-4
            result = (i, j)
            break
        end
    end
    result
end
i0 > 0 || error("No interior interface cell found — check Φ or grid size")
c_central = mesh[i0, j0]
α_central = αs[i0, j0]   # Float64
println("  Representative cell: ($i0, $j0),  α = $(round(α_central, sigdigits=4))")

# 3×3 neighbourhood (same layout as benchmark/reconstruction.jl)
flat_inds       = inds[i0-1:i0+1, j0-1:j0+1][:]
cs              = view(mesh, flat_inds)
αs_local        = view(αs, i0-1:i0+1, j0-1:j0+1)
cmeasures_local = [smeasure(c) for c in cs]

# ── Pre-allocated workspaces ──────────────────────────────────────────────────
workspace           = StaticNgon(c_central)           # for intersect!, shift, reconstruct
shift_workspace     = MVector{30, Float64}(undef)     # vertex shift values in shift()
gl_data             = gausslegendre(16)                # Gauss-Legendre nodes for smeasure(Φ,c)
smeasure_workspaces = (StaticNgon(c_central, 8 * 8),  # refine_edges! output (nref=8 → ≤32 pts)
                       StaticNgon(c_central, 8 * 8))  # ngon_approx! output

# ── Half-space: initial guess and converged solution ─────────────────────────
θ0      = 0.7   # representative angle (radians) — same as benchmark/shift.jl
𝛈₀      = GeometricVOF.angle_to_normal(θ0)
ref_vol = smeasure(c_central) * α_central  # Quantity{Float64, m^2} — fine for shift/PlanarHS
p0      = PlanarHS(𝛈₀, 0u"m")
p_shift = PlanarHS(θ0, ref_vol, c_central; workspace=workspace, shift_workspace=shift_workspace)

# ── donating_region setup ─────────────────────────────────────────────────────
seg     = Segment(c_central.vertices[1], c_central.vertices[2])
dt      = 0.001u"s"
dr_velo = [(0.5u"m/s", 0.0u"m/s"), (0.5u"m/s", 0.0u"m/s")]
dr_out  = StaticNgon(eltype(seg.vertices))

# ── brent_min proxy (analytical minimum at x = 0.3) ─────────────────────────
_bm_fdf(x::Real) = ((x - 0.3)^2, 2(x - 0.3))

println("Setup complete.")

# ──────────────────────────────────────────────────────────────────────────────
# ❷  SECTION 1 — High-level @benchmark
# ──────────────────────────────────────────────────────────────────────────────
section("SECTION 1 — High-level @benchmark")

subsection("reconstruct  (full LVIRA: angle search + shift + cost × neighbours, preallocated)")
b_reconstruct = @benchmark reconstruct(
    $p0, $α_central, $c_central, $αs_local, $cs, $cmeasures_local;
    workspace=$workspace, shift_workspace=$shift_workspace)
show(stdout, MIME"text/plain"(), b_reconstruct)

# ──────────────────────────────────────────────────────────────────────────────
# ❸  SECTION 2 — Mid-level @benchmark
# ──────────────────────────────────────────────────────────────────────────────
section("SECTION 2 — Mid-level @benchmark")

subsection("smeasure(Φ, c)  (level-set initialisation: edge refinement + root-finding + GL quadrature)")
b_smeasure_phi = @benchmark smeasure($Φ, $c_central;
    workspaces=$smeasure_workspaces, gl_data=$gl_data)
show(stdout, MIME"text/plain"(), b_smeasure_phi)

subsection("shift  (volume-constrained root-finding via parabolic bracketing, preallocated)")
b_shift = @benchmark shift($c_central, $𝛈₀, $ref_vol;
    workspace=$workspace, shift_workspace=$shift_workspace)
show(stdout, MIME"text/plain"(), b_shift)

subsection("intersect!  (Sutherland-Hodgman polygon clipping, preallocated workspace)")
b_intersect = @benchmark intersect!($workspace, $c_central, $p_shift)
show(stdout, MIME"text/plain"(), b_intersect)

# ──────────────────────────────────────────────────────────────────────────────
# ❹  SECTION 3 — Low-level @benchmark
# ──────────────────────────────────────────────────────────────────────────────
section("SECTION 3 — Low-level @benchmark")

subsection("lvira_costfun  (per-angle LVIRA cost + derivative over 3×3 neighbourhood)")
b_costfun = @benchmark GeometricVOF.lvira_costfun(
    $p_shift, $cs, $αs_local, $cmeasures_local, $c_central;
    workspace=$workspace)
show(stdout, MIME"text/plain"(), b_costfun)

subsection("donating_region!  (upstream quadrangle/pentagon, mutating, preallocated out)")
b_donating = @benchmark donating_region!($dr_out, $seg, $dr_velo, $dt)
show(stdout, MIME"text/plain"(), b_donating)

subsection("brent_min  (Newton bracket + Brent root, analytical proxy f(x)=(x-0.3)²)")
b_brent = @benchmark GeometricVOF.brent_min($_bm_fdf, 1.0; xatol=√eps(), step_max=0.5, maxiters=25)
show(stdout, MIME"text/plain"(), b_brent)

subsection("smeasure(StaticNgon)  (bare shoelace formula on pre-clipped polygon)")
intersect!(workspace, c_central, p_shift)  # fill workspace with a real polygon
b_smeasure_ngon = @benchmark smeasure($workspace)
show(stdout, MIME"text/plain"(), b_smeasure_ngon)

# ──────────────────────────────────────────────────────────────────────────────
# ❺  SECTION 4 — Allocation analysis
# ──────────────────────────────────────────────────────────────────────────────
section("SECTION 4 — Allocation analysis")

subsection("@allocated — core functions  (target: 0 for preallocated paths)")

alloc_intersect_pre = @allocated intersect!(workspace, c_central, p_shift)
println("  intersect!  (preallocated workspace)                → $alloc_intersect_pre bytes")

alloc_intersect_def = @allocated intersect(c_central, p_shift)
println("  intersect   (default: allocates StaticNgon + Ngon)  → $alloc_intersect_def bytes")

alloc_shift_pre = @allocated shift(c_central, 𝛈₀, ref_vol; workspace=workspace, shift_workspace=shift_workspace)
println("  shift       (preallocated workspace + shift_wsp)    → $alloc_shift_pre bytes")

alloc_shift_def = @allocated shift(c_central, 𝛈₀, ref_vol)
println("  shift       (default: allocates StaticNgon + MVector) → $alloc_shift_def bytes")

alloc_smeasure_ngon = @allocated smeasure(workspace)
println("  smeasure(StaticNgon)  (shoelace, no alloc expected)  → $alloc_smeasure_ngon bytes")

alloc_smeasure_phi_pre = @allocated smeasure(Φ, c_central; workspaces=smeasure_workspaces, gl_data=gl_data)
println("  smeasure(Φ, c)  (preallocated workspaces + gl_data) → $alloc_smeasure_phi_pre bytes")

alloc_smeasure_phi_def = @allocated smeasure(Φ, c_central)
println("  smeasure(Φ, c)  (default: allocates 2×StaticNgon + MVector{N,Bool}) → $alloc_smeasure_phi_def bytes")

alloc_costfun = @allocated GeometricVOF.lvira_costfun(
    p_shift, cs, αs_local, cmeasures_local, c_central; workspace=workspace)
println("  lvira_costfun  (preallocated workspace)             → $alloc_costfun bytes")

alloc_donating = @allocated donating_region!(dr_out, seg, dr_velo, dt)
println("  donating_region!  (mutating, no α constraint)       → $alloc_donating bytes")

alloc_reconstruct = @allocated reconstruct(
    p0, α_central, c_central, αs_local, cs, cmeasures_local;
    workspace=workspace, shift_workspace=shift_workspace)
println("  reconstruct  (preallocated workspace + shift_wsp)   → $alloc_reconstruct bytes")

alloc_reconstruct_def = @allocated reconstruct(p0, α_central, c_central, αs_local, cs, cmeasures_local)
println("  reconstruct  (default: allocates workspaces + cmeasures per call) → $alloc_reconstruct_def bytes")

# ── Profile.Allocs on one reconstruct call ────────────────────────────────────
subsection("Profile.Allocs — top allocation sites in reconstruct(…, preallocated)")
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=1 reconstruct(
    p0, α_central, c_central, αs_local, cs, cmeasures_local;
    workspace=workspace, shift_workspace=shift_workspace)
alloc_results = Profile.Allocs.fetch()
sorted_allocs = sort(alloc_results.allocs, by=a -> -a.size)
println("\n  Top allocation sites by size (bytes):")
println("  ", rpad("Size (B)", 14), rpad("Type", 42), "Location")
println("  ", "-"^82)
for a in first(sorted_allocs, 15)
    loc      = isempty(a.stacktrace) ? "unknown" : string(a.stacktrace[1])
    type_str = rpad(string(a.type), 42)
    println("  ", rpad(string(a.size), 14), type_str, loc)
end

# ──────────────────────────────────────────────────────────────────────────────
# ❻  SECTION 5 — Type stability (JET.jl)
# ──────────────────────────────────────────────────────────────────────────────
section("SECTION 5 — Type stability (JET.jl)")

println("""
Using JET.jl for static analysis:
  @report_call — detects type errors and runtime dispatch in the call tree
  @report_opt  — detects suboptimal code: runtime dispatch and allocation sources

Expected findings:
  • smeasure(Φ, c): MVector{N,Bool} allocation in ngon_approx! (heap alloc per call)
  • smeasure(Φ, c): T == Bool runtime comparison (value-level branch on a type)
  • lvira_costfun: Segment + centroid allocations per neighbour cell
  • reconstruct (default args): StaticNgon + MVector allocation in default kwarg expressions
""")

function jet_report(label, result)
    issues = JET.get_reports(result)
    if isempty(issues)
        println("  ✓  $label — no issues found")
    else
        println("  ✗  $label — $(length(issues)) issue(s):")
        for r in first(issues, 8)
            println("       ", r)
        end
        length(issues) > 8 && println("       … ($(length(issues) - 8) more)")
    end
    return !isempty(issues)
end

subsection("@report_call — detect dispatch / type errors in call trees")

r_intersect = @report_call intersect!(workspace, c_central, p_shift)
has_intersect = jet_report("intersect!(workspace, c, p)", r_intersect)

r_shift = @report_call shift(c_central, 𝛈₀, ref_vol; workspace=workspace, shift_workspace=shift_workspace)
has_shift = jet_report("shift(c, 𝛈, αvol; preallocated)", r_shift)

r_smeasure_ngon = @report_call smeasure(workspace)
has_smeasure_ngon = jet_report("smeasure(StaticNgon)", r_smeasure_ngon)

r_smeasure_phi = @report_call smeasure(Φ, c_central; workspaces=smeasure_workspaces, gl_data=gl_data)
has_smeasure_phi = jet_report("smeasure(Φ, c; preallocated)", r_smeasure_phi)

r_costfun = @report_call GeometricVOF.lvira_costfun(
    p_shift, cs, αs_local, cmeasures_local, c_central; workspace=workspace)
has_costfun = jet_report("lvira_costfun(p, cs, αs, cmeasures, c; preallocated)", r_costfun)

r_donating = @report_call donating_region!(dr_out, seg, dr_velo, dt)
has_donating = jet_report("donating_region!(out, seg, velo, dt)", r_donating)

r_reconstruct = @report_call reconstruct(
    p0, α_central, c_central, αs_local, cs, cmeasures_local;
    workspace=workspace, shift_workspace=shift_workspace)
has_reconstruct = jet_report("reconstruct(p0, α, c, αs, cs, cmeas; preallocated)", r_reconstruct)

subsection("@report_opt — detect runtime dispatch and allocation sources")

o_intersect = @report_opt intersect!(workspace, c_central, p_shift)
jet_report("intersect! (opt)", o_intersect)

o_shift = @report_opt shift(c_central, 𝛈₀, ref_vol; workspace=workspace, shift_workspace=shift_workspace)
jet_report("shift (opt)", o_shift)

o_smeasure_ngon = @report_opt smeasure(workspace)
jet_report("smeasure(StaticNgon) (opt)", o_smeasure_ngon)

o_smeasure_phi = @report_opt smeasure(Φ, c_central; workspaces=smeasure_workspaces, gl_data=gl_data)
jet_report("smeasure(Φ, c) (opt)", o_smeasure_phi)

o_costfun = @report_opt GeometricVOF.lvira_costfun(
    p_shift, cs, αs_local, cmeasures_local, c_central; workspace=workspace)
jet_report("lvira_costfun (opt)", o_costfun)

o_donating = @report_opt donating_region!(dr_out, seg, dr_velo, dt)
jet_report("donating_region! (opt)", o_donating)

o_reconstruct = @report_opt reconstruct(
    p0, α_central, c_central, αs_local, cs, cmeasures_local;
    workspace=workspace, shift_workspace=shift_workspace)
jet_report("reconstruct (opt)", o_reconstruct)

# ──────────────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────────────
section("SUMMARY")
println("""
Allocations (bytes):
  intersect!  (preallocated)           → $alloc_intersect_pre   (target: 0)
  intersect   (default — allocating)   → $alloc_intersect_def
  shift       (preallocated)           → $alloc_shift_pre   (target: 0)
  shift       (default — allocating)   → $alloc_shift_def
  smeasure(StaticNgon)                 → $alloc_smeasure_ngon   (target: 0)
  smeasure(Φ, c)  (preallocated)       → $alloc_smeasure_phi_pre
  smeasure(Φ, c)  (default)            → $alloc_smeasure_phi_def
  lvira_costfun   (preallocated)       → $alloc_costfun
  donating_region!                     → $alloc_donating   (target: 0)
  reconstruct  (preallocated)          → $alloc_reconstruct
  reconstruct  (default)               → $alloc_reconstruct_def

JET issues found:
  intersect!         : $(has_intersect  ? "YES — see above" : "none")
  shift              : $(has_shift      ? "YES — see above" : "none")
  smeasure(StaticNgon): $(has_smeasure_ngon ? "YES — see above" : "none")
  smeasure(Φ, c)     : $(has_smeasure_phi  ? "YES — see above" : "none")
  lvira_costfun      : $(has_costfun    ? "YES — see above" : "none")
  donating_region!   : $(has_donating   ? "YES — see above" : "none")
  reconstruct        : $(has_reconstruct ? "YES — see above" : "none")

Throughput (median):
  reconstruct    (full LVIRA)   : $(BenchmarkTools.median(b_reconstruct))
  smeasure(Φ, c) (initialise)   : $(BenchmarkTools.median(b_smeasure_phi))
  shift                         : $(BenchmarkTools.median(b_shift))
  intersect!                    : $(BenchmarkTools.median(b_intersect))
  lvira_costfun                 : $(BenchmarkTools.median(b_costfun))
  donating_region!              : $(BenchmarkTools.median(b_donating))
  brent_min       (proxy)       : $(BenchmarkTools.median(b_brent))
  smeasure(StaticNgon)          : $(BenchmarkTools.median(b_smeasure_ngon))
""")
