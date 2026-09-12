# Benchmarks

Run the reproducible 2D VOF suite from the repository root:

```sh
julia --project=. benchmark/benchmarks.jl
```

It measures allocating and workspace-based planar/parabolic polygon clipping,
half-space, level-set, and parabolic area evaluation, parabolic moments,
volume-to-plane-shift inversion, and LVIRA, MOF, PMOF, and PLVIRA
reconstruction. Benchmark setup is outside the timed region; compare the
minimum time and allocation estimate on the same machine and Julia version.
