# Benchmarks

Run the reproducible 2D VOF suite from the repository root:

```sh
julia --project=. benchmark/benchmarks.jl
```

It measures allocating and workspace-based polygon clipping, half-space and
level-set area evaluation, volume-to-plane-shift inversion, and both
reconstruction APIs. Benchmark setup is outside the timed region; compare the
minimum time and allocation estimate on the same machine and Julia version.
