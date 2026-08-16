# COPSBenchmark.jl

[![Run tests](https://github.com/madsuite-org/COPSBenchmark.jl/actions/workflows/action.yml/badge.svg)](https://github.com/madsuite-org/COPSBenchmark.jl/actions/workflows/action.yml)
[![Documentation](https://github.com/madsuite-org/COPSBenchmark.jl/actions/workflows/docs.yml/badge.svg)](https://madnlp.github.io/COPSBenchmark.jl/dev/)
[![codecov](https://codecov.io/gh/madsuite-org/COPSBenchmark.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/madsuite-org/COPSBenchmark.jl)

Implementation of the [COPS benchmark](https://www.mcs.anl.gov/~more/cops/) for nonlinear optimization,
with backends for [JuMP](https://github.com/jump-dev/JuMP.jl) and [ExaModels](https://github.com/madsuite-org/ExaModels.jl).

## Quick start

```julia
using COPSBenchmark, JuMP, Ipopt

model = COPSBenchmark.rocket_model(COPSBenchmark.JuMPBackend(), 400)
set_optimizer(model, Ipopt.Optimizer)
set_silent(model)
optimize!(model)
```

```julia
using COPSBenchmark, ExaModels, NLPModelsIpopt

model = COPSBenchmark.rocket_model(COPSBenchmark.ExaModelsBackend(), 400)
result = ipopt(model; print_level=0)
```
