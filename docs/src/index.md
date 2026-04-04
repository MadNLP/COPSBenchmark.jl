# COPSBenchmark.jl

Implementation of the [COPS benchmark](https://www.mcs.anl.gov/~more/cops/) for nonlinear optimization, with backends for [JuMP](https://github.com/jump-dev/JuMP.jl) and [ExaModels](https://github.com/exanauts/ExaModels.jl).

## Installation

```julia
using Pkg
Pkg.add("COPSBenchmark")
```

## Quick start

```julia
using COPSBenchmark, JuMP, Ipopt

# Build and solve the rocket model with JuMP
model = COPSBenchmark.rocket_model(COPSBenchmark.JuMPBackend(), 400)
set_optimizer(model, Ipopt.Optimizer)
set_silent(model)
optimize!(model)
```

```julia
using COPSBenchmark, ExaModels, NLPModelsIpopt

# Build and solve the rocket model with ExaModels
model = COPSBenchmark.rocket_model(COPSBenchmark.ExaModelsBackend(), 400)
result = ipopt(model; print_level=0)
```

## Available models

| Model | Function | Parameters |
|-------|----------|------------|
| Journal bearing | `bearing_model` | `nx, ny` |
| Cam shape | `camshape_model` | `n` |
| Catalyst mixing | `catmix_model` | `nh` |
| Hanging chain | `chain_model` | `n` |
| Channel flow | `channel_model` | `nh` |
| Electrons on sphere | `elec_model` | `np` |
| Gas-oil yield | `gasoil_model` | `nh` |
| Hang glider | `glider_model` | `nh` |
| Marine population | `marine_model` | `nh` |
| Methanol-to-hydrocarbons | `methanol_model` | `nh` |
| Minimal surface | `minsurf_model` | `nx, ny` |
| Isomerization of α-pinene | `pinene_model` | `nh` |
| Largest polygon | `polygon_model` | `n` |
| Robot arm | `robot_model` | `nh` |
| Goddard rocket | `rocket_model` | `nh` |
| Rocket steering | `steering_model` | `nh` |
| Tetrahedral mesh | `tetra_*_model` | — |
| Elastic-plastic torsion | `torsion_model` | `nx, ny` |
| Triangular mesh | `triangle_*_model` | — |
| Dirichlet | `dirichlet_model` | `nh` |
| Henon | `henon_model` | `nh` |
| Lane-Emden | `lane_emden_model` | `nh` |
