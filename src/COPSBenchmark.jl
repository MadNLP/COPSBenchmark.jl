module COPSBenchmark

using Random

abstract type AbstractModelerBackend end
struct JuMPBackend <: AbstractModelerBackend end
struct ExaModelsBackend <: AbstractModelerBackend end

# COPS Instances
function bearing_model end
function camshape_model end
function catmix_model end
function chain_model end
function channel_model end
function clnlbeam_model end
function dirichlet_model end
function elec_model end
function gasoil_model end
function glider_model end
function henon_model end
function lane_emden_model end
function marine_model end
function methanol_model end
function minsurf_model end
function pinene_model end
function polygon_model end
function robot_model end
function rocket_model end
function steering_model end
function tetra_duct12_model end
function tetra_duct15_model end
function tetra_duct20_model end
function tetra_foam5_model end
function tetra_gear_model end
function tetra_hook_model end
function torsion_model end
function triangle_deer_model end
function triangle_pacman_model end
function triangle_turtle_model end

function transition_state_model end

include("pde_utils.jl")

end # module COPSBenchmark
