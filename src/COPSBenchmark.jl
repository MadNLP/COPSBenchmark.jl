module COPSBenchmark

using Random

abstract type AbstractModelerBackend end
struct JuMPBackend <: AbstractModelerBackend end
struct ExaModelsBackend <: AbstractModelerBackend end

# COPS Instances
function bearing_model end
function camshape_model end
function camshape_recipe end
function camshape_args end
function catmix_model end
function catmix_recipe end
function catmix_args end
function chain_model end
function chain_recipe end
function chain_args end
function channel_model end
function channel_recipe end
function channel_args end
function clnlbeam_model end
function dirichlet_model end
function elec_model end
function gasoil_model end
function gasoil_recipe end
function gasoil_args end
function glider_model end
function henon_model end
function lane_emden_model end
function marine_model end
function marine_recipe end
function marine_args end
function methanol_model end
function methanol_recipe end
function methanol_args end
function minsurf_model end
function pinene_model end
function pinene_recipe end
function pinene_args end
function polygon_model end
function robot_model end
function rocket_model end
function rocket_recipe end
function rocket_args end
function steering_model end
function steering_recipe end
function steering_args end
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
