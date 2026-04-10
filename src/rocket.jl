# Goddard Rocket Problem

function rocket_model end

const rocket_h_0  = 1.0
const rocket_v_0  = 0.0
const rocket_m_0  = 1.0
const rocket_g_0  = 1.0
const rocket_T_c  = 3.5
const rocket_h_c  = 500.0
const rocket_v_c  = 620.0
const rocket_m_c  = 0.6
const rocket_c    = 0.5*sqrt(rocket_g_0 * rocket_h_0)
const rocket_m_f  = rocket_m_c * rocket_m_0
const rocket_D_c  = 0.5 * rocket_v_c * (rocket_m_0 / rocket_g_0)
const rocket_T_max = rocket_T_c * rocket_m_0 * rocket_g_0
