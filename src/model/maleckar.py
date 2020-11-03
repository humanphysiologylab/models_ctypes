import numpy as np
import pandas as pd
from numba import njit


size_algebraic = 70
size_states = 30
size_constants = 51


legend_states = pd.DataFrame([
    {'name': 'V', 'component': 'membrane', 'units': 'millivolt'},
    {'name': 'Na_c', 'component': 'cleft_space_ion_concentrations', 'units': 'millimolar'},
    {'name': 'Na_i', 'component': 'intracellular_ion_concentrations', 'units': 'millimolar'},
    {'name': 'm', 'component': 'sodium_current_m_gate', 'units': 'dimensionless'},
    {'name': 'h1', 'component': 'sodium_current_h1_gate', 'units': 'dimensionless'},
    {'name': 'h2', 'component': 'sodium_current_h2_gate', 'units': 'dimensionless'},
    {'name': 'Ca_d', 'component': 'intracellular_ion_concentrations', 'units': 'millimolar'},
    {'name': 'd_L', 'component': 'L_type_Ca_channel_d_L_gate', 'units': 'dimensionless'},
    {'name': 'f_L1', 'component': 'L_type_Ca_channel_f_L1_gate', 'units': 'dimensionless'},
    {'name': 'f_L2', 'component': 'L_type_Ca_channel_f_L2_gate', 'units': 'dimensionless'},
    {'name': 'K_c', 'component': 'cleft_space_ion_concentrations', 'units': 'millimolar'},
    {'name': 'K_i', 'component': 'intracellular_ion_concentrations', 'units': 'millimolar'},
    {'name': 'r', 'component': 'Ca_independent_transient_outward_K_current_r_gate', 'units': 'dimensionless'},
    {'name': 's', 'component': 'Ca_independent_transient_outward_K_current_s_gate', 'units': 'dimensionless'},
    {'name': 'a_ur', 'component': 'ultra_rapid_K_current_aur_gate', 'units': 'dimensionless'},
    {'name': 'i_ur', 'component': 'ultra_rapid_K_current_iur_gate', 'units': 'dimensionless'},
    {'name': 'n', 'component': 'delayed_rectifier_K_currents_n_gate', 'units': 'dimensionless'},
    {'name': 'pa', 'component': 'delayed_rectifier_K_currents_pa_gate', 'units': 'dimensionless'},
    {'name': 'Ca_c', 'component': 'cleft_space_ion_concentrations', 'units': 'millimolar'},
    {'name': 'Ca_i', 'component': 'intracellular_ion_concentrations', 'units': 'millimolar'},
    {'name': 'O_C', 'component': 'intracellular_Ca_buffering', 'units': 'dimensionless'},
    {'name': 'O_TC', 'component': 'intracellular_Ca_buffering', 'units': 'dimensionless'},
    {'name': 'O_TMgC', 'component': 'intracellular_Ca_buffering', 'units': 'dimensionless'},
    {'name': 'O_TMgMg', 'component': 'intracellular_Ca_buffering', 'units': 'dimensionless'},
    {'name': 'O', 'component': 'intracellular_Ca_buffering', 'units': 'dimensionless'},
    {'name': 'Ca_rel', 'component': 'Ca_handling_by_the_SR', 'units': 'millimolar'},
    {'name': 'Ca_up', 'component': 'Ca_handling_by_the_SR', 'units': 'millimolar'},
    {'name': 'O_Calse', 'component': 'Ca_handling_by_the_SR', 'units': 'dimensionless'},
    {'name': 'F1', 'component': 'Ca_handling_by_the_SR', 'units': 'dimensionless'},
    {'name': 'F2', 'component': 'Ca_handling_by_the_SR', 'units': 'dimensionless'}])

legend_constants = pd.DataFrame([
    {'name': 'R', 'component': 'membrane', 'units': 'millijoule_per_mole_kelvin'},
    {'name': 'T', 'component': 'membrane', 'units': 'kelvin'},
    {'name': 'F', 'component': 'membrane', 'units': 'coulomb_per_mole'},
    {'name': 'Cm', 'component': 'membrane', 'units': 'nanoF'},
    {'name': 'stim_offset', 'component': 'membrane', 'units': 'second'},
    {'name': 'stim_period', 'component': 'membrane', 'units': 'second'},
    {'name': 'stim_duration', 'component': 'membrane', 'units': 'second'},
    {'name': 'stim_amplitude', 'component': 'membrane', 'units': 'pA_per_nF'},
    {'name': 'P_Na', 'component': 'sodium_current', 'units': 'nanolitre_per_second'},
    {'name': 'g_Ca_L', 'component': 'L_type_Ca_channel', 'units': 'nanoS'},
    {'name': 'E_Ca_app', 'component': 'L_type_Ca_channel', 'units': 'millivolt'},
    {'name': 'k_Ca', 'component': 'L_type_Ca_channel', 'units': 'millimolar'},
    {'name': 'g_t', 'component': 'Ca_independent_transient_outward_K_current', 'units': 'nanoS'},
    {'name': 'g_kur', 'component': 'ultra_rapid_K_current', 'units': 'nanoS'},
    {'name': 'g_K1', 'component': 'inward_rectifier', 'units': 'nanoS'},
    {'name': 'g_Ks', 'component': 'delayed_rectifier_K_currents', 'units': 'nanoS'},
    {'name': 'g_Kr', 'component': 'delayed_rectifier_K_currents', 'units': 'nanoS'},
    {'name': 'g_B_Na', 'component': 'background_currents', 'units': 'nanoS'},
    {'name': 'g_B_Ca', 'component': 'background_currents', 'units': 'nanoS'},
    {'name': 'K_NaK_K', 'component': 'sodium_potassium_pump', 'units': 'millimolar'},
    {'name': 'i_NaK_max', 'component': 'sodium_potassium_pump', 'units': 'picoA'},
    {'name': 'pow_K_NaK_Na_15', 'component': 'sodium_potassium_pump', 'units': 'millimolar15'},
    {'name': 'i_CaP_max', 'component': 'sarcolemmal_calcium_pump_current', 'units': 'picoA'},
    {'name': 'k_CaP', 'component': 'sarcolemmal_calcium_pump_current', 'units': 'millimolar'},
    {'name': 'K_NaCa', 'component': 'Na_Ca_ion_exchanger_current', 'units': 'picoA_per_millimolar_4'},
    {'name': 'd_NaCa', 'component': 'Na_Ca_ion_exchanger_current', 'units': 'per_millimolar_4'},
    {'name': 'gamma_Na', 'component': 'Na_Ca_ion_exchanger_current', 'units': 'dimensionless'},
    {'name': 'ACh', 'component': 'ACh_dependent_K_current', 'units': 'millimolar'},
    {'name': 'phi_Na_en', 'component': 'intracellular_ion_concentrations', 'units': 'picoA'},
    {'name': 'Vol_i', 'component': 'intracellular_ion_concentrations', 'units': 'nanolitre'},
    {'name': 'Vol_d', 'component': 'intracellular_ion_concentrations', 'units': 'nanolitre'},
    {'name': 'tau_di', 'component': 'intracellular_ion_concentrations', 'units': 'second'},
    {'name': 'Mg_i', 'component': 'intracellular_Ca_buffering', 'units': 'millimolar'},
    {'name': 'Vol_c', 'component': 'cleft_space_ion_concentrations', 'units': 'nanolitre'},
    {'name': 'tau_Na', 'component': 'cleft_space_ion_concentrations', 'units': 'second'},
    {'name': 'tau_K', 'component': 'cleft_space_ion_concentrations', 'units': 'second'},
    {'name': 'tau_Ca', 'component': 'cleft_space_ion_concentrations', 'units': 'second'},
    {'name': 'Na_b', 'component': 'cleft_space_ion_concentrations', 'units': 'millimolar'},
    {'name': 'Ca_b', 'component': 'cleft_space_ion_concentrations', 'units': 'millimolar'},
    {'name': 'K_b', 'component': 'cleft_space_ion_concentrations', 'units': 'millimolar'},
    {'name': 'I_up_max', 'component': 'Ca_handling_by_the_SR', 'units': 'picoA'},
    {'name': 'k_cyca', 'component': 'Ca_handling_by_the_SR', 'units': 'millimolar'},
    {'name': 'k_srca', 'component': 'Ca_handling_by_the_SR', 'units': 'millimolar'},
    {'name': 'k_xcs', 'component': 'Ca_handling_by_the_SR', 'units': 'dimensionless'},
    {'name': 'alpha_rel', 'component': 'Ca_handling_by_the_SR', 'units': 'picoA_per_millimolar'},
    {'name': 'Vol_up', 'component': 'Ca_handling_by_the_SR', 'units': 'nanolitre'},
    {'name': 'Vol_rel', 'component': 'Ca_handling_by_the_SR', 'units': 'nanolitre'},
    {'name': 'r_recov', 'component': 'Ca_handling_by_the_SR', 'units': 'per_second'},
    {'name': 'tau_tr', 'component': 'Ca_handling_by_the_SR', 'units': 'second'},
    {'name': 'k_rel_i', 'component': 'Ca_handling_by_the_SR', 'units': 'millimolar'},
    {'name': 'k_rel_d', 'component': 'Ca_handling_by_the_SR', 'units': 'millimolar'}])

legend_algebraic = pd.DataFrame([
    {'name': 'Q_tot', 'component': 'membrane', 'units': 'millivolt'},
    {'name': 'past', 'component': 'membrane', 'units': 'second'},
    {'name': 'm_factor', 'component': 'sodium_current_m_gate', 'units': 'dimensionless'},
    {'name': 'h_infinity', 'component': 'sodium_current_h1_gate', 'units': 'dimensionless'},
    {'name': 'd_L_infinity', 'component': 'L_type_Ca_channel_d_L_gate', 'units': 'dimensionless'},
    {'name': 'f_L_infinity', 'component': 'L_type_Ca_channel_f_L1_gate', 'units': 'dimensionless'},
    {'name': 'r_infinity', 'component': 'Ca_independent_transient_outward_K_current_r_gate', 'units': 'dimensionless'},
    {'name': 's_infinity', 'component': 'Ca_independent_transient_outward_K_current_s_gate', 'units': 'dimensionless'},
    {'name': 'a_ur_infinity', 'component': 'ultra_rapid_K_current_aur_gate', 'units': 'dimensionless'},
    {'name': 'i_ur_infinity', 'component': 'ultra_rapid_K_current_iur_gate', 'units': 'dimensionless'},
    {'name': 'n_infinity', 'component': 'delayed_rectifier_K_currents_n_gate', 'units': 'dimensionless'},
    {'name': 'p_a_infinity', 'component': 'delayed_rectifier_K_currents_pa_gate', 'units': 'dimensionless'},
    {'name': 'J_O_TMgMg', 'component': 'intracellular_Ca_buffering', 'units': 'per_second'},
    {'name': 'r_Ca_d_term', 'component': 'Ca_handling_by_the_SR', 'units': 'dimensionless'},
    {'name': 'm_infinity', 'component': 'sodium_current_m_gate', 'units': 'dimensionless'},
    {'name': 'h_factor', 'component': 'sodium_current_h1_gate', 'units': 'dimensionless'},
    {'name': 'd_L_factor', 'component': 'L_type_Ca_channel_d_L_gate', 'units': 'dimensionless'},
    {'name': 'f_L_factor', 'component': 'L_type_Ca_channel_f_L1_gate', 'units': 'millivolt'},
    {'name': 'tau_r', 'component': 'Ca_independent_transient_outward_K_current_r_gate', 'units': 'second'},
    {'name': 's_factor', 'component': 'Ca_independent_transient_outward_K_current_s_gate', 'units': 'dimensionless'},
    {'name': 'tau_a_ur', 'component': 'ultra_rapid_K_current_aur_gate', 'units': 'second'},
    {'name': 'tau_i_ur', 'component': 'ultra_rapid_K_current_iur_gate', 'units': 'second'},
    {'name': 'n_factor', 'component': 'delayed_rectifier_K_currents_n_gate', 'units': 'dimensionless'},
    {'name': 'pa_factor', 'component': 'delayed_rectifier_K_currents_pa_gate', 'units': 'dimensionless'},
    {'name': 'i_Stim', 'component': 'membrane', 'units': 'pA_per_nF'},
    {'name': 'r_Ca_i_term', 'component': 'Ca_handling_by_the_SR', 'units': 'dimensionless'},
    {'name': 'tau_m', 'component': 'sodium_current_m_gate', 'units': 'second'},
    {'name': 'tau_h1', 'component': 'sodium_current_h1_gate', 'units': 'second'},
    {'name': 'tau_h2', 'component': 'sodium_current_h2_gate', 'units': 'second'},
    {'name': 'tau_d_L', 'component': 'L_type_Ca_channel_d_L_gate', 'units': 'second'},
    {'name': 'tau_f_L1', 'component': 'L_type_Ca_channel_f_L1_gate', 'units': 'second'},
    {'name': 'tau_f_L2', 'component': 'L_type_Ca_channel_f_L2_gate', 'units': 'second'},
    {'name': 'tau_s', 'component': 'Ca_independent_transient_outward_K_current_s_gate', 'units': 'second'},
    {'name': 'tau_n', 'component': 'delayed_rectifier_K_currents_n_gate', 'units': 'second'},
    {'name': 'tau_pa', 'component': 'delayed_rectifier_K_currents_pa_gate', 'units': 'second'},
    {'name': 'E_Na', 'component': 'sodium_current', 'units': 'millivolt'},
    {'name': 'r_Ca_d_factor', 'component': 'Ca_handling_by_the_SR', 'units': 'dimensionless'},
    {'name': 'i_Na', 'component': 'sodium_current', 'units': 'picoA'},
    {'name': 'r_Ca_i_factor', 'component': 'Ca_handling_by_the_SR', 'units': 'dimensionless'},
    {'name': 'f_Ca', 'component': 'L_type_Ca_channel', 'units': 'dimensionless'},
    {'name': 'r_act', 'component': 'Ca_handling_by_the_SR', 'units': 'per_second'},
    {'name': 'i_Ca_L', 'component': 'L_type_Ca_channel', 'units': 'picoA'},
    {'name': 'r_inact', 'component': 'Ca_handling_by_the_SR', 'units': 'per_second'},
    {'name': 'E_K', 'component': 'Ca_independent_transient_outward_K_current', 'units': 'millivolt'},
    {'name': 'i_t', 'component': 'Ca_independent_transient_outward_K_current', 'units': 'picoA'},
    {'name': 'i_Kur', 'component': 'ultra_rapid_K_current', 'units': 'picoA'},
    {'name': 'i_K1', 'component': 'inward_rectifier', 'units': 'picoA'},
    {'name': 'i_Ks', 'component': 'delayed_rectifier_K_currents', 'units': 'picoA'},
    {'name': 'pip', 'component': 'delayed_rectifier_K_currents_pi_gate', 'units': 'dimensionless'},
    {'name': 'i_Kr', 'component': 'delayed_rectifier_K_currents', 'units': 'picoA'},
    {'name': 'i_B_Na', 'component': 'background_currents', 'units': 'picoA'},
    {'name': 'E_Ca', 'component': 'background_currents', 'units': 'millivolt'},
    {'name': 'i_B_Ca', 'component': 'background_currents', 'units': 'picoA'},
    {'name': 'pow_Na_i_15', 'component': 'sodium_potassium_pump', 'units': 'millimolar15'},
    {'name': 'i_NaK', 'component': 'sodium_potassium_pump', 'units': 'picoA'},
    {'name': 'i_CaP', 'component': 'sarcolemmal_calcium_pump_current', 'units': 'picoA'},
    {'name': 'i_NaCa', 'component': 'Na_Ca_ion_exchanger_current', 'units': 'picoA'},
    {'name': 'i_KACh', 'component': 'ACh_dependent_K_current', 'units': 'picoA'},
    {'name': 'i_di', 'component': 'intracellular_ion_concentrations', 'units': 'picoA'},
    {'name': 'I', 'component': 'membrane', 'units': 'pA_per_nF'},
    {'name': 'J_O_C', 'component': 'intracellular_Ca_buffering', 'units': 'per_second'},
    {'name': 'J_O_TC', 'component': 'intracellular_Ca_buffering', 'units': 'per_second'},
    {'name': 'J_O_TMgC', 'component': 'intracellular_Ca_buffering', 'units': 'per_second'},
    {'name': 'J_O', 'component': 'intracellular_Ca_buffering', 'units': 'per_second'},
    {'name': 'i_rel_f2', 'component': 'Ca_handling_by_the_SR', 'units': 'dimensionless'},
    {'name': 'i_rel_factor', 'component': 'Ca_handling_by_the_SR', 'units': 'dimensionless'},
    {'name': 'i_rel', 'component': 'Ca_handling_by_the_SR', 'units': 'picoA'},
    {'name': 'i_up', 'component': 'Ca_handling_by_the_SR', 'units': 'picoA'},
    {'name': 'i_tr', 'component': 'Ca_handling_by_the_SR', 'units': 'picoA'},
    {'name': 'J_O_Calse', 'component': 'Ca_handling_by_the_SR', 'units': 'per_second'}])
    
    
assert(size_states == legend_states.shape[0])
assert(size_constants == legend_constants.shape[0])
assert(size_algebraic == legend_algebraic.shape[0])

legend = dict(states=legend_states,
              constants=legend_constants,
              algebraic=legend_algebraic)
    

def init_states_constants():
    
    S = np.zeros(size_states)
    C = np.zeros(size_constants)
    
    S[0] = -74.031982
    C[0] = 8314
    C[1] = 306.15
    C[2] = 96487
    C[3] = 50
    C[4] = 0
    C[5] = 1
    C[6] = 0.006
    C[7] = -15
    C[8] = 0.0018
    S[1] = 130.022096
    S[2] = 8.516766
    S[3] = 0.003289
    S[4] = 0.877202
    S[5] = 0.873881
    C[9] = 6.75
    C[10] = 60
    C[11] = 0.025
    S[6] = 7.1e-5
    S[7] = 0.000014
    S[8] = 0.998597
    S[9] = 0.998586
    C[12] = 8.25
    S[10] = 5.560224
    S[11] = 129.485991
    S[12] = 0.001089
    S[13] = 0.948597
    C[13] = 2.25
    S[14] = 0.000367
    S[15] = 0.96729
    C[14] = 3.1
    C[15] = 1
    C[16] = 0.5
    S[16] = 0.004374
    S[17] = 0.000053
    C[17] = 0.060599
    C[18] = 0.078681
    S[18] = 1.815768
    S[19] = 6.5e-5
    C[19] = 1
    C[20] = 68.55
    C[21] = 36.4829
    C[22] = 4
    C[23] = 0.0002
    C[24] = 0.0374842
    C[25] = 0.0003
    C[26] = 0.45
    C[27] = 1e-24
    C[28] = 0
    C[29] = 0.005884
    C[30] = 0.00011768
    C[31] = 0.01
    S[20] = 0.026766
    S[21] = 0.012922
    S[22] = 0.190369
    S[23] = 0.714463
    S[24] = 1.38222
    C[32] = 2.5
    C[33] = 0.000800224
    C[34] = 14.3
    C[35] = 10
    C[36] = 24.7
    C[37] = 130
    C[38] = 1.8
    C[39] = 5.4
    C[40] = 2800
    C[41] = 0.0003
    C[42] = 0.5
    C[43] = 0.4
    C[44] = 200000
    S[25] = 0.632613
    S[26] = 0.649195
    C[45] = 0.0003969
    C[46] = 0.0000441
    C[47] = 0.815
    S[27] = 0.431547
    S[28] = 0.470055
    S[29] = 0.002814
    C[48] = 0.01
    C[49] = 0.0003
    C[50] = 0.003
    return S, C


@njit
def compute_rates_algebraic(t, S, C, R=None, A=None):

    if R is None:
        R = np.zeros_like(S)
    if A is None:
        A = np.zeros(size_algebraic)
    
    A[12] = 2000.00*C[32]*((1.00000-S[22])-S[23])-666.000*S[23]
    R[23] = A[12]
    A[18] = 0.00350000*np.exp(((-S[0]*S[0])/30.0000)/30.0000)+0.00150000
    A[6] = 1.00000/(1.00000+np.exp((S[0]-1.00000)/-11.0000))
    R[12] = (A[6]-S[12])/A[18]
    A[8] = 1.00000/(1.00000+np.exp(-(S[0]+6.00000)/8.60000))
    A[20] = 0.00900000/(1.00000+np.exp((S[0]+5.00000)/12.0000))+0.000500000
    R[14] = (A[8]-S[14])/A[20]
    A[9] = 1.00000/(1.00000+np.exp((S[0]+7.50000)/10.0000))
    A[21] = 0.590000/(1.00000+np.exp((S[0]+60.0000)/10.0000))+3.05000
    R[15] = (A[9]-S[15])/A[21]
    A[14] = 1.00000/(1.00000+np.exp((S[0]+27.1200)/-8.21000))
    A[2] = (S[0]+25.5700)/28.8000
    A[26] = 4.20000e-05*np.exp(-A[2]*A[2])+2.40000e-05
    R[3] = (A[14]-S[3])/A[26]
    A[3] = 1.00000/(1.00000+np.exp((S[0]+63.6000)/5.30000))
    A[15] = 1.00000/(1.00000+np.exp((S[0]+35.1000)/3.20000))
    A[27] = 0.0300000*A[15]+0.000300000
    R[4] = (A[3]-S[4])/A[27]
    A[28] = 0.120000*A[15]+0.00300000
    R[5] = (A[3]-S[5])/A[28]
    A[4] = 1.00000/(1.00000+np.exp((S[0]+9.00000)/-5.80000))
    A[16] = (S[0]+35.0000)/30.0000
    A[29] = 0.00270000*np.exp(-A[16]*A[16])+0.00200000
    R[7] = (A[4]-S[7])/A[29]
    A[5] = 1.00000/(1.00000+np.exp((S[0]+27.4000)/7.10000))
    A[17] = S[0]+40.0000
    A[30] = 0.161000*np.exp(((-A[17]*A[17])/14.4000)/14.4000)+0.0100000
    R[8] = (A[5]-S[8])/A[30]
    A[31] = 1.33230*np.exp(((-A[17]*A[17])/14.2000)/14.2000)+0.0626000
    R[9] = (A[5]-S[9])/A[31]
    A[19] = (S[0]+52.4500)/15.8827
    A[32] = 0.0256350*np.exp(-A[19]*A[19])+0.0141400
    A[7] = 1.00000/(1.00000+np.exp((S[0]+40.5000)/11.5000))
    R[13] = (A[7]-S[13])/A[32]
    A[22] = (S[0]-20.0000)/20.0000
    A[33] = 0.700000+0.400000*np.exp(-A[22]*A[22])
    A[10] = 1.00000/(1.00000+np.exp((S[0]-19.9000)/-12.7000))
    R[16] = (A[10]-S[16])/A[33]
    A[23] = (S[0]+20.1376)/22.1996
    A[34] = 0.0311800+0.217180*np.exp(-A[23]*A[23])
    A[11] = 1.00000/(1.00000+np.exp((S[0]+15.0000)/-6.00000))
    R[17] = (A[11]-S[17])/A[34]
    A[13] = S[6]/(S[6]+C[50])
    A[36] = A[13]*A[13]*A[13]*A[13]
    A[25] = S[19]/(S[19]+C[49])
    A[38] = A[25]*A[25]*A[25]*A[25]
    A[40] = 203.800*(A[38]+A[36])
    R[28] = C[47]*((1.00000-S[28])-S[29])-A[40]*S[28]
    A[42] = 33.9600+339.600*A[38]
    R[29] = A[40]*S[28]-A[42]*S[29]
    A[43] = ((C[0]*C[1])/C[2])*np.log(S[10]/S[11])
    A[44] = C[12]*S[12]*S[13]*(S[0]-A[43])
    A[45] = C[13]*S[14]*S[15]*(S[0]-A[43])
    A[46] = (C[14]*(np.power(S[10]/1.00000, 0.445700))*(S[0]-A[43]))/(1.00000+np.exp((1.50000*((S[0]-A[43])+3.60000)*C[2])/(C[0]*C[1])))
    A[48] = 1.00000/(1.00000+np.exp((S[0]+55.0000)/24.0000))
    A[49] = C[16]*S[17]*A[48]*(S[0]-A[43])
    A[47] = C[15]*S[16]*(S[0]-A[43])
    A[53] = np.power(S[2], 1.50000)
    A[54] = (((((C[20]*S[10])/(S[10]+C[19]))*A[53])/(A[53]+C[21]))*(S[0]+150.000))/(S[0]+200.000)
    A[1] = np.floor(t/C[5])*C[5]
    A[24] = C[7] if ((t-A[1] >= C[4]) and (t-A[1] <= C[4]+C[6])) else 0
    R[11] = -(((A[44]+A[45]+A[46]+A[47]+A[49])-2.00000*A[54])+A[24]*C[3])/(C[29]*C[2])
    R[10] = (C[39]-S[10])/C[35]+((A[44]+A[45]+A[46]+A[47]+A[49])-2.00000*A[54])/(C[33]*C[2])
    A[35] = ((C[0]*C[1])/C[2])*np.log(S[1]/S[2])
    A[37] = (((C[8]*S[3]*S[3]*S[3]*(0.900000*S[4]+0.100000*S[5])*S[1]*S[0]*C[2]*C[2])/(C[0]*C[1]))*(np.exp(((S[0]-A[35])*C[2])/(C[0]*C[1]))-1.00000))/(np.exp((S[0]*C[2])/(C[0]*C[1]))-1.00000)
    A[50] = C[17]*(S[0]-A[35])
    A[56] = (C[24]*(S[2]*S[2]*S[2]*S[18]*np.exp((C[2]*S[0]*C[26])/(C[0]*C[1]))-S[1]*S[1]*S[1]*S[19]*np.exp(((C[26]-1.00000)*S[0]*C[2])/(C[0]*C[1]))))/(1.00000+C[25]*(S[1]*S[1]*S[1]*S[19]+S[2]*S[2]*S[2]*S[18]))
    R[2] = -(A[37]+A[50]+3.00000*A[56]+3.00000*A[54]+C[28])/(C[29]*C[2])
    A[39] = S[6]/(S[6]+C[11])

    #  TODO: make ICaL ( A[41] ) modification clear
    #  A[41] = C[9]*S[7]*(A[39]*S[8]+(1.00000-A[39])*S[9])*(S[0]-C[10])  # C[10] is E_Ca_app
    E_Ca = ((C[0] * C[1]) / (2. * C[2])) * np.log(S[18] / S[6])  # C->R * C->T / (2 * C->F) * log(S->Ca_c / S->Ca_d)
    A[41] = C[9] * S[7] * (A[39] * S[8] + (1. - A[39]) * S[9]) * (S[0] - E_Ca)

    A[51] = ((C[0]*C[1])/(2.00000*C[2]))*np.log(S[18]/S[19])
    A[52] = C[18]*(S[0]-A[51])
    A[55] = (C[22]*S[19])/(S[19]+C[23])
    R[18] = (C[38]-S[18])/C[36]+((A[41]+A[52]+A[55])-2.00000*A[56])/(2.00000*C[33]*C[2])
    R[1] = (C[37]-S[1])/C[34]+(A[37]+A[50]+3.00000*A[56]+3.00000*A[54]+C[28])/(C[33]*C[2])
    A[58] = ((S[6]-S[19])*2.00000*C[30]*C[2])/C[31]
    R[6] = -(A[41]+A[58])/(2.00000*C[30]*C[2])
    A[57] = (10.0000/(1.00000+(9.13652*(np.power(1.00000, 0.477811)))/(np.power(C[27], 0.477811))))*(0.0517000+0.451600/(1.00000+np.exp((S[0]+59.5300)/17.1800)))*(S[0]-A[43])*C[3]
    A[59] = (A[37]+A[41]+A[44]+A[45]+A[46]+A[49]+A[47]+A[50]+A[52]+A[54]+A[55]+A[56]+A[57])/C[3]+A[24]
    R[0] = -A[59]*1000.00
    A[60] = 200000.*S[19]*(1.00000-S[20])-476.000*S[20]
    R[20] = A[60]
    A[61] = 78400.0*S[19]*(1.00000-S[21])-392.000*S[21]
    R[21] = A[61]
    A[62] = 200000.*S[19]*((1.00000-S[22])-S[23])-6.60000*S[22]
    R[22] = A[62]
    A[63] = 0.0800000*A[61]+0.160000*A[62]+0.0450000*A[60]
    R[24] = A[63]
    A[67] = (C[40]*(S[19]/C[41]-(C[43]*C[43]*S[26])/C[42]))/((S[19]+C[41])/C[41]+(C[43]*(S[26]+C[42]))/C[42])
    A[64] = S[29]/(S[29]+0.250000)
    A[65] = A[64]*A[64]
    A[66] = C[44]*A[65]*(S[25]-S[19])
    R[19] = -((A[52]+A[55]+A[67])-(A[58]+A[66]+2.00000*A[56]))/(2.00000*C[29]*C[2])-1.00000*A[63]
    A[68] = ((S[26]-S[25])*2.00000*C[46]*C[2])/C[48]
    R[26] = (A[67]-A[68])/(2.00000*C[45]*C[2])
    A[69] = 480.000*S[25]*(1.00000-S[27])-400.000*S[27]
    R[27] = A[69]
    R[25] = (A[68]-A[66])/(2.00000*C[46]*C[2])-31.0000*A[69]
    return R, A
    
    
def compute_rates(t, S, C, R=None, A=None):
    return compute_rates_algebraic(t, S, C, R, A)[0]
