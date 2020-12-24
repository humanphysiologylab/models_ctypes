#include "./maleckar.h"
#include "math.h"

/*
   There are a total of 70 entries in the algebraic variable array.
   There are a total of 30 entries in each of the rate and state variable arrays.
   There are a total of 51 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (second).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is R in component membrane (millijoule_per_mole_kelvin).
 * CONSTANTS[1] is T in component membrane (kelvin).
 * CONSTANTS[2] is F in component membrane (coulomb_per_mole).
 * CONSTANTS[3] is Cm in component membrane (nanoF).
 * ALGEBRAIC[0] is Q_tot in component membrane (millivolt).
 * ALGEBRAIC[37] is i_Na in component sodium_current (picoA).
 * ALGEBRAIC[41] is i_Ca_L in component L_type_Ca_channel (picoA).
 * ALGEBRAIC[44] is i_t in component Ca_independent_transient_outward_K_current (picoA).
 * ALGEBRAIC[45] is i_Kur in component ultra_rapid_K_current (picoA).
 * ALGEBRAIC[46] is i_K1 in component inward_rectifier (picoA).
 * ALGEBRAIC[49] is i_Kr in component delayed_rectifier_K_currents (picoA).
 * ALGEBRAIC[47] is i_Ks in component delayed_rectifier_K_currents (picoA).
 * ALGEBRAIC[50] is i_B_Na in component background_currents (picoA).
 * ALGEBRAIC[52] is i_B_Ca in component background_currents (picoA).
 * ALGEBRAIC[54] is i_NaK in component sodium_potassium_pump (picoA).
 * ALGEBRAIC[55] is i_CaP in component sarcolemmal_calcium_pump_current (picoA).
 * ALGEBRAIC[56] is i_NaCa in component Na_Ca_ion_exchanger_current (picoA).
 * ALGEBRAIC[57] is i_KACh in component ACh_dependent_K_current (picoA).
 * ALGEBRAIC[59] is I in component membrane (pA_per_nF).
 * ALGEBRAIC[24] is i_Stim in component membrane (pA_per_nF).
 * CONSTANTS[4] is stim_offset in component membrane (second).
 * CONSTANTS[5] is stim_period in component membrane (second).
 * CONSTANTS[6] is stim_duration in component membrane (second).
 * CONSTANTS[7] is stim_amplitude in component membrane (pA_per_nF).
 * ALGEBRAIC[1] is past in component membrane (second).
 * ALGEBRAIC[35] is E_Na in component sodium_current (millivolt).
 * CONSTANTS[8] is P_Na in component sodium_current (nanolitre_per_second).
 * STATES[1] is Na_c in component cleft_space_ion_concentrations (millimolar).
 * STATES[2] is Na_i in component intracellular_ion_concentrations (millimolar).
 * STATES[3] is m in component sodium_current_m_gate (dimensionless).
 * STATES[4] is h1 in component sodium_current_h1_gate (dimensionless).
 * STATES[5] is h2 in component sodium_current_h2_gate (dimensionless).
 * ALGEBRAIC[14] is m_infinity in component sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[2] is m_factor in component sodium_current_m_gate (dimensionless).
 * ALGEBRAIC[26] is tau_m in component sodium_current_m_gate (second).
 * ALGEBRAIC[3] is h_infinity in component sodium_current_h1_gate (dimensionless).
 * ALGEBRAIC[15] is h_factor in component sodium_current_h1_gate (dimensionless).
 * ALGEBRAIC[27] is tau_h1 in component sodium_current_h1_gate (second).
 * ALGEBRAIC[28] is tau_h2 in component sodium_current_h2_gate (second).
 * CONSTANTS[9] is g_Ca_L in component L_type_Ca_channel (nanoS).
 * CONSTANTS[10] is E_Ca_app in component L_type_Ca_channel (millivolt).
 * ALGEBRAIC[39] is f_Ca in component L_type_Ca_channel (dimensionless).
 * CONSTANTS[11] is k_Ca in component L_type_Ca_channel (millimolar).
 * STATES[6] is Ca_d in component intracellular_ion_concentrations (millimolar).
 * STATES[7] is d_L in component L_type_Ca_channel_d_L_gate (dimensionless).
 * STATES[8] is f_L1 in component L_type_Ca_channel_f_L1_gate (dimensionless).
 * STATES[9] is f_L2 in component L_type_Ca_channel_f_L2_gate (dimensionless).
 * ALGEBRAIC[4] is d_L_infinity in component L_type_Ca_channel_d_L_gate (dimensionless).
 * ALGEBRAIC[16] is d_L_factor in component L_type_Ca_channel_d_L_gate (dimensionless).
 * ALGEBRAIC[29] is tau_d_L in component L_type_Ca_channel_d_L_gate (second).
 * ALGEBRAIC[5] is f_L_infinity in component L_type_Ca_channel_f_L1_gate (dimensionless).
 * ALGEBRAIC[17] is f_L_factor in component L_type_Ca_channel_f_L1_gate (millivolt).
 * ALGEBRAIC[30] is tau_f_L1 in component L_type_Ca_channel_f_L1_gate (second).
 * ALGEBRAIC[31] is tau_f_L2 in component L_type_Ca_channel_f_L2_gate (second).
 * ALGEBRAIC[43] is E_K in component Ca_independent_transient_outward_K_current (millivolt).
 * CONSTANTS[12] is g_t in component Ca_independent_transient_outward_K_current (nanoS).
 * STATES[10] is K_c in component cleft_space_ion_concentrations (millimolar).
 * STATES[11] is K_i in component intracellular_ion_concentrations (millimolar).
 * STATES[12] is r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
 * STATES[13] is s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
 * ALGEBRAIC[18] is tau_r in component Ca_independent_transient_outward_K_current_r_gate (second).
 * ALGEBRAIC[6] is r_infinity in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
 * ALGEBRAIC[32] is tau_s in component Ca_independent_transient_outward_K_current_s_gate (second).
 * ALGEBRAIC[7] is s_infinity in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
 * ALGEBRAIC[19] is s_factor in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
 * CONSTANTS[13] is g_kur in component ultra_rapid_K_current (nanoS).
 * STATES[14] is a_ur in component ultra_rapid_K_current_aur_gate (dimensionless).
 * STATES[15] is i_ur in component ultra_rapid_K_current_iur_gate (dimensionless).
 * ALGEBRAIC[8] is a_ur_infinity in component ultra_rapid_K_current_aur_gate (dimensionless).
 * ALGEBRAIC[20] is tau_a_ur in component ultra_rapid_K_current_aur_gate (second).
 * ALGEBRAIC[9] is i_ur_infinity in component ultra_rapid_K_current_iur_gate (dimensionless).
 * ALGEBRAIC[21] is tau_i_ur in component ultra_rapid_K_current_iur_gate (second).
 * CONSTANTS[14] is g_K1 in component inward_rectifier (nanoS).
 * CONSTANTS[15] is g_Ks in component delayed_rectifier_K_currents (nanoS).
 * CONSTANTS[16] is g_Kr in component delayed_rectifier_K_currents (nanoS).
 * STATES[16] is n in component delayed_rectifier_K_currents_n_gate (dimensionless).
 * STATES[17] is pa in component delayed_rectifier_K_currents_pa_gate (dimensionless).
 * ALGEBRAIC[48] is pip in component delayed_rectifier_K_currents_pi_gate (dimensionless).
 * ALGEBRAIC[33] is tau_n in component delayed_rectifier_K_currents_n_gate (second).
 * ALGEBRAIC[10] is n_infinity in component delayed_rectifier_K_currents_n_gate (dimensionless).
 * ALGEBRAIC[22] is n_factor in component delayed_rectifier_K_currents_n_gate (dimensionless).
 * ALGEBRAIC[34] is tau_pa in component delayed_rectifier_K_currents_pa_gate (second).
 * ALGEBRAIC[23] is pa_factor in component delayed_rectifier_K_currents_pa_gate (dimensionless).
 * ALGEBRAIC[11] is p_a_infinity in component delayed_rectifier_K_currents_pa_gate (dimensionless).
 * CONSTANTS[17] is g_B_Na in component background_currents (nanoS).
 * CONSTANTS[18] is g_B_Ca in component background_currents (nanoS).
 * ALGEBRAIC[51] is E_Ca in component background_currents (millivolt).
 * STATES[18] is Ca_c in component cleft_space_ion_concentrations (millimolar).
 * STATES[19] is Ca_i in component intracellular_ion_concentrations (millimolar).
 * CONSTANTS[19] is K_NaK_K in component sodium_potassium_pump (millimolar).
 * CONSTANTS[20] is i_NaK_max in component sodium_potassium_pump (picoA).
 * CONSTANTS[21] is pow_K_NaK_Na_15 in component sodium_potassium_pump (millimolar15).
 * ALGEBRAIC[53] is pow_Na_i_15 in component sodium_potassium_pump (millimolar15).
 * CONSTANTS[22] is i_CaP_max in component sarcolemmal_calcium_pump_current (picoA).
 * CONSTANTS[23] is k_CaP in component sarcolemmal_calcium_pump_current (millimolar).
 * CONSTANTS[24] is K_NaCa in component Na_Ca_ion_exchanger_current (picoA_per_millimolar_4).
 * CONSTANTS[25] is d_NaCa in component Na_Ca_ion_exchanger_current (per_millimolar_4).
 * CONSTANTS[26] is gamma_Na in component Na_Ca_ion_exchanger_current (dimensionless).
 * CONSTANTS[27] is ACh in component ACh_dependent_K_current (millimolar).
 * CONSTANTS[28] is phi_Na_en in component intracellular_ion_concentrations (picoA).
 * CONSTANTS[29] is Vol_i in component intracellular_ion_concentrations (nanolitre).
 * CONSTANTS[30] is Vol_d in component intracellular_ion_concentrations (nanolitre).
 * ALGEBRAIC[58] is i_di in component intracellular_ion_concentrations (picoA).
 * CONSTANTS[31] is tau_di in component intracellular_ion_concentrations (second).
 * ALGEBRAIC[67] is i_up in component Ca_handling_by_the_SR (picoA).
 * ALGEBRAIC[66] is i_rel in component Ca_handling_by_the_SR (picoA).
 * ALGEBRAIC[63] is J_O in component intracellular_Ca_buffering (per_second).
 * STATES[20] is O_C in component intracellular_Ca_buffering (dimensionless).
 * STATES[21] is O_TC in component intracellular_Ca_buffering (dimensionless).
 * STATES[22] is O_TMgC in component intracellular_Ca_buffering (dimensionless).
 * STATES[23] is O_TMgMg in component intracellular_Ca_buffering (dimensionless).
 * STATES[24] is O in component intracellular_Ca_buffering (dimensionless).
 * ALGEBRAIC[60] is J_O_C in component intracellular_Ca_buffering (per_second).
 * ALGEBRAIC[61] is J_O_TC in component intracellular_Ca_buffering (per_second).
 * ALGEBRAIC[62] is J_O_TMgC in component intracellular_Ca_buffering (per_second).
 * ALGEBRAIC[12] is J_O_TMgMg in component intracellular_Ca_buffering (per_second).
 * CONSTANTS[32] is Mg_i in component intracellular_Ca_buffering (millimolar).
 * CONSTANTS[33] is Vol_c in component cleft_space_ion_concentrations (nanolitre).
 * CONSTANTS[34] is tau_Na in component cleft_space_ion_concentrations (second).
 * CONSTANTS[35] is tau_K in component cleft_space_ion_concentrations (second).
 * CONSTANTS[36] is tau_Ca in component cleft_space_ion_concentrations (second).
 * CONSTANTS[37] is Na_b in component cleft_space_ion_concentrations (millimolar).
 * CONSTANTS[38] is Ca_b in component cleft_space_ion_concentrations (millimolar).
 * CONSTANTS[39] is K_b in component cleft_space_ion_concentrations (millimolar).
 * ALGEBRAIC[68] is i_tr in component Ca_handling_by_the_SR (picoA).
 * CONSTANTS[40] is I_up_max in component Ca_handling_by_the_SR (picoA).
 * CONSTANTS[41] is k_cyca in component Ca_handling_by_the_SR (millimolar).
 * CONSTANTS[42] is k_srca in component Ca_handling_by_the_SR (millimolar).
 * CONSTANTS[43] is k_xcs in component Ca_handling_by_the_SR (dimensionless).
 * CONSTANTS[44] is alpha_rel in component Ca_handling_by_the_SR (picoA_per_millimolar).
 * STATES[25] is Ca_rel in component Ca_handling_by_the_SR (millimolar).
 * STATES[26] is Ca_up in component Ca_handling_by_the_SR (millimolar).
 * CONSTANTS[45] is Vol_up in component Ca_handling_by_the_SR (nanolitre).
 * CONSTANTS[46] is Vol_rel in component Ca_handling_by_the_SR (nanolitre).
 * ALGEBRAIC[40] is r_act in component Ca_handling_by_the_SR (per_second).
 * ALGEBRAIC[42] is r_inact in component Ca_handling_by_the_SR (per_second).
 * CONSTANTS[47] is r_recov in component Ca_handling_by_the_SR (per_second).
 * ALGEBRAIC[13] is r_Ca_d_term in component Ca_handling_by_the_SR (dimensionless).
 * ALGEBRAIC[25] is r_Ca_i_term in component Ca_handling_by_the_SR (dimensionless).
 * ALGEBRAIC[36] is r_Ca_d_factor in component Ca_handling_by_the_SR (dimensionless).
 * ALGEBRAIC[38] is r_Ca_i_factor in component Ca_handling_by_the_SR (dimensionless).
 * ALGEBRAIC[64] is i_rel_f2 in component Ca_handling_by_the_SR (dimensionless).
 * ALGEBRAIC[65] is i_rel_factor in component Ca_handling_by_the_SR (dimensionless).
 * STATES[27] is O_Calse in component Ca_handling_by_the_SR (dimensionless).
 * ALGEBRAIC[69] is J_O_Calse in component Ca_handling_by_the_SR (per_second).
 * STATES[28] is F1 in component Ca_handling_by_the_SR (dimensionless).
 * STATES[29] is F2 in component Ca_handling_by_the_SR (dimensionless).
 * CONSTANTS[48] is tau_tr in component Ca_handling_by_the_SR (second).
 * CONSTANTS[49] is k_rel_i in component Ca_handling_by_the_SR (millimolar).
 * CONSTANTS[50] is k_rel_d in component Ca_handling_by_the_SR (millimolar).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[3] is d/dt m in component sodium_current_m_gate (dimensionless).
 * RATES[4] is d/dt h1 in component sodium_current_h1_gate (dimensionless).
 * RATES[5] is d/dt h2 in component sodium_current_h2_gate (dimensionless).
 * RATES[7] is d/dt d_L in component L_type_Ca_channel_d_L_gate (dimensionless).
 * RATES[8] is d/dt f_L1 in component L_type_Ca_channel_f_L1_gate (dimensionless).
 * RATES[9] is d/dt f_L2 in component L_type_Ca_channel_f_L2_gate (dimensionless).
 * RATES[12] is d/dt r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
 * RATES[13] is d/dt s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
 * RATES[14] is d/dt a_ur in component ultra_rapid_K_current_aur_gate (dimensionless).
 * RATES[15] is d/dt i_ur in component ultra_rapid_K_current_iur_gate (dimensionless).
 * RATES[16] is d/dt n in component delayed_rectifier_K_currents_n_gate (dimensionless).
 * RATES[17] is d/dt pa in component delayed_rectifier_K_currents_pa_gate (dimensionless).
 * RATES[11] is d/dt K_i in component intracellular_ion_concentrations (millimolar).
 * RATES[2] is d/dt Na_i in component intracellular_ion_concentrations (millimolar).
 * RATES[19] is d/dt Ca_i in component intracellular_ion_concentrations (millimolar).
 * RATES[6] is d/dt Ca_d in component intracellular_ion_concentrations (millimolar).
 * RATES[20] is d/dt O_C in component intracellular_Ca_buffering (dimensionless).
 * RATES[21] is d/dt O_TC in component intracellular_Ca_buffering (dimensionless).
 * RATES[22] is d/dt O_TMgC in component intracellular_Ca_buffering (dimensionless).
 * RATES[23] is d/dt O_TMgMg in component intracellular_Ca_buffering (dimensionless).
 * RATES[24] is d/dt O in component intracellular_Ca_buffering (dimensionless).
 * RATES[18] is d/dt Ca_c in component cleft_space_ion_concentrations (millimolar).
 * RATES[10] is d/dt K_c in component cleft_space_ion_concentrations (millimolar).
 * RATES[1] is d/dt Na_c in component cleft_space_ion_concentrations (millimolar).
 * RATES[28] is d/dt F1 in component Ca_handling_by_the_SR (dimensionless).
 * RATES[29] is d/dt F2 in component Ca_handling_by_the_SR (dimensionless).
 * RATES[27] is d/dt O_Calse in component Ca_handling_by_the_SR (dimensionless).
 * RATES[26] is d/dt Ca_up in component Ca_handling_by_the_SR (millimolar).
 * RATES[25] is d/dt Ca_rel in component Ca_handling_by_the_SR (millimolar).
 */


 double calculate_ICaL_GHK(double* CONSTANTS, double* RATES,
                           double* STATES, double* ALGEBRAIC) {

     double S_V = STATES[0];

     double A_d_L_infinity = (S_V <= 30.0) ? 1.565608783295298 * exp(-1.56167546 * exp(-0.04160361 * S_V)) : 1.;
     ALGEBRAIC[4]          = A_d_L_infinity;

     double A_tau_d_L      = ALGEBRAIC[29];
     RATES[7]              = (ALGEBRAIC[4] - STATES[7]) / ALGEBRAIC[29]; // R->d_L = (A->d_L_infinity - S->d_L) / A->tau_d_L;

     double A_f_L_infinity = ALGEBRAIC[5];
     double A_tau_f_L1     = ALGEBRAIC[30];
     // R->f_L1 = (A->f_L_infinity - S->f_L1) / A->tau_f_L1;

     double A_tau_f_L2     = ALGEBRAIC[31];
     // R->f_L2 = (A->f_L_infinity - S->f_L2) / A->tau_f_L2;

     double A_f_Ca         = ALGEBRAIC[39];

     double gamma_ca_i = 0.63;
     double gamma_ca_o = 0.63;

     double Ca_d  = STATES[6];
     double Ca_c  = STATES[18];

     double R = CONSTANTS[0];
     double T = CONSTANTS[1];
     double F = CONSTANTS[2];

     double phi_CaL = 4.0 * (S_V * F * F / (R * T)) *
                      (gamma_ca_i * Ca_d * exp(2.0 * (S_V * F / (R * T))) - gamma_ca_o * Ca_c) /
                      (exp(2.0 * (S_V * F / (R * T))) - 1.0);

     double g_Ca_L = CONSTANTS[9];

     double S_d_L  = STATES[7];
     double S_f_L1 = STATES[8];
     double S_f_L2 = STATES[9];

     double ICaL_GHK = 0.004 * (g_Ca_L / 6.75) * phi_CaL * S_d_L *
                       (A_f_Ca * S_f_L1 + (1 - A_f_Ca) * S_f_L2);

     return ICaL_GHK;
}


void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC) {
ALGEBRAIC[12] =  2000.00*CONSTANTS[32]*((1.00000 - STATES[22]) - STATES[23]) -  666.000*STATES[23];
RATES[23] = ALGEBRAIC[12];
ALGEBRAIC[18] =  0.00350000*exp((( - STATES[0]*STATES[0])/30.0000)/30.0000)+0.00150000;
ALGEBRAIC[6] = 1.00000/(1.00000+exp((STATES[0] - 1.00000)/- 11.0000));
RATES[12] = (ALGEBRAIC[6] - STATES[12])/ALGEBRAIC[18];
ALGEBRAIC[8] = 1.00000/(1.00000+exp(- (STATES[0]+6.00000)/8.60000));
ALGEBRAIC[20] = 0.00900000/(1.00000+exp((STATES[0]+5.00000)/12.0000))+0.000500000;
RATES[14] = (ALGEBRAIC[8] - STATES[14])/ALGEBRAIC[20];
ALGEBRAIC[9] = 1.00000/(1.00000+exp((STATES[0]+7.50000)/10.0000));
ALGEBRAIC[21] = 0.590000/(1.00000+exp((STATES[0]+60.0000)/10.0000))+3.05000;
RATES[15] = (ALGEBRAIC[9] - STATES[15])/ALGEBRAIC[21];
ALGEBRAIC[14] = 1.00000/(1.00000+exp((STATES[0]+27.1200)/- 8.21000));
ALGEBRAIC[2] = (STATES[0]+25.5700)/28.8000;
ALGEBRAIC[26] =  4.20000e-05*exp( - ALGEBRAIC[2]*ALGEBRAIC[2])+2.40000e-05;
RATES[3] = (ALGEBRAIC[14] - STATES[3])/ALGEBRAIC[26];
ALGEBRAIC[3] = 1.00000/(1.00000+exp((STATES[0]+63.6000)/5.30000));
ALGEBRAIC[15] = 1.00000/(1.00000+exp((STATES[0]+35.1000)/3.20000));
ALGEBRAIC[27] =  0.0300000*ALGEBRAIC[15]+0.000300000;
RATES[4] = (ALGEBRAIC[3] - STATES[4])/ALGEBRAIC[27];
ALGEBRAIC[28] =  0.120000*ALGEBRAIC[15]+0.00300000;
RATES[5] = (ALGEBRAIC[3] - STATES[5])/ALGEBRAIC[28];
ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[0]+9.00000)/- 5.80000));
ALGEBRAIC[16] = (STATES[0]+35.0000)/30.0000;
ALGEBRAIC[29] =  0.00270000*exp( - ALGEBRAIC[16]*ALGEBRAIC[16])+0.00200000;
RATES[7] = (ALGEBRAIC[4] - STATES[7])/ALGEBRAIC[29];
ALGEBRAIC[5] = 1.00000/(1.00000+exp((STATES[0]+27.4000)/7.10000));
ALGEBRAIC[17] = STATES[0]+40.0000;
ALGEBRAIC[30] =  0.161000*exp((( - ALGEBRAIC[17]*ALGEBRAIC[17])/14.4000)/14.4000)+0.0100000;
RATES[8] = (ALGEBRAIC[5] - STATES[8])/ALGEBRAIC[30];
ALGEBRAIC[31] =  1.33230*exp((( - ALGEBRAIC[17]*ALGEBRAIC[17])/14.2000)/14.2000)+0.0626000;
RATES[9] = (ALGEBRAIC[5] - STATES[9])/ALGEBRAIC[31];
ALGEBRAIC[19] = (STATES[0]+52.4500)/15.8827;
ALGEBRAIC[32] =  0.0256350*exp( - ALGEBRAIC[19]*ALGEBRAIC[19])+0.0141400;
ALGEBRAIC[7] = 1.00000/(1.00000+exp((STATES[0]+40.5000)/11.5000));
RATES[13] = (ALGEBRAIC[7] - STATES[13])/ALGEBRAIC[32];
ALGEBRAIC[22] = (STATES[0] - 20.0000)/20.0000;
ALGEBRAIC[33] = 0.700000+ 0.400000*exp( - ALGEBRAIC[22]*ALGEBRAIC[22]);
ALGEBRAIC[10] = 1.00000/(1.00000+exp((STATES[0] - 19.9000)/- 12.7000));
RATES[16] = (ALGEBRAIC[10] - STATES[16])/ALGEBRAIC[33];
ALGEBRAIC[23] = (STATES[0]+20.1376)/22.1996;
ALGEBRAIC[34] = 0.0311800+ 0.217180*exp( - ALGEBRAIC[23]*ALGEBRAIC[23]);
ALGEBRAIC[11] = 1.00000/(1.00000+exp((STATES[0]+15.0000)/- 6.00000));
RATES[17] = (ALGEBRAIC[11] - STATES[17])/ALGEBRAIC[34];
ALGEBRAIC[13] = STATES[6]/(STATES[6]+CONSTANTS[50]);
ALGEBRAIC[36] =  ALGEBRAIC[13]*ALGEBRAIC[13]*ALGEBRAIC[13]*ALGEBRAIC[13];
ALGEBRAIC[25] = STATES[19]/(STATES[19]+CONSTANTS[49]);
ALGEBRAIC[38] =  ALGEBRAIC[25]*ALGEBRAIC[25]*ALGEBRAIC[25]*ALGEBRAIC[25];
ALGEBRAIC[40] =  203.800*(ALGEBRAIC[38]+ALGEBRAIC[36]);
RATES[28] =  CONSTANTS[47]*((1.00000 - STATES[28]) - STATES[29]) -  ALGEBRAIC[40]*STATES[28];
ALGEBRAIC[42] = 33.9600+ 339.600*ALGEBRAIC[38];
RATES[29] =  ALGEBRAIC[40]*STATES[28] -  ALGEBRAIC[42]*STATES[29];
ALGEBRAIC[43] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(STATES[10]/STATES[11]);
ALGEBRAIC[44] =  CONSTANTS[12]*STATES[12]*STATES[13]*(STATES[0] - ALGEBRAIC[43]);
ALGEBRAIC[45] =  CONSTANTS[13]*STATES[14]*STATES[15]*(STATES[0] - ALGEBRAIC[43]);
ALGEBRAIC[46] = ( CONSTANTS[14]*pow(STATES[10]/1.00000, 0.445700)*(STATES[0] - ALGEBRAIC[43]))/(1.00000+exp(( 1.50000*((STATES[0] - ALGEBRAIC[43])+3.60000)*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])));
ALGEBRAIC[48] = 1.00000/(1.00000+exp((STATES[0]+55.0000)/24.0000));
ALGEBRAIC[49] =  CONSTANTS[16]*STATES[17]*ALGEBRAIC[48]*(STATES[0] - ALGEBRAIC[43]);
ALGEBRAIC[47] =  CONSTANTS[15]*STATES[16]*(STATES[0] - ALGEBRAIC[43]);
ALGEBRAIC[53] = pow(STATES[2], 1.50000);
ALGEBRAIC[54] = ( (( (( CONSTANTS[20]*STATES[10])/(STATES[10]+CONSTANTS[19]))*ALGEBRAIC[53])/(ALGEBRAIC[53]+CONSTANTS[21]))*(STATES[0]+150.000))/(STATES[0]+200.000);
ALGEBRAIC[1] =  floor(VOI/CONSTANTS[5])*CONSTANTS[5];
ALGEBRAIC[24] = CONSTANTS[7]; // (VOI - ALGEBRAIC[1]>=CONSTANTS[4]&&VOI - ALGEBRAIC[1]<=CONSTANTS[4]+CONSTANTS[6] ? CONSTANTS[7] : 0.00000);
RATES[11] = - (((ALGEBRAIC[44]+ALGEBRAIC[45]+ALGEBRAIC[46]+ALGEBRAIC[47]+ALGEBRAIC[49]) -  2.00000*ALGEBRAIC[54])+ ALGEBRAIC[24]*CONSTANTS[3])/( CONSTANTS[29]*CONSTANTS[2]);
RATES[10] = (CONSTANTS[39] - STATES[10])/CONSTANTS[35]+((ALGEBRAIC[44]+ALGEBRAIC[45]+ALGEBRAIC[46]+ALGEBRAIC[47]+ALGEBRAIC[49]) -  2.00000*ALGEBRAIC[54])/( CONSTANTS[33]*CONSTANTS[2]);
ALGEBRAIC[35] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(STATES[1]/STATES[2]);
ALGEBRAIC[37] = ( (( CONSTANTS[8]*STATES[3]*STATES[3]*STATES[3]*( 0.900000*STATES[4]+ 0.100000*STATES[5])*STATES[1]*STATES[0]*CONSTANTS[2]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*(exp(( (STATES[0] - ALGEBRAIC[35])*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - 1.00000))/(exp(( STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - 1.00000);
ALGEBRAIC[50] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[35]);
ALGEBRAIC[56] = ( CONSTANTS[24]*( STATES[2]*STATES[2]*STATES[2]*STATES[18]*exp(( CONSTANTS[2]*STATES[0]*CONSTANTS[26])/( CONSTANTS[0]*CONSTANTS[1])) -  STATES[1]*STATES[1]*STATES[1]*STATES[19]*exp(( (CONSTANTS[26] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))))/(1.00000+ CONSTANTS[25]*( STATES[1]*STATES[1]*STATES[1]*STATES[19]+ STATES[2]*STATES[2]*STATES[2]*STATES[18]));
RATES[2] = - (ALGEBRAIC[37]+ALGEBRAIC[50]+ 3.00000*ALGEBRAIC[56]+ 3.00000*ALGEBRAIC[54]+CONSTANTS[28])/( CONSTANTS[29]*CONSTANTS[2]);
ALGEBRAIC[39] = STATES[6]/(STATES[6]+CONSTANTS[11]);

    //  TODO: make ICaL ( A[41] ) modification clear
    //  ALGEBRAIC[41] =  CONSTANTS[9]*STATES[7]*( ALGEBRAIC[39]*STATES[8]+ (1.00000 - ALGEBRAIC[39])*STATES[9])*(STATES[0] - CONSTANTS[10]);
    //  C[10] is E_Ca_app

    // Nernst
    // double E_Ca = ((CONSTANTS[0] * CONSTANTS[1]) / (2. * CONSTANTS[2])) * log(STATES[18] / STATES[6]);  //  C->R * C->T / (2 * C->F) * log(S->Ca_c / S->Ca_d)
    // ALGEBRAIC[41] = CONSTANTS[9] * STATES[7] * (ALGEBRAIC[39] * STATES[8] + (1. - ALGEBRAIC[39]) * STATES[9]) * (STATES[0] - E_Ca);

    // GHK
    ALGEBRAIC[41] = calculate_ICaL_GHK(CONSTANTS, RATES, STATES, ALGEBRAIC);

ALGEBRAIC[51] =  (( CONSTANTS[0]*CONSTANTS[1])/( 2.00000*CONSTANTS[2]))*log(STATES[18]/STATES[19]);
ALGEBRAIC[52] =  CONSTANTS[18]*(STATES[0] - ALGEBRAIC[51]);
ALGEBRAIC[55] = ( CONSTANTS[22]*STATES[19])/(STATES[19]+CONSTANTS[23]);
RATES[18] = (CONSTANTS[38] - STATES[18])/CONSTANTS[36]+((ALGEBRAIC[41]+ALGEBRAIC[52]+ALGEBRAIC[55]) -  2.00000*ALGEBRAIC[56])/( 2.00000*CONSTANTS[33]*CONSTANTS[2]);
RATES[1] = (CONSTANTS[37] - STATES[1])/CONSTANTS[34]+(ALGEBRAIC[37]+ALGEBRAIC[50]+ 3.00000*ALGEBRAIC[56]+ 3.00000*ALGEBRAIC[54]+CONSTANTS[28])/( CONSTANTS[33]*CONSTANTS[2]);
ALGEBRAIC[58] = ( (STATES[6] - STATES[19])*2.00000*CONSTANTS[30]*CONSTANTS[2])/CONSTANTS[31];
RATES[6] = - (ALGEBRAIC[41]+ALGEBRAIC[58])/( 2.00000*CONSTANTS[30]*CONSTANTS[2]);
ALGEBRAIC[57] =  (10.0000/(1.00000+( 9.13652*pow(1.00000, 0.477811))/pow(CONSTANTS[27], 0.477811)))*(0.0517000+0.451600/(1.00000+exp((STATES[0]+59.5300)/17.1800)))*(STATES[0] - ALGEBRAIC[43])*CONSTANTS[3];
ALGEBRAIC[59] = (ALGEBRAIC[37]+ALGEBRAIC[41]+ALGEBRAIC[44]+ALGEBRAIC[45]+ALGEBRAIC[46]+ALGEBRAIC[49]+ALGEBRAIC[47]+ALGEBRAIC[50]+ALGEBRAIC[52]+ALGEBRAIC[54]+ALGEBRAIC[55]+ALGEBRAIC[56]+ALGEBRAIC[57])/CONSTANTS[3]+ALGEBRAIC[24];
RATES[0] =  - ALGEBRAIC[59]*1000.00;
ALGEBRAIC[60] =  200000.*STATES[19]*(1.00000 - STATES[20]) -  476.000*STATES[20];
RATES[20] = ALGEBRAIC[60];
ALGEBRAIC[61] =  78400.0*STATES[19]*(1.00000 - STATES[21]) -  392.000*STATES[21];
RATES[21] = ALGEBRAIC[61];
ALGEBRAIC[62] =  200000.*STATES[19]*((1.00000 - STATES[22]) - STATES[23]) -  6.60000*STATES[22];
RATES[22] = ALGEBRAIC[62];
ALGEBRAIC[63] =  0.0800000*ALGEBRAIC[61]+ 0.160000*ALGEBRAIC[62]+ 0.0450000*ALGEBRAIC[60];
RATES[24] = ALGEBRAIC[63];
ALGEBRAIC[67] = ( CONSTANTS[40]*(STATES[19]/CONSTANTS[41] - ( CONSTANTS[43]*CONSTANTS[43]*STATES[26])/CONSTANTS[42]))/((STATES[19]+CONSTANTS[41])/CONSTANTS[41]+( CONSTANTS[43]*(STATES[26]+CONSTANTS[42]))/CONSTANTS[42]);
ALGEBRAIC[64] = STATES[29]/(STATES[29]+0.250000);
ALGEBRAIC[65] =  ALGEBRAIC[64]*ALGEBRAIC[64];
ALGEBRAIC[66] =  CONSTANTS[44]*ALGEBRAIC[65]*(STATES[25] - STATES[19]);
RATES[19] = - ((ALGEBRAIC[52]+ALGEBRAIC[55]+ALGEBRAIC[67]) - (ALGEBRAIC[58]+ALGEBRAIC[66]+ 2.00000*ALGEBRAIC[56]))/( 2.00000*CONSTANTS[29]*CONSTANTS[2]) -  1.00000*ALGEBRAIC[63];
ALGEBRAIC[68] = ( (STATES[26] - STATES[25])*2.00000*CONSTANTS[46]*CONSTANTS[2])/CONSTANTS[48];
RATES[26] = (ALGEBRAIC[67] - ALGEBRAIC[68])/( 2.00000*CONSTANTS[45]*CONSTANTS[2]);
ALGEBRAIC[69] =  480.000*STATES[25]*(1.00000 - STATES[27]) -  400.000*STATES[27];
RATES[27] = ALGEBRAIC[69];
RATES[25] = (ALGEBRAIC[68] - ALGEBRAIC[66])/( 2.00000*CONSTANTS[46]*CONSTANTS[2]) -  31.0000*ALGEBRAIC[69];
}


void
initConsts(double* CONSTANTS, double *STATES)
{
STATES[0] = -74.031982;
CONSTANTS[0] = 8314;
CONSTANTS[1] = 306.15;
CONSTANTS[2] = 96487;
CONSTANTS[3] = 50;
CONSTANTS[4] = 0;
CONSTANTS[5] = 1;
CONSTANTS[6] = 0.006;
CONSTANTS[7] = -15;
CONSTANTS[8] = 0.0018;
STATES[1] = 130.022096;
STATES[2] = 8.516766;
STATES[3] = 0.003289;
STATES[4] = 0.877202;
STATES[5] = 0.873881;
CONSTANTS[9] = 6.75;
CONSTANTS[10] = 60;
CONSTANTS[11] = 0.025;
STATES[6] = 7.1e-5;
STATES[7] = 0.000014;
STATES[8] = 0.998597;
STATES[9] = 0.998586;
CONSTANTS[12] = 8.25;
STATES[10] = 5.560224;
STATES[11] = 129.485991;
STATES[12] = 0.001089;
STATES[13] = 0.948597;
CONSTANTS[13] = 2.25;
STATES[14] = 0.000367;
STATES[15] = 0.96729;
CONSTANTS[14] = 3.1;
CONSTANTS[15] = 1;
CONSTANTS[16] = 0.5;
STATES[16] = 0.004374;
STATES[17] = 0.000053;
CONSTANTS[17] = 0.060599;
CONSTANTS[18] = 0.078681;
STATES[18] = 1.815768;
STATES[19] = 6.5e-5;
CONSTANTS[19] = 1;
CONSTANTS[20] = 68.55;
CONSTANTS[21] = 36.4829;
CONSTANTS[22] = 4;
CONSTANTS[23] = 0.0002;
CONSTANTS[24] = 0.0374842;
CONSTANTS[25] = 0.0003;
CONSTANTS[26] = 0.45;
CONSTANTS[27] = 1e-24;
CONSTANTS[28] = 0;
CONSTANTS[29] = 0.005884;
CONSTANTS[30] = 0.00011768;
CONSTANTS[31] = 0.01;
STATES[20] = 0.026766;
STATES[21] = 0.012922;
STATES[22] = 0.190369;
STATES[23] = 0.714463;
STATES[24] = 1.38222;
CONSTANTS[32] = 2.5;
CONSTANTS[33] = 0.000800224;
CONSTANTS[34] = 14.3;
CONSTANTS[35] = 10;
CONSTANTS[36] = 24.7;
CONSTANTS[37] = 130;
CONSTANTS[38] = 1.8;
CONSTANTS[39] = 5.4;
CONSTANTS[40] = 2800;
CONSTANTS[41] = 0.0003;
CONSTANTS[42] = 0.5;
CONSTANTS[43] = 0.4;
CONSTANTS[44] = 200000;
STATES[25] = 0.632613;
STATES[26] = 0.649195;
CONSTANTS[45] = 0.0003969;
CONSTANTS[46] = 0.0000441;
CONSTANTS[47] = 0.815;
STATES[27] = 0.431547;
STATES[28] = 0.470055;
STATES[29] = 0.002814;
CONSTANTS[48] = 0.01;
CONSTANTS[49] = 0.0003;
CONSTANTS[50] = 0.003;
}


void
computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[12] =  2000.00*CONSTANTS[32]*((1.00000 - STATES[22]) - STATES[23]) -  666.000*STATES[23];
ALGEBRAIC[18] =  0.00350000*exp((( - STATES[0]*STATES[0])/30.0000)/30.0000)+0.00150000;
ALGEBRAIC[6] = 1.00000/(1.00000+exp((STATES[0] - 1.00000)/- 11.0000));
ALGEBRAIC[8] = 1.00000/(1.00000+exp(- (STATES[0]+6.00000)/8.60000));
ALGEBRAIC[20] = 0.00900000/(1.00000+exp((STATES[0]+5.00000)/12.0000))+0.000500000;
ALGEBRAIC[9] = 1.00000/(1.00000+exp((STATES[0]+7.50000)/10.0000));
ALGEBRAIC[21] = 0.590000/(1.00000+exp((STATES[0]+60.0000)/10.0000))+3.05000;
ALGEBRAIC[14] = 1.00000/(1.00000+exp((STATES[0]+27.1200)/- 8.21000));
ALGEBRAIC[2] = (STATES[0]+25.5700)/28.8000;
ALGEBRAIC[26] =  4.20000e-05*exp( - ALGEBRAIC[2]*ALGEBRAIC[2])+2.40000e-05;
ALGEBRAIC[3] = 1.00000/(1.00000+exp((STATES[0]+63.6000)/5.30000));
ALGEBRAIC[15] = 1.00000/(1.00000+exp((STATES[0]+35.1000)/3.20000));
ALGEBRAIC[27] =  0.0300000*ALGEBRAIC[15]+0.000300000;
ALGEBRAIC[28] =  0.120000*ALGEBRAIC[15]+0.00300000;
ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[0]+9.00000)/- 5.80000));
ALGEBRAIC[16] = (STATES[0]+35.0000)/30.0000;
ALGEBRAIC[29] =  0.00270000*exp( - ALGEBRAIC[16]*ALGEBRAIC[16])+0.00200000;
ALGEBRAIC[5] = 1.00000/(1.00000+exp((STATES[0]+27.4000)/7.10000));
ALGEBRAIC[17] = STATES[0]+40.0000;
ALGEBRAIC[30] =  0.161000*exp((( - ALGEBRAIC[17]*ALGEBRAIC[17])/14.4000)/14.4000)+0.0100000;
ALGEBRAIC[31] =  1.33230*exp((( - ALGEBRAIC[17]*ALGEBRAIC[17])/14.2000)/14.2000)+0.0626000;
ALGEBRAIC[19] = (STATES[0]+52.4500)/15.8827;
ALGEBRAIC[32] =  0.0256350*exp( - ALGEBRAIC[19]*ALGEBRAIC[19])+0.0141400;
ALGEBRAIC[7] = 1.00000/(1.00000+exp((STATES[0]+40.5000)/11.5000));
ALGEBRAIC[22] = (STATES[0] - 20.0000)/20.0000;
ALGEBRAIC[33] = 0.700000+ 0.400000*exp( - ALGEBRAIC[22]*ALGEBRAIC[22]);
ALGEBRAIC[10] = 1.00000/(1.00000+exp((STATES[0] - 19.9000)/- 12.7000));
ALGEBRAIC[23] = (STATES[0]+20.1376)/22.1996;
ALGEBRAIC[34] = 0.0311800+ 0.217180*exp( - ALGEBRAIC[23]*ALGEBRAIC[23]);
ALGEBRAIC[11] = 1.00000/(1.00000+exp((STATES[0]+15.0000)/- 6.00000));
ALGEBRAIC[13] = STATES[6]/(STATES[6]+CONSTANTS[50]);
ALGEBRAIC[36] =  ALGEBRAIC[13]*ALGEBRAIC[13]*ALGEBRAIC[13]*ALGEBRAIC[13];
ALGEBRAIC[25] = STATES[19]/(STATES[19]+CONSTANTS[49]);
ALGEBRAIC[38] =  ALGEBRAIC[25]*ALGEBRAIC[25]*ALGEBRAIC[25]*ALGEBRAIC[25];
ALGEBRAIC[40] =  203.800*(ALGEBRAIC[38]+ALGEBRAIC[36]);
ALGEBRAIC[42] = 33.9600+ 339.600*ALGEBRAIC[38];
ALGEBRAIC[43] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(STATES[10]/STATES[11]);
ALGEBRAIC[44] =  CONSTANTS[12]*STATES[12]*STATES[13]*(STATES[0] - ALGEBRAIC[43]);
ALGEBRAIC[45] =  CONSTANTS[13]*STATES[14]*STATES[15]*(STATES[0] - ALGEBRAIC[43]);
ALGEBRAIC[46] = ( CONSTANTS[14]*pow(STATES[10]/1.00000, 0.445700)*(STATES[0] - ALGEBRAIC[43]))/(1.00000+exp(( 1.50000*((STATES[0] - ALGEBRAIC[43])+3.60000)*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])));
ALGEBRAIC[48] = 1.00000/(1.00000+exp((STATES[0]+55.0000)/24.0000));
ALGEBRAIC[49] =  CONSTANTS[16]*STATES[17]*ALGEBRAIC[48]*(STATES[0] - ALGEBRAIC[43]);
ALGEBRAIC[47] =  CONSTANTS[15]*STATES[16]*(STATES[0] - ALGEBRAIC[43]);
ALGEBRAIC[53] = pow(STATES[2], 1.50000);
ALGEBRAIC[54] = ( (( (( CONSTANTS[20]*STATES[10])/(STATES[10]+CONSTANTS[19]))*ALGEBRAIC[53])/(ALGEBRAIC[53]+CONSTANTS[21]))*(STATES[0]+150.000))/(STATES[0]+200.000);
ALGEBRAIC[1] =  floor(VOI/CONSTANTS[5])*CONSTANTS[5];
ALGEBRAIC[24] = (VOI - ALGEBRAIC[1]>=CONSTANTS[4]&&VOI - ALGEBRAIC[1]<=CONSTANTS[4]+CONSTANTS[6] ? CONSTANTS[7] : 0.00000);
ALGEBRAIC[35] =  (( CONSTANTS[0]*CONSTANTS[1])/CONSTANTS[2])*log(STATES[1]/STATES[2]);
ALGEBRAIC[37] = ( (( CONSTANTS[8]*STATES[3]*STATES[3]*STATES[3]*( 0.900000*STATES[4]+ 0.100000*STATES[5])*STATES[1]*STATES[0]*CONSTANTS[2]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))*(exp(( (STATES[0] - ALGEBRAIC[35])*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - 1.00000))/(exp(( STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1])) - 1.00000);
ALGEBRAIC[50] =  CONSTANTS[17]*(STATES[0] - ALGEBRAIC[35]);
ALGEBRAIC[56] = ( CONSTANTS[24]*( STATES[2]*STATES[2]*STATES[2]*STATES[18]*exp(( CONSTANTS[2]*STATES[0]*CONSTANTS[26])/( CONSTANTS[0]*CONSTANTS[1])) -  STATES[1]*STATES[1]*STATES[1]*STATES[19]*exp(( (CONSTANTS[26] - 1.00000)*STATES[0]*CONSTANTS[2])/( CONSTANTS[0]*CONSTANTS[1]))))/(1.00000+ CONSTANTS[25]*( STATES[1]*STATES[1]*STATES[1]*STATES[19]+ STATES[2]*STATES[2]*STATES[2]*STATES[18]));
ALGEBRAIC[39] = STATES[6]/(STATES[6]+CONSTANTS[11]);

  //  TODO: make ICaL ( A[41] ) modification clear
  //  ALGEBRAIC[41] =  CONSTANTS[9]*STATES[7]*( ALGEBRAIC[39]*STATES[8]+ (1.00000 - ALGEBRAIC[39])*STATES[9])*(STATES[0] - CONSTANTS[10]);
  //  C[10] is E_Ca_app

  // Nernst
  // double E_Ca = ((CONSTANTS[0] * CONSTANTS[1]) / (2. * CONSTANTS[2])) * log(STATES[18] / STATES[6]);  //  C->R * C->T / (2 * C->F) * log(S->Ca_c / S->Ca_d)
  // ALGEBRAIC[41] = CONSTANTS[9] * STATES[7] * (ALGEBRAIC[39] * STATES[8] + (1. - ALGEBRAIC[39]) * STATES[9]) * (STATES[0] - E_Ca);

  // GHK
  ALGEBRAIC[41] = calculate_ICaL_GHK(CONSTANTS, RATES, STATES, ALGEBRAIC);

ALGEBRAIC[51] =  (( CONSTANTS[0]*CONSTANTS[1])/( 2.00000*CONSTANTS[2]))*log(STATES[18]/STATES[19]);
ALGEBRAIC[52] =  CONSTANTS[18]*(STATES[0] - ALGEBRAIC[51]);
ALGEBRAIC[55] = ( CONSTANTS[22]*STATES[19])/(STATES[19]+CONSTANTS[23]);
ALGEBRAIC[58] = ( (STATES[6] - STATES[19])*2.00000*CONSTANTS[30]*CONSTANTS[2])/CONSTANTS[31];
ALGEBRAIC[57] =  (10.0000/(1.00000+( 9.13652*pow(1.00000, 0.477811))/pow(CONSTANTS[27], 0.477811)))*(0.0517000+0.451600/(1.00000+exp((STATES[0]+59.5300)/17.1800)))*(STATES[0] - ALGEBRAIC[43])*CONSTANTS[3];
ALGEBRAIC[59] = (ALGEBRAIC[37]+ALGEBRAIC[41]+ALGEBRAIC[44]+ALGEBRAIC[45]+ALGEBRAIC[46]+ALGEBRAIC[49]+ALGEBRAIC[47]+ALGEBRAIC[50]+ALGEBRAIC[52]+ALGEBRAIC[54]+ALGEBRAIC[55]+ALGEBRAIC[56]+ALGEBRAIC[57])/CONSTANTS[3]+ALGEBRAIC[24];
ALGEBRAIC[60] =  200000.*STATES[19]*(1.00000 - STATES[20]) -  476.000*STATES[20];
ALGEBRAIC[61] =  78400.0*STATES[19]*(1.00000 - STATES[21]) -  392.000*STATES[21];
ALGEBRAIC[62] =  200000.*STATES[19]*((1.00000 - STATES[22]) - STATES[23]) -  6.60000*STATES[22];
ALGEBRAIC[63] =  0.0800000*ALGEBRAIC[61]+ 0.160000*ALGEBRAIC[62]+ 0.0450000*ALGEBRAIC[60];
ALGEBRAIC[67] = ( CONSTANTS[40]*(STATES[19]/CONSTANTS[41] - ( CONSTANTS[43]*CONSTANTS[43]*STATES[26])/CONSTANTS[42]))/((STATES[19]+CONSTANTS[41])/CONSTANTS[41]+( CONSTANTS[43]*(STATES[26]+CONSTANTS[42]))/CONSTANTS[42]);
ALGEBRAIC[64] = STATES[29]/(STATES[29]+0.250000);
ALGEBRAIC[65] =  ALGEBRAIC[64]*ALGEBRAIC[64];
ALGEBRAIC[66] =  CONSTANTS[44]*ALGEBRAIC[65]*(STATES[25] - STATES[19]);
ALGEBRAIC[68] = ( (STATES[26] - STATES[25])*2.00000*CONSTANTS[46]*CONSTANTS[2])/CONSTANTS[48];
ALGEBRAIC[69] =  480.000*STATES[25]*(1.00000 - STATES[27]) -  400.000*STATES[27];
ALGEBRAIC[0] =  0.0500000*STATES[0];
}
