#include "./ohara_rudy.h"

#include "math.h"

/*
   There are a total of 200 entries in the algebraic variable array.
   There are a total of 49 entries in each of the rate and state variable arrays.
   There are a total of 206 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (millisecond).
 * CONSTANTS[0] is celltype in component environment (dimensionless).
 * CONSTANTS[1] is nao in component extracellular (millimolar).
 * CONSTANTS[2] is cao in component extracellular (millimolar).
 * CONSTANTS[3] is ko in component extracellular (millimolar).
 * CONSTANTS[4] is R in component physical_constants (joule_per_kilomole_kelvin).
 * CONSTANTS[5] is T in component physical_constants (kelvin).
 * CONSTANTS[6] is F in component physical_constants (coulomb_per_mole).
 * CONSTANTS[7] is zna in component physical_constants (dimensionless).
 * CONSTANTS[8] is zca in component physical_constants (dimensionless).
 * CONSTANTS[9] is zk in component physical_constants (dimensionless).
 * CONSTANTS[10] is L in component cell_geometry (centimeter).
 * CONSTANTS[11] is rad in component cell_geometry (centimeter).
 * CONSTANTS[162] is vcell in component cell_geometry (microliter).
 * CONSTANTS[177] is Ageo in component cell_geometry (centimeter_squared).
 * CONSTANTS[183] is Acap in component cell_geometry (centimeter_squared).
 * CONSTANTS[184] is vmyo in component cell_geometry (microliter).
 * CONSTANTS[185] is vnsr in component cell_geometry (microliter).
 * CONSTANTS[186] is vjsr in component cell_geometry (microliter).
 * CONSTANTS[187] is vss in component cell_geometry (microliter).
 * STATES[0] is v in component membrane (millivolt).
 * ALGEBRAIC[12] is vfrt in component membrane (dimensionless).
 * CONSTANTS[169] is ffrt in component membrane (coulomb_per_mole_millivolt).
 * CONSTANTS[149] is frt in component membrane (per_millivolt).
 * ALGEBRAIC[58] is INa in component INa (microA_per_microF).
 * ALGEBRAIC[60] is INaL in component INaL (microA_per_microF).
 * ALGEBRAIC[66] is Ito in component Ito (microA_per_microF).
 * ALGEBRAIC[83] is ICaL in component ICaL (microA_per_microF).
 * ALGEBRAIC[84] is ICaNa in component ICaL (microA_per_microF).
 * ALGEBRAIC[87] is ICaK in component ICaL (microA_per_microF).
 * ALGEBRAIC[90] is IKr in component IKr (microA_per_microF).
 * ALGEBRAIC[96] is IKs in component IKs (microA_per_microF).
 * ALGEBRAIC[98] is IK1 in component IK1 (microA_per_microF).
 * ALGEBRAIC[130] is INaCa_i in component INaCa_i (microA_per_microF).
 * ALGEBRAIC[160] is INaCa_ss in component INaCa_i (microA_per_microF).
 * ALGEBRAIC[179] is INaK in component INaK (microA_per_microF).
 * ALGEBRAIC[185] is INab in component INab (microA_per_microF).
 * ALGEBRAIC[181] is IKb in component IKb (microA_per_microF).
 * ALGEBRAIC[190] is IpCa in component IpCa (microA_per_microF).
 * ALGEBRAIC[189] is ICab in component ICab (microA_per_microF).
 * ALGEBRAIC[0] is Istim in component membrane (microA_per_microF).
 * CONSTANTS[12] is i_Stim_Start in component membrane (millisecond).
 * CONSTANTS[13] is i_Stim_End in component membrane (millisecond).
 * CONSTANTS[14] is i_Stim_Amplitude in component membrane (microA_per_microF).
 * CONSTANTS[15] is i_Stim_Period in component membrane (millisecond).
 * CONSTANTS[16] is i_Stim_PulseDuration in component membrane (millisecond).
 * CONSTANTS[17] is KmCaMK in component CaMK (millimolar).
 * CONSTANTS[18] is aCaMK in component CaMK (per_millimolar_per_millisecond).
 * CONSTANTS[19] is bCaMK in component CaMK (per_millisecond).
 * CONSTANTS[20] is CaMKo in component CaMK (dimensionless).
 * CONSTANTS[21] is KmCaM in component CaMK (millimolar).
 * ALGEBRAIC[36] is CaMKb in component CaMK (millimolar).
 * ALGEBRAIC[42] is CaMKa in component CaMK (millimolar).
 * STATES[1] is CaMKt in component CaMK (millimolar).
 * STATES[2] is cass in component intracellular_ions (millimolar).
 * CONSTANTS[22] is cmdnmax_b in component intracellular_ions (millimolar).
 * CONSTANTS[150] is cmdnmax in component intracellular_ions (millimolar).
 * CONSTANTS[23] is kmcmdn in component intracellular_ions (millimolar).
 * CONSTANTS[24] is trpnmax in component intracellular_ions (millimolar).
 * CONSTANTS[25] is kmtrpn in component intracellular_ions (millimolar).
 * CONSTANTS[26] is BSRmax in component intracellular_ions (millimolar).
 * CONSTANTS[27] is KmBSR in component intracellular_ions (millimolar).
 * CONSTANTS[28] is BSLmax in component intracellular_ions (millimolar).
 * CONSTANTS[29] is KmBSL in component intracellular_ions (millimolar).
 * CONSTANTS[30] is csqnmax in component intracellular_ions (millimolar).
 * CONSTANTS[31] is kmcsqn in component intracellular_ions (millimolar).
 * STATES[3] is nai in component intracellular_ions (millimolar).
 * STATES[4] is nass in component intracellular_ions (millimolar).
 * STATES[5] is ki in component intracellular_ions (millimolar).
 * STATES[6] is kss in component intracellular_ions (millimolar).
 * STATES[7] is cansr in component intracellular_ions (millimolar).
 * STATES[8] is cajsr in component intracellular_ions (millimolar).
 * STATES[9] is cai in component intracellular_ions (millimolar).
 * ALGEBRAIC[187] is JdiffNa in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[191] is Jdiff in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[198] is Jup in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[183] is JdiffK in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[193] is Jrel in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[199] is Jtr in component trans_flux (millimolar_per_millisecond).
 * ALGEBRAIC[44] is Bcai in component intracellular_ions (dimensionless).
 * ALGEBRAIC[48] is Bcajsr in component intracellular_ions (dimensionless).
 * ALGEBRAIC[46] is Bcass in component intracellular_ions (dimensionless).
 * CONSTANTS[32] is cm in component intracellular_ions (microF_per_centimeter_squared).
 * CONSTANTS[33] is PKNa in component reversal_potentials (dimensionless).
 * ALGEBRAIC[50] is ENa in component reversal_potentials (millivolt).
 * ALGEBRAIC[53] is EK in component reversal_potentials (millivolt).
 * ALGEBRAIC[54] is EKs in component reversal_potentials (millivolt).
 * ALGEBRAIC[1] is mss in component INa (dimensionless).
 * ALGEBRAIC[13] is tm in component INa (millisecond).
 * CONSTANTS[34] is mssV1 in component INa (millivolt).
 * CONSTANTS[35] is mssV2 in component INa (millivolt).
 * CONSTANTS[36] is mtV1 in component INa (millivolt).
 * CONSTANTS[37] is mtV2 in component INa (millivolt).
 * CONSTANTS[38] is mtD1 in component INa (dimensionless).
 * CONSTANTS[39] is mtD2 in component INa (dimensionless).
 * CONSTANTS[40] is mtV3 in component INa (millivolt).
 * CONSTANTS[41] is mtV4 in component INa (millivolt).
 * STATES[10] is m in component INa (dimensionless).
 * ALGEBRAIC[2] is hss in component INa (dimensionless).
 * ALGEBRAIC[14] is thf in component INa (millisecond).
 * ALGEBRAIC[15] is ths in component INa (millisecond).
 * CONSTANTS[42] is hssV1 in component INa (millivolt).
 * CONSTANTS[43] is hssV2 in component INa (millivolt).
 * CONSTANTS[151] is Ahs in component INa (dimensionless).
 * CONSTANTS[44] is Ahf in component INa (dimensionless).
 * STATES[11] is hf in component INa (dimensionless).
 * STATES[12] is hs in component INa (dimensionless).
 * ALGEBRAIC[55] is h in component INa (dimensionless).
 * CONSTANTS[45] is GNa in component INa (milliS_per_microF).
 * CONSTANTS[46] is shift_INa_inact in component INa (millivolt).
 * ALGEBRAIC[16] is jss in component INa (dimensionless).
 * ALGEBRAIC[27] is tj in component INa (millisecond).
 * STATES[13] is j in component INa (dimensionless).
 * ALGEBRAIC[28] is hssp in component INa (dimensionless).
 * ALGEBRAIC[37] is thsp in component INa (millisecond).
 * STATES[14] is hsp in component INa (dimensionless).
 * ALGEBRAIC[56] is hp in component INa (dimensionless).
 * ALGEBRAIC[38] is tjp in component INa (millisecond).
 * STATES[15] is jp in component INa (dimensionless).
 * ALGEBRAIC[57] is fINap in component INa (dimensionless).
 * ALGEBRAIC[29] is mLss in component INaL (dimensionless).
 * ALGEBRAIC[39] is tmL in component INaL (millisecond).
 * STATES[16] is mL in component INaL (dimensionless).
 * CONSTANTS[47] is thL in component INaL (millisecond).
 * ALGEBRAIC[3] is hLss in component INaL (dimensionless).
 * STATES[17] is hL in component INaL (dimensionless).
 * ALGEBRAIC[4] is hLssp in component INaL (dimensionless).
 * CONSTANTS[152] is thLp in component INaL (millisecond).
 * STATES[18] is hLp in component INaL (dimensionless).
 * CONSTANTS[48] is GNaL_b in component INaL (milliS_per_microF).
 * CONSTANTS[153] is GNaL in component INaL (milliS_per_microF).
 * ALGEBRAIC[59] is fINaLp in component INaL (dimensionless).
 * CONSTANTS[49] is Gto_b in component Ito (milliS_per_microF).
 * ALGEBRAIC[5] is ass in component Ito (dimensionless).
 * ALGEBRAIC[17] is ta in component Ito (millisecond).
 * STATES[19] is a in component Ito (dimensionless).
 * ALGEBRAIC[6] is iss in component Ito (dimensionless).
 * ALGEBRAIC[18] is delta_epi in component Ito (dimensionless).
 * ALGEBRAIC[30] is tiF_b in component Ito (millisecond).
 * ALGEBRAIC[40] is tiS_b in component Ito (millisecond).
 * ALGEBRAIC[43] is tiF in component Ito (millisecond).
 * ALGEBRAIC[45] is tiS in component Ito (millisecond).
 * ALGEBRAIC[61] is AiF in component Ito (dimensionless).
 * ALGEBRAIC[62] is AiS in component Ito (dimensionless).
 * STATES[20] is iF in component Ito (dimensionless).
 * STATES[21] is iS in component Ito (dimensionless).
 * ALGEBRAIC[63] is i in component Ito (dimensionless).
 * ALGEBRAIC[31] is assp in component Ito (dimensionless).
 * STATES[22] is ap in component Ito (dimensionless).
 * ALGEBRAIC[47] is dti_develop in component Ito (dimensionless).
 * ALGEBRAIC[49] is dti_recover in component Ito (dimensionless).
 * ALGEBRAIC[51] is tiFp in component Ito (millisecond).
 * ALGEBRAIC[52] is tiSp in component Ito (millisecond).
 * STATES[23] is iFp in component Ito (dimensionless).
 * STATES[24] is iSp in component Ito (dimensionless).
 * ALGEBRAIC[64] is ip in component Ito (dimensionless).
 * CONSTANTS[154] is Gto in component Ito (milliS_per_microF).
 * ALGEBRAIC[65] is fItop in component Ito (dimensionless).
 * CONSTANTS[50] is Kmn in component ICaL (millimolar).
 * CONSTANTS[51] is k2n in component ICaL (per_millisecond).
 * CONSTANTS[52] is PCa_b in component ICaL (dimensionless).
 * ALGEBRAIC[7] is dss in component ICaL (dimensionless).
 * STATES[25] is d in component ICaL (dimensionless).
 * ALGEBRAIC[8] is fss in component ICaL (dimensionless).
 * CONSTANTS[155] is Aff in component ICaL (dimensionless).
 * CONSTANTS[170] is Afs in component ICaL (dimensionless).
 * STATES[26] is ff in component ICaL (dimensionless).
 * STATES[27] is fs in component ICaL (dimensionless).
 * ALGEBRAIC[67] is f in component ICaL (dimensionless).
 * ALGEBRAIC[19] is fcass in component ICaL (dimensionless).
 * ALGEBRAIC[68] is Afcaf in component ICaL (dimensionless).
 * ALGEBRAIC[69] is Afcas in component ICaL (dimensionless).
 * STATES[28] is fcaf in component ICaL (dimensionless).
 * STATES[29] is fcas in component ICaL (dimensionless).
 * ALGEBRAIC[70] is fca in component ICaL (dimensionless).
 * STATES[30] is jca in component ICaL (dimensionless).
 * STATES[31] is ffp in component ICaL (dimensionless).
 * ALGEBRAIC[71] is fp in component ICaL (dimensionless).
 * STATES[32] is fcafp in component ICaL (dimensionless).
 * ALGEBRAIC[72] is fcap in component ICaL (dimensionless).
 * ALGEBRAIC[9] is km2n in component ICaL (per_millisecond).
 * ALGEBRAIC[20] is anca in component ICaL (dimensionless).
 * STATES[33] is nca in component ICaL (dimensionless).
 * ALGEBRAIC[75] is PhiCaL in component ICaL (dimensionless).
 * ALGEBRAIC[78] is PhiCaNa in component ICaL (dimensionless).
 * ALGEBRAIC[81] is PhiCaK in component ICaL (dimensionless).
 * CONSTANTS[156] is PCa in component ICaL (dimensionless).
 * CONSTANTS[171] is PCap in component ICaL (dimensionless).
 * CONSTANTS[172] is PCaNa in component ICaL (dimensionless).
 * CONSTANTS[173] is PCaK in component ICaL (dimensionless).
 * CONSTANTS[181] is PCaNap in component ICaL (dimensionless).
 * CONSTANTS[182] is PCaKp in component ICaL (dimensionless).
 * ALGEBRAIC[82] is fICaLp in component ICaL (dimensionless).
 * ALGEBRAIC[21] is td in component ICaL (millisecond).
 * ALGEBRAIC[22] is tff in component ICaL (millisecond).
 * ALGEBRAIC[23] is tfs in component ICaL (millisecond).
 * ALGEBRAIC[32] is tfcaf in component ICaL (millisecond).
 * ALGEBRAIC[33] is tfcas in component ICaL (millisecond).
 * CONSTANTS[157] is tjca in component ICaL (millisecond).
 * ALGEBRAIC[34] is tffp in component ICaL (millisecond).
 * ALGEBRAIC[41] is tfcafp in component ICaL (millisecond).
 * CONSTANTS[158] is v0 in component ICaL (millivolt).
 * ALGEBRAIC[73] is A_1 in component ICaL (dimensionless).
 * CONSTANTS[174] is B_1 in component ICaL (per_millivolt).
 * ALGEBRAIC[74] is U_1 in component ICaL (dimensionless).
 * ALGEBRAIC[76] is A_2 in component ICaL (dimensionless).
 * CONSTANTS[175] is B_2 in component ICaL (per_millivolt).
 * ALGEBRAIC[77] is U_2 in component ICaL (dimensionless).
 * ALGEBRAIC[79] is A_3 in component ICaL (dimensionless).
 * CONSTANTS[176] is B_3 in component ICaL (per_millivolt).
 * ALGEBRAIC[80] is U_3 in component ICaL (dimensionless).
 * CONSTANTS[53] is GKr_b in component IKr (milliS_per_microF).
 * STATES[34] is IC1 in component IKr (dimensionless).
 * STATES[35] is IC2 in component IKr (dimensionless).
 * STATES[36] is C1 in component IKr (dimensionless).
 * STATES[37] is C2 in component IKr (dimensionless).
 * STATES[38] is O in component IKr (dimensionless).
 * STATES[39] is IO in component IKr (dimensionless).
 * STATES[40] is IObound in component IKr (dimensionless).
 * STATES[41] is Obound in component IKr (dimensionless).
 * STATES[42] is Cbound in component IKr (dimensionless).
 * STATES[43] is D in component IKr (dimensionless).
 * CONSTANTS[159] is GKr in component IKr (milliS_per_microF).
 * CONSTANTS[54] is A1 in component IKr (per_millisecond).
 * CONSTANTS[55] is B1 in component IKr (per_millivolt).
 * CONSTANTS[56] is q1 in component IKr (dimensionless).
 * CONSTANTS[57] is A2 in component IKr (per_millisecond).
 * CONSTANTS[58] is B2 in component IKr (per_millivolt).
 * CONSTANTS[59] is q2 in component IKr (dimensionless).
 * CONSTANTS[60] is A3 in component IKr (per_millisecond).
 * CONSTANTS[61] is B3 in component IKr (per_millivolt).
 * CONSTANTS[62] is q3 in component IKr (dimensionless).
 * CONSTANTS[63] is A4 in component IKr (per_millisecond).
 * CONSTANTS[64] is B4 in component IKr (per_millivolt).
 * CONSTANTS[65] is q4 in component IKr (dimensionless).
 * CONSTANTS[66] is A11 in component IKr (per_millisecond).
 * CONSTANTS[67] is B11 in component IKr (per_millivolt).
 * CONSTANTS[68] is q11 in component IKr (dimensionless).
 * CONSTANTS[69] is A21 in component IKr (per_millisecond).
 * CONSTANTS[70] is B21 in component IKr (per_millivolt).
 * CONSTANTS[71] is q21 in component IKr (dimensionless).
 * CONSTANTS[72] is A31 in component IKr (per_millisecond).
 * CONSTANTS[73] is B31 in component IKr (per_millivolt).
 * CONSTANTS[74] is q31 in component IKr (dimensionless).
 * CONSTANTS[75] is A41 in component IKr (per_millisecond).
 * CONSTANTS[76] is B41 in component IKr (per_millivolt).
 * CONSTANTS[77] is q41 in component IKr (dimensionless).
 * CONSTANTS[78] is A51 in component IKr (per_millisecond).
 * CONSTANTS[79] is B51 in component IKr (per_millivolt).
 * CONSTANTS[80] is q51 in component IKr (dimensionless).
 * CONSTANTS[81] is A52 in component IKr (per_millisecond).
 * CONSTANTS[82] is B52 in component IKr (per_millivolt).
 * CONSTANTS[83] is q52 in component IKr (dimensionless).
 * CONSTANTS[84] is A53 in component IKr (per_millisecond).
 * CONSTANTS[85] is B53 in component IKr (per_millivolt).
 * CONSTANTS[86] is q53 in component IKr (dimensionless).
 * CONSTANTS[87] is A61 in component IKr (per_millisecond).
 * CONSTANTS[88] is B61 in component IKr (per_millivolt).
 * CONSTANTS[89] is q61 in component IKr (dimensionless).
 * CONSTANTS[90] is A62 in component IKr (per_millisecond).
 * CONSTANTS[91] is B62 in component IKr (per_millivolt).
 * CONSTANTS[92] is q62 in component IKr (dimensionless).
 * CONSTANTS[93] is A63 in component IKr (per_millisecond).
 * CONSTANTS[94] is B63 in component IKr (per_millivolt).
 * CONSTANTS[95] is q63 in component IKr (dimensionless).
 * CONSTANTS[96] is Kmax in component IKr (dimensionless).
 * CONSTANTS[97] is Ku in component IKr (per_millisecond).
 * CONSTANTS[98] is n in component IKr (dimensionless).
 * CONSTANTS[99] is halfmax in component IKr (dimensionless).
 * CONSTANTS[100] is Kt in component IKr (per_millisecond).
 * CONSTANTS[101] is Vhalf in component IKr (millivolt).
 * CONSTANTS[102] is Temp in component IKr (dimensionless).
 * CONSTANTS[103] is GKs_b in component IKs (milliS_per_microF).
 * CONSTANTS[160] is GKs in component IKs (milliS_per_microF).
 * ALGEBRAIC[10] is xs1ss in component IKs (dimensionless).
 * ALGEBRAIC[24] is xs2ss in component IKs (dimensionless).
 * ALGEBRAIC[25] is txs1 in component IKs (millisecond).
 * CONSTANTS[104] is txs1_max in component IKs (millisecond).
 * STATES[44] is xs1 in component IKs (dimensionless).
 * STATES[45] is xs2 in component IKs (dimensionless).
 * ALGEBRAIC[93] is KsCa in component IKs (dimensionless).
 * ALGEBRAIC[35] is txs2 in component IKs (millisecond).
 * CONSTANTS[161] is GK1 in component IK1 (milliS_per_microF).
 * CONSTANTS[105] is GK1_b in component IK1 (milliS_per_microF).
 * ALGEBRAIC[11] is xk1ss in component IK1 (dimensionless).
 * ALGEBRAIC[26] is txk1 in component IK1 (millisecond).
 * STATES[46] is xk1 in component IK1 (dimensionless).
 * ALGEBRAIC[97] is rk1 in component IK1 (millisecond).
 * CONSTANTS[106] is kna1 in component INaCa_i (per_millisecond).
 * CONSTANTS[107] is kna2 in component INaCa_i (per_millisecond).
 * CONSTANTS[108] is kna3 in component INaCa_i (per_millisecond).
 * CONSTANTS[109] is kasymm in component INaCa_i (dimensionless).
 * CONSTANTS[110] is wna in component INaCa_i (dimensionless).
 * CONSTANTS[111] is wca in component INaCa_i (dimensionless).
 * CONSTANTS[112] is wnaca in component INaCa_i (dimensionless).
 * CONSTANTS[113] is kcaon in component INaCa_i (per_millisecond).
 * CONSTANTS[114] is kcaoff in component INaCa_i (per_millisecond).
 * CONSTANTS[115] is qna in component INaCa_i (dimensionless).
 * CONSTANTS[116] is qca in component INaCa_i (dimensionless).
 * ALGEBRAIC[100] is hna in component INaCa_i (dimensionless).
 * ALGEBRAIC[99] is hca in component INaCa_i (dimensionless).
 * CONSTANTS[117] is KmCaAct in component INaCa_i (millimolar).
 * CONSTANTS[118] is Gncx_b in component INaCa_i (milliS_per_microF).
 * CONSTANTS[194] is Gncx in component INaCa_i (milliS_per_microF).
 * ALGEBRAIC[101] is h1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[102] is h2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[103] is h3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[104] is h4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[105] is h5_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[106] is h6_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[107] is h7_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[108] is h8_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[109] is h9_i in component INaCa_i (dimensionless).
 * CONSTANTS[188] is h10_i in component INaCa_i (dimensionless).
 * CONSTANTS[189] is h11_i in component INaCa_i (dimensionless).
 * CONSTANTS[190] is h12_i in component INaCa_i (dimensionless).
 * CONSTANTS[191] is k1_i in component INaCa_i (dimensionless).
 * CONSTANTS[192] is k2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[110] is k3p_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[111] is k3pp_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[112] is k3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[115] is k4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[113] is k4p_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[114] is k4pp_i in component INaCa_i (dimensionless).
 * CONSTANTS[193] is k5_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[116] is k6_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[117] is k7_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[118] is k8_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[119] is x1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[120] is x2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[121] is x3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[122] is x4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[123] is E1_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[124] is E2_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[125] is E3_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[126] is E4_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[127] is allo_i in component INaCa_i (dimensionless).
 * ALGEBRAIC[128] is JncxNa_i in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[129] is JncxCa_i in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[131] is h1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[132] is h2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[133] is h3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[134] is h4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[135] is h5_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[136] is h6_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[137] is h7_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[138] is h8_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[139] is h9_ss in component INaCa_i (dimensionless).
 * CONSTANTS[195] is h10_ss in component INaCa_i (dimensionless).
 * CONSTANTS[196] is h11_ss in component INaCa_i (dimensionless).
 * CONSTANTS[197] is h12_ss in component INaCa_i (dimensionless).
 * CONSTANTS[198] is k1_ss in component INaCa_i (dimensionless).
 * CONSTANTS[199] is k2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[140] is k3p_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[141] is k3pp_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[142] is k3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[145] is k4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[143] is k4p_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[144] is k4pp_ss in component INaCa_i (dimensionless).
 * CONSTANTS[200] is k5_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[146] is k6_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[147] is k7_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[148] is k8_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[149] is x1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[150] is x2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[151] is x3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[152] is x4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[153] is E1_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[154] is E2_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[155] is E3_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[156] is E4_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[157] is allo_ss in component INaCa_i (dimensionless).
 * ALGEBRAIC[158] is JncxNa_ss in component INaCa_i (millimolar_per_millisecond).
 * ALGEBRAIC[159] is JncxCa_ss in component INaCa_i (millimolar_per_millisecond).
 * CONSTANTS[119] is k1p in component INaK (per_millisecond).
 * CONSTANTS[120] is k1m in component INaK (per_millisecond).
 * CONSTANTS[121] is k2p in component INaK (per_millisecond).
 * CONSTANTS[122] is k2m in component INaK (per_millisecond).
 * CONSTANTS[123] is k3p in component INaK (per_millisecond).
 * CONSTANTS[124] is k3m in component INaK (per_millisecond).
 * CONSTANTS[125] is k4p in component INaK (per_millisecond).
 * CONSTANTS[126] is k4m in component INaK (per_millisecond).
 * CONSTANTS[127] is Knai0 in component INaK (millimolar).
 * CONSTANTS[128] is Knao0 in component INaK (millimolar).
 * CONSTANTS[129] is delta in component INaK (millivolt).
 * CONSTANTS[130] is Kki in component INaK (per_millisecond).
 * CONSTANTS[131] is Kko in component INaK (per_millisecond).
 * CONSTANTS[132] is MgADP in component INaK (millimolar).
 * CONSTANTS[133] is MgATP in component INaK (millimolar).
 * CONSTANTS[134] is Kmgatp in component INaK (millimolar).
 * CONSTANTS[135] is H in component INaK (millimolar).
 * CONSTANTS[136] is eP in component INaK (dimensionless).
 * CONSTANTS[137] is Khp in component INaK (millimolar).
 * CONSTANTS[138] is Knap in component INaK (millimolar).
 * CONSTANTS[139] is Kxkur in component INaK (millimolar).
 * CONSTANTS[140] is Pnak_b in component INaK (milliS_per_microF).
 * CONSTANTS[204] is Pnak in component INaK (milliS_per_microF).
 * ALGEBRAIC[161] is Knai in component INaK (millimolar).
 * ALGEBRAIC[162] is Knao in component INaK (millimolar).
 * ALGEBRAIC[163] is P in component INaK (dimensionless).
 * ALGEBRAIC[164] is a1 in component INaK (dimensionless).
 * CONSTANTS[201] is b1 in component INaK (dimensionless).
 * CONSTANTS[202] is a2 in component INaK (dimensionless).
 * ALGEBRAIC[165] is b2 in component INaK (dimensionless).
 * ALGEBRAIC[166] is a3 in component INaK (dimensionless).
 * ALGEBRAIC[167] is b3 in component INaK (dimensionless).
 * CONSTANTS[203] is a4 in component INaK (dimensionless).
 * ALGEBRAIC[168] is b4 in component INaK (dimensionless).
 * ALGEBRAIC[169] is x1 in component INaK (dimensionless).
 * ALGEBRAIC[170] is x2 in component INaK (dimensionless).
 * ALGEBRAIC[171] is x3 in component INaK (dimensionless).
 * ALGEBRAIC[172] is x4 in component INaK (dimensionless).
 * ALGEBRAIC[173] is E1 in component INaK (dimensionless).
 * ALGEBRAIC[174] is E2 in component INaK (dimensionless).
 * ALGEBRAIC[175] is E3 in component INaK (dimensionless).
 * ALGEBRAIC[176] is E4 in component INaK (dimensionless).
 * ALGEBRAIC[177] is JnakNa in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[178] is JnakK in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[180] is xkb in component IKb (dimensionless).
 * CONSTANTS[141] is GKb_b in component IKb (milliS_per_microF).
 * CONSTANTS[163] is GKb in component IKb (milliS_per_microF).
 * CONSTANTS[142] is PNab in component INab (milliS_per_microF).
 * ALGEBRAIC[182] is A in component INab (microA_per_microF).
 * CONSTANTS[178] is B in component INab (per_millivolt).
 * CONSTANTS[164] is v0 in component INab (millivolt).
 * ALGEBRAIC[184] is U in component INab (dimensionless).
 * CONSTANTS[143] is PCab in component ICab (milliS_per_microF).
 * ALGEBRAIC[186] is A in component ICab (microA_per_microF).
 * CONSTANTS[179] is B in component ICab (per_millivolt).
 * CONSTANTS[165] is v0 in component ICab (millivolt).
 * ALGEBRAIC[188] is U in component ICab (dimensionless).
 * CONSTANTS[144] is GpCa in component IpCa (milliS_per_microF).
 * CONSTANTS[145] is KmCap in component IpCa (millimolar).
 * CONSTANTS[146] is bt in component ryr (millisecond).
 * CONSTANTS[166] is a_rel in component ryr (millisecond).
 * ALGEBRAIC[88] is Jrel_inf in component ryr (dimensionless).
 * ALGEBRAIC[94] is tau_rel in component ryr (millisecond).
 * ALGEBRAIC[89] is Jrel_infp in component ryr (dimensionless).
 * ALGEBRAIC[86] is Jrel_temp in component ryr (dimensionless).
 * ALGEBRAIC[95] is tau_relp in component ryr (millisecond).
 * STATES[47] is Jrelnp in component ryr (dimensionless).
 * STATES[48] is Jrelp in component ryr (dimensionless).
 * CONSTANTS[167] is btp in component ryr (millisecond).
 * CONSTANTS[180] is a_relp in component ryr (millisecond).
 * ALGEBRAIC[85] is Jrel_inf_temp in component ryr (dimensionless).
 * ALGEBRAIC[192] is fJrelp in component ryr (dimensionless).
 * CONSTANTS[147] is Jrel_scaling_factor in component ryr (dimensionless).
 * ALGEBRAIC[91] is tau_rel_temp in component ryr (millisecond).
 * ALGEBRAIC[92] is tau_relp_temp in component ryr (millisecond).
 * CONSTANTS[168] is upScale in component SERCA (dimensionless).
 * ALGEBRAIC[194] is Jupnp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[195] is Jupp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[196] is fJupp in component SERCA (dimensionless).
 * ALGEBRAIC[197] is Jleak in component SERCA (millimolar_per_millisecond).
 * CONSTANTS[148] is Jup_b in component SERCA (dimensionless).
 * RATES[0] is d/dt v in component membrane (millivolt).
 * RATES[1] is d/dt CaMKt in component CaMK (millimolar).
 * RATES[3] is d/dt nai in component intracellular_ions (millimolar).
 * RATES[4] is d/dt nass in component intracellular_ions (millimolar).
 * RATES[5] is d/dt ki in component intracellular_ions (millimolar).
 * RATES[6] is d/dt kss in component intracellular_ions (millimolar).
 * RATES[9] is d/dt cai in component intracellular_ions (millimolar).
 * RATES[2] is d/dt cass in component intracellular_ions (millimolar).
 * RATES[7] is d/dt cansr in component intracellular_ions (millimolar).
 * RATES[8] is d/dt cajsr in component intracellular_ions (millimolar).
 * RATES[10] is d/dt m in component INa (dimensionless).
 * RATES[11] is d/dt hf in component INa (dimensionless).
 * RATES[12] is d/dt hs in component INa (dimensionless).
 * RATES[13] is d/dt j in component INa (dimensionless).
 * RATES[14] is d/dt hsp in component INa (dimensionless).
 * RATES[15] is d/dt jp in component INa (dimensionless).
 * RATES[16] is d/dt mL in component INaL (dimensionless).
 * RATES[17] is d/dt hL in component INaL (dimensionless).
 * RATES[18] is d/dt hLp in component INaL (dimensionless).
 * RATES[19] is d/dt a in component Ito (dimensionless).
 * RATES[20] is d/dt iF in component Ito (dimensionless).
 * RATES[21] is d/dt iS in component Ito (dimensionless).
 * RATES[22] is d/dt ap in component Ito (dimensionless).
 * RATES[23] is d/dt iFp in component Ito (dimensionless).
 * RATES[24] is d/dt iSp in component Ito (dimensionless).
 * RATES[25] is d/dt d in component ICaL (dimensionless).
 * RATES[26] is d/dt ff in component ICaL (dimensionless).
 * RATES[27] is d/dt fs in component ICaL (dimensionless).
 * RATES[28] is d/dt fcaf in component ICaL (dimensionless).
 * RATES[29] is d/dt fcas in component ICaL (dimensionless).
 * RATES[30] is d/dt jca in component ICaL (dimensionless).
 * RATES[31] is d/dt ffp in component ICaL (dimensionless).
 * RATES[32] is d/dt fcafp in component ICaL (dimensionless).
 * RATES[33] is d/dt nca in component ICaL (dimensionless).
 * RATES[34] is d/dt IC1 in component IKr (dimensionless).
 * RATES[35] is d/dt IC2 in component IKr (dimensionless).
 * RATES[36] is d/dt C1 in component IKr (dimensionless).
 * RATES[37] is d/dt C2 in component IKr (dimensionless).
 * RATES[38] is d/dt O in component IKr (dimensionless).
 * RATES[39] is d/dt IO in component IKr (dimensionless).
 * RATES[40] is d/dt IObound in component IKr (dimensionless).
 * RATES[41] is d/dt Obound in component IKr (dimensionless).
 * RATES[42] is d/dt Cbound in component IKr (dimensionless).
 * RATES[43] is d/dt D in component IKr (dimensionless).
 * RATES[44] is d/dt xs1 in component IKs (dimensionless).
 * RATES[45] is d/dt xs2 in component IKs (dimensionless).
 * RATES[46] is d/dt xk1 in component IK1 (dimensionless).
 * RATES[47] is d/dt Jrelnp in component ryr (dimensionless).
 * RATES[48] is d/dt Jrelp in component ryr (dimensionless).
 */
void
initConsts(double* CONSTANTS, double *STATES)
{
CONSTANTS[0] = 0;
CONSTANTS[1] = 140;
CONSTANTS[2] = 1.8;
CONSTANTS[3] = 5.4;
CONSTANTS[4] = 8314;
CONSTANTS[5] = 310;
CONSTANTS[6] = 96485;
CONSTANTS[7] = 1;
CONSTANTS[8] = 2;
CONSTANTS[9] = 1;
CONSTANTS[10] = 0.01;
CONSTANTS[11] = 0.0011;
STATES[0] = -88.00190465;
CONSTANTS[12] = 0; // 10; i_Stim_Start
CONSTANTS[13] = 100000000000000000;
CONSTANTS[14] = -80;
CONSTANTS[15] = 1000; // i_Stim_Period
CONSTANTS[16] = 0.5; // pulse duration
CONSTANTS[17] = 0.15;
CONSTANTS[18] = 0.05;
CONSTANTS[19] = 0.00068;
CONSTANTS[20] = 0.05;
CONSTANTS[21] = 0.0015;
STATES[1] = 0.0125840447;
STATES[2] = 8.49e-05;
CONSTANTS[22] = 0.05;
CONSTANTS[23] = 0.00238;
CONSTANTS[24] = 0.07;
CONSTANTS[25] = 0.0005;
CONSTANTS[26] = 0.047;
CONSTANTS[27] = 0.00087;
CONSTANTS[28] = 1.124;
CONSTANTS[29] = 0.0087;
CONSTANTS[30] = 10;
CONSTANTS[31] = 0.8;
STATES[3] = 7.268004498;
STATES[4] = 7.268089977;
STATES[5] = 144.6555918;
STATES[6] = 144.6555651;
STATES[7] = 1.619574538;
STATES[8] = 1.571234014;
STATES[9] = 8.6e-05;
CONSTANTS[32] = 1;
CONSTANTS[33] = 0.01833;
CONSTANTS[34] = 39.57;
CONSTANTS[35] = 9.871;
CONSTANTS[36] = 11.64;
CONSTANTS[37] = 34.77;
CONSTANTS[38] = 6.765;
CONSTANTS[39] = 8.552;
CONSTANTS[40] = 77.42;
CONSTANTS[41] = 5.955;
STATES[10] = 0.007344121102;
CONSTANTS[42] = 82.9;
CONSTANTS[43] = 6.086;
CONSTANTS[44] = 0.99;
STATES[11] = 0.6981071913;
STATES[12] = 0.6980895801;
CONSTANTS[45] = 75;
CONSTANTS[46] = 0;
STATES[13] = 0.6979908432;
STATES[14] = 0.4549485525;
STATES[15] = 0.6979245865;
STATES[16] = 0.0001882617273;
CONSTANTS[47] = 200;
STATES[17] = 0.5008548855;
STATES[18] = 0.2693065357;
CONSTANTS[48] = 0.019957499999999975;
CONSTANTS[49] = 0.02;
STATES[19] = 0.001001097687;
STATES[20] = 0.9995541745;
STATES[21] = 0.5865061736;
STATES[22] = 0.0005100862934;
STATES[23] = 0.9995541823;
STATES[24] = 0.6393399482;
CONSTANTS[50] = 0.002;
CONSTANTS[51] = 1000;
CONSTANTS[52] = 0.0001007;
STATES[25] = 2.34e-9;
STATES[26] = 0.9999999909;
STATES[27] = 0.9102412777;
STATES[28] = 0.9999999909;
STATES[29] = 0.9998046777;
STATES[30] = 0.9999738312;
STATES[31] = 0.9999999909;
STATES[32] = 0.9999999909;
STATES[33] = 0.002749414044;
CONSTANTS[53] = 0.04658545454545456;
STATES[34] = 0.999637;
STATES[35] = 6.83208e-05;
STATES[36] = 1.80145e-08;
STATES[37] = 8.26619e-05;
STATES[38] = 0.00015551;
STATES[39] = 5.67623e-05;
STATES[40] = 0;
STATES[41] = 0;
STATES[42] = 0;
STATES[43] = 0;
CONSTANTS[54] = 0.0264;
CONSTANTS[55] = 4.631E-05;
CONSTANTS[56] = 4.843;
CONSTANTS[57] = 4.986E-06;
CONSTANTS[58] = -0.004226;
CONSTANTS[59] = 4.23;
CONSTANTS[60] = 0.001214;
CONSTANTS[61] = 0.008516;
CONSTANTS[62] = 4.962;
CONSTANTS[63] = 1.854E-05;
CONSTANTS[64] = -0.04641;
CONSTANTS[65] = 3.769;
CONSTANTS[66] = 0.0007868;
CONSTANTS[67] = 1.535E-08;
CONSTANTS[68] = 4.942;
CONSTANTS[69] = 5.455E-06;
CONSTANTS[70] = -0.1688;
CONSTANTS[71] = 4.156;
CONSTANTS[72] = 0.005509;
CONSTANTS[73] = 7.771E-09;
CONSTANTS[74] = 4.22;
CONSTANTS[75] = 0.001416;
CONSTANTS[76] = -0.02877;
CONSTANTS[77] = 1.459;
CONSTANTS[78] = 0.4492;
CONSTANTS[79] = 0.008595;
CONSTANTS[80] = 5;
CONSTANTS[81] = 0.3181;
CONSTANTS[82] = 3.613E-08;
CONSTANTS[83] = 4.663;
CONSTANTS[84] = 0.149;
CONSTANTS[85] = 0.004668;
CONSTANTS[86] = 2.412;
CONSTANTS[87] = 0.01241;
CONSTANTS[88] = 0.1725;
CONSTANTS[89] = 5.568;
CONSTANTS[90] = 0.3226;
CONSTANTS[91] = -0.0006575;
CONSTANTS[92] = 5;
CONSTANTS[93] = 0.008978;
CONSTANTS[94] = -0.02215;
CONSTANTS[95] = 5.682;
CONSTANTS[96] = 0;
CONSTANTS[97] = 0;
CONSTANTS[98] = 1;
CONSTANTS[99] = 1;
CONSTANTS[100] = 0;
CONSTANTS[101] = 1;
CONSTANTS[102] = 37;
CONSTANTS[103] = 0.006358000000000001;
CONSTANTS[104] = 817.3;
STATES[44] = 0.2707758025;
STATES[45] = 0.0001928503426;
CONSTANTS[105] = 0.3239783999999998;
STATES[46] = 0.9967597594;
CONSTANTS[106] = 15;
CONSTANTS[107] = 5;
CONSTANTS[108] = 88.12;
CONSTANTS[109] = 12.5;
CONSTANTS[110] = 6e4;
CONSTANTS[111] = 6e4;
CONSTANTS[112] = 5e3;
CONSTANTS[113] = 1.5e6;
CONSTANTS[114] = 5e3;
CONSTANTS[115] = 0.5224;
CONSTANTS[116] = 0.167;
CONSTANTS[117] = 150e-6;
CONSTANTS[118] = 0.0008;
CONSTANTS[119] = 949.5;
CONSTANTS[120] = 182.4;
CONSTANTS[121] = 687.2;
CONSTANTS[122] = 39.4;
CONSTANTS[123] = 1899;
CONSTANTS[124] = 79300;
CONSTANTS[125] = 639;
CONSTANTS[126] = 40;
CONSTANTS[127] = 9.073;
CONSTANTS[128] = 27.78;
CONSTANTS[129] = -0.155;
CONSTANTS[130] = 0.5;
CONSTANTS[131] = 0.3582;
CONSTANTS[132] = 0.05;
CONSTANTS[133] = 9.8;
CONSTANTS[134] = 1.698e-7;
CONSTANTS[135] = 1e-7;
CONSTANTS[136] = 4.2;
CONSTANTS[137] = 1.698e-7;
CONSTANTS[138] = 224;
CONSTANTS[139] = 292;
CONSTANTS[140] = 30;
CONSTANTS[141] = 0.003;
CONSTANTS[142] = 3.75e-10;
CONSTANTS[143] = 2.5e-8;
CONSTANTS[144] = 0.0005;
CONSTANTS[145] = 0.0005;
CONSTANTS[146] = 4.75;
STATES[47] = 2.5e-7;
STATES[48] = 3.12e-7;
CONSTANTS[147] = 1.0;
CONSTANTS[148] = 1.0;
CONSTANTS[149] = CONSTANTS[6]/( CONSTANTS[4]*CONSTANTS[5]);
CONSTANTS[150] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[22]*1.30000 : CONSTANTS[22]);
CONSTANTS[151] = 1.00000 - CONSTANTS[44];
CONSTANTS[152] =  3.00000*CONSTANTS[47];
CONSTANTS[153] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[48]*0.600000 : CONSTANTS[48]);
CONSTANTS[154] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[49]*4.00000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[49]*4.00000 : CONSTANTS[49]);
CONSTANTS[155] = 0.600000;
CONSTANTS[156] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[52]*1.20000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[52]*2.50000 : CONSTANTS[52]);
CONSTANTS[157] = 75.0000;
CONSTANTS[158] = 0.00000;
CONSTANTS[159] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[53]*1.30000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[53]*0.800000 : CONSTANTS[53]);
CONSTANTS[160] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[103]*1.40000 : CONSTANTS[103]);
CONSTANTS[161] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[105]*1.20000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[105]*1.30000 : CONSTANTS[105]);
CONSTANTS[162] =  1000.00*3.14000*CONSTANTS[11]*CONSTANTS[11]*CONSTANTS[10];
CONSTANTS[163] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[141]*0.600000 : CONSTANTS[141]);
CONSTANTS[164] = 0.00000;
CONSTANTS[165] = 0.00000;
CONSTANTS[166] =  0.500000*CONSTANTS[146];
CONSTANTS[167] =  1.25000*CONSTANTS[146];
CONSTANTS[168] = (CONSTANTS[0]==1.00000 ? 1.30000 : 1.00000);
CONSTANTS[205] = 0.00000;
CONSTANTS[169] =  CONSTANTS[6]*CONSTANTS[149];
CONSTANTS[170] = 1.00000 - CONSTANTS[155];
CONSTANTS[171] =  1.10000*CONSTANTS[156];
CONSTANTS[172] =  0.00125000*CONSTANTS[156];
CONSTANTS[173] =  0.000357400*CONSTANTS[156];
CONSTANTS[174] =  2.00000*CONSTANTS[149];
CONSTANTS[175] = CONSTANTS[149];
CONSTANTS[176] = CONSTANTS[149];
CONSTANTS[177] =  2.00000*3.14000*CONSTANTS[11]*CONSTANTS[11]+ 2.00000*3.14000*CONSTANTS[11]*CONSTANTS[10];
CONSTANTS[178] = CONSTANTS[149];
CONSTANTS[179] =  2.00000*CONSTANTS[149];
CONSTANTS[180] =  0.500000*CONSTANTS[167];
CONSTANTS[181] =  0.00125000*CONSTANTS[171];
CONSTANTS[182] =  0.000357400*CONSTANTS[171];
CONSTANTS[183] =  2.00000*CONSTANTS[177];
CONSTANTS[184] =  0.680000*CONSTANTS[162];
CONSTANTS[185] =  0.0552000*CONSTANTS[162];
CONSTANTS[186] =  0.00480000*CONSTANTS[162];
CONSTANTS[187] =  0.0200000*CONSTANTS[162];
CONSTANTS[188] = CONSTANTS[109]+1.00000+ (CONSTANTS[1]/CONSTANTS[106])*(1.00000+CONSTANTS[1]/CONSTANTS[107]);
CONSTANTS[189] = ( CONSTANTS[1]*CONSTANTS[1])/( CONSTANTS[188]*CONSTANTS[106]*CONSTANTS[107]);
CONSTANTS[190] = 1.00000/CONSTANTS[188];
CONSTANTS[191] =  CONSTANTS[190]*CONSTANTS[2]*CONSTANTS[113];
CONSTANTS[192] = CONSTANTS[114];
CONSTANTS[193] = CONSTANTS[114];
CONSTANTS[194] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[118]*1.10000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[118]*1.40000 : CONSTANTS[118]);
CONSTANTS[195] = CONSTANTS[109]+1.00000+ (CONSTANTS[1]/CONSTANTS[106])*(1.00000+CONSTANTS[1]/CONSTANTS[107]);
CONSTANTS[196] = ( CONSTANTS[1]*CONSTANTS[1])/( CONSTANTS[195]*CONSTANTS[106]*CONSTANTS[107]);
CONSTANTS[197] = 1.00000/CONSTANTS[195];
CONSTANTS[198] =  CONSTANTS[197]*CONSTANTS[2]*CONSTANTS[113];
CONSTANTS[199] = CONSTANTS[114];
CONSTANTS[200] = CONSTANTS[114];
CONSTANTS[201] =  CONSTANTS[120]*CONSTANTS[132];
CONSTANTS[202] = CONSTANTS[121];
CONSTANTS[203] = (( CONSTANTS[125]*CONSTANTS[133])/CONSTANTS[134])/(1.00000+CONSTANTS[133]/CONSTANTS[134]);
CONSTANTS[204] = (CONSTANTS[0]==1.00000 ?  CONSTANTS[140]*0.900000 : CONSTANTS[0]==2.00000 ?  CONSTANTS[140]*0.700000 : CONSTANTS[140]);
}
void
computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
RATES[43] = CONSTANTS[205];
RATES[34] = (- ( CONSTANTS[66]*exp( CONSTANTS[67]*STATES[0])*STATES[34]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[68]))/10.0000) -  CONSTANTS[69]*exp( CONSTANTS[70]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[71]))/10.0000))+ CONSTANTS[78]*exp( CONSTANTS[79]*STATES[0])*STATES[36]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[80]))/10.0000)) -  CONSTANTS[87]*exp( CONSTANTS[88]*STATES[0])*STATES[34]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[89]))/10.0000);
RATES[35] = ((( CONSTANTS[66]*exp( CONSTANTS[67]*STATES[0])*STATES[34]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[68]))/10.0000) -  CONSTANTS[69]*exp( CONSTANTS[70]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[71]))/10.0000)) - ( CONSTANTS[60]*exp( CONSTANTS[61]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[62]))/10.0000) -  CONSTANTS[63]*exp( CONSTANTS[64]*STATES[0])*STATES[39]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[65]))/10.0000)))+ CONSTANTS[81]*exp( CONSTANTS[82]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[83]))/10.0000)) -  CONSTANTS[90]*exp( CONSTANTS[91]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[92]))/10.0000);
RATES[36] = - ( CONSTANTS[54]*exp( CONSTANTS[55]*STATES[0])*STATES[36]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[56]))/10.0000) -  CONSTANTS[57]*exp( CONSTANTS[58]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[59]))/10.0000)) - ( CONSTANTS[78]*exp( CONSTANTS[79]*STATES[0])*STATES[36]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[80]))/10.0000) -  CONSTANTS[87]*exp( CONSTANTS[88]*STATES[0])*STATES[34]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[89]))/10.0000));
RATES[37] = (( CONSTANTS[54]*exp( CONSTANTS[55]*STATES[0])*STATES[36]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[56]))/10.0000) -  CONSTANTS[57]*exp( CONSTANTS[58]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[59]))/10.0000)) - ( CONSTANTS[72]*exp( CONSTANTS[73]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[74]))/10.0000) -  CONSTANTS[75]*exp( CONSTANTS[76]*STATES[0])*STATES[38]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[77]))/10.0000))) - ( CONSTANTS[81]*exp( CONSTANTS[82]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[83]))/10.0000) -  CONSTANTS[90]*exp( CONSTANTS[91]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[92]))/10.0000));
RATES[38] = (( CONSTANTS[72]*exp( CONSTANTS[73]*STATES[0])*STATES[37]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[74]))/10.0000) -  CONSTANTS[75]*exp( CONSTANTS[76]*STATES[0])*STATES[38]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[77]))/10.0000)) - ( CONSTANTS[84]*exp( CONSTANTS[85]*STATES[0])*STATES[38]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[86]))/10.0000) -  CONSTANTS[93]*exp( CONSTANTS[94]*STATES[0])*STATES[39]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[95]))/10.0000))) - ( (( CONSTANTS[96]*CONSTANTS[97]*exp( CONSTANTS[98]*log(STATES[43])))/(exp( CONSTANTS[98]*log(STATES[43]))+CONSTANTS[99]))*STATES[38] -  CONSTANTS[97]*STATES[41]);
RATES[39] = ((( CONSTANTS[60]*exp( CONSTANTS[61]*STATES[0])*STATES[35]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[62]))/10.0000) -  CONSTANTS[63]*exp( CONSTANTS[64]*STATES[0])*STATES[39]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[65]))/10.0000))+ CONSTANTS[84]*exp( CONSTANTS[85]*STATES[0])*STATES[38]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[86]))/10.0000)) -  CONSTANTS[93]*exp( CONSTANTS[94]*STATES[0])*STATES[39]*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[95]))/10.0000)) - ( (( CONSTANTS[96]*CONSTANTS[97]*exp( CONSTANTS[98]*log(STATES[43])))/(exp( CONSTANTS[98]*log(STATES[43]))+CONSTANTS[99]))*STATES[39] -  (( CONSTANTS[97]*CONSTANTS[84]*exp( CONSTANTS[85]*STATES[0])*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[86]))/10.0000))/( CONSTANTS[93]*exp( CONSTANTS[94]*STATES[0])*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[95]))/10.0000)))*STATES[40]);
RATES[40] = (( (( CONSTANTS[96]*CONSTANTS[97]*exp( CONSTANTS[98]*log(STATES[43])))/(exp( CONSTANTS[98]*log(STATES[43]))+CONSTANTS[99]))*STATES[39] -  (( CONSTANTS[97]*CONSTANTS[84]*exp( CONSTANTS[85]*STATES[0])*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[86]))/10.0000))/( CONSTANTS[93]*exp( CONSTANTS[94]*STATES[0])*exp(( (CONSTANTS[102] - 20.0000)*log(CONSTANTS[95]))/10.0000)))*STATES[40])+ (CONSTANTS[100]/(1.00000+exp(- (STATES[0] - CONSTANTS[101])/6.78900)))*STATES[42]) -  CONSTANTS[100]*STATES[40];
RATES[41] = (( (( CONSTANTS[96]*CONSTANTS[97]*exp( CONSTANTS[98]*log(STATES[43])))/(exp( CONSTANTS[98]*log(STATES[43]))+CONSTANTS[99]))*STATES[38] -  CONSTANTS[97]*STATES[41])+ (CONSTANTS[100]/(1.00000+exp(- (STATES[0] - CONSTANTS[101])/6.78900)))*STATES[42]) -  CONSTANTS[100]*STATES[41];
RATES[42] = - ( (CONSTANTS[100]/(1.00000+exp(- (STATES[0] - CONSTANTS[101])/6.78900)))*STATES[42] -  CONSTANTS[100]*STATES[41]) - ( (CONSTANTS[100]/(1.00000+exp(- (STATES[0] - CONSTANTS[101])/6.78900)))*STATES[42] -  CONSTANTS[100]*STATES[40]);
ALGEBRAIC[3] = 1.00000/(1.00000+exp((STATES[0]+87.6100)/7.48800));
RATES[17] = (ALGEBRAIC[3] - STATES[17])/CONSTANTS[47];
ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[0]+93.8100)/7.48800));
RATES[18] = (ALGEBRAIC[4] - STATES[18])/CONSTANTS[152];
ALGEBRAIC[1] = 1.00000/(1.00000+exp(- (STATES[0]+CONSTANTS[34])/CONSTANTS[35]));
ALGEBRAIC[13] = 1.00000/( CONSTANTS[38]*exp((STATES[0]+CONSTANTS[36])/CONSTANTS[37])+ CONSTANTS[39]*exp(- (STATES[0]+CONSTANTS[40])/CONSTANTS[41]));
RATES[10] = (ALGEBRAIC[1] - STATES[10])/ALGEBRAIC[13];
ALGEBRAIC[2] = 1.00000/(1.00000+exp(((STATES[0]+CONSTANTS[42]) - CONSTANTS[46])/CONSTANTS[43]));
ALGEBRAIC[14] = 1.00000/( 1.43200e-05*exp(- ((STATES[0]+1.19600) - CONSTANTS[46])/6.28500)+ 6.14900*exp(((STATES[0]+0.509600) - CONSTANTS[46])/20.2700));
RATES[11] = (ALGEBRAIC[2] - STATES[11])/ALGEBRAIC[14];
ALGEBRAIC[15] = 1.00000/( 0.00979400*exp(- ((STATES[0]+17.9500) - CONSTANTS[46])/28.0500)+ 0.334300*exp(((STATES[0]+5.73000) - CONSTANTS[46])/56.6600));
RATES[12] = (ALGEBRAIC[2] - STATES[12])/ALGEBRAIC[15];
ALGEBRAIC[5] = 1.00000/(1.00000+exp(- (STATES[0] - 14.3400)/14.8200));
ALGEBRAIC[17] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- (STATES[0] - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[0]+100.000)/29.3814)));
RATES[19] = (ALGEBRAIC[5] - STATES[19])/ALGEBRAIC[17];
ALGEBRAIC[7] = 1.00000/(1.00000+exp(- (STATES[0]+3.94000)/4.23000));
ALGEBRAIC[21] = 0.600000+1.00000/(exp( - 0.0500000*(STATES[0]+6.00000))+exp( 0.0900000*(STATES[0]+14.0000)));
RATES[25] = (ALGEBRAIC[7] - STATES[25])/ALGEBRAIC[21];
ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+19.5800)/3.69600));
ALGEBRAIC[22] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[0]+20.0000)/10.0000)+ 0.00450000*exp((STATES[0]+20.0000)/10.0000));
RATES[26] = (ALGEBRAIC[8] - STATES[26])/ALGEBRAIC[22];
ALGEBRAIC[23] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[0]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[0]+5.00000)/6.00000));
RATES[27] = (ALGEBRAIC[8] - STATES[27])/ALGEBRAIC[23];
ALGEBRAIC[19] = ALGEBRAIC[8];
RATES[30] = (ALGEBRAIC[19] - STATES[30])/CONSTANTS[157];
ALGEBRAIC[9] =  STATES[30]*1.00000;
ALGEBRAIC[20] = 1.00000/(CONSTANTS[51]/ALGEBRAIC[9]+pow(1.00000+CONSTANTS[50]/STATES[2], 4.00000));
RATES[33] =  ALGEBRAIC[20]*CONSTANTS[51] -  STATES[33]*ALGEBRAIC[9];
ALGEBRAIC[10] = 1.00000/(1.00000+exp(- (STATES[0]+11.6000)/8.93200));
ALGEBRAIC[25] = CONSTANTS[104]+1.00000/( 0.000232600*exp((STATES[0]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[0]+210.000)/230.000));
RATES[44] = (ALGEBRAIC[10] - STATES[44])/ALGEBRAIC[25];
ALGEBRAIC[11] = 1.00000/(1.00000+exp(- (STATES[0]+ 2.55380*CONSTANTS[3]+144.590)/( 1.56920*CONSTANTS[3]+3.81150)));
ALGEBRAIC[26] = 122.200/(exp(- (STATES[0]+127.200)/20.3600)+exp((STATES[0]+236.800)/69.3300));
RATES[46] = (ALGEBRAIC[11] - STATES[46])/ALGEBRAIC[26];
ALGEBRAIC[36] = ( CONSTANTS[20]*(1.00000 - STATES[1]))/(1.00000+CONSTANTS[21]/STATES[2]);
RATES[1] =  CONSTANTS[18]*ALGEBRAIC[36]*(ALGEBRAIC[36]+STATES[1]) -  CONSTANTS[19]*STATES[1];
ALGEBRAIC[16] = ALGEBRAIC[2];
ALGEBRAIC[27] = 2.03800+1.00000/( 0.0213600*exp(- ((STATES[0]+100.600) - CONSTANTS[46])/8.28100)+ 0.305200*exp(((STATES[0]+0.994100) - CONSTANTS[46])/38.4500));
RATES[13] = (ALGEBRAIC[16] - STATES[13])/ALGEBRAIC[27];
ALGEBRAIC[31] = 1.00000/(1.00000+exp(- (STATES[0] - 24.3400)/14.8200));
RATES[22] = (ALGEBRAIC[31] - STATES[22])/ALGEBRAIC[17];
ALGEBRAIC[32] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[0] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[0] - 4.00000)/7.00000));
RATES[28] = (ALGEBRAIC[19] - STATES[28])/ALGEBRAIC[32];
ALGEBRAIC[33] = 100.000+1.00000/( 0.000120000*exp(- STATES[0]/3.00000)+ 0.000120000*exp(STATES[0]/7.00000));
RATES[29] = (ALGEBRAIC[19] - STATES[29])/ALGEBRAIC[33];
ALGEBRAIC[34] =  2.50000*ALGEBRAIC[22];
RATES[31] = (ALGEBRAIC[8] - STATES[31])/ALGEBRAIC[34];
ALGEBRAIC[24] = ALGEBRAIC[10];
ALGEBRAIC[35] = 1.00000/( 0.0100000*exp((STATES[0] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[0]+66.5400)/31.0000));
RATES[45] = (ALGEBRAIC[24] - STATES[45])/ALGEBRAIC[35];
ALGEBRAIC[28] = 1.00000/(1.00000+exp(((STATES[0]+89.1000) - CONSTANTS[46])/6.08600));
ALGEBRAIC[37] =  3.00000*ALGEBRAIC[15];
RATES[14] = (ALGEBRAIC[28] - STATES[14])/ALGEBRAIC[37];
ALGEBRAIC[38] =  1.46000*ALGEBRAIC[27];
RATES[15] = (ALGEBRAIC[16] - STATES[15])/ALGEBRAIC[38];
ALGEBRAIC[29] = 1.00000/(1.00000+exp(- (STATES[0]+42.8500)/5.26400));
ALGEBRAIC[39] = ALGEBRAIC[13];
RATES[16] = (ALGEBRAIC[29] - STATES[16])/ALGEBRAIC[39];
ALGEBRAIC[41] =  2.50000*ALGEBRAIC[32];
RATES[32] = (ALGEBRAIC[19] - STATES[32])/ALGEBRAIC[41];
ALGEBRAIC[6] = 1.00000/(1.00000+exp((STATES[0]+43.9400)/5.71100));
ALGEBRAIC[18] = (CONSTANTS[0]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[0]+70.0000)/5.00000)) : 1.00000);
ALGEBRAIC[30] = 4.56200+1.00000/( 0.393300*exp(- (STATES[0]+100.000)/100.000)+ 0.0800400*exp((STATES[0]+50.0000)/16.5900));
ALGEBRAIC[43] =  ALGEBRAIC[30]*ALGEBRAIC[18];
RATES[20] = (ALGEBRAIC[6] - STATES[20])/ALGEBRAIC[43];
ALGEBRAIC[40] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[0]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[0]+114.100)/8.07900));
ALGEBRAIC[45] =  ALGEBRAIC[40]*ALGEBRAIC[18];
RATES[21] = (ALGEBRAIC[6] - STATES[21])/ALGEBRAIC[45];
ALGEBRAIC[47] = 1.35400+0.000100000/(exp((STATES[0] - 167.400)/15.8900)+exp(- (STATES[0] - 12.2300)/0.215400));
ALGEBRAIC[49] = 1.00000 - 0.500000/(1.00000+exp((STATES[0]+70.0000)/20.0000));
ALGEBRAIC[51] =  ALGEBRAIC[47]*ALGEBRAIC[49]*ALGEBRAIC[43];
RATES[23] = (ALGEBRAIC[6] - STATES[23])/ALGEBRAIC[51];
ALGEBRAIC[52] =  ALGEBRAIC[47]*ALGEBRAIC[49]*ALGEBRAIC[45];
RATES[24] = (ALGEBRAIC[6] - STATES[24])/ALGEBRAIC[52];
ALGEBRAIC[67] =  CONSTANTS[155]*STATES[26]+ CONSTANTS[170]*STATES[27];
ALGEBRAIC[68] = 0.300000+0.600000/(1.00000+exp((STATES[0] - 10.0000)/10.0000));
ALGEBRAIC[69] = 1.00000 - ALGEBRAIC[68];
ALGEBRAIC[70] =  ALGEBRAIC[68]*STATES[28]+ ALGEBRAIC[69]*STATES[29];
ALGEBRAIC[71] =  CONSTANTS[155]*STATES[31]+ CONSTANTS[170]*STATES[27];
ALGEBRAIC[72] =  ALGEBRAIC[68]*STATES[32]+ ALGEBRAIC[69]*STATES[29];
ALGEBRAIC[12] =  STATES[0]*CONSTANTS[149];
ALGEBRAIC[73] = ( 4.00000*CONSTANTS[169]*( STATES[2]*exp( 2.00000*ALGEBRAIC[12]) -  0.341000*CONSTANTS[2]))/CONSTANTS[174];
ALGEBRAIC[74] =  CONSTANTS[174]*(STATES[0] - CONSTANTS[158]);
ALGEBRAIC[75] = (- 1.00000e-07<=ALGEBRAIC[74]&&ALGEBRAIC[74]<=1.00000e-07 ?  ALGEBRAIC[73]*(1.00000 -  0.500000*ALGEBRAIC[74]) : ( ALGEBRAIC[73]*ALGEBRAIC[74])/(exp(ALGEBRAIC[74]) - 1.00000));
ALGEBRAIC[42] = ALGEBRAIC[36]+STATES[1];
ALGEBRAIC[82] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[83] =  (1.00000 - ALGEBRAIC[82])*CONSTANTS[156]*ALGEBRAIC[75]*STATES[25]*( ALGEBRAIC[67]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[70]*STATES[33])+ ALGEBRAIC[82]*CONSTANTS[171]*ALGEBRAIC[75]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[72]*STATES[33]);
ALGEBRAIC[85] = ( CONSTANTS[166]*- ALGEBRAIC[83])/(1.00000+ 1.00000*pow(1.50000/STATES[8], 8.00000));
ALGEBRAIC[88] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[85]*1.70000 : ALGEBRAIC[85]);
ALGEBRAIC[91] = CONSTANTS[146]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[94] = (ALGEBRAIC[91]<0.00100000 ? 0.00100000 : ALGEBRAIC[91]);
RATES[47] = (ALGEBRAIC[88] - STATES[47])/ALGEBRAIC[94];
ALGEBRAIC[86] = ( CONSTANTS[180]*- ALGEBRAIC[83])/(1.00000+pow(1.50000/STATES[8], 8.00000));
ALGEBRAIC[89] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[86]*1.70000 : ALGEBRAIC[86]);
ALGEBRAIC[92] = CONSTANTS[167]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[95] = (ALGEBRAIC[92]<0.00100000 ? 0.00100000 : ALGEBRAIC[92]);
RATES[48] = (ALGEBRAIC[89] - STATES[48])/ALGEBRAIC[95];
ALGEBRAIC[53] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log(CONSTANTS[3]/STATES[5]);
ALGEBRAIC[61] = 1.00000/(1.00000+exp((STATES[0] - 213.600)/151.200));
ALGEBRAIC[62] = 1.00000 - ALGEBRAIC[61];
ALGEBRAIC[63] =  ALGEBRAIC[61]*STATES[20]+ ALGEBRAIC[62]*STATES[21];
ALGEBRAIC[64] =  ALGEBRAIC[61]*STATES[23]+ ALGEBRAIC[62]*STATES[24];
ALGEBRAIC[65] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[66] =  CONSTANTS[154]*(STATES[0] - ALGEBRAIC[53])*( (1.00000 - ALGEBRAIC[65])*STATES[19]*ALGEBRAIC[63]+ ALGEBRAIC[65]*STATES[22]*ALGEBRAIC[64]);
ALGEBRAIC[90] =  CONSTANTS[159]* pow((CONSTANTS[3]/5.40000), 1.0 / 2)*STATES[38]*(STATES[0] - ALGEBRAIC[53]);
ALGEBRAIC[54] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log((CONSTANTS[3]+ CONSTANTS[33]*CONSTANTS[1])/(STATES[5]+ CONSTANTS[33]*STATES[3]));
ALGEBRAIC[93] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[9], 1.40000));
ALGEBRAIC[96] =  CONSTANTS[160]*ALGEBRAIC[93]*STATES[44]*STATES[45]*(STATES[0] - ALGEBRAIC[54]);
ALGEBRAIC[97] = 1.00000/(1.00000+exp(((STATES[0]+105.800) -  2.60000*CONSTANTS[3])/9.49300));
ALGEBRAIC[98] =  CONSTANTS[161]* pow(CONSTANTS[3], 1.0 / 2)*ALGEBRAIC[97]*STATES[46]*(STATES[0] - ALGEBRAIC[53]);
ALGEBRAIC[162] =  CONSTANTS[128]*exp(( (1.00000 - CONSTANTS[129])*STATES[0]*CONSTANTS[6])/( 3.00000*CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[166] = ( CONSTANTS[123]*pow(CONSTANTS[3]/CONSTANTS[131], 2.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[162], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[131], 2.00000)) - 1.00000);
ALGEBRAIC[163] = CONSTANTS[136]/(1.00000+CONSTANTS[135]/CONSTANTS[137]+STATES[3]/CONSTANTS[138]+STATES[5]/CONSTANTS[139]);
ALGEBRAIC[167] = ( CONSTANTS[124]*ALGEBRAIC[163]*CONSTANTS[135])/(1.00000+CONSTANTS[133]/CONSTANTS[134]);
ALGEBRAIC[161] =  CONSTANTS[127]*exp(( CONSTANTS[129]*STATES[0]*CONSTANTS[6])/( 3.00000*CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[164] = ( CONSTANTS[119]*pow(STATES[3]/ALGEBRAIC[161], 3.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[161], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[130], 2.00000)) - 1.00000);
ALGEBRAIC[165] = ( CONSTANTS[122]*pow(CONSTANTS[1]/ALGEBRAIC[162], 3.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[162], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[131], 2.00000)) - 1.00000);
ALGEBRAIC[168] = ( CONSTANTS[126]*pow(STATES[5]/CONSTANTS[130], 2.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[161], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[130], 2.00000)) - 1.00000);
ALGEBRAIC[169] =  CONSTANTS[203]*ALGEBRAIC[164]*CONSTANTS[202]+ ALGEBRAIC[165]*ALGEBRAIC[168]*ALGEBRAIC[167]+ CONSTANTS[202]*ALGEBRAIC[168]*ALGEBRAIC[167]+ ALGEBRAIC[167]*ALGEBRAIC[164]*CONSTANTS[202];
ALGEBRAIC[170] =  ALGEBRAIC[165]*CONSTANTS[201]*ALGEBRAIC[168]+ ALGEBRAIC[164]*CONSTANTS[202]*ALGEBRAIC[166]+ ALGEBRAIC[166]*CONSTANTS[201]*ALGEBRAIC[168]+ CONSTANTS[202]*ALGEBRAIC[166]*ALGEBRAIC[168];
ALGEBRAIC[171] =  CONSTANTS[202]*ALGEBRAIC[166]*CONSTANTS[203]+ ALGEBRAIC[167]*ALGEBRAIC[165]*CONSTANTS[201]+ ALGEBRAIC[165]*CONSTANTS[201]*CONSTANTS[203]+ ALGEBRAIC[166]*CONSTANTS[203]*CONSTANTS[201];
ALGEBRAIC[172] =  ALGEBRAIC[168]*ALGEBRAIC[167]*ALGEBRAIC[165]+ ALGEBRAIC[166]*CONSTANTS[203]*ALGEBRAIC[164]+ ALGEBRAIC[165]*CONSTANTS[203]*ALGEBRAIC[164]+ ALGEBRAIC[167]*ALGEBRAIC[165]*ALGEBRAIC[164];
ALGEBRAIC[173] = ALGEBRAIC[169]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[174] = ALGEBRAIC[170]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[177] =  3.00000*( ALGEBRAIC[173]*ALGEBRAIC[166] -  ALGEBRAIC[174]*ALGEBRAIC[167]);
ALGEBRAIC[175] = ALGEBRAIC[171]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[176] = ALGEBRAIC[172]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[178] =  2.00000*( ALGEBRAIC[176]*CONSTANTS[201] -  ALGEBRAIC[175]*ALGEBRAIC[164]);
ALGEBRAIC[179] =  CONSTANTS[204]*( CONSTANTS[7]*ALGEBRAIC[177]+ CONSTANTS[9]*ALGEBRAIC[178]);
ALGEBRAIC[180] = 1.00000/(1.00000+exp(- (STATES[0] - 14.4800)/18.3400));
ALGEBRAIC[181] =  CONSTANTS[163]*ALGEBRAIC[180]*(STATES[0] - ALGEBRAIC[53]);
ALGEBRAIC[0] = (VOI>=CONSTANTS[12]&&VOI<=CONSTANTS[13]&&(VOI - CONSTANTS[12]) -  floor((VOI - CONSTANTS[12])/CONSTANTS[15])*CONSTANTS[15]<=CONSTANTS[16] ? CONSTANTS[14] : 0.00000);
ALGEBRAIC[183] = (STATES[6] - STATES[5])/2.00000;
RATES[5] = ( - ((ALGEBRAIC[66]+ALGEBRAIC[90]+ALGEBRAIC[96]+ALGEBRAIC[98]+ALGEBRAIC[181]+ALGEBRAIC[0]) -  2.00000*ALGEBRAIC[179])*CONSTANTS[32]*CONSTANTS[183])/( CONSTANTS[6]*CONSTANTS[184])+( ALGEBRAIC[183]*CONSTANTS[187])/CONSTANTS[184];
ALGEBRAIC[79] = ( 0.750000*CONSTANTS[169]*( STATES[6]*exp(ALGEBRAIC[12]) - CONSTANTS[3]))/CONSTANTS[176];
ALGEBRAIC[80] =  CONSTANTS[176]*(STATES[0] - CONSTANTS[158]);
ALGEBRAIC[81] = (- 1.00000e-07<=ALGEBRAIC[80]&&ALGEBRAIC[80]<=1.00000e-07 ?  ALGEBRAIC[79]*(1.00000 -  0.500000*ALGEBRAIC[80]) : ( ALGEBRAIC[79]*ALGEBRAIC[80])/(exp(ALGEBRAIC[80]) - 1.00000));
ALGEBRAIC[87] =  (1.00000 - ALGEBRAIC[82])*CONSTANTS[173]*ALGEBRAIC[81]*STATES[25]*( ALGEBRAIC[67]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[70]*STATES[33])+ ALGEBRAIC[82]*CONSTANTS[182]*ALGEBRAIC[81]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[72]*STATES[33]);
RATES[6] = ( - ALGEBRAIC[87]*CONSTANTS[32]*CONSTANTS[183])/( CONSTANTS[6]*CONSTANTS[187]) - ALGEBRAIC[183];
ALGEBRAIC[50] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log(CONSTANTS[1]/STATES[3]);
ALGEBRAIC[55] =  CONSTANTS[44]*STATES[11]+ CONSTANTS[151]*STATES[12];
ALGEBRAIC[56] =  CONSTANTS[44]*STATES[11]+ CONSTANTS[151]*STATES[14];
ALGEBRAIC[57] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[58] =  CONSTANTS[45]*(STATES[0] - ALGEBRAIC[50])*pow(STATES[10], 3.00000)*( (1.00000 - ALGEBRAIC[57])*ALGEBRAIC[55]*STATES[13]+ ALGEBRAIC[57]*ALGEBRAIC[56]*STATES[15]);
ALGEBRAIC[59] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[60] =  CONSTANTS[153]*(STATES[0] - ALGEBRAIC[50])*STATES[16]*( (1.00000 - ALGEBRAIC[59])*STATES[17]+ ALGEBRAIC[59]*STATES[18]);
ALGEBRAIC[127] = 1.00000/(1.00000+pow(CONSTANTS[117]/STATES[9], 2.00000));
ALGEBRAIC[100] = exp(( CONSTANTS[115]*STATES[0]*CONSTANTS[6])/( CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[107] = 1.00000+ (CONSTANTS[1]/CONSTANTS[108])*(1.00000+1.00000/ALGEBRAIC[100]);
ALGEBRAIC[108] = CONSTANTS[1]/( CONSTANTS[108]*ALGEBRAIC[100]*ALGEBRAIC[107]);
ALGEBRAIC[111] =  ALGEBRAIC[108]*CONSTANTS[112];
ALGEBRAIC[101] = 1.00000+ (STATES[3]/CONSTANTS[108])*(1.00000+ALGEBRAIC[100]);
ALGEBRAIC[102] = ( STATES[3]*ALGEBRAIC[100])/( CONSTANTS[108]*ALGEBRAIC[101]);
ALGEBRAIC[114] =  ALGEBRAIC[102]*CONSTANTS[112];
ALGEBRAIC[104] = 1.00000+ (STATES[3]/CONSTANTS[106])*(1.00000+STATES[3]/CONSTANTS[107]);
ALGEBRAIC[105] = ( STATES[3]*STATES[3])/( ALGEBRAIC[104]*CONSTANTS[106]*CONSTANTS[107]);
ALGEBRAIC[117] =  ALGEBRAIC[105]*ALGEBRAIC[102]*CONSTANTS[110];
ALGEBRAIC[118] =  ALGEBRAIC[108]*CONSTANTS[189]*CONSTANTS[110];
ALGEBRAIC[109] = 1.00000/ALGEBRAIC[107];
ALGEBRAIC[110] =  ALGEBRAIC[109]*CONSTANTS[111];
ALGEBRAIC[112] = ALGEBRAIC[110]+ALGEBRAIC[111];
ALGEBRAIC[99] = exp(( CONSTANTS[116]*STATES[0]*CONSTANTS[6])/( CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[103] = 1.00000/ALGEBRAIC[101];
ALGEBRAIC[113] = ( ALGEBRAIC[103]*CONSTANTS[111])/ALGEBRAIC[99];
ALGEBRAIC[115] = ALGEBRAIC[113]+ALGEBRAIC[114];
ALGEBRAIC[106] = 1.00000/ALGEBRAIC[104];
ALGEBRAIC[116] =  ALGEBRAIC[106]*STATES[9]*CONSTANTS[113];
ALGEBRAIC[119] =  CONSTANTS[192]*ALGEBRAIC[115]*(ALGEBRAIC[117]+ALGEBRAIC[116])+ CONSTANTS[193]*ALGEBRAIC[117]*(CONSTANTS[192]+ALGEBRAIC[112]);
ALGEBRAIC[120] =  CONSTANTS[191]*ALGEBRAIC[117]*(ALGEBRAIC[115]+CONSTANTS[193])+ ALGEBRAIC[115]*ALGEBRAIC[116]*(CONSTANTS[191]+ALGEBRAIC[118]);
ALGEBRAIC[121] =  CONSTANTS[191]*ALGEBRAIC[112]*(ALGEBRAIC[117]+ALGEBRAIC[116])+ ALGEBRAIC[118]*ALGEBRAIC[116]*(CONSTANTS[192]+ALGEBRAIC[112]);
ALGEBRAIC[122] =  CONSTANTS[192]*ALGEBRAIC[118]*(ALGEBRAIC[115]+CONSTANTS[193])+ ALGEBRAIC[112]*CONSTANTS[193]*(CONSTANTS[191]+ALGEBRAIC[118]);
ALGEBRAIC[123] = ALGEBRAIC[119]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[124] = ALGEBRAIC[120]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[125] = ALGEBRAIC[121]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[126] = ALGEBRAIC[122]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[128] = ( 3.00000*( ALGEBRAIC[126]*ALGEBRAIC[117] -  ALGEBRAIC[123]*ALGEBRAIC[118])+ ALGEBRAIC[125]*ALGEBRAIC[114]) -  ALGEBRAIC[124]*ALGEBRAIC[111];
ALGEBRAIC[129] =  ALGEBRAIC[124]*CONSTANTS[192] -  ALGEBRAIC[123]*CONSTANTS[191];
ALGEBRAIC[130] =  0.800000*CONSTANTS[194]*ALGEBRAIC[127]*( CONSTANTS[7]*ALGEBRAIC[128]+ CONSTANTS[8]*ALGEBRAIC[129]);
ALGEBRAIC[182] = ( CONSTANTS[142]*CONSTANTS[169]*( STATES[3]*exp(ALGEBRAIC[12]) - CONSTANTS[1]))/CONSTANTS[178];
ALGEBRAIC[184] =  CONSTANTS[178]*(STATES[0] - CONSTANTS[164]);
ALGEBRAIC[185] = (- 1.00000e-07<=ALGEBRAIC[184]&&ALGEBRAIC[184]<=1.00000e-07 ?  ALGEBRAIC[182]*(1.00000 -  0.500000*ALGEBRAIC[184]) : ( ALGEBRAIC[182]*ALGEBRAIC[184])/(exp(ALGEBRAIC[184]) - 1.00000));
ALGEBRAIC[187] = (STATES[4] - STATES[3])/2.00000;
RATES[3] = ( - (ALGEBRAIC[58]+ALGEBRAIC[60]+ 3.00000*ALGEBRAIC[130]+ 3.00000*ALGEBRAIC[179]+ALGEBRAIC[185])*CONSTANTS[183]*CONSTANTS[32])/( CONSTANTS[6]*CONSTANTS[184])+( ALGEBRAIC[187]*CONSTANTS[187])/CONSTANTS[184];
ALGEBRAIC[76] = ( 0.750000*CONSTANTS[169]*( STATES[4]*exp(ALGEBRAIC[12]) - CONSTANTS[1]))/CONSTANTS[175];
ALGEBRAIC[77] =  CONSTANTS[175]*(STATES[0] - CONSTANTS[158]);
ALGEBRAIC[78] = (- 1.00000e-07<=ALGEBRAIC[77]&&ALGEBRAIC[77]<=1.00000e-07 ?  ALGEBRAIC[76]*(1.00000 -  0.500000*ALGEBRAIC[77]) : ( ALGEBRAIC[76]*ALGEBRAIC[77])/(exp(ALGEBRAIC[77]) - 1.00000));
ALGEBRAIC[84] =  (1.00000 - ALGEBRAIC[82])*CONSTANTS[172]*ALGEBRAIC[78]*STATES[25]*( ALGEBRAIC[67]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[70]*STATES[33])+ ALGEBRAIC[82]*CONSTANTS[181]*ALGEBRAIC[78]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[72]*STATES[33]);
ALGEBRAIC[157] = 1.00000/(1.00000+pow(CONSTANTS[117]/STATES[2], 2.00000));
ALGEBRAIC[137] = 1.00000+ (CONSTANTS[1]/CONSTANTS[108])*(1.00000+1.00000/ALGEBRAIC[100]);
ALGEBRAIC[138] = CONSTANTS[1]/( CONSTANTS[108]*ALGEBRAIC[100]*ALGEBRAIC[137]);
ALGEBRAIC[141] =  ALGEBRAIC[138]*CONSTANTS[112];
ALGEBRAIC[131] = 1.00000+ (STATES[4]/CONSTANTS[108])*(1.00000+ALGEBRAIC[100]);
ALGEBRAIC[132] = ( STATES[4]*ALGEBRAIC[100])/( CONSTANTS[108]*ALGEBRAIC[131]);
ALGEBRAIC[144] =  ALGEBRAIC[132]*CONSTANTS[112];
ALGEBRAIC[134] = 1.00000+ (STATES[4]/CONSTANTS[106])*(1.00000+STATES[4]/CONSTANTS[107]);
ALGEBRAIC[135] = ( STATES[4]*STATES[4])/( ALGEBRAIC[134]*CONSTANTS[106]*CONSTANTS[107]);
ALGEBRAIC[147] =  ALGEBRAIC[135]*ALGEBRAIC[132]*CONSTANTS[110];
ALGEBRAIC[148] =  ALGEBRAIC[138]*CONSTANTS[196]*CONSTANTS[110];
ALGEBRAIC[139] = 1.00000/ALGEBRAIC[137];
ALGEBRAIC[140] =  ALGEBRAIC[139]*CONSTANTS[111];
ALGEBRAIC[142] = ALGEBRAIC[140]+ALGEBRAIC[141];
ALGEBRAIC[133] = 1.00000/ALGEBRAIC[131];
ALGEBRAIC[143] = ( ALGEBRAIC[133]*CONSTANTS[111])/ALGEBRAIC[99];
ALGEBRAIC[145] = ALGEBRAIC[143]+ALGEBRAIC[144];
ALGEBRAIC[136] = 1.00000/ALGEBRAIC[134];
ALGEBRAIC[146] =  ALGEBRAIC[136]*STATES[2]*CONSTANTS[113];
ALGEBRAIC[149] =  CONSTANTS[199]*ALGEBRAIC[145]*(ALGEBRAIC[147]+ALGEBRAIC[146])+ CONSTANTS[200]*ALGEBRAIC[147]*(CONSTANTS[199]+ALGEBRAIC[142]);
ALGEBRAIC[150] =  CONSTANTS[198]*ALGEBRAIC[147]*(ALGEBRAIC[145]+CONSTANTS[200])+ ALGEBRAIC[145]*ALGEBRAIC[146]*(CONSTANTS[198]+ALGEBRAIC[148]);
ALGEBRAIC[151] =  CONSTANTS[198]*ALGEBRAIC[142]*(ALGEBRAIC[147]+ALGEBRAIC[146])+ ALGEBRAIC[148]*ALGEBRAIC[146]*(CONSTANTS[199]+ALGEBRAIC[142]);
ALGEBRAIC[152] =  CONSTANTS[199]*ALGEBRAIC[148]*(ALGEBRAIC[145]+CONSTANTS[200])+ ALGEBRAIC[142]*CONSTANTS[200]*(CONSTANTS[198]+ALGEBRAIC[148]);
ALGEBRAIC[153] = ALGEBRAIC[149]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[154] = ALGEBRAIC[150]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[155] = ALGEBRAIC[151]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[156] = ALGEBRAIC[152]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[158] = ( 3.00000*( ALGEBRAIC[156]*ALGEBRAIC[147] -  ALGEBRAIC[153]*ALGEBRAIC[148])+ ALGEBRAIC[155]*ALGEBRAIC[144]) -  ALGEBRAIC[154]*ALGEBRAIC[141];
ALGEBRAIC[159] =  ALGEBRAIC[154]*CONSTANTS[199] -  ALGEBRAIC[153]*CONSTANTS[198];
ALGEBRAIC[160] =  0.200000*CONSTANTS[194]*ALGEBRAIC[157]*( CONSTANTS[7]*ALGEBRAIC[158]+ CONSTANTS[8]*ALGEBRAIC[159]);
RATES[4] = ( - (ALGEBRAIC[84]+ 3.00000*ALGEBRAIC[160])*CONSTANTS[32]*CONSTANTS[183])/( CONSTANTS[6]*CONSTANTS[187]) - ALGEBRAIC[187];
ALGEBRAIC[190] = ( CONSTANTS[144]*STATES[9])/(CONSTANTS[145]+STATES[9]);
ALGEBRAIC[186] = ( CONSTANTS[143]*4.00000*CONSTANTS[169]*( STATES[9]*exp( 2.00000*ALGEBRAIC[12]) -  0.341000*CONSTANTS[2]))/CONSTANTS[179];
ALGEBRAIC[188] =  CONSTANTS[179]*(STATES[0] - CONSTANTS[165]);
ALGEBRAIC[189] = (- 1.00000e-07<=ALGEBRAIC[188]&&ALGEBRAIC[188]<=1.00000e-07 ?  ALGEBRAIC[186]*(1.00000 -  0.500000*ALGEBRAIC[188]) : ( ALGEBRAIC[186]*ALGEBRAIC[188])/(exp(ALGEBRAIC[188]) - 1.00000));
RATES[0] = - (ALGEBRAIC[58]+ALGEBRAIC[60]+ALGEBRAIC[66]+ALGEBRAIC[83]+ALGEBRAIC[84]+ALGEBRAIC[87]+ALGEBRAIC[90]+ALGEBRAIC[96]+ALGEBRAIC[98]+ALGEBRAIC[130]+ALGEBRAIC[160]+ALGEBRAIC[179]+ALGEBRAIC[185]+ALGEBRAIC[181]+ALGEBRAIC[190]+ALGEBRAIC[189]+ALGEBRAIC[0]);
ALGEBRAIC[191] = (STATES[2] - STATES[9])/0.200000;
ALGEBRAIC[192] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[193] =  CONSTANTS[147]*( (1.00000 - ALGEBRAIC[192])*STATES[47]+ ALGEBRAIC[192]*STATES[48]);
ALGEBRAIC[46] = 1.00000/(1.00000+( CONSTANTS[26]*CONSTANTS[27])/pow(CONSTANTS[27]+STATES[2], 2.00000)+( CONSTANTS[28]*CONSTANTS[29])/pow(CONSTANTS[29]+STATES[2], 2.00000));
RATES[2] =  ALGEBRAIC[46]*((( - (ALGEBRAIC[83] -  2.00000*ALGEBRAIC[160])*CONSTANTS[32]*CONSTANTS[183])/( 2.00000*CONSTANTS[6]*CONSTANTS[187])+( ALGEBRAIC[193]*CONSTANTS[186])/CONSTANTS[187]) - ALGEBRAIC[191]);
ALGEBRAIC[194] = ( CONSTANTS[168]*0.00437500*STATES[9])/(STATES[9]+0.000920000);
ALGEBRAIC[195] = ( CONSTANTS[168]*2.75000*0.00437500*STATES[9])/((STATES[9]+0.000920000) - 0.000170000);
ALGEBRAIC[196] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[197] = ( 0.00393750*STATES[7])/15.0000;
ALGEBRAIC[198] =  CONSTANTS[148]*(( (1.00000 - ALGEBRAIC[196])*ALGEBRAIC[194]+ ALGEBRAIC[196]*ALGEBRAIC[195]) - ALGEBRAIC[197]);
ALGEBRAIC[44] = 1.00000/(1.00000+( CONSTANTS[150]*CONSTANTS[23])/pow(CONSTANTS[23]+STATES[9], 2.00000)+( CONSTANTS[24]*CONSTANTS[25])/pow(CONSTANTS[25]+STATES[9], 2.00000));
RATES[9] =  ALGEBRAIC[44]*((( - ((ALGEBRAIC[190]+ALGEBRAIC[189]) -  2.00000*ALGEBRAIC[130])*CONSTANTS[32]*CONSTANTS[183])/( 2.00000*CONSTANTS[6]*CONSTANTS[184]) - ( ALGEBRAIC[198]*CONSTANTS[185])/CONSTANTS[184])+( ALGEBRAIC[191]*CONSTANTS[187])/CONSTANTS[184]);
ALGEBRAIC[199] = (STATES[7] - STATES[8])/100.000;
RATES[7] = ALGEBRAIC[198] - ( ALGEBRAIC[199]*CONSTANTS[186])/CONSTANTS[185];
ALGEBRAIC[48] = 1.00000/(1.00000+( CONSTANTS[30]*CONSTANTS[31])/pow(CONSTANTS[31]+STATES[8], 2.00000));
RATES[8] =  ALGEBRAIC[48]*(ALGEBRAIC[199] - ALGEBRAIC[193]);
}
void
computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[3] = 1.00000/(1.00000+exp((STATES[0]+87.6100)/7.48800));
ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[0]+93.8100)/7.48800));
ALGEBRAIC[1] = 1.00000/(1.00000+exp(- (STATES[0]+CONSTANTS[34])/CONSTANTS[35]));
ALGEBRAIC[13] = 1.00000/( CONSTANTS[38]*exp((STATES[0]+CONSTANTS[36])/CONSTANTS[37])+ CONSTANTS[39]*exp(- (STATES[0]+CONSTANTS[40])/CONSTANTS[41]));
ALGEBRAIC[2] = 1.00000/(1.00000+exp(((STATES[0]+CONSTANTS[42]) - CONSTANTS[46])/CONSTANTS[43]));
ALGEBRAIC[14] = 1.00000/( 1.43200e-05*exp(- ((STATES[0]+1.19600) - CONSTANTS[46])/6.28500)+ 6.14900*exp(((STATES[0]+0.509600) - CONSTANTS[46])/20.2700));
ALGEBRAIC[15] = 1.00000/( 0.00979400*exp(- ((STATES[0]+17.9500) - CONSTANTS[46])/28.0500)+ 0.334300*exp(((STATES[0]+5.73000) - CONSTANTS[46])/56.6600));
ALGEBRAIC[5] = 1.00000/(1.00000+exp(- (STATES[0] - 14.3400)/14.8200));
ALGEBRAIC[17] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- (STATES[0] - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[0]+100.000)/29.3814)));
ALGEBRAIC[7] = 1.00000/(1.00000+exp(- (STATES[0]+3.94000)/4.23000));
ALGEBRAIC[21] = 0.600000+1.00000/(exp( - 0.0500000*(STATES[0]+6.00000))+exp( 0.0900000*(STATES[0]+14.0000)));
ALGEBRAIC[8] = 1.00000/(1.00000+exp((STATES[0]+19.5800)/3.69600));
ALGEBRAIC[22] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[0]+20.0000)/10.0000)+ 0.00450000*exp((STATES[0]+20.0000)/10.0000));
ALGEBRAIC[23] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[0]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[0]+5.00000)/6.00000));
ALGEBRAIC[19] = ALGEBRAIC[8];
ALGEBRAIC[9] =  STATES[30]*1.00000;
ALGEBRAIC[20] = 1.00000/(CONSTANTS[51]/ALGEBRAIC[9]+pow(1.00000+CONSTANTS[50]/STATES[2], 4.00000));
ALGEBRAIC[10] = 1.00000/(1.00000+exp(- (STATES[0]+11.6000)/8.93200));
ALGEBRAIC[25] = CONSTANTS[104]+1.00000/( 0.000232600*exp((STATES[0]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[0]+210.000)/230.000));
ALGEBRAIC[11] = 1.00000/(1.00000+exp(- (STATES[0]+ 2.55380*CONSTANTS[3]+144.590)/( 1.56920*CONSTANTS[3]+3.81150)));
ALGEBRAIC[26] = 122.200/(exp(- (STATES[0]+127.200)/20.3600)+exp((STATES[0]+236.800)/69.3300));
ALGEBRAIC[36] = ( CONSTANTS[20]*(1.00000 - STATES[1]))/(1.00000+CONSTANTS[21]/STATES[2]);
ALGEBRAIC[16] = ALGEBRAIC[2];
ALGEBRAIC[27] = 2.03800+1.00000/( 0.0213600*exp(- ((STATES[0]+100.600) - CONSTANTS[46])/8.28100)+ 0.305200*exp(((STATES[0]+0.994100) - CONSTANTS[46])/38.4500));
ALGEBRAIC[31] = 1.00000/(1.00000+exp(- (STATES[0] - 24.3400)/14.8200));
ALGEBRAIC[32] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[0] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[0] - 4.00000)/7.00000));
ALGEBRAIC[33] = 100.000+1.00000/( 0.000120000*exp(- STATES[0]/3.00000)+ 0.000120000*exp(STATES[0]/7.00000));
ALGEBRAIC[34] =  2.50000*ALGEBRAIC[22];
ALGEBRAIC[24] = ALGEBRAIC[10];
ALGEBRAIC[35] = 1.00000/( 0.0100000*exp((STATES[0] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[0]+66.5400)/31.0000));
ALGEBRAIC[28] = 1.00000/(1.00000+exp(((STATES[0]+89.1000) - CONSTANTS[46])/6.08600));
ALGEBRAIC[37] =  3.00000*ALGEBRAIC[15];
ALGEBRAIC[38] =  1.46000*ALGEBRAIC[27];
ALGEBRAIC[29] = 1.00000/(1.00000+exp(- (STATES[0]+42.8500)/5.26400));
ALGEBRAIC[39] = ALGEBRAIC[13];
ALGEBRAIC[41] =  2.50000*ALGEBRAIC[32];
ALGEBRAIC[6] = 1.00000/(1.00000+exp((STATES[0]+43.9400)/5.71100));
ALGEBRAIC[18] = (CONSTANTS[0]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[0]+70.0000)/5.00000)) : 1.00000);
ALGEBRAIC[30] = 4.56200+1.00000/( 0.393300*exp(- (STATES[0]+100.000)/100.000)+ 0.0800400*exp((STATES[0]+50.0000)/16.5900));
ALGEBRAIC[43] =  ALGEBRAIC[30]*ALGEBRAIC[18];
ALGEBRAIC[40] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[0]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[0]+114.100)/8.07900));
ALGEBRAIC[45] =  ALGEBRAIC[40]*ALGEBRAIC[18];
ALGEBRAIC[47] = 1.35400+0.000100000/(exp((STATES[0] - 167.400)/15.8900)+exp(- (STATES[0] - 12.2300)/0.215400));
ALGEBRAIC[49] = 1.00000 - 0.500000/(1.00000+exp((STATES[0]+70.0000)/20.0000));
ALGEBRAIC[51] =  ALGEBRAIC[47]*ALGEBRAIC[49]*ALGEBRAIC[43];
ALGEBRAIC[52] =  ALGEBRAIC[47]*ALGEBRAIC[49]*ALGEBRAIC[45];
ALGEBRAIC[67] =  CONSTANTS[155]*STATES[26]+ CONSTANTS[170]*STATES[27];
ALGEBRAIC[68] = 0.300000+0.600000/(1.00000+exp((STATES[0] - 10.0000)/10.0000));
ALGEBRAIC[69] = 1.00000 - ALGEBRAIC[68];
ALGEBRAIC[70] =  ALGEBRAIC[68]*STATES[28]+ ALGEBRAIC[69]*STATES[29];
ALGEBRAIC[71] =  CONSTANTS[155]*STATES[31]+ CONSTANTS[170]*STATES[27];
ALGEBRAIC[72] =  ALGEBRAIC[68]*STATES[32]+ ALGEBRAIC[69]*STATES[29];
ALGEBRAIC[12] =  STATES[0]*CONSTANTS[149];
ALGEBRAIC[73] = ( 4.00000*CONSTANTS[169]*( STATES[2]*exp( 2.00000*ALGEBRAIC[12]) -  0.341000*CONSTANTS[2]))/CONSTANTS[174];
ALGEBRAIC[74] =  CONSTANTS[174]*(STATES[0] - CONSTANTS[158]);
ALGEBRAIC[75] = (- 1.00000e-07<=ALGEBRAIC[74]&&ALGEBRAIC[74]<=1.00000e-07 ?  ALGEBRAIC[73]*(1.00000 -  0.500000*ALGEBRAIC[74]) : ( ALGEBRAIC[73]*ALGEBRAIC[74])/(exp(ALGEBRAIC[74]) - 1.00000));
ALGEBRAIC[42] = ALGEBRAIC[36]+STATES[1];
ALGEBRAIC[82] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[83] =  (1.00000 - ALGEBRAIC[82])*CONSTANTS[156]*ALGEBRAIC[75]*STATES[25]*( ALGEBRAIC[67]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[70]*STATES[33])+ ALGEBRAIC[82]*CONSTANTS[171]*ALGEBRAIC[75]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[72]*STATES[33]);
ALGEBRAIC[85] = ( CONSTANTS[166]*- ALGEBRAIC[83])/(1.00000+ 1.00000*pow(1.50000/STATES[8], 8.00000));
ALGEBRAIC[88] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[85]*1.70000 : ALGEBRAIC[85]);
ALGEBRAIC[91] = CONSTANTS[146]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[94] = (ALGEBRAIC[91]<0.00100000 ? 0.00100000 : ALGEBRAIC[91]);
ALGEBRAIC[86] = ( CONSTANTS[180]*- ALGEBRAIC[83])/(1.00000+pow(1.50000/STATES[8], 8.00000));
ALGEBRAIC[89] = (CONSTANTS[0]==2.00000 ?  ALGEBRAIC[86]*1.70000 : ALGEBRAIC[86]);
ALGEBRAIC[92] = CONSTANTS[167]/(1.00000+0.0123000/STATES[8]);
ALGEBRAIC[95] = (ALGEBRAIC[92]<0.00100000 ? 0.00100000 : ALGEBRAIC[92]);
ALGEBRAIC[53] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log(CONSTANTS[3]/STATES[5]);
ALGEBRAIC[61] = 1.00000/(1.00000+exp((STATES[0] - 213.600)/151.200));
ALGEBRAIC[62] = 1.00000 - ALGEBRAIC[61];
ALGEBRAIC[63] =  ALGEBRAIC[61]*STATES[20]+ ALGEBRAIC[62]*STATES[21];
ALGEBRAIC[64] =  ALGEBRAIC[61]*STATES[23]+ ALGEBRAIC[62]*STATES[24];
ALGEBRAIC[65] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[66] =  CONSTANTS[154]*(STATES[0] - ALGEBRAIC[53])*( (1.00000 - ALGEBRAIC[65])*STATES[19]*ALGEBRAIC[63]+ ALGEBRAIC[65]*STATES[22]*ALGEBRAIC[64]);
ALGEBRAIC[90] =  CONSTANTS[159]* pow((CONSTANTS[3]/5.40000), 1.0 / 2)*STATES[38]*(STATES[0] - ALGEBRAIC[53]);
ALGEBRAIC[54] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log((CONSTANTS[3]+ CONSTANTS[33]*CONSTANTS[1])/(STATES[5]+ CONSTANTS[33]*STATES[3]));
ALGEBRAIC[93] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[9], 1.40000));
ALGEBRAIC[96] =  CONSTANTS[160]*ALGEBRAIC[93]*STATES[44]*STATES[45]*(STATES[0] - ALGEBRAIC[54]);
ALGEBRAIC[97] = 1.00000/(1.00000+exp(((STATES[0]+105.800) -  2.60000*CONSTANTS[3])/9.49300));
ALGEBRAIC[98] =  CONSTANTS[161]* pow(CONSTANTS[3], 1.0 / 2)*ALGEBRAIC[97]*STATES[46]*(STATES[0] - ALGEBRAIC[53]);
ALGEBRAIC[162] =  CONSTANTS[128]*exp(( (1.00000 - CONSTANTS[129])*STATES[0]*CONSTANTS[6])/( 3.00000*CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[166] = ( CONSTANTS[123]*pow(CONSTANTS[3]/CONSTANTS[131], 2.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[162], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[131], 2.00000)) - 1.00000);
ALGEBRAIC[163] = CONSTANTS[136]/(1.00000+CONSTANTS[135]/CONSTANTS[137]+STATES[3]/CONSTANTS[138]+STATES[5]/CONSTANTS[139]);
ALGEBRAIC[167] = ( CONSTANTS[124]*ALGEBRAIC[163]*CONSTANTS[135])/(1.00000+CONSTANTS[133]/CONSTANTS[134]);
ALGEBRAIC[161] =  CONSTANTS[127]*exp(( CONSTANTS[129]*STATES[0]*CONSTANTS[6])/( 3.00000*CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[164] = ( CONSTANTS[119]*pow(STATES[3]/ALGEBRAIC[161], 3.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[161], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[130], 2.00000)) - 1.00000);
ALGEBRAIC[165] = ( CONSTANTS[122]*pow(CONSTANTS[1]/ALGEBRAIC[162], 3.00000))/((pow(1.00000+CONSTANTS[1]/ALGEBRAIC[162], 3.00000)+pow(1.00000+CONSTANTS[3]/CONSTANTS[131], 2.00000)) - 1.00000);
ALGEBRAIC[168] = ( CONSTANTS[126]*pow(STATES[5]/CONSTANTS[130], 2.00000))/((pow(1.00000+STATES[3]/ALGEBRAIC[161], 3.00000)+pow(1.00000+STATES[5]/CONSTANTS[130], 2.00000)) - 1.00000);
ALGEBRAIC[169] =  CONSTANTS[203]*ALGEBRAIC[164]*CONSTANTS[202]+ ALGEBRAIC[165]*ALGEBRAIC[168]*ALGEBRAIC[167]+ CONSTANTS[202]*ALGEBRAIC[168]*ALGEBRAIC[167]+ ALGEBRAIC[167]*ALGEBRAIC[164]*CONSTANTS[202];
ALGEBRAIC[170] =  ALGEBRAIC[165]*CONSTANTS[201]*ALGEBRAIC[168]+ ALGEBRAIC[164]*CONSTANTS[202]*ALGEBRAIC[166]+ ALGEBRAIC[166]*CONSTANTS[201]*ALGEBRAIC[168]+ CONSTANTS[202]*ALGEBRAIC[166]*ALGEBRAIC[168];
ALGEBRAIC[171] =  CONSTANTS[202]*ALGEBRAIC[166]*CONSTANTS[203]+ ALGEBRAIC[167]*ALGEBRAIC[165]*CONSTANTS[201]+ ALGEBRAIC[165]*CONSTANTS[201]*CONSTANTS[203]+ ALGEBRAIC[166]*CONSTANTS[203]*CONSTANTS[201];
ALGEBRAIC[172] =  ALGEBRAIC[168]*ALGEBRAIC[167]*ALGEBRAIC[165]+ ALGEBRAIC[166]*CONSTANTS[203]*ALGEBRAIC[164]+ ALGEBRAIC[165]*CONSTANTS[203]*ALGEBRAIC[164]+ ALGEBRAIC[167]*ALGEBRAIC[165]*ALGEBRAIC[164];
ALGEBRAIC[173] = ALGEBRAIC[169]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[174] = ALGEBRAIC[170]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[177] =  3.00000*( ALGEBRAIC[173]*ALGEBRAIC[166] -  ALGEBRAIC[174]*ALGEBRAIC[167]);
ALGEBRAIC[175] = ALGEBRAIC[171]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[176] = ALGEBRAIC[172]/(ALGEBRAIC[169]+ALGEBRAIC[170]+ALGEBRAIC[171]+ALGEBRAIC[172]);
ALGEBRAIC[178] =  2.00000*( ALGEBRAIC[176]*CONSTANTS[201] -  ALGEBRAIC[175]*ALGEBRAIC[164]);
ALGEBRAIC[179] =  CONSTANTS[204]*( CONSTANTS[7]*ALGEBRAIC[177]+ CONSTANTS[9]*ALGEBRAIC[178]);
ALGEBRAIC[180] = 1.00000/(1.00000+exp(- (STATES[0] - 14.4800)/18.3400));
ALGEBRAIC[181] =  CONSTANTS[163]*ALGEBRAIC[180]*(STATES[0] - ALGEBRAIC[53]);
ALGEBRAIC[0] = (VOI>=CONSTANTS[12]&&VOI<=CONSTANTS[13]&&(VOI - CONSTANTS[12]) -  floor((VOI - CONSTANTS[12])/CONSTANTS[15])*CONSTANTS[15]<=CONSTANTS[16] ? CONSTANTS[14] : 0.00000);
ALGEBRAIC[183] = (STATES[6] - STATES[5])/2.00000;
ALGEBRAIC[79] = ( 0.750000*CONSTANTS[169]*( STATES[6]*exp(ALGEBRAIC[12]) - CONSTANTS[3]))/CONSTANTS[176];
ALGEBRAIC[80] =  CONSTANTS[176]*(STATES[0] - CONSTANTS[158]);
ALGEBRAIC[81] = (- 1.00000e-07<=ALGEBRAIC[80]&&ALGEBRAIC[80]<=1.00000e-07 ?  ALGEBRAIC[79]*(1.00000 -  0.500000*ALGEBRAIC[80]) : ( ALGEBRAIC[79]*ALGEBRAIC[80])/(exp(ALGEBRAIC[80]) - 1.00000));
ALGEBRAIC[87] =  (1.00000 - ALGEBRAIC[82])*CONSTANTS[173]*ALGEBRAIC[81]*STATES[25]*( ALGEBRAIC[67]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[70]*STATES[33])+ ALGEBRAIC[82]*CONSTANTS[182]*ALGEBRAIC[81]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[72]*STATES[33]);
ALGEBRAIC[50] =  (( CONSTANTS[4]*CONSTANTS[5])/CONSTANTS[6])*log(CONSTANTS[1]/STATES[3]);
ALGEBRAIC[55] =  CONSTANTS[44]*STATES[11]+ CONSTANTS[151]*STATES[12];
ALGEBRAIC[56] =  CONSTANTS[44]*STATES[11]+ CONSTANTS[151]*STATES[14];
ALGEBRAIC[57] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[58] =  CONSTANTS[45]*(STATES[0] - ALGEBRAIC[50])*pow(STATES[10], 3.00000)*( (1.00000 - ALGEBRAIC[57])*ALGEBRAIC[55]*STATES[13]+ ALGEBRAIC[57]*ALGEBRAIC[56]*STATES[15]);
ALGEBRAIC[59] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[60] =  CONSTANTS[153]*(STATES[0] - ALGEBRAIC[50])*STATES[16]*( (1.00000 - ALGEBRAIC[59])*STATES[17]+ ALGEBRAIC[59]*STATES[18]);
ALGEBRAIC[127] = 1.00000/(1.00000+pow(CONSTANTS[117]/STATES[9], 2.00000));
ALGEBRAIC[100] = exp(( CONSTANTS[115]*STATES[0]*CONSTANTS[6])/( CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[107] = 1.00000+ (CONSTANTS[1]/CONSTANTS[108])*(1.00000+1.00000/ALGEBRAIC[100]);
ALGEBRAIC[108] = CONSTANTS[1]/( CONSTANTS[108]*ALGEBRAIC[100]*ALGEBRAIC[107]);
ALGEBRAIC[111] =  ALGEBRAIC[108]*CONSTANTS[112];
ALGEBRAIC[101] = 1.00000+ (STATES[3]/CONSTANTS[108])*(1.00000+ALGEBRAIC[100]);
ALGEBRAIC[102] = ( STATES[3]*ALGEBRAIC[100])/( CONSTANTS[108]*ALGEBRAIC[101]);
ALGEBRAIC[114] =  ALGEBRAIC[102]*CONSTANTS[112];
ALGEBRAIC[104] = 1.00000+ (STATES[3]/CONSTANTS[106])*(1.00000+STATES[3]/CONSTANTS[107]);
ALGEBRAIC[105] = ( STATES[3]*STATES[3])/( ALGEBRAIC[104]*CONSTANTS[106]*CONSTANTS[107]);
ALGEBRAIC[117] =  ALGEBRAIC[105]*ALGEBRAIC[102]*CONSTANTS[110];
ALGEBRAIC[118] =  ALGEBRAIC[108]*CONSTANTS[189]*CONSTANTS[110];
ALGEBRAIC[109] = 1.00000/ALGEBRAIC[107];
ALGEBRAIC[110] =  ALGEBRAIC[109]*CONSTANTS[111];
ALGEBRAIC[112] = ALGEBRAIC[110]+ALGEBRAIC[111];
ALGEBRAIC[99] = exp(( CONSTANTS[116]*STATES[0]*CONSTANTS[6])/( CONSTANTS[4]*CONSTANTS[5]));
ALGEBRAIC[103] = 1.00000/ALGEBRAIC[101];
ALGEBRAIC[113] = ( ALGEBRAIC[103]*CONSTANTS[111])/ALGEBRAIC[99];
ALGEBRAIC[115] = ALGEBRAIC[113]+ALGEBRAIC[114];
ALGEBRAIC[106] = 1.00000/ALGEBRAIC[104];
ALGEBRAIC[116] =  ALGEBRAIC[106]*STATES[9]*CONSTANTS[113];
ALGEBRAIC[119] =  CONSTANTS[192]*ALGEBRAIC[115]*(ALGEBRAIC[117]+ALGEBRAIC[116])+ CONSTANTS[193]*ALGEBRAIC[117]*(CONSTANTS[192]+ALGEBRAIC[112]);
ALGEBRAIC[120] =  CONSTANTS[191]*ALGEBRAIC[117]*(ALGEBRAIC[115]+CONSTANTS[193])+ ALGEBRAIC[115]*ALGEBRAIC[116]*(CONSTANTS[191]+ALGEBRAIC[118]);
ALGEBRAIC[121] =  CONSTANTS[191]*ALGEBRAIC[112]*(ALGEBRAIC[117]+ALGEBRAIC[116])+ ALGEBRAIC[118]*ALGEBRAIC[116]*(CONSTANTS[192]+ALGEBRAIC[112]);
ALGEBRAIC[122] =  CONSTANTS[192]*ALGEBRAIC[118]*(ALGEBRAIC[115]+CONSTANTS[193])+ ALGEBRAIC[112]*CONSTANTS[193]*(CONSTANTS[191]+ALGEBRAIC[118]);
ALGEBRAIC[123] = ALGEBRAIC[119]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[124] = ALGEBRAIC[120]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[125] = ALGEBRAIC[121]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[126] = ALGEBRAIC[122]/(ALGEBRAIC[119]+ALGEBRAIC[120]+ALGEBRAIC[121]+ALGEBRAIC[122]);
ALGEBRAIC[128] = ( 3.00000*( ALGEBRAIC[126]*ALGEBRAIC[117] -  ALGEBRAIC[123]*ALGEBRAIC[118])+ ALGEBRAIC[125]*ALGEBRAIC[114]) -  ALGEBRAIC[124]*ALGEBRAIC[111];
ALGEBRAIC[129] =  ALGEBRAIC[124]*CONSTANTS[192] -  ALGEBRAIC[123]*CONSTANTS[191];
ALGEBRAIC[130] =  0.800000*CONSTANTS[194]*ALGEBRAIC[127]*( CONSTANTS[7]*ALGEBRAIC[128]+ CONSTANTS[8]*ALGEBRAIC[129]);
ALGEBRAIC[182] = ( CONSTANTS[142]*CONSTANTS[169]*( STATES[3]*exp(ALGEBRAIC[12]) - CONSTANTS[1]))/CONSTANTS[178];
ALGEBRAIC[184] =  CONSTANTS[178]*(STATES[0] - CONSTANTS[164]);
ALGEBRAIC[185] = (- 1.00000e-07<=ALGEBRAIC[184]&&ALGEBRAIC[184]<=1.00000e-07 ?  ALGEBRAIC[182]*(1.00000 -  0.500000*ALGEBRAIC[184]) : ( ALGEBRAIC[182]*ALGEBRAIC[184])/(exp(ALGEBRAIC[184]) - 1.00000));
ALGEBRAIC[187] = (STATES[4] - STATES[3])/2.00000;
ALGEBRAIC[76] = ( 0.750000*CONSTANTS[169]*( STATES[4]*exp(ALGEBRAIC[12]) - CONSTANTS[1]))/CONSTANTS[175];
ALGEBRAIC[77] =  CONSTANTS[175]*(STATES[0] - CONSTANTS[158]);
ALGEBRAIC[78] = (- 1.00000e-07<=ALGEBRAIC[77]&&ALGEBRAIC[77]<=1.00000e-07 ?  ALGEBRAIC[76]*(1.00000 -  0.500000*ALGEBRAIC[77]) : ( ALGEBRAIC[76]*ALGEBRAIC[77])/(exp(ALGEBRAIC[77]) - 1.00000));
ALGEBRAIC[84] =  (1.00000 - ALGEBRAIC[82])*CONSTANTS[172]*ALGEBRAIC[78]*STATES[25]*( ALGEBRAIC[67]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[70]*STATES[33])+ ALGEBRAIC[82]*CONSTANTS[181]*ALGEBRAIC[78]*STATES[25]*( ALGEBRAIC[71]*(1.00000 - STATES[33])+ STATES[30]*ALGEBRAIC[72]*STATES[33]);
ALGEBRAIC[157] = 1.00000/(1.00000+pow(CONSTANTS[117]/STATES[2], 2.00000));
ALGEBRAIC[137] = 1.00000+ (CONSTANTS[1]/CONSTANTS[108])*(1.00000+1.00000/ALGEBRAIC[100]);
ALGEBRAIC[138] = CONSTANTS[1]/( CONSTANTS[108]*ALGEBRAIC[100]*ALGEBRAIC[137]);
ALGEBRAIC[141] =  ALGEBRAIC[138]*CONSTANTS[112];
ALGEBRAIC[131] = 1.00000+ (STATES[4]/CONSTANTS[108])*(1.00000+ALGEBRAIC[100]);
ALGEBRAIC[132] = ( STATES[4]*ALGEBRAIC[100])/( CONSTANTS[108]*ALGEBRAIC[131]);
ALGEBRAIC[144] =  ALGEBRAIC[132]*CONSTANTS[112];
ALGEBRAIC[134] = 1.00000+ (STATES[4]/CONSTANTS[106])*(1.00000+STATES[4]/CONSTANTS[107]);
ALGEBRAIC[135] = ( STATES[4]*STATES[4])/( ALGEBRAIC[134]*CONSTANTS[106]*CONSTANTS[107]);
ALGEBRAIC[147] =  ALGEBRAIC[135]*ALGEBRAIC[132]*CONSTANTS[110];
ALGEBRAIC[148] =  ALGEBRAIC[138]*CONSTANTS[196]*CONSTANTS[110];
ALGEBRAIC[139] = 1.00000/ALGEBRAIC[137];
ALGEBRAIC[140] =  ALGEBRAIC[139]*CONSTANTS[111];
ALGEBRAIC[142] = ALGEBRAIC[140]+ALGEBRAIC[141];
ALGEBRAIC[133] = 1.00000/ALGEBRAIC[131];
ALGEBRAIC[143] = ( ALGEBRAIC[133]*CONSTANTS[111])/ALGEBRAIC[99];
ALGEBRAIC[145] = ALGEBRAIC[143]+ALGEBRAIC[144];
ALGEBRAIC[136] = 1.00000/ALGEBRAIC[134];
ALGEBRAIC[146] =  ALGEBRAIC[136]*STATES[2]*CONSTANTS[113];
ALGEBRAIC[149] =  CONSTANTS[199]*ALGEBRAIC[145]*(ALGEBRAIC[147]+ALGEBRAIC[146])+ CONSTANTS[200]*ALGEBRAIC[147]*(CONSTANTS[199]+ALGEBRAIC[142]);
ALGEBRAIC[150] =  CONSTANTS[198]*ALGEBRAIC[147]*(ALGEBRAIC[145]+CONSTANTS[200])+ ALGEBRAIC[145]*ALGEBRAIC[146]*(CONSTANTS[198]+ALGEBRAIC[148]);
ALGEBRAIC[151] =  CONSTANTS[198]*ALGEBRAIC[142]*(ALGEBRAIC[147]+ALGEBRAIC[146])+ ALGEBRAIC[148]*ALGEBRAIC[146]*(CONSTANTS[199]+ALGEBRAIC[142]);
ALGEBRAIC[152] =  CONSTANTS[199]*ALGEBRAIC[148]*(ALGEBRAIC[145]+CONSTANTS[200])+ ALGEBRAIC[142]*CONSTANTS[200]*(CONSTANTS[198]+ALGEBRAIC[148]);
ALGEBRAIC[153] = ALGEBRAIC[149]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[154] = ALGEBRAIC[150]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[155] = ALGEBRAIC[151]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[156] = ALGEBRAIC[152]/(ALGEBRAIC[149]+ALGEBRAIC[150]+ALGEBRAIC[151]+ALGEBRAIC[152]);
ALGEBRAIC[158] = ( 3.00000*( ALGEBRAIC[156]*ALGEBRAIC[147] -  ALGEBRAIC[153]*ALGEBRAIC[148])+ ALGEBRAIC[155]*ALGEBRAIC[144]) -  ALGEBRAIC[154]*ALGEBRAIC[141];
ALGEBRAIC[159] =  ALGEBRAIC[154]*CONSTANTS[199] -  ALGEBRAIC[153]*CONSTANTS[198];
ALGEBRAIC[160] =  0.200000*CONSTANTS[194]*ALGEBRAIC[157]*( CONSTANTS[7]*ALGEBRAIC[158]+ CONSTANTS[8]*ALGEBRAIC[159]);
ALGEBRAIC[190] = ( CONSTANTS[144]*STATES[9])/(CONSTANTS[145]+STATES[9]);
ALGEBRAIC[186] = ( CONSTANTS[143]*4.00000*CONSTANTS[169]*( STATES[9]*exp( 2.00000*ALGEBRAIC[12]) -  0.341000*CONSTANTS[2]))/CONSTANTS[179];
ALGEBRAIC[188] =  CONSTANTS[179]*(STATES[0] - CONSTANTS[165]);
ALGEBRAIC[189] = (- 1.00000e-07<=ALGEBRAIC[188]&&ALGEBRAIC[188]<=1.00000e-07 ?  ALGEBRAIC[186]*(1.00000 -  0.500000*ALGEBRAIC[188]) : ( ALGEBRAIC[186]*ALGEBRAIC[188])/(exp(ALGEBRAIC[188]) - 1.00000));
ALGEBRAIC[191] = (STATES[2] - STATES[9])/0.200000;
ALGEBRAIC[192] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[193] =  CONSTANTS[147]*( (1.00000 - ALGEBRAIC[192])*STATES[47]+ ALGEBRAIC[192]*STATES[48]);
ALGEBRAIC[46] = 1.00000/(1.00000+( CONSTANTS[26]*CONSTANTS[27])/pow(CONSTANTS[27]+STATES[2], 2.00000)+( CONSTANTS[28]*CONSTANTS[29])/pow(CONSTANTS[29]+STATES[2], 2.00000));
ALGEBRAIC[194] = ( CONSTANTS[168]*0.00437500*STATES[9])/(STATES[9]+0.000920000);
ALGEBRAIC[195] = ( CONSTANTS[168]*2.75000*0.00437500*STATES[9])/((STATES[9]+0.000920000) - 0.000170000);
ALGEBRAIC[196] = 1.00000/(1.00000+CONSTANTS[17]/ALGEBRAIC[42]);
ALGEBRAIC[197] = ( 0.00393750*STATES[7])/15.0000;
ALGEBRAIC[198] =  CONSTANTS[148]*(( (1.00000 - ALGEBRAIC[196])*ALGEBRAIC[194]+ ALGEBRAIC[196]*ALGEBRAIC[195]) - ALGEBRAIC[197]);
ALGEBRAIC[44] = 1.00000/(1.00000+( CONSTANTS[150]*CONSTANTS[23])/pow(CONSTANTS[23]+STATES[9], 2.00000)+( CONSTANTS[24]*CONSTANTS[25])/pow(CONSTANTS[25]+STATES[9], 2.00000));
ALGEBRAIC[199] = (STATES[7] - STATES[8])/100.000;
ALGEBRAIC[48] = 1.00000/(1.00000+( CONSTANTS[30]*CONSTANTS[31])/pow(CONSTANTS[31]+STATES[8], 2.00000));
}

