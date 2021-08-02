#include "./koivumaki.h"
#include "math.h"
#include "stdio.h"


void _calc_means(double *STATES, const double *CONSTANTS)
{
    double CaSR = 0, VSR_sum = 0;
    double Cai = 0, V_sum = 0;
    double fluo = 0;

    for (int i = 0; i < 4; ++i) {

        double VSRi = CONSTANTS[13 + i];
        double CaSRi = STATES[0 + i];
        CaSR += CaSRi * VSRi;
        VSR_sum += VSRi;

        double Vnonjuncti = CONSTANTS[7 + i];
        double Caii = STATES[4 + i];
        double fluo_i = STATES[43 + i];
        Cai += Caii * Vnonjuncti;
        fluo += fluo_i * Vnonjuncti;
        V_sum += Vnonjuncti;

    }

    double Cass = STATES[8];
    double fluo_ss = STATES[47];
    double Vss = CONSTANTS[1];
    Cai += Cass * Vss;
    fluo += fluo_ss * Vss;
    V_sum += Vss;

    CaSR /= VSR_sum;
    Cai /= V_sum;
    fluo /= V_sum;

    STATES[48] = CaSR;
    STATES[49] = Cai;
    STATES[50] = fluo;
}


void initialize_states_default(double *STATES, const double *CONSTANTS) {
    STATES[0]	=	0.6189225;	//  CaSR1
    STATES[1]	=	0.6076289;	//  CaSR2
    STATES[2]	=	0.5905266;	//  CaSR3
    STATES[3]	=	0.5738108;	//  CaSR4
    STATES[4]	=	1.35496500000000013e-04;	//  Cai1
    STATES[5]	=	1.38142100000000014e-04;	//  Cai2
    STATES[6]	=	1.44208699999999994e-04;	//  Cai3
    STATES[7]	=	1.56184399999999995e-04;	//  Cai4
    STATES[8]	=	1.61937700000000013e-04;	//  Cass
    STATES[9]	=	1.06091699999999996e-05;	//  d
    STATES[10]	=	0.9988566;	//  f1
    STATES[11]	=	0.9988624;	//  f2
    STATES[12]	=	0.9744374;	//  fca
    STATES[13]	=	5.62066499999999969e-02;	//  y
    STATES[14]	=	4.18941700000000008e-05;	//  pa
    STATES[15]	=	4.10975100000000003e-03;	//  n
    STATES[16]	=	3.11170299999999984e-04;	//  ikur_r
    STATES[17]	=	0.9751094;	//  ikur_s
    STATES[18]	=	0.90391;	//  h1
    STATES[19]	=	0.9039673;	//  h2
    STATES[20]	=	2.77581199999999990e-03;	//  m
    STATES[21]	=	9.59425800000000026e-04;	//  it_r
    STATES[22]	=	0.954338;	//  it_s
    STATES[23]	=	-75.42786;	//  V
    STATES[24]	=	134.6313;	//  Ki
    STATES[25]	=	0.1925362;	//  ryr_a1
    STATES[26]	=	0.2010345;	//  ryr_a2
    STATES[27]	=	0.2163122;	//  ryr_a3
    STATES[28]	=	0.2455297;	//  ryr_ass
    STATES[29]	=	0.9993722;	//  c1
    STATES[30]	=	0.9995086;	//  c2
    STATES[31]	=	0.9995604;	//  c3
    STATES[32]	=	0.9999717;	//  css
    STATES[33]	=	9.47851400000000044e-05;	//  o1
    STATES[34]	=	7.76550300000000031e-05;	//  o2
    STATES[35]	=	5.67494700000000006e-05;	//  o3
    STATES[36]	=	3.97509699999999973e-05;	//  oss
    STATES[37]	=	4.63856499999999988e-03 / CONSTANTS[53];	//  serca_a1
    STATES[38]	=	4.51207800000000010e-03 / CONSTANTS[53];	//  serca_a2
    STATES[39]	=	4.32640899999999981e-03 / CONSTANTS[53];	//  serca_a3
    STATES[40]	=	4.25044500000000026e-03 / CONSTANTS[53];	//  serca_ass
    STATES[41]	=	9.28686;	//  Nai
    STATES[42]	=	8.691504;	//  Nass

    // Fluo-3
    const double fluo_tot = CONSTANTS[89];
    const double k_on = CONSTANTS[90], k_off = CONSTANTS[91];

    for (int i = 0; i < 5; ++i) {

      int i_fluo = 43 + i;
      int i_Cai  = 4 + i;

      double Cai  = STATES[i_Cai];

      STATES[i_fluo] = Cai * fluo_tot / (k_off / k_on + Cai);

    }

    _calc_means(STATES, CONSTANTS);

}


void initialize_constants_default(double *CONSTANTS) {
    CONSTANTS[0]	=	0.05;	//  Cm
    CONSTANTS[1]	=	4.99231999999999966e-05;	//  Vss
    CONSTANTS[2]	=	1.625;	//  dx
    CONSTANTS[3]	=	122.051;	//  lcell
    CONSTANTS[4]	=	3.14159265358979312e+00;	//  cell_pi
    CONSTANTS[5]	=	6.5;	//  rjunct
    CONSTANTS[6]	=	(CONSTANTS[4] * CONSTANTS[5] * 2.0 * CONSTANTS[3] * 0.5);	//  Aj_nj
    CONSTANTS[7]	=	((pow(1.0 * CONSTANTS[2], 2.0) - pow(0.0 * CONSTANTS[2], 2.0)) * CONSTANTS[4] * CONSTANTS[3] * 0.5 * 1e-06);	//  Vnonjunct1
    CONSTANTS[8]	=	((pow(2.0 * CONSTANTS[2], 2.0) - pow(1.0 * CONSTANTS[2], 2.0)) * CONSTANTS[4] * CONSTANTS[3] * 0.5 * 1e-06);	//  Vnonjunct2
    CONSTANTS[9]	=	((pow(3.0 * CONSTANTS[2], 2.0) - pow(2.0 * CONSTANTS[2], 2.0)) * CONSTANTS[4] * CONSTANTS[3] * 0.5 * 1e-06);	//  Vnonjunct3
    CONSTANTS[10]	=	((pow(4.0 * CONSTANTS[2], 2.0) - pow(3.0 * CONSTANTS[2], 2.0)) * CONSTANTS[4] * CONSTANTS[3] * 0.5 * 1e-06);	//  Vnonjunct4
    CONSTANTS[11]	=	(0.02 / 2.0 + CONSTANTS[2] / 2.0);	//  xj_nj
    CONSTANTS[12]	=	(0.02 / 2.0 + 2.0 * CONSTANTS[2]);	//  xj_nj_Nai
    CONSTANTS[13]	=	(0.05 * CONSTANTS[7] / 2.0 * 0.9);	//  VSR1
    CONSTANTS[14]	=	(0.05 * CONSTANTS[8] / 2.0 * 0.9);	//  VSR2
    CONSTANTS[15]	=	(0.05 * CONSTANTS[9] / 2.0 * 0.9);	//  VSR3
    CONSTANTS[16]	=	(0.05 * CONSTANTS[10] / 2.0 * 0.9);	//  VSR4
    CONSTANTS[17]	=	(CONSTANTS[7] + CONSTANTS[8] + CONSTANTS[9] + CONSTANTS[10]);	//  Vnonjunct_Nai
    CONSTANTS[18]	=	(CONSTANTS[17] + CONSTANTS[1]);	//  Vcytosol
    CONSTANTS[19]	=	1.8;	//  Cao
    CONSTANTS[20]	=	5.4;	//  Ko
    CONSTANTS[21]	=	130.0;	//  Nao
    CONSTANTS[22]	=	60.0;	//  ECa_app
    CONSTANTS[23]	=	25.3125;	//  gCaL
    CONSTANTS[24]	=	0.002;	//  ical_fca_tau
    CONSTANTS[25]	=	0.001;	//  kCa
    CONSTANTS[26]	=	2.0;	//  kCan
    CONSTANTS[27]	=	2.0;	//  ICaPmax
    CONSTANTS[28]	=	0.0005;	//  kCaP
    CONSTANTS[29]	=	96487.0;	//  F
    CONSTANTS[30]	=	8314.0;	//  R
    CONSTANTS[31]	=	306.15;	//  T
    CONSTANTS[32]	=	(CONSTANTS[30] * CONSTANTS[31] / CONSTANTS[29]);	//  RTF
    CONSTANTS[33]	=	(1.0 / CONSTANTS[32]);	//  1 / RTF
    CONSTANTS[34]	=	((-2500.0));	//  amplitude
    CONSTANTS[35]	=	0.0003;	//  dNaCa
    CONSTANTS[36]	=	1.0;	//  fCaNCX
    CONSTANTS[37]	=	0.45;	//  gam
    CONSTANTS[38]	=	0.0084;	//  kNaCa
    CONSTANTS[39]	=	70.8253;	//  INaKmax
    CONSTANTS[40]	=	1.0;	//  kNaKK
    CONSTANTS[41]	=	11.0;	//  kNaKNa
    CONSTANTS[42]	=	(1.0 * CONSTANTS[7]);	//  nu1
    CONSTANTS[43]	=	(1.0 * CONSTANTS[8]);	//  nu2
    CONSTANTS[44]	=	(1.0 * CONSTANTS[9]);	//  nu3
    CONSTANTS[45]	=	(625.0 * CONSTANTS[1]);	//  nuss
    CONSTANTS[46]	=	0.01875;	//  tau_act
    CONSTANTS[47]	=	0.005;	//  tau_actss
    CONSTANTS[48]	=	1.0;	//  tau_adapt
    CONSTANTS[49]	=	0.0875;	//  tau_inact
    CONSTANTS[50]	=	0.015;	//  tau_inactss
    CONSTANTS[51]	=	0.00025;	//  SERCAKmf
    CONSTANTS[52]	=	1.8;	//  SERCAKmr
    CONSTANTS[53]	=	0.04;	//  cpumps
    CONSTANTS[54]	=	7.5;	//  k4
    CONSTANTS[55]	=	(pow(1000.0, 2.0) * CONSTANTS[54]);	//  k1
    CONSTANTS[56]	=	(CONSTANTS[54] / pow(CONSTANTS[52], 2.0));	//  k3
    CONSTANTS[57]	=	(CONSTANTS[55] * pow(CONSTANTS[51], 2.0));	//  k2
    CONSTANTS[58]	=	0.0952;	//  gCab
    CONSTANTS[59]	=	1.0;	//  gIf
    CONSTANTS[60]	=	(3.825 * 0.9);	//  gK1
    CONSTANTS[61]	=	0.5;	//  gKr
    CONSTANTS[62]	=	1.0;	//  gKs
    CONSTANTS[63]	=	(0.89 * 2.75);	//  gKur
    CONSTANTS[64]	=	0.0018;	//  PNa
    CONSTANTS[65]	=	0.060599;	//  gNab
    CONSTANTS[66]	=	(1.09 * 7.5);	//  gt
    CONSTANTS[67]	=	0.024;	//  BCa
    CONSTANTS[68]	=	6.7;	//  CSQN
    CONSTANTS[69]	=	780.0;	//  DCa
    CONSTANTS[70]	=	25.0;	//  DCaBm
    CONSTANTS[71]	=	44.0;	//  DCaSR
    CONSTANTS[72]	=	0.00238;	//  KdBCa
    CONSTANTS[73]	=	0.8;	//  KdCSQN
    CONSTANTS[74]	=	0.013;	//  KdSLhigh
    CONSTANTS[75]	=	1.1;	//  KdSLlow
    CONSTANTS[76]	=	13.0;	//  SLhigh
    CONSTANTS[77]	=	165.0;	//  SLlow
    CONSTANTS[78]	=	0.006;	//  kSRleak
    CONSTANTS[79]	=	(0.49 * 2.31);	//  BNa
    CONSTANTS[80]	=	0.12;	//  DNa
    CONSTANTS[81]	=	10.0;	//  KdBNa

    CONSTANTS[82]   =   1.0;    //  STIM_LEVEL
    CONSTANTS[83]   =   0.001;  //  STIM_DURATION
    CONSTANTS[84]   =   0.0;    //  STIM_OFFSET
    CONSTANTS[85]   =   1.0;    //  STIM_PERIOD

    CONSTANTS[86]   =   1.0;    //  Jrel_multiplier
    CONSTANTS[87]   =   1.0;    //  J_SERCASR_multiplier
    CONSTANTS[88]   =   1.0;    //  J_bulkSERCA_multiplier

    CONSTANTS[89]   =   0.0;    //  fluo_tot
    CONSTANTS[90]   =   236000; //  fluo_k_on
    CONSTANTS[91]   =   175;    //  fluo_k_off

    CONSTANTS[92]   =   4.00E-15;   //  G_seal
    CONSTANTS[93]   =   0.01;      // gKb

    CONSTANTS[94]   = 4.75462;  // K_mobility
    CONSTANTS[95]   = 3.24155;  // Na_mobility

    CONSTANTS[96]   = 4e-15;    // K_ghk_scaler_seal
    CONSTANTS[97]   = 4e-15;    // K_ghk_scaler_membrane
    CONSTANTS[98]   = 3.33e-15; // Na_ghk_scaler_seal
    CONSTANTS[99]   = 3.33e-15; // Na_ghk_scaler_membrane

    CONSTANTS[100]  = 48;       // Cl_i, mM
    CONSTANTS[101]  = 150;      // Cl_o, mM
    CONSTANTS[102]  = 4.93639;  // Cl_mobility

    CONSTANTS[103]  = 92;       // Aspartate_i, mM
    CONSTANTS[104]  = 0;        // Aspartate_o, mM
    CONSTANTS[105]  = 1.42639;  // Aspartate_mobility

    CONSTANTS[106]  = 1.92499;  // Ca_mobility
    CONSTANTS[107]  = 1e-15;  // Ca_ghk_scaler_membrane

}


double calc_ghk(int z, double u, double Ci, double Co, double *STATES, double *CONSTANTS) {

    double F = CONSTANTS[29];
    double R = CONSTANTS[30];
    double T = CONSTANTS[31];
    // double F_RT = CONSTANTS[33];

    double V = STATES[23];
    double zFVRT = z * V * F / R / T;

    double scaler = 1; // gamma / h;
    double P = u * R * T * scaler;

    double J = P * zFVRT * (Co - Ci * exp(zFVRT)) / (1 - exp(zFVRT));
    double I = z * F * J;
    return I;
}


void calc_ghk_K(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES) {
    int z_K = +1;
    double u_K = CONSTANTS[94];
    double Ci_K = STATES[24], Co_K = CONSTANTS[20];
    double ghk_K = calc_ghk(z_K, u_K, Ci_K, Co_K, STATES, CONSTANTS);
    ALGEBRAIC[109] = ghk_K;
}


void calc_ghk_Na(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES) {
    int z_Na = +1;
    double u_Na = CONSTANTS[95];
    double Ci_Na = STATES[41], Co_Na = CONSTANTS[21];
    double ghk_Na = calc_ghk(z_Na, u_Na, Ci_Na, Co_Na, STATES, CONSTANTS);
    ALGEBRAIC[110] = ghk_Na;
}


void calc_ghk_Cl(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES) {
    int z_Cl = -1;
    double u_Cl = CONSTANTS[102];
    double Ci_Cl = CONSTANTS[100], Co_Cl = CONSTANTS[101];
    double ghk_Cl = calc_ghk(z_Cl, u_Cl, Ci_Cl, Co_Cl, STATES, CONSTANTS);
    ALGEBRAIC[111] = ghk_Cl;
}


void calc_ghk_Aspartate(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES) {
    int z_Aspartate = -1;
    double u_Aspartate = CONSTANTS[105];
    double Ci_Aspartate = CONSTANTS[103], Co_Aspartate = CONSTANTS[104];
    double ghk_Aspartate = calc_ghk(z_Aspartate, u_Aspartate, Ci_Aspartate, Co_Aspartate, STATES, CONSTANTS);
    ALGEBRAIC[112] = ghk_Aspartate;
}

void calc_ghk_Ca(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES) {
    int z_Ca = +2;
    double u_Ca = CONSTANTS[106];
    double Ci_Ca = STATES[8], Co_Ca = CONSTANTS[19];
    double ghk_Ca = calc_ghk(z_Ca, u_Ca, Ci_Ca, Co_Ca, STATES, CONSTANTS);
    ALGEBRAIC[113] = ghk_Ca;
}


/* Components */
void calc_calcium(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[2] = (-ALGEBRAIC[84]) /*- ALGEBRAIC[86] */ -  - ALGEBRAIC[95] + 2.0 * ALGEBRAIC[62];
    ALGEBRAIC[2] -= ALGEBRAIC[86] + CONSTANTS[92] * ALGEBRAIC[113];  // ICab + G_seal * ICa_ghk

    ALGEBRAIC[24] = CONSTANTS[78] * (STATES[0] - STATES[4]) * CONSTANTS[7];
    ALGEBRAIC[23] = CONSTANTS[78] * (STATES[1] - STATES[5]) * CONSTANTS[8];
    ALGEBRAIC[36] = CONSTANTS[78] * (STATES[2] - STATES[6]) * CONSTANTS[9];
    ALGEBRAIC[21] = CONSTANTS[78] * (STATES[3] - STATES[8]) * CONSTANTS[1];
    ALGEBRAIC[82] = CONSTANTS[69] * CONSTANTS[6] / CONSTANTS[11] * (STATES[8] - STATES[7]) * 1e-06;
    ALGEBRAIC[3] = 1.0 / (1.0 + CONSTANTS[68] * CONSTANTS[73] / pow(STATES[0] + CONSTANTS[73], 2.0));
    ALGEBRAIC[7] = 1.0 / (1.0 + CONSTANTS[68] * CONSTANTS[73] / pow(STATES[1] + CONSTANTS[73], 2.0));
    ALGEBRAIC[5] = 1.0 / (1.0 + CONSTANTS[68] * CONSTANTS[73] / pow(STATES[2] + CONSTANTS[73], 2.0));
    ALGEBRAIC[0] = 1.0 / (1.0 + CONSTANTS[68] * CONSTANTS[73] / pow(STATES[3] + CONSTANTS[73], 2.0));
    ALGEBRAIC[6] = CONSTANTS[67] * CONSTANTS[72] / pow(STATES[4] + CONSTANTS[72], 2.0);
    ALGEBRAIC[4] = CONSTANTS[67] * CONSTANTS[72] / pow(STATES[5] + CONSTANTS[72], 2.0);
    ALGEBRAIC[1] = CONSTANTS[67] * CONSTANTS[72] / pow(STATES[6] + CONSTANTS[72], 2.0);
    ALGEBRAIC[8] = CONSTANTS[67] * CONSTANTS[72] / pow(STATES[7] + CONSTANTS[72], 2.0);
    ALGEBRAIC[13] = 1.0 / (1.0 + CONSTANTS[77] * CONSTANTS[75] / pow(STATES[8] + CONSTANTS[75], 2.0) + CONSTANTS[76] * CONSTANTS[74] / pow(STATES[8] + CONSTANTS[74], 2.0) + CONSTANTS[67] * CONSTANTS[72] / pow(STATES[8] + CONSTANTS[72], 2.0));
    ALGEBRAIC[94] = (-ALGEBRAIC[17]) + ALGEBRAIC[24] + ALGEBRAIC[66];
    ALGEBRAIC[89] = (-ALGEBRAIC[15]) + ALGEBRAIC[23] + ALGEBRAIC[64];
    ALGEBRAIC[91] = (-ALGEBRAIC[16]) + ALGEBRAIC[36] + ALGEBRAIC[65];
    ALGEBRAIC[85] = ALGEBRAIC[82];
    ALGEBRAIC[63] = (-ALGEBRAIC[82]) + ALGEBRAIC[21] - ALGEBRAIC[14] + ALGEBRAIC[58];
    ALGEBRAIC[54] = ALGEBRAIC[30] - ALGEBRAIC[24] - ALGEBRAIC[66];
    ALGEBRAIC[59] = ALGEBRAIC[33] - ALGEBRAIC[23] - ALGEBRAIC[64];
    ALGEBRAIC[50] = ALGEBRAIC[34] - ALGEBRAIC[36] - ALGEBRAIC[65];
    ALGEBRAIC[53] = ALGEBRAIC[19] - ALGEBRAIC[21] - ALGEBRAIC[58];
    ALGEBRAIC[11] = 1.0 / (1.0 + ALGEBRAIC[6]);
    ALGEBRAIC[9] = 1.0 / (1.0 + ALGEBRAIC[4]);
    ALGEBRAIC[10] = 1.0 / (1.0 + ALGEBRAIC[1]);
    ALGEBRAIC[12] = 1.0 / (1.0 + ALGEBRAIC[8]);
    RATES[0] = ALGEBRAIC[3] * CONSTANTS[71] * ((STATES[1] - 2.0 * STATES[0] + STATES[0]) / pow(CONSTANTS[2], 2.0) + (STATES[1] - STATES[0]) / (2.0 * 1.0 * pow(CONSTANTS[2], 2.0))) + ALGEBRAIC[54] / CONSTANTS[13] * ALGEBRAIC[3];
    RATES[1] = ALGEBRAIC[7] * CONSTANTS[71] * ((STATES[2] - 2.0 * STATES[1] + STATES[0]) / pow(CONSTANTS[2], 2.0) + (STATES[2] - STATES[0]) / (2.0 * 2.0 * pow(CONSTANTS[2], 2.0))) + ALGEBRAIC[59] / CONSTANTS[14] * ALGEBRAIC[7];
    RATES[2] = ALGEBRAIC[5] * CONSTANTS[71] * ((STATES[3] - 2.0 * STATES[2] + STATES[1]) / pow(CONSTANTS[2], 2.0) + (STATES[3] - STATES[1]) / (2.0 * 3.0 * pow(CONSTANTS[2], 2.0))) + ALGEBRAIC[50] / CONSTANTS[15] * ALGEBRAIC[5];
    RATES[3] = ALGEBRAIC[0] * CONSTANTS[71] * ((STATES[3] - 2.0 * STATES[3] + STATES[2]) / pow(CONSTANTS[2], 2.0) + (STATES[3] - STATES[2]) / (2.0 * 4.0 * pow(CONSTANTS[2], 2.0))) + ALGEBRAIC[53] / CONSTANTS[16] * ALGEBRAIC[0];
    RATES[4] = ALGEBRAIC[11] * (CONSTANTS[69] + ALGEBRAIC[6] * CONSTANTS[70]) * ((STATES[5] - 2.0 * STATES[4] + STATES[4]) / pow(CONSTANTS[2], 2.0) + (STATES[5] - STATES[4]) / (2.0 * 1.0 * pow(CONSTANTS[2], 2.0))) - 2.0 * ALGEBRAIC[11] * ALGEBRAIC[6] * CONSTANTS[70] / (CONSTANTS[72] + STATES[4]) * pow((STATES[5] - STATES[4]) / (2.0 * CONSTANTS[2]), 2.0) + ALGEBRAIC[94] / CONSTANTS[7] * ALGEBRAIC[11];
    RATES[5] = ALGEBRAIC[9] * (CONSTANTS[69] + ALGEBRAIC[4] * CONSTANTS[70]) * ((STATES[6] - 2.0 * STATES[5] + STATES[4]) / pow(CONSTANTS[2], 2.0) + (STATES[6] - STATES[4]) / (2.0 * 2.0 * pow(CONSTANTS[2], 2.0))) - 2.0 * ALGEBRAIC[9] * ALGEBRAIC[4] * CONSTANTS[70] / (CONSTANTS[72] + STATES[5]) * pow((STATES[6] - STATES[4]) / (2.0 * CONSTANTS[2]), 2.0) + ALGEBRAIC[89] / CONSTANTS[8] * ALGEBRAIC[9];
    RATES[6] = ALGEBRAIC[10] * (CONSTANTS[69] + ALGEBRAIC[1] * CONSTANTS[70]) * ((STATES[7] - 2.0 * STATES[6] + STATES[5]) / pow(CONSTANTS[2], 2.0) + (STATES[7] - STATES[5]) / (2.0 * 3.0 * pow(CONSTANTS[2], 2.0))) - 2.0 * ALGEBRAIC[10] * ALGEBRAIC[1] * CONSTANTS[70] / (CONSTANTS[72] + STATES[6]) * pow((STATES[7] - STATES[5]) / (2.0 * CONSTANTS[2]), 2.0) + ALGEBRAIC[91] / CONSTANTS[9] * ALGEBRAIC[10];
    RATES[7] = ALGEBRAIC[12] * (CONSTANTS[69] + ALGEBRAIC[8] * CONSTANTS[70]) * ((STATES[7] - 2.0 * STATES[7] + STATES[6]) / pow(CONSTANTS[2], 2.0) + (STATES[7] - STATES[6]) / (2.0 * 4.0 * pow(CONSTANTS[2], 2.0))) - 2.0 * ALGEBRAIC[12] * ALGEBRAIC[8] * CONSTANTS[70] / (CONSTANTS[72] + STATES[7]) * pow((STATES[7] - STATES[6]) / (2.0 * CONSTANTS[2]), 2.0) + ALGEBRAIC[85] / CONSTANTS[10] * ALGEBRAIC[12];
    RATES[8] = ALGEBRAIC[13] * (ALGEBRAIC[63] / CONSTANTS[1] + ALGEBRAIC[2] / (2.0 * CONSTANTS[1] * CONSTANTS[29]));
}

void calc_icab(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[86] = CONSTANTS[58] * (STATES[23] - ALGEBRAIC[98]);
    // double ghk_Ca = ALGEBRAIC[113];
    // double Ca_ghk_scaler_membrane = CONSTANTS[107];
    // ALGEBRAIC[86] = Ca_ghk_scaler_membrane * ghk_Ca;
}

void calc_ical(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[80] = 1.0 / (1.0 + exp((STATES[23] + 27.4) / 7.1));
    ALGEBRAIC[37] = 1.0 / (1.0 + exp((STATES[23] + 9.0) / (-5.8)));
    ALGEBRAIC[29] = 0.0027 * exp((-pow((STATES[23] + 35.0) / 30.0, 2.0))) + 0.002;
    ALGEBRAIC[22] = 0.98698 * exp((-pow((STATES[23] + 30.16047) / 7.09396, 2.0))) + 0.04275 / (1.0 + exp((STATES[23] - 51.61555) / (-80.61331))) + 0.03576 / (1.0 + exp((STATES[23] + 29.57272) / 13.21758)) - 0.00821;
    ALGEBRAIC[20] = 1.3323 * exp((-pow((STATES[23] + 40.0) / 14.2, 2.0))) + 0.0626;
    RATES[9] = (ALGEBRAIC[37] - STATES[9]) / ALGEBRAIC[29];
    RATES[10] = (ALGEBRAIC[80] - STATES[10]) / ALGEBRAIC[22];
    RATES[11] = (ALGEBRAIC[80] - STATES[11]) / ALGEBRAIC[20];

    //  TODO: make ICaL ( A[41] ) modification clear
    double E_Ca = CONSTANTS[32] / 2 * log(CONSTANTS[19] / STATES[8]);  //  RT/F * log(Ca_c / Ca_d)
    // ALGEBRAIC[84] = CONSTANTS[23] * STATES[9] * STATES[12] * STATES[10] * STATES[11] * (STATES[23] - CONSTANTS[22]); // ICaL
    ALGEBRAIC[84] = CONSTANTS[23] * STATES[9] * STATES[12] * STATES[10] * STATES[11] * (STATES[23] - E_Ca); // ICaL
    ALGEBRAIC[18] = 1.0 - 1.0 / (1.0 + pow(CONSTANTS[25] / STATES[8], CONSTANTS[26]));
    RATES[12] = (ALGEBRAIC[18] - STATES[12]) / CONSTANTS[24];
}

void calc_icap(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[95] = CONSTANTS[27] * STATES[8] / (CONSTANTS[28] + STATES[8]);
}

void calc_if(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[45] = 1.0 / (1.0 + exp((STATES[23] + 97.82874) / 12.48025));
    ALGEBRAIC[43] = 1.0 / (0.00332 * exp((-STATES[23]) / 16.54103) + 23.71839 * exp(STATES[23] / 16.54103));
    RATES[13] = (ALGEBRAIC[45] - STATES[13]) / ALGEBRAIC[43];
    ALGEBRAIC[102] = CONSTANTS[59] * STATES[13] * ((1.0 - 0.2677) * (STATES[23] - ALGEBRAIC[105]));
    ALGEBRAIC[90] = CONSTANTS[59] * STATES[13] * (0.2677 * (STATES[23] - ALGEBRAIC[97]));
    ALGEBRAIC[106] = ALGEBRAIC[102] + ALGEBRAIC[90];
}

void calc_ik1(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[103] = CONSTANTS[60] * pow(CONSTANTS[20] * 1.0, 0.4457) * (STATES[23] - ALGEBRAIC[105]) / (1.0 + exp(1.5 * (STATES[23] - ALGEBRAIC[105] + 3.6) * CONSTANTS[33]));
}

void calc_ikr(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[32] = 1.0 / (1.0 + exp((STATES[23] + 15.0) / (-6.0)));
    ALGEBRAIC[35] = 0.21718 * exp((-pow((STATES[23] + 20.1376) / 22.1996, 2.0))) + 0.03118;
    ALGEBRAIC[52] = 1.0 / (1.0 + exp((STATES[23] + 55.0) / 24.0));
    RATES[14] = (ALGEBRAIC[32] - STATES[14]) / ALGEBRAIC[35];
    ALGEBRAIC[101] = CONSTANTS[61] * STATES[14] * ALGEBRAIC[52] * (STATES[23] - ALGEBRAIC[105]);
}

void calc_iks(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[41] = 1.0 / (1.0 + exp((STATES[23] - 19.9) / (-12.7)));
    ALGEBRAIC[39] = 0.4 * exp((-pow((STATES[23] - 20.0) / 20.0, 2.0))) + 0.7;
    RATES[15] = (ALGEBRAIC[41] - STATES[15]) / ALGEBRAIC[39];
    ALGEBRAIC[100] = CONSTANTS[62] * STATES[15] * (STATES[23] - ALGEBRAIC[105]);
}

void calc_ikur(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[25] = 1.0 / (1.0 + exp((STATES[23] + 6.0) / (-8.6)));
    ALGEBRAIC[28] = 0.009 / (1.0 + exp((STATES[23] + 5.0) / 12.0)) + 0.0005;
    ALGEBRAIC[27] = 1.0 / (1.0 + exp((STATES[23] + 7.5) / 10.0));
    ALGEBRAIC[31] = 0.59 / (1.0 + exp((STATES[23] + 60.0) / 10.0)) + 3.05;
    RATES[16] = (ALGEBRAIC[25] - STATES[16]) / ALGEBRAIC[28];
    RATES[17] = (ALGEBRAIC[27] - STATES[17]) / ALGEBRAIC[31];
    ALGEBRAIC[87] = CONSTANTS[63] * STATES[16] * STATES[17] * (STATES[23] - ALGEBRAIC[105]);
}

void calc_ina(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[76] = 1.0 / (1.0 + exp((STATES[23] + 63.6) / 5.3));
    ALGEBRAIC[26] = 0.03 / (1.0 + exp((STATES[23] + 35.1) / 3.2)) + 0.0003;
    ALGEBRAIC[38] = 0.12 / (1.0 + exp((STATES[23] + 35.1) / 3.2)) + 0.003;
    ALGEBRAIC[42] = 1.0 / (1.0 + exp((STATES[23] + 27.12) / (-8.21)));
    ALGEBRAIC[40] = 4.2e-05 * exp((-pow((STATES[23] + 25.57) / 28.8, 2.0))) + 2.4e-05;
    RATES[18] = (ALGEBRAIC[76] - STATES[18]) / ALGEBRAIC[26];
    RATES[19] = (ALGEBRAIC[76] - STATES[19]) / ALGEBRAIC[38];
    RATES[20] = (ALGEBRAIC[42] - STATES[20]) / ALGEBRAIC[40];
    ALGEBRAIC[99] = CONSTANTS[64] * pow(STATES[20], 3.0) * (0.9 * STATES[18] + 0.1 * STATES[19]) * CONSTANTS[21] * STATES[23] * CONSTANTS[29] * CONSTANTS[33] * (exp((STATES[23] - ALGEBRAIC[97]) * CONSTANTS[33]) - 1.0) / (exp(STATES[23] * CONSTANTS[33]) - 1.0);
}

void calc_inab(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    // ALGEBRAIC[92] = CONSTANTS[65] * (STATES[23] - ALGEBRAIC[97]);
    double ghk_Na = ALGEBRAIC[110];
    double Na_ghk_scaler_membrane = CONSTANTS[99];
    ALGEBRAIC[92] = Na_ghk_scaler_membrane * ghk_Na;
}

void calc_inaca(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[62] = CONSTANTS[38] * (exp(CONSTANTS[37] * STATES[23] * CONSTANTS[33]) * pow(STATES[42], 3.0) * CONSTANTS[19] - exp((CONSTANTS[37] - 1.0) * STATES[23] * CONSTANTS[33]) * pow(CONSTANTS[21], 3.0) * STATES[8] * CONSTANTS[36]) / (1.0 + CONSTANTS[35] * (pow(CONSTANTS[21], 3.0) * STATES[8] * CONSTANTS[36] + pow(STATES[42], 3.0) * CONSTANTS[19]));
}

void calc_inak(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[55] = pow(STATES[42] * 1.0, 1.5);
    ALGEBRAIC[88] = CONSTANTS[39] * CONSTANTS[20] / (CONSTANTS[20] + CONSTANTS[40]) * ALGEBRAIC[55] / (ALGEBRAIC[55] + pow(CONSTANTS[41], 1.5)) * (STATES[23] + 150.0) / (STATES[23] + 200.0);
}

void calc_it(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[46] = 1.0 / (1.0 + exp((STATES[23] - 1.0) / (-11.0)));
    ALGEBRAIC[48] = 0.0035 * exp((-pow((STATES[23] + 0.0) / 30.0, 2.0))) + 0.0015;
    ALGEBRAIC[44] = 1.0 / (1.0 + exp((STATES[23] + 40.5) / 11.5));
    ALGEBRAIC[49] = 0.025635 * exp((-pow((STATES[23] + 52.45) / 15.8827, 2.0))) + 0.01414;
    RATES[21] = (ALGEBRAIC[46] - STATES[21]) / ALGEBRAIC[48];
    RATES[22] = (ALGEBRAIC[44] - STATES[22]) / ALGEBRAIC[49];
    ALGEBRAIC[104] = CONSTANTS[66] * STATES[21] * STATES[22] * (STATES[23] - ALGEBRAIC[105]);
}

void calc_membrane(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[73] = ALGEBRAIC[99] + ALGEBRAIC[84] + ALGEBRAIC[104] + ALGEBRAIC[87] + ALGEBRAIC[103] + ALGEBRAIC[101] + ALGEBRAIC[100] + ALGEBRAIC[92] + ALGEBRAIC[86] + ALGEBRAIC[88] + ALGEBRAIC[95] + ALGEBRAIC[62] + ALGEBRAIC[106];
    ALGEBRAIC[73] += ALGEBRAIC[107]; // I_Iseal
    RATES[23] = (-(ALGEBRAIC[73] + ALGEBRAIC[56])) / CONSTANTS[0];
}

void calc_nernst(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[98] = CONSTANTS[32] * log(CONSTANTS[19] / STATES[8]) / 2.0;
    ALGEBRAIC[105] = CONSTANTS[32] * log(CONSTANTS[20] / STATES[24]);
    ALGEBRAIC[97] = CONSTANTS[32] * log(CONSTANTS[21] / STATES[42]);
}

void calc_potassium(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[71] = ALGEBRAIC[104] + ALGEBRAIC[87] + ALGEBRAIC[103] + ALGEBRAIC[101] + ALGEBRAIC[100] - 2.0 * ALGEBRAIC[88] + ALGEBRAIC[102] + ALGEBRAIC[56];
    ALGEBRAIC[71] += ALGEBRAIC[108] + CONSTANTS[92] * ALGEBRAIC[109];  // IKb + G_seal * IK_ghk
    RATES[24] = (-ALGEBRAIC[71]) / (CONSTANTS[18] * CONSTANTS[29]);
}

void calc_ryr(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[72] = 1.0 - 1.0 / (1.0 + exp((STATES[0] - 0.3) / 0.1));
    ALGEBRAIC[83] = 1.0 - 1.0 / (1.0 + exp((STATES[1] - 0.3) / 0.1));
    ALGEBRAIC[77] = 1.0 - 1.0 / (1.0 + exp((STATES[2] - 0.3) / 0.1));
    ALGEBRAIC[61] = 1.0 - 1.0 / (1.0 + exp((STATES[3] - 0.3) / 0.1));
    ALGEBRAIC[74] = 0.505 - 0.427 / (1.0 + exp((STATES[4] * 1000.0 - 0.29) / 0.082));
    ALGEBRAIC[69] = 0.505 - 0.427 / (1.0 + exp((STATES[5] * 1000.0 - 0.29) / 0.082));
    ALGEBRAIC[67] = 0.505 - 0.427 / (1.0 + exp((STATES[6] * 1000.0 - 0.29) / 0.082));
    ALGEBRAIC[60] = 0.505 - 0.427 / (1.0 + exp((STATES[8] * 1000.0 - 0.29) / 0.082));
    ALGEBRAIC[81] = 1.0 / (1.0 + exp((STATES[4] * 1000.0 - (STATES[25] + 0.02)) / 0.01));
    ALGEBRAIC[75] = 1.0 / (1.0 + exp((STATES[5] * 1000.0 - (STATES[26] + 0.02)) / 0.01));
    ALGEBRAIC[70] = 1.0 / (1.0 + exp((STATES[6] * 1000.0 - (STATES[27] + 0.02)) / 0.01));
    ALGEBRAIC[57] = 1.0 / (1.0 + exp((STATES[8] * 1000.0 - (STATES[28] + 0.02)) / 0.01));
    ALGEBRAIC[79] = 1.0 - 1.0 / (1.0 + exp((STATES[4] * 1000.0 - (STATES[25] + 0.22)) / 0.03));
    ALGEBRAIC[78] = 1.0 - 1.0 / (1.0 + exp((STATES[5] * 1000.0 - (STATES[26] + 0.22)) / 0.03));
    ALGEBRAIC[68] = 1.0 - 1.0 / (1.0 + exp((STATES[6] * 1000.0 - (STATES[27] + 0.22)) / 0.03));
    ALGEBRAIC[51] = 1.0 - 1.0 / (1.0 + exp((STATES[8] * 1000.0 - (STATES[28] + 0.22)) / 0.03));
    RATES[25] = (ALGEBRAIC[74] - STATES[25]) / CONSTANTS[48];
    RATES[26] = (ALGEBRAIC[69] - STATES[26]) / CONSTANTS[48];
    RATES[27] = (ALGEBRAIC[67] - STATES[27]) / CONSTANTS[48];
    RATES[28] = (ALGEBRAIC[60] - STATES[28]) / CONSTANTS[48];
    RATES[29] = (ALGEBRAIC[81] - STATES[29]) / CONSTANTS[49];
    RATES[30] = (ALGEBRAIC[75] - STATES[30]) / CONSTANTS[49];
    RATES[31] = (ALGEBRAIC[70] - STATES[31]) / CONSTANTS[49];
    RATES[32] = (ALGEBRAIC[57] - STATES[32]) / CONSTANTS[50];
    RATES[33] = (ALGEBRAIC[79] - STATES[33]) / CONSTANTS[46];
    RATES[34] = (ALGEBRAIC[78] - STATES[34]) / CONSTANTS[46];
    RATES[35] = (ALGEBRAIC[68] - STATES[35]) / CONSTANTS[46];
    RATES[36] = (ALGEBRAIC[51] - STATES[36]) / CONSTANTS[47];
    ALGEBRAIC[66] = CONSTANTS[86] * CONSTANTS[42] * STATES[33] * STATES[29] * ALGEBRAIC[72] * (STATES[0] - STATES[4]);
    ALGEBRAIC[64] = CONSTANTS[86] * CONSTANTS[43] * STATES[34] * STATES[30] * ALGEBRAIC[83] * (STATES[1] - STATES[5]);
    ALGEBRAIC[65] = CONSTANTS[86] * CONSTANTS[44] * STATES[35] * STATES[31] * ALGEBRAIC[77] * (STATES[2] - STATES[6]);
    ALGEBRAIC[58] = CONSTANTS[86] * CONSTANTS[45] * STATES[36] * STATES[32] * ALGEBRAIC[61] * (STATES[3] - STATES[8]);
}

void calc_serca(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[30] = CONSTANTS[87] * ((-CONSTANTS[56]) * pow(STATES[0], 2.0) * (CONSTANTS[53] - CONSTANTS[53] * STATES[37]) + CONSTANTS[54] * CONSTANTS[53] * STATES[37]) * CONSTANTS[7] * 2.0;
    ALGEBRAIC[33] = CONSTANTS[87] * ((-CONSTANTS[56]) * pow(STATES[1], 2.0) * (CONSTANTS[53] - CONSTANTS[53] * STATES[38]) + CONSTANTS[54] * CONSTANTS[53] * STATES[38]) * CONSTANTS[8] * 2.0;
    ALGEBRAIC[34] = CONSTANTS[87] * ((-CONSTANTS[56]) * pow(STATES[2], 2.0) * (CONSTANTS[53] - CONSTANTS[53] * STATES[39]) + CONSTANTS[54] * CONSTANTS[53] * STATES[39]) * CONSTANTS[9] * 2.0;
    ALGEBRAIC[19] = CONSTANTS[87] * ((-CONSTANTS[56]) * pow(STATES[3], 2.0) * (CONSTANTS[53] - CONSTANTS[53] * STATES[40]) + CONSTANTS[54] * CONSTANTS[53] * STATES[40]) * CONSTANTS[1] * 2.0;

    ALGEBRAIC[17] = CONSTANTS[88] * (CONSTANTS[55] * pow(STATES[4], 2.0) * (CONSTANTS[53] - CONSTANTS[53] * STATES[37]) - CONSTANTS[57] * CONSTANTS[53] * STATES[37]) * CONSTANTS[7] * 2.0;
    ALGEBRAIC[15] = CONSTANTS[88] * (CONSTANTS[55] * pow(STATES[5], 2.0) * (CONSTANTS[53] - CONSTANTS[53] * STATES[38]) - CONSTANTS[57] * CONSTANTS[53] * STATES[38]) * CONSTANTS[8] * 2.0;
    ALGEBRAIC[16] = CONSTANTS[88] * (CONSTANTS[55] * pow(STATES[6], 2.0) * (CONSTANTS[53] - CONSTANTS[53] * STATES[39]) - CONSTANTS[57] * CONSTANTS[53] * STATES[39]) * CONSTANTS[9] * 2.0;
    ALGEBRAIC[14] = CONSTANTS[88] * (CONSTANTS[55] * pow(STATES[8], 2.0) * (CONSTANTS[53] - CONSTANTS[53] * STATES[40]) - CONSTANTS[57] * CONSTANTS[53] * STATES[40]) * CONSTANTS[1] * 2.0;

    RATES[37] = 0.5 * ((-ALGEBRAIC[30]) + ALGEBRAIC[17]) / CONSTANTS[7] / CONSTANTS[53];
    RATES[38] = 0.5 * ((-ALGEBRAIC[33]) + ALGEBRAIC[15]) / CONSTANTS[8] / CONSTANTS[53];
    RATES[39] = 0.5 * ((-ALGEBRAIC[34]) + ALGEBRAIC[16]) / CONSTANTS[9] / CONSTANTS[53];
    RATES[40] = 0.5 * ((-ALGEBRAIC[19]) + ALGEBRAIC[14]) / CONSTANTS[1] / CONSTANTS[53];
}

void calc_sodium(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    ALGEBRAIC[93] = ALGEBRAIC[99] + 3.0 * ALGEBRAIC[88] + 3.0 * ALGEBRAIC[62] + ALGEBRAIC[90];
    ALGEBRAIC[93] += ALGEBRAIC[92] + CONSTANTS[92] * ALGEBRAIC[110];  // INab + G_seal * INa_ghk
    ALGEBRAIC[96] = CONSTANTS[80] * CONSTANTS[6] / CONSTANTS[12] * (STATES[42] - STATES[41]) * 1e-06;
    ALGEBRAIC[47] = 1.0 / (1.0 + CONSTANTS[79] * CONSTANTS[81] / pow(STATES[42] + CONSTANTS[81], 2.0));
    RATES[41] = ALGEBRAIC[96] / CONSTANTS[17];
    RATES[42] = ALGEBRAIC[47] * ((-ALGEBRAIC[96]) / CONSTANTS[1] - ALGEBRAIC[93] / (CONSTANTS[1] * CONSTANTS[29]));
}

void calc_stimulus(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    // ALGEBRAIC[56] = fmod(time, CONSTANTS[85]) - CONSTANTS[84];
    // if ((ALGEBRAIC[56] >= 0) && (ALGEBRAIC[56] < CONSTANTS[83])) {
    //     ALGEBRAIC[56] = CONSTANTS[82];
    // } else {
    //     ALGEBRAIC[56] = 0.0;
    // }
    // ALGEBRAIC[56] *= CONSTANTS[34];
    // double past = floor(time / CONSTANTS[85]) * CONSTANTS[85];
    // ALGEBRAIC[56]  = ((time - past >= CONSTANTS[84]) && (time - past <= CONSTANTS[84] + CONSTANTS[83])) ? CONSTANTS[82] * CONSTANTS[34] : 0.;

    // ALGEBRAIC[56] = CONSTANTS[82] * CONSTANTS[34]; // STIM LEVEL * amplitude

    // Dirty hack
    // ALGEBRAIC[56]  = ((time - CONSTANTS[85] >= CONSTANTS[84]) && (time - CONSTANTS[85] <= CONSTANTS[84] + CONSTANTS[83])) ? CONSTANTS[82] * CONSTANTS[34] : 0.;
    // printf("Stimulation: %f\n", ALGEBRAIC[56]);
}


void calc_fluo(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    const double fluo_tot = CONSTANTS[89];
    const double k_on = CONSTANTS[90], k_off = CONSTANTS[91];

    for (int i = 0; i < 5; ++i) {

        int i_fluo = 43 + i;
        int i_Cai  = 4 + i;

        double fluo = STATES[i_fluo];
        double Cai  = STATES[i_Cai];

        RATES[i_fluo] = k_on * Cai * (fluo_tot - fluo) - k_off * fluo;
        RATES[i_Cai] -= RATES[i_fluo];

    }

}


void calc_means(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    RATES[48] = 0;
    RATES[49] = 0;
    RATES[50] = 0;
}


void calc_iseal(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    // ALGEBRAIC[107] = CONSTANTS[92] * STATES[23];  // I_Iseal = G_seal * V

    // I_seal = G_seal * (IK_ghk + INa_ghk + ICl_ghk + IAspartate_ghk + ICa_ghk)
    ALGEBRAIC[107] = CONSTANTS[92] * (ALGEBRAIC[109] + ALGEBRAIC[110] + ALGEBRAIC[111] + ALGEBRAIC[112] + ALGEBRAIC[113]);
}


void calc_ikb(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES)
{
    /* IKb = gKb * (V - EK) */
    double ghk_K = ALGEBRAIC[109];
    double K_ghk_scaler_membrane = CONSTANTS[97];
    ALGEBRAIC[108] = K_ghk_scaler_membrane * ghk_K; // CONSTANTS[93] * (STATES[23] - ALGEBRAIC[105]);
}


void compute_rates_algebraic(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES) {

    calc_ghk_K(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_ghk_Na(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_ghk_Cl(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_ghk_Aspartate(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_ghk_Ca(time, STATES, CONSTANTS, ALGEBRAIC, RATES);

    calc_iseal(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_ikb(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_inab(time, STATES, CONSTANTS, ALGEBRAIC, RATES);

    calc_ical(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_icap(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_stimulus(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_inaca(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_inak(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_nernst(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_ryr(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_serca(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_icab(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_if(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_ik1(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_ikr(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_iks(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_ikur(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_ina(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_it(time, STATES, CONSTANTS, ALGEBRAIC, RATES);

    calc_calcium(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_membrane(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_potassium(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_sodium(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_fluo(time, STATES, CONSTANTS, ALGEBRAIC, RATES);
    calc_means(time, STATES, CONSTANTS, ALGEBRAIC, RATES);

}
