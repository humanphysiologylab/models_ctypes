#include <math.h>

void initialize_states_default(double *STATES) {
    STATES[0] =  7.02128101897185673e-04;
    STATES[1] =  3.94923428392655786e-03;
    STATES[2] =  1.35538532457244482e-01;
    STATES[3] =  1.03674364292988680e-01;
    STATES[4] =  1.90759804527589089e-01;
    STATES[5] =  1.35640688636079511e-02;
    STATES[6] =  2.14063418881809235e-02;
    STATES[7] =  4.45327242854324807e-03;
    STATES[8] =  1.27856586024588575e-01;
    STATES[9] =  5.69999505293381902e-03;
    STATES[10] =  1.83143535034222225e-02;
    STATES[11] =  2.10808768153058460e-04;
    STATES[12] =  3.25814677291117296e-04;
    STATES[13] =  2.33018340557575125e-04;
    STATES[14] =  3.61396062660070427e+00;
    STATES[15] =  7.88607791910409195e-01;
    STATES[16] =  9.15153381546177336e+00;
    STATES[17] =  9.15182798281732346e+00;
    STATES[18] =  5.02305826642838293e-01;
    STATES[19] =  1.13337536953687845e+00;
    STATES[20] = -7.34336366728778671e+01;
    STATES[21] =  2.16850216379767157e-05;
    STATES[22] =  9.98384427312367095e-01;
    STATES[23] =  4.49572164109603364e-02;
    STATES[24] =  3.28512098597005947e-02;
    STATES[25] = 120.0;
    STATES[26] =  1.31290096227093382e-03;
    STATES[27] =  7.49436760722081534e-03;
    STATES[28] =  9.15199678386256998e+00;
    STATES[29] =  3.93548562883350357e-04;
    STATES[30] =  9.58234428284286399e-01;
    STATES[31] =  3.15482710277587786e-01;
    STATES[32] =  2.48034071360795916e-01;
    STATES[33] =  1.89326933812916480e-02;
    STATES[34] =  3.79829335413739144e-02;
    STATES[35] =  1.01974216400706526e-02;
    STATES[36] =  1.37939236359928058e-03;
    STATES[37] =  9.45874848392074696e-01;
    STATES[38] =  5.01323282772066123e-07;
    STATES[39] =  2.01567245823636694e-06;
    STATES[40] =  8.00819151705148946e-01;
}

void initialize_constants_default(double* CONSTANTS) {
    CONSTANTS[0] = 0.024;
    CONSTANTS[1] = 0.0171;
    CONSTANTS[2] = 0.14;
    CONSTANTS[3] = 0.07;
    CONSTANTS[4] = 0.14;
    CONSTANTS[5] = 0.238;
    CONSTANTS[6] = 0.00046;
    CONSTANTS[7] = 5.7e-05;
    CONSTANTS[8] = 0.03;
    CONSTANTS[9] = 1.3;
    CONSTANTS[10] = 0.06;
    CONSTANTS[11] = 3.2e-05;
    CONSTANTS[12] = 0.00333;
    CONSTANTS[13] = 34.0;
    CONSTANTS[14] = 13.8;
    CONSTANTS[15] = 0.0157;
    CONSTANTS[16] = 100.0;
    CONSTANTS[17] = 100.0;
    CONSTANTS[18] = 100.0;
    CONSTANTS[19] = 2.37;
    CONSTANTS[20] = 0.003;
    CONSTANTS[21] = 32.7;
    CONSTANTS[22] = 1.0;
    CONSTANTS[23] = 0.0;
    CONSTANTS[24] = 7.561;
    CONSTANTS[25] = 1.65;
    CONSTANTS[26] = 0.001;
    CONSTANTS[27] = 0.0001;
    CONSTANTS[28] =  3.72425607984805052e-12;
    CONSTANTS[29] = 1.1e-10;
    CONSTANTS[30] =  8.24130542277896849e-13;
    CONSTANTS[31] = 96485.0;
    CONSTANTS[32] = 65.0;
    CONSTANTS[33] = 100.0;
    CONSTANTS[34] = 0.0;
    CONSTANTS[35] = 0.0;
    CONSTANTS[36] =  1.83127823220607955e-14;
    CONSTANTS[37] =  1.63862792221979433e-12;
    CONSTANTS[38] = 100.0;
    CONSTANTS[39] = 10.25;
    CONSTANTS[40] =  3.14159265358979312e+00;
    CONSTANTS[41] =  6.06430000000000033e-04;
    CONSTANTS[42] = 0.11;
    CONSTANTS[43] = 1.8;
    CONSTANTS[44] = 1.0;
    CONSTANTS[45] = 0.0;
    CONSTANTS[46] = 0.0;
    CONSTANTS[47] = 0.00027;
    CONSTANTS[48] = 1.35e-07;
    CONSTANTS[49] = 7.5e-09;
    CONSTANTS[50] = 0.9;
    CONSTANTS[51] = 1.8;
    CONSTANTS[52] = 5.4;
    CONSTANTS[53] = 140.0;
    CONSTANTS[54] = 0.009;
    CONSTANTS[55] = 0.0548;
    CONSTANTS[56] = 0.1;
    CONSTANTS[57] = 0.0525;
    CONSTANTS[58] = 0.002;
    CONSTANTS[59] = 0.035;
    CONSTANTS[60] = 0.0035;
    CONSTANTS[61] = 0.01833;
    CONSTANTS[62] = 0.045;
    CONSTANTS[63] = 23.0;
    CONSTANTS[64] = 0.000597;
    CONSTANTS[65] = 3.15;
    CONSTANTS[66] = 0.000384;
    CONSTANTS[67] = 0.00359;
    CONSTANTS[68] = 1.3;
    CONSTANTS[69] = 12.29;
    CONSTANTS[70] = 87.5;
    CONSTANTS[71] = 1.57;
    CONSTANTS[72] = 0.27;
    CONSTANTS[73] = 0.35;
    CONSTANTS[74] = 1.26;
    CONSTANTS[75] = 1.5;
    CONSTANTS[76] = 0.0025;
    CONSTANTS[77] = 600.0;
    CONSTANTS[78] = 15.0;
    CONSTANTS[79] = 150.0;
    CONSTANTS[80] = 0.0471;
    CONSTANTS[81] = 0.0005;
    CONSTANTS[82] = 2.35;
    CONSTANTS[83] = 0.165;
    CONSTANTS[84] = 8314.0;
    CONSTANTS[85] = 310.0;
    CONSTANTS[86] = 5.348e-06;
    CONSTANTS[87] = 1.7;
    CONSTANTS[88] = 15.0;
    CONSTANTS[89] = 1.0;
    CONSTANTS[90] = 2.6;
    CONSTANTS[91] = 0.0053114;
    CONSTANTS[92] = 0.45;
    CONSTANTS[93] = 1.787;
    CONSTANTS[94] = 0.5;
    CONSTANTS[95] = 0.005;
    CONSTANTS[96] = 0.06;
    CONSTANTS[97] = 25.0;
    CONSTANTS[98] = -12.5;
    CONSTANTS[99] = 5.0;
    CONSTANTS[100] = 50.0;
    CONSTANTS[101] = 1000.0;
    CONSTANTS[102] =  (1.00000+ 0.500000*CONSTANTS[23])*0.0196000;
    CONSTANTS[103] =  ( ( CONSTANTS[40]*pow(CONSTANTS[39], 2.00000))*CONSTANTS[38])*1.00000e-15;
    CONSTANTS[104] = pow( CONSTANTS[81]*1.00000, 1.60000);
    CONSTANTS[105] =  ( (1.00000+ 0.500000*CONSTANTS[23])*(1.00000 -  0.500000*CONSTANTS[34]))*CONSTANTS[47];
    CONSTANTS[106] =  ( (1.00000+ 0.500000*CONSTANTS[23])*(1.00000 -  0.500000*CONSTANTS[34]))*CONSTANTS[48];
    CONSTANTS[107] =  ( (1.00000+ 0.500000*CONSTANTS[23])*(1.00000 -  0.500000*CONSTANTS[34]))*CONSTANTS[49];
    CONSTANTS[108] = 1.00000 - CONSTANTS[50];
    CONSTANTS[109] = (CONSTANTS[31]/CONSTANTS[84])/CONSTANTS[85];
    CONSTANTS[110] =  ( (1.00000+CONSTANTS[34])* pow((CONSTANTS[52]/5.40000), 1.0 / 2))*CONSTANTS[57];
    CONSTANTS[111] =  CONSTANTS[59]* pow((CONSTANTS[52]/5.40000), 1.0 / 2);
    CONSTANTS[112] =  ((1.00000+CONSTANTS[34])+ 2.00000*CONSTANTS[23])*CONSTANTS[60];
    CONSTANTS[113] =  ((1.00000+CONSTANTS[34])+ 2.00000*CONSTANTS[23])*CONSTANTS[60];
    CONSTANTS[114] =  ( ( (1.00000 -  0.500000*CONSTANTS[34])*(1.00000+ 2.00000*CONSTANTS[23]))*(1.00000+ 0.200000*CONSTANTS[35]))*CONSTANTS[62];
    CONSTANTS[115] =  CONSTANTS[63]*(1.00000 -  0.100000*CONSTANTS[34]);
    CONSTANTS[116] =  (1.00000+ 0.400000*CONSTANTS[34])*CONSTANTS[65];
    CONSTANTS[117] =  11.0000*(1.00000 -  0.250000*CONSTANTS[23]);
    CONSTANTS[118] = (exp(CONSTANTS[53]/67.3000) - 1.00000)/7.00000;
    CONSTANTS[119] =  CONSTANTS[76]*CONSTANTS[34];
    CONSTANTS[120] =  (1.00000 -  0.700000*CONSTANTS[34])*CONSTANTS[83];
    CONSTANTS[121] =  (2.50000 -  1.25000*CONSTANTS[23])*0.000246000;
    CONSTANTS[122] =  ((10.0000+ 20.0000*CONSTANTS[34])+ ( 10.0000*CONSTANTS[23])*(1.00000 - CONSTANTS[34]))*1.00000;
    CONSTANTS[123] =  ( 0.0539000*0.0100000)*CONSTANTS[103];
    CONSTANTS[124] = pow( CONSTANTS[81]*1.00000, 1.60000);
    CONSTANTS[125] =  (1.00000/CONSTANTS[109])*log(CONSTANTS[78]/CONSTANTS[79]);
    CONSTANTS[126] =  0.650000*CONSTANTS[103];
    CONSTANTS[127] = 1.00000 - CONSTANTS[42];
    CONSTANTS[128] =  (( 0.00165000*CONSTANTS[126])/CONSTANTS[123])*0.100000;
    CONSTANTS[129] =  (( 0.00460000*CONSTANTS[126])/CONSTANTS[123])*0.100000;
    CONSTANTS[130] =  0.0200000*CONSTANTS[103];
    CONSTANTS[131] =  0.0350000*CONSTANTS[103];
    CONSTANTS[132] = (CONSTANTS[85] - 310.000)/10.0000;
    CONSTANTS[133] =  (CONSTANTS[126]/CONSTANTS[130])*0.0134000;
    CONSTANTS[134] =  (CONSTANTS[126]/CONSTANTS[130])*0.0374000;
    CONSTANTS[135] =  (CONSTANTS[126]/CONSTANTS[131])*0.140000;

    for (int i = 0; i < 18; ++i) { // scalers
        CONSTANTS[136 + i] = 1;
    }

}