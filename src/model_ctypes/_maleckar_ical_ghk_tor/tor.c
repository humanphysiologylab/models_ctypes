#include "./tor.h"
#include "math.h"


inline double calculate_dss(double v) {
    return (v < 31.4978) ? 1.0763 * exp(-1.0070 * exp(-0.0829 * v)) : 1; // ToR_ORd
}

inline double calculate_fss(double v) {
    return 1.0 / (1.0 + exp((v + 19.58) / 3.696));
}

inline double calculate_tff(double v) {
    return 7.0 + 1.0 / (0.0045 * exp(-(v + 20.0) / 10.0) +
                        0.0045 * exp((v + 20.0) / 10.0));
}

inline double calculate_tfs(double v) {
    return 1000.0 + 1.0 / (0.000035 * exp(-(v + 5.0) / 4.0) +
                           0.000035 * exp((v + 5.0) / 6.0));
}

inline double calculate_tfcaf(double v) {
    return 7.0 + 1.0 / (0.04 * exp(-(v - 4.0) / 7.0) +
                        0.04 * exp((v - 4.0) / 7.0));
}

inline double calculate_tfcas(double v) {
    return 100.0 + 1.0 / (0.00012 * exp(-v / 3.0) + 0.00012 * exp(v / 7.0));
}

inline double calculate_Afcaf(double v) {
    return 0.3 + 0.6 / (1.0 + exp((v - 10.0) / 10.0));
}

inline double calculate_td(double v) {
    return 0.6 + 1.0 / (exp(-0.05 * (v + 6.0)) + exp(0.09 * (v + 14.0)));
}

double calculate_ical_tor(double* CONSTANTS, double* RATES,
                          double* STATES, double* ALGEBRAIC) {

    double S_V = STATES[0];
    double dss = calculate_dss(S_V); // dimensionless
    double fss = calculate_fss(S_V); // dimensionless
    double td = calculate_td(S_V) /*ms*/ / 1000.; // sec
    double tff = calculate_tff(S_V) /*ms*/ / 1000.; // sec
    double tfs = calculate_tfs(S_V) /*ms*/ / 1000.; // sec
    double tfcaf = calculate_tfcaf(S_V) /*ms*/ / 1000.; // sec
    double tfcas = calculate_tfcas(S_V) /*ms*/ / 1000.; // sec
    double Afcaf = calculate_Afcaf(S_V); // dimensionless

    double S_d_ord = STATES[30];
    double S_ff    = STATES[31];
    double S_fs    = STATES[32];
    double S_fcaf  = STATES[33];
    double S_fcas  = STATES[34];
    double S_jca   = STATES[35];
    double S_ffp   = STATES[36];
    double S_fcafp = STATES[37];
    double S_nca   = STATES[38];

    double Ca_d  = STATES[6];
    double Ca_c  = STATES[18];

    double R = CONSTANTS[0];
    double T = CONSTANTS[1];
    double F = CONSTANTS[2];


    double Aff    = 0.6;
    double Afs    = 1.0 - Aff;
    double f      = Aff * S_ff + Afs * S_fs;
    double fcass  = fss;
    double Afcas  = 1.0 - Afcaf;
    double fca    = Afcaf * S_fcaf + Afcas * S_fcas;
    double tjca   = 75.0 /*ms*/ / 1000.; // sec
    double tffp   = 2.5 * tff;
    double tfcafp = 2.5 * tfcaf;
    double Kmn    = 0.002; // mM
    double k2n    = 1000.0  /*1/ms*/ * 1000.; // sec
    double km2n   = S_jca * 1.0  /1/*ms*/ * 1000.; // sec
    double anca   = 1.0 / (k2n / km2n + pow(1.0 + Kmn / Ca_d, 4.0));

    double gamma_ca_i = 0.63, gamma_ca_o = 0.63;
    double phi_CaL = 4.0 * (S_V * F * F / (R * T)) *
                     (gamma_ca_i * Ca_d * exp(2.0 * (S_V * F / (R * T))) - gamma_ca_o * Ca_c) /
                     (exp(2.0 * (S_V * F / (R * T))) - 1.0);

    // Rush-Larsen method was used in ORd
    // Here I use simple Euler (like everywhere else)
    RATES[30] = (dss   - S_d_ord) / td;     // NextSt->d = Sc->dss - (Sc->dss - S->d_ord) * exp(-Par->dt / td);
    RATES[31] = (fss   - S_ff   ) / tff;    // NextSt->ff = Sc->fss - (Sc->fss - S->ff) * exp(-Par->dt / tff);
    RATES[32] = (fss   - S_fs   ) / tfs;    // NextSt->fs = Sc->fss - (Sc->fss - S->fs) * exp(-Par->dt / tfs);
    RATES[33] = (fcass - S_fcaf ) / tfcaf;  // NextSt->fcaf = fcass - (fcass - S->fcaf) * exp(-Par->dt / tfcaf);
    RATES[34] = (fcass - S_fcas   ) / tfcas;  // NextSt->fcas = fcass - (fcass - S->fcas) * exp(-Par->dt / tfcas);
    RATES[35] = (fcass - S_jca  ) / tjca;   // NextSt->jca = fcass - (fcass - S->jca) * exp(-Par->dt / tjca);
    RATES[36] = (fss   - S_ffp  ) / tffp;   // NextSt->ffp = Sc->fss - (Sc->fss - S->ffp) * exp(-Par->dt / tffp);
    RATES[37] = (fcass - S_fcafp) / tfcafp; // NextSt->fcafp = fcass - (fcass - S->fcafp) * exp(-Par->dt / tfcafp);
    RATES[38] = anca * k2n - S_nca * km2n;  // NextSt->nca = anca * k2n / km2n - (anca * k2n / km2n - S->nca) * exp(-km2n * Par->dt);


    double Cm      = CONSTANTS[3];
    double pca_tor = CONSTANTS[51];

    double ical_tor_max = pca_tor * phi_CaL; // Col / m^3
    double ical_tor     = ical_tor_max * S_d_ord * (f  * (1.0 - S_nca) + S_jca * fca  * S_nca);

    ical_tor = ical_tor * Cm /*nF*/; // in ORd model [ICaL] = [uA/uF]

    return ical_tor;
}
