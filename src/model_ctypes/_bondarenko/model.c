#include "math.h"

void compute_rates_algebraic(double VOI,
                             const double* STATES, const double* CONSTANTS,
                             double* ALGEBRAIC, double* RATES)
//void fun(common& cmn, int const& /* neq */, double const& /* time */, double* y, double* ydot, double const& istim)
{//C
  //C Cell Geometry Parameters
  //C
  if(VOI <= 1.0f)
  {
      ALGEBRAIC[113] = 80;
  }else{
      ALGEBRAIC[113] = 0;
  }
  //ALGEBRAIC[114] = VOI;
 ALGEBRAIC[6]/*acap*/= 0.5931e-4f;
 double vcell = 22.1200e-6f;
 ALGEBRAIC[7]/*vmyo*/ = 13.9600e-6f;
 double vcyt = 13.9600e-6f;
 ALGEBRAIC[8]/*vjsr*/= 0.1349e-6f;
 ALGEBRAIC[9]/*vnsr*/= 2.5523e-6f;
 ALGEBRAIC[10]/*vss*/= 1.4850e-9f * 2.0f;
  double vcav = 0.02f * vcell;
  double vecav = 0.04f * vcell;
  //C
  //C Standard Ionic Concentrations
  //C
  double cko = 5400.0f;
  double cnao = 140000.0f;
  double ccao = 1800.0f;
  //C
  //C SR Parameters
  //C
  double icasmax = 2.0f;
  double t1 = 0.04f;
  double t2 = -0.1f / icasmax;
  double v1 = 4.5f * 0.09f;
  double v2 = 1.74e-5f * 3.0f;
  double v3 = 0.306f * 2.63f;
  //C      Kmup  =  0.5
  double kmup = (0.38f * (1 - STATES[104]) + 0.25f * STATES[104]) * 0.815f;
  double ttr = 20.0f;
  double txfer = 8.0f * 1.3f;
  double kap = 0.006075f;
  double kam = 0.07125f;
  double kbp = 0.00405f;
  double kbm = 0.965f;
  double kcp = 0.009f;
  double kcm = 0.0008f;
  double kapp = kap * 5.0f;
  double kamp = kam * 3.0f;
  double kbpp = kbp * 5.0f;
  double kbmp = kbm * 3.0f;
  double kcpp = kcp * 50.0f;
  double kcmp = kcm * 30.0f;
  double f_ryr = 0.001f;
  //C
  //C L-type Ca2+ Channel Parameters
  //C
  double kpcmax = 0.11662f * 2.0f;
  double kpchalf = 10.0f;
  double kpcb = 0.0005f * 4.0f * 1.2f;
  double kpcf = 2.5f * 8.0f * 2.0f;
  double kco = 1.0f;
  double kcop = 4.0f;
  double koc = 1.0f;
  double sfc4i1 = 0.01f;
  double sfc4i2 = 0.002f;
  double sfc4i3 = 1.0f;
  double sfoi1 = 1.0f;
  double sfoi2 = 0.001f;
  double sfi1i3 = 0.001f;
  double sfi2i3 = 1.0f;
  double sfica = 0.1127f;
  double sficap = 0.2351f;
  //C
  //C T-type Ca2+ Channel Parameters
  //C
  double k_catpv = 9.5f * exp(STATES[1] / 15.0f);
  double k_catmv = 0.008f * exp(-1.0f * STATES[1] / 18.0f) * 4.0f;
  double k_catpo = 4.0f;
  double k_catmo = 0.025f * exp(-1.0f * STATES[1] / 34.0f);
  double k_catpi = 0.070f;
  double k_catmi = 0.0014f * 0.7f;
  double f_cat = 0.2f;
  double h_cat = 0.5f;
  double g_cat = 0.012725f;
  double g_catp = 0.02545f;
  double ficatcav = 0.68f;
  //C
  //C Buffering Parameters
  //C
  double cltrpntot = 70.0f;
  double chtrpntot = 140.0f;
  double khtrpnp = 0.00237f;
  double khtrpnm = 0.000032f;
  double kltrpnp = 0.0327f;
  //C      kltrpnm   =     0.0196
  double kltrpnm = 0.0196f * (1 - STATES[105]) + 0.0294f * STATES[105];
  double ccmdntot = 50.0f;
  double ccsqntot = 15000.0f * 0.6f;
  double kmcmdn = 0.238f;
  double kmcsqn = 800.0f;
  //C
  //C Membrane Current Parameters
  //C
 ALGEBRAIC[11]/*cm*/= 1.0f;
 ALGEBRAIC[12]/*f*/= 96.5f;
  double t = 298.0f;
  double r = 8.314f;
  double factor = r * t / ALGEBRAIC[12]/*f*/;
  double ifactor =ALGEBRAIC[12]/*f*// (r * t);
  double gna = 14.4f;
  double gnap = 18.0f;
  double knaca = 275.0f * 1.0f;
  double kmna = 87500.0f;
  double kmca = 1380.0f;
  double ksat = 0.27f;
  double nu = 0.35f;
  double inakmax = 4.0f * 0.5f;
  //C      Kmnai   = 21000.0
  double kmnai = 18800.0f * (1 - STATES[129]) + 13600.0f * STATES[129];
  double kmko = 1500.0f;
  double sigma = (exp(cnao / 67300.0f) - 1.0f) / 7.0f;
  double ipcamax = 0.051f * 1.5f * 2.0f;
  double kmpca = 0.5f;
  double gcab = 0.000284f * 0.75f;
  double gnab = 0.0063f * 1.0f;
  double dito = 7.0f;
  //C      Gks     =     0.00575
  double gks = 0.0f;
  double gkto = 0.14067f;
  double gktop = 0.14067f * (1.0f - 0.387f);
  double gkur1 = 0.0f;
  double gkur2 = 0.05766f;
  double gkur2p = 0.07496f;
  double gkur3 = 0.0428f;
  double gkur4 = 0.05766f;
  //C
  //C HERG current parameters
  //C
  double kf = 0.023761f;
  double kb = 0.036778f;
  double gkr = 0.078f * 10.0f;
  //C
  //C Calcium activated chloride current
  //C
  double gclca = 10.0f;
  double poclcamax = 0.2f;
  double kmclca = 10.0f;
  double ecl = -40.0f;
  //C
  //C Small-conductance calcium-activated K+ current
  //C
  double k_cakf1 = 0.2f * STATES[2] * 0.6f;
  double k_cakf2 = 0.16f * STATES[2] * 0.6f;
  double k_cakf3 = 0.08f * STATES[2] * 0.6f;
  double k_cakb1 = 0.08f;
  double k_cakb2 = 0.08f;
  double k_cakb3 = 0.2f;
  double k_cakfo1 = 0.16f;
  double k_cakbo1 = 1.0f;
  double k_cakfo2 = 1.2f;
  double k_cakbo2 = 0.1f;
  double g_cak = (1.1e-6f * STATES[1] * STATES[1] * STATES[1] + 1.0e-4f * STATES[1] * STATES[
    1] + 0.0041f * STATES[1] + 0.1872f) * 4.0f;
  //C
  //C Beta1-adrenergic receptor module
  //C
  //C      L       =   0.0
  double rb1tot = 0.0103f;
  double fcavb1 = 0.01f;
  double fecavb1 = 0.5f;
  double fcytb1 = 1 - fcavb1 - fecavb1;
  double gstot = 2.054f;
  double fcavgs = 0.4f;
  double fecavgs = 0.4f;
  double fcytgs = 1 - fcavgs - fecavgs;
  double kb1l = 0.567f;
  double kb1h = 0.0617f;
  double kb1c = 2.86f;
  double kpkap = 8.1e-7f;
  double kpkam = 0.25f * kpkap;
  double kbarkp = 0.3f * kpkap;
  double kbarkm = kpkam;
  double factgsgi = 0.04f;
  double kact1gs = 4.9e-3f;
  double kact2gs = 2.6e-4f;
  double khydgs = 8.0e-4f;
  double kreasgs = 1.2f;
  //C
  //C Beta2-adrenergic receptor module
  //C
  double pi = 3.14159026f;
  double rb2tot = 0.0053f;
  double fcavb2 = 0.99f;
  double fecavb2 = 1 - fcavb2;
  double gitot = 10.086f;
  double fcavgi = 0.99f;
  double fecavgi = 1 - fcavgi;
  double kb2l = 1.053f;
  double kb2h = 0.0118f;
  double kb2c = 5.86f;
  double kb2f = 0.0189f;
  double kb2a = 28.79f;
  double kact1gi = 4.0e-3f * 0.5f;
  double kact2gi = 5.0e-5f;
  double khydgi = khydgs;
  double kreasgi = kreasgs;
  //C
  //C Adenylyl Cyclase module
  //C
  double kmatp = 340.0f;
  double atp = 5000.0f;
  double actot = 0.02622f;
  double fac56ac47 = 0.74f;
  double fcavac56 = 0.0875f;
  double fecavac47 = 0.1648f;
  double kac56mgsa = 0.0852f;
  double hac56gsa = 1.357f;
  double vac56gbg = 1.430f;
  double kac56mgsbg = 0.003793f;
  double hac56gsbg = 1.0842f;
  double ac56bas = 0.0377f;
  double af56 = 0.0511335f;
  double kac47mgsa = 0.05008f;
  double hac47gsa = 1.1657f;
  double vac47gbg = 1.35f;
  double kac47mgsbg = 0.004466f;
  double hac47gsbg = 0.870f;
  double ac47bas = 0.04725f;
  double af47 = 0.009283f;
  double kac56mgsgi = 0.482f;
  double kac56mgi = 0.0465f;
  double hac56gsgi = 0.662f;
  double vac56gsgi = 0.857f;
  //C
  //C Phosphodiesterase Module
  //C
  double ibmx = 1.0e-20f;
  double hibmxpde2 = 1.0f;
  double kibmxmpde2 = 29.5f;
  double hibmxpde3 = 1.0f;
  double kibmxmpde3 = 5.1f;
  double hibmxpde4 = 1.0f;
  double kibmxmpde4 = 16.2f;
  double kfpdep = 1.96e-5f;
  double kbpdep = 1.02e-5f;
  double deltakpde34 = 3.0f;
  double kpde2 = 0.020f;
  double kmpde2 = 33.0f;
  double kpde3 = 0.0025f;
  double kmpde3 = 0.44f;
  double kpde4 = 0.0035f;
  double kmpde4 = 1.4f;
  double kpde8 = 0.0f;
  double kmpde8 = 0.15f;
  //double fpdepart = 0.2f;
  //double rpartpde23 = 0.570f;
  //double rpartpde34 = 0.748f;
  double pde2tot = 0.034610f;
  double pde3tot = 0.010346f;
  double pde4tot = 0.026687f;
  double pde8tot = 0.016805f;
  double fcavpde2 = 0.06608f;
  double fecavpde2 = 2 * fcavpde2;
  double fcytpde2 = 1 - fcavpde2 - fecavpde2;
  double fcavpde3 = 0.29814f;
  double fecavpde3 = 0.0f;
  double fcytpde3 = 1 - fcavpde3 - fecavpde3;
  double fcavpde4 = 0.05366f;
  double fecavpde4 = 2 * fcavpde4;
  double fcytpde4 = 1 - fcavpde4 - fecavpde4;
  double fcavpde8 = 0.06608f;
  double fecavpde8 = 2 * fcavpde8;
  double fcytpde8 = 1 - fcavpde8 - fecavpde8;
  //C
  //C cAMP-PKA module
  //C
  double pkatot = 0.5176f;
  double fcavpka = 0.08f;
  double fecavpka = 0.20f;
  double fcytpka = 1 - fcavpka - fecavpka;
  double pkitot = 2 * 0.2f * pkatot;
  double fcavpki = fcavpka;
  double fecavpki = fecavpka;
  double fcytpki = fcytpka;
  double kpkaif1 = 0.0056f;
  double kpkai1 = 2.9f;
  double kpkaif2 = kpkaif1;
  double kpkai2 = 2.9f;
  //C      kpkaif3   = 7.99e-4
  double kpkaif3 = 0.0026f;
  double kpkai3 = 1.3f;
  double kpkaiif1 = kpkaif1;
  double kpkaii1 = 2.5f;
  double kpkaiif2 = kpkaif1;
  double kpkaii2 = 2.5f;
  //C      kpkaiif3  = 0.00511
  double kpkaiif3 = kpkaif3;
  double kpkaii3 = kpkai3;
  double kpkif = 0.050f;
  double kpki = 2.6e-4f;
  //C
  //C Protein Phosphatase module
  //C
  double pp1cyttot = 0.2f;
  double inhib1cyttot = 0.08543f;
  double kinhib1 = 1.0e-3f;
  double kpkainhib1 = 1.08f;
  double kmpkainhib1 = 1.5f;
  double kpp2ainhib1pp2acyt = 0.00308f;
  double kmpp2ainhib1 = 0.001f;
  //C
  //C Caveolae
  //C
 ALGEBRAIC[51]/*rb1cavtot*/= fcavb1 * rb1tot * vcell / vcav;
 ALGEBRAIC[93]/*rb2cavtot*/= fcavb2 * rb2tot * vcell / vcav;
 ALGEBRAIC[38]/*gscavabg*/= fcavgs * gstot * vcell / vcav - STATES[47] - STATES[49];
 ALGEBRAIC[91]/*gicavabg*/= fcavgi * gitot * vcell / vcav - STATES[157] - STATES[158];
  double rb1cavnptot =ALGEBRAIC[51]/*rb1cavtot*/- STATES[45] - STATES[46];
  double rb2cavnptot =ALGEBRAIC[93]/*rb2cavtot*/- STATES[155] - STATES[156];
 ALGEBRAIC[99]/*acavb2i*/= (kb2l + ALGEBRAIC[80]/*l*/) * (kb2f + ALGEBRAIC[80]/*l*/) / kb2l;
 ALGEBRAIC[100]/*bcavb2i*/=ALGEBRAIC[91]/*gicavabg*/* (kb2f + ALGEBRAIC[80]/*l*/) - STATES[155] * (kb2f + ALGEBRAIC[80]/*l*/) + kb2a *
    kb2f * (1.0f +ALGEBRAIC[80]/*l*// kb2l);
 ALGEBRAIC[101]/*ccavb2i*/= -1.0f * STATES[155] * kb2a * kb2f;

 ALGEBRAIC[97]/*rb2cavpkaf*/= (-ALGEBRAIC[100]/*bcavb2i*/ + sqrt(pow(ALGEBRAIC[100]/*bcavb2i*/,2) - 4.0f *
   ALGEBRAIC[99]/*acavb2i*/* ALGEBRAIC[101]/*ccavb2i*/)) / (2.0f * ALGEBRAIC[99]/*acavb2i*/);

 ALGEBRAIC[98]/*gicavf*/=ALGEBRAIC[91]/*gicavabg*// (1.0f +ALGEBRAIC[97]/*rb2cavpkaf*/* (1.0f / kb2a +ALGEBRAIC[80]/*l*// (kb2a * kb2f)));
  //double lrb2cavpka =ALGEBRAIC[80]/*l*/*ALGEBRAIC[97]/*rb2cavpkaf*// kb2l;
 ALGEBRAIC[95]/*rb2gicavpka*/=ALGEBRAIC[97]/*rb2cavpkaf*/*ALGEBRAIC[98]/*gicavf*// kb2a;
 ALGEBRAIC[96]/*lrb2gicavpka*/=ALGEBRAIC[80]/*l*/*ALGEBRAIC[97]/*rb2cavpkaf*/*ALGEBRAIC[98]/*gicavf*// (kb2a * kb2f);
  double acavb2s = (kb1h + ALGEBRAIC[80]/*l*/) * (kb2h + ALGEBRAIC[80]/*l*/);
  double bcavb2s = (ALGEBRAIC[80]/*l*/ + kb1h) * (ALGEBRAIC[80]/*l*/ + kb2h) * (rb1cavnptot +
    rb2cavnptot) + (kb1c * kb1h +ALGEBRAIC[80]/*l*/* kb1c * kb1h / kb1l) * (ALGEBRAIC[80]/*l*/ +
    kb2h) + (kb2c * kb2h +ALGEBRAIC[80]/*l*/* kb2c * kb2h / kb2l) * (ALGEBRAIC[80]/*l*/ + kb1h) -
   ALGEBRAIC[38]/*gscavabg*/* (ALGEBRAIC[80]/*l*/ + kb1h) * (ALGEBRAIC[80]/*l*/ + kb2h);
  double ccavb2s = (ALGEBRAIC[80]/*l*/ + kb1h) * (kb2c * kb2h +ALGEBRAIC[80]/*l*/* kb2c * kb2h /
    kb2l) * (rb1cavnptot - ALGEBRAIC[38]/*gscavabg*/) + (ALGEBRAIC[80]/*l*/ + kb2h) * (kb1c * kb1h +
   ALGEBRAIC[80]/*l*/* kb1c * kb1h / kb1l) * (rb2cavnptot - ALGEBRAIC[38]/*gscavabg*/) + (kb1c *
    kb1h +ALGEBRAIC[80]/*l*/* kb1c * kb1h / kb1l) * (kb2c * kb2h +ALGEBRAIC[80]/*l*/* kb2c * kb2h /
    kb2l);
  double dcavb2s = -ALGEBRAIC[38]/*gscavabg*/ * (kb1c * kb1h +ALGEBRAIC[80]/*l*/* kb1c * kb1h /
    kb1l) * (kb2c * kb2h +ALGEBRAIC[80]/*l*/* kb2c * kb2h / kb2l);
  //C
  //C Cubic equation solution
  //C
  double pcav = bcavb2s / acavb2s;
  double qcav = ccavb2s / acavb2s;
  double rcav = dcavb2s / acavb2s;
  double acav = (3.0f * qcav - pcav * pcav) / 3.0f;
  double bcav = (2.0f * pcav * pcav * pcav - 9.0f * pcav * qcav +
    27.0f * rcav) / 27.0f;
  double dcav = acav * acav * acav / 27.0f + bcav * bcav / 4.0f;
  double mcav = pow((-bcav * 0.5f + sqrt(dcav)), (1.0f / 3.0f));
  double ncav = pow((-bcav * 0.5f - sqrt(dcav)), (1.0f / 3.0f));
  //C
  double y1cav;
  double y2cav;
  double y3cav;
  if (dcav > 0.0f) {
    y1cav = mcav + ncav;
    y2cav = 0.0f;
    y3cav = 0.0f;
  }
  if (dcav == 0.0f) {
    y1cav = mcav + ncav;
    y2cav = (mcav + ncav) / (-2.0f);
    y3cav = (mcav + ncav) / (-2.0f);
  }
  double ficav;
  if (dcav < 0.0f) {
    if (bcav > 0.0f) {
      ficav = acos(-sqrt((bcav * bcav / 4.0f) / (-acav *
        acav * acav / 27.0f)));
    }
    else {
      ficav = acos(sqrt((bcav * bcav / 4.0f) / (-acav *
        acav * acav / 27.0f)));
    }
    y1cav = 2.0f * sqrt(-acav / 3.0f) * cos(ficav);
    y2cav = 2.0f * sqrt(-acav / 3.0f) * cos(ficav + 2.0f * pi / 3.0f);
    y3cav = 2.0f * sqrt(-acav / 3.0f) * cos(ficav + 4.0f * pi / 3.0f);
  }
  //C
  double z1cav = y1cav - pcav / 3.0f;
  double z2cav = y2cav - pcav / 3.0f;
  double z3cav = y3cav - pcav / 3.0f;
  double gscavf = fmax(z1cav, fmax(z2cav, z3cav));
  double rb1cavnpf = rb1cavnptot / (1.0f +ALGEBRAIC[80]/*l*// kb1l + gscavf * (ALGEBRAIC[80]/*l*/ / (
    kb1c * kb1h) + 1.0f / kb1c));
  double rb2cavnpf = rb2cavnptot / (1.0f +ALGEBRAIC[80]/*l*// kb2l + gscavf * (ALGEBRAIC[80]/*l*/ / (
    kb2c * kb2h) + 1.0f / kb2c));
  double lrb1cavnp =ALGEBRAIC[80]/*l*/* rb1cavnpf / kb1l;
  double rb1gscavnp = rb1cavnpf * gscavf / kb1c;
  double lrb1gscavnp =ALGEBRAIC[80]/*l*/* rb1cavnpf * gscavf / (kb1c * kb1h);
  double lrb2cavnp =ALGEBRAIC[80]/*l*/* rb2cavnpf / kb2l;
  double rb2gscavnp = rb2cavnpf * gscavf / kb2c;
  double lrb2gscavnp =ALGEBRAIC[80]/*l*/* rb2cavnpf * gscavf / (kb2c * kb2h);
  double ac56cav = fcavac56 * fac56ac47 * actot * (vcell / vcav);
  double kcavac56 = af56 * (ac56bas + pow(STATES[47], hac56gsa) / (
    kac56mgsa + pow(STATES[47], hac56gsa))) * (1 + vac56gbg * pow(STATES[48],
    hac56gsbg) / (kac56mgsbg + pow(STATES[48], hac56gsbg))) * (1.0f -
    (1.0f - vac56gsgi * pow(STATES[47], hac56gsgi) / (kac56mgsgi +
    pow(STATES[47], hac56gsgi))) * STATES[157] / (kac56mgi + STATES[157]));
 ALGEBRAIC[54]/*pde2cavtot*/= (1.0f - pow(ibmx, hibmxpde2) / (kibmxmpde2 + pow(ibmx,
    hibmxpde2))) * fcavpde2 * pde2tot * (vcell / vcav);
 ALGEBRAIC[55]/*pde3cavtot*/= (1.0f - pow(ibmx, hibmxpde3) / (kibmxmpde3 + pow(ibmx,
    hibmxpde3))) * fcavpde3 * pde3tot * (vcell / vcav);
 ALGEBRAIC[56]/*pde4cavtot*/= (1.0f - pow(ibmx, hibmxpde4) / (kibmxmpde4 + pow(ibmx,
    hibmxpde4))) * fcavpde4 * pde4tot * (vcell / vcav);
 ALGEBRAIC[63]/*pde8cavtot*/= fcavpde8 * pde8tot * (vcell / vcav);
 ALGEBRAIC[71]/*pkacav*/= fcavpka * pkatot * (vcell / vcav);
 ALGEBRAIC[35]/*rccavf*/= 2.0f *ALGEBRAIC[71]/*pkacav*/- STATES[80] - STATES[81] - STATES[82];
 ALGEBRAIC[74]/*pkicavf*/= fcavpki * pkitot * (vcell / vcav) - STATES[84];
  double kpkaiib1 = kpkaiif1 * kpkaii1;
  double kpkaiib2 = kpkaiif2 * kpkaii2;
  double kpkaiib3 = kpkaiif3 / kpkaii3;
  double kpkib = kpkif * kpki;
  //C
  //C Extracaveolae
  //C
 ALGEBRAIC[52]/*rb1ecavtot*/= fecavb1 * rb1tot * vcell / vecav;
 ALGEBRAIC[94]/*rb2ecavtot*/= fecavb2 * rb2tot * vcell / vecav;
 ALGEBRAIC[39]/*gsecavabg*/= fecavgs * gstot * vcell / vecav - STATES[52] - STATES[54];
 ALGEBRAIC[92]/*giecavabg*/= fecavgi * gitot * vcell / vecav - STATES[161] - STATES[162];
  double rb1ecavnptot =ALGEBRAIC[52]/*rb1ecavtot*/- STATES[50] - STATES[51];
  double rb2ecavnptot =ALGEBRAIC[94]/*rb2ecavtot*/- STATES[159] - STATES[160];
 ALGEBRAIC[106]/*aecavb2i*/= (kb2l + ALGEBRAIC[80]/*l*/) * (kb2f + ALGEBRAIC[80]/*l*/) / kb2l;
 ALGEBRAIC[107]/*becavb2i*/=ALGEBRAIC[92]/*giecavabg*/* (kb2f + ALGEBRAIC[80]/*l*/) - STATES[159] * (kb2f + ALGEBRAIC[80]/*l*/) + kb2a *
    kb2f * (1.0f +ALGEBRAIC[80]/*l*// kb2l);
 ALGEBRAIC[108]/*cecavb2i*/= -1.0f * STATES[159] * kb2a * kb2f;
 ALGEBRAIC[104]/*rb2ecavpkaf*/= (-ALGEBRAIC[107]/*becavb2i*/ + sqrt(pow(ALGEBRAIC[107]/*becavb2i*/,2) - 4.0f *
   ALGEBRAIC[106]/*aecavb2i*/* ALGEBRAIC[108]/*cecavb2i*/)) / (2.0f * ALGEBRAIC[106]/*aecavb2i*/);
 ALGEBRAIC[105]/*giecavf*/=ALGEBRAIC[92]/*giecavabg*// (1.0f +ALGEBRAIC[104]/*rb2ecavpkaf*/* (1.0f / kb2a +ALGEBRAIC[80]/*l*// (
    kb2a * kb2f)));
  //double lrb2ecavpka =ALGEBRAIC[80]/*l*/*ALGEBRAIC[104]/*rb2ecavpkaf*// kb2l;
 ALGEBRAIC[102]/*rb2giecavpka*/=ALGEBRAIC[104]/*rb2ecavpkaf*/*ALGEBRAIC[105]/*giecavf*// kb2a;
 ALGEBRAIC[103]/*lrb2giecavpka*/=ALGEBRAIC[80]/*l*/*ALGEBRAIC[104]/*rb2ecavpkaf*/*ALGEBRAIC[105]/*giecavf*// (kb2a * kb2f);
  double aecavb2s = (kb1h + ALGEBRAIC[80]/*l*/) * (kb2h + ALGEBRAIC[80]/*l*/);
  double becavb2s = (ALGEBRAIC[80]/*l*/ + kb1h) * (ALGEBRAIC[80]/*l*/ + kb2h) * (rb1ecavnptot +
    rb2ecavnptot) + (kb1c * kb1h +ALGEBRAIC[80]/*l*/* kb1c * kb1h / kb1l) * (ALGEBRAIC[80]/*l*/ +
    kb2h) + (kb2c * kb2h +ALGEBRAIC[80]/*l*/* kb2c * kb2h / kb2l) * (ALGEBRAIC[80]/*l*/ + kb1h) -
   ALGEBRAIC[39]/*gsecavabg*/* (ALGEBRAIC[80]/*l*/ + kb1h) * (ALGEBRAIC[80]/*l*/ + kb2h);
  double cecavb2s = (ALGEBRAIC[80]/*l*/ + kb1h) * (kb2c * kb2h +ALGEBRAIC[80]/*l*/* kb2c * kb2h /
    kb2l) * (rb1ecavnptot - ALGEBRAIC[39]/*gsecavabg*/) + (ALGEBRAIC[80]/*l*/ + kb2h) * (kb1c * kb1h +
   ALGEBRAIC[80]/*l*/* kb1c * kb1h / kb1l) * (rb2ecavnptot - ALGEBRAIC[39]/*gsecavabg*/) + (kb1c *
    kb1h +ALGEBRAIC[80]/*l*/* kb1c * kb1h / kb1l) * (kb2c * kb2h +ALGEBRAIC[80]/*l*/* kb2c * kb2h /
    kb2l);
  double decavb2s = -ALGEBRAIC[39]/*gsecavabg*/ * (kb1c * kb1h +ALGEBRAIC[80]/*l*/* kb1c * kb1h /
    kb1l) * (kb2c * kb2h +ALGEBRAIC[80]/*l*/* kb2c * kb2h / kb2l);
  //C
  //C Cubic equation solution
  //C
  double pecav = becavb2s / aecavb2s;
  double qecav = cecavb2s / aecavb2s;
  double recav = decavb2s / aecavb2s;
  double aecav = (3.0f * qecav - pecav * pecav) / 3.0f;
  double becav = (2.0f * pecav * pecav * pecav - 9.0f * pecav *
    qecav + 27.0f * recav) / 27.0f;
  double decav = aecav * aecav * aecav / 27.0f + becav * becav / 4.0f;
  double mecav = pow((-becav * 0.5f + sqrt(decav)), (1.0f / 3.0f));
  double necav = pow((-becav * 0.5f - sqrt(decav)), (1.0f / 3.0f));
  //C
  double y1ecav;
  double y2ecav;
  double y3ecav;
  if (decav > 0.0f) {
    y1ecav = mecav + necav;
    y2ecav = 0.0f;
    y3ecav = 0.0f;
  }
  if (decav == 0.0f) {
    y1ecav = mecav + necav;
    y2ecav = (mecav + necav) / (-2.0f);
    y3ecav = (mecav + necav) / (-2.0f);
  }
  double fiecav;
  if (decav < 0.0f) {
    if (becav > 0.0f) {
      fiecav = acos(-sqrt((becav * becav / 4.0f) / (
        -aecav * aecav * aecav / 27.0f)));
    }
    else {
      fiecav = acos(sqrt((becav * becav / 4.0f) / (-aecav *
        aecav * aecav / 27.0f)));
    }
    y1ecav = 2.0f * sqrt(-aecav / 3.0f) * cos(fiecav);
    y2ecav = 2.0f * sqrt(-aecav / 3.0f) * cos(fiecav +
      2.0f * pi / 3.0f);
    y3ecav = 2.0f * sqrt(-aecav / 3.0f) * cos(fiecav +
      4.0f * pi / 3.0f);
  }
  //C
  double z1ecav = y1ecav - pecav / 3.0f;
  double z2ecav = y2ecav - pecav / 3.0f;
  double z3ecav = y3ecav - pecav / 3.0f;
  double gsecavf = fmax(z1ecav, fmax(z2ecav, z3ecav));
  double rb1ecavnpf = rb1ecavnptot / (1.0f +ALGEBRAIC[80]/*l*// kb1l + gsecavf * (ALGEBRAIC[80]/*l*/ /
    (kb1c * kb1h) + 1.0f / kb1c));
  double rb2ecavnpf = rb2ecavnptot / (1.0f +ALGEBRAIC[80]/*l*// kb2l + gsecavf * (ALGEBRAIC[80]/*l*/ /
    (kb2c * kb2h) + 1.0f / kb2c));
  double lrb1ecavnp =ALGEBRAIC[80]/*l*/* rb1ecavnpf / kb1l;
  double rb1gsecavnp = rb1ecavnpf * gsecavf / kb1c;
  double lrb1gsecavnp =ALGEBRAIC[80]/*l*/* rb1ecavnpf * gsecavf / (kb1c * kb1h);
  double lrb2ecavnp =ALGEBRAIC[80]/*l*/* rb2ecavnpf / kb2l;
  double rb2gsecavnp = rb2ecavnpf * gsecavf / kb2c;
  double lrb2gsecavnp =ALGEBRAIC[80]/*l*/* rb2ecavnpf * gsecavf / (kb2c * kb2h);
  double ac47ecav = fecavac47 * (1.0f - fac56ac47) * actot * vcell / vecav;
  double kecavac47 = af47 * (ac47bas + pow(STATES[52], hac47gsa) / (
    kac47mgsa + pow(STATES[52], hac47gsa))) * (1.0f + vac47gbg *
    pow(STATES[53], hac47gsbg) / (kac47mgsbg + pow(STATES[53],
    hac47gsbg)));
 ALGEBRAIC[57]/*pde2ecavtot*/= (1.0f - (pow(ibmx, hibmxpde2) / (kibmxmpde2 +
    pow(ibmx, hibmxpde2)))) * fecavpde2 * pde2tot * (vcell /
    vecav);
 ALGEBRAIC[58]/*pde3ecavtot*/= (1.0f - (pow(ibmx, hibmxpde3) / (kibmxmpde3 +
    pow(ibmx, hibmxpde3)))) * fecavpde3 * pde3tot * (vcell /
    vecav);
 ALGEBRAIC[59]/*pde4ecavtot*/= (1.0f - (pow(ibmx, hibmxpde4) / (kibmxmpde4 +
    pow(ibmx, hibmxpde4)))) * fecavpde4 * pde4tot * (vcell /
    vecav);
 ALGEBRAIC[64]/*pde8ecavtot*/= fecavpde8 * pde8tot * (vcell / vecav);
 ALGEBRAIC[72]/*pkaecav*/= fecavpka * pkatot * (vcell / vecav);
 ALGEBRAIC[36]/*rcecavf*/= 2.0f *ALGEBRAIC[72]/*pkaecav*/- STATES[86] - STATES[87] - STATES[88];
 ALGEBRAIC[75]/*pkiecavf*/= fecavpki * pkitot * (vcell / vecav) - STATES[90];
  kpkaiib1 = kpkaiif1 * kpkaii1;
  kpkaiib2 = kpkaiif2 * kpkaii2;
  kpkaiib3 = kpkaiif3 / kpkaii3;
  kpkib = kpkif * kpki;
  //C
  //C Cytosol
  //C
    ALGEBRAIC[53]/*rb1cyttot*/ = fcytb1 * rb1tot * vcell / vcyt;
    ALGEBRAIC[40]/*gscytabg*/ = fcytgs * gstot * vcell / vcyt - STATES[57] - STATES[59];
  double rb1cytnptot = ALGEBRAIC[53]/*rb1cyttot*/ - STATES[55] - STATES[56];
  double acytb1 = (kb1l + ALGEBRAIC[80]/*l*/) * (kb1h + ALGEBRAIC[80]/*l*/) / kb1l;
  double bcytb1 = ALGEBRAIC[40]/*gscytabg*/ * (kb1h + ALGEBRAIC[80]/*l*/) - rb1cytnptot * (kb1h + ALGEBRAIC[80]/*l*/) +
    kb1c * kb1h * (1.0f +ALGEBRAIC[80]/*l*// kb1l);
  double ccytb1 = -rb1cytnptot * kb1c * kb1h;
  double rb1cytnpf = (-bcytb1 + sqrt(pow(bcytb1,2) - 4.0f *
    acytb1 * ccytb1)) / (2.0f * acytb1);
  double gscytf = ALGEBRAIC[40]/*gscytabg*/ / (1.0f + rb1cytnpf * (1.0f / kb1c +ALGEBRAIC[80]/*l*// (
    kb1c * kb1h)));
  double lrb1cytnp =ALGEBRAIC[80]/*l*/* rb1cytnpf / kb1l;
  double rb1gscytnp = rb1cytnpf * gscytf / kb1c;
  double lrb1gscytnp =ALGEBRAIC[80]/*l*/* rb1cytnpf * gscytf / (kb1c * kb1h);
  double ac56cyt = (1.0f - fcavac56) * fac56ac47 * actot * (vcell / vcyt);
  double kcytac56 = af56 * (ac56bas + pow(STATES[57], hac56gsa) / (
    kac56mgsa + pow(STATES[57], hac56gsa))) * (1.0f + vac56gbg *
    pow(STATES[58], hac56gsbg) / (kac56mgsbg + pow(STATES[58],
    hac56gsbg)));
  double ac47cyt = (1.0f - fecavac47) * (1.0f - fac56ac47) * actot *
    vcell / vcyt;
  double kcytac47 = af47 * (ac47bas + pow(STATES[57], hac47gsa) / (
    kac47mgsa + pow(STATES[57], hac47gsa))) * (1.0f + vac47gbg *
    pow(STATES[58], hac47gsbg) / (kac47mgsbg + pow(STATES[58],
    hac47gsbg)));
    ALGEBRAIC[60]/*pde2cyttot*/ = (1.0f - (pow(ibmx, hibmxpde2) / (kibmxmpde2 +
    pow(ibmx, hibmxpde2)))) * fcytpde2 * pde2tot * (vcell /
    vcyt);
    ALGEBRAIC[61]/*pde3cyttot*/ = (1.0f - (pow(ibmx, hibmxpde3) / (kibmxmpde3 +
    pow(ibmx, hibmxpde3)))) * fcytpde3 * pde3tot * (vcell /
    vcyt);
    ALGEBRAIC[62]/*pde4cyttot*/ = (1.0f - (pow(ibmx, hibmxpde4) / (kibmxmpde4 +
    pow(ibmx, hibmxpde4)))) * fcytpde4 * pde4tot * (vcell /
    vcyt);
    ALGEBRAIC[65]/*pde8cyttot*/ = fcytpde8 * pde8tot * (vcell / vcyt);
    ALGEBRAIC[73]/*pkacyt*/ = fcytpka * pkatot * (vcell / vcyt);
    ALGEBRAIC[37]/*rccytf*/= 2.0f * ALGEBRAIC[73]/*pkacyt*/ - STATES[92] - STATES[93] - STATES[94];
    ALGEBRAIC[76]/*pkicytf*/ = fcytpki * pkitot * (vcell / vcyt) - STATES[96];
  double kpkaib1 = kpkaif1 * kpkai1;
  double kpkaib2 = kpkaif2 * kpkai2;
  double kpkaib3 = kpkaif3 / kpkai3;
  kpkib = kpkif * kpki;
  //C
  //C Protein phosphatase module
  //C
  double inhib1cytf = inhib1cyttot - STATES[100];
  double ainh1 = 1.0f;
  double kinh1 = 0.0f;//uninitialized in m9da.f
  double binh1 = kinh1 + pp1cyttot - STATES[100];
  double cinh1 = -STATES[100] * kinhib1;
  double inhib1cytp = -binh1 / 2.0f + sqrt(pow(binh1,2) -
    4.0f * ainh1 * cinh1) / 2.0f;
    ALGEBRAIC[89]/*pp1cytf*/ = pp1cyttot * kinhib1 / (kinhib1 + inhib1cytp);
  //C
  //C cAMP fluxes
  //C
  double jcavecav = 5.00e-12f;
  double jcavcyt = 7.500e-11f;
  double jecavcyt = 9.00e-12f;
  //C
  //C PLB module
  //C
  double kplbpka = 1.08917e-4f;
  double kplbmpka = 4.90970e-4f;
  double kplbpp1 = 4.41956e-5f;
  double kplbmpp1 = 1.69376e-2f;
  //C
  //C TnI module
  //C
  double pp2acyt = 0.0607843f;
  double ktnipka = 2.47254e-5f;
  double ktnimpka = 2.71430e-5f;
  double ktnipp2a = 8.65898e-5f;
  double ktnimpp2a = 8.01420e-1f;
  //C
  //C ICaL module
  //C
  double pp1cav = 0.1f;
  double pp2acav = 0.1f;
  double ppcav = pp1cav + pp2acav;
  double icaltot = 0.0273f;
 ALGEBRAIC[90]/*ficalcav*/= 0.2f;
  double icalcav =ALGEBRAIC[90]/*ficalcav*/* icaltot * vcell / vcav;
  double icalecav = (1.0f - ALGEBRAIC[90]/*ficalcav*/) * icaltot * vcell / vecav;
  double kicalpka = 2.0e-5f * 0.87f;
  double kicalpp = 1.55e-7f * 1.5f;
  double kicalmpka = 0.5f;
  double kicalmpp = 0.2f;
  //C
  //C ICaT module
  //C
  double kicatpka = 1.74e-5f * 0.5f * 20.0f;
  double kicatpp = 2.33e-7f * 200.0f;
  double kicatmpka = 5.0f;
  double kicatmpp = 0.1f;
  //C
  //C ICaK module
  //C
  double kicakpka = 1.74e-5f * 0.085f;
  double kicakpp = 2.33e-7f * 1.25f;
  double kicakmpka = 2.5f;
  double kicakmpp = 0.01f;
  //C
  //C RyR module
  //C
  double ryrtot = 0.1993f;
  double ryrecav = ryrtot * vcell / vecav;
  double kryrpka = 3.85e-5f * 1.5f;
  double kryrpp = 1.925e-4f * 1.5f;
  double kryrmpka = 0.5f;
  double kryrmpp = 0.05f;
  //C
  //C INa module
  //C
  //double inatot = 0.0f;
  //double inacav = inatot * vcell / vcav;
  double kinapka = 6.8400e-6f;
  double kinapp = 1.9804e-5f;
  double kinampka = 5.49415e-3f;
  double kinampp = 0.393025f;
  //C
  //C INaK module
  //C
  double kinakpka = 3.053e-6f;
  double kinakpp = 1.8491e-5f;
  double kinakmpka = 1.1001e-3f;
  double kinakmpp = 5.7392f;
  //C
  //C IKur module
  //C
  double pp1ecav = 0.1f;
  double kikurpka = 6.9537e-6f;
  double kikurpp = 3.1700e-5f;
  double kikurmpka = 0.138115f;
  double kikurmpp = 0.23310f;
  //C
  //C IK1 module
  //C
  double kik1pka = 1.9954e-5f;
  double kik1pp = 9.0968e-5f;
  double kik1mpka = 0.027623f;
  double kik1mpp = 0.023310f;
  //C
  //C IKto module
  //C
  double kiktopka = 4.38983e-5f;
  double kiktopp = 9.09678e-5f;
  double kiktompka = 0.27623f;
  double kiktompp = 0.23310f;
  //C
    ALGEBRAIC[41]/*campcell*/ = (STATES[97] * vcav + STATES[98] * vecav + STATES[99] * vcyt) / vcell;
    ALGEBRAIC[42]/*ccell*/ = (STATES[83] * vcav + STATES[89] * vecav + STATES[95] * vcyt) / vcell;
    ALGEBRAIC[43]/*acactcell*/ = (RATES[60] * vcav + RATES[61] * vecav + (RATES[62] +
    RATES[63]) * vcyt) / vcell;
 ALGEBRAIC[66]/*pdeact2*/= (RATES[66] * vcav + RATES[71] * vecav + RATES[76] * vcyt) / vcell;
 ALGEBRAIC[67]/*pdeact3*/= (RATES[67] * vcav + RATES[72] * vecav + RATES[77] * vcyt) / vcell;
 ALGEBRAIC[68]/*pdeact4*/= (RATES[68] * vcav + RATES[73]* vecav + RATES[78] * vcyt) / vcell;
 ALGEBRAIC[69]/*pdeact8*/= (RATES[101] * vcav + RATES[102] * vecav + RATES[103] * vcyt) / vcell;
 ALGEBRAIC[84]/*pde2mem*/= (RATES[66] * vcav + RATES[71] * vecav) / vcell;
 ALGEBRAIC[85]/*pde3mem*/= (RATES[67] * vcav + RATES[72] * vecav) / vcell;
 ALGEBRAIC[86]/*pde4mem*/= (RATES[68] * vcav + RATES[73] * vecav) / vcell;
 ALGEBRAIC[87]/*pde8mem*/= (RATES[101] * vcav + RATES[102] * vecav) / vcell;
    ALGEBRAIC[44]/*pdeactcell*/ =ALGEBRAIC[66]/*pdeact2*/+ALGEBRAIC[67]/*pdeact3*/+ALGEBRAIC[68]/*pdeact4*/+ ALGEBRAIC[69]/*pdeact8*/;
    ALGEBRAIC[88]/*pdemem*/ =ALGEBRAIC[84]/*pde2mem*/+ALGEBRAIC[85]/*pde3mem*/+ALGEBRAIC[86]/*pde4mem*/+ ALGEBRAIC[87]/*pde8mem*/;
    ALGEBRAIC[70]/*pkaactcell*/ = (RATES[83] * vcav + RATES[89] * vecav + RATES[95] *
    vcyt) / vcell;
 ALGEBRAIC[77]/*rb1ppka*/= (STATES[45] * vcav + STATES[50] * vecav + STATES[55] * vcyt) / vcell;
 ALGEBRAIC[78]/*rb1pbark*/= (STATES[46] * vcav + STATES[51] * vecav + STATES[56] * vcyt) / vcell;
    ALGEBRAIC[79]/*rb1ptot*/ =ALGEBRAIC[77]/*rb1ppka*/+ ALGEBRAIC[78]/*rb1pbark*/;
 ALGEBRAIC[81]/*rb2ppka*/= (STATES[155] * vcav + STATES[159] * vecav) / vcell;
 ALGEBRAIC[82]/*rb2pbark*/= (STATES[156] * vcav + STATES[160] * vecav) / vcell;
    ALGEBRAIC[83]/*rb2ptot*/ =ALGEBRAIC[81]/*rb2ppka*/+ ALGEBRAIC[82]/*rb2pbark*/;
  //C
  //C***********************************************************************
  //C            Calcium fluxes
  //C***********************************************************************
  //C      Pryr     = exp((STATES[1]-5.0]*(STATES[1]-5.0]/(-648.))
 ALGEBRAIC[0]/*jrelt*/= v1 * (STATES[18] + STATES[19] + STATES[118] + STATES[119]) * (STATES[6] - STATES[7]) * STATES[44];
 ALGEBRAIC[1]/*jtrt*/= (STATES[3] - STATES[6]) / ttr;
 ALGEBRAIC[2]/*jxfert*/= (STATES[7] - STATES[2]) / txfer;
 ALGEBRAIC[3]/*jleak*/= v2 * (STATES[3] - STATES[2]);
 ALGEBRAIC[4]/*jup*/= v3 * STATES[2] * STATES[2] / (kmup * kmup + STATES[2] * STATES[2]);
  double temp12 = khtrpnp * STATES[2] * (chtrpntot - STATES[5]) - khtrpnm * STATES[5];
  double temp13 = kltrpnp * STATES[2] * (cltrpntot - STATES[4]) - kltrpnm * STATES[4];
 ALGEBRAIC[5]/*jtrpn*/= temp12 + temp13;
  //C***********************************************************************
  //C             Ionic currents
  //C***********************************************************************
  //C
  //C STATES[1]   : Membrane potential (mV)
  //C
  //double v = STATES[1];
  //C
  //C Icab   : Calcium background current
  //C
  double ecan = 0.5f * factor * log(ccao / STATES[2]);
 ALGEBRAIC[17]/*icab*/= gcab * (STATES[1] - ecan);
  //C
  //C Ipca   : Calcium pump current
  //C
 ALGEBRAIC[19]/*ipca*/= ipcamax * STATES[2] * STATES[2] / (kmpca * kmpca + STATES[2] * STATES[2]);
  //C
  //C Inaca  : Na-Ca exchange current
  //C
  double tempa1 = knaca / (kmna * kmna * kmna + cnao * cnao * cnao);
  double tempa2 = 1.0f / (kmca + ccao);
  double tempa3 = 1.0f / (1.0f + ksat * exp((nu - 1.0f) * STATES[1] * ifactor));
  double tempa41 = exp(nu * STATES[1] * ifactor) * STATES[23] * STATES[23] * STATES[23] * ccao;
  double tempa42 = exp((nu - 1.0f) * STATES[1] * ifactor) * cnao *
    cnao * cnao * STATES[2];
  double tempa4 = tempa41 - 2.0f * tempa42;
 ALGEBRAIC[18]/*inaca*/= tempa1 * tempa2 * tempa3 * tempa4;
  //C
  //C Ica    : L-type calcium current
  //C
 ALGEBRAIC[23]/*icasc*/=ALGEBRAIC[90]/*ficalcav*/* (sfica * STATES[8] + sficap * STATES[107]) * (STATES[1] - 52.0f);
 ALGEBRAIC[24]/*icase*/= (1.0f - ALGEBRAIC[90]/*ficalcav*/) * (sfica * STATES[137] + sficap * STATES[146]) * (STATES[
    1] - 52.0f);
 ALGEBRAIC[16]/*icas*/=ALGEBRAIC[23]/*icasc*/+ ALGEBRAIC[24]/*icase*/;
  //C
  //C Icat   : T-type calcium current
  //C
 ALGEBRAIC[110]/*icatc*/= ficatcav * (g_cat * STATES[168] * (1 - STATES[175]) + g_catp * STATES[
    168] * STATES[175]) * (STATES[1] - 35.0f);
 ALGEBRAIC[111]/*icats*/= (1.0f - ficatcav) * (g_cat * STATES[168] * (1 - STATES[176]) +
    g_catp * STATES[168] * STATES[176]) * (STATES[1] - 35.0f);
 ALGEBRAIC[109]/*icat*/=ALGEBRAIC[110]/*icatc*/+ ALGEBRAIC[111]/*icats*/;
  //C
  //C Ina    : Na fast current (Luo and Rudy, 1994]
  //C
  double ena = factor * log((0.9f * cnao + 0.1f * cko) / (0.9f *
    STATES[23] + 0.1f * STATES[24]));
 ALGEBRAIC[20]/*ina*/= (gna * STATES[37] + gnap * STATES[123]) * (STATES[1] - ena);
  //C
  //C Inab   : Na background current
  //C
 ALGEBRAIC[21]/*inab*/= gnab * (STATES[1] - ena);
  //C
  //C Inak    : Na-K exchange current
  //C
  double tempa5 = 1.0f / (1.0f + pow((kmnai / STATES[23]), 3.0f));
  double tempa6 = cko / (cko + kmko);
  double tempa11 = 1.0f + 0.1245f * exp(-0.1f * STATES[1] * ifactor) +
    0.0365f * sigma * exp(-1.0f * STATES[1] * ifactor);
 ALGEBRAIC[25]/*inak*/= inakmax * tempa5 * tempa6 / tempa11;
  //C
  //C Ikto    : transient outward current (Liu and Rasmusson, 1997]
  //C
  double ek = factor * log(cko / STATES[24]);
 ALGEBRAIC[26]/*ikto*/= gkto * STATES[25] * STATES[25] * STATES[25] * STATES[26] * (STATES[1] - ek) * (1 - STATES[
    134]) + gktop * STATES[135] * STATES[135] * STATES[135] * STATES[136] * (STATES[1] - ek) *
    STATES[134];
  //C
  //C Ik1    : Time independent K+ current (Rasmusson et al. 1990]
  //C
  //C      tempa7 = 0.3397*cKo/(cKo+210.0]
  //C      tempa8 = 1.0+exp(0.0448*(STATES[1]-EK))
  //C      tempa7p = 0.11394*cKo/(cKo+210.0]
  //C      tempa8p = 1.0+exp(0.0298*(STATES[1]-EK))
  //C      Ik1 = (tempa7*STATES[133]/tempa8
  //C     *    + tempa7*[1-STATES[133])/tempa8]*(STATES[1]-EK)
  //C
  double tempa7 = 1.02f / (1.0f + exp(0.2385f * (STATES[1] - ek - 59.215f)));
  double tempa8 = (0.8f * exp(0.08032f * (STATES[1] - ek + 5.476f)) +
    exp(0.06175f * (STATES[1] - ek - 594.31f))) / (1.0f + exp(
    -0.5143f * (STATES[1] - ek + 4.753f)));
  double gk1 = 1.7f * 0.27f * sqrt(cko / 5400.0f);
  double gk1p = 1.7f * 0.27f * sqrt(cko / 5400.0f);
 ALGEBRAIC[27]/*ik1*/= (gk1 * STATES[133] * tempa7 / (tempa7 + tempa8) + gk1p * (1 - STATES[
    133]) * tempa7 / (tempa7 + tempa8)) * (STATES[1] - ek);
  //C
  //C Iks    : Delayed rectifier K+ current (Rasmusson et al. 1990]
  //C
 ALGEBRAIC[28]/*iks*/= gks * STATES[27] * STATES[27] * (STATES[1] - ek);
  //C
  //C Ikur   : Ultra-rapidly activating delayed rectifier Ikur
  //C                  (Zhou et al., 1998]
  //C
 ALGEBRAIC[31]/*ikur1*/= gkur1 * STATES[28] * STATES[29] * (STATES[1] - ek);
 ALGEBRAIC[32]/*ikur2*/= (gkur2 * STATES[35] * STATES[36] * STATES[132] + gkur2p * STATES[130] * STATES[
    131] * (1 - STATES[132])) * (STATES[1] - ek);
 ALGEBRAIC[33]/*ikur3*/= gkur3 * STATES[43] * (STATES[1] - ek);
 ALGEBRAIC[34]/*ikur4*/= gkur4 * STATES[184] * STATES[185] * (STATES[1] - ek);
 ALGEBRAIC[29]/*ikur*/=ALGEBRAIC[31]/*ikur1*/+ALGEBRAIC[32]/*ikur2*/+ALGEBRAIC[33]/*ikur3*/+ ALGEBRAIC[34]/*ikur4*/;
  //C
  //C Ikr    : HERG current (Wang et al., 1997]
  //C
  double ekr = factor * log((0.98f * cko + 0.02f * cnao) / (
    0.98f * STATES[24] + 0.02f * STATES[23]));
 ALGEBRAIC[30]/*ikr*/= gkr * STATES[33] * (STATES[1] - ekr);
  //C
  //C Iclca  : calcium-activated chloride current (Xu et al., 2002]
  //C
  double tempa9 = poclcamax / (1 + exp((46.7f - STATES[1]) / 7.8f));
  double tempa10 = STATES[2] / (STATES[2] + kmclca);
 ALGEBRAIC[22]/*iclca*/= gclca * tempa9 * tempa10 * (STATES[1] - ecl);
  //C
  //C ICaK : small-conductance calcium-activated K+ current
  //C
 ALGEBRAIC[112]/*icak*/= g_cak * (STATES[181] + STATES[182]) * (1 - STATES[183]);
  //C
  //C***********************************************************************
  //C
  //C STATES[1]  membrane potential
  //C
  double sum_i =ALGEBRAIC[16]/*icas*/+ALGEBRAIC[17]/*icab*/+ALGEBRAIC[18]/*inaca*/+ALGEBRAIC[19]/*ipca*/+ALGEBRAIC[20]/*ina*/+ALGEBRAIC[21]/*inab*/+ALGEBRAIC[22]/*iclca*/+
   ALGEBRAIC[25]/*inak*/+ALGEBRAIC[26]/*ikto*/+ALGEBRAIC[27]/*ik1*/+ALGEBRAIC[28]/*iks*/+ALGEBRAIC[29]/*ikur*/+ALGEBRAIC[30]/*ikr*/- ALGEBRAIC[113]/*istim*/ +ALGEBRAIC[109]/*icat*/+ ALGEBRAIC[112]/*icak*/;
  RATES[1] = -sum_i / ALGEBRAIC[11]/*cm*/;
  //C
  //C STATES[2]  intracellular calcium Cai
  //C
 ALGEBRAIC[13]/*bi*/= 1.0f / (1 + ccmdntot * kmcmdn / ((kmcmdn + STATES[2]) * (kmcmdn + STATES[2])));
  double temp1 = 0.5f * (ALGEBRAIC[17]/*icab*/ - 2.0f *ALGEBRAIC[18]/*inaca*/+ALGEBRAIC[19]/*ipca*/+ALGEBRAIC[23]/*icasc*/+ ALGEBRAIC[109]/*icat*/) *
   ALGEBRAIC[6]/*acap*// (ALGEBRAIC[7]/*vmyo*/ * ALGEBRAIC[12]/*f*/);
  RATES[2] =ALGEBRAIC[13]/*bi*/* (ALGEBRAIC[3]/*jleak*/ +ALGEBRAIC[2]/*jxfert*/-ALGEBRAIC[4]/*jup*/-ALGEBRAIC[5]/*jtrpn*/- temp1);
  //C
  //C STATES[3]  network SR calcium Cansr
  //C
  RATES[3] = (ALGEBRAIC[4]/*jup*/ - ALGEBRAIC[3]/*jleak*/) * ALGEBRAIC[7]/*vmyo*/ /ALGEBRAIC[9]/*vnsr*/-ALGEBRAIC[1]/*jtrt*/*ALGEBRAIC[8]/*vjsr*// ALGEBRAIC[9]/*vnsr*/;
  //C
  //C STATES[4]  cLTRPNca
  //C
  RATES[4] = kltrpnp * STATES[2] * (cltrpntot - STATES[4]) - kltrpnm * STATES[4];
  //C
  //C STATES[5]  cHTRPNca
  //C
  RATES[5] = khtrpnp * STATES[2] * (chtrpntot - STATES[5]) - khtrpnm * STATES[5];
  //C
  //C     Partial contributions
  //C
  double alpha = 0.4f * exp((STATES[1] + 15.0f) / 15.0f);
  double alphap = 0.4f * exp((STATES[1] + 15.0f + 20.0f) / 15.0f);
  double beta = 0.13f * exp(-1.0f * (STATES[1] + 15.0f) / 18.0f);
  //C
  double gammac = kpcmax * STATES[2] / (kpchalf + STATES[2]);
  //C
  //C STATES[6] cCajsr junction SR calcium concentartion
  //C
  double temp2 = (kmcsqn + STATES[6]) * (kmcsqn + STATES[6]);
 ALGEBRAIC[15]/*bjsr*/= 1.0f / (1.0f + ccsqntot * kmcsqn / temp2);
  RATES[6] =ALGEBRAIC[15]/*bjsr*/* (ALGEBRAIC[1]/*jtrt*/ - ALGEBRAIC[0]/*jrelt*/);
  //C
  //C STATES[7] cCass  subspace calcium concentration
  //C
  double temp3 = (kmcmdn + STATES[7]) * (kmcmdn + STATES[7]);
 ALGEBRAIC[14]/*bss*/= 1.0f / (1.0f + ccmdntot * kmcmdn / temp3);
  double temp4 =ALGEBRAIC[0]/*jrelt*/*ALGEBRAIC[8]/*vjsr*//ALGEBRAIC[10]/*vss*/-ALGEBRAIC[2]/*jxfert*/* ALGEBRAIC[7]/*vmyo*/ / ALGEBRAIC[10]/*vss*/;
  RATES[7] =ALGEBRAIC[14]/*bss*/* (temp4 -ALGEBRAIC[24]/*icase*/*ALGEBRAIC[6]/*acap*// (2.0f *ALGEBRAIC[10]/*vss*/* ALGEBRAIC[12]/*f*/));
  //C
  //C               ICaLcav
  //C
  //C STATES[8]     O  Ca channel variable
  //C
  RATES[8] = kco * STATES[106] - koc * STATES[8] + sfoi1 * kpcb * STATES[13] -
    sfoi1 * gammac * STATES[8] + sfoi2 * alpha * STATES[14] - sfoi2 * kpcf * STATES[
    8] - kicalpka * STATES[83] * STATES[8] / (kicalmpka + icalcav * STATES[8]) + (
    alpha / alphap) * kicalpp * ppcav * STATES[107] / (kicalmpp +
    icalcav * STATES[107]);
  //C
  //C STATES[9]    C1  Ca channel variable
  //C
  //C      RATES[9] = beta*STATES[10]-4.0*alpha*STATES[9]
  //C
  //C STATES[10]   C2  Ca channel variable
  //C
  RATES[10] = 4.0f * alpha * STATES[9] - beta * STATES[10] + 2.0f * beta * STATES[
    11] - 3.0f * alpha * STATES[10] - kicalpka * STATES[83] * STATES[10] / (
    kicalmpka + icalcav * STATES[10]) + (alphap * alphap * kcop / (alpha *
    alpha * kco)) * kicalpp * ppcav * STATES[109] / (kicalmpp + icalcav *
    STATES[109]);
  //C
  //C STATES[11]   C3  Ca channel variable
  //C
  RATES[11] = 3.0f * alpha * STATES[10] - 2.0f * beta * STATES[11] + 3.0f *
    beta * STATES[12] - 2.0f * alpha * STATES[11] - kicalpka * STATES[83] * STATES[11] / (
    kicalmpka + icalcav * STATES[11]) + (alphap * kcop / (alpha * kco)) *
    kicalpp * ppcav * STATES[110] / (kicalmpp + icalcav * STATES[110]);
  //C
  //C STATES[12]   C4  Ca channel variable
  //C
  RATES[12] = 2.0f * alpha * STATES[11] - 3.0f * beta * STATES[12] + 4.0f *
    beta * STATES[106] - alpha * STATES[12] + 4.0f * sfc4i1 * kpcb * beta * STATES[
    13] - sfc4i1 * alpha * gammac * (kco / koc) * STATES[12] + 4.0f *
    sfc4i2 * beta * STATES[14] - sfc4i2 * kpcf * (kco / koc) * STATES[12] +
    4.0f * sfc4i3 * beta * kpcb * STATES[15] - sfc4i3 * gammac * kpcf * (
    kco / koc) * STATES[12] - kicalpka * STATES[83] * STATES[12] / (kicalmpka +
    icalcav * STATES[12]) + (kcop / kco) * kicalpp * ppcav * STATES[111] / (
    kicalmpp + icalcav * STATES[111]);
  //C
  //C STATES[106] CP Ca channel variable
  //C
  RATES[106] = alpha * STATES[12] - 4.0f * beta * STATES[106] + koc * STATES[8] -
    kco * STATES[106] - kicalpka * STATES[83] * STATES[106] / (kicalmpka + icalcav *
    STATES[106]) + (alpha * kcop / (alphap * kco)) * kicalpp * ppcav * STATES[
    115] / (kicalmpp + icalcav * STATES[115]);
  //C
  //C STATES[13]   I1  Ca channel variable
  //C
  RATES[13] = sfoi1 * gammac * STATES[8] - sfoi1 * kpcb * STATES[13] + sfi1i3 *
    alpha * STATES[15] - sfi1i3 * kpcf * STATES[13] + sfc4i1 * alpha * gammac *
    (kco / koc) * STATES[12] - 4.0f * sfc4i1 * beta * kpcb * STATES[13] -
    kicalpka * STATES[83] * STATES[13] / (kicalmpka + icalcav * STATES[13]) + (
    alpha / alphap) * kicalpp * ppcav * STATES[112] / (kicalmpp +
    icalcav * STATES[112]);
  //C
  //C STATES[14]   I2  Ca channel variable
  //C
  RATES[14] = sfoi2 * kpcf * STATES[8] - sfoi2 * alpha * STATES[14] + sfi2i3 *
    kpcb * STATES[15] - sfi2i3 * gammac * STATES[14] + sfc4i2 * kpcf * (kco /
    koc) * STATES[12] - 4.0f * sfc4i2 * beta * STATES[14] - kicalpka * STATES[83] * STATES[
    14] / (kicalmpka + icalcav * STATES[14]) + kicalpp * ppcav * STATES[113] / (
    kicalmpp + icalcav * STATES[113]);
  //C
  //C STATES[15]   I3  Ca channel variable
  //C
  RATES[15] = sfi1i3 * kpcf * STATES[13] - sfi1i3 * alpha * STATES[15] +
    sfi2i3 * gammac * STATES[14] - sfi2i3 * kpcb * STATES[15] + sfc4i3 *
    gammac * kpcf * (kco / koc) * STATES[12] - 4.0f * sfc4i3 * beta *
    kpcb * STATES[15] - kicalpka * STATES[83] * STATES[15] / (kicalmpka + icalcav *
    STATES[15]) + kicalpp * ppcav * STATES[114] / (kicalmpp + icalcav * STATES[114]);
  //C
  //C STATES[107] OP Ca channel variable
  //C
  RATES[107] = kcop * STATES[115] - koc * STATES[107] + sfoi1 * kpcb * STATES[112] -
    sfoi1 * gammac * STATES[107] + sfoi2 * alphap * STATES[113] - sfoi2 *
    kpcf * STATES[107] + kicalpka * STATES[83] * STATES[8] / (kicalmpka + icalcav *
    STATES[8]) - (alpha / alphap) * kicalpp * ppcav * STATES[107] / (kicalmpp +
    icalcav * STATES[107]);
  //C
  //C STATES[108] C1P Ca channel variable
  //C
  RATES[108] = beta * STATES[109] - 4.0f * alphap * STATES[108] + kicalpka * STATES[
    83] * STATES[9] / (kicalmpka + icalcav * STATES[9]) - (alphap * alphap *
    alphap * kcop / (alpha * alpha * alpha * kco)) * kicalpp *
    ppcav * STATES[108] / (kicalmpp + icalcav * STATES[108]);
  //C
  //C STATES[109] C2P Ca channel variable
  //C
  RATES[109] = 4.0f * alphap * STATES[108] - beta * STATES[109] + 2.0f * beta *
    STATES[110] - 3.0f * alphap * STATES[109] + kicalpka * STATES[83] * STATES[10] / (
    kicalmpka + icalcav * STATES[10]) - (alphap * alphap * kcop / (alpha *
    alpha * kco)) * kicalpp * ppcav * STATES[109] / (kicalmpp + icalcav *
    STATES[109]);
  //C
  //C STATES[110] C3P Ca channel variable
  //C
  RATES[110] = 3.0f * alphap * STATES[109] - 2.0f * beta * STATES[110] + 3.0f *
    beta * STATES[111] - 2.0f * alphap * STATES[110] + kicalpka * STATES[83] * STATES[
    11] / (kicalmpka + icalcav * STATES[11]) - (alphap * kcop / (alpha *
    kco)) * kicalpp * ppcav * STATES[110] / (kicalmpp + icalcav * STATES[110]);
  //C
  //C STATES[111] C4P Ca channel variable
  //C
  RATES[111] = 2.0f * alphap * STATES[110] - 3.0f * beta * STATES[111] + 4.0f *
    beta * STATES[115] - alphap * STATES[111] + 4.0f * sfc4i1 * kpcb * beta * STATES[
    112] - sfc4i1 * alphap * gammac * (kcop / koc) * STATES[111] + 4.0f *
    sfc4i2 * beta * STATES[113] - sfc4i2 * kpcf * (kcop / koc) * STATES[111] +
    4.0f * sfc4i3 * beta * kpcb * STATES[114] - sfc4i3 * gammac * kpcf * (
    kcop / koc) * STATES[111] + kicalpka * STATES[83] * STATES[12] / (kicalmpka +
    icalcav * STATES[12]) - (kcop / kco) * kicalpp * ppcav * STATES[111] / (
    kicalmpp + icalcav * STATES[111]);
  //C
  //C STATES[115] CPP Ca channel variable
  //C
  RATES[115] = alphap * STATES[111] - 4.0f * beta * STATES[115] + koc * STATES[107] -
    kcop * STATES[115] + kicalpka * STATES[83] * STATES[106] / (kicalmpka +
    icalcav * STATES[106]) - (alpha * kcop / (alphap * kco)) * kicalpp *
    ppcav * STATES[115] / (kicalmpp + icalcav * STATES[115]);
  //C
  //C STATES[112] I1P Ca channel variable
  //C
  RATES[112] = sfoi1 * gammac * STATES[107] - sfoi1 * kpcb * STATES[112] +
    sfi1i3 * alphap * STATES[114] - sfi1i3 * kpcf * STATES[112] + sfc4i1 *
    alphap * gammac * (kcop / koc) * STATES[111] - 4.0f * sfc4i1 * beta *
    kpcb * STATES[112] + kicalpka * STATES[83] * STATES[13] / (kicalmpka + icalcav *
    STATES[13]) - (alpha / alphap) * kicalpp * ppcav * STATES[112] / (
    kicalmpp + icalcav * STATES[112]);
  //C
  //C STATES[113] I2P Ca channel variable
  //C
  RATES[113] = sfoi2 * kpcf * STATES[107] - sfoi2 * alphap * STATES[113] + sfi2i3 *
    kpcb * STATES[114] - sfi2i3 * gammac * STATES[113] + sfc4i2 * kpcf * (kcop /
    koc) * STATES[111] - 4.0f * sfc4i2 * beta * STATES[113] + kicalpka * STATES[83] * STATES[
    14] / (kicalmpka + icalcav * STATES[14]) - kicalpp * ppcav * STATES[113] / (
    kicalmpp + icalcav * STATES[113]);
  //C
  //C STATES[114] I3P Ca channel variable
  //C
  RATES[114] = sfi1i3 * kpcf * STATES[112] - sfi1i3 * alphap * STATES[114] +
    sfi2i3 * gammac * STATES[113] - sfi2i3 * kpcb * STATES[114] + sfc4i3 *
    gammac * kpcf * (kcop / koc) * STATES[111] - 4.0f * sfc4i3 * beta *
    kpcb * STATES[114] + kicalpka * STATES[83] * STATES[15] / (kicalmpka + icalcav *
    STATES[15]) - kicalpp * ppcav * STATES[114] / (kicalmpp + icalcav * STATES[114]);
  //C
  //C     Partial contributions
  //C
  double gammae = kpcmax * STATES[7] / (kpchalf + STATES[7]);
  //C
  //C                ICaLecav
  //C
  //C STATES[137]   O  Ca channel variable
  //C
  RATES[137] = kco * STATES[142] - koc * STATES[137] + sfoi1 * kpcb * STATES[143] -
    sfoi1 * gammae * STATES[137] + sfoi2 * alpha * STATES[144] - sfoi2 * kpcf *
    STATES[137] - kicalpka * STATES[89] * STATES[137] / (kicalmpka + icalecav * STATES[
    137]) + (alpha / alphap) * kicalpp * pp1ecav * STATES[146] / (
    kicalmpp + icalecav * STATES[146]);
  //C
  //C STATES[138]   C1  Ca channel variable
  //C
  //C      RATES[138] = beta*STATES[139]-4.0*alpha*STATES[138]
  //C
  //C STATES[139]   C2  Ca channel variable
  //C
  RATES[139] = 4.0f * alpha * STATES[138] - beta * STATES[139] + 2.0f * beta * STATES[
    140] - 3.0f * alpha * STATES[139] - kicalpka * STATES[89] * STATES[139] / (
    kicalmpka + icalecav * STATES[139]) + (alphap * alphap * kcop / (
    alpha * alpha * kco)) * kicalpp * pp1ecav * STATES[148] / (kicalmpp +
    icalecav * STATES[148]);
  //C
  //C STATES[140]   C3  Ca channel variable
  //C
  RATES[140] = 3.0f * alpha * STATES[139] - 2.0f * beta * STATES[140] + 3.0f *
    beta * STATES[141] - 2.0f * alpha * STATES[140] - kicalpka * STATES[89] * STATES[
    140] / (kicalmpka + icalecav * STATES[140]) + (alphap * kcop / (
    alpha * kco)) * kicalpp * pp1ecav * STATES[149] / (kicalmpp +
    icalecav * STATES[149]);
  //C
  //C STATES[141]   C4  Ca channel variable
  //C
  RATES[141] = 2.0f * alpha * STATES[140] - 3.0f * beta * STATES[141] + 4.0f *
    beta * STATES[142] - alpha * STATES[141] + 4.0f * sfc4i1 * kpcb * beta * STATES[
    143] - sfc4i1 * alpha * gammae * (kco / koc) * STATES[141] + 4.0f *
    sfc4i2 * beta * STATES[144] - sfc4i2 * kpcf * (kco / koc) * STATES[141] +
    4.0f * sfc4i3 * beta * kpcb * STATES[145] - sfc4i3 * gammae * kpcf * (
    kco / koc) * STATES[141] - kicalpka * STATES[89] * STATES[141] / (kicalmpka +
    icalecav * STATES[141]) + (kcop / kco) * kicalpp * pp1ecav * STATES[150] / (
    kicalmpp + icalecav * STATES[150]);
  //C
  //C STATES[142] CP Ca channel variable
  //C
  RATES[142] = alpha * STATES[141] - 4.0f * beta * STATES[142] + koc * STATES[137] -
    kco * STATES[142] - kicalpka * STATES[89] * STATES[142] / (kicalmpka +
    icalecav * STATES[142]) + (alpha * kcop / (alphap * kco)) * kicalpp *
    pp1ecav * STATES[151] / (kicalmpp + icalecav * STATES[151]);
  //C
  //C STATES[143]   I1  Ca channel variable
  //C
  RATES[143] = sfoi1 * gammae * STATES[137] - sfoi1 * kpcb * STATES[143] +
    sfi1i3 * alpha * STATES[145] - sfi1i3 * kpcf * STATES[143] + sfc4i1 *
    alpha * gammae * (kco / koc) * STATES[141] - 4.0f * sfc4i1 * beta *
    kpcb * STATES[143] - kicalpka * STATES[89] * STATES[143] / (kicalmpka +
    icalecav * STATES[143]) + (alpha / alphap) * kicalpp * pp1ecav * STATES[
    152] / (kicalmpp + icalecav * STATES[152]);
  //C
  //C STATES[144]   I2  Ca channel variable
  //C
  RATES[144] = sfoi2 * kpcf * STATES[137] - sfoi2 * alpha * STATES[144] +
    sfi2i3 * kpcb * STATES[145] - sfi2i3 * gammae * STATES[144] + sfc4i2 *
    kpcf * (kco / koc) * STATES[141] - 4.0f * sfc4i2 * beta * STATES[144] -
    kicalpka * STATES[89] * STATES[144] / (kicalmpka + icalecav * STATES[144]) +
    kicalpp * pp1ecav * STATES[153] / (kicalmpp + icalecav * STATES[153]);
  //C
  //C STATES[145]   I3  Ca channel variable
  //C
  RATES[145] = sfi1i3 * kpcf * STATES[143] - sfi1i3 * alpha * STATES[145] +
    sfi2i3 * gammae * STATES[144] - sfi2i3 * kpcb * STATES[145] + sfc4i3 *
    gammae * kpcf * (kco / koc) * STATES[141] - 4.0f * sfc4i3 * beta *
    kpcb * STATES[145] - kicalpka * STATES[89] * STATES[145] / (kicalmpka +
    icalecav * STATES[145]) + kicalpp * pp1ecav * STATES[154] / (kicalmpp +
    icalecav * STATES[154]);
  //C
  //C STATES[146] OP Ca channel variable
  //C
  RATES[146] = kcop * STATES[151] - koc * STATES[146] + sfoi1 * kpcb * STATES[152] -
    sfoi1 * gammae * STATES[146] + sfoi2 * alphap * STATES[153] - sfoi2 * kpcf *
    STATES[146] + kicalpka * STATES[89] * STATES[137] / (kicalmpka + icalecav * STATES[
    137]) - (alpha / alphap) * kicalpp * pp1ecav * STATES[146] / (
    kicalmpp + icalecav * STATES[146]);
  //C
  //C STATES[147] C1P Ca channel variable
  //C
  RATES[147] = beta * STATES[148] - 4.0f * alphap * STATES[147] + kicalpka * STATES[
    89] * STATES[138] / (kicalmpka + icalecav * STATES[138]) - (alphap *
    alphap * alphap * kcop / (alpha * alpha * alpha * kco)) *
    kicalpp * pp1ecav * STATES[147] / (kicalmpp + icalecav * STATES[147]);
  //C
  //C STATES[148] C2P Ca channel variable
  //C
  RATES[148] = 4.0f * alphap * STATES[147] - beta * STATES[148] + 2.0f * beta *
    STATES[149] - 3.0f * alphap * STATES[148] + kicalpka * STATES[89] * STATES[139] / (
    kicalmpka + icalecav * STATES[139]) - (alphap * alphap * kcop / (
    alpha * alpha * kco)) * kicalpp * pp1ecav * STATES[148] / (kicalmpp +
    icalecav * STATES[148]);
  //C
  //C STATES[149] C3P Ca channel variable
  //C
  RATES[149] = 3.0f * alphap * STATES[148] - 2.0f * beta * STATES[149] + 3.0f *
    beta * STATES[150] - 2.0f * alphap * STATES[149] + kicalpka * STATES[89] * STATES[
    140] / (kicalmpka + icalecav * STATES[140]) - (alphap * kcop / (
    alpha * kco)) * kicalpp * pp1ecav * STATES[149] / (kicalmpp +
    icalecav * STATES[149]);
  //C
  //C STATES[150] C4P Ca channel variable
  //C
  RATES[150] = 2.0f * alphap * STATES[149] - 3.0f * beta * STATES[150] + 4.0f *
    beta * STATES[151] - alphap * STATES[150] + 4.0f * sfc4i1 * kpcb * beta * STATES[
    152] - sfc4i1 * alphap * gammae * (kcop / koc) * STATES[150] + 4.0f *
    sfc4i2 * beta * STATES[153] - sfc4i2 * kpcf * (kcop / koc) * STATES[150] +
    4.0f * sfc4i3 * beta * kpcb * STATES[154] - sfc4i3 * gammae * kpcf * (
    kcop / koc) * STATES[150] + kicalpka * STATES[89] * STATES[141] / (kicalmpka +
    icalecav * STATES[141]) - (kcop / kco) * kicalpp * pp1ecav * STATES[150] / (
    kicalmpp + icalecav * STATES[150]);
  //C
  //C STATES[151] CPP Ca channel variable
  //C
  RATES[151] = alphap * STATES[150] - 4.0f * beta * STATES[151] + koc * STATES[146] -
    kcop * STATES[151] + kicalpka * STATES[89] * STATES[142] / (kicalmpka +
    icalecav * STATES[142]) - (alpha * kcop / (alphap * kco)) * kicalpp *
    pp1ecav * STATES[151] / (kicalmpp + icalecav * STATES[151]);
  //C
  //C STATES[152] I1P Ca channel variable
  //C
  RATES[152] = sfoi1 * gammae * STATES[146] - sfoi1 * kpcb * STATES[152] +
    sfi1i3 * alphap * STATES[154] - sfi1i3 * kpcf * STATES[152] + sfc4i1 *
    alphap * gammae * (kcop / koc) * STATES[150] - 4.0f * sfc4i1 * beta *
    kpcb * STATES[152] + kicalpka * STATES[89] * STATES[143] / (kicalmpka +
    icalecav * STATES[143]) - (alpha / alphap) * kicalpp * pp1ecav * STATES[
    152] / (kicalmpp + icalecav * STATES[152]);
  //C
  //C STATES[153] I2P Ca channel variable
  //C
  RATES[153] = sfoi2 * kpcf * STATES[146] - sfoi2 * alphap * STATES[153] +
    sfi2i3 * kpcb * STATES[154] - sfi2i3 * gammae * STATES[153] + sfc4i2 *
    kpcf * (kcop / koc) * STATES[150] - 4.0f * sfc4i2 * beta * STATES[153] +
    kicalpka * STATES[89] * STATES[144] / (kicalmpka + icalecav * STATES[144]) -
    kicalpp * pp1ecav * STATES[153] / (kicalmpp + icalecav * STATES[153]);
  //C
  //C STATES[154] I3P Ca channel variable
  //C
  RATES[154] = sfi1i3 * kpcf * STATES[152] - sfi1i3 * alphap * STATES[154] +
    sfi2i3 * gammae * STATES[153] - sfi2i3 * kpcb * STATES[154] + sfc4i3 *
    gammae * kpcf * (kcop / koc) * STATES[150] - 4.0f * sfc4i3 * beta *
    kpcb * STATES[154] + kicalpka * STATES[89] * STATES[145] / (kicalmpka +
    icalecav * STATES[145]) - kicalpp * pp1ecav * STATES[154] / (kicalmpp +
    icalecav * STATES[154]);
  //C
  //C STATES[16]-STATES[19]  RyR channel states
  //C
  double temp5 = STATES[7] * STATES[7] * STATES[7];
  double temp6 = STATES[7] * STATES[7] * STATES[7] * STATES[7];
  //C
  //C STATES[16] RyR C1 channel variable
  //C
  RATES[16] = -kap * temp6 * STATES[16] + kam * STATES[18] - kryrpka * STATES[89] * STATES[
    16] / (kryrmpka + ryrecav * STATES[16]) + kryrpp * pp1ecav * STATES[116] / (
    kryrmpp + ryrecav * STATES[116]);
  //C
  //C STATES[17] RyR C2 channel variable
  //C
  RATES[17] = kcp * STATES[18] - kcm * STATES[17] - kryrpka * STATES[89] * STATES[17] * f_ryr / (
    kryrmpka + ryrecav * STATES[17]) + ((kap * kamp * kcp * kcmp) / (kapp *
    kam * kcpp * kcm)) * f_ryr * kryrpp * pp1ecav * STATES[117] / (kryrmpp +
    ryrecav * STATES[117]);
  //C
  //C STATES[18] RyR O1 channel variable
  //C
  RATES[18] = kap * temp6 * STATES[16] - kam * STATES[18] - kbp * temp5 * STATES[
    18] + kbm * STATES[19] - kcp * STATES[18] + kcm * STATES[17] - kryrpka * STATES[89] *
    STATES[18] * f_ryr / (kryrmpka + ryrecav * STATES[18]) + ((kap * kamp) / (
    kapp * kam)) * f_ryr * kryrpp * pp1ecav * STATES[118] / (kryrmpp +
    ryrecav * STATES[118]);
  //C
  //C STATES[19] RyR O2 channel variable
  //C
  RATES[19] = kbp * temp5 * STATES[18] - kbm * STATES[19] - kryrpka * STATES[89] * STATES[
    19] * f_ryr / (kryrmpka + ryrecav * STATES[19]) + ((kap * kamp * kbp *
    kbmp) / (kapp * kam * kbpp * kbm)) * f_ryr * kryrpp * pp1ecav * STATES[
    119] / (kryrmpp + ryrecav * STATES[119]);
  //C
  //C STATES[116] RyR C1P channel variable
  //C
  RATES[116] = -kapp * temp6 * STATES[116] + kamp * STATES[118] + kryrpka * STATES[
    89] * STATES[16] / (kryrmpka + ryrecav * STATES[16]) - kryrpp * pp1ecav * STATES[
    116] / (kryrmpp + ryrecav * STATES[116]);
  //C
  //C STATES[117] RyR C2P channel variable
  //C
  RATES[117] = kcpp * STATES[118] - kcmp * STATES[117] + kryrpka * STATES[89] * STATES[
    17] * f_ryr / (kryrmpka + ryrecav * STATES[17]) - ((kap * kamp * kcp *
    kcmp) / (kapp * kam * kcpp * kcm)) * f_ryr * kryrpp * pp1ecav * STATES[
    117] / (kryrmpp + ryrecav * STATES[117]);
  //C
  //C STATES[118] RyR O1P channel variable
  //C
  RATES[118] = kapp * temp6 * STATES[116] - kamp * STATES[118] - kbpp * temp5 * STATES[
    118] + kbmp * STATES[119] - kcpp * STATES[118] + kcmp * STATES[117] + kryrpka * STATES[
    89] * STATES[18] * f_ryr / (kryrmpka + ryrecav * STATES[18]) - ((kap *
    kamp) / (kapp * kam)) * f_ryr * kryrpp * pp1ecav * STATES[118] / (
    kryrmpp + ryrecav * STATES[118]);
  //C
  //C STATES[119] RyR O2P channel variable
  //C
  RATES[119] = kbpp * temp5 * STATES[118] - kbmp * STATES[119] + kryrpka * STATES[
    89] * STATES[19] * f_ryr / (kryrmpka + ryrecav * STATES[19]) - ((kap *
    kamp * kbp * kbmp) / (kapp * kam * kbpp * kbm)) * f_ryr *
    kryrpp * pp1ecav * STATES[119] / (kryrmpp + ryrecav * STATES[119]);
  //C
  //C  STATES[20]-STATES[22],STATES[37]-STATES[42] Na fast current (Clancy-Rudy, Circulation, 2002]
  //C
  double va = STATES[1] - 2.5f;
  double vi = STATES[1] + 7.0f;
  double alp11 = 3.802f / (0.1027f * exp(-va / 17.0f) + 0.20f *
    exp(-va / 150.f));
  double alp12 = 3.802f / (0.1027f * exp(-va / 15.0f) + 0.23f *
    exp(-va / 150.f));
  double alp13 = 3.802f / (0.1027f * exp(-va / 12.0f) + 0.25f *
    exp(-va / 150.f));
  double bet11 = 0.1917f * exp(-va / 20.3f);
  double bet12 = 0.20f * exp(-(va - 5.0f) / 20.3f);
  double bet13 = 0.22f * exp(-(va - 10.0f) / 20.3f);
  double alp3 = 7.0e-7f * exp(-vi / 7.7f);
  double bet3 = (0.0084f + 0.00002f * vi);
  double alp2 = 1.0f / (0.188495f * exp(-vi / 16.6f) + 0.393956f);
  double bet2 = alp13 * alp2 * alp3 / (bet13 * bet3);
  double alp4 = alp2 / 100.f;
  double bet4 = alp3;
  double alp5 = alp2 / 9.5e4f;
  double bet5 = alp3 / 50.0f;
  double vap = STATES[1] - 2.5f;
  //double vip = STATES[1] + 7.0f;
  double alp11p = 3.802f / (0.1027f * exp(-vap / 17.0f) + 0.20f *
    exp(-vap / 150.f));
  double alp12p = 3.802f / (0.1027f * exp(-vap / 15.0f) + 0.23f *
    exp(-vap / 150.f));
  double alp13p = 3.802f / (0.1027f * exp(-vap / 12.0f) + 0.25f *
    exp(-vap / 150.f));
  double bet2p = alp13p * alp2 * alp3 / (bet13 * bet3);
  //C
  //C STATES[20] INa CNa3 channel variable
  //C
  //C      RATES[20] = bet11*STATES[21] - alp11*STATES[20]
  //C     *         + alp3*STATES[42]  - bet3*STATES[20]
  //C     * - kINAPKA*STATES[83]*STATES[20]/(KINamPKA+STATES[20])
  //C     * + (alp11p*alp12p/(alp11*alp12])*kINaPP*PPcav*STATES[120]
  //C     * /(KINamPP+STATES[120])
  //C
  //C STATES[21] INa CNa2 channel variable
  //C
  RATES[21] = alp11 * STATES[20] - bet11 * STATES[21] + bet12 * STATES[22] - alp12 * STATES[
    21] + alp3 * STATES[41] - bet3 * STATES[21] - kinapka * STATES[83] * STATES[21] / (
    kinampka + STATES[21]) + (alp12p / alp12) * kinapp * ppcav * STATES[121] / (
    kinampp + STATES[121]);
  //C
  //C STATES[22] INa CNa1 channel variable
  //C
  RATES[22] = alp12 * STATES[21] - bet12 * STATES[22] + bet13 * STATES[37] - alp13 *
    STATES[22] + alp3 * STATES[38] - bet3 * STATES[22] - kinapka * STATES[83] * STATES[22] / (
    kinampka + STATES[22]) + kinapp * ppcav * STATES[122] / (kinampp + STATES[122]);
  //C
  //C STATES[37] INa ONa channel variable
  //C
  RATES[37] = alp13 * STATES[22] - bet13 * STATES[37] + bet2 * STATES[38] - alp2 * STATES[
    37] - kinapka * STATES[83] * STATES[37] / (kinampka + STATES[37]) + (alp13 /
    alp13p) * kinapp * ppcav * STATES[123] / (kinampp + STATES[123]);
  //C
  //C STATES[38] INa IFNa channel variable
  //C
  RATES[38] = alp2 * STATES[37] - bet2 * STATES[38] + bet3 * STATES[22] - alp3 * STATES[38] +
    bet4 * STATES[39] - alp4 * STATES[38] + alp12 * STATES[41] - bet12 * STATES[38] -
    kinapka * STATES[83] * STATES[38] / (kinampka + STATES[38]) + kinapp * ppcav * STATES[
    124] / (kinampp + STATES[124]);
  //C
  //C STATES[39] INa I1Na channel variable
  //C
  RATES[39] = alp4 * STATES[38] - bet4 * STATES[39] + bet5 * STATES[40] - alp5 * STATES[
    39] - kinapka * STATES[83] * STATES[39] / (kinampka + STATES[39]) + kinapp *
    ppcav * STATES[125] / (kinampp + STATES[125]);
  //C
  //C STATES[40] INa I2Na channel variable
  //C
  RATES[40] = alp5 * STATES[39] - bet5 * STATES[40] - kinapka * STATES[83] * STATES[40] / (
    kinampka + STATES[40]) + kinapp * ppcav * STATES[126] / (kinampp + STATES[126]);
  //C
  //C STATES[41] INa ICNa2 channel variable
  //C
  RATES[41] = alp11 * STATES[42] - bet11 * STATES[41] + bet12 * STATES[38] - alp12 * STATES[
    41] + bet3 * STATES[21] - alp3 * STATES[41] - kinapka * STATES[83] * STATES[41] / (
    kinampka + STATES[41]) + (alp12p / alp12) * kinapp * ppcav * STATES[127] / (
    kinampp + STATES[127]);
  //C
  //C STATES[42] INa ICNa3 channel variable
  //C
  RATES[42] = bet11 * STATES[41] - alp11 * STATES[42] + bet3 * STATES[20] - alp3 * STATES[
    42] - kinapka * STATES[83] * STATES[42] / (kinampka + STATES[42]) + (alp11p *
    alp12p / (alp11 * alp12)) * kinapp * ppcav * STATES[128] / (kinampp +
    STATES[128]);
  //C
  //C STATES[120] INa CNa3p channel variable
  //C
  RATES[120] = bet11 * STATES[121] - alp11p * STATES[120] + alp3 * STATES[128] -
    bet3 * STATES[120] + kinapka * STATES[83] * STATES[20] / (kinampka + STATES[20]) - (
    alp11p * alp12p / (alp11 * alp12)) * kinapp * ppcav * STATES[120] / (
    kinampp + STATES[120]);
  //C
  //C STATES[121] INa CNa2p channel variable
  //C
  RATES[121] = alp11p * STATES[120] - bet11 * STATES[121] + bet12 * STATES[122] - alp12p *
    STATES[121] + alp3 * STATES[127] - bet3 * STATES[121] + kinapka * STATES[83] * STATES[21] / (
    kinampka + STATES[21]) - (alp12p / alp12) * kinapp * ppcav * STATES[121] / (
    kinampp + STATES[121]);
  //C
  //C STATES[122] INa CNa1p channel variable
  //C
  RATES[122] = alp12p * STATES[121] - bet12 * STATES[122] + bet13 * STATES[123] -
    alp13p * STATES[122] + alp3 * STATES[124] - bet3 * STATES[122] + kinapka * STATES[
    83] * STATES[22] / (kinampka + STATES[22]) - kinapp * ppcav * STATES[122] / (
    kinampp + STATES[122]);
  //C
  //C STATES[123] INa ONap channel variable
  //C
  RATES[123] = alp13p * STATES[122] - bet13 * STATES[123] + bet2p * STATES[124] -
    alp2 * STATES[123] + kinapka * STATES[83] * STATES[37] / (kinampka + STATES[37]) - (
    alp13 / alp13p) * kinapp * ppcav * STATES[123] / (kinampp + STATES[123]);
  //C
  //C STATES[124] INa IFNap channel variable
  //C
  RATES[124] = alp2 * STATES[123] - bet2p * STATES[124] + bet3 * STATES[122] - alp3 *
    STATES[124] + bet4 * STATES[125] - alp4 * STATES[124] + alp12p * STATES[127] -
    bet12 * STATES[124] + kinapka * STATES[83] * STATES[38] / (kinampka + STATES[38]) -
    kinapp * ppcav * STATES[124] / (kinampp + STATES[124]);
  //C
  //C STATES[125] INa I1Nap channel variable
  //C
  RATES[125] = alp4 * STATES[124] - bet4 * STATES[125] + bet5 * STATES[126] - alp5 *
    STATES[125] + kinapka * STATES[83] * STATES[39] / (kinampka + STATES[39]) - kinapp *
    ppcav * STATES[125] / (kinampp + STATES[125]);
  //C
  //C STATES[126] INa I2Nap channel variable
  //C
  RATES[126] = alp5 * STATES[125] - bet5 * STATES[126] + kinapka * STATES[83] * STATES[
    40] / (kinampka + STATES[40]) - kinapp * ppcav * STATES[126] / (kinampp + STATES[
    126]);
  //C
  //C STATES[127] INa ICNa2p channel variable
  //C
  RATES[127] = alp11p * STATES[128] - bet11 * STATES[127] + bet12 * STATES[124] - alp12p *
    STATES[127] + bet3 * STATES[121] - alp3 * STATES[127] + kinapka * STATES[83] * STATES[41] / (
    kinampka + STATES[41]) - (alp12p / alp12) * kinapp * ppcav * STATES[127] / (
    kinampp + STATES[127]);
  //C
  //C STATES[128] INa ICNa3p channel variable
  //C
  RATES[128] = bet11 * STATES[127] - alp11p * STATES[128] + bet3 * STATES[120] -
    alp3 * STATES[128] + kinapka * STATES[83] * STATES[42] / (kinampka + STATES[42]) - (
    alp11p * alp12p / (alp11 * alp12)) * kinapp * ppcav * STATES[128] / (
    kinampp + STATES[128]);
  //C
  //C  STATES[23]  Na intracellular concentration
  //C
  RATES[23] = -1.0f * (ALGEBRAIC[20]/*ina*/ +ALGEBRAIC[21]/*inab*/+ 3.0f * (ALGEBRAIC[18]/*inaca*/ + ALGEBRAIC[25]/*inak*/)) *ALGEBRAIC[6]/*acap*// (ALGEBRAIC[7]/*vmyo*/ * ALGEBRAIC[12]/*f*/);
  //C
  //C  STATES[24]  K  intracellular concentration
  //C
  RATES[24] = -1.0f * (ALGEBRAIC[26]/*ikto*/ +ALGEBRAIC[27]/*ik1*/+ALGEBRAIC[28]/*iks*/+ALGEBRAIC[29]/*ikur*/+ALGEBRAIC[30]/*ikr*/+ALGEBRAIC[112]/*icak*/- 2.0f *
   ALGEBRAIC[25]/*inak*/- ALGEBRAIC[113]/*istim*/) *ALGEBRAIC[6]/*acap*// (ALGEBRAIC[7]/*vmyo*/ * ALGEBRAIC[12]/*f*/);
  //C
  //C  STATES[25],STATES[26]  ato and ito gating variables for Ikto
  //C
  double alp25 = 0.04516f * exp(0.03577f * (STATES[1] + 33.0f)) * 4.0f;
  double bet25 = 0.0989f * exp(-0.06237f * (STATES[1] + 33.0f)) * 4.0f;
  double temp7 = 0.0019f * exp((STATES[1] + 15.5f) / (-1.0f * dito));
  double temp8 = 0.067083f * exp((STATES[1] + 15.5f + 20.0f) / (-1.0f * dito));
  double alp26 = 0.08f * temp7 / (1.0f + temp8);
  double temp9 = 0.0019f * exp((STATES[1] + 15.5f + 20.0f) / dito);
  double temp10 = 0.051335f * exp((STATES[1] + 15.5f + 20.0f) / dito);
  double bet26 = 0.5f * temp9 / (1.0f + temp10);
  //C
  RATES[25] = alp25 * (1.0f - STATES[25]) - bet25 * STATES[25];
  RATES[26] = alp26 * (1.0f - STATES[26]) - bet26 * STATES[26];
  //C
  //C  STATES[135],STATES[136] ato and ito gating variables for Ikto,phosphorylated
  //C
  double alp25p = 0.04516f * exp(0.03577f * (STATES[1] + 30.0f - 13.0f)) * 4.0f;
  double bet25p = 0.0989f * exp(-0.06237f * (STATES[1] + 30.0f - 13.0f)) * 4.0f;
  double temp7p = 0.0019f * exp((STATES[1] + 13.5f - 6.0f) / (-1.0f * dito));
  double temp8p = 0.067083f * exp((STATES[1] + 13.5f + 20.0f - 6.0f) /
    (-1.0f * dito));
  double alp26p = 0.08f * temp7p / (1.0f + temp8p);
  double temp9p = 0.0019f * exp((STATES[1] + 13.5f + 20.0f - 6.0f) / dito);
  double temp10p = 0.051335f * exp((STATES[1] + 13.5f + 20.0f - 6.0f) / dito);
  double bet26p = 0.5f * temp9p / (1.0f + temp10p);
  //C
  RATES[135] = alp25p * (1.0f - STATES[135]) - bet25p * STATES[135];
  RATES[136] = alp26p * (1.0f - STATES[136]) - bet26p * STATES[136];
  //C
  //C  STATES[27]  nks gating variable for Iks
  //C
  double temp11 = 0.00001444f * (STATES[1] + 26.5f);
  double alp27 = temp11 / (1.0f - exp(-0.128f * (STATES[1] + 26.5f))) / 3.f;
  double bet27 = 0.000286f * exp(-0.038f * (STATES[1] + 26.5f)) / 3.f;
  //C
  RATES[27] = alp27 * (1.0f - STATES[27]) - bet27 * STATES[27];
  //C
  //C  STATES[28], STATES[29]  aur and iur gating variables for Ikur1
  //C
  double ass = 1.0f / (1.0f + exp((-22.5f - STATES[1]) / 7.7f));
  double iss1 = 1.0f / (1.0f + exp((45.2f + STATES[1]) / 5.7f));
  //C      taua1 = [0.493*exp(-0.0629*STATES[1])+2.058]
  double taua1 = 6.1f / (exp(0.0629f * (STATES[1] + 40.0f)) + exp(
    -0.0629f * (STATES[1] + 40.0f))) + 2.058f;
  double taui1 = 270.f + 1050 / (1.0f + exp((45.2f + STATES[1]) / 5.7f));
  //C
  RATES[28] = (ass - STATES[28]) / taua1;
  RATES[29] = (iss1 - STATES[29]) / taui1;
  //C
  //C  STATES[35], STATES[36]  aur and iur gating variables for Ikur2
  //C
  double iss2 = 1.0f / (1.0f + exp((45.2f + STATES[1]) / 5.7f));
  //C      taua2 = [0.493*exp(-0.0629*STATES[1])+2.058]
  double taua2 = 6.1f / (exp(0.0629f * (STATES[1] + 40.0f)) + exp(
    -0.0629f * (STATES[1] + 40.0f))) + 2.058f;
  double taui2 = 803.f - 18.f / (1.0f + exp((45.2f + STATES[1]) / 5.7f));
  //C
  RATES[35] = (ass - STATES[35]) / taua2;
  RATES[36] = (iss2 - STATES[36]) / taui2;
  RATES[130] = (ass - STATES[130]) / taua2;
  RATES[131] = (iss2 - STATES[131]) / taui2;
  RATES[132] = kikurpp * pp1ecav * (1 - STATES[132]) / (kikurmpp + (1 - STATES[
    132])) - kikurpka * STATES[89] * STATES[132] / (kikurmpka + STATES[132]);
  //C
  //C  STATES[43] aur gating variable for Ikur3
  //C
  //C      taua3 = [39.3*exp(-0.0862*STATES[1])+13.17]
  double taua3 = 1235.5f / (exp(0.0862f * (STATES[1] + 40.0f)) +
    exp(-0.0862f * (STATES[1] + 40.0f))) + 13.17f;
  RATES[43] = (ass - STATES[43]) / taua3;
  //C
  //C  STATES[184], STATES[185] aur and iur gating variables for Ikur4
  //C
  double taui4 = 5334.f - 4912.f / (1.0f + exp((45.2f + STATES[1]) / 5.7f));
  RATES[184] = (ass - STATES[184]) / taua2;
  RATES[185] = (iss2 - STATES[185]) / taui4;
  //C
  //C STATES[30]-STATES[34] HERG channel state variables
  //C
  double ala0 = 0.022348f * exp(0.01176f * STATES[1]);
  double bea0 = 0.047002f * exp(-0.0631f * STATES[1]);
  double ala1 = 0.013733f * exp(0.038198f * STATES[1]);
  double bea1 = 0.0000689f * exp(-0.04178f * STATES[1]);
  double ali = 0.090821f * exp(0.023391f * (STATES[1] + 5.0f));
  double bei = 0.006497f * exp(-0.03268f * (STATES[1] + 5.0f));
  //C
  RATES[30] = bea0 * STATES[31] - ala0 * STATES[30];
  RATES[31] = ala0 * STATES[30] - bea0 * STATES[31] + kb * STATES[32] - kf * STATES[31];
  RATES[32] = kf * STATES[31] - kb * STATES[32] + bea1 * STATES[33] - ala1 * STATES[32];
  RATES[33] = ala1 * STATES[32] - bea1 * STATES[33] + bei * STATES[34] - ali * STATES[33];
  RATES[34] = ali * STATES[33] - bei * STATES[34];
  //C
  //C   STATES[44] Pryr Ryanodine receptor modulation factor
  //C
  RATES[44] = -t1 * STATES[44] + t2 *ALGEBRAIC[16]/*icas*/* exp((STATES[1] + 5.0f) * (STATES[
    1] + 5.0f) / (-648.f));
  //C
  //C   Caveolae domain
  //C
  //C   STATES[45] beta1 tot concentration phosphorylated by PKA caveolae
  //C
  RATES[45] = kpkap * STATES[83] * rb1cavnptot - kpkam * STATES[45];
  //C
  //C   STATES[46] beta1 tot concentration phosphorylated by BARK caveolae
  //C
  RATES[46] = kbarkp * (lrb1cavnp + lrb1gscavnp) - kbarkm * STATES[46];
  //C
  //C   STATES[155] beta2 tot concentration phosphorylated by PKA caveolae
  //C
  RATES[155] = kpkap * STATES[83] * rb2cavnptot - kpkam * STATES[155];
  //C
  //C   STATES[156] beta2 tot concentration phosphorylated by BARK caveolae
  //C
  RATES[156] = kbarkp * (lrb2cavnp + lrb2gscavnp) - kbarkm * STATES[156];
  //C
  //C   STATES[47] Gs-alpha with GTP caveolae
  //C
  RATES[47] = kact2gs * rb1gscavnp + factgsgi * kact2gs * rb2gscavnp +
    kact1gs * lrb1gscavnp + factgsgi * kact1gs * lrb2gscavnp -
    khydgs * STATES[47];
  //C
  //C   STATES[48]  G Beta-Gamma caveolae
  //C
  RATES[48] = kact2gs * rb1gscavnp + factgsgi * kact2gs * rb2gscavnp +
    kact1gs * lrb1gscavnp + factgsgi * kact1gs * lrb2gscavnp + kact2gi *
   ALGEBRAIC[95]/*rb2gicavpka*/+ kact1gi *ALGEBRAIC[96]/*lrb2gicavpka*/- kreasgs * STATES[48] * STATES[49] -
    kreasgi * STATES[48] * STATES[158];
  //C
  //C   STATES[49]  Gs-alpha with GDP caveolae
  //C
  RATES[49] = khydgs * STATES[47] - kreasgs * STATES[48] * STATES[49];
  //C
  //C   STATES[157] Gi-alpha with GTP caveolae
  //C
  RATES[157] = kact2gi *ALGEBRAIC[95]/*rb2gicavpka*/+ kact1gi *ALGEBRAIC[96]/*lrb2gicavpka*/- khydgi * STATES[157];
  //C
  //C   STATES[158]  Gi-alpha with GDP caveolae
  //C
  RATES[158] = khydgi * STATES[157] - kreasgi * STATES[48] * STATES[158];
  //C
  //C   STATES[50] beta1 tot concentration phosphorylated by PKA extracaveolae
  //C
  RATES[50] = kpkap * STATES[89] * rb1ecavnptot - kpkam * STATES[50];
  //C
  //C   STATES[51] beta1 tot concentration phosphorylated by BARK extracaveolae
  //C
  RATES[51] = kbarkp * (lrb1ecavnp + lrb1gsecavnp) - kbarkm * STATES[51];
  //C
  //C   STATES[159] beta2 tot concentration phosphorylated by PKA extracaveolae
  //C
  RATES[159] = kpkap * STATES[89] * rb2ecavnptot - kpkam * STATES[159];
  //C
  //C   STATES[160] beta2 tot concentration phosphorylated by BARK extracaveolae
  //C
  RATES[160] = kbarkp * (lrb2ecavnp + lrb2gsecavnp) - kbarkm * STATES[160];
  //C
  //C   STATES[52] Gs-alpha with GTP extracaveolae
  //C
  RATES[52] = kact2gs * rb1gsecavnp + factgsgi * kact2gs * rb2gsecavnp +
    kact1gs * lrb1gsecavnp + factgsgi * kact1gs * lrb2gsecavnp -
    khydgs * STATES[52];
  //C
  //C   STATES[53]  G Beta-Gamma extracaveolae
  //C
  RATES[53] = kact2gs * rb1gsecavnp + factgsgi * kact2gs *
    rb2gsecavnp + kact1gs * lrb1gsecavnp + factgsgi * kact1gs *
    lrb2gsecavnp + kact2gi *ALGEBRAIC[102]/*rb2giecavpka*/+ kact1gi *ALGEBRAIC[103]/*lrb2giecavpka*/-
    kreasgs * STATES[53] * STATES[54] - kreasgi * STATES[53] * STATES[162];
  //C
  //C   STATES[54]  Gs-alpha with GDP extracaveolae
  //C
  RATES[54] = khydgs * STATES[52] - kreasgs * STATES[53] * STATES[54];
  //C
  //C   STATES[161] Gi-alpha with GTP extracaveolae
  //C
  RATES[161] = kact2gi *ALGEBRAIC[102]/*rb2giecavpka*/+ kact1gi *ALGEBRAIC[103]/*lrb2giecavpka*/-
    khydgi * STATES[161];
  //C
  //C   STATES[162]  Gi-alpha with GDP extracaveolae
  //C
  RATES[162] = khydgi * STATES[161] - kreasgi * STATES[53] * STATES[162];
  //C
  //C   STATES[55] beta1 tot concentration phosphorylated by PKA cytosol
  //C
  RATES[55] = kpkap * STATES[95] * rb1cytnptot - kpkam * STATES[55];
  //C
  //C   STATES[56] beta1 tot concentration phosphorylated by BARK
  //C
  RATES[56] = kbarkp * (lrb1cytnp + lrb1gscytnp) - kbarkm * STATES[56];
  //C
  //C   STATES[57] Gs-alpha with GTP
  //C
  RATES[57] = kact2gs * rb1gscytnp + kact1gs * lrb1gscytnp - khydgs * STATES[57];
  //C
  //C   STATES[58]  Gs Beta-Gamma
  //C
  RATES[58] = kact2gs * rb1gscytnp + kact1gs * lrb1gscytnp - kreasgs *
    STATES[58] * STATES[59];
  //C
  //C   STATES[59]  Gs-alpha with GDP
  //C
  RATES[59] = khydgs * STATES[57] - kreasgs * STATES[58] * STATES[59];
  //C
  //C   STATES[60] cAMP from AC56 in ceveolae
  //C
  RATES[60] = kcavac56 * ac56cav * atp / (kmatp + atp);
  //C
  //C   STATES[61] cAMP from AC47 in extracaveolae
  //C
  RATES[61] = kecavac47 * ac47ecav * atp / (kmatp + atp);
  //C
  //C   STATES[62] cAMP from AC56 in cytosol
  //C
  RATES[62] = kcytac56 * ac56cyt * atp / (kmatp + atp);
  //C
  //C   STATES[63] cAMP from AC47 in cytosol
  //C
  RATES[63] = kcytac47 * ac47cyt * atp / (kmatp + atp);
  //C
  //C   STATES[64] PDE3 caveolar phosphorylated
  //C
  RATES[64] = kfpdep * STATES[83] * (ALGEBRAIC[55]/*pde3cavtot*/ - STATES[64]) - kbpdep * STATES[64];
  //C
  //C   STATES[65] PDE4 caveolar phosphorylated
  //C
  RATES[65] = kfpdep * STATES[83] * (ALGEBRAIC[56]/*pde4cavtot*/ - STATES[65]) - kbpdep * STATES[65];
  //C
  //C   STATES[66] cAMP change due to PDE2 caveolar domain
  //C
  RATES[66] = (kpde2 *ALGEBRAIC[54]/*pde2cavtot*/* STATES[97]) / (kmpde2 + STATES[97]);
  //C
  //C   STATES[67] cAMP change due to PDE3 caveolar domain
  //C
  RATES[67] = (kpde3 * (ALGEBRAIC[55]/*pde3cavtot*/ - STATES[64]) * STATES[97] + deltakpde34 *
    kpde3 * STATES[64] * STATES[97]) / (kmpde3 + STATES[97]);
  //C
  //C   STATES[68] cAMP change due to PDE4 caveolar domain
  //C
  RATES[68] = (kpde4 * (ALGEBRAIC[56]/*pde4cavtot*/ - STATES[65]) * STATES[97] + deltakpde34 *
    kpde4 * STATES[65] * STATES[97]) / (kmpde4 + STATES[97]);
  //C
  //C   STATES[69] PDE3 extracaveolar phosphorylated
  //C
  //C      RATES[69] = kfpdep*STATES[89]*(PDE3ecavtot-STATES[69])-kbpdep*STATES[69]
  RATES[69] = 0.0f;
  //C
  //C   STATES[70] PDE4 extracaveolar phosphorylated
  //C
  RATES[70] = kfpdep * STATES[89] * (ALGEBRAIC[59]/*pde4ecavtot*/ - STATES[70]) - kbpdep * STATES[70];
  //C
  //C   STATES[71] cAMP change due to PDE2 extracaveolar domain
  //C
  RATES[71] = (kpde2 *ALGEBRAIC[57]/*pde2ecavtot*/* STATES[98]) / (kmpde2 + STATES[98]);
  //C
  //C   STATES[72] cAMP change due to PDE3 extracaveolar domain
  //C
  RATES[72] = (kpde3 * (ALGEBRAIC[58]/*pde3ecavtot*/ - STATES[69]) * STATES[98] + deltakpde34 *
    kpde3 * STATES[69] * STATES[98]) / (kmpde3 + STATES[98]);
  //C
  //C   STATES[73] cAMP change due to PDE4 extracaveolar domain
  //C
  RATES[73] = (kpde4 * (ALGEBRAIC[59]/*pde4ecavtot*/ - STATES[70]) * STATES[98] + deltakpde34 *
    kpde4 * STATES[70] * STATES[98]) / (kmpde4 + STATES[98]);
  //C
  //C   STATES[74] PDE3 cytosol phosphorylated
  //C
  RATES[74] = kfpdep * STATES[95] * (ALGEBRAIC[61]/*pde3cyttot*/ - STATES[74]) - kbpdep * STATES[74];
  //C
  //C   STATES[75] PDE4 cytosol phosphorylated
  //C
  RATES[75] = kfpdep * STATES[95] * (ALGEBRAIC[62]/*pde4cyttot*/ - STATES[75]) - kbpdep * STATES[75];
  //C
  //C   STATES[76] cAMP change due to PDE2 cytosol domain
  //C
  RATES[76] = (kpde2 * ALGEBRAIC[60]/*pde2cyttot*/ * STATES[99]) / (kmpde2 + STATES[99]);
  //C
  //C   STATES[77] cAMP change due to PDE3 cytosol domain
  //C
  RATES[77] = (kpde3 * (ALGEBRAIC[61]/*pde3cyttot*/ - STATES[74]) * STATES[99] + deltakpde34 *
    kpde3 * STATES[74] * STATES[99]) / (kmpde3 + STATES[99]);
  //C
  //C   STATES[78] cAMP change due to PDE4 cytosol domain
  //C
  RATES[78] = (kpde4 * (ALGEBRAIC[62]/*pde4cyttot*/ - STATES[75]) * STATES[99] + deltakpde34 *
    kpde4 * STATES[75] * STATES[99]) / (kmpde4 + STATES[99]);
  //C
  //C   STATES[79] cAMP binding to PKA caveolar
  //C
  RATES[79] = -kpkaiif1 *ALGEBRAIC[35]/*rccavf*/* STATES[97] + kpkaiib1 * STATES[80] -
    kpkaiif2 * STATES[80] * STATES[97] + kpkaiib2 * STATES[81];
  //C      write[51,44] time,-kpkaiif1*RCcavf*STATES[97],kpkaiib1*STATES[80],
  //C     *                  -kpkaiif2*STATES[80]*STATES[97],kpkaiib2*STATES[81]
  //C
  //C   STATES[80] RC with bound one cAMP in caveolae
  //C
  RATES[80] = kpkaiif1 *ALGEBRAIC[35]/*rccavf*/* STATES[97] - kpkaiib1 * STATES[80] -
    kpkaiif2 * STATES[80] * STATES[97] + kpkaiib2 * STATES[81];
  //C
  //C   STATES[81] RC with bound two cAMP in caveolae
  //C
  RATES[81] = kpkaiif2 * STATES[80] * STATES[97] - (kpkaiib2 + kpkaiif3) * STATES[
    81] + kpkaiib3 * STATES[82] * STATES[83];
  //C
  //C   STATES[82] R with bound two cAMP in caveolae
  //C
  RATES[82] = kpkaiif3 * STATES[81] - kpkaiib3 * STATES[82] * STATES[83];
  //C
  //C   STATES[83] concentration of catalytic subunit of PKA in caveolae
  //C
  RATES[83] = kpkaiif3 * STATES[81] - kpkaiib3 * STATES[82] * STATES[83] + kpkib * STATES[
    84] - kpkif *ALGEBRAIC[74]/*pkicavf*/* STATES[83];
  //C
  //C   STATES[84] concentration of PKI bound to catalytic subunit of PKA cav
  //C
  RATES[84] = -kpkib * STATES[84] + kpkif *ALGEBRAIC[74]/*pkicavf*/* STATES[83];
  //C
  //C   Extracaveolar domain
  //C
  //C   STATES[85] cAMP binding to PKA extracaveolar
  //C
  RATES[85] = -kpkaiif1 *ALGEBRAIC[36]/*rcecavf*/* STATES[98] + kpkaiib1 * STATES[86] -
    kpkaiif2 * STATES[86] * STATES[98] + kpkaiib2 * STATES[87];
  //C
  //C   STATES[86] RC with bound one cAMP in extracaveolae
  //C
  RATES[86] = kpkaiif1 *ALGEBRAIC[36]/*rcecavf*/* STATES[98] - kpkaiib1 * STATES[86] -
    kpkaiif2 * STATES[86] * STATES[98] + kpkaiib2 * STATES[87];
  //C
  //C   STATES[87] RC with bound two cAMP in extracaveolae
  //C
  RATES[87] = kpkaiif2 * STATES[86] * STATES[98] - (kpkaiib2 + kpkaiif3) * STATES[
    87] + kpkaiib3 * STATES[88] * STATES[89];
  //C
  //C   STATES[88] R with bound two cAMP in extracaveolae
  //C
  RATES[88] = kpkaiif3 * STATES[87] - kpkaiib3 * STATES[88] * STATES[89];
  //C
  //C   STATES[89] concentration of catalytic subunit of PKA in extracaveolae
  //C
  RATES[89] = kpkaiif3 * STATES[87] - kpkaiib3 * STATES[88] * STATES[89] + kpkib * STATES[
    90] - kpkif *ALGEBRAIC[75]/*pkiecavf*/* STATES[89];
  //C
  //C   STATES[90] concentration of PKI bound to catalytic subunit of PKA ecav
  //C
  RATES[90] = -kpkib * STATES[90] + kpkif *ALGEBRAIC[75]/*pkiecavf*/* STATES[89];
  //C
  //C   Cytosol domain
  //C
  //C   STATES[91] cAMP binding to PKA cytosol
  //C
  RATES[91] = -kpkaif1 * ALGEBRAIC[37]/*rccytf*/ * STATES[99] + kpkaib1 * STATES[92] - kpkaif2 *
    STATES[92] * STATES[99] + kpkaib2 * STATES[93];
  //C
  //C   STATES[92] RC with bound one cAMP in cytosol
  //C
  RATES[92] = kpkaif1 * ALGEBRAIC[37]/*rccytf*/ * STATES[99] - kpkaib1 * STATES[92] - kpkaif2 * STATES[
    92] * STATES[99] + kpkaib2 * STATES[93];
  //C
  //C   STATES[93] RC with bound two cAMP in cytosol
  //C
  RATES[93] = kpkaif2 * STATES[92] * STATES[99] - (kpkaib2 + kpkaif3) * STATES[93] +
    kpkaib3 * STATES[94] * STATES[95];
  //C
  //C   STATES[94] R with bound two cAMP in cytosol
  //C
  RATES[94] = kpkaif3 * STATES[93] - kpkaib3 * STATES[94] * STATES[95];
  //C
  //C   STATES[95] concentration of catalytic subunit of PKA in cytosol
  //C
  RATES[95] = kpkaif3 * STATES[93] - kpkaib3 * STATES[94] * STATES[95] + kpkib * STATES[
    96] - kpkif * ALGEBRAIC[76]/*pkicytf*/ * STATES[95];
  //C
  //C   STATES[96] concentration of PKI bound to catalytic subunit of PKA cyt
  //C
  RATES[96] = -kpkib * STATES[96] + kpkif * ALGEBRAIC[76]/*pkicytf*/ * STATES[95];
  //C
  //C   STATES[101] cAMP change due to PDE8 caveolar domain
  //C
  RATES[101] = (kpde8 *ALGEBRAIC[63]/*pde8cavtot*/* STATES[97]) / (kmpde8 + STATES[97]);
  //C
  //C   STATES[97] total change of cAMP in ceveolae
  //C
  RATES[97] = RATES[79] + RATES[60] - RATES[66] - RATES[67] - RATES[68] - RATES[
    101] - jcavecav * (STATES[97] - STATES[98]) / vcav - jcavcyt * (STATES[97] - STATES[
    99]) / vcav;
    ALGEBRAIC[45]/*fluxcavecav*/ = jcavecav * (STATES[97] - STATES[98]);
    ALGEBRAIC[46]/*fluxcavcyt*/ = jcavcyt * (STATES[97] - STATES[99]);
  //C
  //C   STATES[102] cAMP change due to PDE8 extracaveolar domain
  //C
  RATES[102] = (kpde8 *ALGEBRAIC[64]/*pde8ecavtot*/* STATES[98]) / (kmpde8 + STATES[98]);
  //C
  //C   STATES[98] total change of cAMP in extracaveolae
  //C
  RATES[98] = RATES[85] + RATES[61] - RATES[71] - RATES[72] - RATES[73] - RATES[
    102] - jcavecav * (STATES[98] - STATES[97]) / vecav - jecavcyt * (STATES[98] - STATES[
    99]) / vecav;
    ALGEBRAIC[47]/*fluxecavcav*/= jcavecav * (STATES[98] - STATES[97]);
    ALGEBRAIC[48]/*fluxecavcyt*/ = jecavcyt * (STATES[98] - STATES[99]);
  //C
  //C   STATES[103] cAMP change due to PDE8 cytosolic domain
  //C
  RATES[103] = (kpde8 * ALGEBRAIC[65]/*pde8cyttot*/ * STATES[99]) / (kmpde8 + STATES[99]);
  //C
  //C   STATES[99] total change of cAMP in cytosol
  //C
  RATES[99] = RATES[91] + RATES[62] + RATES[63] - RATES[76] - RATES[77] -
    RATES[78] - RATES[103] - jcavcyt * (STATES[99] - STATES[97]) / vcyt -
    jecavcyt * (STATES[99] - STATES[98]) / vcyt;
    ALGEBRAIC[49]/*fluxcytcav*/ = jcavcyt * (STATES[99] - STATES[97]);
    ALGEBRAIC[50]/*fluxcytecav*/ = jecavcyt * (STATES[99] - STATES[98]);
  //C
  //C    STATES[100] Inhibitor 1 cytosol phosphorylated
  //C
  RATES[100] = kpkainhib1 * STATES[95] * inhib1cytf / (kmpkainhib1 +
    inhib1cytf) - kpp2ainhib1pp2acyt * STATES[100] / (kmpp2ainhib1 + STATES[
    100]);
  //C
  //C    STATES[104] PLB fraction phosphorylated
  //C
  RATES[104] = kplbpka * STATES[95] * (1 - STATES[104]) / (kplbmpka + (1 - STATES[
    104])) - kplbpp1 * ALGEBRAIC[89]/*pp1cytf*/ * STATES[104] / (kplbmpp1 + STATES[104]);
  //C
  //C    STATES[105] TnI fraction phosphorylated
  //C
  RATES[105] = ktnipka * STATES[95] * (1 - STATES[105]) / (ktnimpka + (1 - STATES[
    105])) - ktnipp2a * pp2acyt * STATES[105] / (ktnimpp2a + STATES[105]);
  //C
  //C    STATES[129] INaK fraction phosphorylated
  //C
  RATES[129] = kinakpka * STATES[83] * (1 - STATES[129]) / (kinakmpka + (1 - STATES[
    129])) - kinakpp * ppcav * STATES[129] / (kinakmpp + STATES[129]);
  //C
  //C    STATES[133] IK1 fraction nonphophorylated
  //C
  RATES[133] = kik1pp * pp1ecav * (1 - STATES[133]) / (kik1mpp + (1 - STATES[
    133])) - kik1pka * STATES[89] * STATES[133] / (kik1mpka + STATES[133]);
  //C
  //C    STATES[134] IKto fraction phophorylated
  //C
  RATES[134] = kiktopka * STATES[89] * (1 - STATES[134]) / (kiktompka + (1 - STATES[
    134])) - kiktopp * pp1ecav * STATES[134] / (kiktompp + STATES[134]);
  //C
  //C    STATES[163] ICaT C0 channel variable
  //C
  RATES[163] = k_catmv * STATES[164] - 4.0f * k_catpv * STATES[163] + k_catmi *
    STATES[169] / pow(h_cat,3) - k_catpi * STATES[163] * pow(f_cat,3);
  //C
  //C    STATES[164] ICaT C1 channel variable
  //C
  RATES[164] = 4.0f * k_catpv * STATES[163] - k_catmv * STATES[164] + 2.0f *
    k_catmv * STATES[165] - 3.0f * k_catpv * STATES[164] + k_catmi * STATES[170] /
    pow(h_cat,2) - k_catpi * STATES[164] * pow(f_cat,2);
  //C
  //C    STATES[165] ICaT C2 channel variable
  //C
  RATES[165] = 3.0f * k_catpv * STATES[164] - 2.0f * k_catmv * STATES[165] +
    3.0f * k_catmv * STATES[166] - 2.0f * k_catpv * STATES[165] + k_catmi * STATES[
    171] / h_cat - k_catpi * STATES[165] * f_cat;
  //C
  //C    STATES[166] ICaT C3 channel variable
  //C
  RATES[166] = 2.0f * k_catpv * STATES[165] - 3.0f * k_catmv * STATES[166] + 4.0f *
    k_catmv * STATES[167] - k_catpv * STATES[166] + k_catmi * STATES[172] - k_catpi * STATES[
    166];
  //C
  //C    STATES[167] ICaT C4 channel variable
  //C
  RATES[167] = k_catpv * STATES[166] - 4.0f * k_catmv * STATES[167] + k_catmo *
    STATES[168] - k_catpo * STATES[167] + k_catmi * STATES[173] - k_catpi * STATES[167];
  //C
  //C    STATES[168] ICaT O channel variable
  //C
  RATES[168] = k_catpo * STATES[167] - k_catmo * STATES[168] + k_catmi * STATES[
    174] - k_catpi * STATES[168];
  //C
  //C    STATES[169] ICaT I0 channel variable
  //C
  RATES[169] = k_catmv * STATES[170] * h_cat - 4.0f * k_catpv * STATES[169] /
    f_cat + k_catpi * STATES[163] * pow(f_cat,3) - k_catmi * STATES[169] /
    pow(h_cat,3);
  //C
  //C    STATES[170] ICaT I1 channel variable
  //C
  RATES[170] = 4.0f * k_catpv * STATES[169] / f_cat - k_catmv * STATES[170] *
    h_cat + 2.0f * k_catmv * STATES[171] * h_cat - 3.0f * k_catpv * STATES[170] /
    f_cat + k_catpi * STATES[164] * pow(f_cat,2) - k_catmi * STATES[170] /
    pow(h_cat,2);
  //C
  //C    STATES[171] ICaT I2 channel variable
  //C
  RATES[171] = 3.0f * k_catpv * STATES[170] / f_cat - 2.0f * k_catmv * STATES[
    171] * h_cat + 3.0f * k_catmv * STATES[172] * h_cat - 2.0f * k_catpv *
    STATES[171] / f_cat + k_catpi * STATES[165] * f_cat - k_catmi * STATES[171] /
    h_cat;
  //C
  //C    STATES[172] ICaT I3 channel variable
  //C
  RATES[172] = 2.0f * k_catpv * STATES[171] / f_cat - 3.0f * k_catmv * STATES[
    172] * h_cat + 4.0f * k_catmv * STATES[173] - k_catpv * STATES[172] +
    k_catpi * STATES[166] - k_catmi * STATES[172];
  //C
  //C     STATES[173] ICaT I4 channel variable
  //C
  RATES[173] = k_catpv * STATES[172] - 4.0f * k_catmv * STATES[173] + k_catmo *
    STATES[174] - k_catpo * STATES[173] + k_catpi * STATES[167] - k_catmi * STATES[173];
  //C
  //C    STATES[174] ICaT IO channel variable
  //C
  RATES[174] = k_catpo * STATES[173] - k_catmo * STATES[174] + k_catpi * STATES[
    168] - k_catmi * STATES[174];
  //C
  //C    STATES[175] ICaT fraction phophorylated caveolae
  //C
  RATES[175] = kicatpka * STATES[83] * (1 - STATES[175]) / (kicatmpka + (1 - STATES[
    175])) - kicatpp * ppcav * STATES[175] / (kicatmpp + STATES[175]);
  //C
  //C    STATES[176] ICaT fraction phophorylated cytosol
  //C
  RATES[176] = kicatpka * STATES[95] * (1 - STATES[176]) / (kicatmpka + (1 - STATES[
    176])) - kicatpp * (pp2acyt + ALGEBRAIC[89]/*pp1cytf*/) * STATES[176] / (kicatmpp + STATES[
    176]);
  //C
  //C  Small-conductance calcium activated K+ channel
  //C
  RATES[177] = k_cakb1 * STATES[178] - k_cakf1 * STATES[177];
  RATES[178] = k_cakf1 * STATES[177] - k_cakb1 * STATES[178] + k_cakb2 * STATES[
    179] - k_cakf2 * STATES[178];
  RATES[179] = k_cakf2 * STATES[178] - k_cakb2 * STATES[179] + k_cakb3 * STATES[
    180] - k_cakf3 * STATES[179] + k_cakbo1 * STATES[181] - k_cakfo1 * STATES[179];
  RATES[180] = k_cakf3 * STATES[179] - k_cakb3 * STATES[180] + k_cakbo2 * STATES[
    182] - k_cakfo2 * STATES[180];
  RATES[181] = k_cakfo1 * STATES[179] - k_cakbo1 * STATES[181];
  RATES[182] = k_cakfo2 * STATES[180] - k_cakbo2 * STATES[182];
  //C
  //C    STATES[183] ICaK fraction phophorylated caveolae
  //C
  RATES[183] = kicakpka * STATES[83] * (1 - STATES[183]) / (kicakmpka + (1 - STATES[
    183])) - kicakpp * ppcav * STATES[183] / (kicakmpp + STATES[183]);
  //C
}



  void initialize_states_default(double *STATES) {
      //time = 0.0f;
      STATES[1] = -0.8088310412e+02d;
      STATES[2] = 0.9279026900e-01d;
      STATES[3] = 0.1620528775e+04d;
      STATES[4] = 0.8053999501e+01d;
      STATES[5] = 0.1222160781e+03d;
      STATES[6] = 0.1620528775e+04d;
      STATES[7] = 0.9279026946e-01d;
      STATES[8] = 0.8851691228e-12d;
      STATES[9] = 0.9615600724e+00d;
      STATES[10] = 0.3767463428e-02d;
      STATES[11] = 0.5535450499e-05d;
      STATES[12] = 0.3614717015e-08d;
      STATES[13] = 0.7908511386e-12d;
      STATES[14] = 0.7153680173e-08d;
      STATES[15] = 0.6391662360e-08d;
      STATES[16] = 0.9967529680e+00d;
      STATES[17] = 0.7087660992e-04d;
      STATES[18] = 0.6300248702e-05d;
      STATES[19] = 0.2112480482e-10d;
      STATES[20] = 0.4628585970e+00d;
      STATES[21] = 0.1062737990e-01d;
      STATES[22] = 0.9600858225e-04d;
      STATES[23] = 0.1376322742e+05d;
      STATES[24] = 0.1418267743e+06d;
      STATES[25] = 0.4139033547e-02d;
      STATES[26] = 0.9999623535e+00d;
      STATES[27] = 0.3298367341e-03d;
      STATES[28] = 0.5091689794e-03d;
      STATES[29] = 0.9980927689e+00d;
      STATES[30] = 0.9979010383e+00d;
      STATES[31] = 0.1113258849e-02d;
      STATES[32] = 0.7192382272e-03d;
      STATES[33] = 0.2223437325e-03d;
      STATES[34] = 0.4412087393e-04d;
      STATES[35] = 0.5091689794e-03d;
      STATES[36] = 0.9980927689e+00d;
      STATES[37] = 0.1552405441e-06d;
      STATES[38] = 0.6461071366e-04d;
      STATES[39] = 0.3762774987e-05d;
      STATES[40] = 0.8212756509e-08d;
      STATES[41] = 0.7151890060e-02d;
      STATES[42] = 0.3114893963e+00d;
      STATES[43] = 0.5091689794e-03d;
      STATES[44] = 0.3779115232e-12d;
      STATES[45] = 0.1193078173e-02d;
      STATES[46] = 0.3630488755e-61d;
      STATES[47] = 0.3490326981e-02d;
      STATES[48] = 0.5883268221e-02d;
      STATES[49] = 0.3955088510e-03d;
      STATES[50] = 0.4450890836e-01d;
      STATES[51] = 0.3630488755e-61d;
      STATES[52] = 0.2392020281e-01d;
      STATES[53] = 0.6409630097e-02d;
      STATES[54] = 0.2487944176e-02d;
      STATES[55] = 0.1722672905e-02d;
      STATES[56] = 0.3630488755e-61d;
      STATES[57] = 0.3771728607e-03d;
      STATES[58] = 0.7236469606e-03d;
      STATES[59] = 0.3474740964e-03d;
      STATES[60] = 0.1300635182e+06d;
      STATES[61] = 0.5077713406e+05d;
      STATES[62] = 0.2353396779e+05d;
      STATES[63] = 0.2193421483e+04d;
      STATES[64] = 0.1951294807e-01d;
      STATES[65] = 0.9059006924e-02d;
      STATES[66] = 0.1097267542e+05d;
      STATES[67] = 0.9355062972e+05d;
      STATES[68] = 0.2873599432e+05d;
      STATES[69] = 0.6773138136e-90d;
      STATES[70] = 0.1449457584e-01d;
      STATES[71] = 0.1340498892e+05d;
      STATES[72] = 0.2242397227e+00d;
      STATES[73] = 0.3690967042e+05d;
      STATES[74] = 0.1340729197e-02d;
      STATES[75] = 0.4134187505e-02d;
      STATES[76] = 0.4430125145e+04d;
      STATES[77] = 0.7006267909e+04d;
      STATES[78] = 0.1421393882e+05d;
      STATES[79] = 0.7736146378e+01d;
      STATES[80] = 0.3812957050e+00d;
      STATES[81] = 0.5158797692e-01d;
      STATES[82] = 0.8896934822e+00d;
      STATES[83] = 0.7537918546e-01d;
      STATES[84] = 0.8253133062e+00d;
      STATES[85] = 0.6836347022e+01d;
      STATES[86] = 0.6205890724e+00d;
      STATES[87] = 0.1172793205e+00d;
      STATES[88] = 0.1154255173e+01d;
      STATES[89] = 0.1320878782e+00d;
      STATES[90] = 0.1033166329e+01d;
      STATES[91] = 0.9274466720e+01d;
      STATES[92] = 0.1086997100e+00d;
      STATES[93] = 0.1546752825e-01d;
      STATES[94] = 0.2929517738e+00d;
      STATES[95] = 0.6863855598e-01d;
      STATES[96] = 0.2353122184e+00d;
      STATES[97] = 0.3382412669e+00d;
      STATES[98] = 0.4724516017e+00d;
      STATES[99] = 0.4126582483e+00d;
      STATES[100] = 0.2320312876e-01d;
      STATES[101] = 0.6026039427e+00d;
      STATES[102] = 0.6026039427e+00d;
      STATES[103] = 0.2596146734e+00d;
      STATES[104] = 0.3675832948e+00d;
      STATES[105] = 0.3813640007e+00d;
      STATES[106] = 0.8851705616e-12d;
      STATES[107] = 0.2605069110e-10d;
      STATES[108] = 0.3415627435e-01d;
      STATES[109] = 0.5076951530e-03d;
      STATES[110] = 0.2829869165e-05d;
      STATES[111] = 0.7010470517e-08d;
      STATES[112] = 0.2327573109e-10d;
      STATES[113] = 0.5549734835e-07d;
      STATES[114] = 0.4958553899e-07d;
      STATES[115] = 0.6512671925e-11d;
      STATES[116] = 0.3169194559e-02d;
      STATES[117] = 0.6261242517e-06d;
      STATES[118] = 0.3339118148e-07d;
      STATES[119] = 0.1866017795e-12d;
      STATES[120] = 0.1213360465e+00d;
      STATES[121] = 0.2785949016e-02d;
      STATES[122] = 0.2516918244e-04d;
      STATES[123] = 0.4069800161e-07d;
      STATES[124] = 0.1694073185e-04d;
      STATES[125] = 0.1028703688e-05d;
      STATES[126] = 0.6604024237e-08d;
      STATES[127] = 0.1875148561e-02d;
      STATES[128] = 0.8166786118e-01d;
      STATES[129] = 0.3801202195e+00d;
      STATES[130] = 0.5091689794e-03d;
      STATES[131] = 0.9980927689e+00d;
      STATES[132] = 0.9214774521e+00d;
      STATES[133] = 0.9908502903e+00d;
      STATES[134] = 0.2087704319e+00d;
      STATES[135] = 0.8637307739e-03d;
      STATES[136] = 0.9999881755e+00d;
      STATES[137] = 0.8187695769e-12d;
      STATES[138] = 0.8894320487e+00d;
      STATES[139] = 0.3484859972e-02d;
      STATES[140] = 0.5120227116e-05d;
      STATES[141] = 0.3343570256e-08d;
      STATES[142] = 0.8187721513e-12d;
      STATES[143] = 0.7315039113e-12d;
      STATES[144] = 0.6616944790e-08d;
      STATES[145] = 0.5912117549e-08d;
      STATES[146] = 0.8046441685e-10d;
      STATES[147] = 0.1055007176e+00d;
      STATES[148] = 0.1568150421e-02d;
      STATES[149] = 0.8740796912e-05d;
      STATES[150] = 0.2165368505e-07d;
      STATES[151] = 0.2011610269e-10d;
      STATES[152] = 0.7189314575e-10d;
      STATES[153] = 0.1714179939e-06d;
      STATES[154] = 0.1531578593e-06d;
      STATES[155] = 0.6077748720e-01d;
      STATES[156] = 0.5796341475e-38d;
      STATES[157] = 0.3591463241e-02d;
      STATES[158] = 0.4069691772e-03d;
      STATES[159] = 0.4580528431e-03d;
      STATES[160] = 0.5796341475e-39d;
      STATES[161] = 0.2305390265e-05d;
      STATES[162] = 0.2397840156e-06d;
      STATES[163] = 0.8422801321e+00d;
      STATES[164] = 0.5090725598e-01d;
      STATES[165] = 0.1153809439e-02d;
      STATES[166] = 0.1162267165e-04d;
      STATES[167] = 0.4390450831e-07d;
      STATES[168] = 0.6508479153e-06d;
      STATES[169] = 0.6016286658e-01d;
      STATES[170] = 0.3636232570e-01d;
      STATES[171] = 0.8241495995e-02d;
      STATES[172] = 0.8301908322e-03d;
      STATES[173] = 0.3136036308e-05d;
      STATES[174] = 0.4648913681e-04d;
      STATES[175] = 0.2965524523e-01d;
      STATES[176] = 0.2144114378e-01d;
      STATES[177] = 0.8608445698e+00d;
      STATES[178] = 0.1198169988e+00d;
      STATES[179] = 0.1334142186e-01d;
      STATES[180] = 0.2971089895e-03d;
      STATES[181] = 0.2134627497e-02d;
      STATES[182] = 0.3565307875e-02d;
      STATES[183] = 0.1184324373e-01d;
      STATES[184] = 0.5091689794e-03d;
      STATES[185] = 0.9980927689e+00d;

  }
