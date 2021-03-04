#define _GNU_SOURCE
#include <stdio.h>
#include <fenv.h>
#include "math.h"

void fun(double time, double* y, double* ydot, double* params)
{
    y[9] = 1.0f - (y[8] + y[10] + y[11] + y[12] + y[13] + y[14] + y[15] + y[106] + y[107] + y[108] + y[109] + y[110] + y[111] + y[112] + y[113] + y[114] + y[115]);
    y[138] = 1.0f - (y[137] + y[139] + y[140] + y[141] + y[142] + y[143] + y[144] + y[145] + y[146] + y[147] + y[148] + y[149] + y[150] + y[151] + y[152] + y[153] + y[154]);
    y[16] = 1.0-(y[17]+y[18]+y[19]+y[116]+y[117]+y[118]+y[119]);
    y[20] = 1.0f - (y[21] + y[22] + y[37] + y[38] + y[39] + y[40] + y[41] + y[42] + y[120] + y[121] + y[122] + y[123] + y[124] + y[125] + y[126] + y[127] + y[128]);
    y[30] = 1.0f - (y[31] + y[32] + y[33] + y[34]);
    double CL = params[0];
    int stims_passed = (int)floor(time/CL);
    double r_t = time - stims_passed*CL;
    double istim = 0.0;
    if(r_t <= 1.0f)
    {
        istim = 80.0;
    }else{
        istim = 0.0;
    }

  // COMMON a1
  double jrelt, jtrt, jxfert, jleak, jup, jtrpn;

  // COMMON a2
  double acap, vmyo, vjsr, vnsr, vss, cm, f, bi, bss, bjsr;

  // COMMON a10
  double icas, icab, inaca, ipca, ina, inab, iclca, icasc, icase;

  // COMMON a11
  double inak, ikto, ik1, iks, ikur, ikr, ikur1, ikur2, ikur3, ikur4;

  // COMMON b1
  double rccavf, rcecavf, rccytf, gscavabg, gsecavabg, gscytabg;

  // COMMON b5
  double rb1cavtot, rb1ecavtot, rb1cyttot;

  // COMMON b6
  double pde2cavtot, pde3cavtot, pde4cavtot, pde2ecavtot,
  pde3ecavtot, pde4ecavtot, pde2cyttot, pde3cyttot,
  pde4cyttot, pde8cavtot, pde8ecavtot, pde8cyttot;
  // COMMON b7
  double pdeact2, pdeact3, pdeact4, pdeact8;
  // COMMON b8
  double pkacav, pkaecav, pkacyt,
  pkicavf, pkiecavf, pkicytf;
  // COMMON b9
  double rb1ppka, rb1pbark, l, rb2ppka,rb2pbark;
  // COMMON b10
  double pde2mem, pde3mem, pde4mem, pde8mem;
  // COMMON b11
  double pp1cytf, ficalcav;
  // COMMON b12
  double gicavabg, giecavabg, rb2cavtot, rb2ecavtot;
  // COMMON c1
  double rb2gicavpka, lrb2gicavpka, rb2cavpkaf,
  gicavf, acavb2i, bcavb2i, ccavb2i;
  // COMMON c2
  double rb2giecavpka, lrb2giecavpka, rb2ecavpkaf,
  giecavf, aecavb2i, becavb2i, cecavb2i;
  // COMMON d1
  double icat, icatc, icats, icak;

  l = 0;
  //
  //C
  //C****************************************************************
  //C
  //C Cell Geometry Parameters
  //C
  acap = 0.5931e-4f;
  double vcell = 22.1200e-6f;
  vmyo = 13.9600e-6f;
  double vcyt = 13.9600e-6f;
  vjsr = 0.1349e-6f;
  vnsr = 2.5523e-6f;
  vss = 1.4850e-9f * 2.0f;
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
  double v1 = 4.5f * 0.09f;//jrelt
  v1 = v1 * params[13];
  double v2 = 1.74e-5f * 3.0f;
  v2 = v2 * params[14];
  double v3 = 0.306f * 2.63f;
  v3 = v3 * params[15];
  //C      Kmup  =  0.5
  double kmup = (0.38f * (1 - y[104]) + 0.25f * y[104]) * 0.815f;
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
  double k_catpv = 9.5f * exp(y[1] / 15.0f);
  double k_catmv = 0.008f * exp(-1.0f * y[1] / 18.0f) * 4.0f;
  double k_catpo = 4.0f;
  double k_catmo = 0.025f * exp(-1.0f * y[1] / 34.0f);
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
  double kltrpnm = 0.0196f * (1 - y[105]) + 0.0294f * y[105];
  double ccmdntot = 50.0f;
  double ccsqntot = 15000.0f * 0.6f;
  double kmcmdn = 0.238f;
  double kmcsqn = 800.0f;
  //C
  //C Membrane Current Parameters
  //C
  cm = 1.0f;
  f = 96.5f;
  double t = 298.0f;
  double r = 8.314f;
  double factor = r * t / f;
  double ifactor = f / (r * t);
  double gna = 14.4f; //Sodium
  gna = gna * params[1];
  //double gnap = 18.0f; //Sodium phosphorylated
  double gnap = gna * 1.25;
  //Sodium- Calcium exchanger
  double knaca = 275.0f * 1.0f;
  knaca = knaca * params[2];
  double kmna = 87500.0f;
  double kmca = 1380.0f;
  double ksat = 0.27f;
  double nu = 0.35f;
  //Sodium-Potassium exchanger
  double inakmax = 4.0f * 0.5f;
  inakmax = inakmax * params[3];
  //C      Kmnai   = 21000.0
  double kmnai = 18800.0f * (1 - y[129]) + 13600.0f * y[129];
  double kmko = 1500.0f;
  double sigma = (exp(cnao / 67300.0f) - 1.0f) / 7.0f;
  //Calcium pump
  double ipcamax = 0.051f * 1.5f * 2.0f;
  ipcamax = ipcamax * params[4];
  double kmpca = 0.5f;
  //Calcium background
  double gcab = 0.000284f * 0.75f;
  gcab = gcab * params[5];
  //Sodium background
  double gnab = 0.0063f * 1.0f;
  gnab = gnab * params[6];
  //Ikto??
  double dito = 7.0f;
  //C      Gks     =     0.00575
  //Iks
  double gks = 0.0f;
  gks = gks * params[7];
  //Ikto
  double gkto = 0.14067f;
  gkto = gkto * params[8];
  //double gktop = 0.14067f * (1.0f - 0.387f);??
  double gktop = gkto * 0.613f;
  //Ikur
  double gkur1 = 0.0f;
  double gkur2 = 0.05766f;
  double gkur2p = 0.07496f;
  double gkur3 = 0.0428f;
  double gkur4 = 0.05766f;
  //C
  //C HERG current parameters (IKr)
  //C
  double kf = 0.023761f;
  double kb = 0.036778f;
  double gkr = 0.078f * 10.0f;
  gkr = gkr * params[9];
  //C
  //C Calcium activated chloride current
  //C
  double gclca = 10.0f;
  gclca = gclca * params[10];
  double poclcamax = 0.2f;
  double kmclca = 10.0f;
  double ecl = -40.0f;
  //C
  //C Small-conductance calcium-activated K+ current
  //C
  double k_cakf1 = 0.2f * y[2] * 0.6f;
  double k_cakf2 = 0.16f * y[2] * 0.6f;
  double k_cakf3 = 0.08f * y[2] * 0.6f;
  double k_cakb1 = 0.08f;
  double k_cakb2 = 0.08f;
  double k_cakb3 = 0.2f;
  double k_cakfo1 = 0.16f;
  double k_cakbo1 = 1.0f;
  double k_cakfo2 = 1.2f;
  double k_cakbo2 = 0.1f;
  double g_cak_factor = params[16];
  double g_cak = (1.1e-6f * y[1] * y[1] * y[1] + 1.0e-4f * y[1] * y[
    1] + 0.0041f * y[1] + 0.1872f) * 4.0f * g_cak_factor;

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
  double fpdepart = 0.2f;
  double rpartpde23 = 0.570f;
  double rpartpde34 = 0.748f;
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
  rb1cavtot = fcavb1 * rb1tot * vcell / vcav;
  rb2cavtot = fcavb2 * rb2tot * vcell / vcav;
  gscavabg = fcavgs * gstot * vcell / vcav - y[47] - y[49];
  gicavabg = fcavgi * gitot * vcell / vcav - y[157] - y[158];
  double rb1cavnptot = rb1cavtot - y[45] - y[46];
  double rb2cavnptot = rb2cavtot - y[155] - y[156];
  acavb2i = (kb2l + l) * (kb2f + l) / kb2l;
  bcavb2i = gicavabg * (kb2f + l) - y[155] * (kb2f + l) + kb2a *
    kb2f * (1.0f + l / kb2l);
  ccavb2i = -1.0f * y[155] * kb2a * kb2f;

  rb2cavpkaf = (-bcavb2i + sqrt(pow(bcavb2i,2) - 4.0f *
    acavb2i * ccavb2i)) / (2.0f * acavb2i);

  gicavf = gicavabg / (1.0f + rb2cavpkaf * (1.0f / kb2a + l / (kb2a * kb2f)));
  double lrb2cavpka = l * rb2cavpkaf / kb2l;
  rb2gicavpka = rb2cavpkaf * gicavf / kb2a;
  lrb2gicavpka = l * rb2cavpkaf * gicavf / (kb2a * kb2f);
  double acavb2s = (kb1h + l) * (kb2h + l);
  double bcavb2s = (l + kb1h) * (l + kb2h) * (rb1cavnptot +
    rb2cavnptot) + (kb1c * kb1h + l * kb1c * kb1h / kb1l) * (l +
    kb2h) + (kb2c * kb2h + l * kb2c * kb2h / kb2l) * (l + kb1h) -
    gscavabg * (l + kb1h) * (l + kb2h);
  double ccavb2s = (l + kb1h) * (kb2c * kb2h + l * kb2c * kb2h /
    kb2l) * (rb1cavnptot - gscavabg) + (l + kb2h) * (kb1c * kb1h +
    l * kb1c * kb1h / kb1l) * (rb2cavnptot - gscavabg) + (kb1c *
    kb1h + l * kb1c * kb1h / kb1l) * (kb2c * kb2h + l * kb2c * kb2h /
    kb2l);
  double dcavb2s = -gscavabg * (kb1c * kb1h + l * kb1c * kb1h /
    kb1l) * (kb2c * kb2h + l * kb2c * kb2h / kb2l);
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
  double rb1cavnpf = rb1cavnptot / (1.0f + l / kb1l + gscavf * (l / (
    kb1c * kb1h) + 1.0f / kb1c));
  double rb2cavnpf = rb2cavnptot / (1.0f + l / kb2l + gscavf * (l / (
    kb2c * kb2h) + 1.0f / kb2c));
  double lrb1cavnp = l * rb1cavnpf / kb1l;
  double rb1gscavnp = rb1cavnpf * gscavf / kb1c;
  double lrb1gscavnp = l * rb1cavnpf * gscavf / (kb1c * kb1h);
  double lrb2cavnp = l * rb2cavnpf / kb2l;
  double rb2gscavnp = rb2cavnpf * gscavf / kb2c;
  double lrb2gscavnp = l * rb2cavnpf * gscavf / (kb2c * kb2h);
  double ac56cav = fcavac56 * fac56ac47 * actot * (vcell / vcav);
  double kcavac56 = af56 * (ac56bas + pow(y[47], hac56gsa) / (
    kac56mgsa + pow(y[47], hac56gsa))) * (1 + vac56gbg * pow(y[48],
    hac56gsbg) / (kac56mgsbg + pow(y[48], hac56gsbg))) * (1.0f -
    (1.0f - vac56gsgi * pow(y[47], hac56gsgi) / (kac56mgsgi +
    pow(y[47], hac56gsgi))) * y[157] / (kac56mgi + y[157]));
  pde2cavtot = (1.0f - pow(ibmx, hibmxpde2) / (kibmxmpde2 + pow(ibmx,
    hibmxpde2))) * fcavpde2 * pde2tot * (vcell / vcav);
  pde3cavtot = (1.0f - pow(ibmx, hibmxpde3) / (kibmxmpde3 + pow(ibmx,
    hibmxpde3))) * fcavpde3 * pde3tot * (vcell / vcav);
  pde4cavtot = (1.0f - pow(ibmx, hibmxpde4) / (kibmxmpde4 + pow(ibmx,
    hibmxpde4))) * fcavpde4 * pde4tot * (vcell / vcav);
  pde8cavtot = fcavpde8 * pde8tot * (vcell / vcav);
  pkacav = fcavpka * pkatot * (vcell / vcav);
  rccavf = 2.0f * pkacav - y[80] - y[81] - y[82];
  pkicavf = fcavpki * pkitot * (vcell / vcav) - y[84];
  double kpkaiib1 = kpkaiif1 * kpkaii1;
  double kpkaiib2 = kpkaiif2 * kpkaii2;
  double kpkaiib3 = kpkaiif3 / kpkaii3;
  double kpkib = kpkif * kpki;
  //C
  //C Extracaveolae
  //C
  rb1ecavtot = fecavb1 * rb1tot * vcell / vecav;
  rb2ecavtot = fecavb2 * rb2tot * vcell / vecav;
  gsecavabg = fecavgs * gstot * vcell / vecav - y[52] - y[54];
  giecavabg = fecavgi * gitot * vcell / vecav - y[161] - y[162];
  double rb1ecavnptot = rb1ecavtot - y[50] - y[51];
  double rb2ecavnptot = rb2ecavtot - y[159] - y[160];
  aecavb2i = (kb2l + l) * (kb2f + l) / kb2l;
  becavb2i = giecavabg * (kb2f + l) - y[159] * (kb2f + l) + kb2a *
    kb2f * (1.0f + l / kb2l);
  cecavb2i = -1.0f * y[159] * kb2a * kb2f;
  rb2ecavpkaf = (-becavb2i + sqrt(pow(becavb2i,2) - 4.0f *
    aecavb2i * cecavb2i)) / (2.0f * aecavb2i);
  giecavf = giecavabg / (1.0f + rb2ecavpkaf * (1.0f / kb2a + l / (
    kb2a * kb2f)));
  double lrb2ecavpka = l * rb2ecavpkaf / kb2l;
  rb2giecavpka = rb2ecavpkaf * giecavf / kb2a;
  lrb2giecavpka = l * rb2ecavpkaf * giecavf / (kb2a * kb2f);
  double aecavb2s = (kb1h + l) * (kb2h + l);
  double becavb2s = (l + kb1h) * (l + kb2h) * (rb1ecavnptot +
    rb2ecavnptot) + (kb1c * kb1h + l * kb1c * kb1h / kb1l) * (l +
    kb2h) + (kb2c * kb2h + l * kb2c * kb2h / kb2l) * (l + kb1h) -
    gsecavabg * (l + kb1h) * (l + kb2h);
  double cecavb2s = (l + kb1h) * (kb2c * kb2h + l * kb2c * kb2h /
    kb2l) * (rb1ecavnptot - gsecavabg) + (l + kb2h) * (kb1c * kb1h +
    l * kb1c * kb1h / kb1l) * (rb2ecavnptot - gsecavabg) + (kb1c *
    kb1h + l * kb1c * kb1h / kb1l) * (kb2c * kb2h + l * kb2c * kb2h /
    kb2l);
  double decavb2s = -gsecavabg * (kb1c * kb1h + l * kb1c * kb1h /
    kb1l) * (kb2c * kb2h + l * kb2c * kb2h / kb2l);
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
  double rb1ecavnpf = rb1ecavnptot / (1.0f + l / kb1l + gsecavf * (l /
    (kb1c * kb1h) + 1.0f / kb1c));
  double rb2ecavnpf = rb2ecavnptot / (1.0f + l / kb2l + gsecavf * (l /
    (kb2c * kb2h) + 1.0f / kb2c));
  double lrb1ecavnp = l * rb1ecavnpf / kb1l;
  double rb1gsecavnp = rb1ecavnpf * gsecavf / kb1c;
  double lrb1gsecavnp = l * rb1ecavnpf * gsecavf / (kb1c * kb1h);
  double lrb2ecavnp = l * rb2ecavnpf / kb2l;
  double rb2gsecavnp = rb2ecavnpf * gsecavf / kb2c;
  double lrb2gsecavnp = l * rb2ecavnpf * gsecavf / (kb2c * kb2h);
  double ac47ecav = fecavac47 * (1.0f - fac56ac47) * actot * vcell / vecav;
  double kecavac47 = af47 * (ac47bas + pow(y[52], hac47gsa) / (
    kac47mgsa + pow(y[52], hac47gsa))) * (1.0f + vac47gbg *
    pow(y[53], hac47gsbg) / (kac47mgsbg + pow(y[53],
    hac47gsbg)));
  pde2ecavtot = (1.0f - (pow(ibmx, hibmxpde2) / (kibmxmpde2 +
    pow(ibmx, hibmxpde2)))) * fecavpde2 * pde2tot * (vcell /
    vecav);
  pde3ecavtot = (1.0f - (pow(ibmx, hibmxpde3) / (kibmxmpde3 +
    pow(ibmx, hibmxpde3)))) * fecavpde3 * pde3tot * (vcell /
    vecav);
  pde4ecavtot = (1.0f - (pow(ibmx, hibmxpde4) / (kibmxmpde4 +
    pow(ibmx, hibmxpde4)))) * fecavpde4 * pde4tot * (vcell /
    vecav);
  pde8ecavtot = fecavpde8 * pde8tot * (vcell / vecav);
  pkaecav = fecavpka * pkatot * (vcell / vecav);
  rcecavf = 2.0f * pkaecav - y[86] - y[87] - y[88];
  pkiecavf = fecavpki * pkitot * (vcell / vecav) - y[90];
  kpkaiib1 = kpkaiif1 * kpkaii1;
  kpkaiib2 = kpkaiif2 * kpkaii2;
  kpkaiib3 = kpkaiif3 / kpkaii3;
  kpkib = kpkif * kpki;
  //C
  //C Cytosol
  //C
  rb1cyttot = fcytb1 * rb1tot * vcell / vcyt;
  gscytabg = fcytgs * gstot * vcell / vcyt - y[57] - y[59];
  double rb1cytnptot = rb1cyttot - y[55] - y[56];
  double acytb1 = (kb1l + l) * (kb1h + l) / kb1l;
  double bcytb1 = gscytabg * (kb1h + l) - rb1cytnptot * (kb1h + l) +
    kb1c * kb1h * (1.0f + l / kb1l);
  double ccytb1 = -rb1cytnptot * kb1c * kb1h;
  double rb1cytnpf = (-bcytb1 + sqrt(pow(bcytb1,2) - 4.0f *
    acytb1 * ccytb1)) / (2.0f * acytb1);
  double gscytf = gscytabg / (1.0f + rb1cytnpf * (1.0f / kb1c + l / (
    kb1c * kb1h)));
  double lrb1cytnp = l * rb1cytnpf / kb1l;
  double rb1gscytnp = rb1cytnpf * gscytf / kb1c;
  double lrb1gscytnp = l * rb1cytnpf * gscytf / (kb1c * kb1h);
  double ac56cyt = (1.0f - fcavac56) * fac56ac47 * actot * (vcell / vcyt);
  double kcytac56 = af56 * (ac56bas + pow(y[57], hac56gsa) / (
    kac56mgsa + pow(y[57], hac56gsa))) * (1.0f + vac56gbg *
    pow(y[58], hac56gsbg) / (kac56mgsbg + pow(y[58],
    hac56gsbg)));
  double ac47cyt = (1.0f - fecavac47) * (1.0f - fac56ac47) * actot *
    vcell / vcyt;
  double kcytac47 = af47 * (ac47bas + pow(y[57], hac47gsa) / (
    kac47mgsa + pow(y[57], hac47gsa))) * (1.0f + vac47gbg *
    pow(y[58], hac47gsbg) / (kac47mgsbg + pow(y[58],
    hac47gsbg)));
  pde2cyttot = (1.0f - (pow(ibmx, hibmxpde2) / (kibmxmpde2 +
    pow(ibmx, hibmxpde2)))) * fcytpde2 * pde2tot * (vcell /
    vcyt);
  pde3cyttot = (1.0f - (pow(ibmx, hibmxpde3) / (kibmxmpde3 +
    pow(ibmx, hibmxpde3)))) * fcytpde3 * pde3tot * (vcell /
    vcyt);
  pde4cyttot = (1.0f - (pow(ibmx, hibmxpde4) / (kibmxmpde4 +
    pow(ibmx, hibmxpde4)))) * fcytpde4 * pde4tot * (vcell /
    vcyt);
  pde8cyttot = fcytpde8 * pde8tot * (vcell / vcyt);
  pkacyt = fcytpka * pkatot * (vcell / vcyt);
  rccytf = 2.0f * pkacyt - y[92] - y[93] - y[94];
  pkicytf = fcytpki * pkitot * (vcell / vcyt) - y[96];
  double kpkaib1 = kpkaif1 * kpkai1;
  double kpkaib2 = kpkaif2 * kpkai2;
  double kpkaib3 = kpkaif3 / kpkai3;
  kpkib = kpkif * kpki;
  //C
  //C Protein phosphatase module
  //C
  double inhib1cytf = inhib1cyttot - y[100];
  double ainh1 = 1.0f;
  double kinh1 = 0.0f; //uninitialized in m9da.f
  double binh1 = kinh1 + pp1cyttot - y[100];
  double cinh1 = -y[100] * kinhib1;
  double inhib1cytp = -binh1 / 2.0f + sqrt(pow(binh1,2) -
    4.0f * ainh1 * cinh1) / 2.0f;
  pp1cytf = pp1cyttot * kinhib1 / (kinhib1 + inhib1cytp);
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
  ficalcav = 0.2f;
  double icalcav = ficalcav * icaltot * vcell / vcav;
  double icalecav = (1.0f - ficalcav) * icaltot * vcell / vecav;
  double kicalpka = 2.0e-5f * 0.87f;
  double kicalpp = 1.55e-7f * 1.5f;
  double kicalmpka = 0.5f;
  double kicalmpp = 0.2f;
  double CaLFactor = params[11];
  //C
  //C ICaT module
  //C
  double kicatpka = 1.74e-5f * 0.5f * 20.0f;
  double kicatpp = 2.33e-7f * 200.0f;
  double kicatmpka = 5.0f;
  double kicatmpp = 0.1f;
  double CaTFactor = params[12];
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
  double inatot = 0.0f;
  double inacav = inatot * vcell / vcav;
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
  double campcell = (y[97] * vcav + y[98] * vecav + y[99] * vcyt) / vcell;
  double ccell = (y[83] * vcav + y[89] * vecav + y[95] * vcyt) / vcell;
  double acactcell = (ydot[60] * vcav + ydot[61] * vecav + (ydot[62] +
    ydot[63]) * vcyt) / vcell;
  pdeact2 = (ydot[66] * vcav + ydot[71] * vecav + ydot[76] * vcyt) / vcell;
  pdeact3 = (ydot[67] * vcav + ydot[72] * vecav + ydot[77] * vcyt) / vcell;
  pdeact4 = (ydot[68] * vcav + ydot[73]* vecav + ydot[78] * vcyt) / vcell;
  pdeact8 = (ydot[101] * vcav + ydot[102] * vecav + ydot[103] * vcyt) / vcell;
  pde2mem = (ydot[66] * vcav + ydot[71] * vecav) / vcell;
  pde3mem = (ydot[67] * vcav + ydot[72] * vecav) / vcell;
  pde4mem = (ydot[68] * vcav + ydot[73] * vecav) / vcell;
  pde8mem = (ydot[101] * vcav + ydot[102] * vecav) / vcell;
  double pdeactcell = pdeact2 + pdeact3 + pdeact4 + pdeact8;
  double pdemem = pde2mem + pde3mem + pde4mem + pde8mem;
  double pkaactcell = (ydot[83] * vcav + ydot[89] * vecav + ydot[95] *
    vcyt) / vcell;
  rb1ppka = (y[45] * vcav + y[50] * vecav + y[55] * vcyt) / vcell;
  rb1pbark = (y[46] * vcav + y[51] * vecav + y[56] * vcyt) / vcell;
  double rb1ptot = rb1ppka + rb1pbark;
  rb2ppka = (y[155] * vcav + y[159] * vecav) / vcell;
  rb2pbark = (y[156] * vcav + y[160] * vecav) / vcell;
  double rb2ptot = rb2ppka + rb2pbark;
  //C
  //C***********************************************************************
  //C            Calcium fluxes
  //C***********************************************************************
  //C      Pryr     = exp((y[1]-5.0]*(y[1]-5.0]/(-648.))
  jrelt = v1 * (y[18] + y[19] + y[118] + y[119]) * (y[6] - y[7]) * y[44];
  jtrt = (y[3] - y[6]) / ttr;
  jxfert = (y[7] - y[2]) / txfer;
  jleak = v2 * (y[3] - y[2]);
  jup = v3 * y[2] * y[2] / (kmup * kmup + y[2] * y[2]);
  double temp12 = khtrpnp * y[2] * (chtrpntot - y[5]) - khtrpnm * y[5];
  double temp13 = kltrpnp * y[2] * (cltrpntot - y[4]) - kltrpnm * y[4];
  jtrpn = temp12 + temp13;
  //C***********************************************************************
  //C             Ionic currents
  //C***********************************************************************
  //C
  //C y[1]   : Membrane potential (mV)
  //C
  double v = y[1];
  //C
  //C Icab   : Calcium background current
  //C
  double ecan = 0.5f * factor * log(ccao / y[2]);
  icab = gcab * (y[1] - ecan);
  //C
  //C Ipca   : Calcium pump current
  //C
  ipca = ipcamax * y[2] * y[2] / (kmpca * kmpca + y[2] * y[2]);
  //C
  //C Inaca  : Na-Ca exchange current
  //C
  double tempa1 = knaca / (kmna * kmna * kmna + cnao * cnao * cnao);
  double tempa2 = 1.0f / (kmca + ccao);
  double tempa3 = 1.0f / (1.0f + ksat * exp((nu - 1.0f) * y[1] * ifactor));
  double tempa41 = exp(nu * y[1] * ifactor) * y[23] * y[23] * y[23] * ccao;
  double tempa42 = exp((nu - 1.0f) * y[1] * ifactor) * cnao *
    cnao * cnao * y[2];
  double tempa4 = tempa41 - 2.0f * tempa42;
  inaca = tempa1 * tempa2 * tempa3 * tempa4;
  //C
  //C Ica    : L-type calcium current
  //C
  icasc = ficalcav * (sfica * y[8] + sficap * y[107]) * (y[1] - 52.0f);
  icase = (1.0f - ficalcav) * (sfica * y[137] + sficap * y[146]) * (y[
    1] - 52.0f);
  icas = CaLFactor * (icasc + icase);
  //C
  //C Icat   : T-type calcium current
  //C
  icatc = ficatcav * (g_cat * y[168] * (1 - y[175]) + g_catp * y[
    168] * y[175]) * (y[1] - 35.0f);
  icats = (1.0f - ficatcav) * (g_cat * y[168] * (1 - y[176]) +
    g_catp * y[168] * y[176]) * (y[1] - 35.0f);
  icat = CaTFactor * (icatc + icats);
  //C
  //C Ina    : Na fast current (Luo and Rudy, 1994]
  //C
  double ena = factor * log((0.9f * cnao + 0.1f * cko) / (0.9f *
    y[23] + 0.1f * y[24]));
  ina = (gna * y[37] + gnap * y[123]) * (y[1] - ena);
  //C
  //C Inab   : Na background current
  //C
  inab = gnab * (y[1] - ena);
  //C
  //C Inak    : Na-K exchange current
  //C
  double tempa5 = 1.0f / (1.0f + pow((kmnai / y[23]), 3.0f));
  double tempa6 = cko / (cko + kmko);
  double tempa11 = 1.0f + 0.1245f * exp(-0.1f * y[1] * ifactor) +
    0.0365f * sigma * exp(-1.0f * y[1] * ifactor);
  inak = inakmax * tempa5 * tempa6 / tempa11;
  //C
  //C Ikto    : transient outward current (Liu and Rasmusson, 1997]
  //C
  double ek = factor * log(cko / y[24]);
  ikto = gkto * y[25] * y[25] * y[25] * y[26] * (y[1] - ek) * (1 - y[
    134]) + gktop * y[135] * y[135] * y[135] * y[136] * (y[1] - ek) *
    y[134];
  //C
  //C Ik1    : Time independent K+ current (Rasmusson et al. 1990]
  //C
  //C      tempa7 = 0.3397*cKo/(cKo+210.0]
  //C      tempa8 = 1.0+exp(0.0448*(y[1]-EK))
  //C      tempa7p = 0.11394*cKo/(cKo+210.0]
  //C      tempa8p = 1.0+exp(0.0298*(y[1]-EK))
  //C      Ik1 = (tempa7*y[133]/tempa8
  //C     *    + tempa7*[1-y[133])/tempa8]*(y[1]-EK)
  //C
  double tempa7 = 1.02f / (1.0f + exp(0.2385f * (y[1] - ek - 59.215f)));
  double tempa8 = (0.8f * exp(0.08032f * (y[1] - ek + 5.476f)) +
    exp(0.06175f * (y[1] - ek - 594.31f))) / (1.0f + exp(
    -0.5143f * (y[1] - ek + 4.753f)));
  double gk1 = 1.7f * 0.27f * sqrt(cko / 5400.0f);
  double gk1p = 1.7f * 0.27f * sqrt(cko / 5400.0f);
  ik1 = (gk1 * y[133] * tempa7 / (tempa7 + tempa8) + gk1p * (1 - y[
    133]) * tempa7 / (tempa7 + tempa8)) * (y[1] - ek);
  //C
  //C Iks    : Delayed rectifier K+ current (Rasmusson et al. 1990]
  //C
  iks = gks * y[27] * y[27] * (y[1] - ek);
  //C
  //C Ikur   : Ultra-rapidly activating delayed rectifier Ikur
  //C                  (Zhou et al., 1998]
  //C
  ikur1 = gkur1 * y[28] * y[29] * (y[1] - ek);
  ikur2 = (gkur2 * y[35] * y[36] * y[132] + gkur2p * y[130] * y[
    131] * (1 - y[132])) * (y[1] - ek);
  ikur3 = gkur3 * y[43] * (y[1] - ek);
  ikur4 = gkur4 * y[184] * y[185] * (y[1] - ek);
  ikur = ikur1 + ikur2 + ikur3 + ikur4;
  //C
  //C Ikr    : HERG current (Wang et al., 1997]
  //C
  double ekr = factor * log((0.98f * cko + 0.02f * cnao) / (
    0.98f * y[24] + 0.02f * y[23]));
  ikr = gkr * y[33] * (y[1] - ekr);
  //C
  //C Iclca  : calcium-activated chloride current (Xu et al., 2002]
  //C
  double tempa9 = poclcamax / (1 + exp((46.7f - y[1]) / 7.8f));
  double tempa10 = y[2] / (y[2] + kmclca);
  iclca = gclca * tempa9 * tempa10 * (y[1] - ecl);
  //C
  //C ICaK : small-conductance calcium-activated K+ current
  //C
  icak = g_cak * (y[181] + y[182]) * (1 - y[183]);
  //C
  //C***********************************************************************
  //C
  //C y[1]  membrane potential
  //C
  double sum_i = icas + icab + inaca + ipca + ina + inab + iclca +
    inak + ikto + ik1 + iks + ikur + ikr - istim + icat + icak;
  ydot[1] = -sum_i / cm;
  //C
  //C y[2]  intracellular calcium Cai
  //C
  bi = 1.0f / (1 + ccmdntot * kmcmdn / ((kmcmdn + y[2]) * (kmcmdn + y[2])));
  double temp1 = 0.5f * (icab - 2.0f * inaca + ipca + icasc + icat) *
    acap / (vmyo * f);
  ydot[2] = bi * (jleak + jxfert - jup - jtrpn - temp1);
  //C
  //C y[3]  network SR calcium Cansr
  //C
  ydot[3] = (jup - jleak) * vmyo / vnsr - jtrt * vjsr / vnsr;
  //C
  //C y[4]  cLTRPNca
  //C
  ydot[4] = kltrpnp * y[2] * (cltrpntot - y[4]) - kltrpnm * y[4];
  //C
  //C y[5]  cHTRPNca
  //C
  ydot[5] = khtrpnp * y[2] * (chtrpntot - y[5]) - khtrpnm * y[5];
  //C
  //C     Partial contributions
  //C
  double alpha = 0.4f * exp((y[1] + 15.0f) / 15.0f);
  double alphap = 0.4f * exp((y[1] + 15.0f + 20.0f) / 15.0f);
  double beta = 0.13f * exp(-1.0f * (y[1] + 15.0f) / 18.0f);
  //C
  double gammac = kpcmax * y[2] / (kpchalf + y[2]);
  //C
  //C y[6] cCajsr junction SR calcium concentartion
  //C
  double temp2 = (kmcsqn + y[6]) * (kmcsqn + y[6]);
  bjsr = 1.0f / (1.0f + ccsqntot * kmcsqn / temp2);
  ydot[6] = bjsr * (jtrt - jrelt);
  //C
  //C y[7] cCass  subspace calcium concentration
  //C
  double temp3 = (kmcmdn + y[7]) * (kmcmdn + y[7]);
  bss = 1.0f / (1.0f + ccmdntot * kmcmdn / temp3);
  double temp4 = jrelt * vjsr / vss - jxfert * vmyo / vss;
  ydot[7] = bss * (temp4 - icase * acap / (2.0f * vss * f));
  //C
  //C               ICaLcav
  //C
  //C y[8]     O  Ca channel variable
  //C
  ydot[8] = kco * y[106] - koc * y[8] + sfoi1 * kpcb * y[13] -
    sfoi1 * gammac * y[8] + sfoi2 * alpha * y[14] - sfoi2 * kpcf * y[
    8] - kicalpka * y[83] * y[8] / (kicalmpka + icalcav * y[8]) + (
    alpha / alphap) * kicalpp * ppcav * y[107] / (kicalmpp +
    icalcav * y[107]);
  //C
  //C y[9]    C1  Ca channel variable
  //C
  //C      ydot[9] = beta*y[10]-4.0*alpha*y[9]
  //C
  //C y[10]   C2  Ca channel variable
  //C
  ydot[10] = 4.0f * alpha * y[9] - beta * y[10] + 2.0f * beta * y[
    11] - 3.0f * alpha * y[10] - kicalpka * y[83] * y[10] / (
    kicalmpka + icalcav * y[10]) + (alphap * alphap * kcop / (alpha *
    alpha * kco)) * kicalpp * ppcav * y[109] / (kicalmpp + icalcav *
    y[109]);
  //C
  //C y[11]   C3  Ca channel variable
  //C
  ydot[11] = 3.0f * alpha * y[10] - 2.0f * beta * y[11] + 3.0f *
    beta * y[12] - 2.0f * alpha * y[11] - kicalpka * y[83] * y[11] / (
    kicalmpka + icalcav * y[11]) + (alphap * kcop / (alpha * kco)) *
    kicalpp * ppcav * y[110] / (kicalmpp + icalcav * y[110]);
  //C
  //C y[12]   C4  Ca channel variable
  //C
  ydot[12] = 2.0f * alpha * y[11] - 3.0f * beta * y[12] + 4.0f *
    beta * y[106] - alpha * y[12] + 4.0f * sfc4i1 * kpcb * beta * y[
    13] - sfc4i1 * alpha * gammac * (kco / koc) * y[12] + 4.0f *
    sfc4i2 * beta * y[14] - sfc4i2 * kpcf * (kco / koc) * y[12] +
    4.0f * sfc4i3 * beta * kpcb * y[15] - sfc4i3 * gammac * kpcf * (
    kco / koc) * y[12] - kicalpka * y[83] * y[12] / (kicalmpka +
    icalcav * y[12]) + (kcop / kco) * kicalpp * ppcav * y[111] / (
    kicalmpp + icalcav * y[111]);
  //C
  //C y[106] CP Ca channel variable
  //C
  ydot[106] = alpha * y[12] - 4.0f * beta * y[106] + koc * y[8] -
    kco * y[106] - kicalpka * y[83] * y[106] / (kicalmpka + icalcav *
    y[106]) + (alpha * kcop / (alphap * kco)) * kicalpp * ppcav * y[
    115] / (kicalmpp + icalcav * y[115]);
  //C
  //C y[13]   I1  Ca channel variable
  //C
  ydot[13] = sfoi1 * gammac * y[8] - sfoi1 * kpcb * y[13] + sfi1i3 *
    alpha * y[15] - sfi1i3 * kpcf * y[13] + sfc4i1 * alpha * gammac *
    (kco / koc) * y[12] - 4.0f * sfc4i1 * beta * kpcb * y[13] -
    kicalpka * y[83] * y[13] / (kicalmpka + icalcav * y[13]) + (
    alpha / alphap) * kicalpp * ppcav * y[112] / (kicalmpp +
    icalcav * y[112]);
  //C
  //C y[14]   I2  Ca channel variable
  //C
  ydot[14] = sfoi2 * kpcf * y[8] - sfoi2 * alpha * y[14] + sfi2i3 *
    kpcb * y[15] - sfi2i3 * gammac * y[14] + sfc4i2 * kpcf * (kco /
    koc) * y[12] - 4.0f * sfc4i2 * beta * y[14] - kicalpka * y[83] * y[
    14] / (kicalmpka + icalcav * y[14]) + kicalpp * ppcav * y[113] / (
    kicalmpp + icalcav * y[113]);
  //C
  //C y[15]   I3  Ca channel variable
  //C
  ydot[15] = sfi1i3 * kpcf * y[13] - sfi1i3 * alpha * y[15] +
    sfi2i3 * gammac * y[14] - sfi2i3 * kpcb * y[15] + sfc4i3 *
    gammac * kpcf * (kco / koc) * y[12] - 4.0f * sfc4i3 * beta *
    kpcb * y[15] - kicalpka * y[83] * y[15] / (kicalmpka + icalcav *
    y[15]) + kicalpp * ppcav * y[114] / (kicalmpp + icalcav * y[114]);
  //C
  //C y[107] OP Ca channel variable
  //C
  ydot[107] = kcop * y[115] - koc * y[107] + sfoi1 * kpcb * y[112] -
    sfoi1 * gammac * y[107] + sfoi2 * alphap * y[113] - sfoi2 *
    kpcf * y[107] + kicalpka * y[83] * y[8] / (kicalmpka + icalcav *
    y[8]) - (alpha / alphap) * kicalpp * ppcav * y[107] / (kicalmpp +
    icalcav * y[107]);
  //C
  //C y[108] C1P Ca channel variable
  //C
  ydot[108] = beta * y[109] - 4.0f * alphap * y[108] + kicalpka * y[
    83] * y[9] / (kicalmpka + icalcav * y[9]) - (alphap * alphap *
    alphap * kcop / (alpha * alpha * alpha * kco)) * kicalpp *
    ppcav * y[108] / (kicalmpp + icalcav * y[108]);
  //C
  //C y[109] C2P Ca channel variable
  //C
  ydot[109] = 4.0f * alphap * y[108] - beta * y[109] + 2.0f * beta *
    y[110] - 3.0f * alphap * y[109] + kicalpka * y[83] * y[10] / (
    kicalmpka + icalcav * y[10]) - (alphap * alphap * kcop / (alpha *
    alpha * kco)) * kicalpp * ppcav * y[109] / (kicalmpp + icalcav *
    y[109]);
  //C
  //C y[110] C3P Ca channel variable
  //C
  ydot[110] = 3.0f * alphap * y[109] - 2.0f * beta * y[110] + 3.0f *
    beta * y[111] - 2.0f * alphap * y[110] + kicalpka * y[83] * y[
    11] / (kicalmpka + icalcav * y[11]) - (alphap * kcop / (alpha *
    kco)) * kicalpp * ppcav * y[110] / (kicalmpp + icalcav * y[110]);
  //C
  //C y[111] C4P Ca channel variable
  //C
  ydot[111] = 2.0f * alphap * y[110] - 3.0f * beta * y[111] + 4.0f *
    beta * y[115] - alphap * y[111] + 4.0f * sfc4i1 * kpcb * beta * y[
    112] - sfc4i1 * alphap * gammac * (kcop / koc) * y[111] + 4.0f *
    sfc4i2 * beta * y[113] - sfc4i2 * kpcf * (kcop / koc) * y[111] +
    4.0f * sfc4i3 * beta * kpcb * y[114] - sfc4i3 * gammac * kpcf * (
    kcop / koc) * y[111] + kicalpka * y[83] * y[12] / (kicalmpka +
    icalcav * y[12]) - (kcop / kco) * kicalpp * ppcav * y[111] / (
    kicalmpp + icalcav * y[111]);
  //C
  //C y[115] CPP Ca channel variable
  //C
  ydot[115] = alphap * y[111] - 4.0f * beta * y[115] + koc * y[107] -
    kcop * y[115] + kicalpka * y[83] * y[106] / (kicalmpka +
    icalcav * y[106]) - (alpha * kcop / (alphap * kco)) * kicalpp *
    ppcav * y[115] / (kicalmpp + icalcav * y[115]);
  //C
  //C y[112] I1P Ca channel variable
  //C
  ydot[112] = sfoi1 * gammac * y[107] - sfoi1 * kpcb * y[112] +
    sfi1i3 * alphap * y[114] - sfi1i3 * kpcf * y[112] + sfc4i1 *
    alphap * gammac * (kcop / koc) * y[111] - 4.0f * sfc4i1 * beta *
    kpcb * y[112] + kicalpka * y[83] * y[13] / (kicalmpka + icalcav *
    y[13]) - (alpha / alphap) * kicalpp * ppcav * y[112] / (
    kicalmpp + icalcav * y[112]);
  //C
  //C y[113] I2P Ca channel variable
  //C
  ydot[113] = sfoi2 * kpcf * y[107] - sfoi2 * alphap * y[113] + sfi2i3 *
    kpcb * y[114] - sfi2i3 * gammac * y[113] + sfc4i2 * kpcf * (kcop /
    koc) * y[111] - 4.0f * sfc4i2 * beta * y[113] + kicalpka * y[83] * y[
    14] / (kicalmpka + icalcav * y[14]) - kicalpp * ppcav * y[113] / (
    kicalmpp + icalcav * y[113]);
  //C
  //C y[114] I3P Ca channel variable
  //C
  ydot[114] = sfi1i3 * kpcf * y[112] - sfi1i3 * alphap * y[114] +
    sfi2i3 * gammac * y[113] - sfi2i3 * kpcb * y[114] + sfc4i3 *
    gammac * kpcf * (kcop / koc) * y[111] - 4.0f * sfc4i3 * beta *
    kpcb * y[114] + kicalpka * y[83] * y[15] / (kicalmpka + icalcav *
    y[15]) - kicalpp * ppcav * y[114] / (kicalmpp + icalcav * y[114]);
  //C
  //C     Partial contributions
  //C
  double gammae = kpcmax * y[7] / (kpchalf + y[7]);
  //C
  //C                ICaLecav
  //C
  //C y[137]   O  Ca channel variable
  //C
  ydot[137] = kco * y[142] - koc * y[137] + sfoi1 * kpcb * y[143] -
    sfoi1 * gammae * y[137] + sfoi2 * alpha * y[144] - sfoi2 * kpcf *
    y[137] - kicalpka * y[89] * y[137] / (kicalmpka + icalecav * y[
    137]) + (alpha / alphap) * kicalpp * pp1ecav * y[146] / (
    kicalmpp + icalecav * y[146]);
  //C
  //C y[138]   C1  Ca channel variable
  //C
  //C      ydot[138] = beta*y[139]-4.0*alpha*y[138]
  //C
  //C y[139]   C2  Ca channel variable
  //C
  ydot[139] = 4.0f * alpha * y[138] - beta * y[139] + 2.0f * beta * y[
    140] - 3.0f * alpha * y[139] - kicalpka * y[89] * y[139] / (
    kicalmpka + icalecav * y[139]) + (alphap * alphap * kcop / (
    alpha * alpha * kco)) * kicalpp * pp1ecav * y[148] / (kicalmpp +
    icalecav * y[148]);
  //C
  //C y[140]   C3  Ca channel variable
  //C
  ydot[140] = 3.0f * alpha * y[139] - 2.0f * beta * y[140] + 3.0f *
    beta * y[141] - 2.0f * alpha * y[140] - kicalpka * y[89] * y[
    140] / (kicalmpka + icalecav * y[140]) + (alphap * kcop / (
    alpha * kco)) * kicalpp * pp1ecav * y[149] / (kicalmpp +
    icalecav * y[149]);
  //C
  //C y[141]   C4  Ca channel variable
  //C
  ydot[141] = 2.0f * alpha * y[140] - 3.0f * beta * y[141] + 4.0f *
    beta * y[142] - alpha * y[141] + 4.0f * sfc4i1 * kpcb * beta * y[
    143] - sfc4i1 * alpha * gammae * (kco / koc) * y[141] + 4.0f *
    sfc4i2 * beta * y[144] - sfc4i2 * kpcf * (kco / koc) * y[141] +
    4.0f * sfc4i3 * beta * kpcb * y[145] - sfc4i3 * gammae * kpcf * (
    kco / koc) * y[141] - kicalpka * y[89] * y[141] / (kicalmpka +
    icalecav * y[141]) + (kcop / kco) * kicalpp * pp1ecav * y[150] / (
    kicalmpp + icalecav * y[150]);
  //C
  //C y[142] CP Ca channel variable
  //C
  ydot[142] = alpha * y[141] - 4.0f * beta * y[142] + koc * y[137] -
    kco * y[142] - kicalpka * y[89] * y[142] / (kicalmpka +
    icalecav * y[142]) + (alpha * kcop / (alphap * kco)) * kicalpp *
    pp1ecav * y[151] / (kicalmpp + icalecav * y[151]);
  //C
  //C y[143]   I1  Ca channel variable
  //C
  ydot[143] = sfoi1 * gammae * y[137] - sfoi1 * kpcb * y[143] +
    sfi1i3 * alpha * y[145] - sfi1i3 * kpcf * y[143] + sfc4i1 *
    alpha * gammae * (kco / koc) * y[141] - 4.0f * sfc4i1 * beta *
    kpcb * y[143] - kicalpka * y[89] * y[143] / (kicalmpka +
    icalecav * y[143]) + (alpha / alphap) * kicalpp * pp1ecav * y[
    152] / (kicalmpp + icalecav * y[152]);
  //C
  //C y[144]   I2  Ca channel variable
  //C
  ydot[144] = sfoi2 * kpcf * y[137] - sfoi2 * alpha * y[144] +
    sfi2i3 * kpcb * y[145] - sfi2i3 * gammae * y[144] + sfc4i2 *
    kpcf * (kco / koc) * y[141] - 4.0f * sfc4i2 * beta * y[144] -
    kicalpka * y[89] * y[144] / (kicalmpka + icalecav * y[144]) +
    kicalpp * pp1ecav * y[153] / (kicalmpp + icalecav * y[153]);
  //C
  //C y[145]   I3  Ca channel variable
  //C
  ydot[145] = sfi1i3 * kpcf * y[143] - sfi1i3 * alpha * y[145] +
    sfi2i3 * gammae * y[144] - sfi2i3 * kpcb * y[145] + sfc4i3 *
    gammae * kpcf * (kco / koc) * y[141] - 4.0f * sfc4i3 * beta *
    kpcb * y[145] - kicalpka * y[89] * y[145] / (kicalmpka +
    icalecav * y[145]) + kicalpp * pp1ecav * y[154] / (kicalmpp +
    icalecav * y[154]);
  //C
  //C y[146] OP Ca channel variable
  //C
  ydot[146] = kcop * y[151] - koc * y[146] + sfoi1 * kpcb * y[152] -
    sfoi1 * gammae * y[146] + sfoi2 * alphap * y[153] - sfoi2 * kpcf *
    y[146] + kicalpka * y[89] * y[137] / (kicalmpka + icalecav * y[
    137]) - (alpha / alphap) * kicalpp * pp1ecav * y[146] / (
    kicalmpp + icalecav * y[146]);
  //C
  //C y[147] C1P Ca channel variable
  //C
  ydot[147] = beta * y[148] - 4.0f * alphap * y[147] + kicalpka * y[
    89] * y[138] / (kicalmpka + icalecav * y[138]) - (alphap *
    alphap * alphap * kcop / (alpha * alpha * alpha * kco)) *
    kicalpp * pp1ecav * y[147] / (kicalmpp + icalecav * y[147]);
  //C
  //C y[148] C2P Ca channel variable
  //C
  ydot[148] = 4.0f * alphap * y[147] - beta * y[148] + 2.0f * beta *
    y[149] - 3.0f * alphap * y[148] + kicalpka * y[89] * y[139] / (
    kicalmpka + icalecav * y[139]) - (alphap * alphap * kcop / (
    alpha * alpha * kco)) * kicalpp * pp1ecav * y[148] / (kicalmpp +
    icalecav * y[148]);
  //C
  //C y[149] C3P Ca channel variable
  //C
  ydot[149] = 3.0f * alphap * y[148] - 2.0f * beta * y[149] + 3.0f *
    beta * y[150] - 2.0f * alphap * y[149] + kicalpka * y[89] * y[
    140] / (kicalmpka + icalecav * y[140]) - (alphap * kcop / (
    alpha * kco)) * kicalpp * pp1ecav * y[149] / (kicalmpp +
    icalecav * y[149]);
  //C
  //C y[150] C4P Ca channel variable
  //C
  ydot[150] = 2.0f * alphap * y[149] - 3.0f * beta * y[150] + 4.0f *
    beta * y[151] - alphap * y[150] + 4.0f * sfc4i1 * kpcb * beta * y[
    152] - sfc4i1 * alphap * gammae * (kcop / koc) * y[150] + 4.0f *
    sfc4i2 * beta * y[153] - sfc4i2 * kpcf * (kcop / koc) * y[150] +
    4.0f * sfc4i3 * beta * kpcb * y[154] - sfc4i3 * gammae * kpcf * (
    kcop / koc) * y[150] + kicalpka * y[89] * y[141] / (kicalmpka +
    icalecav * y[141]) - (kcop / kco) * kicalpp * pp1ecav * y[150] / (
    kicalmpp + icalecav * y[150]);
  //C
  //C y[151] CPP Ca channel variable
  //C
  ydot[151] = alphap * y[150] - 4.0f * beta * y[151] + koc * y[146] -
    kcop * y[151] + kicalpka * y[89] * y[142] / (kicalmpka +
    icalecav * y[142]) - (alpha * kcop / (alphap * kco)) * kicalpp *
    pp1ecav * y[151] / (kicalmpp + icalecav * y[151]);
  //C
  //C y[152] I1P Ca channel variable
  //C
  ydot[152] = sfoi1 * gammae * y[146] - sfoi1 * kpcb * y[152] +
    sfi1i3 * alphap * y[154] - sfi1i3 * kpcf * y[152] + sfc4i1 *
    alphap * gammae * (kcop / koc) * y[150] - 4.0f * sfc4i1 * beta *
    kpcb * y[152] + kicalpka * y[89] * y[143] / (kicalmpka +
    icalecav * y[143]) - (alpha / alphap) * kicalpp * pp1ecav * y[
    152] / (kicalmpp + icalecav * y[152]);
  //C
  //C y[153] I2P Ca channel variable
  //C
  ydot[153] = sfoi2 * kpcf * y[146] - sfoi2 * alphap * y[153] +
    sfi2i3 * kpcb * y[154] - sfi2i3 * gammae * y[153] + sfc4i2 *
    kpcf * (kcop / koc) * y[150] - 4.0f * sfc4i2 * beta * y[153] +
    kicalpka * y[89] * y[144] / (kicalmpka + icalecav * y[144]) -
    kicalpp * pp1ecav * y[153] / (kicalmpp + icalecav * y[153]);
  //C
  //C y[154] I3P Ca channel variable
  //C
  ydot[154] = sfi1i3 * kpcf * y[152] - sfi1i3 * alphap * y[154] +
    sfi2i3 * gammae * y[153] - sfi2i3 * kpcb * y[154] + sfc4i3 *
    gammae * kpcf * (kcop / koc) * y[150] - 4.0f * sfc4i3 * beta *
    kpcb * y[154] + kicalpka * y[89] * y[145] / (kicalmpka +
    icalecav * y[145]) - kicalpp * pp1ecav * y[154] / (kicalmpp +
    icalecav * y[154]);
  //C
  //C y[16]-y[19]  RyR channel states
  //C
  double temp5 = y[7] * y[7] * y[7];
  double temp6 = y[7] * y[7] * y[7] * y[7];
  //C
  //C y[16] RyR C1 channel variable
  //C
  ydot[16] = -kap * temp6 * y[16] + kam * y[18] - kryrpka * y[89] * y[
    16] / (kryrmpka + ryrecav * y[16]) + kryrpp * pp1ecav * y[116] / (
    kryrmpp + ryrecav * y[116]);
  //C
  //C y[17] RyR C2 channel variable
  //C
  ydot[17] = kcp * y[18] - kcm * y[17] - kryrpka * y[89] * y[17] * f_ryr / (
    kryrmpka + ryrecav * y[17]) + ((kap * kamp * kcp * kcmp) / (kapp *
    kam * kcpp * kcm)) * f_ryr * kryrpp * pp1ecav * y[117] / (kryrmpp +
    ryrecav * y[117]);
  //C
  //C y[18] RyR O1 channel variable
  //C
  ydot[18] = kap * temp6 * y[16] - kam * y[18] - kbp * temp5 * y[
    18] + kbm * y[19] - kcp * y[18] + kcm * y[17] - kryrpka * y[89] *
    y[18] * f_ryr / (kryrmpka + ryrecav * y[18]) + ((kap * kamp) / (
    kapp * kam)) * f_ryr * kryrpp * pp1ecav * y[118] / (kryrmpp +
    ryrecav * y[118]);
  //C
  //C y[19] RyR O2 channel variable
  //C
  ydot[19] = kbp * temp5 * y[18] - kbm * y[19] - kryrpka * y[89] * y[
    19] * f_ryr / (kryrmpka + ryrecav * y[19]) + ((kap * kamp * kbp *
    kbmp) / (kapp * kam * kbpp * kbm)) * f_ryr * kryrpp * pp1ecav * y[
    119] / (kryrmpp + ryrecav * y[119]);
  //C
  //C y[116] RyR C1P channel variable
  //C
  ydot[116] = -kapp * temp6 * y[116] + kamp * y[118] + kryrpka * y[
    89] * y[16] / (kryrmpka + ryrecav * y[16]) - kryrpp * pp1ecav * y[
    116] / (kryrmpp + ryrecav * y[116]);
  //C
  //C y[117] RyR C2P channel variable
  //C
  ydot[117] = kcpp * y[118] - kcmp * y[117] + kryrpka * y[89] * y[
    17] * f_ryr / (kryrmpka + ryrecav * y[17]) - ((kap * kamp * kcp *
    kcmp) / (kapp * kam * kcpp * kcm)) * f_ryr * kryrpp * pp1ecav * y[
    117] / (kryrmpp + ryrecav * y[117]);
  //C
  //C y[118] RyR O1P channel variable
  //C
  ydot[118] = kapp * temp6 * y[116] - kamp * y[118] - kbpp * temp5 * y[
    118] + kbmp * y[119] - kcpp * y[118] + kcmp * y[117] + kryrpka * y[
    89] * y[18] * f_ryr / (kryrmpka + ryrecav * y[18]) - ((kap *
    kamp) / (kapp * kam)) * f_ryr * kryrpp * pp1ecav * y[118] / (
    kryrmpp + ryrecav * y[118]);
  //C
  //C y[119] RyR O2P channel variable
  //C
  ydot[119] = kbpp * temp5 * y[118] - kbmp * y[119] + kryrpka * y[
    89] * y[19] * f_ryr / (kryrmpka + ryrecav * y[19]) - ((kap *
    kamp * kbp * kbmp) / (kapp * kam * kbpp * kbm)) * f_ryr *
    kryrpp * pp1ecav * y[119] / (kryrmpp + ryrecav * y[119]);
  //C
  //C  y[20]-y[22],y[37]-y[42] Na fast current (Clancy-Rudy, Circulation, 2002]
  //C
  double va = y[1] - 2.5f;
  double vi = y[1] + 7.0f;
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
  double vap = y[1] - 2.5f;
  double vip = y[1] + 7.0f;
  double alp11p = 3.802f / (0.1027f * exp(-vap / 17.0f) + 0.20f *
    exp(-vap / 150.f));
  double alp12p = 3.802f / (0.1027f * exp(-vap / 15.0f) + 0.23f *
    exp(-vap / 150.f));
  double alp13p = 3.802f / (0.1027f * exp(-vap / 12.0f) + 0.25f *
    exp(-vap / 150.f));
  double bet2p = alp13p * alp2 * alp3 / (bet13 * bet3);
  //C
  //C y[20] INa CNa3 channel variable
  //C
  //C      ydot[20] = bet11*y[21] - alp11*y[20]
  //C     *         + alp3*y[42]  - bet3*y[20]
  //C     * - kINAPKA*y[83]*y[20]/(KINamPKA+y[20])
  //C     * + (alp11p*alp12p/(alp11*alp12])*kINaPP*PPcav*y[120]
  //C     * /(KINamPP+y[120])
  //C
  //C y[21] INa CNa2 channel variable
  //C
  ydot[21] = alp11 * y[20] - bet11 * y[21] + bet12 * y[22] - alp12 * y[
    21] + alp3 * y[41] - bet3 * y[21] - kinapka * y[83] * y[21] / (
    kinampka + y[21]) + (alp12p / alp12) * kinapp * ppcav * y[121] / (
    kinampp + y[121]);
  //C
  //C y[22] INa CNa1 channel variable
  //C
  ydot[22] = alp12 * y[21] - bet12 * y[22] + bet13 * y[37] - alp13 *
    y[22] + alp3 * y[38] - bet3 * y[22] - kinapka * y[83] * y[22] / (
    kinampka + y[22]) + kinapp * ppcav * y[122] / (kinampp + y[122]);
  //C
  //C y[37] INa ONa channel variable
  //C
  ydot[37] = alp13 * y[22] - bet13 * y[37] + bet2 * y[38] - alp2 * y[
    37] - kinapka * y[83] * y[37] / (kinampka + y[37]) + (alp13 /
    alp13p) * kinapp * ppcav * y[123] / (kinampp + y[123]);
  //C
  //C y[38] INa IFNa channel variable
  //C
  ydot[38] = alp2 * y[37] - bet2 * y[38] + bet3 * y[22] - alp3 * y[38] +
    bet4 * y[39] - alp4 * y[38] + alp12 * y[41] - bet12 * y[38] -
    kinapka * y[83] * y[38] / (kinampka + y[38]) + kinapp * ppcav * y[
    124] / (kinampp + y[124]);
  //C
  //C y[39] INa I1Na channel variable
  //C
  ydot[39] = alp4 * y[38] - bet4 * y[39] + bet5 * y[40] - alp5 * y[
    39] - kinapka * y[83] * y[39] / (kinampka + y[39]) + kinapp *
    ppcav * y[125] / (kinampp + y[125]);
  //C
  //C y[40] INa I2Na channel variable
  //C
  ydot[40] = alp5 * y[39] - bet5 * y[40] - kinapka * y[83] * y[40] / (
    kinampka + y[40]) + kinapp * ppcav * y[126] / (kinampp + y[126]);
  //C
  //C y[41] INa ICNa2 channel variable
  //C
  ydot[41] = alp11 * y[42] - bet11 * y[41] + bet12 * y[38] - alp12 * y[
    41] + bet3 * y[21] - alp3 * y[41] - kinapka * y[83] * y[41] / (
    kinampka + y[41]) + (alp12p / alp12) * kinapp * ppcav * y[127] / (
    kinampp + y[127]);
  //C
  //C y[42] INa ICNa3 channel variable
  //C
  ydot[42] = bet11 * y[41] - alp11 * y[42] + bet3 * y[20] - alp3 * y[
    42] - kinapka * y[83] * y[42] / (kinampka + y[42]) + (alp11p *
    alp12p / (alp11 * alp12)) * kinapp * ppcav * y[128] / (kinampp +
    y[128]);
  //C
  //C y[120] INa CNa3p channel variable
  //C
  ydot[120] = bet11 * y[121] - alp11p * y[120] + alp3 * y[128] -
    bet3 * y[120] + kinapka * y[83] * y[20] / (kinampka + y[20]) - (
    alp11p * alp12p / (alp11 * alp12)) * kinapp * ppcav * y[120] / (
    kinampp + y[120]);
  //C
  //C y[121] INa CNa2p channel variable
  //C
  ydot[121] = alp11p * y[120] - bet11 * y[121] + bet12 * y[122] - alp12p *
    y[121] + alp3 * y[127] - bet3 * y[121] + kinapka * y[83] * y[21] / (
    kinampka + y[21]) - (alp12p / alp12) * kinapp * ppcav * y[121] / (
    kinampp + y[121]);
  //C
  //C y[122] INa CNa1p channel variable
  //C
  ydot[122] = alp12p * y[121] - bet12 * y[122] + bet13 * y[123] -
    alp13p * y[122] + alp3 * y[124] - bet3 * y[122] + kinapka * y[
    83] * y[22] / (kinampka + y[22]) - kinapp * ppcav * y[122] / (
    kinampp + y[122]);
  //C
  //C y[123] INa ONap channel variable
  //C
  ydot[123] = alp13p * y[122] - bet13 * y[123] + bet2p * y[124] -
    alp2 * y[123] + kinapka * y[83] * y[37] / (kinampka + y[37]) - (
    alp13 / alp13p) * kinapp * ppcav * y[123] / (kinampp + y[123]);
  //C
  //C y[124] INa IFNap channel variable
  //C
  ydot[124] = alp2 * y[123] - bet2p * y[124] + bet3 * y[122] - alp3 *
    y[124] + bet4 * y[125] - alp4 * y[124] + alp12p * y[127] -
    bet12 * y[124] + kinapka * y[83] * y[38] / (kinampka + y[38]) -
    kinapp * ppcav * y[124] / (kinampp + y[124]);
  //C
  //C y[125] INa I1Nap channel variable
  //C
  ydot[125] = alp4 * y[124] - bet4 * y[125] + bet5 * y[126] - alp5 *
    y[125] + kinapka * y[83] * y[39] / (kinampka + y[39]) - kinapp *
    ppcav * y[125] / (kinampp + y[125]);
  //C
  //C y[126] INa I2Nap channel variable
  //C
  ydot[126] = alp5 * y[125] - bet5 * y[126] + kinapka * y[83] * y[
    40] / (kinampka + y[40]) - kinapp * ppcav * y[126] / (kinampp + y[
    126]);
  //C
  //C y[127] INa ICNa2p channel variable
  //C
  ydot[127] = alp11p * y[128] - bet11 * y[127] + bet12 * y[124] - alp12p *
    y[127] + bet3 * y[121] - alp3 * y[127] + kinapka * y[83] * y[41] / (
    kinampka + y[41]) - (alp12p / alp12) * kinapp * ppcav * y[127] / (
    kinampp + y[127]);
  //C
  //C y[128] INa ICNa3p channel variable
  //C
  ydot[128] = bet11 * y[127] - alp11p * y[128] + bet3 * y[120] -
    alp3 * y[128] + kinapka * y[83] * y[42] / (kinampka + y[42]) - (
    alp11p * alp12p / (alp11 * alp12)) * kinapp * ppcav * y[128] / (
    kinampp + y[128]);
  //C
  //C  y[23]  Na intracellular concentration
  //C
  ydot[23] = -1.0f * (ina + inab + 3.0f * (inaca + inak)) * acap / (vmyo * f);
  //C
  //C  y[24]  K  intracellular concentration
  //C
  ydot[24] = -1.0f * (ikto + ik1 + iks + ikur + ikr + icak - 2.0f *
    inak - istim) * acap / (vmyo * f);
  //C
  //C  y[25],y[26]  ato and ito gating variables for Ikto
  //C
  double alp25 = 0.04516f * exp(0.03577f * (y[1] + 33.0f)) * 4.0f;
  double bet25 = 0.0989f * exp(-0.06237f * (y[1] + 33.0f)) * 4.0f;
  double temp7 = 0.0019f * exp((y[1] + 15.5f) / (-1.0f * dito));
  double temp8 = 0.067083f * exp((y[1] + 15.5f + 20.0f) / (-1.0f * dito));
  double alp26 = 0.08f * temp7 / (1.0f + temp8);
  double temp9 = 0.0019f * exp((y[1] + 15.5f + 20.0f) / dito);
  double temp10 = 0.051335f * exp((y[1] + 15.5f + 20.0f) / dito);
  double bet26 = 0.5f * temp9 / (1.0f + temp10);
  //C
  ydot[25] = alp25 * (1.0f - y[25]) - bet25 * y[25];
  ydot[26] = alp26 * (1.0f - y[26]) - bet26 * y[26];
  //C
  //C  y[135],y[136] ato and ito gating variables for Ikto,phosphorylated
  //C
  double alp25p = 0.04516f * exp(0.03577f * (y[1] + 30.0f - 13.0f)) * 4.0f;
  double bet25p = 0.0989f * exp(-0.06237f * (y[1] + 30.0f - 13.0f)) * 4.0f;
  double temp7p = 0.0019f * exp((y[1] + 13.5f - 6.0f) / (-1.0f * dito));
  double temp8p = 0.067083f * exp((y[1] + 13.5f + 20.0f - 6.0f) /
    (-1.0f * dito));
  double alp26p = 0.08f * temp7p / (1.0f + temp8p);
  double temp9p = 0.0019f * exp((y[1] + 13.5f + 20.0f - 6.0f) / dito);
  double temp10p = 0.051335f * exp((y[1] + 13.5f + 20.0f - 6.0f) / dito);
  double bet26p = 0.5f * temp9p / (1.0f + temp10p);
  //C
  ydot[135] = alp25p * (1.0f - y[135]) - bet25p * y[135];
  ydot[136] = alp26p * (1.0f - y[136]) - bet26p * y[136];
  //C
  //C  y[27]  nks gating variable for Iks
  //C
  double temp11 = 0.00001444f * (y[1] + 26.5f);
  double alp27 = temp11 / (1.0f - exp(-0.128f * (y[1] + 26.5f))) / 3.f;
  double bet27 = 0.000286f * exp(-0.038f * (y[1] + 26.5f)) / 3.f;
  //C
  ydot[27] = alp27 * (1.0f - y[27]) - bet27 * y[27];
  //C
  //C  y[28], y[29]  aur and iur gating variables for Ikur1
  //C
  double ass = 1.0f / (1.0f + exp((-22.5f - y[1]) / 7.7f));
  double iss1 = 1.0f / (1.0f + exp((45.2f + y[1]) / 5.7f));
  //C      taua1 = [0.493*exp(-0.0629*y[1])+2.058]
  double taua1 = 6.1f / (exp(0.0629f * (y[1] + 40.0f)) + exp(
    -0.0629f * (y[1] + 40.0f))) + 2.058f;
  double taui1 = 270.f + 1050 / (1.0f + exp((45.2f + y[1]) / 5.7f));
  //C
  ydot[28] = (ass - y[28]) / taua1;
  ydot[29] = (iss1 - y[29]) / taui1;
  //C
  //C  y[35], y[36]  aur and iur gating variables for Ikur2
  //C
  double iss2 = 1.0f / (1.0f + exp((45.2f + y[1]) / 5.7f));
  //C      taua2 = [0.493*exp(-0.0629*y[1])+2.058]
  double taua2 = 6.1f / (exp(0.0629f * (y[1] + 40.0f)) + exp(
    -0.0629f * (y[1] + 40.0f))) + 2.058f;
  double taui2 = 803.f - 18.f / (1.0f + exp((45.2f + y[1]) / 5.7f));
  //C
  ydot[35] = (ass - y[35]) / taua2;
  ydot[36] = (iss2 - y[36]) / taui2;
  ydot[130] = (ass - y[130]) / taua2;
  ydot[131] = (iss2 - y[131]) / taui2;
  ydot[132] = kikurpp * pp1ecav * (1 - y[132]) / (kikurmpp + (1 - y[
    132])) - kikurpka * y[89] * y[132] / (kikurmpka + y[132]);
  //C
  //C  y[43] aur gating variable for Ikur3
  //C
  //C      taua3 = [39.3*exp(-0.0862*y[1])+13.17]
  double taua3 = 1235.5f / (exp(0.0862f * (y[1] + 40.0f)) +
    exp(-0.0862f * (y[1] + 40.0f))) + 13.17f;
  ydot[43] = (ass - y[43]) / taua3;
  //C
  //C  y[184], y[185] aur and iur gating variables for Ikur4
  //C
  double taui4 = 5334.f - 4912.f / (1.0f + exp((45.2f + y[1]) / 5.7f));
  ydot[184] = (ass - y[184]) / taua2;
  ydot[185] = (iss2 - y[185]) / taui4;
  //C
  //C y[30]-y[34] HERG channel state variables
  //C
  double ala0 = 0.022348f * exp(0.01176f * y[1]);
  double bea0 = 0.047002f * exp(-0.0631f * y[1]);
  double ala1 = 0.013733f * exp(0.038198f * y[1]);
  double bea1 = 0.0000689f * exp(-0.04178f * y[1]);
  double ali = 0.090821f * exp(0.023391f * (y[1] + 5.0f));
  double bei = 0.006497f * exp(-0.03268f * (y[1] + 5.0f));
  //C
  ydot[30] = bea0 * y[31] - ala0 * y[30];
  ydot[31] = ala0 * y[30] - bea0 * y[31] + kb * y[32] - kf * y[31];
  ydot[32] = kf * y[31] - kb * y[32] + bea1 * y[33] - ala1 * y[32];
  ydot[33] = ala1 * y[32] - bea1 * y[33] + bei * y[34] - ali * y[33];
  ydot[34] = ali * y[33] - bei * y[34];
  //C
  //C   y[44] Pryr Ryanodine receptor modulation factor
  //C
  ydot[44] = -t1 * y[44] + t2 * icas * exp((y[1] + 5.0f) * (y[
    1] + 5.0f) / (-648.f));
  //C
  //C   Caveolae domain
  //C
  //C   y[45] beta1 tot concentration phosphorylated by PKA caveolae
  //C
  ydot[45] = kpkap * y[83] * rb1cavnptot - kpkam * y[45];
  //C
  //C   y[46] beta1 tot concentration phosphorylated by BARK caveolae
  //C
  ydot[46] = kbarkp * (lrb1cavnp + lrb1gscavnp) - kbarkm * y[46];
  //C
  //C   y[155] beta2 tot concentration phosphorylated by PKA caveolae
  //C
  ydot[155] = kpkap * y[83] * rb2cavnptot - kpkam * y[155];
  //C
  //C   y[156] beta2 tot concentration phosphorylated by BARK caveolae
  //C
  ydot[156] = kbarkp * (lrb2cavnp + lrb2gscavnp) - kbarkm * y[156];
  //C
  //C   y[47] Gs-alpha with GTP caveolae
  //C
  ydot[47] = kact2gs * rb1gscavnp + factgsgi * kact2gs * rb2gscavnp +
    kact1gs * lrb1gscavnp + factgsgi * kact1gs * lrb2gscavnp -
    khydgs * y[47];
  //C
  //C   y[48]  G Beta-Gamma caveolae
  //C
  ydot[48] = kact2gs * rb1gscavnp + factgsgi * kact2gs * rb2gscavnp +
    kact1gs * lrb1gscavnp + factgsgi * kact1gs * lrb2gscavnp + kact2gi *
    rb2gicavpka + kact1gi * lrb2gicavpka - kreasgs * y[48] * y[49] -
    kreasgi * y[48] * y[158];
  //C
  //C   y[49]  Gs-alpha with GDP caveolae
  //C
  ydot[49] = khydgs * y[47] - kreasgs * y[48] * y[49];
  //C
  //C   y[157] Gi-alpha with GTP caveolae
  //C
  ydot[157] = kact2gi * rb2gicavpka + kact1gi * lrb2gicavpka - khydgi * y[157];
  //C
  //C   y[158]  Gi-alpha with GDP caveolae
  //C
  ydot[158] = khydgi * y[157] - kreasgi * y[48] * y[158];
  //C
  //C   y[50] beta1 tot concentration phosphorylated by PKA extracaveolae
  //C
  ydot[50] = kpkap * y[89] * rb1ecavnptot - kpkam * y[50];
  //C
  //C   y[51] beta1 tot concentration phosphorylated by BARK extracaveolae
  //C
  ydot[51] = kbarkp * (lrb1ecavnp + lrb1gsecavnp) - kbarkm * y[51];
  //C
  //C   y[159] beta2 tot concentration phosphorylated by PKA extracaveolae
  //C
  ydot[159] = kpkap * y[89] * rb2ecavnptot - kpkam * y[159];
  //C
  //C   y[160] beta2 tot concentration phosphorylated by BARK extracaveolae
  //C
  ydot[160] = kbarkp * (lrb2ecavnp + lrb2gsecavnp) - kbarkm * y[160];
  //C
  //C   y[52] Gs-alpha with GTP extracaveolae
  //C
  ydot[52] = kact2gs * rb1gsecavnp + factgsgi * kact2gs * rb2gsecavnp +
    kact1gs * lrb1gsecavnp + factgsgi * kact1gs * lrb2gsecavnp -
    khydgs * y[52];
  //C
  //C   y[53]  G Beta-Gamma extracaveolae
  //C
  ydot[53] = kact2gs * rb1gsecavnp + factgsgi * kact2gs *
    rb2gsecavnp + kact1gs * lrb1gsecavnp + factgsgi * kact1gs *
    lrb2gsecavnp + kact2gi * rb2giecavpka + kact1gi * lrb2giecavpka -
    kreasgs * y[53] * y[54] - kreasgi * y[53] * y[162];
  //C
  //C   y[54]  Gs-alpha with GDP extracaveolae
  //C
  ydot[54] = khydgs * y[52] - kreasgs * y[53] * y[54];
  //C
  //C   y[161] Gi-alpha with GTP extracaveolae
  //C
  ydot[161] = kact2gi * rb2giecavpka + kact1gi * lrb2giecavpka -
    khydgi * y[161];
  //C
  //C   y[162]  Gi-alpha with GDP extracaveolae
  //C
  ydot[162] = khydgi * y[161] - kreasgi * y[53] * y[162];
  //C
  //C   y[55] beta1 tot concentration phosphorylated by PKA cytosol
  //C
  ydot[55] = kpkap * y[95] * rb1cytnptot - kpkam * y[55];
  //C
  //C   y[56] beta1 tot concentration phosphorylated by BARK
  //C
  ydot[56] = kbarkp * (lrb1cytnp + lrb1gscytnp) - kbarkm * y[56];
  //C
  //C   y[57] Gs-alpha with GTP
  //C
  ydot[57] = kact2gs * rb1gscytnp + kact1gs * lrb1gscytnp - khydgs * y[57];
  //C
  //C   y[58]  Gs Beta-Gamma
  //C
  ydot[58] = kact2gs * rb1gscytnp + kact1gs * lrb1gscytnp - kreasgs *
    y[58] * y[59];
  //C
  //C   y[59]  Gs-alpha with GDP
  //C
  ydot[59] = khydgs * y[57] - kreasgs * y[58] * y[59];
  //C
  //C   y[60] cAMP from AC56 in ceveolae
  //C
  ydot[60] = kcavac56 * ac56cav * atp / (kmatp + atp);
  //C
  //C   y[61] cAMP from AC47 in extracaveolae
  //C
  ydot[61] = kecavac47 * ac47ecav * atp / (kmatp + atp);
  //C
  //C   y[62] cAMP from AC56 in cytosol
  //C
  ydot[62] = kcytac56 * ac56cyt * atp / (kmatp + atp);
  //C
  //C   y[63] cAMP from AC47 in cytosol
  //C
  ydot[63] = kcytac47 * ac47cyt * atp / (kmatp + atp);
  //C
  //C   y[64] PDE3 caveolar phosphorylated
  //C
  ydot[64] = kfpdep * y[83] * (pde3cavtot - y[64]) - kbpdep * y[64];
  //C
  //C   y[65] PDE4 caveolar phosphorylated
  //C
  ydot[65] = kfpdep * y[83] * (pde4cavtot - y[65]) - kbpdep * y[65];
  //C
  //C   y[66] cAMP change due to PDE2 caveolar domain
  //C
  ydot[66] = (kpde2 * pde2cavtot * y[97]) / (kmpde2 + y[97]);
  //C
  //C   y[67] cAMP change due to PDE3 caveolar domain
  //C
  ydot[67] = (kpde3 * (pde3cavtot - y[64]) * y[97] + deltakpde34 *
    kpde3 * y[64] * y[97]) / (kmpde3 + y[97]);
  //C
  //C   y[68] cAMP change due to PDE4 caveolar domain
  //C
  ydot[68] = (kpde4 * (pde4cavtot - y[65]) * y[97] + deltakpde34 *
    kpde4 * y[65] * y[97]) / (kmpde4 + y[97]);
  //C
  //C   y[69] PDE3 extracaveolar phosphorylated
  //C
  //C      ydot[69] = kfpdep*y[89]*(PDE3ecavtot-y[69])-kbpdep*y[69]
  ydot[69] = 0.0f;
  //C
  //C   y[70] PDE4 extracaveolar phosphorylated
  //C
  ydot[70] = kfpdep * y[89] * (pde4ecavtot - y[70]) - kbpdep * y[70];
  //C
  //C   y[71] cAMP change due to PDE2 extracaveolar domain
  //C
  ydot[71] = (kpde2 * pde2ecavtot * y[98]) / (kmpde2 + y[98]);
  //C
  //C   y[72] cAMP change due to PDE3 extracaveolar domain
  //C
  ydot[72] = (kpde3 * (pde3ecavtot - y[69]) * y[98] + deltakpde34 *
    kpde3 * y[69] * y[98]) / (kmpde3 + y[98]);
  //C
  //C   y[73] cAMP change due to PDE4 extracaveolar domain
  //C
  ydot[73] = (kpde4 * (pde4ecavtot - y[70]) * y[98] + deltakpde34 *
    kpde4 * y[70] * y[98]) / (kmpde4 + y[98]);
  //C
  //C   y[74] PDE3 cytosol phosphorylated
  //C
  ydot[74] = kfpdep * y[95] * (pde3cyttot - y[74]) - kbpdep * y[74];
  //C
  //C   y[75] PDE4 cytosol phosphorylated
  //C
  ydot[75] = kfpdep * y[95] * (pde4cyttot - y[75]) - kbpdep * y[75];
  //C
  //C   y[76] cAMP change due to PDE2 cytosol domain
  //C
  ydot[76] = (kpde2 * pde2cyttot * y[99]) / (kmpde2 + y[99]);
  //C
  //C   y[77] cAMP change due to PDE3 cytosol domain
  //C
  ydot[77] = (kpde3 * (pde3cyttot - y[74]) * y[99] + deltakpde34 *
    kpde3 * y[74] * y[99]) / (kmpde3 + y[99]);
  //C
  //C   y[78] cAMP change due to PDE4 cytosol domain
  //C
  ydot[78] = (kpde4 * (pde4cyttot - y[75]) * y[99] + deltakpde34 *
    kpde4 * y[75] * y[99]) / (kmpde4 + y[99]);
  //C
  //C   y[79] cAMP binding to PKA caveolar
  //C
  ydot[79] = -kpkaiif1 * rccavf * y[97] + kpkaiib1 * y[80] -
    kpkaiif2 * y[80] * y[97] + kpkaiib2 * y[81];
  //C      write[51,44] time,-kpkaiif1*RCcavf*y[97],kpkaiib1*y[80],
  //C     *                  -kpkaiif2*y[80]*y[97],kpkaiib2*y[81]
  //C
  //C   y[80] RC with bound one cAMP in caveolae
  //C
  ydot[80] = kpkaiif1 * rccavf * y[97] - kpkaiib1 * y[80] -
    kpkaiif2 * y[80] * y[97] + kpkaiib2 * y[81];
  //C
  //C   y[81] RC with bound two cAMP in caveolae
  //C
  ydot[81] = kpkaiif2 * y[80] * y[97] - (kpkaiib2 + kpkaiif3) * y[
    81] + kpkaiib3 * y[82] * y[83];
  //C
  //C   y[82] R with bound two cAMP in caveolae
  //C
  ydot[82] = kpkaiif3 * y[81] - kpkaiib3 * y[82] * y[83];
  //C
  //C   y[83] concentration of catalytic subunit of PKA in caveolae
  //C
  ydot[83] = kpkaiif3 * y[81] - kpkaiib3 * y[82] * y[83] + kpkib * y[
    84] - kpkif * pkicavf * y[83];
  //C
  //C   y[84] concentration of PKI bound to catalytic subunit of PKA cav
  //C
  ydot[84] = -kpkib * y[84] + kpkif * pkicavf * y[83];
  //C
  //C   Extracaveolar domain
  //C
  //C   y[85] cAMP binding to PKA extracaveolar
  //C
  ydot[85] = -kpkaiif1 * rcecavf * y[98] + kpkaiib1 * y[86] -
    kpkaiif2 * y[86] * y[98] + kpkaiib2 * y[87];
  //C
  //C   y[86] RC with bound one cAMP in extracaveolae
  //C
  ydot[86] = kpkaiif1 * rcecavf * y[98] - kpkaiib1 * y[86] -
    kpkaiif2 * y[86] * y[98] + kpkaiib2 * y[87];
  //C
  //C   y[87] RC with bound two cAMP in extracaveolae
  //C
  ydot[87] = kpkaiif2 * y[86] * y[98] - (kpkaiib2 + kpkaiif3) * y[
    87] + kpkaiib3 * y[88] * y[89];
  //C
  //C   y[88] R with bound two cAMP in extracaveolae
  //C
  ydot[88] = kpkaiif3 * y[87] - kpkaiib3 * y[88] * y[89];
  //C
  //C   y[89] concentration of catalytic subunit of PKA in extracaveolae
  //C
  ydot[89] = kpkaiif3 * y[87] - kpkaiib3 * y[88] * y[89] + kpkib * y[
    90] - kpkif * pkiecavf * y[89];
  //C
  //C   y[90] concentration of PKI bound to catalytic subunit of PKA ecav
  //C
  ydot[90] = -kpkib * y[90] + kpkif * pkiecavf * y[89];
  //C
  //C   Cytosol domain
  //C
  //C   y[91] cAMP binding to PKA cytosol
  //C
  ydot[91] = -kpkaif1 * rccytf * y[99] + kpkaib1 * y[92] - kpkaif2 *
    y[92] * y[99] + kpkaib2 * y[93];
  //C
  //C   y[92] RC with bound one cAMP in cytosol
  //C
  ydot[92] = kpkaif1 * rccytf * y[99] - kpkaib1 * y[92] - kpkaif2 * y[
    92] * y[99] + kpkaib2 * y[93];
  //C
  //C   y[93] RC with bound two cAMP in cytosol
  //C
  ydot[93] = kpkaif2 * y[92] * y[99] - (kpkaib2 + kpkaif3) * y[93] +
    kpkaib3 * y[94] * y[95];
  //C
  //C   y[94] R with bound two cAMP in cytosol
  //C
  ydot[94] = kpkaif3 * y[93] - kpkaib3 * y[94] * y[95];
  //C
  //C   y[95] concentration of catalytic subunit of PKA in cytosol
  //C
  ydot[95] = kpkaif3 * y[93] - kpkaib3 * y[94] * y[95] + kpkib * y[
    96] - kpkif * pkicytf * y[95];
  //C
  //C   y[96] concentration of PKI bound to catalytic subunit of PKA cyt
  //C
  ydot[96] = -kpkib * y[96] + kpkif * pkicytf * y[95];
  //C
  //C   y[101] cAMP change due to PDE8 caveolar domain
  //C
  ydot[101] = (kpde8 * pde8cavtot * y[97]) / (kmpde8 + y[97]);
  //C
  //C   y[97] total change of cAMP in ceveolae
  //C
  ydot[97] = ydot[79] + ydot[60] - ydot[66] - ydot[67] - ydot[68] - ydot[
    101] - jcavecav * (y[97] - y[98]) / vcav - jcavcyt * (y[97] - y[
    99]) / vcav;
  double fluxcavecav = jcavecav * (y[97] - y[98]);
  double fluxcavcyt = jcavcyt * (y[97] - y[99]);
  //C
  //C   y[102] cAMP change due to PDE8 extracaveolar domain
  //C
  ydot[102] = (kpde8 * pde8ecavtot * y[98]) / (kmpde8 + y[98]);
  //C
  //C   y[98] total change of cAMP in extracaveolae
  //C
  ydot[98] = ydot[85] + ydot[61] - ydot[71] - ydot[72] - ydot[73] - ydot[
    102] - jcavecav * (y[98] - y[97]) / vecav - jecavcyt * (y[98] - y[
    99]) / vecav;
  double fluxecavcav = jcavecav * (y[98] - y[97]);
  double fluxecavcyt = jecavcyt * (y[98] - y[99]);
  //C
  //C   y[103] cAMP change due to PDE8 cytosolic domain
  //C
  ydot[103] = (kpde8 * pde8cyttot * y[99]) / (kmpde8 + y[99]);
  //C
  //C   y[99] total change of cAMP in cytosol
  //C
  ydot[99] = ydot[91] + ydot[62] + ydot[63] - ydot[76] - ydot[77] -
    ydot[78] - ydot[103] - jcavcyt * (y[99] - y[97]) / vcyt -
    jecavcyt * (y[99] - y[98]) / vcyt;
  double fluxcytcav = jcavcyt * (y[99] - y[97]);
  double fluxcytecav = jecavcyt * (y[99] - y[98]);
  //C
  //C    y[100] Inhibitor 1 cytosol phosphorylated
  //C
  ydot[100] = kpkainhib1 * y[95] * inhib1cytf / (kmpkainhib1 +
    inhib1cytf) - kpp2ainhib1pp2acyt * y[100] / (kmpp2ainhib1 + y[
    100]);
  //C
  //C    y[104] PLB fraction phosphorylated
  //C
  ydot[104] = kplbpka * y[95] * (1 - y[104]) / (kplbmpka + (1 - y[
    104])) - kplbpp1 * pp1cytf * y[104] / (kplbmpp1 + y[104]);
  //C
  //C    y[105] TnI fraction phosphorylated
  //C
  ydot[105] = ktnipka * y[95] * (1 - y[105]) / (ktnimpka + (1 - y[
    105])) - ktnipp2a * pp2acyt * y[105] / (ktnimpp2a + y[105]);
  //C
  //C    y[129] INaK fraction phosphorylated
  //C
  ydot[129] = kinakpka * y[83] * (1 - y[129]) / (kinakmpka + (1 - y[
    129])) - kinakpp * ppcav * y[129] / (kinakmpp + y[129]);
  //C
  //C    y[133] IK1 fraction nonphophorylated
  //C
  ydot[133] = kik1pp * pp1ecav * (1 - y[133]) / (kik1mpp + (1 - y[
    133])) - kik1pka * y[89] * y[133] / (kik1mpka + y[133]);
  //C
  //C    y[134] IKto fraction phophorylated
  //C
  ydot[134] = kiktopka * y[89] * (1 - y[134]) / (kiktompka + (1 - y[
    134])) - kiktopp * pp1ecav * y[134] / (kiktompp + y[134]);
  //C
  //C    y[163] ICaT C0 channel variable
  //C
  ydot[163] = k_catmv * y[164] - 4.0f * k_catpv * y[163] + k_catmi *
    y[169] / pow(h_cat,3) - k_catpi * y[163] * pow(f_cat,3);
  //C
  //C    y[164] ICaT C1 channel variable
  //C
  ydot[164] = 4.0f * k_catpv * y[163] - k_catmv * y[164] + 2.0f *
    k_catmv * y[165] - 3.0f * k_catpv * y[164] + k_catmi * y[170] /
    pow(h_cat,2) - k_catpi * y[164] * pow(f_cat,2);
  //C
  //C    y[165] ICaT C2 channel variable
  //C
  ydot[165] = 3.0f * k_catpv * y[164] - 2.0f * k_catmv * y[165] +
    3.0f * k_catmv * y[166] - 2.0f * k_catpv * y[165] + k_catmi * y[
    171] / h_cat - k_catpi * y[165] * f_cat;
  //C
  //C    y[166] ICaT C3 channel variable
  //C
  ydot[166] = 2.0f * k_catpv * y[165] - 3.0f * k_catmv * y[166] + 4.0f *
    k_catmv * y[167] - k_catpv * y[166] + k_catmi * y[172] - k_catpi * y[
    166];
  //C
  //C    y[167] ICaT C4 channel variable
  //C
  ydot[167] = k_catpv * y[166] - 4.0f * k_catmv * y[167] + k_catmo *
    y[168] - k_catpo * y[167] + k_catmi * y[173] - k_catpi * y[167];
  //C
  //C    y[168] ICaT O channel variable
  //C
  ydot[168] = k_catpo * y[167] - k_catmo * y[168] + k_catmi * y[
    174] - k_catpi * y[168];
  //C
  //C    y[169] ICaT I0 channel variable
  //C
  ydot[169] = k_catmv * y[170] * h_cat - 4.0f * k_catpv * y[169] /
    f_cat + k_catpi * y[163] * pow(f_cat,3) - k_catmi * y[169] /
    pow(h_cat,3);
  //C
  //C    y[170] ICaT I1 channel variable
  //C
  ydot[170] = 4.0f * k_catpv * y[169] / f_cat - k_catmv * y[170] *
    h_cat + 2.0f * k_catmv * y[171] * h_cat - 3.0f * k_catpv * y[170] /
    f_cat + k_catpi * y[164] * pow(f_cat,2) - k_catmi * y[170] /
    pow(h_cat,2);
  //C
  //C    y[171] ICaT I2 channel variable
  //C
  ydot[171] = 3.0f * k_catpv * y[170] / f_cat - 2.0f * k_catmv * y[
    171] * h_cat + 3.0f * k_catmv * y[172] * h_cat - 2.0f * k_catpv *
    y[171] / f_cat + k_catpi * y[165] * f_cat - k_catmi * y[171] /
    h_cat;
  //C
  //C    y[172] ICaT I3 channel variable
  //C
  ydot[172] = 2.0f * k_catpv * y[171] / f_cat - 3.0f * k_catmv * y[
    172] * h_cat + 4.0f * k_catmv * y[173] - k_catpv * y[172] +
    k_catpi * y[166] - k_catmi * y[172];
  //C
  //C     y[173] ICaT I4 channel variable
  //C
  ydot[173] = k_catpv * y[172] - 4.0f * k_catmv * y[173] + k_catmo *
    y[174] - k_catpo * y[173] + k_catpi * y[167] - k_catmi * y[173];
  //C
  //C    y[174] ICaT IO channel variable
  //C
  ydot[174] = k_catpo * y[173] - k_catmo * y[174] + k_catpi * y[
    168] - k_catmi * y[174];
  //C
  //C    y[175] ICaT fraction phophorylated caveolae
  //C
  ydot[175] = kicatpka * y[83] * (1 - y[175]) / (kicatmpka + (1 - y[
    175])) - kicatpp * ppcav * y[175] / (kicatmpp + y[175]);
  //C
  //C    y[176] ICaT fraction phophorylated cytosol
  //C
  ydot[176] = kicatpka * y[95] * (1 - y[176]) / (kicatmpka + (1 - y[
    176])) - kicatpp * (pp2acyt + pp1cytf) * y[176] / (kicatmpp + y[
    176]);
  //C
  //C  Small-conductance calcium activated K+ channel
  //C
  ydot[177] = k_cakb1 * y[178] - k_cakf1 * y[177];
  ydot[178] = k_cakf1 * y[177] - k_cakb1 * y[178] + k_cakb2 * y[
    179] - k_cakf2 * y[178];
  ydot[179] = k_cakf2 * y[178] - k_cakb2 * y[179] + k_cakb3 * y[
    180] - k_cakf3 * y[179] + k_cakbo1 * y[181] - k_cakfo1 * y[179];
  ydot[180] = k_cakf3 * y[179] - k_cakb3 * y[180] + k_cakbo2 * y[
    182] - k_cakfo2 * y[180];
  ydot[181] = k_cakfo1 * y[179] - k_cakbo1 * y[181];
  ydot[182] = k_cakfo2 * y[180] - k_cakbo2 * y[182];
  //C
  //C    y[183] ICaK fraction phophorylated caveolae
  //C
  ydot[183] = kicakpka * y[83] * (1 - y[183]) / (kicakmpka + (1 - y[
    183])) - kicakpp * ppcav * y[183] / (kicakmpp + y[183]);
  //C


}

