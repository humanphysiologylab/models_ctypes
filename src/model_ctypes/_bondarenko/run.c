#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../liblsoda/src/common.h"
#include "../../liblsoda/src/lsoda.h"
#include "./bondarenko.h"

int rhs(double t, double *y, double *ydot, void *data) {

    double *C = (double *)data, *A = ((double *)data) + C_SIZE;

    // TODO: fix fun()
    for (int i = 0; i < S_SIZE; ++i) {
	ydot[i] = 0;
    };

    fun(t, y, ydot, A, C);
    return 0;
}


int euler(double *t, double *y, void *data,
          double dt, double t_out) {

    double ydot[S_SIZE];

    while (*t < t_out) {

        // printf("Euler step: %f %f; ", *t, y[23]);

        rhs(*t, y, ydot, data);

        if (*t + dt <= t_out) {
            *t += dt;
        } else {
            dt = t_out - *t;
            *t = t_out;
        }

        for (int i = 0; i < S_SIZE; ++i) {
            y[i] += dt * ydot[i];
        }
    }

    return 0;
}


int run(double *S, double *C, int n_beats, double t_sampling, double tol, double *output,
        /* can be None -> */ double *output_A, double *output_t, double *stim_protocol) {


    double data[C_SIZE + A_SIZE];
    double *A = data + C_SIZE;

    for (int i = 0; i < C_SIZE; ++i) {
        data[i] = C[i];
    }

    double  stim_period = C[0];
    int     n_samples   = stim_period / t_sampling;

    double  atol[S_SIZE], rtol[S_SIZE];
    int     neq = S_SIZE;

    struct lsoda_opt_t opt = {0};
    opt.ixpr    = 0;
    opt.rtol    = rtol;
    opt.atol    = atol;
    opt.itask   = 1;

    double atol_mult[] = {/*s_0*/ 1e-06,      /*V*/ 1.0,          /*Cai*/ 0.01,       /*Cansr*/ 100.0,    /*cLTRPNca*/ 1.0,
/*cHTRPNca*/ 10.0,  /*cCajsr*/ 100.0,   /*cCass*/ 0.01,     /*O_CaL_cav*/ 1e-06, /*C1_CaL_cav*/ 0.1,
/*C2_CaL_cav*/ 1e-06, /*C3_CaL_cav*/ 1e-06, /*C4_CaL_cav*/ 1e-06, /*I1_CaL_cav*/ 1e-06, /*I2_CaL_cav*/ 1e-06,
/*I3_CaL_cav*/ 1e-06, /*C1_RyR*/ 1e-06,   /*C2_RyR*/ 1e-05,   /*O1_RyR*/ 1e-06,   /*O2_RyR*/ 1e-06,
/*CNa3_INa*/ 0.1,   /*CNa2_INa*/ 1e-06, /*CNa1_INa*/ 1e-06, /*Na_intra*/ 1000.0, /*K_intra*/ 10000.0,
/*ato_Ikto*/ 0.001, /*ito_Ikto*/ 0.1,   /*nks_Iks*/ 0.0001, /*aur_Ikur1*/ 0.0001, /*iur_Ikur1*/ 0.1,
/*HERG_v1*/ 0.1,    /*HERG_v2*/ 0.0001, /*HERG_v3*/ 0.0001, /*HERG_v4*/ 1e-05,  /*HERG_v5*/ 1e-05,
/*aur_Ikur2*/ 0.0001, /*iur_Ikur2*/ 0.1,  /*ONa_INa*/ 1e-06,  /*IFNa_ONa*/ 1e-05, /*I1Na_INa*/ 1e-06,
/*I2Na_INa*/ 1e-06, /*ICNa2_INa*/ 0.0001, /*ICNa3_INa*/ 1e-06, /*aur_Ikur3*/ 0.0001, /*Pryr*/ 1e-06,
/*beta1_tot_p_PKA_cav*/ 0.0001, /*beta1_tot_p_BARK_cav*/ 1e-06, /*Gs_alpa_GTP_cav*/ 0.001, /*G_Beta_Gamma_cav*/ 0.001, /*Gs_alpha_GDP_cav*/ 0.0001,
/*s_50*/ 0.01,      /*s_51*/ 1e-06,     /*s_52*/ 0.001,     /*s_53*/ 0.001,     /*s_54*/ 0.0001,
/*s_55*/ 0.0001,    /*s_56*/ 1e-06,     /*s_57*/ 0.0001,    /*s_58*/ 0.0001,    /*s_59*/ 0.0001,
/*s_60*/ 10000.0,   /*s_61*/ 10000.0,   /*s_62*/ 1000.0,    /*s_63*/ 100.0,     /*s_64*/ 0.001,
/*s_65*/ 0.001,     /*s_66*/ 1000.0,    /*s_67*/ 10000.0,   /*s_68*/ 1000.0,    /*s_69*/ 1e-06,
/*s_70*/ 0.001,     /*s_71*/ 1000.0,    /*s_72*/ 0.01,      /*s_73*/ 10000.0,   /*s_74*/ 0.0001,
/*s_75*/ 0.001,     /*s_76*/ 1000.0,    /*s_77*/ 1000.0,    /*s_78*/ 1000.0,    /*s_79*/ 1.0,
/*s_80*/ 0.1,       /*s_81*/ 0.01,      /*s_82*/ 0.1,       /*s_83*/ 0.01,      /*s_84*/ 0.1,
/*s_85*/ 1.0,       /*s_86*/ 0.1,       /*s_87*/ 0.01,      /*s_88*/ 0.1,       /*s_89*/ 0.01,
/*s_90*/ 0.1,       /*s_91*/ 1.0,       /*s_92*/ 0.01,      /*s_93*/ 0.001,     /*s_94*/ 0.01,
/*s_95*/ 0.01,      /*s_96*/ 0.01,      /*s_97*/ 0.1,       /*s_98*/ 0.1,       /*s_99*/ 0.1,
/*s_100*/ 0.001,    /*s_101*/ 0.1,      /*s_102*/ 0.1,      /*s_103*/ 0.01,     /*s_104*/ 0.1,
/*s_105*/ 0.1,      /*CP_CaL_cav*/ 1e-06, /*OP_CaL_cav*/ 1e-06, /*C1P_CaL_cav*/ 1e-06, /*C2P_CaL_cav*/ 1e-06,
/*C3P_CaL_cav*/ 1e-06, /*C4P_CaL_cav*/ 1e-06, /*I1P_CaL_cav*/ 1e-06, /*I2P_CaL_cav*/ 1e-06, /*I3P_CaL_cav*/ 1e-06,
/*CPP_CaL_cav*/ 1e-06, /*C1P_RyR*/ 1e-06,  /*C2P_RyR*/ 1e-06,  /*O1P_RyR*/ 1e-06,  /*O2P_RyR*/ 1e-06,
/*CNa3p_INa*/ 1e-06, /*CNa2p_INa*/ 1e-06, /*CNa1p_INa*/ 1e-06, /*ONap_INa*/ 1e-06, /*IFNap_INa*/ 1e-06,
/*I1Nap_INa*/ 1e-06, /*I2Nap_INa*/ 1e-06, /*ICNa2p_INa*/ 0.0001, /*ICNa3p_INa*/ 1e-06, /*s_129*/ 0.1,
/*s_130*/ 0.0001,   /*s_131*/ 0.1,      /*s_132*/ 0.1,      /*s_133*/ 0.1,      /*s_134*/ 0.01,
/*s_135*/ 0.0001,   /*s_136*/ 0.1,      /*O_CaL_ecav*/ 1e-06, /*C1_CaL_ecav*/ 0.1, /*C2_CaL_ecav*/ 1e-06,
/*C3_CaL_ecav*/ 1e-06, /*C4_CaL_ecav*/ 1e-06, /*CP_CaL_ecav*/ 1e-06, /*I1_CaL_ecav*/ 1e-06, /*I2_CaL_ecav*/ 1e-06,
/*I3_CaL_ecav*/ 1e-06, /*OP_CaL_ecav*/ 1e-06, /*C1P_CaL_ecav*/ 1e-06, /*C2P_CaL_ecav*/ 1e-06, /*C3P_CaL_ecav*/ 1e-06,
/*C4P_CaL_ecav*/ 1e-06, /*CPP_CaL_ecav*/ 1e-06, /*I1P_CaL_ecav*/ 1e-06, /*I2P_CaL_ecav*/ 1e-06, /*I3P_CaL_ecav*/ 1e-06,
/*s_155*/ 0.01,     /*s_156*/ 1e-06,    /*Gi_alpha_GTP_cav*/ 0.001, /*Gi_alpha_GDP_cav*/ 0.0001, /*beta2_tot_p_PKA_ecav*/ 0.0001,
/*beta2_tot_p_BARK_ecav*/ 1e-06, /*s_161*/ 1e-06,    /*s_162*/ 1e-06,    /*s_163*/ 1e-06,    /*s_164*/ 1e-06,
/*s_165*/ 1e-06,    /*s_166*/ 1e-06,    /*s_167*/ 1e-06,    /*s_168*/ 1e-06,    /*s_169*/ 1e-06,
/*s_170*/ 1e-06,    /*s_171*/ 1e-06,    /*s_172*/ 1e-06,    /*s_173*/ 1e-06,    /*s_174*/ 1e-05,
/*s_175*/ 0.001,    /*s_176*/ 0.001,    /*s_177*/ 0.1,      /*s_178*/ 0.01,     /*s_179*/ 0.001,
/*s_180*/ 1e-05,    /*s_181*/ 0.0001,   /*s_182*/ 0.001,    /*s_183*/ 0.001,    /*aur_Ikur4*/ 0.0001,
/*iur_Ikur4*/ 0.1};

    for (int i = 0; i < S_SIZE; ++i) {
        rtol[i] = tol;
        atol[i] = atol_mult[i];
    }

    double  t               = 0;
    double  t_out           = 0;
    int     i_out_global    = 1;
    int     ctx_state       = 0;

    memcpy(output, S, S_SIZE * sizeof(double));

    if (output_A != 0) {
        double R[S_SIZE];
        rhs(0, S, R, data);
        memcpy(output_A, A, A_SIZE * sizeof(double));
    }

    if (output_t != 0) {
        memcpy(output_t, &t, 1 * sizeof(double));
    }

    for (int i_beat = 0; i_beat < n_beats; ++i_beat) {

        int i_out_local = 1;

    //        for (; i_out_local <= n_samples_per_stim; i_out_local++, i_out_global++) {
    //            t_out = i_out_global * t_sampling;
    //            euler(&t, S, data, 1e-5, t_out);
    //            memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));
    //            printf("E: %f %f\n", t, S[23]);
    //        }

        struct lsoda_context_t ctx = {
            .function = rhs,
            .neq = neq,
            .data = data,
            .state = 1,
        };
        lsoda_prepare(&ctx, &opt);

        for (; i_out_local <= n_samples; i_out_local++, i_out_global++) {

            t_out = i_out_global * t_sampling;

            // TODO: move max amplitude to CONSTANTS
            if (stim_protocol != 0) {
                A[5] = stim_protocol[i_out_local];
            } else { // if stim_protocol is None in the python code
                A[5] = (i_out_local == 1) * 80;
            }

            lsoda(&ctx, S, &t, t_out);

            memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));

            if (output_A != 0) {
                memcpy(output_A + i_out_global * A_SIZE, A, A_SIZE * sizeof(double));
            }
            if (output_t != 0) {
                memcpy(output_t + i_out_global, &t, 1 * sizeof(double));
            }

            if (ctx.state != 2) {
                return ctx.state;
            }
            // printf("L: %f %f\n", t, S[23]);

        }

        // printf("beat %d completed\n", i_beat);

        ctx_state = ctx.state;
        lsoda_free(&ctx);

    }

    return ctx_state;
}
