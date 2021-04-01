#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../liblsoda/src/common.h"
#include "../../liblsoda/src/lsoda.h"
#include "./model.h"


int rhs(double t, double *y, double *ydot, void *data) {

    double *C = (double *)data, *A = ((double *)data) + C_SIZE;
    computeRates(t, C, ydot, y, A);
    return 0;
}


int run(double *S, double *C, int n_beats, double t_sampling, double tol,
        double *output/*, double *output_A, double *output_t*/) {

    double data[C_SIZE + A_SIZE];
    for (int i = 0; i < C_SIZE; ++i) {
        data[i] = C[i];
    }

    double  stim_period = C[101];
    int     n_samples   = stim_period / t_sampling;

    double  atol[S_SIZE], rtol[S_SIZE];
	int     neq = S_SIZE;

	struct lsoda_opt_t opt = {0};
    opt.ixpr    = 0;
    opt.rtol    = rtol;
    opt.atol    = atol;
    opt.itask   = 1;

    double atol_mult[] = {/*CaM*/ 0.0001,     /*Myoc*/ 0.001,     /*Myom*/ 0.001,     /*SLH_jn*/ 0.001,   /*SLH_sl*/ 0.001,
                          /*SLL_jn*/ 0.001,   /*SLL_sl*/ 0.001,   /*SRB*/ 0.001,      /*TnCHc*/ 0.001,    /*TnCHm*/ 0.001,
                          /*TnCL*/ 0.001,     /*Ca_i*/ 1e-05,     /*Ca_jn*/ 0.0001,   /*Ca_sl*/ 1e-05,    /*NaB_jn*/ 0.001,
                          /*NaB_sl*/ 0.001,   /*Na_jn*/ 0.001,    /*Na_sl*/ 0.001,    /*Ca_sr*/ 0.001,    /*Csqn*/ 0.001,
                          /*V*/ 0.001,        /*d*/ 1e-06,        /*f*/ 0.001,        /*fCaB_jn*/ 0.001,  /*fCaB_sl*/ 0.001,
                          /*K_i*/ 0.001,      /*xr*/ 0.0001,      /*xs*/ 0.001,       /*Na_i*/ 0.001,     /*ikur_r*/ 0.0001,
                          /*s*/ 0.001,        /*h*/ 1e-06,        /*j*/ 1e-06,        /*m*/ 0.001,        /*hl*/ 0.001,
                          /*ml*/ 0.001,       /*x*/ 0.0001,       /*y*/ 0.001,        /*i*/ 1e-06,        /*o*/ 1e-06,
                          /*ryr_r*/ 0.001,    /*fluo_i*/1e-6,     /*fluo_jn*/1e-6,    /*fluo_ss*/1e-6,    /*Cai_mean*/ 1e-6,
                          /*fluo_mean*/ 1e-6};

    for (int i = 0; i < S_SIZE; ++i) {
        rtol[i] = tol;
        atol[i] = atol_mult[i];
    }

    double  t               = 0;
    double  t_out           = 0;
    int     i_out_global    = 1;
    int     ctx_state       = 0;

    _calc_means(S, data);

    memcpy(output, S, S_SIZE * sizeof(double));

    /*
    double *A = data + C_SIZE;
    double R[S_SIZE];
    rhs(0, S, R, data);
    memcpy(output_A, A, A_SIZE * sizeof(double));

    memcpy(output_t, &t, 1 * sizeof(double));
    */

    for (int i_beat = 0; i_beat < n_beats; ++i_beat) {

        int i_out_local = 1;

        struct lsoda_context_t ctx = {
            .function = rhs,
            .neq = neq,
            .data = data,
            .state = 1,
        };
        lsoda_prepare(&ctx, &opt);

        for (; i_out_local <= n_samples; i_out_local++, i_out_global++) {

            t_out = i_out_global * t_sampling;

            data[C_SIZE + 118] = (i_out_local == 1); // ALGEBRAIC[118] = pace;
            lsoda(&ctx, S, &t, t_out);

            _calc_means(S, data);

            memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));
            //memcpy(output_A + i_out_global * A_SIZE, A, A_SIZE * sizeof(double));
            //memcpy(output_t + i_out_global, &t, 1 * sizeof(double));

            if (ctx.state != 2) {
                return ctx.state;
            }

	    }

        ctx_state = ctx.state;
        lsoda_free(&ctx);

	}

	return ctx_state;
}


int run_algebraic(double *S, double *C, int n_beats, double t_sampling, double tol,
                  double *output, double *output_A, double *output_t) {

    double data[C_SIZE + A_SIZE];
    double *A = data + C_SIZE;
    for (int i = 0; i < C_SIZE; ++i) {
        data[i] = C[i];
    }

    double  stim_period = C[101];
    int     n_samples   = stim_period / t_sampling;

    double  atol[S_SIZE], rtol[S_SIZE];
	int     neq = S_SIZE;

	struct lsoda_opt_t opt = {0};
    opt.ixpr    = 0;
    opt.rtol    = rtol;
    opt.atol    = atol;
    opt.itask   = 1;

    double atol_mult[] = {/*CaM*/ 0.0001,     /*Myoc*/ 0.001,     /*Myom*/ 0.001,     /*SLH_jn*/ 0.001,   /*SLH_sl*/ 0.001,
                          /*SLL_jn*/ 0.001,   /*SLL_sl*/ 0.001,   /*SRB*/ 0.001,      /*TnCHc*/ 0.001,    /*TnCHm*/ 0.001,
                          /*TnCL*/ 0.001,     /*Ca_i*/ 1e-05,     /*Ca_jn*/ 0.0001,   /*Ca_sl*/ 1e-05,    /*NaB_jn*/ 0.001,
                          /*NaB_sl*/ 0.001,   /*Na_jn*/ 0.001,    /*Na_sl*/ 0.001,    /*Ca_sr*/ 0.001,    /*Csqn*/ 0.001,
                          /*V*/ 0.001,        /*d*/ 1e-06,        /*f*/ 0.001,        /*fCaB_jn*/ 0.001,  /*fCaB_sl*/ 0.001,
                          /*K_i*/ 0.001,      /*xr*/ 0.0001,      /*xs*/ 0.001,       /*Na_i*/ 0.001,     /*ikur_r*/ 0.0001,
                          /*s*/ 0.001,        /*h*/ 1e-06,        /*j*/ 1e-06,        /*m*/ 0.001,        /*hl*/ 0.001,
                          /*ml*/ 0.001,       /*x*/ 0.0001,       /*y*/ 0.001,        /*i*/ 1e-06,        /*o*/ 1e-06,
                          /*ryr_r*/ 0.001,    /*fluo_i*/1e-6,     /*fluo_jn*/1e-6,    /*fluo_ss*/1e-6,    /*Cai_mean*/ 1e-6,
                          /*fluo_mean*/ 1e-6};

    for (int i = 0; i < S_SIZE; ++i) {
        rtol[i] = tol;
        atol[i] = atol_mult[i];
    }

    double  t               = 0;
    double  t_out           = 0;
    int     i_out_global    = 1;
    int     ctx_state       = 0;

    _calc_means(S, data);
    computeVariables(t, C, S, A);

    memcpy(output, S, S_SIZE * sizeof(double));
    memcpy(output_A, A, A_SIZE * sizeof(double));
    memcpy(output_t, &t, 1 * sizeof(double));

    for (int i_beat = 0; i_beat < n_beats; ++i_beat) {

        int i_out_local = 1;

        struct lsoda_context_t ctx = {
            .function = rhs,
            .neq = neq,
            .data = data,
            .state = 1,
        };
        lsoda_prepare(&ctx, &opt);

        for (; i_out_local <= n_samples; i_out_local++, i_out_global++) {

            t_out = i_out_global * t_sampling;

            data[C_SIZE + 118] = (i_out_local == 1); // ALGEBRAIC[118] = pace;
            lsoda(&ctx, S, &t, t_out);

            _calc_means(S, data);
            computeVariables(t, C, S, A);

            memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));
            memcpy(output_A + i_out_global * A_SIZE, A, A_SIZE * sizeof(double));
            memcpy(output_t + i_out_global, &t, 1 * sizeof(double));

            if (ctx.state != 2) {
                return ctx.state;
            }

	    }

        ctx_state = ctx.state;
        lsoda_free(&ctx);

	}

	return ctx_state;
}
