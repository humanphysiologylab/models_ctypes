#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../liblsoda/src/common.h"
#include "../../liblsoda/src/lsoda.h"
#include "./gattoni.h"

int rhs(double t, double *y, double *ydot, void *data) {

    double *C = (double *)data, *A = ((double *)data) + C_SIZE;
    computeRates(t, C, ydot, y, A);
    return 0;
}


int run(double *S, double *C,
        int n_beats, double t_sampling, double tol,
        double *output,
        /* can be None -> */ double *output_A) {

    double data[C_SIZE + A_SIZE];
    double *A = data + C_SIZE;

    for (int i = 0; i < C_SIZE; ++i) {
        data[i] = C[i];
    }

    double  stim_period = C[4];
    int     n_samples   = stim_period / t_sampling;

    double  atol[S_SIZE], rtol[S_SIZE];
    int     neq = S_SIZE;

    struct lsoda_opt_t opt = {0};
    opt.ixpr    = 0;
    opt.rtol    = rtol;
    opt.atol    = atol;
    opt.itask   = 1;

    double atol_mult[] = {/*V*/ 0.001,        /*Na_i*/ 0.001,     /*m*/ 0.0001,       /*h*/ 1e-06,        /*j*/ 0.001,
                          /*K_i*/ 0.001,      /*r*/ 0.0001,       /*s*/ 0.001,        /*s_slow*/ 0.001,   /*r_ss*/ 0.0001,
                          /*s_ss*/ 0.001,     /*y*/ 0.0001,       /*Ca_i*/ 1e-05,     /*Ca_SR*/ 0.001,    /*z_1*/ 0.001,
                          /*z_2*/ 0.001,      /*z_3*/ 0.001,      /*TRPN*/ 0.001};

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
        computeVariables(t, data, S, A);
        memcpy(output_A, A, A_SIZE * sizeof(double));
    }

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

            lsoda(&ctx, S, &t, t_out);

            memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));

            if (output_A != 0) {
                computeVariables(t, data, S, A);
                memcpy(output_A + i_out_global * A_SIZE, A, A_SIZE * sizeof(double));
            }

            if (ctx.state != 2) {
                return ctx.state;
            }
        }

        ctx_state = ctx.state;
        lsoda_free(&ctx);

    }

    return ctx_state;
}
