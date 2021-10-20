#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../liblsoda/src/common.h"
#include "../../liblsoda/src/lsoda.h"
#include "./koivumaki.h"

int rhs(double t, double *y, double *ydot, void *data) {

    double *C = (double *)data, *A = ((double *)data) + C_SIZE;
    for (int i = 0; i < S_SIZE; ++i) {
        ydot[i] = 0.;  // for safety
    }
    compute_rates_algebraic(t, y, C, A, ydot);
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
           /* can be None -> */ double *output_A, double *output_t, double *stim_protocol_Ist, double *stim_protocol_t) {

    double CaSR = S[48];
    _calc_means(S, C);
    for (int i = 0; i < 4; ++i) {
      S[i] *= CaSR / S[48]; // CaSRi *= CaSR_desired / CaSR_real;
    }

    double data[C_SIZE + A_SIZE];
    double *A = data + C_SIZE;

    for (int i = 0; i < C_SIZE; ++i) {
        if (i < C_SIZE) {
            data[i] = C[i];
        } else {
            data[i] = 0.;
        }
    }

    double  stim_period = C[85];

    double  atol[S_SIZE], rtol[S_SIZE];
    int     neq = S_SIZE;

    struct lsoda_opt_t opt = {0};
    opt.ixpr    = 0;
    opt.rtol    = rtol;
    opt.atol    = atol;
    opt.itask   = 1;

    double atol_mult[] = {/*CaSR1*/ 0.001,    /*CaSR2*/ 0.001,    /*CaSR3*/ 0.001,    /*CaSR4*/ 0.001,    /*Cai1*/ 1e-05,
                          /*Cai2*/ 1e-05,     /*Cai3*/ 1e-05,     /*Cai4*/ 1e-05,     /*Cass*/ 1e-05,     /*d*/ 1e-06,
                          /*f1*/ 0.001,       /*f2*/ 0.001,       /*fca*/ 0.001,      /*y*/ 0.0001,       /*pa*/ 1e-06,
                          /*n*/ 0.001,        /*ikur_r*/ 1e-05,   /*ikur_s*/ 0.001,   /*h1*/ 1e-06,       /*h2*/ 1e-05,
                          /*m*/ 0.0001,       /*it_r*/ 0.0001,    /*it_s*/ 0.001,     /*V*/ 0.0001,       /*Ki*/ 0.001,
                          /*ryr_a1*/ 0.001,   /*ryr_a2*/ 0.001,   /*ryr_a3*/ 0.001,   /*ryr_ass*/ 0.001,  /*c1*/ 0.001,
                          /*c2*/ 0.001,       /*c3*/ 0.001,       /*css*/ 1e-06,      /*o1*/ 1e-06,       /*o2*/ 1e-06,
                          /*o3*/ 1e-06,       /*oss*/ 1e-06,      /*serca_a1*/ 1e-6, /*serca_a2*/ 1e-6, /*serca_a3*/ 1e-6,
                          /*serca_ass*/ 1e-6, /*Nai*/ 0.001,      /*Nass*/ 0.001,     /*fluo_1*/ 1e-05,   /*fluo_2*/ 1e-05,
                          /*fluo_3*/ 1e-05,   /*fluo_4*/ 1e-05,   /*fluo_ss*/ 1e-05,  /*CaSR*/ 0.001,     /*Cai*/ 1e-05,
                          /*fluo*/ 1e-05};

    for (int i = 0; i < S_SIZE; ++i) {
        rtol[i] = 1e-6; // tol;
        atol[i] = 1e-9; // atol_mult[i];
    }

    double  t               = 0;
    double  t_out           = 0;
    int     i_out_global    = 1;
    int     ctx_state       = 0;


    double dt_euler = tol; // 5e-6;

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

        int i_out_local = 1, i_protocol = 1;

        double t_start_local = stim_period * i_beat;
        double t_end_local = stim_period * (i_beat + 1);
        double t_end_local_euler = stim_period * i_beat + 0.025;  // s = 25 ms


        /* ----------------------*/
        /* - - - - EULER - - - - */

        while (t < t_end_local_euler) {

            t_out = i_out_global * t_sampling;

            double t_protocol = (stim_protocol_t != 0) ? stim_protocol_t[i_protocol] + t_start_local : t_out;

            if (stim_protocol_Ist != 0) {
                A[56] = -stim_protocol_Ist[i_protocol] * C[82] /* STIM_LEVEL */;
            } else {
                C[82] = (i_out_local * t_sampling <= C[83] /* STIM_DURATION */) ? 1 : 0;
                A[56] = C[82] * C[34]; // STIM LEVEL * amplitude
            }

            double t_nearest = (t_out < t_protocol) ? t_out : t_protocol;

            euler(&t, S, data, dt_euler, t_nearest);
            _calc_means(S, data);

            if (t == t_out) {

                memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));

                if (output_A != 0) {
                    memcpy(output_A + i_out_global * A_SIZE, A, A_SIZE * sizeof(double));
                }

                if (output_t != 0) {
                    memcpy(output_t + i_out_global, &t, 1 * sizeof(double));
                }

                i_out_local++; i_out_global++;
            }

            if (t == t_protocol) {
                i_protocol++;
            }

        }

        /* ----------------------*/
        /* - - - - LSODA - - - - */

        struct lsoda_context_t ctx = {
            .function = rhs,
            .neq = neq,
            .data = data,
            .state = 1,
        };
        lsoda_prepare(&ctx, &opt);

        while (t < t_end_local) {

            // t_out = i_out_local * t_sampling + t_start_local;
            t_out = i_out_global * t_sampling;

            double t_protocol = (stim_protocol_t != 0) ? stim_protocol_t[i_protocol] + t_start_local : t_out;

            if (stim_protocol_Ist != 0) {
                A[56] = -stim_protocol_Ist[i_protocol] * C[82] /* STIM_LEVEL */;
            } else {
                C[82] = (i_out_local * t_sampling <= C[83] /* STIM_DURATION */) ? 1 : 0;
                A[56] = C[82] * C[34]; // STIM LEVEL * amplitude
            }

            double t_nearest = (t_out < t_protocol) ? t_out : t_protocol;

            lsoda(&ctx, S, &t, t_nearest);
            _calc_means(S, data);

            if (t == t_out) {

                memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));

                if (output_A != 0) {
                    memcpy(output_A + i_out_global * A_SIZE, A, A_SIZE * sizeof(double));
                }

                if (output_t != 0) {
                    memcpy(output_t + i_out_global, &t, 1 * sizeof(double));
                }

                i_out_local++; i_out_global++;
            }

            if (t == t_protocol) {
                i_protocol++;
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
