#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../liblsoda/src/common.h"
#include "../../liblsoda/src/lsoda.h"
#include "./model.h"

#define USE_LSODA 1

int rhs(double t, double *y, double *ydot, void *data) {

    double *C = (double *)data, *A = ((double *)data) + C_SIZE;
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


int run(double *S, double *C, int n_beats, double t_sampling, double tol,
        double *output/*, double *output_A, double *output_t*/) {

  double data[C_SIZE + A_SIZE];
  for (int i = 0; i < C_SIZE; ++i) {
    data[i] = C[i];
  }

  double  stim_period = C[87];
  int     n_samples   = stim_period / t_sampling;
  //  int     n_samples_per_stim   = 0; // ceil((C[84] + C[83]) / t_sampling); // from start till the end of the stimulation

  double  atol[S_SIZE], rtol[S_SIZE];
  int     neq = S_SIZE;

  struct lsoda_opt_t opt = {0};
  opt.ixpr    = 0;
  opt.rtol    = rtol;
  opt.atol    = atol;
  opt.itask   = 1;

  for (int i = 0; i < S_SIZE; ++i) {
    rtol[i] = tol;
    atol[i] = tol;
    if ((i == 2) || (i == 21)) {
      atol[i] *= 1e-4;  // Cai, O
    }
    if ((i == 15) || (i == 20) || (i == 22)) {
      atol[i] *= 1e-2;  // Xf, R, I
    }
  }

  double  t               = 0;
  double  t_out           = 0;
  int     i_out_global    = 1;
  int     ctx_state       = 0;

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

#if USE_LSODA == 0

    for (; i_out_local <= n_samples; i_out_local++, i_out_global++) {
      t_out = i_out_global * t_sampling;
      euler(&t, S, data, 1e-3, t_out);
      memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));
    }

#else // USE_LSODA == 0

    struct lsoda_context_t ctx = {
      .function = rhs,
      .neq = neq,
      .data = data,
      .state = 1,
    };
    lsoda_prepare(&ctx, &opt);

    for (; i_out_local <= n_samples; i_out_local++, i_out_global++) {

      t_out = i_out_global * t_sampling;

      // double stim_shift = data[83];
      // data[82] = (i_out_local > data[83]) && (i_out_local <= data[83] + 5); //  const double i_stim_PulseDuration = 5;   // milisecond (in stim_mode)

      lsoda(&ctx, S, &t, t_out);

      memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));
      //memcpy(output_A + i_out_global * A_SIZE, A, A_SIZE * sizeof(double));
      //memcpy(output_t + i_out_global, &t, 1 * sizeof(double));

      if (ctx.state != 2) {
        return ctx.state;
      }
      // printf("L: %f %f\n", t, S[23]);

    }

    // printf("beat %d completed\n", i_beat);

    ctx_state = ctx.state;
    lsoda_free(&ctx);

#endif // USE_LSODA

}

    return ctx_state;
}
