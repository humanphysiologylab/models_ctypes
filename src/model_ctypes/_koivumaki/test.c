#include <stdio.h>
#include "koivumaki.h"
#include "run.h"

// int run(double *S, double *C, int n_beats, double t_sampling, double tol, double *output,
//            /* can be None -> */ double *output_A, double *output_t, double *stim_protocol_Ist, double *stim_protocol_t)

int main() {

    double S[S_SIZE];
    double C[C_SIZE];

    FILE *f;
    f = fopen("legend_states.bin", "rb");
    fread(S, sizeof *S, S_SIZE, f);
    fclose(f);
    for (int i = 0; i < S_SIZE; ++i){
        printf("S[%d] = %e\n", i, S[i]);
    }

    f = fopen("legend_constants.bin", "rb");
    fread(C, sizeof *C, C_SIZE, f);
    fclose(f);
    for (int i = 0; i < C_SIZE; ++i){
        printf("C[%d] = %e\n", i, C[i]);
    }

    int OUTPUT_SIZE = 10001 * S_SIZE;
    double output[OUTPUT_SIZE];  // hahaha
    // f = fopen("output_zeros.bin", "rb");
    // fread(output, sizeof *output, OUTPUT_SIZE, f);
    // fclose(f);

    int n_beats = 10;
    double t_sampling = 1e-3;
    double tol = 1e-4;

    printf("Prepared to run\n");

    int status = run(S, C, n_beats, t_sampling, tol, output, 0, 0, 0, 0);

    f = fopen("output.bin", "wb");
    fwrite(output, sizeof(output), 1, f);
    fclose(f);

    printf("status %d\nDONE\n", status);

    return 0;
}
