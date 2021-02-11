#ifndef _MODEL_H_
#define _MODEL_H_

#define S_SIZE 25
#define C_SIZE 91
#define A_SIZE 20

void initialize_states_default(double *STATES);

void initialize_constants_default(double *CONSTANTS);

void compute_rates_algebraic(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES);

#endif // _MODEL_H_
