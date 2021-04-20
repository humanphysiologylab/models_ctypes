#ifndef _KOIVUMAKI_H_
#define _KOIVUMAKI_H_

#define S_SIZE 51
#define C_SIZE 93
#define A_SIZE 107

void initialize_states_default(double *STATES, const double *CONSTANTS);

void initialize_constants_default(double *CONSTANTS);

void compute_rates_algebraic(const double time, double *STATES, double *CONSTANTS, double *ALGEBRAIC, double *RATES);

void _calc_means(double *STATES, const double *CONSTANTS);

#endif // _KOIVUMAKI_H_
