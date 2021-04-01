#ifndef _MODEL_H_
#define _MODEL_H_

#define S_SIZE 46 // (41 + 3 /*fluo*/ + 2 /*means*/)
#define C_SIZE 157 // (136 + 18 /*scalers*/ + 3 /*fluo*/)
#define A_SIZE 120

void initialize_states_default(double *STATES, const double *CONSTANTS);

void initialize_constants_default(double* CONSTANTS);

void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

void computeVariables(double VOI, double* CONSTANTS, double* STATES, double* ALGEBRAIC);

void _calc_means(double *STATES, const double *CONSTANTS);


#endif // _MODEL_H_
