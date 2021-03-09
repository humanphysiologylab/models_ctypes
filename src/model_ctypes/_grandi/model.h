#ifndef _MODEL_H_
#define _MODEL_H_

#define S_SIZE 41
#define C_SIZE 136
#define A_SIZE 120

void initialize_states_default(double *STATES);

void initialize_constants_default(double* CONSTANTS);

void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

void computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);


#endif // _MODEL_H_
