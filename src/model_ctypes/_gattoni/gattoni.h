#ifndef _GATTONI_H_
#define _GATTONI_H_

#define S_SIZE 18
#define C_SIZE 71
#define A_SIZE 85

void initialize_states_default(double* STATES);

void initialize_constants_default(double* CONSTANTS);

void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

void computeVariables(double VOI, double* CONSTANTS, double* STATES, double* ALGEBRAIC);

#endif // _GATTONI_H_
