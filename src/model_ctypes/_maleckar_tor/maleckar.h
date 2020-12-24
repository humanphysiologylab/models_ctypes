#ifndef _MALECKAR_H_
#define _MALECKAR_H_

// #define S_SIZE 30 // 39
// #define C_SIZE 51 // 52
// #define A_SIZE 70 // 71

#define S_SIZE 39
#define C_SIZE 52
#define A_SIZE 71

void initialize_states_default(double *STATES);

void initialize_constants_default(double *CONSTANTS);

void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

void computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

#endif // _MALECKAR_H_
