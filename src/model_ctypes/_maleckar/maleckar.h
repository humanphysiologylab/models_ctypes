#ifndef _MALECKAR_H_
#define _MALECKAR_H_

#define S_SIZE 30
#define C_SIZE 51
#define A_SIZE 70

void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

void initConsts(double* CONSTANTS, double *STATES);

void computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

#endif // _MALECKAR_H_
