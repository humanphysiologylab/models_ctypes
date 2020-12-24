#ifndef _ORD_H_
#define _ORD_H_

#define S_SIZE 49
#define C_SIZE 206
#define A_SIZE 200

void initConsts(double* CONSTANTS, double *STATES);

void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

void computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

#endif // _ORD_H_

