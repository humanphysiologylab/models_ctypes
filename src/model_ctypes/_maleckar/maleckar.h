#ifndef _MALECKAR_H_
#define _MALECKAR_H_

#define S_SIZE 34
#define C_SIZE 54
#define A_SIZE 70

void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

void initConsts(double* CONSTANTS, double *STATES);

void computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC);

void computeVariablesVec(const int states_input_length, double* constants, double* states_input, double* algebraic_output);

void _calc_means(double *STATES, const double *CONSTANTS);

#endif // _MALECKAR_H_
