#ifndef _INITIALIZE_H_   // include guard
#define _INITIALIZE_H_

void Memory_Initialization_AtStart(void);

void Memory_VerifyMalloc(void);

void Global_Array_Initialization_AtStart(void);

void Reset_Global_Arrays(void);

void Initial_Conditions_Simple(void);

void Initial_Conditions_FromFile(void);

void Initial_Conditions_BreakBonds(void);

float Temperature_Function(int mode, long nGen);

void Calculate_Rot_Bias(float CurrentTemp);

int *Array_Create_1D_int(const int xDim, const char* ArrName);

long *Array_Create_1D_long(const int xDim, const char* ArrName);

float *Array_Create_1D_float(const int xDim, const char* ArrName);

long double *Array_Create_1D_longdouble(const int xDim, const char* ArrName);

int **Array_Create_2D_int(const int xDim, const int yDim, const char* ArrName);

long double **Array_Create_2D_longdouble(const int xDim, const int yDim, const char* ArrName);

#endif // _INITIALIZE_H_
