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

int *Create1DInt(const int xDim, const char* ArrName);

long *Create1DLong(const int xDim, const char* ArrName);

float *Create1DFloat(const int xDim, const char* ArrName);

long double *Create1DLongdouble(const int xDim, const char* ArrName);

int **Create2DInt(const int xDim, const int yDim, const char* ArrName);

long double **Create2DLongdouble(const int xDim, const int yDim, const char* ArrName);

#endif // _INITIALIZE_H_
