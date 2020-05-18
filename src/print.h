#ifndef _PRINT_H_   // include guard
#define _PRINT_H_

#include <stdio.h>
#include <stdlib.h>

long TrajArr_Index(const int beadID, const int nFrameNumber, const int beadProp);

void Print_LogToScreen(long nGen, int run_it);

void Write_MCMove(char *filename, long nGen, float fMCTemp);

void Write_Energy(char *filename, long nGen);

void Print_Key(void);

void HandleTrajectory(char *fileStruct, const int run_it, const long nGen);

void Write_Trajectory(char *filename, long nGen);

void Save_Trajectory(const long nGen, const long curFrame);

void Write_Saved_Trajectory(char *filename, const int run_it);

void Write_TopFile(char *filename);

void Write_SysProp(char *filename);

void Write_TotalSysProp(char *filename, int run_it);

void Print_Data(long nGen, int run_it);

void Copy_Data(int run_it);


#endif // _PRINT_H_
