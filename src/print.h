#ifndef _PRINT_H_ // include guard
#define _PRINT_H_

#include "global.h"

long TrajArr_Index(const int beadID, const int nFrameNumber, const int beadProp);

void Print_LogToScreen(long nGen, int run_it);

void FileIO_WriteTo_MCMoveFile(char* filename, const long nGen, const float fMCTemp);

void PrintToScreen_KeyFile(void);

void PrintToScreen_SystemEnergy(void);

void PrintToScreen_AcceptanceRatios(void);

void PrintToScreen_Log_Thermalization(const long nGen);

void PrintToScreen_Log_FullRun(const long nGen, const int run_cycle);

char ForPrinting_GetReportState(const long nGen, const long thisReport);

void DataPrinting_Thermalization(const long nGen);

void HandleTrajectory(char* fileStruct, const int run_it, const long nGen);

void Write_Trajectory(const char* filename, const long nGen);

void Save_Trajectory(const long nGen, const long curFrame);

void Write_Saved_Trajectory(char* filename, const int run_it);

void Write_SysProp(char* filename);

void Print_Data(const long nGen, const int run_it);

void Copy_Data(int run_it);

void CopyData_RDF(const int run_it);

void CopyData_COMDen(const int run_it);

void CopyData_Clus(const int run_it);

void FileIO_CreateFile(const char* fileName);

void FileIO_CreateRunningDataFiles(void);

void FileIO_WriteTo_EnergyFile(char* filename, long nGen);

void FileIO_WriteTo_TopFile(const char* filename);

void FileIO_Write_TotalSysProp(const int run_it);

#endif // _PRINT_H_
