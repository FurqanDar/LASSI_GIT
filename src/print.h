#ifndef _PRINT_H_ // include guard
#define _PRINT_H_

#include "global.h"

long TrajArr_Index(const int beadID, const int nFrameNumber, const int beadProp);

void Print_LogToScreen(long nGen, int run_it);

void ScreenIO_Print_KeyFile(void);

void ScreenIO_Print_SystemEnergy(void);

void ScreenIO_Print_AcceptanceRatios(void);

void ScreenIO_Print_Log_Thermalization(const long nGen);

void ScreenIO_Print_Log_FullRun(const long nGen, const int run_cycle);

char ForPrinting_GetReportState(const long nGen, const long thisReport);

void DataPrinting_Thermalization(const long nGen);

void DataPrinting_DuringRunCycles(const long nGen, const int run_it);

void DataAnalysis_DuringRunCycles(const long nGen, const int run_it);

void Save_Trajectory(const long nGen, const long curFrame);

void Write_Saved_Trajectory(char* filename, const int run_it);

void Write_SysProp(char* filename);

void Print_Data(const long nGen, const int run_it);

void CopyData_All(const int run_it);

void CopyData_RDF(const int run_it);

void CopyData_COMDen(const int run_it);

void CopyData_Clus(const int run_it);

void FileIO_CreateFile(const char* fileName);

void FileIO_CreateRunningDataFiles(void);

void FileIO_WriteTo_MCMoveFile(const char* filename, const long nGen, const float fMCTemp);

void FileIO_PreCycle_Init(const int run_it);

void FileIO_AppendEnergyTo_EnergyFile(const char* fileNameStr, const long nGen);

void FileIO_WriteTo_TopFile(const char* filename);

void FileIO_HandleTrajectory(const char* fileNameStr, const int run_it, const long nGen);

void FileIO_AppendTrajFrame_ToFile(const char* filename, const long nGen);

void FileIO_WriteRestart_ForThermalization(void);

void FileIO_WriteRestart_ForRun(const int run_it);

void FileIO_Write_TotalSysProp(const int run_it);

#endif // _PRINT_H_
