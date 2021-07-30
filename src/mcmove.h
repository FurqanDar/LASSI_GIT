#ifndef _MCMOVE_H_ // include guard
#define _MCMOVE_H_

#include "global.h"

int MC_Step(float fMCTemp);

int MC_Step_Equil(float fMCTemp);

int Move_Rot(int beadID, float MyTemp);

int Move_Local(int beadID, float MyTemp);

int Move_Snake(int chainID, float MyTemp);

int Move_Trans(int chainID, float MyTemp);

int Move_Clus_Network(float MyTemp);

int Move_SmallClus_Network(int chainID, float MyTemp);

int Move_DbPvt(int beadID);

int Move_CoLocal(int beadID, float MyTemp);

int Move_MultiLocal(int beadID, float MyTemp);

int Move_Pivot(int chainID, float MyTemp);

int Move_BranchedRot(int chainID, float MyTemp);

int Move_SmallClus_Proximity(int chainID);

int Move_Local_Equil(int beadID, float MyTemp);

int Move_Snake_Equil(int chainID, float MyTemp);

int Move_Trans_Equil(int chainID, float MyTemp);

int Move_MultiLocal_Equil(int beadID, float MyTemp);

int Move_Pivot_Equil(int chainID, float MyTemp);

int Move_BranchedRot_Equil(int chainID, float MyTemp);

int Check_ChainDisp(const int chainID, const int *vec_disp);

int Check_MoveBeadTo(int *newPos);

void OP_DispChain(int chainID, const int *movR);

void OP_DispChain_ForTrans(int chainID, const int *movR);

void OP_RestoreChain(int chainID);

void OP_RestoreChain_ForTrans(int chainID);

void OP_CopyBead(int* copy_arr, const int* orig_arr);

void OP_CopyBeadsToOld(const int firstB, const int lastB);

void OP_RestoreBeadsFromOld(const int firstB, const int lastB);

void OP_RestoreChain_ForSnake(const int fB, const int lB);

void OP_MoveBeadTo(const int beadID, const int *newPos);

void OP_Inv_MoveBeadTo(int beadID);

void OP_MoveBeadTo_ForMTLocal(int beadID, const int *newPos);

void OP_SwapBeads(int bead1, int bead2);

void OP_Rotation(int PivotM, int beadID, int *tmpR);

void Rot_X_90(int beadID, const int tmpR[]);

void Rot_X_180(int beadID, const int tmpR[]);

void Rot_X_270(int beadID, const int tmpR[]);

void Rot_Y_90(int beadID, const int tmpR[]);

void Rot_Y_180(int beadID, const int tmpR[]);

void Rot_Y_270(int beadID, const int tmpR[]);

void Rot_Z_90(int beadID, const int tmpR[]);

void Rot_Z_180(int beadID, const int tmpR[]);

void Rot_Z_270(int beadID, const int tmpR[]);

void OP_ShuffleRotIndecies(void);

void OP_ShuffleArray(const int arr_size, int *dum_arr);

int Check_RotStates_wNeighList(int const beadID, int const resi, const int *neighList, int const neighNum);

int Check_RotStatesOld(int const beadID, int const resi, float const MyTemp);

int Check_RotStatesNew(int const beadID, int const resi, float const MyTemp);

void OP_NormalizeRotState(int beadVal, int CandNums);

int OP_PickRotState(int CandNums);

lLDub OP_GenMHValue(lLDub fRos, lLDub bRos, lLDub Delta_En, lLDub Cur_Temp);

#endif // _MCMOVE_H_
