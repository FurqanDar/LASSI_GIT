#ifndef _CLUSTER_H_ // include guard
#define _CLUSTER_H_

#include "global.h"

int Clus_Network_ChainCluster_General(int const chainID);

int Clus_Network_ChainCluster_ForTotal(int const chainID);

int Clus_Network_LimitedCluster(int const chainID);

void Clus_Network_Distribution_Avg(void);

void Clus_Network_Distribution_MolWise_Avg(void);

int Clus_Network_SecondLargestCluster(void);

void Clus_Network_TotalAnalysis(void);

void Clus_Network_MolWise_LargestClusters(void);

int Clus_Proximity_ChainCluster_ForTotal_All(int const chainID);

int Clus_Proximity_ChainCluster_ForTotal_IntOnly(int const chainID);

void Clus_Proximity_TotalAnalysis(void);

int Clus_Proximity_SecondLargestCluster(void);

int Clus_Proximity_LimitedCluster_IntOnly(int const chainID);

int Clus_Proximity_LimitedCluster_IntOnly_Check(int const chainID, int const* OldList);

int Clus_Proximity_LimitedCluster_All(int const chainID, int* clusList);

int Clus_Proximity_LimitedCluster_All_Check(int const chainID, int const* OldList, int* NewList);

void Clus_Proximity_Distribution_Avg(void);

void Clus_Proximity_Distribution_IntOnly_MolWise_Avg(void);

void Clus_Proximity_Distribution_All_MolWise_Avg(void);

void Clus_Proximity_IntOnly_MolWise_LargestClusters(void);

void Clus_Perform_Analysis(void);

int Clus_Perform_ChainCluster_ForTotal(int const chainID);

void Clus_Perform_MolWise_LargestClusters(void);

void Clus_Find_LargestClusters(void);

int ClusUtil_GetOvlpNeighborBeads_ForBead(const int beadID, int* restrict neighList);

int ClusUtil_GetOvlpNeighborChains_ForBead(const int beadID, int* restrict neighList);

int ClusUtil_AddOvlpCluster_OfBead(const int beadID, char* restrict nTotClusTable, int* restrict clusList,
                                   int* restrict clusSize);

int ClusUtil_AddOvlpCluster_OfChain(const int chainID, char* restrict nTotClusTable, int* restrict clusList,
                                    int* restrict clusSize);

int Clus_Ovlp_OfChain(const int chainID, char* restrict nTotClusTable, int* restrict clusList);

#endif // _CLUSTER_H_
