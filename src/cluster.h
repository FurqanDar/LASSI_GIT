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

int Clus_Proximity_LimitedCluster_All(int const chainID);

int Clus_Proximity_LimitedCluster_All_Check(int const chainID, int const* OldList);

void Clus_Proximity_Distribution_Avg(void);

void Clus_Proximity_Distribution_IntOnly_MolWise_Avg(void);

void Clus_Proximity_Distribution_All_MolWise_Avg(void);

void Clus_Proximity_IntOnly_MolWise_LargestClusters(void);

void Clus_Perform_Analysis(void);

int Clus_Perform_ChainCluster_ForTotal(int const chainID);

void Clus_Perform_MolWise_LargestClusters(void);

void Clus_Find_LargestClusters(void);

#endif // _CLUSTER_H_
