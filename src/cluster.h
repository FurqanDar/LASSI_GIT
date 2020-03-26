#ifndef _CLUSTER_H_   // include guard
#define _CLUSTER_H_


int Clus_Netwrok_ChainCluster_General(int chainID);

int Clus_Network_ChainCluster_ForTotal(int chainID);

int Clus_Network_LimitedCluster(int chainID);

void Clus_Network_Distribution_Avg(void);

void Clus_Network_Distribution_MolWise_Avg(void);

int Clus_Network_SecondLargestCluster(void);

void Clus_Network_TotalAnalysis(void);

void Clus_Network_MolWise_LargestCluster(void);

int Clus_Proximity_ChainCluster_ForTotal(int chainID);

void Clus_Proximity_TotalAnalysis(void);

int Clus_Proximity_SecondLargestCluster(void);

int Clus_Proximity_LimitedCluster(int const chainID);

int Clus_Proximity_LimitedCluster_Check(int const chainID, int const *OldList);

void Clus_Proximity_Distribution_Avg(void);

void Clus_Proximity_Distribution_MolWise_Avg(void);

void Clus_Proximity_MolWise_LargestCluster(void);

#endif // _CLUSTER_H_
