#ifndef _CLUSTER_H_   // include guard
#define _CLUSTER_H_


int Clus_ChainNetwork_General(int chainID);

int Clus_ChainNetwork_ForTotal(int chainID);

int Clus_LimitedNetworkCluster(int chainID);

void Clus_Distribution_Avg(void);

void Clus_DistributionMolWise_Avg(void);

int Clus_SecondLargestCluster(void);

void Clus_TotalAnalysis(void);

void Clus_MolWiseLargestCluster(void);

int Clus_ChainProximity_ForTotal(int chainID);

int Clus_LimitedProximityCluster(int const chainID);

int Clus_LimitedProximityCluster_Check(int const chainID, int const *OldList);

#endif // _CLUSTER_H_
