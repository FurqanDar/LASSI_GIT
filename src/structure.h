#ifndef _STRUCTURE_H_   // include guard
#define _STRUCTURE_H_

#include "global.h"

int Lat_Ind_FromCoords(int i, int j, int k);

int Lat_Ind_FromVec(int *xArr);

int Lat_Ind_OfBead(int beadID);

float Dist_BeadToBead(int bead1, int bead2);

float Dist_PointTotPoint_Float(float *f1, float *f2);

float Dist_PointToPoint(int *f1, int *f2);

float Dist_VecMag(const int *f1);

int Check_LinkerConstraint(int beadID, int *tmpR);

int Check_MTLinkerConstraint(int beadID, int (*tmpR)[POS_MAX]);

int Check_System_Structure(void);

void RDF_ComponentWise_Avg(void);

int RDF_ComponentIndex(const int i, const int j);

int RDFArr_Index(const int run_cycle, const int rdf_comp, const int x_pos);

int RadDenArr_Index(const int run_cycle, const int chain_type, const int x_pos);

int MolClusArr_Index(const int run_cycle, const int chain_type, const int clus_size);

void Calc_SystemCenterOfMass(int *tmpR);

void Calc_SystemCenterOfMass_OfMolType(int *tmpR, int const thisType);

void Calc_SystemCenterOfMass_WithoutMolType(int *tmpR, int const thisType);

void RadDen_MolTypeWise_Avg(void);

void GyrTensor_GyrRad_Avg(void);

#endif // _STRUCTURE_H_
