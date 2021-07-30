#ifndef _ENERGY_H_ // include guard
#define _ENERGY_H_

#include <stdio.h>
#include <stdlib.h>

void Energy_Total_System(void);

float Energy_Anisotropic(const int beadID);

float Energy_Anisotropic_Self(const int beadID);

float Energy_Anisotropic_For_Chain(const int beadID);

float Energy_Anisotropic_Contiguous_Range(const int beadID, const int smallest_bead, const int largest_bead);

float Energy_Anisotropic_With_List(const int beadID, const int *bead_list, const int list_size);

float Energy_Iso_Ovlp(int const beadType1, int const beadType2, float const xDis);

float Energy_Iso_Cont(int const beadType1, int const beadType2, float const xDis);

float Energy_Iso_fSol(int const beadType);

float Energy_OfOvlp_wNeighList(int const beadID, const int *neighList, int const neighNum);

float Energy_OfCont_wNeighList(int const beadID, const int *neighList, int const neighNum);

float Energy_ofSol_wNeighList(const int *neighList, int const neighNum);

float Energy_ofPairs_wNeighList(int const beadID, const int *neighList, int const neighNum);

float Energy_Isotropic(const int beadID);

float Energy_Isotropic_Self(const int beadID);

float Energy_Isotropic_For_Chain(const int beadID);

float Energy_Isotropic_Contiguous_Range(const int beadID, const int smallest_bead, const int largest_bead);

float Energy_Isotropic_With_List(const int beadID, const int *bead_list, const int list_size);

float Energy_Of_Chain(const int chainID);

float Energy_Of_Chain_Self(const int chainID);

float Energy_InitPotential(const int beadID);

void Energy_Iso_ForLocal(const int beadID, const int resi, const int* r_pos0,
                         long double *oldEn, long double *newEn,
                         int *ovlp_num, int *cont_num, int *ovlp_neighs, int *cont_neighs);

#endif // _ENERGY_H_
