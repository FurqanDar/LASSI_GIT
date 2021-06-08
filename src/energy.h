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

float Energy_OVLP(const float fEn, const float xDis);

float Energy_CONT(const float fEn, const float xDis);

float Energy_Isotropic(const int beadID);

float Energy_Isotropic_Self(const int beadID);

float Energy_Isotropic_For_Chain(const int beadID);

float Energy_Isotropic_Contiguous_Range(const int beadID, const int smallest_bead, const int largest_bead);

float Energy_Isotropic_With_List(const int beadID, const int *bead_list, const int list_size);

float Energy_Of_Chain(const int chainID);

float Energy_Of_Chain_Self(const int chainID);

float Energy_InitPotential(const int beadID);

#endif // _ENERGY_H_
