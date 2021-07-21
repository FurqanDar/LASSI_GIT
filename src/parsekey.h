#ifndef _PARSEKEY_H_   // include guard
#define _PARSEKEY_H_

#include <stdio.h>
#include <string.h>

int Parse_Keyfile(char *filename);

int Parse_EnergyFile(char *strEnFile);

void Parse_StructureFile(char *filename);

void Parse_StructureFile_CalcBeadsAndChains(const char* filename,
                                            int* n_bead_num, int* n_chain_num,
                                            int* n_bead_types, int* n_chain_types);

#endif // _PARSEKEY_H_
