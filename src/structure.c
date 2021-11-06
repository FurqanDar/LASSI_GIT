#include "structure.h"
#include "cluster.h"

/// Lat_Ind_FromCoords - helper function to get the correct 1D index of this position
/// \param i
/// \param j
/// \param k
/// \return the 1D index location for (i,j,k) position
int Lat_Ind_FromCoords(const int i, const int j, const int k)
{ // Lattice index from 3D to 1D array
    return i + nBoxSize[POS_X] * (j + nBoxSize[POS_Y] * k);
}

/// Lat_Ind_FromVec - returns the 1D index given the array xArr
/// \param xArr
/// \return
int Lat_Ind_FromVec(const int* xArr)
{ // Just the vector form of the above function. Easier to read sometimes.
    return xArr[POS_X] + nBoxSize[POS_X] * (xArr[POS_Y] + nBoxSize[POS_Y] * xArr[POS_Z]);
}

/// Lat_Ind_OfBead - returns the 1D index of this bead's location
/// \param beadID
/// \return
int Lat_Ind_OfBead(const int beadID)
{
    return bead_info[beadID][POS_X] +
           nBoxSize[POS_X] * (bead_info[beadID][POS_Y] + nBoxSize[POS_Y] * bead_info[beadID][POS_Z]);
}

// Note that the distance functions account for periodic boundaries

/// Dist_PointTotPoint_Float - euclidean distance between the vectors (arrays) f1 and f2 where f1 and f2 are floats.
/// \param f1
/// \param f2
/// \return
float Dist_PointToPoint_Float(const float* f1, const float* f2)
{
    float d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = fabsf(f1[i] - f2[i]);
            d[i] = d[i] > (float) nBoxSize[i] / 2. ? (float) nBoxSize[i] - d[i] : d[i];
        }

    return sqrtf(d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]);
}

/// Dist_PointToPoint - euclidean distance between the vectors (arrays) f1 and f2 where f1 and f2 are integer points.
/// \param f1
/// \param f2
/// \return
float Dist_PointToPoint(const int* f1, const int* f2)
{
    int d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = abs(f1[i] - f2[i]);
            d[i] = d[i] > nBoxSize[i] / 2 ? nBoxSize[i] - d[i] : d[i];
        }

    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToPoint - euclidean distance between beadID and the vector f1.
/// \param beadID
/// \param f1
/// \return
float Dist_BeadToPoint(const int beadID, const int* f1)
{
    int d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = abs(bead_info[beadID][i] - f1[i]);
            d[i] = d[i] > nBoxSize[i] / 2 ? nBoxSize[i] - d[i] : d[i];
        }
    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToPoint_Double - euclidean distance between beadID and the double array f1.
/// \param beadID
/// \param f1
/// \return
float Dist_BeadToPoint_Double(const int beadID, const lDub* f1)
{
    lDub d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = fabs((lDub) bead_info[beadID][i] - f1[i]);
            d[i] = d[i] > (lDub) nBoxSize[i] / 2. ? (lDub) nBoxSize[i] - d[i] : d[i];
        }
    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToPoint_Float - euclidean distance between beadID and the float array f1.
/// \param beadID
/// \param f1
/// \return
float Dist_BeadToPoint_Float(const int beadID, const float* f1)
{
    float d[POS_MAX];
    int i;
    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = fabsf((float) bead_info[beadID][i] - f1[i]);
            d[i] = d[i] > (float) nBoxSize[i] / 2. ? (float) nBoxSize[i] - d[i] : d[i];
        }
    return sqrtf((d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Dist_BeadToBead - euclidean distance between the two beads.
/// \param n1
/// \param n2
/// \return
float Dist_BeadToBead(const int n1, const int n2)
{
    lInt d[POS_MAX];
    lInt i;

    for (i = 0; i < POS_MAX; i++)
        {
            d[i] = abs(bead_info[n1][i] - bead_info[n2][i]);
            d[i] = d[i] > nBoxSize[i] / 2 ? nBoxSize[i] - d[i] : d[i];
        }

    return sqrtf((float) (d[POS_X] * d[POS_X] + d[POS_Y] * d[POS_Y] + d[POS_Z] * d[POS_Z]));
}

/// Check_System_Structure - performs a intensive and extensive sanity check. Check:
/// 1. All distances between beads in a molecule are legal.
/// 2. All bonds are symmetric.
/// 3. Bonds are only between possibly interacting beads.
/// 4. Beads' locations correspond with actual lattice location.
/// 5. Bonds are not between beads that are too far apart.
/// \return 0 if everything is okay, beadID+1 if failed.
int Check_System_Structure(void)
{
    int i, j;          // Looping variables
    int idx;           // Internal iterators for covalent bonds.
    int tmpR[POS_MAX]; // Just a vector to store coordinates
    int bondPart;
    for (i = 0; i < tot_beads; i++)
        {
            idx      = 0;
            bondPart = topo_info[i][idx];
            for (j = 0; j < POS_MAX; j++)
                {
                    tmpR[j] = bead_info[i][j];
                }
            if (naTotLattice[Lat_Ind_FromVec(tmpR)] == -1)
                {
                    printf("Lattice Position for bead %d is empty! Chain: %d\n", i, bead_info[i][BEAD_CHAINID]);
                    return i + 1;
                }
            if (i - naTotLattice[Lat_Ind_FromVec(tmpR)] != 0)
                { // This means there is a mismatch between where the bead is and where the lattice thinks the bead is
                    printf("Bead position and lattice value not the same. Crashing\t\t");
                    printf("B1:%d B2:%d\t C1:%d C2:%d\n", i, naTotLattice[Lat_Ind_FromVec(tmpR)],
                           bead_info[i][BEAD_CHAINID], bead_info[naTotLattice[Lat_Ind_FromVec(tmpR)]][BEAD_CHAINID]);
                    return i + 1;
                }
            while (topo_info[i][idx] != -1 && idx < MAX_BONDS)
                {
                    bondPart = topo_info[i][idx];
                    if (Dist_BeadToBead(i, bondPart) > LINKER_RSCALE * linker_len[i][idx])
                        {
                            printf("Bad beads! %d\t(%d %d %d)\t\tTopo:(%d %d %d)\t\tLinkers:(%.5f\t%.5f\t%.5f)\n", i,
                                   bead_info[i][0], bead_info[i][1], bead_info[i][2], topo_info[i][0], topo_info[i][1],
                                   topo_info[i][2], (float) linker_len[i][0], (float) linker_len[i][1],
                                   (float) linker_len[i][2]);
                            printf("\t\t\t\t\t-------------------->\t\t%f\tSHOULD BE\t%f\n",
                                   Dist_BeadToBead(i, bondPart), LINKER_RSCALE * (float) linker_len[i][idx]);
                            printf("Bad beads! %d\t(%d %d %d)\t\tTopo:(%d %d %d)\t\tLinkers:(%.5f\t%.5f\t%.5f)\n\n",
                                   bondPart, bead_info[bondPart][0], bead_info[bondPart][1], bead_info[bondPart][2],
                                   topo_info[bondPart][0], topo_info[bondPart][1], topo_info[bondPart][2],
                                   (float) linker_len[bondPart][0], (float) linker_len[bondPart][1],
                                   (float) linker_len[bondPart][2]);
                            return i + 1;
                        }
                    idx++;
                }
            if (bead_info[i][BEAD_FACE] != -1)
                {
                    if (nBeadTypeIsSticker[bead_info[i][BEAD_TYPE]] == 0)
                        {
                            printf("This bead -- %d -- should not have a bond. Crashing.\n", i);
                            return i + 1;
                        }
                    if (bead_info[i][BEAD_FACE] == i)
                        {
                            printf("Self bonded.\n");
                            return i + 1;
                        }
                    if (i != bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE])
                        {
                            printf("Bad bond!\n\t%d %d %d %f\nCrashing.\n", i, bead_info[i][BEAD_FACE],
                                   bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE],
                                   fEnergy[bead_info[i][BEAD_TYPE]][bead_info[bead_info[i][BEAD_FACE]][BEAD_TYPE]]
                                          [E_SC_SC]);
                            return i + 1;
                        }
                    if (Dist_BeadToBead(i, bead_info[i][BEAD_FACE]) > LINKER_RSCALE)
                        {
                            printf("Bad bond! Distance is wrong\n\t%d %d %d\n Distance is %f. Crashing.\n", i,
                                   bead_info[i][BEAD_FACE], bead_info[bead_info[i][BEAD_FACE]][BEAD_FACE],
                                   Dist_BeadToBead(i, bead_info[i][BEAD_FACE]));
                            return i + 1;
                        }
                }
        }
    return 0;
}

/// Dist_Vec3n - non periodic boundary euclidean magnitude of vector
/// \param f1: The array where indicies 0,1 and 3 correspond to x y and z.
/// \return
float Dist_Vec3n(const int* f1)
{ // Outputs the magnitude of the vector
    return sqrtf((float) (f1[0] * f1[0] + f1[1] * f1[1] + f1[2] * f1[2]));
}

/// Dist_VecMagSq - non periodic boundary euclidean magnitude of vector
/// \param f1: The array where indicies 0,1 and 3 correspond to x y and z.
/// \return Square of magnitude: x^2 + y^2 + z^2
int Dist_VecMagSq(const int* f1)
{ // Outputs the magnitude of the vector
    return (f1[0] * f1[0] + f1[1] * f1[1] + f1[2] * f1[2]);
}

/// GyrTensor_ClusterSpecific - calculates the total Gyration Tensor for a given cluster
/// \param ClusSize - the total size of the cluster.
/// \param ClusIndex - the index on naCluster where the cluster is stored.
/// THIS IS VERY OLD AND HASN'T BEEN LOOKED AT IN A WHILE
/// TODO: Update this for the new version
void GyrTensor_ClusterSpecific(int ClusSize, int ClusIndex)
{
    // Calculate the components of the gyration tensor for a given cluster.
    // ClusSize is the size of the cluster -- obviously -- whereas ClusIndex tell us
    // where in naCluster the chain indecies are located. naCluster[ClusIndex][0-ClusSize] is all the chainID's I need
    // for the calculation
    // Remember that the Gyration Tensor is a 3x3 symmetric object so we only need 6 numbers.
    int i, k, j, j2; // Basic indecies for loops
    for (i = 0; i < 7; i++)
        {
            fGyrTensor[i] = 0.;
        }                          // Initializing
    int firstB, lastB;             // Tracks the first and last bead of the given chain
    float tot_COM[POS_MAX] = {0.}; // This is where we shall store the COM of the cluster.
    int NumRes             = 0;    // Tracks how many residues are in this cluster

    // The only thing one needs to be careful about is to take PBC into account; the rest is tedium.
    // Use good old differential geometry to map each coordinate to two new coordinates, and unpack at the end.
    float theta[POS_MAX] = {0.};
    float zeta[POS_MAX]  = {0.}; // Extra coordinates for goodness
    float dumArg         = 0.;   // Just a dummy variable to be more efficient
    float dumArg2        = 0.;   // Another one
    // printf("Starting with COM\n");
    // Calculating the COM
    for (i = 0; i < ClusSize; i++)
        {
            firstB = chain_info[naCluster[ClusIndex][i]][CHAIN_START];
            lastB  = firstB + chain_info[naCluster[ClusIndex][i]][CHAIN_LENGTH];
            // printf("%d %d\n", firstB, lastB);
            // Just easier to track each chain like this
            for (k = firstB; k < lastB; k++)
                {
                    NumRes++; // Adding a residue to the total
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = 2. * M_PI * ((float) bead_info[k][j] / (float) nBoxSize[j]);
                            theta[j] += cosf(dumArg);
                            zeta[j] += sinf(dumArg); // Since I am taking the average just keep adding
                        }
                }
            // printf("Next bead\n");
        }
    // printf("Done with COM\n");
    // Calculating average, and then COM
    for (j = 0; j < POS_MAX; j++)
        {
            theta[j]   = theta[j] / (float) NumRes;
            zeta[j]    = zeta[j] / (float) NumRes;
            tot_COM[j] = atan2f(-theta[j], -zeta[j]) + M_PI;
            tot_COM[j] = nBoxSize[j] * (tot_COM[j] / 2. / M_PI);
        }

    // Using the COM to calculate the Gyration Tensor
    //  GyrTen_{ij} = 1/N sum_1^N (r_i-com_i)(r_j-com_j) so just take the sums and divide at the end
    for (i = 0; i < ClusSize; i++)
        {
            firstB = chain_info[naCluster[ClusIndex][i]][CHAIN_START];
            lastB  = firstB + chain_info[naCluster[ClusIndex][i]][CHAIN_LENGTH];
            // Just easier to track each chain like this
            for (k = firstB; k < lastB; k++)
                {
                    for (j = 0; j < POS_MAX; j++)
                        {
                            for (j2 = j; j2 < POS_MAX; j2++)
                                {
                                    dumArg = fabsf((float) bead_info[k][j] - tot_COM[j]) <
                                                     (float) nBoxSize[j] - fabsf((float) bead_info[k][j] - tot_COM[j]) ?
                                                 fabsf((float) bead_info[k][j] - tot_COM[j]) :
                                                 (float) nBoxSize[j] - fabsf((float) bead_info[k][j] - tot_COM[j]);
                                    dumArg2 =
                                        fabsf((float) bead_info[k][j2] - tot_COM[j2]) <
                                                (float) nBoxSize[j2] - fabsf((float) bead_info[k][j2] - tot_COM[j2]) ?
                                            fabsf((float) bead_info[k][j2] - tot_COM[j2]) :
                                            (float) nBoxSize[j2] - fabsf((float) bead_info[k][j2] - tot_COM[j2]);
                                    fGyrTensor[j + 3 * (j2 - j)] += dumArg * dumArg2;
                                    // 0 = xx; 1 = yy; 2 = zz; 3 = xy; 6 = xz; 4 = yz; Need smarter indexing
                                }
                        }
                }
        }

    // Calculating the average;
    for (i = 0; i < 7; i++)
        {
            fGyrTensor[i] /= (float) NumRes;
        } // printf("%f\n",fGyrTensor[i]);}printf("\n");

    // exit(1);
}

/// GyrTensor_GyrRad - calculates the total Gyration Tensor of the system.
/// THIS IS VERY OLD AND HASN'T BEEN LOOKED AT IN A WHILE
/// TODO: Update this for the new version
void GyrTensor_GyrRad(void)
{                    // Calculates the gyration tensor for the whole system
    int i, k, j, j2; // Basic indecies for loops
    for (i = 0; i < 7; i++)
        {
            fGyrTensor[i] = 0.f;
        }                          // Initializing
    float tot_COM[POS_MAX] = {0.f}; // This is where we shall store the COM of the cluster.

    // The only thing one needs to be careful about is to take PBC into account; the rest is tedium.
    // Use good old differential geometry to map each coordinate to two new coordinates, and unpack at the end.
    float theta[POS_MAX] = {0.f};
    float zeta[POS_MAX]  = {0.f}; // Extra coordinates for goodness
    float dumArg         = 0.f;   // Just a dummy variable to be more efficient
    float dumArg2        = 0.f;   // Another one

    for (i = 0; i < tot_beads; i++)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    tot_COM[j] += bead_info[i][j];
                }
        }

    // Calculating average, and then COM
    for (j = 0; j < POS_MAX; j++)
        {
            tot_COM[j] /= (float) tot_beads;
        }
    // printf("\n");
    // Using the COM to calculate the Gyration Tensor
    //  GyrTen_{ij} = 1/N sum_1^N (r_i-com_i)(r_j-com_j) so just take the sums and divide at the end
    for (i = 0; i < tot_beads; i++)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    dumArg = (float) bead_info[i][j] - tot_COM[j];
                    for (j2 = j; j2 < POS_MAX; j2++)
                        {
                            dumArg2 = (float) bead_info[i][j2] - tot_COM[j2];
                            fGyrTensor[j + 3 * (j2 - j)] += dumArg * dumArg2;
                            // printf("%.2f\t", dumArg2);
                            // 0 = xx; 1 = yy; 2 = zz; 3 = xy; 6 = xz; 4 = yz; Need smarter indexing
                        }
                } // printf("\n");
        }
    // exit(1);
    // Calculating the average;
    for (i = 0; i < 7; i++)
        {
            fGyrTensor[i] /= (float) tot_beads;
        } // printf("%f\n",fGyrTensor[i]);}printf("\n");
}

/// GyrTensor_GyrRad_Avg - calculates the total radius of gyration of the system, while not being smart about the
/// periodic boundaries. This is used as a proxy to detect phase separation, but is a relic of the old formalism. The
/// RDF should be used in general. Although this can be used without the need for a non-interacting prior.
void GyrTensor_GyrRad_Avg(void)
{
    /*
    Only calculates the diagonals of the gyration tensor, and calculates the sum of the
    diagonals. Remember that Rg^2 = Tr(GyrTen) so we only need to calculate the diagonals, and then then sum.
    I shall borrow most of the code from above, and so read GyrTensor_ClusterSpecific for what's happening here.
    */
    int i, j; // Loop indecies
    lDub tot_COM[POS_MAX] = {0.};

    Calc_SystemCenterOfMass(tot_COM);
    int dumArg[POS_MAX] = {0};
    for (i = 0; i < 7; i++)
        { // Initializing to 0
            fGyrTensor[i] = 0.f;
        }

    for (i = 0; i < tot_beads; i++)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    dumArg[j] = abs(bead_info[i][j] - (int) tot_COM[j]);
                    dumArg[j] = dumArg[j] > nBoxSize[i] / 2 ? nBoxSize[i] - dumArg[j] : dumArg[j];
                    fGyrTensor[j] += (float) (dumArg[j] * dumArg[j]);
                }
        }

    // Adding to the total fSysGyrRad to be averaged at the end.
    fSysGyrRad += sqrtf((fGyrTensor[0] + fGyrTensor[1] + fGyrTensor[2]) / (float) tot_beads);
    nTotGyrRadCounter++; // Remembering that we have calculated the radius; for final averaging.
}

/// RDF_ComponentIndex - 1D index for the symmetric g_{ij} matrix given i and j.
/// \param i
/// \param j
/// \return The index of the array
/// The way it is set up, the indexing goes through the diagonal and then 0-1, 0-2, ... 0-N, 1-2, ... 1-N and so on
int RDF_ComponentIndex(const int i, const int j)
{
    if (i == j)
        {
            return i + 1;
        }
    else if (i < j)
        {
            return nBeadTypes + j - (i * (3 + i - 2 * nBeadTypes)) / 2;
        }
    else
        {
            return nBeadTypes + i - (j * (3 + j - 2 * nBeadTypes)) / 2;
        }
}

/// RDFArr_Index - 1D index for ld_TOTRDF_Arr which is used to globally store the different RDFs
/// \param run_cycle
/// \param rdf_comp
/// \param x_pos
/// \return 1D index for the totalRDFArray
int RDFArr_Index(const int run_cycle, const int rdf_comp, const int x_pos)
{
    return x_pos + nRDF_TotBins * (rdf_comp + nRDF_TotComps * run_cycle);
}

int RadDen_ComponentIndex(const int i, const int j)
{
    if (i < 0)
        {
            return j;
        }
    else
        {
            return tot_chain_types + j + tot_chain_types * i;
        }
}

int RadDenArr_Index(const int run_cycle, const int rad_comp, const int x_pos)
{
    return x_pos + nRDF_TotBins * (rad_comp + nRadDen_TotComps * run_cycle);
}

int MolClusArr_Index(const int run_cycle, const int chain_type, const int clus_size)
{
    return clus_size + tot_chains * (chain_type + tot_chain_types * run_cycle);
}

/// RDF_ComponentWise_Avg - calculates the pair-distribution of the system where every bead acts as the center
/// of a radial histogram of pairs. Note that dr = 1/4 lattice units.
void RDF_ComponentWise_Avg(void)
{
    /*
    Calculates the RDF and adds it all up so that it can be averaged out at the end of the run.
    */
    float x; // For distance
    int i, j, k;
    int resi, resj;
    int myBin = 0;
    int array_pos;

    // Calculating where, and how many, pairs exist
    for (i = 0; i < tot_beads; i++)
        {
            resi = bead_info[i][BEAD_TYPE];
            for (j = i + 1; j < tot_beads; j++)
                {
                    resj = bead_info[j][BEAD_TYPE];
                    x    = Dist_BeadToBead(i, j);
                    // Note that Dist_BeadToBead(i,j) automatically ensures no distance is greater than (L/2)*sqrt(3)
                    myBin = (int) floor(4. * x);                 // I am assuming for now that dr=1/4
                    ldRDF_Arr[RDFArr_Index(0, 0, myBin)] += 2.0; // Adding a pair to that bin
                    array_pos = RDF_ComponentIndex(resi, resj);
                    ldRDF_Arr[RDFArr_Index(0, array_pos, myBin)] += 2.0;
                }
        }
    nRDFCounter++;
}

/// Check_LinkerConstraint - if I move beadID to tmpR, do I still satisfy the linker lengths for beadID?
/// For the proposed location, we loop over all bonded partners of beadID and check if the distance is within the
/// linker_len distance.
/// \param beadID
/// \param tmpR
/// \return 1 means all is good, 0 means bad.
int Check_LinkerConstraint(const int beadID, const int* tmpR)
{
    // Check if the proposed new location for beadID is such that all the linkers are unbroken.
    int idx;         // Iterator to loop over bond Partners
    int bondPartner; // It is what it is.
    idx         = 0;
    bondPartner = topo_info[beadID][idx]; // Initializing the two.
    while (idx < MAX_BONDS && topo_info[beadID][idx] != -1)
        { // Keep going till we run out of partners
            bondPartner = topo_info[beadID][idx];
            if (Dist_PointToPoint(bead_info[bondPartner], tmpR) > LINKER_RSCALE * (float) linker_len[beadID][idx])
                {
                    return 0; // This means that we have broken one of the linkers.
                }
            idx++;
        }
    return 1; // This means that all linker constraints are satisfied.
}

/// Check_MTLinkerConstraint_OLD - if I move beadID and all of it's covalent bonded beads to the locations stored in
/// tmpR[][], would all the linker constraints be satisfied?
/// \param beadID
/// \param tmpR
/// \return
int Check_MTLinkerConstraint_OLD(int beadID, int (*tmpR)[POS_MAX])
{

    int curID = beadID;
    int idx, bPart;
    int j;
    int topIt = 0;
    int canI  = 1;

    while (curID != -1)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    bead_info[curID][j] = tmpR[topIt][j]; // Moving
                }
            curID = topo_info[beadID][topIt++];
        }

    curID = beadID;
    topIt = 0;
    while (curID != -1 && canI == 1)
        {
            idx   = 0;
            bPart = topo_info[curID][idx];
            while (bPart != -1 && idx < MAX_BONDS)
                {
                    if (Dist_BeadToBead(curID, bPart) > LINKER_RSCALE * (float) linker_len[curID][idx])
                        {
                            canI = 0;
                            break;
                        }
                    bPart = topo_info[curID][++idx];
                }
            curID = topo_info[beadID][topIt++];
        }
    curID = beadID;
    topIt = 0;
    while (curID != -1)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    bead_info[curID][j] = old_bead[curID][j]; // Moving back
                }
            curID = topo_info[beadID][topIt++];
        }

    return canI;
}

/// Check_LinkerConstraints_ForBeadList - given that all the beads in beadList are at the new
/// positions already, as part of the MTLocal move, we iterate over every bead and check if
/// any linker constraints are broken. If so, return 0. If not, return 1 for success.
/// \param beadID
/// \param tmpR
/// \return
int Check_LinkerConstraints_ForBeadList(const int listSize, const int* beadList)
{

    int i, j, bead1, bead2;
    int dumBonds[MAX_BONDS + 1];
    int dumNum;
    float xDis;
    for (i = 0; i < listSize; i++)
        {
            bead1  = beadList[i];
            dumNum = OP_GetTopoBonds(bead1, dumBonds);
            for (j = 1; j < dumNum; j++)
                { // No need for self-distance = 0.
                    bead2 = dumBonds[j];
                    xDis  = Dist_BeadToBead(bead1, bead2);
                    if (xDis > LINKER_RSCALE * (float) linker_len[bead1][j - 1])
                        {
                            return 0;
                        }
                }
        }

    return 1;
}

/// Check_BeadID_InList: Iterate over the list of beads to see if the bead is in the list.
/// If so, return 1. Otherwise return 0.
/// \param thisBeadID
/// \param listSize
/// \param beadList
/// \return
int Check_BeadID_InList(const int thisBeadID, const int listSize, const int beadList[MAX_BONDS + 1])
{
    int i;
    for (i = 0; i < listSize; i++)
        {
            if (thisBeadID == beadList[i])
                {
                    return 1;
                }
        }
    return 0;
}

void Calc_SystemCenterOfMass(lDub* tmpR)
{

    int i, j;                     // Iterators
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lDub) nBoxSize[j];
        }

    for (i = 0; i < tot_beads; i++)
        {
            for (j = 0; j < POS_MAX; j++)
                {
                    dumArg = dumConst[j] * (lDub) bead_info[i][j];
                    zeta[j] += sin(dumArg);
                    xi[j] += cos(dumArg);
                }
        }

    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lDub) tot_beads;
                    zeta[j] /= (lDub) tot_beads;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
            else
                {
                    tot_COM[j] = 0.;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = tot_COM[j];
        }
}

void Calc_CenterOfMass_OfCluster(lDub* tmpR, const int cluster_size, const int* ClusList)
{
    // This version measures the COM of a cluster of size c
    //  cluster size, given the molecule ID's in naList.
    // The COM from this is not necessarily the COM of the system as a whole.
    int thisMol, i, j, k;         // Iterators
    int fB, lB;                   // Keep track of first and last beads of a given molecule
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored
    int bead_total_now    = 0;

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lLDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lLDub) nBoxSize[j];
        }

    for (k = 0; k < cluster_size; k++)
        {
            thisMol = ClusList[k];
            fB      = chain_info[thisMol][CHAIN_START];
            lB      = fB + chain_info[thisMol][CHAIN_LENGTH];
            for (i = fB; i < lB; i++)
                {
                    bead_total_now++;
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = dumConst[j] * (lDub) bead_info[i][j];
                            zeta[j] += sin(dumArg);
                            xi[j] += cos(dumArg);
                        }
                }
        }
    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lLDub) bead_total_now;
                    zeta[j] /= (lLDub) bead_total_now;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = tot_COM[j];
        }
}

void Calc_SystemCenterOfMass_OfMolType(lDub* tmpR, const int thisType)
{
    // This version measures the COM of only type thisType
    // The COM from this is not necessarily the COM of the system as a whole.
    int thisMol, i, j, k;         // Iterators
    int fB, lB;                   // Keep track of first and last beads of a given molecule
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored
    int bead_total_now    = 0;

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lDub) nBoxSize[j];
        }

    for (k = 0; k < tot_chains; k++)
        {
            thisMol = chain_info[k][CHAIN_TYPE];
            if (thisMol != thisType)
                {
                    continue;
                }
            fB = chain_info[k][CHAIN_START];
            lB = fB + chain_info[k][CHAIN_LENGTH];
            for (i = fB; i < lB; i++)
                {
                    bead_total_now++;
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = dumConst[j] * (lDub) bead_info[i][j];
                            zeta[j] += sin(dumArg);
                            xi[j] += cos(dumArg);
                        }
                }
        }
    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lDub) bead_total_now;
                    zeta[j] /= (lDub) bead_total_now;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
        }

    /*for (j=0; j<POS_MAX; j++){
        tmpR[j] = ((int) tot_COM[j]);
    }*/
    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = tot_COM[j];
        }
}

void Calc_SystemCenterOfMass_WithoutMolType(lDub* tmpR, const int thisType)
{
    // This version measures the COM of everything except type thisType
    // The COM from this is not necessarily the COM of the system as a whole.
    int thisMol, i, j, k;         // Iterators
    int fB, lB;                   // Keep track of first and last beads of a given molecule
    lDub tot_COM[POS_MAX] = {0.}; // This is where the COM will be stored
    int bead_total_now    = 0;

    lDub zeta[POS_MAX] = {0.};
    lDub xi[POS_MAX]   = {0.};
    lDub dumArg        = 0.;

    lDub dumConst[POS_MAX] = {0.};
    for (j = 0; j < POS_MAX; j++)
        {
            dumConst[j] = 2. * M_PI / (lDub) nBoxSize[j];
        }

    for (k = 0; k < tot_chains; k++)
        {
            thisMol = chain_info[k][CHAIN_TYPE];
            if (thisMol == thisType)
                {
                    continue;
                }
            fB = chain_info[k][CHAIN_START];
            lB = fB + chain_info[k][CHAIN_LENGTH];
            for (i = fB; i < lB; i++)
                {
                    bead_total_now++;
                    for (j = 0; j < POS_MAX; j++)
                        {
                            dumArg = dumConst[j] * (lDub) bead_info[i][j];
                            zeta[j] += sin(dumArg);
                            xi[j] += cos(dumArg);
                        }
                }
        }
    int nCheck[POS_MAX] = {0};

    for (j = 0; j < POS_MAX; j++)
        {
            if ((zeta[j] == 0.) && (xi[j] == 0.))
                { // If both 0, then undefined, so skip.
                    nCheck[j] = 0;
                }
            else
                {
                    nCheck[j] = 1;
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            if (nCheck[j] == 1)
                {
                    xi[j] /= (lDub) bead_total_now;
                    zeta[j] /= (lDub) bead_total_now;
                    tot_COM[j] = atan2(-zeta[j], -xi[j]) + M_PI;
                    tot_COM[j] /= dumConst[j];
                }
        }

    for (j = 0; j < POS_MAX; j++)
        {
            tmpR[j] = tot_COM[j];
        }
}

void RadDen_Avg_MolTypeWise_FromSysCen(void)
{

    int i;        // Iterator for loop
    int thisType; // Tracks the type of the chain
    lDub sysCOM[POS_MAX] = {0.};
    int myBin;
    float xDis = 0.; // Tracks the distance between the COM and the specific bead

    Calc_SystemCenterOfMass(sysCOM);
    // printf("SYS COM = (%d, %d, %d) \n", sysCOM[0], sysCOM[1], sysCOM[2]);
    for (i = 0; i < tot_beads; i++)
        {
            thisType = bead_info[i][BEAD_CHAINID];
            thisType = chain_info[thisType][CHAIN_TYPE];
            xDis     = Dist_BeadToPoint_Double(i, sysCOM);
            myBin    = (int) floor(4. * xDis);
            ldRadDen_Arr[RadDenArr_Index(0, thisType, myBin)] += 1.0;
        }
    nRadDenCounter++;
}

void RadDen_Avg_MolTypeWise_FromMolTypeCen_Old_CorrectVersion(void)
{

    int i, j;     // Iterators for loop
    int thisType; // Tracks the type of the chain
    int thisComp; // Tracks which component of ldRadDen
    lDub typeCOM[POS_MAX] = {0};
    int COM_int[POS_MAX]  = {0.};
    int myBin;
    float xDis   = 0.; // Tracks the distance between the COM and the specific bead
    int cur_type = 0;

    Calc_SystemCenterOfMass(typeCOM);
    for (j = 0; j < POS_MAX; j++)
        {
            COM_int[j] = (int) typeCOM[j];
        }
    for (i = 0; i < tot_beads; i++)
        {
            thisType = bead_info[i][BEAD_CHAINID];
            thisType = chain_info[thisType][CHAIN_TYPE];
            thisComp = RadDen_ComponentIndex(-1, thisType);
            xDis     = Dist_BeadToPoint(i, COM_int);
            myBin    = (int) (4. * xDis);
            ldRadDen_Arr[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
        }

    for (cur_type = 0; cur_type < tot_chain_types; cur_type++)
        { // Go through each molType's center
            Calc_SystemCenterOfMass_OfMolType(typeCOM, cur_type);
            for (j = 0; j < POS_MAX; j++)
                {
                    COM_int[j] = (int) typeCOM[j];
                }
            for (i = 0; i < tot_beads; i++)
                {
                    thisType = bead_info[i][BEAD_CHAINID];
                    thisType = chain_info[thisType][CHAIN_TYPE];
                    thisComp = RadDen_ComponentIndex(cur_type, thisType);
                    xDis     = Dist_BeadToPoint(i, COM_int);
                    myBin    = (int) (4. * xDis);
                    ldRadDen_Arr[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
                }
        }

    /*int tmpBead;
    int cur_POS[POS_MAX] = {0.};
    int cur_DIS[POS_MAX] = {0.};
    int radRange;
    radRange = nBoxSize[0];
    thisComp = RadDen_ComponentIndex(2,2);
    for(j=0; j<POS_MAX; j++){
        COM_int[j] = 0;
    }
    for (i = 0; i < radRange; i++) {
        //cur_POS[POS_X] = (COM_int[POS_X] + i) % nBoxSize[POS_X];
        cur_DIS[POS_X] = i > nBoxSize[0]/2 ? nBoxSize[0] - i: i;
        cur_DIS[POS_X] *= cur_DIS[POS_X];
        for (j = 0; j < radRange; j++) {
            //cur_POS[POS_Y] = (COM_int[POS_Y] + j ) % nBoxSize[POS_Y];
            cur_DIS[POS_Y] = j > nBoxSize[0]/2 ? nBoxSize[0] - j: j;
            cur_DIS[POS_Y] *= cur_DIS[POS_Y];
            for (k = 0; k < radRange; k++) {
                //cur_POS[POS_Z] = (COM_int[POS_Z] + k) % nBoxSize[POS_Z];
                cur_DIS[POS_Z] = k > nBoxSize[0]/2 ? nBoxSize[0] - k: k;
                cur_DIS[POS_Z] *= cur_DIS[POS_Z];
                xDis = sqrtf((float)(cur_DIS[POS_X]+cur_DIS[POS_Y]+cur_DIS[POS_Z] ));
                myBin = (int) (4.*xDis);
                ldRadDen_Arr[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
            }
        }
    }
    */

    nRadDenCounter++;
}

void RadDen_Avg_MolTypeWise_FromMolTypeCen(void)
{

    int i, j, k;  // Iterators for loop
    int thisType; // Tracks the type of the chain
    int thisComp; // Tracks which component of ldRadDen
    lDub typeCOM[POS_MAX] = {0.};
    int COM_int[POS_MAX]  = {0};
    int myBin;
    float xDis   = 0.; // Tracks the distance between the COM and the specific bead
    int cur_type = 0;
    int fB, lB;
    int clus_size = 0;
    int clus_id   = -1;
    int* clus_id_list; // Tracks the cluster ID's after the clustering analysis.
    clus_id_list = malloc((tot_chain_types + 1) * sizeof(lInt));
    for (i = 0; i <= tot_chain_types; i++)
        {
            clus_id_list[i] = -1;
        }

    // Perform clustering analysis
    Clus_Perform_MolWise_LargestClusters();
    // Remember the cluster ID's
    for (i = 0; i <= tot_chain_types; i++)
        {
            clus_id_list[i] = naList[i];
        }

    for (cur_type = 0; cur_type <= tot_chain_types; cur_type++)
        {
            clus_id   = clus_id_list[cur_type];
            clus_size = naCluster[clus_id][0];
            // printf("(%d %d) ", clus_id, clus_size);
            for (i = 0; i < clus_size; i++)
                {
                    naList[i] = naCluster[clus_id][i + 1];
                }

            Calc_CenterOfMass_OfCluster(typeCOM, clus_size, naList);
            for (j = 0; j < POS_MAX; j++)
                {
                    COM_int[j] = (int) typeCOM[j];
                }

            for (k = 0; k < clus_size; k++)
                {
                    fB = chain_info[naList[k]][CHAIN_START];
                    lB = fB + chain_info[naList[k]][CHAIN_LENGTH];
                    for (i = fB; i < lB; i++)
                        {
                            thisType = bead_info[i][BEAD_CHAINID];
                            thisType = chain_info[thisType][CHAIN_TYPE];
                            thisComp = RadDen_ComponentIndex(cur_type - 1, thisType);
                            xDis     = Dist_BeadToPoint(i, COM_int);
                            myBin    = (int) (4. * xDis);
                            ldRadDen_Arr[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
                        }
                }
        }

    for (cur_type = 0; cur_type <= tot_chain_types; cur_type++)
        {
            if (cur_type == 0)
                {
                    Calc_SystemCenterOfMass(typeCOM);
                }
            else
                {
                    Calc_SystemCenterOfMass_OfMolType(typeCOM, cur_type - 1);
                }
            for (j = 0; j < POS_MAX; j++)
                {
                    COM_int[j] = (int) typeCOM[j];
                }

            for (i = 0; i < tot_beads; i++)
                {
                    thisType = bead_info[i][BEAD_CHAINID];
                    thisType = chain_info[thisType][CHAIN_TYPE];
                    thisComp = nRadDen_CompShift + RadDen_ComponentIndex(cur_type - 1, thisType);
                    xDis     = Dist_BeadToPoint(i, COM_int);
                    myBin    = (int) (4. * xDis);
                    ldRadDen_Arr[RadDenArr_Index(0, thisComp, myBin)] += 1.0;
                }
        }

    // printf("\n");

    nRadDenCounter++;

    free(clus_id_list);
}

void Calculate_Distances_For_Radius(float* thisList, const int nRad)
{
    int r_disp[POS_MAX] = {0};

    int dum_iter = 0;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            thisList[dum_iter++] = Dist_Vec3n(r_disp);
                        }
                }
        }
    //    printf("%d", dum_iter);
    //    exit(1);
}

/// NeighborSearch_ForOvlp: Given the beadID, which should be at startVec, we calculate all the neighbors
/// within the +-1 cube. We return the number of neighbors.
/// \param beadID
/// \param startVec
/// \param neighList
/// \return Number of neighbors.
int NeighborSearch_ForOvlp(int const beadID, const int* startVec, int* neighList)
{
    int neigh_num = 0;

    int tmpBead;
    int r_disp[POS_MAX] = {0};
    int r_chck[POS_MAX] = {0};
    const int nRad      = 1;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_X);
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Y);
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Z);
                            tmpBead = naTotLattice[Lat_Ind_FromVec(r_chck)];
                            if (tmpBead != -1 && tmpBead != beadID)
                                {
                                    neighList[neigh_num++] = tmpBead;
                                }
                        }
                }
        }
    neighList[neigh_num] = -1;
    return neigh_num;
}

/// NeighborSearch_ForOvlp: Given the beadID, which should be at startVec, we calculate all the neighbors
/// within the +-1 cube. We return the number of neighbors.
/// \param beadID
/// \param startVec
/// \param neighList
/// \return Number of neighbors.
int NeighborSearch_ForCont(int const beadID, const int* startVec, int* contList, int* ovlpList, int* ovlpNum)
{
    int neigh_num = 0;

    int tmpBead;
    int r_disp[POS_MAX] = {0};
    int r_chck[POS_MAX] = {0};
    const int nRad      = LARGEST_RADIUS;

    *ovlpNum = 0;

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_X);
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Y);
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC_ofComp(r_chck, startVec, r_disp, POS_Z);
                            tmpBead = naTotLattice[Lat_Ind_FromVec(r_chck)];
                            if (tmpBead != -1 && tmpBead != beadID)
                                {
                                    if ((abs(r_disp[0]) <= 1) && (abs(r_disp[1]) <= 1) && (abs(r_disp[2]) <= 1))
                                        {
                                            ovlpList[*ovlpNum] = tmpBead;
                                            *ovlpNum           = *ovlpNum + 1;
                                        }
                                    contList[neigh_num++] = tmpBead;
                                }
                        }
                }
        }
    contList[neigh_num] = -1;
    ovlpList[*ovlpNum]  = -1;
    return neigh_num;
}

/// NeighborSearch_EmptySitesAround_wRad: Given the starting point, and radius, we calculate how many sites
/// around this point are empty.
/// \param startVec
/// \param nRad
/// \return Number of empty lattice-sites around startVec
int NeighborSearch_EmptySitesAround_wRad(const int* startVec, const int nRad)
{
    int empty_num = 0;
    int tmpBead;
    int r_disp[POS_MAX]   = {0};
    int r_search[POS_MAX] = {0};

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC(r_search, startVec, r_disp);
                            tmpBead = naTotLattice[Lat_Ind_FromVec(r_search)];
                            if (tmpBead == -1)
                                {
                                    empty_num++;
                                }
                        }
                }
        }
    return empty_num;
}

/// NeighborSearch_SolSitesAround: number of solvation sites around this point. Finds the number of empty lattice
/// sites in the 26 closest sites.
/// \param startVec
/// \return Empty sites around this point.
int NeighborSearch_SolSitesAround(const int* startVec)
{
    return NeighborSearch_EmptySitesAround_wRad(startVec, 1);
}

int NeighborSearch_AroundPoint_UptoIndex(const int beadID, const int* startVec, const int upTO, int* neighList)
{
    int i;
    int neigh_num = 0;
    int tmpBead;

    int r_search[POS_MAX] = {0};

    int r_all[100][POS_MAX] = {0};

    for (i = 0; i < upTO; ++i)
        {
            LatPos_add_wPBC(r_search, startVec, r_all[i]);
            tmpBead = naTotLattice[Lat_Ind_FromVec(r_search)];
            if (tmpBead != -1)
                {
                    neighList[neigh_num++] = tmpBead;
                }
        }

    return neigh_num;
}

int NeighborSearch_AroundPoint_wRad_IgnBead(const int beadID, const int* startVec, const int nRad, int* neighList)
{

    int neigh_num = 0;
    int tmpBead;
    int r_disp[POS_MAX]   = {0};
    int r_search[POS_MAX] = {0};

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC(r_search, startVec, r_disp);
                            tmpBead = naTotLattice[Lat_Ind_FromVec(r_search)];
                            if (tmpBead != -1 && tmpBead != beadID)
                                {
                                    neighList[neigh_num++] = tmpBead;
                                }
                        }
                }
        }
    return neigh_num;
}

int NeighborSearch_AroundPoint_wRad_wDists(const int beadID, const int* startVec, const int nRad, int* neighList,
                                           float* distList)
{

    int neigh_num = 0;
    int tmpBead;
    int r_disp[POS_MAX] = {0};

    int r_search[POS_MAX] = {0};

    for (r_disp[0] = -nRad; r_disp[0] <= nRad; r_disp[0]++)
        {
            for (r_disp[1] = -nRad; r_disp[1] <= nRad; r_disp[1]++)
                {
                    for (r_disp[2] = -nRad; r_disp[2] <= nRad; r_disp[2]++)
                        {
                            LatPos_add_wPBC(r_search, startVec, r_disp);
                            tmpBead = naTotLattice[Lat_Ind_FromVec(r_search)];
                            if (tmpBead != -1 && tmpBead != beadID)
                                {
                                    neighList[neigh_num] = tmpBead;
                                    distList[neigh_num]  = Dist_Vec3n(r_disp);
                                    neigh_num++;
                                }
                        }
                }
        }
    return neigh_num;
}

/// LatPos_copy: Copies the position elements to copy_vec.
/// \param copy_vec
/// \param in_vec
void LatPos_copy(int* copy_vec, const int* in_vec)
{
    copy_vec[0] = in_vec[0];
    copy_vec[1] = in_vec[1];
    copy_vec[2] = in_vec[2];
}

/// LatPos_gen_rand_wRad. Given a radius, we generate a random set of coordinates between -Radius and Radius
/// for each dimension. The result is stored in outVec
/// \param outVec
/// \param nRadius
void LatPos_gen_rand_wRad(int* outVec, const int nRadius)
{
    int j;
    int radUp = 2 * nRadius + 1;
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = rand() % radUp;
        }
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] -= nRadius;
        }
}

/// LatPos_add_wPBC: Given firVec and secVec, we add the two vectors, and take care of PBCs, storing in outVec.
/// The implementation assumes that the sum of any coordinate can not be larger than twice the lattice-size,
/// because there is no explicit modulo operation.
/// \param outVec
/// \param firVec
/// \param secVec
void LatPos_add_wPBC(int* outVec, const int* firVec, const int* secVec)
{
    short j;
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = firVec[j] + secVec[j];
        }
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = outVec[j] < 0 ? outVec[j] + nBoxSize[j] : outVec[j];
        }
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = outVec[j] >= nBoxSize[j] ? outVec[j] - nBoxSize[j] : outVec[j];
        }
}

/// LatPos_add_wPBC_ofComp: Given firVec and secVec, we add the compNum component of the two vectors,
/// and take care of PBCs, storing in outVec.
/// The implementation assumes that the sum of any coordinate can not be larger than twice the lattice-size,
/// because there is no explicit modulo operation.
/// \param outVec
/// \param firVec
/// \param secVec
inline void LatPos_add_wPBC_ofComp(int* outVec, const int* firVec, const int* secVec, const int compNum)
{
    outVec[compNum] = firVec[compNum] + secVec[compNum];
    outVec[compNum] = outVec[compNum] < 0 ? outVec[compNum] + nBoxSize[compNum] : outVec[compNum];
    outVec[compNum] = outVec[compNum] >= nBoxSize[compNum] ? outVec[compNum] - nBoxSize[compNum] : outVec[compNum];
}

/// LatPos_add_noPBC: Given the two arrays firVec and secVec, we store the sum of the two arrays in outVec.
/// \param outVec
/// \param firVec
/// \param secVec
void LatPos_add_noPBC(int* outVec, const int* firVec, const int* secVec)
{
    short j;
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = firVec[j] + secVec[j];
        }
}

/// BeadPos_sub_wPBC: Given firVec and secVec, we subtract second from first, and take care of PBCs, storing in outVec.
/// The implementation assumes that the sum of any coordinate can not be larger than the lattice size
/// because there is no explicit modulo operation.
/// This version is specifically for bead-positions, or for inter-bead calculations.
/// \param outVec
/// \param firVec
/// \param secVec
void BeadPos_sub_wPBC(int* outVec, const int* firVec, const int* secVec)
{
    short j;
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = firVec[j] - secVec[j];
        }
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = outVec[j] < -nBoxSize[j] / 2 ? outVec[j] + nBoxSize[j] / 2 : outVec[j];
        }
    for (j = 0; j < POS_MAX; j++)
        {
            outVec[j] = outVec[j] >= nBoxSize[j] / 2 ? outVec[j] - nBoxSize[j] / 2 : outVec[j];
        }
}

/// OP_GetTopoBonds - including beadID, get all the beads that are covalently bonded to beadID.
/// \param beadID
/// \param dum_list
/// \return Number_of_bonds+1. Can be used to loop over list which includes beadID.
int OP_GetTopoBonds(const int beadID, int* dum_list)
{
    int top_it = 0;
    int curID  = beadID;

    while (curID != -1)
        {
            dum_list[top_it] = curID;
            curID            = topo_info[beadID][top_it++];
        }

    return top_it;
}

int Vec3n_DotProd(const int* vec_1, const int* vec_2)
{
    return vec_1[0] * vec_2[0] + vec_1[1] * vec_2[1] + vec_1[2] * vec_2[2];
}

/// Vec3n_CosTheta - calculates the Cos(theta) between the vectors v1 and v2.
/// Cos(theta) = v_1 . v_2 / |v_1| / |v_2|
/// \param v1
/// \param v2
/// \return
float Vec3n_CosTheta(const int* v1, const int* v2)
{
    const int dotProd       = Vec3n_DotProd(v1, v2);
    const int vec1MagSq     = Vec3n_DotProd(v1, v1);
    const int vec2MagSq     = Vec3n_DotProd(v2, v2);
    const float dumDenom    = (float) vec1MagSq * (float) vec2MagSq;
    const float dumCosTheta = (float) dotProd / sqrtf(dumDenom);

    return dumCosTheta;
}

float Vec3n_AngleBetweenVecs(const int* vec_1, const int* vec_2)
{

    const float dumAng = Vec3n_CosTheta(vec_1, vec_2);
    return acosf(dumAng);
}

float BeadPos_CosThetaOfBeads(const int bead1, const int bead2)
{

    const int r_pos1[POS_MAX] = {bead_info[bead1][0], bead_info[bead1][1], bead_info[bead1][2]};

    const int r_pos2[POS_MAX] = {bead_info[bead2][0], bead_info[bead2][1], bead_info[bead2][2]};

    return Vec3n_CosTheta(r_pos1, r_pos2);
}
