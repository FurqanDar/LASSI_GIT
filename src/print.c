#include "print.h"
#include "cluster.h"
#include "energy.h"
#include "global.h"
#include "structure.h"

/// PrintToScreen_EnergyMatrix - fancy function to print symmetric matrices with labels
/// \param strTitle
/// \param nSeqEn
/// \param fArray
/// \param param
void PrintToScreen_EnergyMatrix(char* strTitle, int nSeqEn, float fArray[MAX_AA][MAX_AA][MAX_E], int param)
{
    const int nLen   = nSeqEn;
    const int outLen = (nLen + 3) * 5;
    int i, j;

    char sSectionHead[512];
    memset(sSectionHead, '-', outLen);
    sSectionHead[outLen] = NULL;

    printf("%s\n", sSectionHead);

    printf("| %-5s ", strTitle);
    for (i = 0; i < nLen; i++)
        {
            printf("%4d  ", i);
        }
    printf("%3s\n", "|");

    for (i = 0; i < nLen; i++)
        {
            printf("| %-5d ", i);
            for (j = 0; j < nLen; j++)
                {
                    printf("%5.2f ", fArray[i][j][param]);
                }
            printf("%3s\n", "|");
        }

    printf("%s\n", sSectionHead);
}

/// TrajArr_Index
/// \param beadID
/// \param nFrameNumber
/// \param beadProp
/// \return
long TrajArr_Index(const int beadID, const int nFrameNumber, const int beadProp)
{
    return beadProp + BEADINFO_MAX * (beadID + tot_beads * nFrameNumber);
}

/// Write_ClusterDist - write the cluster histogram to a separate file. NOT USED
/// ANYMORE, but might find a use later Always appends to the file for this run.
/// Stopped using it because the IO load was slowing things down at our cluster
/// Would be a good way to gather proper statistics on clusters as the runs went
/// on \param filename \param nGen
void Write_ClusterDist(char* filename, long nGen)
{
    FILE* fp;
    int i;
    if (nGen == -1)
        {
            fp = fopen(filename, "w"); // overwrite
        }
    else
        {
            fp = fopen(filename, "a");
        }

    if (nGen == -1)
        { // title
            fprintf(fp, "#Step followed by histogram\n");
        }
    else
        {
            fprintf(fp, "#%ld\n", nGen);
            for (i = 0; i <= tot_chains; i++)
                {
                    fprintf(fp, "%ld\t", naClusHistList[i]);
                }
            fprintf(fp, "\n");
        }

    fclose(fp);
}

/// Write_GyrTen - write the total gyration tensor to a file. Remember that we
/// only have 6 things in a symmetric 3x3 tensor. Always appends to the file for
/// this run. Stopped using it because the IO load was slowing things down at
/// our cluster Would be a good way to gather proper statistics on the shapes
/// and such as the runs went on \param filename \param nGen
void Write_GyrTen(char* filename, long nGen)
{
    FILE* fp;
    int i;
    if (nGen == -1)
        {
            fp = fopen(filename, "w"); // overwrite
        }
    else
        {
            fp = fopen(filename, "a");
        }

    if (nGen == -1)
        { // title
            fprintf(fp, "#Step followed by Gyration Tensor Vals: xx yy zz xy yz xz\n");
        }
    else
        {
            fprintf(fp, "#%ld\n", nGen);
            for (i = 0; i < 7; i++)
                {
                    if (i == 5)
                        {
                            continue;
                        } // Not smart enough to make a better indexing function
                    fprintf(fp, "%f\t", fGyrTensor[i]);
                }
            fprintf(fp, "\n");
        }

    fclose(fp);
}

/// FileIO_Write_MCMoveHeader - write the header for the MC move acceptance file.
/// \param fileName
void FileIO_Write_MCMoveHeader(const char* fileName)
{
    FILE* fp = fopen(fileName, "a");
    fprintf(fp, "#Steps, and Temp, are followed by (rej, acc) for each MC Move.\n");
    fprintf(fp, "#(nStep, T)             "
                "STROT                   "
                "LOCAL                   "
                "COLOCAL                 "
                "MTLOCAL                 "
                "SNAKE                   "
                "TRANS                   "
                "SM_CLST                 "
                "CLST                    "
                "PIVOT                   "
                "BRROT                   "
                "DBPVT                   "
                "PR_CLST                 "
                "\n");
    fclose(fp);
}

/// FileIO_WriteTo_MCMoveFile - writes the acceptance and rejection ratios of the various moves. Just keeps appending to
/// that file. Can be used to track the 'dynamics' during a simulation. ALSO: Zeros out all the acceptance arrays.
/// \param filename
/// \param nGen
/// \param fMCTemp
void FileIO_WriteTo_MCMoveFile(const char* filename, const long nGen, float const fMCTemp)
{
    FILE* fp = fopen(filename, "a+");
    int i;                                           // Iterator
    fprintf(fp, "%-10ld  %-10.2f  ", nGen, fMCTemp); // Step and Temp
    for (i = 1; i < MAX_MV; i++)
        {
            fprintf(fp, "%-10ld  %-10ld  ", naMCAccepMat[0][i], naMCAccepMat[1][i]);
        }
    fprintf(fp, "\n");

    for (i = 1; i < MAX_MV; i++)
        {
            naMCAccepMat[0][i] = 0;
            naMCAccepMat[1][i] = 0;
            // This way the print function will zero out the matrix every time
            // we print to a file!
        }

    fclose(fp);
}

/// Print_LogToScreen - prints current status of the run. Overall move
/// acceptance ratios, and energies. \param nGen \param run_it
void Print_LogToScreen(long nGen, const int run_it)
{
    ScreenIO_Print_Log_FullRun(nGen, run_it);
}

/// FileIO_Write_EnergyHeader - write the header for the energy file.
/// \param fileName
void FileIO_Write_EnergyHeader(const char* fileName)
{
    FILE* fp = fopen(fileName, "a+");
    fprintf(fp, "#           "
                "STEP        "
                "TOT         "
                "OVLP        "
                "CONT        "
                "SC-SC       "
                "FSOL        "
                "T_IND       "
                "STIFF       "
                "\n");
    fclose(fp);
}

/// FileIO_AppendEnergyTo_EnergyFile - write the decomposed energy of the system to a file.
/// Just appends to the file for this run.
/// \param fileNameStr
/// \param nGen
void FileIO_AppendEnergyTo_EnergyFile(const char* fileNameStr, const long nGen)
{
    FILE* fp = fopen(fileNameStr, "a");

    int i;
    fprintf(fp, "%-10ld", nGen);
    for (i = 0; i < (MAX_E); i++)
        {
            fprintf(fp, "  %-10.2e", faCurrEn[i]);
        }
    fprintf(fp, "\n");

    fclose(fp);
}

/// FileIO_HandleTrajectory -
/// \param fileNameStr
/// \param run_it
/// \param nGen
void FileIO_HandleTrajectory(const char* fileNameStr, const int run_it, const long nGen)
{
    if (nTrajMode == 1)
        {
            Save_Trajectory(nGen, nTrajCurFrame);
            nTrajCurFrame++;
        }
    else
        {
            if (run_it >= 0)
                {
                    FileIO_AppendTrajFrame_ToFile(fileNameStr, nGen + nMCStepsForTherm + run_it * nMCStepsPerCycle);
                }
            else
                {
                    FileIO_AppendTrajFrame_ToFile(fileNameStr, nGen);
                }
        }
}

/// FileIO_AppendTrajFrame_ToFile -  write a LAMMPS style trajectory file for the system.
/// Appends to the file for this run.
/// \param filename
/// \param nGen
void FileIO_AppendTrajFrame_ToFile(const char* filename, const long nGen)
{
    // Writes the trajectory in LAMMPS format. To be viewed with VMD (hopefully). Read the LAMMPS dump documentation for
    // the actual format of the file
    FILE* fp = fopen(filename, "a"); // We will append to the file made for this run;

    int i; // Looping index
    fprintf(fp, "ITEM: TIMESTEP\n");
    fprintf(fp, "%ld\n", nGen); // Timestep

    fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fp, "%ld\n", tot_beads); // Total atom number

    fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");                               // BCs are always periodic for now
    fprintf(fp, "0 %d\n0 %d\n0 %d\n", nBoxSize[0], nBoxSize[1], nBoxSize[2]); // Box dimensions

    fprintf(fp, "ITEM: ATOMS id type mol x y z bP\n"); // What we are printing

    for (i = 0; i < tot_beads; i++)
        {
            fprintf(fp, "%d %d %d %d %d %d %d\n", i, bead_info[i][BEAD_TYPE], bead_info[i][BEAD_CHAINID],
                    bead_info[i][POS_X], bead_info[i][POS_Y], bead_info[i][POS_Z], bead_info[i][BEAD_FACE]);
        }

    fclose(fp);
}

/// Save_Trajectory -  saves the current position onto the total array
/// \param nGen
/// \param curFrame
void Save_Trajectory(const long nGen, const long curFrame)
{

    int i, j;
    long ar_idx;
    for (i = 0; i < tot_beads; i++)
        {
            for (j = 0; j < BEADINFO_MAX; j++)
                {
                    ar_idx                = TrajArr_Index(i, curFrame, j);
                    n_TOTTRAJ_Arr[ar_idx] = bead_info[i][j];
                }
        }
}

/// Write_Saved_Trajectory -  write a LAMMPS style trajectory file for the
/// system. This funcion writes all the stored frames for this temperature cycle
/// \param filename
/// \param nGen
void Write_Saved_Trajectory(char* filename, const int run_it)
{
    // Writes the trajectory in LAMMPS format. To be viewed with VMD
    // (Hopefully). Read the LAMMPS dump documentation for the actual formate of
    // the file
    FILE* fp;
    fp = fopen(filename, "a"); // We will append to the file made for this run

    int i, j, k; // Looping index
    int ar_idx;

    for (i = 0; i < nTrajCurFrame; i++)
        {
            fprintf(fp, "ITEM: TIMESTEP\n");
            if (run_it >= 0)
                {
                    fprintf(fp, "%ld\n",
                            nMCStepsForTherm + run_it * nMCStepsPerCycle + i * nReport[REPORT_CONFIG]); // Timestep
                }
            else
                {
                    fprintf(fp, "%ld\n", i * nReport[REPORT_CONFIG]); // Timestep
                }

            fprintf(fp, "ITEM: NUMBER OF ATOMS\n");
            fprintf(fp, "%ld\n", tot_beads); // Total atom number

            fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");                               // BCs are always periodic for now
            fprintf(fp, "0 %d\n0 %d\n0 %d\n", nBoxSize[0], nBoxSize[1], nBoxSize[2]); // Box dimensions

            fprintf(fp, "ITEM: ATOMS id type mol x y z bP\n"); // What we are printing

            for (j = 0; j < tot_beads; j++)
                {
                    fprintf(fp, "%d", j);
                    ar_idx = TrajArr_Index(j, i, BEAD_TYPE);
                    fprintf(fp, " %d", n_TOTTRAJ_Arr[ar_idx]);
                    ar_idx = TrajArr_Index(j, i, BEAD_CHAINID);
                    fprintf(fp, " %d", n_TOTTRAJ_Arr[ar_idx]);
                    ar_idx = TrajArr_Index(j, i, POS_X);
                    fprintf(fp, " %d", n_TOTTRAJ_Arr[ar_idx]);
                    ar_idx = TrajArr_Index(j, i, POS_Y);
                    fprintf(fp, " %d", n_TOTTRAJ_Arr[ar_idx]);
                    ar_idx = TrajArr_Index(j, i, POS_Z);
                    fprintf(fp, " %d", n_TOTTRAJ_Arr[ar_idx]);
                    ar_idx = TrajArr_Index(j, i, BEAD_FACE);
                    fprintf(fp, " %d", n_TOTTRAJ_Arr[ar_idx]);
                    fprintf(fp, "\n");
                }
        }

    fclose(fp);
}

/// PrintToScreen_AllEnergyMatrices
void PrintToScreen_AllEnergyMatrices(void)
{

    const char lBrace[] = "<======      ";
    const char rBrace[] = "      ======>";

    printf("%s Energy Matrices %s\n", lBrace, rBrace);

    PrintToScreen_EnergyMatrix("OVLP ", nBeadTypes, fEnergy, E_OVLP);

    PrintToScreen_EnergyMatrix("CONT ", nBeadTypes, fEnergy, E_CONT);

    PrintToScreen_EnergyMatrix("SC_SC", nBeadTypes, fEnergy, E_SC_SC);

    PrintToScreen_EnergyMatrix("FSOL ", nBeadTypes, fEnergy, E_F_SOL);

    PrintToScreen_EnergyMatrix("STIFF", nBeadTypes, fEnergy, E_STIFF);
}

/// PrintToScreen_MCMoveFreqs
void PrintToScreen_MCMoveFreqs(void)
{
    char* MoveName[MAX_MV];
    MoveName[MV_PIVOT]      = "Pivot           ";
    MoveName[MV_DBPVT]      = "Double Pivot    ";
    MoveName[MV_CLSTR]      = "Larger Cluster  ";
    MoveName[MV_SMCLSTR]    = "Smaller Cluster ";
    MoveName[MV_STROT]      = "Face Change     ";
    MoveName[MV_LOCAL]      = "Local           ";
    MoveName[MV_COLOCAL]    = "Co-local        ";
    MoveName[MV_MTLOCAL]    = "Shake           ";
    MoveName[MV_BRROT]      = "Rotate Branched ";
    MoveName[MV_SNAKE]      = "Slithering Snake";
    MoveName[MV_TRANS]      = "Translation     ";
    MoveName[MV_PR_SMCLSTR] = "Pr. Smal Cluster";

    float freqMin = 1e10f;
    int i;
    for (i = MV_NULL + 1; i < MAX_MV; i++)
        {
            if (freqMin >= fMCFreq[i] && fMCFreq[i] != 0.)
                {
                    freqMin = fMCFreq[i];
                }
        }

    for (i = 0; i < 30; i++)
        {
            printf("-");
        }
    printf("\n");
    printf("| MC Move Frequencies%9s\n", "|");
    for (i = MV_NULL + 1; i < MAX_MV; i++)
        {
            printf("| %-16s %5d %5s\n", MoveName[i], (int) ceilf(fMCFreq[i] / freqMin), "|");
        }
    for (i = 0; i < 30; i++)
        {
            printf("-");
        }
}

/// ScreenIO_Print_KeyFile - print the keyfile that was read in to the screen
void ScreenIO_Print_KeyFile(void)
{ // should be output-dependent (stdout, stderr, other files)

    int i;
    const char lBrace[] = "<======      ";
    const char rBrace[] = "      ======>";
    printf("%s System Settings %s\n", lBrace, rBrace);
    printf("Number of Bead Types = %d\n", nBeadTypes);
    printf("Number of Beads      = %ld\n", tot_beads);
    printf("Number of Chains     = %ld\n", tot_chains);
    printf("Number of Components = %ld\n", tot_chain_types);
    printf("Box Sizes            = %3d %3d %3d\n", nBoxSize[0], nBoxSize[1], nBoxSize[2]);
    printf("Monomer density      = %1.1e\n",
           (float) tot_beads / (float) nBoxSize[0] / (float) nBoxSize[1] / (float) nBoxSize[2]);
    printf("\n");

    PrintToScreen_AllEnergyMatrices();

    printf("\n");

    printf("%s Linker Info %s\n", lBrace, rBrace);
    printf("Linker length             = %.2f\n", fLinkerLength);
    printf("Linker spring constant    = %.2f\n", fLinkerSprCon);
    printf("Linker equilibrium length = %.2f\n", fLinkerEqLen);
    printf("\n");

    printf("%s MC Setup %s\n", lBrace, rBrace);
    printf("Temperature Inverted           = %d\n", nTemp_inv);
    printf("MC Temperatures: (First, Last) = (%.2f, %.2f)\n", fKT, fKT + (float) (nTot_CycleNum - 1) * fdelta_temp);
    printf("Temperature Mode               = %d\n", nAnnealing_Mode);
    printf("Indent Mode                    = %d\n", nInitialPotential_Mode);
    printf("Rotational Bias Mode           = %d\n", RotBias_Mode);
    printf("Number of MC Cycles            = %d\n", nTot_CycleNum);
    printf("Number of MC Steps/Cycle       = %e\n", (float) nMCStepsPerCycle);
    printf("Thermalizing Temperature       = %.2f\n", fPreKT);
    printf("Number of Thermalizing Steps   = %e\n", (float) nMCStepsForTherm);
    printf("RNG Seed                       = %d\n", nRNG_Seed);
    printf("Clustering Mode                = %d\n", nClusteringMode);

    PrintToScreen_MCMoveFreqs();
    printf("\n");
}

/// ScreenIO_Print_SystemEnergy
void ScreenIO_Print_SystemEnergy(void)
{

    int i;

    char sSectionHead[32];
    memset(sSectionHead, '-', 17);
    sSectionHead[17] = NULL;

    printf("%s\n", sSectionHead);
    printf("Energies\n");
    printf("Tot  : %-10.2e |\n", faCurrEn[E_TOT]);
    printf("Ovlp : %-10.2e |\n", faCurrEn[E_OVLP]);
    printf("Cont : %-10.2e |\n", faCurrEn[E_CONT]);
    printf("Aniso: %-10.2e |\n", faCurrEn[E_SC_SC]);
    printf("FSol : %-10.2e |\n", faCurrEn[E_F_SOL]);
    printf("Stiff: %-10.2e |\n", faCurrEn[E_STIFF]);
    printf("%s\n", sSectionHead);
}

/// void ScreenIO_Print_AcceptanceRatios
void ScreenIO_Print_AcceptanceRatios(void)
{

    char* MoveName[MAX_MV];
    MoveName[MV_PIVOT]      = "Pivot       ";
    MoveName[MV_DBPVT]      = "Double Pivot";
    MoveName[MV_CLSTR]      = "La Cluster  ";
    MoveName[MV_SMCLSTR]    = "Sm Cluster  ";
    MoveName[MV_STROT]      = "Face Change ";
    MoveName[MV_LOCAL]      = "Local       ";
    MoveName[MV_COLOCAL]    = "Co-local    ";
    MoveName[MV_MTLOCAL]    = "Multi Local ";
    MoveName[MV_BRROT]      = "Rot. Br.    ";
    MoveName[MV_SNAKE]      = "Sli. Snake  ";
    MoveName[MV_TRANS]      = "Translation ";
    MoveName[MV_PR_SMCLSTR] = "Pr. Sm. Cls.";

    int i, j;
    float fAccRatio = 0.f;
    lLong nMoveSum  = 0;

    char sSectionHead[32];
    memset(sSectionHead, '-', 22);
    sSectionHead[22] = NULL;

    printf("%s\n", sSectionHead);
    printf("Acceptance Ratios:\n");
    for (i = 1; i < MAX_MV; i++)
        {
            nMoveSum = naMCAccepMat[0][i] + naMCAccepMat[1][i];
            if (nMoveSum)
                {
                    fAccRatio = 100.f * (float) naMCAccepMat[1][i] / (float) nMoveSum;
                    printf("%-12s: %-7.2f|\n", MoveName[i], fAccRatio);
                }
            else
                {
                    printf("%-12s: %-7s|\n", MoveName[i], "NA");
                }
        }
    printf("%s\n", sSectionHead);
}

/// ScreenIO_Print_Log_Thermalization - print the log to the screen.
void ScreenIO_Print_Log_Thermalization(const long nGen)
{
    printf("Run Cycle: Thermalization\n");
    printf("Step     : %8.3e\n", (float) nGen);
    printf("MC Temp  : %8.3e\n", fCuTemp);

    if (nTotClusCounter > 0)
        {
            printf("Perc Phi : %8.3f\n",
                   ((float) nLargestClusterRightNow) / ((float) nTotClusCounter) / (float) tot_chains);
        }

    ScreenIO_Print_SystemEnergy();
    ScreenIO_Print_AcceptanceRatios();
}

/// ScreenIO_Print_Log_FullRun - print the log to the screen.
void ScreenIO_Print_Log_FullRun(const long nGen, const int run_cycle)
{
    printf("Run Cycle: %d\n", run_cycle);
    printf("Step     : %8.3e\n", (float) nGen);
    printf("MC Temp  : %8.3e\n", fCuTemp);

    if (nTotClusCounter > 0)
        {
            printf("Perc Phi : %8.3f\n",
                   ((float) nLargestClusterRightNow) / ((float) nTotClusCounter) / (float) tot_chains);
        }

    ScreenIO_Print_SystemEnergy();
    ScreenIO_Print_AcceptanceRatios();
}

/// Write_RDF_ComponentWise - old implementation of printing the RDF, component
/// by component. Always appends to the file for this run. Stopped using it
/// because the IO load was slowing things down at our cluster Would be a good
/// way to gather proper statistics on clusters as the runs went on. TODO:
/// update this to work with the new indexing
/// \param filename
/// \param nGen
void Write_RDF_ComponentWise(char* filename, long nGen)
{
    FILE* fp;
    if (nGen == -1)
        {
            fp = fopen(filename, "w"); // overwrite
        }
    else
        {
            fp = fopen(filename, "a");
        }
    int i, k;
    if (nGen == -1)
        { // title
            fprintf(fp, "#Row-by-row split RDF. dr = 1/4\n");
        }
    else
        {
            for (k = 0; k < nRDF_TotComps; k++)
                {
                    for (i = 0; i < nRDF_TotBins; i++)
                        {
                            fprintf(fp, "%.5Lf\t", ldRDF_Arr[RDFArr_Index(0, k, i)]);
                        }
                    fprintf(fp, "\n");
                }
        }

    fclose(fp);
}

/// FileIO_WriteTo_TopFile - write a LAMMPS DATA file, which can be read in VMD to store
/// topological information for the system Only writes the topology once at the
/// beginning of the runs. \param filename
void FileIO_WriteTo_TopFile(const char* filename)
{
    /*
    Writes a topology file for VMD. Since the trajectory is saved in the LAMMPS
    format, the topology file is also in the LAMMPS format. The format for this
    file is one which is used as input data for LAMMPS. This function copies the
    format given in read_data for LAMMPS:
    https://lammps.sandia.gov/doc/read_data.html
    */
    FILE* fp;
    fp = fopen(filename, "w"); // Just overwrite over if last file exists.
    int i, j, k;               // Loop iterators.
    int numBonds;              // Used to count total number of bonds!
    printf("Writing the topology file!\n");
    fprintf(fp, "LAMMPS Description\n");      // The file must start with this.
    fprintf(fp, "\n");                        // Empty line.
    fprintf(fp, "\t%ld\tatoms\n", tot_beads); // Listing total number of atoms
    numBonds = 0;
    for (i = 0; i < tot_beads; i++)
        {
            for (j = 0; j < MAX_BONDS; j++)
                {
                    if (topo_info[i][j] != -1)
                        { // There is a bond between i and topo_info[i][j]
                            numBonds++;
                        }
                }
        }

    fprintf(fp, "\t%d\tbonds\n", numBonds);      // Listing number of bonds
    fprintf(fp, "\t0\tangles\n");                // Systems don't have angle depenece yet
    fprintf(fp, "\t0\tdihedrals\n");             // Systems don't have dihedrals  yet
    fprintf(fp, "\t0\timpropers\n");             // Systems don't have imporopers yet
    fprintf(fp, "\n");                           // Empty line.
    fprintf(fp, "%d\tatom types\n", nBeadTypes); // This many bead-types
    fprintf(fp, "1\tbond types\n");              // System can have multiple bond-lengths
    fprintf(fp, "0\tangle types\n");             // Systems don't have any angular forces yet
    fprintf(fp, "\n");                           // Empty line.
    fprintf(fp, " 0 %d xlo xhi\n", nBoxSize[0]);
    fprintf(fp, " 0 %d ylo yhi\n", nBoxSize[1]);
    fprintf(fp, " 0 %d zlo zhi\n", nBoxSize[2]);
    fprintf(fp, "\n");       // Empty line.
    fprintf(fp, "Masses\n"); // These don't really mean anything
    fprintf(fp, "\n");       // Empty line.
    for (i = 0; i < nBeadTypes; i++)
        {
            fprintf(fp, "%d\t1.0\n", i);
        }

    fprintf(fp, "\n");             // Empty line.
    fprintf(fp, "Atoms # bond\n"); // Signifying the beginning of atom coordinates.
    fprintf(fp, "\n");             // Empty line.
    for (i = 0; i < tot_beads; i++)
        {
            fprintf(fp, "%d %d %d %d %d %d\n", i, bead_info[i][BEAD_CHAINID], bead_info[i][BEAD_TYPE],
                    bead_info[i][POS_X], bead_info[i][POS_Y], bead_info[i][POS_Z]);
        }                   // Done with the coordinates
    fprintf(fp, "\n");      // Empty line.
    fprintf(fp, "Bonds\n"); // Signifying the beginning of atom coordinates.
    fprintf(fp, "\n");      // Empty line.
    k = 0;                  // This guy counts bondIDs
    for (i = 0; i < tot_beads; i++)
        {
            for (j = 0; j < MAX_BONDS; j++)
                {
                    if (topo_info[i][j] != -1)
                        { // There is a bond between i and topo_info[i][j]
                            fprintf(fp, "%d %d %d %d\n", k, linker_len[i][j] / (int) fLinkerLength, i, topo_info[i][j]);
                            k++;
                        }
                }
        }              // Done with the coordinates
    fprintf(fp, "\n"); // Empty line.

    fclose(fp);
}

/// Write_SysProp - writes the RDF, CLUS and GyrTen avg rrays to a file. OLD implementation that should not be used yet
/// Also, TODO: update this to work with the new indexing
/// \param filename
void Write_SysProp(char* filename)
{
    FILE* fp;
    int i, j;
    fp = fopen(filename, "w");
    fprintf(fp, "#Contains various averaged quatities.\n");
    // Gyration Radius
    fprintf(fp, "#Total Rg\n%f\t%f\n#Cluster Hist\n", fSysGyrRad / (float) nTotGyrRadCounter, (float) nBoxSize[0] / 2.);
    // Cluster Histogram
    fprintf(fp, "%f\t", (float) nLargestClusterRightNow / (float) nTotClusCounter);
    for (i = 1; i <= tot_chains; i++)
        {
            fprintf(fp, "%f\t", (float) naClusHistList[i] / (float) nTotClusCounter);
        }
    // Split RDFs
    fprintf(fp, "\n#Split RDFs. ALL-ALL; DIAGONALS and then from 0 onwards \n");
    for (j = 0; j < nRDF_TotComps; j++)
        {
            for (i = 0; i < nRDF_TotBins; i++)
                {
                    // fprintf(fp, "%LE\t", ldRDF_ARR[j][i] / (float)nRDFCounter);
                    fprintf(fp, "%LE\t", ldRDF_Arr[RDFArr_Index(0, j, i)] / (float) nRDFCounter);
                }
            fprintf(fp, "\n");
        }
    fprintf(fp, "\n#Done");
}

/// FileIO_WriteTo_RDFTotFile - writes the RDFTot file.
/// \param run_it Number of temperature cycles.
void FileIO_WriteTo_RDFTotFile(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_RDF.dat",
            strReportPrefix); // Name Of the RDF Files
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Split RDFs. ALL-ALL; DIAGONALS and then from 0-0, 0-1, and onwards \n");
    for (i = 0; i < run_it; i++)
        {
            fprintf(fp, "#Run_Cycle = %d\n", i);
            for (j = 0; j < nRDF_TotComps; j++)
                {
                    for (k = 0; k < nRDF_TotBins; k++)
                        {
                            fprintf(fp, "%LE\t", ld_TOTRDF_Arr[RDFArr_Index(i, j, k)]);
                        }
                    fprintf(fp, "\n");
                }
        }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_COMDenFile - writes the COMDenTot file.
/// \param run_it Number of temperature cycles.
void FileIO_WriteTo_COMDenFile(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_COMDen.dat",
            strReportPrefix); // Name Of the COM Density Distribution
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Density distribution from the COM outwards. Order is "
                "ChainType.\n");
    for (i = 0; i < run_it; i++)
        {
            fprintf(fp, "#Run_Cycle = %d\n", i);
            for (j = 0; j < nRadDen_TotComps; j++)
                {
                    for (k = 0; k < nRDF_TotBins; k++)
                        {
                            fprintf(fp, "%LE\t", ld_TOTRadDen_Arr[RadDenArr_Index(i, j, k)]);
                        }
                    fprintf(fp, "\n");
                }
        }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_ClusFile - writes the CLUS file.
/// \param run_it Number of temperature cycles.
void FileIO_WriteTo_ClusFile(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_CLUS.dat",
            strReportPrefix); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Cluster Histograms: 1st column is largest cluster, and then clusters of size 1, 2, and so on\n");
    for (i = 0; i < run_it; i++)
        {
            fprintf(fp, "#Run_Cycle = %d\n", i);
            for (k = 0; k <= tot_chains; k++)
                {
                    fprintf(fp, "%LE\t", ld_TOTCLUS_Arr[i][k]);
                }
            fprintf(fp, "\n");
        }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_MolClus - writes the MolClus file.
/// \param run_it Number of temperature cycles.
void FileIO_WriteTo_MolClusFile(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_MolClus.dat",
            strReportPrefix); // Name Of the COM Density Distribution
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Cluster Histograms: 1st column is largest cluster, and then clusters of size 1, 2, and so on. Each "
                "row is a different moltype\n");
    for (i = 0; i < run_it; i++)
        {
            fprintf(fp, "#Run_Cycle = %d\n", i);
            for (j = 0; j < tot_chain_types; j++)
                {
                    for (k = 0; k < tot_chains; k++)
                        {
                            fprintf(fp, "%LE\t", ld_TOTMOLCLUS_Arr[MolClusArr_Index(i, j, k)]);
                        }
                    fprintf(fp, "\n");
                }
        }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_WriteTo_GyrRad - writes the GR file.
/// \param run_it Number of temperature cycles.
void FileIO_WriteTo_GyrRadFile(const int run_it)
{
    FILE* fp;
    int i, j, k;
    char dumFile[512];
    sprintf(dumFile, "%s_GR.dat",
            strReportPrefix); // Name Of the ClusterHistogram Files
    fp = fopen(dumFile, "w");
    fprintf(fp, "#Cluster Histograms: 1st column is largest cluster, and then "
                "clusters of size 1, 2, and so on\n");
    for (i = 0; i < run_it; i++)
        {
            fprintf(fp, "#Run_Cycle = %d\n", i);
            fprintf(fp, "%LE\t", ld_TOTRg_Arr[i][0]);
            fprintf(fp, "%LE\n", ld_TOTRg_Arr[i][1]);
        }
    fprintf(fp, "#Done");
    fclose(fp);
}

/// FileIO_Write_TotalSysProp - writes the RDF, CLUS and GyrTen arrays to their
/// respective files. These arrays store the data over the course of an ENTIRE
/// run, or over all the temp cycles.
/// \param run_it
void FileIO_Write_TotalSysProp(const int run_it)
{
    /* This function writes one large file with all the averaged values from
     * each run_cycle. run_it let's the function know how many cycles to write.
     * As opposed to Write_SysProp, each of the relevant averaged parameters
     * will be in their appropriately named files. The naming convention shall
     * be: $REPORTNAME_*.dat where * will be GR, CLUS and RDF for the different
     * measured quantities. Within each file, each run_cycle shall be written
     * one-by-one in the style of Write_SysProp
     */
    if (nReport[REPORT_RDFTOT] != 0)
        {
            FileIO_WriteTo_RDFTotFile(run_it);
        }

    if (nReport[REPORT_COMDEN] != 0)
        {
            FileIO_WriteTo_COMDenFile(run_it);
        }

    if (nReport[REPORT_NETWORK] != 0)
        {
            FileIO_WriteTo_ClusFile(run_it);
            FileIO_WriteTo_MolClusFile(run_it);
            FileIO_WriteTo_GyrRadFile(run_it);
        }
}

/// FileIO_CreateFile - Creates a new overwritten file with the given name.
/// \param fileName
void FileIO_CreateFile(const char* fileName)
{
    FILE* fp = fopen(fileName, "w+");
    fclose(fp);
}

/// FileIO_CreateRunningDataFiles - Creates the necessary files that are continuously written to over the course
/// of a simulation.
void FileIO_CreateRunningDataFiles(void)
{
    // Trajectory
    if (nReport[REPORT_CONFIG])
        {
            sprintf(fileTraj, "%s_topo.lammpstrj", strReportPrefix); // Name of the topology file
            FileIO_WriteTo_TopFile(fileTraj);                        // Write the topology file. Only need to write once
            if (nTrajMode != 1)
                {
                    sprintf(fileTraj, "%s_trj.lammpstrj", strReportPrefix); // Naming convention for trajectory files.
                    FileIO_CreateFile(fileTraj); // This opens a new trajectory file; each run_it will have its own
                }
        }

    // Energy
    if (nReport[REPORT_ENERGY])
        {
            sprintf(fileEnergy, "%s_energy.dat", strReportPrefix);
            FileIO_CreateFile(fileEnergy); // Open a new energy file; each run_it will have its own
            FileIO_Write_EnergyHeader(fileEnergy);
        }

    // MC Move Acceptance
    if (nReport[REPORT_MCMOVE])
        {
            sprintf(fileMCMove, "%s_mcmove.dat", strReportPrefix);
            FileIO_CreateFile(fileMCMove); // Open a new MCInfo file; each run_it will have its own
            FileIO_Write_MCMoveHeader(fileMCMove);
        }
}

/// ForPrinting_GetReportState. Calculates if the current step is a multiple of the reporting frequency
/// for this particular report.
/// Convenient wrapper around some modulo arithmetic.
/// \param nGen
/// \param thisReport
/// \return
char ForPrinting_GetReportState(const long nGen, const long thisReport)
{
    char dum_log = (char) ((nGen % thisReport) == 0);
    dum_log      = dum_log ? 1 : 0;
    return dum_log;
}

/// DataPrinting_Thermalization - helper function for printing out data during the thermalization cycle.
/// No data analysis is performed during the thermalization procedure.
/// This function decides if, given the MCStep, it is time to print the following:
/// 1. The Log - to the screen.
/// 2. Trajectory - to the file (or saved if in total format)
/// 3. Energy - to the file.
/// 4. MCMove - to the file.
/// \param nGen
void DataPrinting_Thermalization(const long nGen)
{

    char cFlagForEnCal = 0;
    char cLogFlag      = 0;
    char cEnergyFlag   = 0;
    char cAccFlag      = 0;
    char cConfigFlag   = 0;

    if (nReport[REPORT_LOG])
        {
            cLogFlag = ForPrinting_GetReportState(nGen, nReport[REPORT_LOG]);
            if (cLogFlag)
                {
                    // TODO: I think this whole business can be abstracted away as well.
                    if (Check_System_Structure())
                        {
                            fprintf(stderr, "Molecular structure is inconsistent with initial "
                                            "structure.\nCRASHING\n\n");
                            exit(1);
                        }
                    Energy_Total_System();
                    cFlagForEnCal = 1;
                    ScreenIO_Print_Log_Thermalization(nGen);
                }
        }

    if (nGen)
        {
            if (nReport[REPORT_CONFIG])
                {
                    cConfigFlag = ForPrinting_GetReportState(nGen, nReport[REPORT_CONFIG]);
                    if (cConfigFlag)
                        {
                            FileIO_HandleTrajectory(fileTraj, -1, nGen);
                        }
                }

            if (nReport[REPORT_ENERGY])
                {
                    // DO ENERGY SHIT
                    cEnergyFlag = ForPrinting_GetReportState(nGen, nReport[REPORT_ENERGY]);
                    if (cEnergyFlag)
                        {
                            if (! cFlagForEnCal)
                                {
                                    Energy_Total_System();
                                    cFlagForEnCal = 1;
                                }
                            FileIO_AppendEnergyTo_EnergyFile(fileEnergy, nGen);
                        }
                }

            if (nReport[REPORT_MCMOVE])
                {
                    // DO MC_ACC SHIT
                    cAccFlag = ForPrinting_GetReportState(nGen, nReport[REPORT_MCMOVE]);
                    if (cAccFlag)
                        {
                            FileIO_WriteTo_MCMoveFile(fileMCMove, nGen, fCuTemp);
                        }
                }
        }
}

/// DataPrinting_DuringRunCycles
/// \param nGen
/// \param run_it
void DataPrinting_DuringRunCycles(const long nGen, const int run_it)
{
    char cFlagForEnCal = 0;
    char cLogFlag      = 0;
    char cEnergyFlag   = 0;
    char cAccFlag      = 0;
    char cConfigFlag   = 0;

    if (nReport[REPORT_LOG])
        {
            cLogFlag = ForPrinting_GetReportState(nGen, nReport[REPORT_LOG]);
            if (cLogFlag)
                {
                    // TODO: I think this whole business can be abstracted away as well.
                    if (Check_System_Structure())
                        {
                            fprintf(stderr, "Molecular structure is inconsistent with initial "
                                            "structure.\nCRASHING\n\n");
                            exit(1);
                        }
                    Energy_Total_System();
                    cFlagForEnCal = 1;
                    ScreenIO_Print_Log_FullRun(nGen, run_it);
                }
        }

    if (nGen)
        {
            if (nReport[REPORT_CONFIG])
                {
                    cConfigFlag = ForPrinting_GetReportState(nGen, nReport[REPORT_CONFIG]);
                    if (cConfigFlag)
                        {
                            sprintf(fileTraj, "%s_trj.lammpstrj", strReportPrefix);
                            FileIO_HandleTrajectory(fileTraj, run_it, nGen);
                        }
                }

            if (nReport[REPORT_ENERGY])
                {
                    // DO ENERGY SHIT
                    cEnergyFlag = ForPrinting_GetReportState(nGen, nReport[REPORT_ENERGY]);
                    if (cEnergyFlag)
                        {
                            if (! cFlagForEnCal)
                                {
                                    Energy_Total_System();
                                    cFlagForEnCal = 1;
                                }
                            FileIO_AppendEnergyTo_EnergyFile(fileEnergy, nGen);
                        }
                }

            if (nReport[REPORT_MCMOVE])
                {
                    // DO MC_ACC SHIT
                    cAccFlag = ForPrinting_GetReportState(nGen, nReport[REPORT_MCMOVE]);
                    if (cAccFlag)
                        {
                            FileIO_WriteTo_MCMoveFile(fileMCMove, nGen, fCuTemp);
                        }
                }
        }
}

/// DataAnalysis_DuringRunCycles
/// \param nGen
/// \param run_it
void DataAnalysis_DuringRunCycles(const long nGen, const int run_it)
{

    //    char cRDF_flag  = 0;
    //    char cCOM_flag  = 0;
    //    char cCLUS_flag = 0;

    if (nReport[REPORT_RDFTOT])
        { // SysProp is printed outside of this function in main.c, lol
            if (nGen % nReport[REPORT_RDFTOT] == 0)
                {
                    RDF_ComponentWise_Avg();
                }
        }
    if (nReport[REPORT_COMDEN])
        { // SysProp is printed outside of this function in main.c, lol
            if (nGen % nReport[REPORT_COMDEN] == 0)
                {
                    RadDen_Avg_MolTypeWise_FromMolTypeCen();
                }
        }
    if (nReport[REPORT_NETWORK])
        { // SysProp is printed outside of this function in main.c, lol
            if (nGen % nReport[REPORT_NETWORK] == 0)
                {
                    Clus_Perform_Analysis();
                    GyrTensor_GyrRad_Avg();
                }
        }
}

/// FileIO_PreCycle_Init
/// \param run_it
void FileIO_PreCycle_Init(const int run_it)
{
    if (nReport[REPORT_CONFIG])
        {
            if (nTrajMode != 1)
                {
                    sprintf(fileEnergy, "%s_trj.dat", strReportPrefix);
                }
        }

    if (nReport[REPORT_ENERGY] != 0)
        {
            sprintf(fileEnergy, "%s_%d_energy.dat", strReportPrefix, run_it);
            FileIO_CreateFile(fileEnergy); // Open a new Energy file; each run_it will have its own
        }
    if (nReport[REPORT_MCMOVE] != 0)
        {
            sprintf(fileMCMove, "%s_%d_mcmove.dat", strReportPrefix, run_it);
            FileIO_CreateFile(fileMCMove); // Open a new MCInfo file; each run_it will have its own
        }
}
/// FileIO_WriteRestart_ForRun
/// \param run_it
void FileIO_WriteRestart_ForRun(const int run_it)
{
    sprintf(fileTraj, "%s_%d_restart.lammpstrj", strReportPrefix,
            run_it); // Naming convention for trajectory files.
    FileIO_CreateFile(fileTraj);
    FileIO_AppendTrajFrame_ToFile(fileTraj, nMCStepsForTherm + (run_it + 1) * nMCStepsPerCycle);
}

/// FileIO_WriteRestart_ForThermalization
void FileIO_WriteRestart_ForThermalization(void)
{
    if (nMCStepsForTherm)
        {
            sprintf(fileTraj, "%s_EQ_restart.lammpstrj", strReportPrefix); // Naming convention for trajectory files.
            FileIO_CreateFile(fileTraj);
            FileIO_AppendTrajFrame_ToFile(fileTraj, nMCStepsForTherm);
        }
}

/// CopyData_All - copies data from run_it specific data arrays to the overall
/// global data arrays that are printed later. Also note that the averaging, or
/// dividing by the frequency of acquisitions, occurs here.
/// \param run_it
void CopyData_All(const int run_it)
{
    if (nReport[REPORT_RDFTOT])
        {
            CopyData_RDF(run_it);
        }
    if (nReport[REPORT_COMDEN])
        {
            CopyData_COMDen(run_it);
        }
    if (nReport[REPORT_NETWORK])
        {
            CopyData_Clus(run_it);
        }
}

/// CopyData_RDF: Copies ldRDF_Arr into ld_TOTRDF_Arr
/// \param run_it: which run cycle we are on. Should be >= 0.
void CopyData_RDF(const int run_it)
{
    int i, j;
    for (i = 0; i < nRDF_TotComps; i++)
        {
            for (j = 0; j < nRDF_TotBins; j++)
                {
                    ld_TOTRDF_Arr[RDFArr_Index(run_it, i, j)] =
                        ldRDF_Arr[RDFArr_Index(0, i, j)] / (long double) nRDFCounter;
                }
        }
}

/// CopyData_RDF: Copies ldRadDen_Arr into ld_TOTRadDen_Arr
/// \param run_it: which run cycle we are on. Should be >= 0.
void CopyData_COMDen(const int run_it)
{
    int i, j;
    for (i = 0; i < nRadDen_TotComps; i++)
        {
            for (j = 0; j < nRDF_TotBins; j++)
                {
                    ld_TOTRadDen_Arr[RadDenArr_Index(run_it, i, j)] =
                        ldRadDen_Arr[RadDenArr_Index(0, i, j)] / (long double) nRadDenCounter;
                }
        }
}

/// CopyData_Clus - Copies ldMOLClUS_Arr into ld_TOTMOLCLUS_Arr. Also copies naClusHistList into ld_TOTCLUS_Arr.
/// And stores the Gyration radius.
/// \param run_it: which run cycle we are on. Should be >= 0.
void CopyData_Clus(const int run_it)
{
    int i, j;
    ld_TOTCLUS_Arr[run_it][0] = (long double) nLargestClusterRightNow / (long double) nTotClusCounter;
    for (i = 1; i <= tot_chains; i++)
        {
            ld_TOTCLUS_Arr[run_it][i] = (long double) naClusHistList[i] / (long double) nTotClusCounter;
        }
    for (i = 0; i < tot_chains; i++)
        {
            for (j = 0; j < tot_chain_types; j++)
                {
                    ld_TOTMOLCLUS_Arr[MolClusArr_Index(run_it, j, i)] =
                        ldMOLCLUS_Arr[MolClusArr_Index(0, j, i)] / (long double) nTotClusCounter;
                }
        }
    ld_TOTRg_Arr[run_it][0] = (long double) fSysGyrRad / (long double) nTotGyrRadCounter;
    ld_TOTRg_Arr[run_it][1] = (long double) nBoxSize[0] / 2.;
}
