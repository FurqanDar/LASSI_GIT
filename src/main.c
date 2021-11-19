#include "cluster.h"
#include "energy.h"
#include "global.h"
#include "initialize.h"
#include "mcmove.h"
#include "parsekey.h"
#include "print.h"
#include "structure.h"

int main(int argc, char* argv[])
{
    // Read in the system commands
    char keyfile[100];
    if (argc == 3 && strcmp(argv[1], "-k") == 0)
        {
            strcpy(keyfile, argv[2]);
        }
    else
        {
            strcpy(keyfile, "param.key");
        }
    // Text formatting helpers.
    char lBrace[] = "<======      ";
    char rBrace[] = "      ======>";
    // Read the param file!
    int errorcode;
    errorcode = Parse_Keyfile(keyfile);
    if (errorcode == 0)
        {
            printf("Key file %s was successfully parsed.\n\n", keyfile);
            srand(nRNG_Seed);
            ScreenIO_Print_KeyFile();
        }
    else
        {
            printf("ERROR: unable to parse key file %s.\n", keyfile);
            exit(1);
        }

    printf("%s Initializing %s\n", lBrace, rBrace);

    clock_t tStart = clock();

    // Allocating memory and initializing all global arrays.
    Memory_Initialization_AtStart();
    Global_Array_Initialization_AtStart();

    // Reading in the structure file, and figuring out initial conditions.
    if (bReadConf == 0)
        {
            Initial_Conditions_Simple();
        }
    else if (bReadConf == 1)
        {
            Initial_Conditions_FromFile();
            if (nMCStepsForTherm > 0)
                { // Remember that all thermalizing variants of
                  // MC moves cannot handle bonds.
                    Initial_Conditions_BreakBonds();
                }
        }
    else
        {
            printf("Wrong initial conditions condition given. Crashing!\n");
            exit(1);
        }

    // Performing a sanity check to see if all the beads and structures are
    // correct.
    if (Check_System_Structure() == 0)
        {
            printf("Check structure sanity: OK\n");
        }
    else
        {
            printf("ERROR: wrong structure.\n");
            exit(1);
        }

    clock_t tEnd        = clock();
    double elapsed_time = (double) (tEnd - tStart) / (double) CLOCKS_PER_SEC;
    printf("Initialization done. %.2f sec elapsed.\n", elapsed_time);

    tStart = clock();

    long nGen;
    int nMCInfo; // MC accept/reject

    printf("%s Beginning MC Simulation %s\n", lBrace, rBrace);
    printf("******************************************\n\n");
    printf("_____________________\n");
    printf("Thermalizing system.\n");
    printf("---------------------\n\n");
    // Thermalizing the system.
    fCuTemp = fPreKT;
    //FILE Initialization
    FileIO_CreateRunningDataFiles();

    for (nGen = 0; nGen < nMCStepsForTherm; nGen++)
        { // Intentionally not performing any data acquisition in the thermalization phase.
            nMCInfo = MC_Step_Equil(fCuTemp);
            //        printf("(%d,%d)\n", nMCInfo / 12, nMCInfo % 2);
            DataPrinting_Thermalization(nGen);
        }
    /*
     * Post-thermalization
     * */
    FileIO_WriteRestart_ForThermalization();


    printf("____________________________\n");
    printf("System has been thermalized!\n");
    printf("----------------------------\n\n");
    /*
    The system has thermalized!
    */
    int run_cycle;
    // Going through the MC cycles.
    for (run_cycle = 0; run_cycle < nTot_CycleNum; run_cycle++)
        {
            /*
             * Pre run-cycle specific initialization.
             */
            fKT = fKT_Cycle[run_cycle];
            Calculate_Rot_Bias(fKT);
            FileIO_PreCycle_Init(run_cycle);
            for (nGen = 0; nGen < nMCStepsPerCycle / 2; nGen++)
                {
                    fCuTemp = nAnnealing_Mode == -1 ? fKT : Temperature_Function(nAnnealing_Mode, nGen);
                    nMCInfo = MC_Step(fCuTemp);
                    //            printf("(%d,%d)\n", nMCInfo / 12, nMCInfo % 2);
                    DataPrinting_DuringRunCycles(nGen, run_cycle);
                }

            for (nGen = nMCStepsPerCycle / 2; nGen <= nMCStepsPerCycle; nGen++)
                {
                    fCuTemp = nAnnealing_Mode == -1 ? fKT : Temperature_Function(nAnnealing_Mode, nGen);
                    nMCInfo = MC_Step(fCuTemp);
                    //            printf("(%d,%d)\n", nMCInfo / 12, nMCInfo % 2);
                    DataPrinting_DuringRunCycles(nGen, run_cycle);
                    DataAnalysis_DuringRunCycles(nGen, run_cycle);
                }
            /*
             * Post run-cycle specific cleanup.
             */
            nAnnealing_Mode = -1;
            FileIO_WriteRestart_ForRun(run_cycle);
            CopyData_All(run_cycle);
            Reset_Global_Arrays();
        }

    // Writing everything
    FileIO_Write_TotalSysProp(run_cycle);

    tEnd         = clock();
    elapsed_time = (double) (tEnd - tStart) / (double) CLOCKS_PER_SEC;

    double elapsed_minutes = elapsed_time / 60.;
    double elapsed_hours   = elapsed_minutes / 60.;
    double elapsed_days    = elapsed_hours / 24.;

    printf("____________________________________\n");
    printf("Simulation finished! Timings:\n");
    printf("%.2lf mins.\n", elapsed_minutes);
    printf("%.2lf hours.\n", elapsed_hours);
    printf("%.2lf days.\n", elapsed_days);
    printf("------------------------------------\n");
    printf("%s ENDING %s\n", lBrace, rBrace);
    printf("******************************************\n");

    return 0;
}
