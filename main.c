/**
* Author : Ali CHAOUCHE
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <errno.h>

#include "./headers/typeDeclaration.h"
#include "./headers/communesFunctions.h"
#include "./headers/compilationConditionnelle.h"

#include "./headers/genetic_alg_BE.h"
#include "./headers/genetic_alg_SVTC.h"
#include "./headers/genetic_alg_FC.h"
#include "./headers/genetic_alg_DVTC.h"
#include "./headers/genetic_alg_PMP.h"
#include "./headers/Graph.h"
#include "./headers/prog_params.h"
#include "headers/migration.h"


///############################################
// Global Variables.
///############################################
int sommeTotalFlux = 0,
        nbrNoeuds = 0,
        nbr_constraint = 4,
        nbrArretes = 0,
        fixedNbrArretes = 0,
        nbrNoeudsInCluster = 0,
        fluxMatrix[tailleMax][tailleMax],
        neighborMatrix[tailleMax][tailleMax],
        buckUpFluxMatrix[tailleMax][tailleMax],
        *fluxVector,
        max_clusters = 0,
        min_clusters = 2,
        min_sizeCluster = 2,
        max_sizeCluster = 0,
        bestSolution = 0,
        nbrApparition = 0,
        iteration = 0,
        nbrCohabitationConstraintes = 0,
        nbrNonCohabitationConstraintes = 0,
        nbrRun = 1,
        migration = 0;

///***************************************07102017
float MBF=0, ART=0, AES=0, SDBF=0, BF[40];
///***********************************************
edge *edgeVector,*fixedEdgeVector,cohabitationConstraintes[100],nonCohabitationConstraintes[100];
neighbors *neighborsVector;
///***********************************************
partitionBE *populationBE1 = NULL,*populationBE2 = NULL, *solutionDominanteBE = NULL,*bestSolutionOverRunsBE = NULL;
partitionDVTC *populationDVTC1 = NULL,*populationDVTC2 = NULL, *solutionDominanteDVTC = NULL,*bestSolutionOverRunsDVTC = NULL;
partitionFC *populationFC1 = NULL,*populationFC2 = NULL, *solutionDominanteFC = NULL,*bestSolutionOverRunsFC = NULL;
partitionSVTC *populationSVTC1 = NULL,*populationSVTC2 = NULL, *solutionDominanteSVTC = NULL,*bestSolutionOverRunsSVTC = NULL;
partitionVAE *populationVAE1 = NULL,*populationVAE2 = NULL, *solutionDominanteVAE = NULL,*bestSolutionOverRunsVAE = NULL;
partitionPMP *populationPMP1 = NULL,*populationPMP2 = NULL, *solutionDominantePMP = NULL,*bestSolutionOverRunsPMP = NULL;

///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
/// DECLARATION DES VARIABLES FILES POUR NOS FICHIERS
///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
char cheminBestSolutionBE[150],  cheminOptimalSolutionBE[150], cheminFichierDetailsProbleme[150];
FILE *bestSolutionsOverIterationBE, *optimalSolutionFileBE, *FichierDetailsProbleme;

char cheminBestSolutionDVTC[150],  cheminOptimalSolutionDVTC[150];
FILE *bestSolutionsOverIterationDVTC, *optimalSolutionFileDVTC;

char cheminBestSolutionPMP[150],  cheminOptimalSolutionPMP[150];
FILE *bestSolutionsOverIterationPMP, *optimalSolutionFilePMP;

char cheminBestSolutionPGA[150],  cheminOptimalSolutionPGA[150];
FILE *bestSolutionsOverIterationPGA, *optimalSolutionFilePGA;

int main()
{
    printf("==================================================================\n");
    printf("          Solving the GPP using parallel genetic algorithm        \n");
    printf("==================================================================\n");
    srand(time(NULL));
    int i, numeroGraphe = 1;
    memoryAllocationForPopulations();
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///iiii  SOME STATISTICS ON THE GRAPHS WEIGHTS iiiiiiiiiiiiiii
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    float moyenneDesPoids=0, variancePoids=0, standardDeviation=0;
    int minWeight, maxWeight;

    ///#######La boucle des graphes #########
    for(numeroGraphe = 1 ; numeroGraphe <= __NBR_GRAPHS__; numeroGraphe++) {

        char cheminGraph[100], cheminInputFileConstrainte[100];
        sprintf(cheminGraph, "../data/graphs/graph%d.txt", numeroGraphe);
        sprintf(cheminInputFileConstrainte, "../data/constraintes/constraintsFileGraph%d.txt", numeroGraphe);
        FILE *graphFile = fopen(cheminGraph, "r");

        if (graphFile == NULL) {
            printf("graph number %d \n", numeroGraphe);
            perror("Impossible d'ouvrir le fichier fluxMatrix en Lecture\n");
            exit(EXIT_FAILURE);
        }

        FILE *inputFileConstrainte = fopen(cheminInputFileConstrainte, "r");
        if (inputFileConstrainte == NULL) {
            perror("Impossible d'ouvrir le fichier inputFileConstrainte en Lecture\n");
            exit(EXIT_FAILURE);
        }

        ///########## lecture de la matrice des flux and nbrNoueds  ########//
        readFlowMatrixAndNbrNoeudsFromFlowMatrixFile(graphFile);
        if (!isFlowMatrixSymetric()) {
            printf("flow matrix not symetric \n");
            exit(-4);
        }

        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        ///iiiii STATISTICS ON THE GRAPH WEIGHTS iiiiiiiiiiiiiiiiiiiii
        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        moyenneDesPoids = 0, variancePoids = 0, standardDeviation = 0;

        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        /// CALCULATE THE MIN AND THE MAX WEIGHT IN THE GRAPH
        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

        maxMinWeight(&minWeight, &maxWeight);

        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        /// READ THE CONSTRAINTS FROM THE CONSTRAINTS FILE
        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        fscanf(inputFileConstrainte, "%d", &max_sizeCluster);
        min_sizeCluster = 2;
        //nbr_constraint++;

        ///###### COHABITATION CONSTRAINT #######

        fscanf(inputFileConstrainte, "%d", &nbrCohabitationConstraintes);
        int nbrCC = nbrCohabitationConstraintes;
        if (nbrCohabitationConstraintes != 0) {
            //nbr_constraint++;
            int nd, na;
            for (i = 0; i < nbrCohabitationConstraintes; i++) {
                fscanf(inputFileConstrainte, "%d %d", &nd, &na);
                if (fluxMatrix[nd][na] > 0) {
                    cohabitationConstraintes[i].nouedDepart = nd;
                    cohabitationConstraintes[i].nouedArrive = na;
                } else {
                    fluxMatrix[nd][na] = 0;
                    cohabitationConstraintes[i].nouedDepart = nd;
                    cohabitationConstraintes[i].nouedArrive = na;
                }
            }
        }

        ///###### NON COHABITATION CONSTRAINT #######
        fscanf(inputFileConstrainte, "%d", &nbrNonCohabitationConstraintes);
        if (nbrNonCohabitationConstraintes != 0) {
            //nbr_constraint++;
            for (i = 0; i < nbrNonCohabitationConstraintes; i++) {
                fscanf(inputFileConstrainte, "%d %d", &nonCohabitationConstraintes[i].nouedDepart,
                       &nonCohabitationConstraintes[i].nouedArrive);
            }

        }

        ///############ Connecter le graphe #############
        //TO DO: CHECK THE GRAPH'S CONNECTEDNESS
        //fixedNbrArretes = calculerNombreArrete();
        graphConnexion();

        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        /// TOTAL NUMBER OF EDGES AND THEIR ARRAY #
        /// WE CALL IT fixedNbrArretes EVEN IF IT IS CALCULATED UPON
        /// THE fluxMatrix BECAUSE WE WILL USE ANOTHER VARIABLE
        /// nbrAretes WHERE WE DELETE SOME EGES ACCORDING TO THEIR FLOW
        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        nbrArretes = fixedNbrArretes = calculerNombreArrete();

        edgeVector = (edge*)malloc(fixedNbrArretes*sizeof(edge));
        if(edgeVector == NULL) perror("allocation memory failed for edge vector\n");
        calculateEdgeVector();
        //displayEdgeVector();

        ///######### Nombre max des clusters #########
        max_clusters = nbrNoeuds;
        //max_clusters = ceil((float) (nbrNoeuds) / (float)(max_sizeCluster));
        ///######### Calculer la somme totale des poids ##########
        sommeTotalFlux = calculSommeTotalFlux();
        //printf("sommeTotalFlux = %d\n",sommeTotalFlux);

        ///######### Weigth's informations ##################
        moyenneDesPoids = 0, variancePoids = 0, standardDeviation = 0;
        moyenneDesPoids = (float) sommeTotalFlux / (float) fixedNbrArretes;
        //printf("%.4f", moyenneDesPoids);

        for (int i = 0; i < fixedNbrArretes; i++) {
            //variancePoids += fixedEdgeVector[i].poids;
            variancePoids += pow((edgeVector[i].poids - moyenneDesPoids), 2);
        }
        variancePoids = variancePoids / (float) fixedNbrArretes;
        standardDeviation = sqrt(variancePoids);
        //printf("moyenneDesPoids = %.2f variancePoids = %.2f\n",moyenneDesPoids,variancePoids);

        MBF = 0;
        ART = 0;
        AES = 0;
        SDBF = 0;
        metrics metricsDVTC, metricsBE, metricsPMP;
        metricsDVTC.runTime = metricsDVTC.ES = 0;
        metricsBE.runTime = metricsBE.ES  = 0;
        metricsPMP.runTime = metricsPMP.ES  = 0;
        ///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHE
        /// To be checked : 25/20/2024
        ///#### OPEN OUTPUT FILES HERE TO ADD SOME ADDITIONAL INFO ####
        openFilesForOptimalAndBestSolutions(numeroGraphe);

        ///#### Add the header for the files (Columns names) ####
        createFilesHeader();

        //##### Calculate the matrix of neighbors
        calculateNeighborhoodMatrix();
        ///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK


        while (nbrRun <= maxRun) {
            printf("\nParallel Genetic Algorithm Run N = %d\n",nbrRun);
            printf("=================================\n");
            metricsDVTC.runTime = metricsDVTC.ES = 0;
            metricsBE.runTime = metricsBE.ES = 0;
            metricsPMP.runTime = metricsPMP.ES = 0;

            for (migration  = 0; migration < 7; migration++){
                //binaryEncoding(numeroGraphe, &metricsBE);
                //vertexToClusterEncoding(numeroGraphe, &metricsDVTC);
                pMedianEncodingAffectationParArretes(numeroGraphe, &metricsPMP);
#if pGA
                migrationFromDVTCtoBE();
                migrationFromDVTCtoPMP();
                migrationFromBEtoDVTC();
                migrationFromBEtoPMP();
                migrationFromPMPtoDVTC();
                migrationFromPMPtoBE();

#endif
            }
            //==========================================================================================================
#if DVTC
            writeOptimalSolutionInFileDVTC(solutionDominanteDVTC,optimalSolutionFileDVTC,nbrRun,metricsDVTC);
            printf("\nDVTC : cut size = %d\t run time = %.2f\t Evaluation to a solution = %.2f contraintes violees = %d\n",
                   solutionDominanteDVTC->coutCoupe, metricsDVTC.runTime, metricsDVTC.ES, solutionDominanteDVTC->contrainteViole);
            for(int i=0; i<4; i++)
                printf("%d\t",solutionDominanteDVTC->constraintVector[i]);
            printf("\n");
#endif
            //==========================================================================================================
#if BE
            writeOptimalSolutionInFileBE(solutionDominanteBE,optimalSolutionFileBE,nbrRun,metricsBE);
            printf("BE : cut size = %d\t run time = %.2f\t Evaluation to a solution = %.2f contraintes violees = %d\n",
                   solutionDominanteBE->coutCoupe, metricsBE.runTime, metricsBE.ES, solutionDominanteBE->contrainteViole);
#endif
            //==========================================================================================================
#if PMP
            writeOptimalSolutionInFilePMP(solutionDominantePMP,optimalSolutionFilePMP,nbrRun, metricsPMP);
            printf("\nPMP : cut size = %d\t run time = %.2f\t Evaluation to a solution = %.2f contraintes violees = %d\n",
                   solutionDominantePMP->coutCoupe, metricsPMP.runTime, metricsPMP.ES, solutionDominantePMP->contrainteViole);
#endif
            //==========================================================================================================

#if pGA
       /// iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
       /// compare optimale solutions and the best
       /// one will be taken as the optimal solution for the pGA
       /// iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
       if(solutionDominantePMP->coutCoupeNormalise >= solutionDominanteBE->coutCoupeNormalise &&
        solutionDominantePMP->coutCoupeNormalise >= solutionDominanteDVTC->coutCoupeNormalise) {
           printf("PMP win\n");
           writeOptimalSolutionInFilePMP(solutionDominantePMP,optimalSolutionFilePGA,nbrRun, metricsPMP);
       }
       else if(solutionDominanteBE->coutCoupeNormalise >= solutionDominantePMP->coutCoupeNormalise &&
        solutionDominanteBE->coutCoupeNormalise >= solutionDominanteDVTC->coutCoupeNormalise) {
           printf("BE win\n");
           writeOptimalSolutionInFileBE(solutionDominanteBE,optimalSolutionFilePGA,nbrRun, metricsBE);
       }

       else {
           printf("DVTC win\n");
           writeOptimalSolutionInFileDVTC(solutionDominanteDVTC,optimalSolutionFilePGA,nbrRun, metricsDVTC);
       }

#endif
            nbrRun++;
        }

        fclose(bestSolutionsOverIterationDVTC);
        fclose(optimalSolutionFileDVTC);
        fclose(bestSolutionsOverIterationBE);
        fclose(optimalSolutionFileBE);
        fclose(bestSolutionsOverIterationPMP);
        fclose(optimalSolutionFilePMP);
        fclose(optimalSolutionFilePGA);

        printf("\n#######################################################################################\n");
        printf("YOU CAN FIND THE RESULTS IN THE FOLLWING FAILE :\n");
        printf("best solutions over generation in : data/results/encodingScheme\n");
        printf("\n#######################################################################################\n\n");
        nbrRun = 1;


        fclose(graphFile);
        fclose(inputFileConstrainte);
    }
    return 1;
}
