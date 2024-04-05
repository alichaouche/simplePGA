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
#include "./headers/pMediansEncoding.h"
#include "./headers/Graph.h"
#include "./headers/prog_params.h"


///############################################
// Global Variables.
///############################################
int sommeTotalFlux = 0,
        nbrNoeuds = 0,
        nbr_constraint = 0,
        nbrArretes = 0,
        fixedNbrArretes = 0,
        nbrNoeudsInCluster = 0,
        fluxMatrix[tailleMax][tailleMax],
        neighborMatrix[1000][1000],
        buckUpFluxMatrix[tailleMax][tailleMax],
        *fluxVector,
        max_clusters = 0,
        min_clusters = 0,
        min_sizeCluster = 0,
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
///***********************************************
partitionBE *populationBE1 = NULL,*populationBE2 = NULL, *solutionDominanteBE = NULL,*bestSolutionOverRunsBE = NULL;
partitionDVTC *populationDVTC1 = NULL,*populationDVTC2 = NULL, *solutionDominanteDVTC = NULL,*bestSolutionOverRunsDVTC = NULL;
partitionFC *populationFC1 = NULL,*populationFC2 = NULL, *solutionDominanteFC = NULL,*bestSolutionOverRunsFC = NULL;
partitionSVTC *populationSVTC1 = NULL,*populationSVTC2 = NULL, *solutionDominanteSVTC = NULL,*bestSolutionOverRunsSVTC = NULL;
partitionVAE *populationVAE1 = NULL,*populationVAE2 = NULL, *solutionDominanteVAE = NULL,*bestSolutionOverRunsVAE = NULL;
partitionVAE *populationPMP1 = NULL,*populationPMP2 = NULL, *solutionDominantePMP = NULL,*bestSolutionOverRunsPMP = NULL;

///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
/// DECLARATION DES VARIABLES FILES POUR NOS FICHIERS
///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
char cheminBestSolutionBE[150],  cheminOptimalSolutionBE[150], cheminFichierDetailsProbleme[150];
FILE *bestSolutionsOverIterationBE, *optimalSolutionFileBE, *FichierDetailsProbleme;

char cheminBestSolutionDVTC[150],  cheminOptimalSolutionDVTC[150], cheminFichierDetailsProbleme[150];
FILE *bestSolutionsOverIterationDVTC, *optimalSolutionFileDVTC, *FichierDetailsProbleme;

#define f 1
int main()
{
    printf("===================================================================================\n");
    printf("          application de l'algorithme genetique pour le probleme GPP               \n");
    printf("    utilisation de différents types de codage pour la representation des solution  \n");
    printf("===================================================================================\n");


    int i,j,nbrGeneration=100, numeroGraphe = 11;
/*
    printf("Veuillez introduire le numero de graphe \n");
    scanf("%d",&numeroGraphe);
*/
    char cheminGraph[100], cheminInputFileConstrainte[100];
    FILE *graphFile, *inputFileConstrainte;


    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///iiii  SOME STATISTICS ON THE GRAPHS WEIGHTS iiiiiiiiiiiiiii
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    float moyenneDesPoids=0, variancePoids=0, standardDeviation=0;
    int minWeight, maxWeight;

    ///#######La boucle des graphes ##########
    //TODO (ALI): Graphs with CC = 6,7,14,20,21,23,32,34,36,37
    for( ; numeroGraphe <= __NBR_GRAPHS__; numeroGraphe++) {

        sprintf(cheminGraph, "../data/graphs/graph%d.txt", numeroGraphe);
        sprintf(cheminInputFileConstrainte, "../data/constraintes/constraintsFileGraph%d.txt", numeroGraphe);

        graphFile = fopen(cheminGraph, "r");

        if (graphFile == NULL) {
            //printf("graph number %d \n", numeroGraphe);
            perror("Impossible d'ouvrir le fichier fluxMatrix en Lecture\n");
            exit(EXIT_FAILURE);
        }

        inputFileConstrainte = fopen(cheminInputFileConstrainte, "r");
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
        nbr_constraint++;

        ///###### COHABITATION CONSTRAINT #######

        fscanf(inputFileConstrainte, "%d", &nbrCohabitationConstraintes);
        int nbrCC = nbrCohabitationConstraintes;
        if (nbrCohabitationConstraintes != 0) {
            nbr_constraint++;
            int nd, na;
            for (i = 0; i < nbrCohabitationConstraintes; i++) {
                fscanf(inputFileConstrainte, "%d %d", &nd, &na);
                if (fluxMatrix[nd][na] > 0) {
                    cohabitationConstraintes[i].nouedDepart = nd;
                    cohabitationConstraintes[i].nouedArrive = na;
                } else {
                    printf("%d %d should not cohabitate\n", nd, na);
                    nbrCC--;
                }
            }
        }
        ///###### NON COHABITATION CONSTRAINT #######
        fscanf(inputFileConstrainte, "%d", &nbrNonCohabitationConstraintes);
        if (nbrNonCohabitationConstraintes != 0) {
            nbr_constraint++;
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
        metrics metricsDVTC, metricsBE;
        metricsDVTC.runTime = metricsDVTC.ES = 0;
        metricsBE.runTime = metricsBE.ES = 0;

        while (nbrRun <= maxRun) {
            printf("\nParallel Genetic Algorithm Run N = %d\n",nbrRun);
            printf("=================================\n");
            metricsDVTC.runTime = metricsDVTC.ES = 0;
            metricsBE.runTime = metricsBE.ES = 0;

            for (migration  = 0; migration < 10; migration++){
                binaryEncoding(numeroGraphe, &metricsBE);
                vertexToClusterEncoding(numeroGraphe, &metricsDVTC);
                migrationFromDVTCtoBE();
                migrationFromBEtoDVTC();
            }
            writeOptimalSolutionInFileDVTC(solutionDominanteDVTC,optimalSolutionFileDVTC,nbrRun,metricsDVTC);

            printf("\nDVTC : cut size = %d\t run time = %.2f\t Evaluation to a solution = %.2f contraintes violees = %d\n",
                   solutionDominanteDVTC->coutCoupe, metricsDVTC.runTime, metricsDVTC.ES, solutionDominanteDVTC->contrainteViole);

            writeOptimalSolutionInFileBE(solutionDominanteBE,optimalSolutionFileBE,nbrRun,
                                         metricsBE,
                                         moyenneDesPoids,variancePoids,
                                         standardDeviation);
            printf("BE : cut size = %d\t run time = %.2f\t Evaluation to a solution = %.2f contraintes violees = %d\n",
                   solutionDominanteBE->coutCoupe, metricsBE.runTime, metricsBE.ES, solutionDominanteBE->contrainteViole);

            nbrRun++;
        }

        fclose(bestSolutionsOverIterationDVTC);
        fclose(optimalSolutionFileDVTC);
        fclose(bestSolutionsOverIterationBE);
        fclose(optimalSolutionFileBE);


        printf("\n#######################################################################################\n");
        printf("YOU CAN FIND THE RESULTS IN THE FOLLWING FAILE :\n");
        printf("best solutions over generation in : data/results/encodingScheme\n");
        printf("\n#######################################################################################\n\n");
        nbrRun = 1;


        fclose(graphFile);
        fclose(inputFileConstrainte);
    }
    return 0;
}
