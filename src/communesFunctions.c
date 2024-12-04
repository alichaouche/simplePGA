//
// Created by Ali on 01-Apr-24.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "../headers/typeDeclaration.h"
#include "../headers/communesFunctions.h"
#include "../headers/genetic_alg_BE.h"
#include "../headers/prog_params.h"

#define tailleMax 10000
#define dynamiqueAllocation 0

///############################################
// Global Variables.
///############################################
extern int sommeTotalFlux,
        nbrNoeuds,
        nbr_constraint,
        nbrArretes,
        fixedNbrArretes,
        nbrNoeudsInCluster,
        fluxMatrix[tailleMax][tailleMax],
        buckUpFluxMatrix[tailleMax][tailleMax],
        *fluxVector,
        max_clusters,
        min_clusters,
        min_sizeCluster,
        max_sizeCluster,
        bestSolution,
        nbrApparition,
        iteration,
        nbrCohabitationConstraintes,
        nbrNonCohabitationConstraintes,
        nbrRun;

extern float MBF, ART, AES, SDBF, BF[40];

extern edge *edgeVector,*fixedEdgeVector,cohabitationConstraintes[100],nonCohabitationConstraintes[100];
///***********************************************
extern char cheminBestSolutionBE[150],  cheminOptimalSolutionBE[150],
cheminBestSolutionDVTC[150],  cheminOptimalSolutionDVTC[150],
cheminBestSolutionPMP[150],  cheminOptimalSolutionPMP[150],
cheminBestSolutionPGA[150],  cheminOptimalSolutionPGA[150];

extern FILE *bestSolutionsOverIterationBE, *optimalSolutionFileBE,
*bestSolutionsOverIterationDVTC, *optimalSolutionFileDVTC,
*bestSolutionsOverIterationPMP, *optimalSolutionFilePMP,
*bestSolutionsOverIterationPGA, *optimalSolutionFilePGA;
///***********************************************
extern partitionBE *populationBE1,*populationBE2, *solutionDominanteBE,*bestSolutionOverRunsBE;
extern partitionDVTC *populationDVTC1, *populationDVTC2, *solutionDominanteDVTC, *bestSolutionOverRunsDVTC;
extern partitionFC *populationFC1, *populationFC2, *solutionDominanteFC, *bestSolutionOverRunsFC;
extern partitionSVTC *populationSVTC1, *populationSVTC2, *solutionDominanteSVTC, *bestSolutionOverRunsSVTC;
extern partitionVAE *populationVAE1, *populationVAE2, *solutionDominanteVAE, *bestSolutionOverRunsVAE;
extern partitionPMP *populationPMP1, *populationPMP2, *solutionDominantePMP, *bestSolutionOverRunsPMP;

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Memory allocation for populations
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void memoryAllocationForPopulations() {
    ///BE :
    populationBE1 = (partitionBE*)malloc(taillePopulation*sizeof (partitionBE));
    if(populationBE1==NULL) printf("memory allocation failed for the populationBE1 \n");

    populationBE2 = (partitionBE*)malloc(taillePopulation*sizeof(partitionBE));
    if(populationBE2==NULL) printf("memory allocation failed for the populationBE2 \n");

    solutionDominanteBE = (partitionBE*)malloc(sizeof(partitionBE));
    if(solutionDominanteBE==NULL) printf("memory allocation failed for the solutionDominanteBE \n");

    ///DVTC :
    populationDVTC1 = (partitionDVTC*)malloc(taillePopulation*sizeof (partitionDVTC));
    if(populationDVTC1==NULL) printf("memory allocation failed for the populationDVTC1 \n");

    populationDVTC2 = (partitionDVTC*)malloc(taillePopulation*sizeof(partitionDVTC));
    if(populationDVTC2==NULL) printf("memory allocation failed for the populationDVTC2 \n");

    solutionDominanteDVTC = (partitionDVTC*)malloc(sizeof(partitionDVTC));
    if(solutionDominanteDVTC==NULL) printf("memory allocation failed for the solutionDominanteDVTC \n");

    ///PMP :
    populationPMP1 = (partitionPMP*)malloc(taillePopulation*sizeof (partitionPMP));
    if(populationPMP1==NULL) printf("memory allocation failed for the populationPMP1 \n");

    populationPMP2 = (partitionPMP*)malloc(taillePopulation*sizeof(partitionPMP));
    if(populationPMP2==NULL) printf("memory allocation failed for the populationPMP2 \n");

    solutionDominantePMP = (partitionPMP*)malloc(sizeof(partitionPMP));
    if(solutionDominantePMP==NULL) printf("memory allocation failed for the solutionDominantePMP \n");

}


///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Generate random values in a specific range
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int rnd(int v1, int v2)
{
    if(v1 == v2)
    {
        printf("la fonction rnd : v1 == v2 could not generate a random value");
    }
    return (rand()%(v2 - v1)+v1);;
}
///END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Generate random values in the range [0,1]
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
double drand48ForWindows(int v1, int v2)
{
    return (double)(rnd(v1,v2)) /100.0;
}
///END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Generate random values in a specific range
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void mon_sleep( int nbr_seconds )
{
    clock_t goal;
    goal = ( nbr_seconds * CLOCKS_PER_SEC ) +  clock();
    while( goal > clock() )
    {
        ; /* loop */
    }

}
///END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Used by the qsort function for an ascendin sorting
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int compareTrieCroissant (void const *a, void const *b)
{
    /* definir des pointeurs type's et initialise's
       avec les parametres */
    int pa = *(const int*)a;
    int pb = *(const int*)b;

    /* evaluer et retourner l'etat de l'evaluation (tri d�scroissant) */
    return  pa-pb;
}
///END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Used by the qsort function for an decending sorting
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int compareTrieDecroissant (void const *a, void const *b)
{
    /* definir des pointeurs type's et initialise's
       avec les parametres */
    int pa = *(const int*)a;
    int pb = *(const int*)b;
    /* evaluer et retourner l'etat de l'evaluation (tri d�scroissant) */
    return  pb-pa;
}
///END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Save all the informations of the graph in a file
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void writeDetailsProblemInFile(FILE* detailsFile, int seuil, edge deletedEdges[],int nbrDeletedEdge)
{
    int i,j;

    ///###################  GRAPH BEFORE MODIFICATION ########################
    fprintf(detailsFile,"THE ORIGINAL GRAPH\n");
    for(i=0; i<nbrNoeuds; i++)
    {
        for(j=0; j<nbrNoeuds; j++)
        {
            if(buckUpFluxMatrix[i][j] == -1)
                fprintf(detailsFile,"%d;",0);
            else
                fprintf(detailsFile,"%d;",buckUpFluxMatrix[i][j]);
        }
        fprintf(detailsFile,"\n");
    }

    ///###################  GRAPH AFTER MODIFICATION ########################
    fprintf(detailsFile,"THE NEW GRAPH AFTER DELETING SOMME EDGES \n");
    for(i=0; i<nbrNoeuds; i++)
    {
        for(j=0; j<nbrNoeuds; j++)
        {
            if(fluxMatrix[i][j] == -1)
                fprintf(detailsFile,"%d;",0);
            else
                fprintf(detailsFile,"%d;",fluxMatrix[i][j]);
        }
        fprintf(detailsFile,"\n");
    }
    ///###################  INFOS ABOUT THE GRAPH ##############################
    fprintf(detailsFile,"LES SOMMETS;LES ARETES;SOMME DES PIODS\n");
    fprintf(detailsFile,"%d;%d;%d\n\n", nbrNoeuds, nbrArretes, sommeTotalFlux);

    ///###################  COHABITATIONS CONSTRAINTS #####################
    if(nbrCohabitationConstraintes != 0)
    {
        fprintf(detailsFile,"CONTRAINTES DE COHABITATION \n");
        fprintf(detailsFile,"SOMMET 1; SOMMET 2\n");
        for(i=0; i<nbrCohabitationConstraintes; i++)
            fprintf(detailsFile,"%d;%d\n",cohabitationConstraintes[i].nouedDepart,cohabitationConstraintes[i].nouedArrive);
    }
    fprintf(detailsFile,"\n\n");

    ///###################  NON COHABITATIONS CONTRAINTS ######################
    if(nbrNonCohabitationConstraintes != 0)
    {
        fprintf(detailsFile,"CONTRAINTES DE NON COHABITATION : \n");
        fprintf(detailsFile,"SOMMET 1; SOMMET 2\n");
        for(i=0; i<nbrNonCohabitationConstraintes; i++)
            fprintf(detailsFile,"%d;%d\n",nonCohabitationConstraintes[i].nouedDepart,nonCohabitationConstraintes[i].nouedArrive);
    }

    ///###################  DELETED EDGES ############################
    fprintf(detailsFile,"ARETES SUPPRIMEES \n");
    fprintf(detailsFile,"ARETES;POIDS\n");
    for(i=0; i<nbrDeletedEdge; i++)
        fprintf(detailsFile,"{%d,%d};%d\n",deletedEdges[i].nouedDepart,deletedEdges[i].nouedArrive,deletedEdges[i].poids);

}
///END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC_END_FUNC



///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Get the right tuples of vertices for CC and NCC
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void getCohabitationCouples(int maxWeight){
    int nbrNCC = 0, nbrCC = 0;
    for(int i=0; i<nbrNoeuds-1;i++){
        for(int j = i+1; j< nbrNoeuds; j++){
            if (fluxMatrix[i][j] == maxWeight && nbrCC < 5){
                printf("CC = %d %d\n",i,j);
                nbrCC++;
                i++;
            }
            else if(fluxMatrix[i][j] < 0 && nbrNCC < 5 && i > nbrNoeuds/2){
                printf("NCC = %d %d\n",i,j);
                nbrNCC++;
                i++;
            }

        }
    }
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Create files headers
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void createFilesHeader(){
    //###################  BEST SOLUTION PER ITERATION #################
    fprintf(bestSolutionsOverIterationBE,"Run num,Iteration num,Dominante Solution ID,Intra Cut Size,"
            "Inter Cut Size,Normalised Cut Size,Fitness,Nbr Violated Constraints,Numbre of clusters,Genotype,Phenotype \n");

    ///###################  OPTIMAL SOLUTION PER RUN ####################
    fprintf(optimalSolutionFileBE,"nbrRun,intra cut size,inter cut size,total flow,run time,ES,nbr vertices,nbr edges,nbr violated constraints,nbr-Clusters,M-S-C,Coh Con, Non Coh Con,Partition\n");

    //###################  BEST SOLUTION PER ITERATION #################
    fprintf(bestSolutionsOverIterationDVTC,"Run num,Iteration num,Dominante Solution ID,Intra Cut Size,"
            "Inter Cut Size,Normalised Cut Size,Fitness,Nbr Violated Constraints,Numbre of clusters,Genotype,Phenotype \n");

    ///###################  OPTIMAL SOLUTION PER RUN ####################
    fprintf(optimalSolutionFileDVTC,    "nbrRun,intra cut size,inter cut size,total flow,run time,ES,nbr vertices,nbr edges,nbr violated constraints,nbr-Clusters,M-S-C,Coh Con, Non Coh Con,Partition\n");

    //###################  BEST SOLUTION PER ITERATION #################
    fprintf(bestSolutionsOverIterationPMP,"Run num,Iteration num,Dominante Solution ID,Intra Cut Size,"
            "Inter Cut Size,Normalised Cut Size,Fitness,Nbr Violated Constraints,Numbre of clusters,Genotype,Phenotype \n");

    ///###################  OPTIMAL SOLUTION PER RUN ####################
    fprintf(optimalSolutionFilePMP,    "nbrRun,intra cut size,inter cut size,total flow,run time,ES,nbr vertices,nbr edges,nbr violated constraints,nbr-Clusters,M-S-C,Coh Con, Non Coh Con,Partition\n");

    ///###################  OPTIMAL SOLUTION PER RUN ####################
    fprintf(optimalSolutionFilePGA,"nbrRun,intra cut size,inter cut size,total flow,run time,ES,nbr vertices,nbr edges,nbr violated constraints,nbr-Clusters,M-S-C,Coh Con, Non Coh Con,Partition\n");
}


void openFilesForOptimalAndBestSolutions(int numeroGraphe) {
    sprintf(cheminBestSolutionBE,"../data/results/BE/bestSolutions/bestSolutionBE%d.csv",numeroGraphe);
    sprintf(cheminOptimalSolutionBE,"../data/results/BE/OptimalSolution/OptimalSolutionBE%d.csv",numeroGraphe);
    bestSolutionsOverIterationBE=fopen(cheminBestSolutionBE,"w");
    optimalSolutionFileBE=fopen(cheminOptimalSolutionBE,"w");
    if(bestSolutionsOverIterationBE == NULL)
    {
        printf("An error occurred: This is a custom error message\n");
        exit(EXIT_FAILURE);
    }

    if(optimalSolutionFileBE== NULL)
    {
        perror("Impossible d'ouvrir le fichier OptimalSolutionBE en ecriture\n");
            exit(EXIT_FAILURE);
    }
    sprintf(cheminBestSolutionDVTC,"../data/results/DVTC/bestSolutions/bestSolutionDVTC%d.csv",numeroGraphe);
    sprintf(cheminOptimalSolutionDVTC,"../data/results/DVTC/OptimalSolution/OptimalSolutionDVTC%d.csv",numeroGraphe);
    bestSolutionsOverIterationDVTC=fopen(cheminBestSolutionDVTC,"w");
    optimalSolutionFileDVTC=fopen(cheminOptimalSolutionDVTC,"w");
    if(bestSolutionsOverIterationDVTC == NULL)
    {
        printf("An error occurred: This is a custom error message\n");
        exit(EXIT_FAILURE);
    }

        if(optimalSolutionFileDVTC== NULL)
        {
            perror("Impossible d'ouvrir le fichier OptimalSolutionBE en ecriture\n");
            exit(EXIT_FAILURE);
        }

        sprintf(cheminBestSolutionPMP,"../data/results/PMP/bestSolutions/bestSolutionPMP%d.csv",numeroGraphe);
        sprintf(cheminOptimalSolutionPMP,"../data/results/PMP/OptimalSolution/OptimalSolutionPMP%d.csv",numeroGraphe);
        bestSolutionsOverIterationPMP=fopen(cheminBestSolutionPMP,"w");
        optimalSolutionFilePMP=fopen(cheminOptimalSolutionPMP,"w");
        if(bestSolutionsOverIterationPMP == NULL)
        {
            printf("An error occurred when opening PMP file: This is a custom error message\n");
            exit(EXIT_FAILURE);
        }

        if(optimalSolutionFilePMP== NULL)
        {
            perror("Impossible d'ouvrir le fichier OptimalSolutionPMP en ecriture\n");
            exit(EXIT_FAILURE);
        }
    /// IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    /// When dealing with pGA, we dont need to store best solutions overIteration
    /// Instead, for each run we compare optimal solutions and the best one is
    /// assigned to the pGA
    /// IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
        sprintf(cheminOptimalSolutionPGA,"../data/results/pGA/OptimalSolutionPGA%d.csv",numeroGraphe);
        optimalSolutionFilePGA=fopen(cheminOptimalSolutionPGA,"w");
        if(optimalSolutionFilePGA== NULL)
        {
            perror("Impossible d'ouvrir le fichier OptimalSolutionPGA en ecriture\n");
            exit(EXIT_FAILURE);
        }

}

