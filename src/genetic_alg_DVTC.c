#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../headers/typeDeclaration.h"
#include "../headers/communesFunctions.h"
#include "../headers/genetic_alg_DVTC.h"
#include "../headers/compilationConditionnelle.h"
#include "../headers/prog_params.h"

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
        nbrRun,
        migration;

extern float MBF, ART, AES, SDBF, BF[40];

extern edge *edgeVector,*fixedEdgeVector,cohabitationConstraintes[100],nonCohabitationConstraintes[100];
///***********************************************
extern partitionDVTC *populationDVTC1,*populationDVTC2, *solutionDominanteDVTC,*bestSolutionOverRunsDVTC;
extern char cheminBestSolutionDVTC[150],  cheminOptimalSolutionDVTC[150], cheminFichierDetailsProbleme[150];
extern FILE *bestSolutionsOverIterationDVTC, *optimalSolutionFileDVTC, *FichierDetailsProbleme;
extern partitionBE *solutionDominanteBE;

///***********************************************
///generation of the initial population
///***********************************************

void generatePopulationDVTC(partitionDVTC* populationDVTC, int indiceFirstElt)
{
    int i,j;
    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        for(j=0; j<nbrNoeuds; j++)
        {
            /// on va essayer de limiter le nombre de clusters pour voir ce que �a va donner
            ///(populationDVTC+i)->genotype[j] = rnd(0,nbrParties);
            (populationDVTC+i)->genotype[j] = rnd(0,nbrNoeuds); /// [0;nbrNoeuds-1]
        }
        (populationDVTC+i)->id =i;
    }
}
///**************************************************************************************
/// affichage des population
void affichePopulationDVTC(partitionDVTC* populationDVTC)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        printf("\n id = %d \t cout de coupe = %d \t fitness = %0.4f \t",
               (populationDVTC+i)->id,(populationDVTC+i)->coutCoupe, (populationDVTC+i)->fitness);
        printf("le cout de coupe normalise = %d \n",(populationDVTC+i)->coutCoupeNormalise);
        printf("la partition est  = \t");
        for(j=0; j<nbrNoeuds; j++)
        {
            printf("%d ",(populationDVTC+i)->genotype[j]);
        }
        printf("\n");
        printf("le nombre des clusters est %d \n",(populationDVTC+i)->nbrCluster);
        printf("le nombre des contraintes violees est  = %d \n",(populationDVTC+i)->contrainteViole);
        printf("affichage des tailles des clusters \n");
        for(j=0; j<(populationDVTC+i)->nbrCluster; j++)
        {
            printf("%d ",(populationDVTC+i)->clustersSize[j]);
        }
        printf("\n============================================================\n");

    }
}
///************************************************************************************************************
/// calculer le cout de coupe de chaque partition ainsi que son fitness
///***************************************************************************************
/// on va appliquer la formule suivante pour le calcul des couts de coupe
/// T'(C) = B - T(C)            <=> la minimisation devient une maximisation
/// Z(C)=T'(C)+(u-v(C)) �B      <=>  calculer la nouvelle fitness par rapport aux contrainte
/// T'(C) : cout de coupe (maximiser les flux INTRA), u ; nombre de contrainte totale
/// v(C): nombre de contrainte viol�,  B: la somme totale des flux
///***************************************************************************************
void calculCoutCoupeEtFitnessDVTC(partitionDVTC* populationDVTC)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationDVTC+i)->coutCoupe = 0; /// r�initialisation de la valriable apr�s chaque it�ration (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                /// on peut avoir deux noeuds du m�me cluster mais qui sont par relier par une arr�tes
                if((populationDVTC+i)->genotype[j] == (populationDVTC+i)->genotype[k] && fluxMatrix[j][k] >= 0)
                {
                    (populationDVTC+i)->coutCoupe = (populationDVTC+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationDVTC+i)->coutCoupeNormalise = (populationDVTC+i)->coutCoupe + ((nbr_constraint - (populationDVTC+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationDVTC+i)->coutCoupeNormalise;
    }
#if scaling


    ///======================================================================================================
    /// Scaling Fintnss : sigma scaling Melanie mitchelle r�f�rence Goldberg | le 26/11/2015
    ///pour r�gler le probl�me des valeurs n�gatives des expected values ==> utiliser sigma truncation de goldberg
    /// g(f) = f + (moyenne - c * sigma )
    /// 1 <= c <= 5, g�n�ralement initilis� � 2
    ///======================================================================================================
    float varianceCoutDeCoupeNormalise = 0, moyenneCoutDeCoupeNormalise, sigma, sommeTotalOfExpectdValues = 0;
    int c= 2;

    moyenneCoutDeCoupeNormalise = sommeTotalCoutCoupe / taillePopulation;
    for(i=0; i<taillePopulation; i++)
    {
        varianceCoutDeCoupeNormalise =
            varianceCoutDeCoupeNormalise + pow(((populationDVTC+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajout�
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationDVTC+i)->expectedValue = (populationDVTC+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationDVTC+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationDVTC+i)->fitness = (float)((populationDVTC+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/vertexToClusterEncoding/expectedValueDVTC.txt","w");
    for(i=0; i<taillePopulation ; i++)
    {
/**
        if(iteration == 1)
        {
            fprintf(file,"%d \t %d \t %0.2f \n",i,(populationDVTC+i)->coutCoupeNormalise,(populationDVTC+i)->expectedValue);
        }
*/
        (populationDVTC+i)->fitness = (float)((populationDVTC+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling
}
///**************************************************************************************
void naturalSelectionDVTC(partitionDVTC* populationDVTC1,partitionDVTC* populationDVTC2)
{
    int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<elitismeRate; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationDVTC1+maxFitness)->coutCoupeNormalise < (populationDVTC1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationDVTC2+i) = *(populationDVTC1+maxFitness);
        ///printf("(populationDVTC2+%d)->coutCoupeNormalise = %d \n",i,(populationDVTC2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationDVTC2, taillePopulation, sizeof *populationDVTC2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationDVTC2,tmpPopulation);

    for(i=elitismeRate; i<taillePopulation; i++)
    {
#if Windows
        lotterie= drand48ForWindows(0,1001);
#else
        lotterie= drand48();
#endif // Windows
        sommeFitness = 0;

        for(j=0; j<taillePopulation; j++)
        {
            sommeFitness = sommeFitness + (populationDVTC1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationDVTC2+i) = *(populationDVTC1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationDVTC2+i) = *(populationDVTC1+j);
            }
        }
    }
    ///free(tmpPopulation);

}
///**************************************************************************************
void crossOverDVTC(partitionDVTC* populationDVTC1, partitionDVTC* populationDVTC2)
{
    int i = 0,j =0;
    int choixInd1, choixInd2, choixLocus;
    ///***********************************************************************************
    for(i=0; i<elitismeRate; i++)
    {
        *(populationDVTC2+i) = *(populationDVTC1+i);
        (populationDVTC2+i)->id = i;
    }
    ///***********************************************************************************

#if NEW_CROSSOVER
    while (i < taillePopulation-regeneration)
    {
        ///choixInd1 = rnd(elitismeRate,taillePopulation);
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            /// ce traitement �vite l'appariement d'un element avec lui m�me/
            ///choixInd2 = rnd(elitismeRate,taillePopulation);
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );
        choixLocus = rnd(0,nbrNoeuds-1); /// l'intervalle sera entre 1 et nbrNoeuds-2 (la borne superieur -1)
        for(j=0; j<=choixLocus; j++)
        {
            (populationDVTC2+i)->genotype[j] = (populationDVTC1+choixInd1)->genotype[j];
            (populationDVTC2+i+1)->genotype[j] = (populationDVTC1+choixInd2)->genotype[j];
        }
        ///**********************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)
        {
            (populationDVTC2+i)->genotype[j] = (populationDVTC1+choixInd2)->genotype[j];
            (populationDVTC2+i+1)->genotype[j] = (populationDVTC1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationDVTC2+i)->id = i;
        (populationDVTC2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }
    generatePopulationDVTC( populationDVTC2, taillePopulation-regeneration);
#else
    while (i < taillePopulation)
    {
        ///choixInd1 = rnd(elitismeRate,taillePopulation);
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            /// ce traitement �vite l'appariement d'un element avec lui m�me/
            ///choixInd2 = rnd(elitismeRate,taillePopulation);
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );
        choixLocus = rnd(0,nbrNoeuds-1); /// l'intervalle sera entre 1 et nbrNoeuds-2 (la borne superieur -1)
        for(j=0; j<=choixLocus; j++)
        {
            (populationDVTC2+i)->genotype[j] = (populationDVTC1+choixInd1)->genotype[j];
            (populationDVTC2+i+1)->genotype[j] = (populationDVTC1+choixInd2)->genotype[j];
        }
        ///**********************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)
        {
            (populationDVTC2+i)->genotype[j] = (populationDVTC1+choixInd2)->genotype[j];
            (populationDVTC2+i+1)->genotype[j] = (populationDVTC1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationDVTC2+i)->id = i;
        (populationDVTC2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }

#endif
}
///**************************************************************************************
int findTheBestSolutionDVTC(partitionDVTC *populationDVTC)
{
    /// la fitness est calcul� � partir des couts de coupe normalis�s
    float maxFitness = populationDVTC->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationDVTC+i)->fitness)
        {
            maxFitness = (populationDVTC+i)->fitness;
            indice = i;
        }
    }
    return indice;
}
///**************************************************************************************
void mutationDVTC(partitionDVTC* populationDVTC)
{

    int i,numeroGeneMute=-1;
    double applicationOfMutation;
    for(i=elitismeRate; i<taillePopulation; i++)
    {
#if Windows
        applicationOfMutation = drand48ForWindows(0,1001);
#else
        applicationOfMutation = drand48();
#endif // Windows


        if(applicationOfMutation <= tauxMutation)
        {
            /// choisir un g�ne al�atoirement pour lui appliquer la mutation
            numeroGeneMute = rnd(0,nbrNoeuds);
            (populationDVTC+i)->genotype[numeroGeneMute] = rnd(0,max_clusters);

        }
    }

}
///**************************************************************************************
float testerLaSommeDesFitnessDVTC(partitionDVTC* populationDVTC)
{
#if mouchard
    printf("testerLaSommeDesFitness ...\n");
#endif // mouchard
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationDVTC+i)->fitness;
    }
    return sommeFitness;
}
///************************************************************************************************************
void displayTheBestSolutionDVTC(partitionDVTC* solutionDominante)
{
    int i;

    printf("\nid = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
           solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    printf("le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte sur le nombre de clusters : OK \n");
            else printf("contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte sur la taille des clusters : OK \n");
            else printf("contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte de cohabitation : OK \n");
            else printf("contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte Non Cohabitation : OK \n");
            else printf("contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            printf("Cette contrainte n'est pas prise en charge par le syst�me \n");

        }

    }
}

///**********************************************************************************
void writeSolutionInFileDVTC(partitionDVTC *populationDVTC, FILE *outputFilePop, int iteration)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///�criture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.4f\t%d\t%d\t%d\t",iteration,(populationDVTC+i)->id,
                (populationDVTC+i)->coutCoupe,(populationDVTC+i)->fitness,(populationDVTC+i)->coutCoupeNormalise,
                (populationDVTC+i)->contrainteViole,(populationDVTC+i)->nbrCluster);
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationDVTC+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t");

        for(j=0; j<nbr_constraint; j++)
        {
            fprintf(outputFilePop,"%d ",(populationDVTC+i)->constraintVector[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    fprintf(outputFilePop,"\n\n");

}
///**************************************************************************************
void writeBestSolutionInFileDVTC(partitionDVTC *solutionDominante, FILE *bestSolutionsOverIterationDVTC,int iteration)
{
    int i;
    fprintf(bestSolutionsOverIterationDVTC,"%d,%d,%d,%d,%d,%d,%0.2f,%d,%d,",nbrRun,iteration,solutionDominante->id,
            solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),solutionDominante->coutCoupeNormalise,
            solutionDominante->fitness, solutionDominante->contrainteViole,solutionDominante->nbrCluster);

    for(i=0; i<nbrNoeuds; i++)
    {
        fprintf(bestSolutionsOverIterationDVTC,"%d ",solutionDominante->genotype[i]);
    }
    fprintf(bestSolutionsOverIterationDVTC,"\n");
}
///======================================================================================
void vertexToClusterEncoding(int numeroGraphe, metrics *metricsDVTC)
{

    int bestSolution,iteration=1 ,bestSolutionIteration,ES = 0;
    clock_t t1,t2;
    double temps;
    if(nbrRun == 1){
        ///###################  OPTIMAL SOLUTION PER RUN ####################
        sprintf(cheminOptimalSolutionDVTC, "../data/results/DVTC/OptimalSolution/OptimalSolutionDVTC%d.csv", numeroGraphe);
        optimalSolutionFileDVTC = fopen(cheminOptimalSolutionDVTC, "w");

        if (optimalSolutionFileDVTC == NULL) {
            perror("Impossible d'ouvrir le fichier OptimalSolutionDVTC en ecriture\n");
            exit(EXIT_FAILURE);
        }

        fprintf(optimalSolutionFileDVTC, "Run N,"
                                         "intra cut size , inter cut size, Total flows,"
                                         "Run time,ES,"
                                         "Nbr vertices , Nbr Edges,Nbr Violated constraintes, Cont T_M_C,Cont Cohab, Cont NON Coha, Partition \n");
    }

    if(migration == 0){
        ///############ ALLOCATION DE LA MEMOIRE POUR LES VARIABLES POPULATION ##########
        if(populationDVTC1 == NULL) populationDVTC1 = (partitionDVTC*)malloc(taillePopulation*sizeof (partitionDVTC));
        if(populationDVTC1==NULL) printf("memory allocation failed for the populationDVTC1 \n");

        if(populationDVTC2 == NULL) populationDVTC2 = (partitionDVTC*)malloc(taillePopulation*sizeof(partitionDVTC));
        if(populationDVTC2==NULL) printf("memory allocation failed for the populationDVTC2 \n");

        if(solutionDominanteDVTC == NULL)solutionDominanteDVTC = (partitionDVTC*)malloc(sizeof(partitionDVTC));
        if(solutionDominanteDVTC==NULL) printf("memory allocation failed for the solutionDominanteDVTC \n");

        if(bestSolutionOverRunsDVTC == NULL)bestSolutionOverRunsDVTC = (partitionDVTC*)malloc(sizeof(partitionDVTC));
        if(bestSolutionOverRunsDVTC==NULL) printf("memory allocation failed for the bestSolutionOverRuns \n");

        ///#### OPEN OUTPUT FILES HERE TO ADD SOME ADDITIONAL INFO ####
        sprintf(cheminBestSolutionDVTC, "../data/results/DVTC/bestSolutions/bestSolutionDVTC%d.csv", numeroGraphe);

        bestSolutionsOverIterationDVTC = fopen(cheminBestSolutionDVTC, "w");

        if (bestSolutionsOverIterationDVTC == NULL) {
            perror("Impossible d'ouvrir le fichier bestSolutionDVTC.csv en ecriture\n");
            exit(EXIT_FAILURE);
        }

        ///###################  BEST SOLUTION PER ITERATION #################
        fprintf(bestSolutionsOverIterationDVTC,
                "Numéro de l'exécution,Numéro de l'itération,ID Solution Dominante,Cout de coupe intra,"
                "Cout de Coupe inter,Cout de coupe Normalisé,Fitness,Nbr contraintes violées,Nombre de cluster,Genotype,Phenotype \n");

        ///####################################################################
        /// The best solution of the previous run is kept in the new generation
        ///####################################################################
        t1=clock(); ///START THE CHRONOMETER
        generatePopulationDVTC(populationDVTC1,0);
        checkContrainstAndFitnessPenalizationDVTC(populationDVTC1);
        calculCoutCoupeEtFitnessDVTC(populationDVTC1);

        bestSolution=findTheBestSolutionDVTC(populationDVTC1);
        *solutionDominanteDVTC=*(populationDVTC1+bestSolution);
        bestSolutionIteration= 1;
        nbrApparition=1;
        writeBestSolutionInFileDVTC(solutionDominanteDVTC,bestSolutionsOverIterationDVTC,(10*migration+iteration));
    }
    else{
        t1=clock(); ///START THE CHRONOMETER
    }
    for(iteration = 1 ; iteration <= 10; iteration++) {
        //do{
        ///=====================================================================================================
        printf(".");

        /// application des trois operateurs de bases de genetiques
        naturalSelectionDVTC(populationDVTC1, populationDVTC2);
        crossOverDVTC(populationDVTC2, populationDVTC1); /// nbrNoeuds repr�sente la taille des genotype
        mutationDVTC(populationDVTC1);
        checkContrainstAndFitnessPenalizationDVTC(populationDVTC1);
        calculCoutCoupeEtFitnessDVTC(populationDVTC1);
        bestSolution = findTheBestSolutionDVTC(populationDVTC1);
        if ((populationDVTC1 + bestSolution)->coutCoupeNormalise > solutionDominanteDVTC->coutCoupeNormalise)
            /** le 23/12/2015 : pour avoir la meilleur solution avec le moins de contrainte viol�es*/
            /// && (populationSDVTC1+bestSolution)->contrainteViole   < solutionDominanteDVTC->contrainteViole)
        {
            *(solutionDominanteDVTC) = *(populationDVTC1 + bestSolution);
            nbrApparition = 1;
            bestSolutionIteration = iteration;
        } else if ((populationDVTC1 + bestSolution)->coutCoupeNormalise <= solutionDominanteDVTC->coutCoupeNormalise) {
            nbrApparition++;
        }
        //fprintf(bestSolutionsOverIterationDVTC,"\n");
        writeBestSolutionInFileDVTC(solutionDominanteDVTC,bestSolutionsOverIterationDVTC,bestSolutionIteration);

        solutionsReproductionDVTC(populationDVTC1);
        //iteration++;
        //}while(nbrApparition < max_steady_generation);
    }

    t2=clock();
    //printf("\n");
    //displayTheBestSolutionDVTC(solutionDominanteDVTC);
    metricsDVTC->runTime += (float)(t2-t1)/(float)CLOCKS_PER_SEC;
    if(bestSolutionIteration>=2){
        metricsDVTC->ES += ((taillePopulation-elitismeRate)*(bestSolutionIteration-2))+(taillePopulation+ solutionDominanteDVTC->id -elitismeRate+1);
    }
    else {
        metricsDVTC->ES += solutionDominanteDVTC->id +1;
    }

}

///**************************************************************************************
void checkContrainstAndFitnessPenalizationDVTC(partitionDVTC *populationDVTC)
{
    int i,j,k;
    int tailleClusters;

    for(i=0; i<taillePopulation; i++)
    {
        (populationDVTC+i)->nbrCluster=0; /// le nombre de cluster doit etre initialis� � 1 pour prendre en cons�diration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes viol�es
        (populationDVTC+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationDVTC+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************

        (populationDVTC+i)->nbrCluster=0;
        for(j=0; j<nbrNoeuds; j++)
        {
            tailleClusters =0;
            for(k=0; k<nbrNoeuds; k++){
                if ((populationDVTC+i)->genotype[k] == j){
                    (populationDVTC+i)->clusters[(populationDVTC+i)->nbrCluster].clusterNoueds[tailleClusters]=k;
                    tailleClusters++;
                }
            }
            if (tailleClusters !=0){
                (populationDVTC+i)->clusters[(populationDVTC+i)->nbrCluster].clusterId = (populationDVTC+i)->nbrCluster;
                (populationDVTC+i)->clusters[(populationDVTC+i)->nbrCluster].clusterSize = tailleClusters;
                (populationDVTC+i)->nbrCluster++;
            }
        }
        /*-----------------------------------------------------------------*/


        for(j=0; j < (populationDVTC+i)->nbrCluster; j++)
        {
            if((populationDVTC+i)->clusters[j].clusterSize  > max_sizeCluster)
            {
                (populationDVTC+i)->constraintVector[0]++;

            }
        }
        if((populationDVTC+i)->constraintVector[0] != 0)
        {
            (populationDVTC+i)->contrainteViole++;
        }
    }

    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstrainteDVTC(populationDVTC);
    }

}
///**************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteDVTC(partitionDVTC *populationDVTC)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
        /// Cette tache est d�j� r�alis� au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
        (populationDVTC+i)->constraintVector[2]=0;
        (populationDVTC+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationDVTC+i)->genotype[noeud1]!= (populationDVTC+i)->genotype[noeud2])
                {
                    (populationDVTC+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes viol�es
            if((populationDVTC+i)->constraintVector[2]!=0)
            {
                (populationDVTC+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;
                if((populationDVTC+i)->genotype[noeud1]== (populationDVTC+i)->genotype[noeud2])
                {
                    (populationDVTC+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes viol�es
            if((populationDVTC+i)->constraintVector[3]!=0)
            {
                (populationDVTC+i)->contrainteViole++;
            }
        }


    }

}

///********************************************************************************************
void writeOptimalSolutionInFileDVTC(partitionDVTC *solutionDominante,FILE* optimalSolutionFileBE,
                                  int nbrRun, metrics metricsDVTC)
{
    int i,j,k;
    fprintf(optimalSolutionFileBE,"%d,",nbrRun);
    fprintf(optimalSolutionFileBE,"%d,%d,%d,%.3f,%.3f,", solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,metricsDVTC.runTime, metricsDVTC.ES);
    fprintf(optimalSolutionFileBE,"%d,%d,%d,",nbrNoeuds,nbrArretes,solutionDominante->contrainteViole);
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///RAJOUTER LES CONTRAINTES VIOLEES
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    fprintf(optimalSolutionFileDVTC,"%d,%d,%d,",solutionDominante->constraintVector[0],
            solutionDominante->constraintVector[1],
            solutionDominante->constraintVector[2]);

    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///ECRIRE LA PARTITION CORRESPONDANTE A LA MEILLEUR SOLUTION
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK

    int countNoued = 0;
    for(j=0; j<solutionDominante->nbrCluster; j++){
        countNoued += solutionDominante->clusters[j].clusterSize;
    }
    if(countNoued < nbrNoeuds){
        printf("writing solutions  \n");
        printf("error : nbrNoued not conform ! %d != %d\n", nbrNoeuds, countNoued);
        exit(countNoued);
    }

    ///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK
    for(j=0; j < solutionDominante->nbrCluster ; j++)
    {
        fprintf(optimalSolutionFileBE,"{");
        for(k=0;k< solutionDominante->clusters[j].clusterSize; k++){
            fprintf(optimalSolutionFileBE,"%d ",solutionDominante->clusters[j].clusterNoueds[k]);
        }
        fprintf(optimalSolutionFileBE,"}\t");
    }
    fprintf(optimalSolutionFileBE,"\n");
}
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
///FFFFFFFFFFF compareCroissantFitnessDVTC FFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
int compareCroissantFitnessDVTC (void const *a, void const *b)
{

    partitionDVTC const *pa = a;
    partitionDVTC const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}

///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
///FFFFFFFFFFF compareDecroissantFitnessDVTC FFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
int compareDecroissantFitnessDVTC (void const *a, void const *b)
{

    partitionDVTC const *pa = a;
    partitionDVTC const *pb = b;
    return  pa->coutCoupeNormalise - pb->coutCoupeNormalise;
}
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
///FFFFFFFFFFF solutionsReproduction FFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
void solutionsReproductionDVTC(partitionDVTC* populationDVTC) {

    qsort(populationDVTC, taillePopulation, sizeof *populationDVTC, compareDecroissantFitnessDVTC);
    /// REPRODUCTION OF NEW INDIVIDUALS
    for (int i = taillePopulation - 10; i < taillePopulation; i++) {
        /// THE SOLUTION ID
        (populationDVTC + i)->id = i;
        for(int j=0; j<nbrNoeuds; j++)
        {
            (populationDVTC+i)->genotype[j] = rnd(0,max_clusters); /// [0;nbrNoeuds-1]
        }

    }
}

///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
///FFFFFFFFFFF solutionsReproduction FFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
void migrationFromBEtoDVTC(){
    qsort(populationDVTC1, taillePopulation, sizeof *populationDVTC1, compareDecroissantFitnessDVTC);
    for(int i = 0; i <nbrNoeuds; i++)
        (populationDVTC1+taillePopulation-1)->genotype[i] = solutionDominanteBE->phenotype[i];

    (populationDVTC1+taillePopulation-1)->id = taillePopulation - 1;
    (populationDVTC1+taillePopulation-1)->coutCoupe = solutionDominanteBE->coutCoupe;
    (populationDVTC1+taillePopulation-1)->coutCoupeNormalise = solutionDominanteBE->coutCoupeNormalise;
    (populationDVTC1+taillePopulation-1)->contrainteViole = solutionDominanteBE->contrainteViole;
    (populationDVTC1+taillePopulation-1)->nbrCluster = solutionDominanteBE->nbrCluster;

    for(int i = 0; i <(populationDVTC1+taillePopulation-1)->nbrCluster; i++)
        (populationDVTC1+taillePopulation-1)->clusters[i] = solutionDominanteBE->clusters[i];

    for(int i = 0; i < nbr_constraint; i++)
        (populationDVTC1+taillePopulation-1)->constraintVector[i] = solutionDominanteBE->constraintVector[i];
    int sommeTotalCoutCoupe = 0;
    for(int i = 0; i < taillePopulation; i++)
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationDVTC1+i)->coutCoupeNormalise;

    float varianceCoutDeCoupeNormalise = 0, moyenneCoutDeCoupeNormalise, sigma, sommeTotalOfExpectdValues = 0;
    int c= 2;

    moyenneCoutDeCoupeNormalise = sommeTotalCoutCoupe / taillePopulation;
    for(int i=0; i<taillePopulation; i++)
    {
        varianceCoutDeCoupeNormalise =
        varianceCoutDeCoupeNormalise + pow(((populationDVTC1+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);

    /// calcul de sigma truncation
    (populationDVTC1+taillePopulation-1)->expectedValue = (populationDVTC1+taillePopulation-1)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);

    for(int i=0; i< taillePopulation;i++)
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationDVTC1+i)->expectedValue;

    for(int i=0; i<taillePopulation ; i++)
    (populationDVTC1+i)->fitness = (float)((populationDVTC1+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);

}


