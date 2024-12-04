#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "../headers/prog_params.h"
#include "../headers/typeDeclaration.h"
#include "../headers/communesFunctions.h"
#include "../headers/genetic_alg_PMP.h"
#include "../headers/Graph.h"
#include "../headers/compilationConditionnelle.h"
#define affectationArretes 0

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
        neighborMatrix[tailleMax][tailleMax],
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
extern neighbors *neighborsVector;

///***********************************************
extern partitionPMP *populationPMP1,*populationPMP2, *solutionDominantePMP,*bestSolutionOverRunsPMP;
extern partitionBE *solutionDominanteBE;
extern partitionDVTC *solutionDominanteDVTC;
///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
/// DECLARATION DES VARIABLES FILES POUR NOS FICHIERS
///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
extern char cheminBestSolutionPMP[150],  cheminOptimalSolutionPMP[150];
extern FILE *bestSolutionsOverIterationPMP, *optimalSolutionFilePMP;

///***************************************************************************
///                     les corps des fonctions
///***************************************************************************

///***************************************************************************
///generation of the first population
///***************************************************************************
void generatePopulationPMP(partitionPMP* populationPMP, int indiceFirstElt)
{
    int i,j,indice;
    ///********************************************
    /// generating zero vertices : medians
    ///********************************************
    for(i=indiceFirstElt; i < taillePopulation; i++)
    {
        (populationPMP+i)->id =i;
        /// affecter la valeur 0 pour tous les noeuds de la solution
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationPMP+i)->genotype[j] = 0;
        }

        //### nbr Medians and None Medians
        (populationPMP+i)->medians = max_clusters;
        // (populationPMP+i)->medians =rnd(0,max_clusters); // max_clusters = ceil(nbrNoeuds / max_sizeClusters)
        (populationPMP+i)->nonMedians = nbrNoeuds - (populationPMP+i)->medians;

        // 1 : is a median node | 0 : is a non median vertex
        for(j=0; j<(populationPMP+i)->medians; j++)
        {
            do
            {
                indice = rnd(0,nbrNoeuds);
            }
            while((populationPMP+i)->genotype[indice] == 1);
            (populationPMP+i)->genotype[indice] = 1;
        }
    }
}

///****************************************************************
/// get the list of medians (1) and non medians vertices (0)
///****************************************************************
void calculerMediansNonMediansVertices(partitionPMP* populationPMP)
{
    int i,j,testValueGenotype;
    //printf("calculerMediansNonMediansVertices : ...\n");
    for(i=0; i<taillePopulation; i++)
    {
        (populationPMP+i)->medians = 0;
        (populationPMP+i)->nonMedians = 0;
        for(j=0; j<nbrNoeuds; j++)
        {
            testValueGenotype = (populationPMP+i)->genotype[j];
            switch (testValueGenotype)
            {
            case 0 :
                (populationPMP+i)->nonMediansVertices[(populationPMP+i)->nonMedians] = j;
                (populationPMP+i)->nonMedians++;
                break;
            case 1 :
                (populationPMP+i)->mediansVertices[(populationPMP+i)->medians] = j;
                (populationPMP+i)->medians++;
                break;
            default :
                printf("Cette valeur n\'est pas prise en cons�diration \n");
                break;
            }
        }
    }
}

///****************************************************************
/// calculate the cut size and fitness value for each solution
///****************************************************************
void cutSizeAndFitnessPMP(partitionPMP* populationPMP)
{
    //printf("cutSizeAndFitnessPMP\n");
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationPMP+i)->coutCoupe = 0; /// r�initialisation de la valriable apr�s chaque it�ration (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                if((populationPMP+i)->phenotype[j] == (populationPMP+i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (populationPMP+i)->coutCoupe = (populationPMP+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationPMP+i)->coutCoupeNormalise = (populationPMP+i)->coutCoupe + ((nbr_constraint - (populationPMP+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationPMP+i)->coutCoupeNormalise;
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
            varianceCoutDeCoupeNormalise + pow(((populationPMP+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajout�
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationPMP+i)->expectedValue = (populationPMP+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationPMP+i)->expectedValue;
    }


    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationPMP+i)->fitness = (float)((populationPMP+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }
#else

    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/pMedianEncoding/expectedValuePMP.txt","w");
    for(i=0; i<taillePopulation ; i++)
    {
        (populationPMP+i)->fitness = (float)((populationPMP+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling

}

///****************************************************************
/// naturalSelection function using roulette wheel
///****************************************************************
void naturalSelectionPMP(partitionPMP* populationPMP1,partitionPMP* populationPMP2)
{
    //printf("naturalSelectionPMP\n");
    int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<elitismeRate; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationPMP1+maxFitness)->coutCoupeNormalise < (populationPMP1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationPMP2+i) = *(populationPMP1+maxFitness);
        ///printf("(populationPMP2+%d)->coutCoupeNormalise = %d \n",i,(populationPMP2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationPMP2, taillePopulation, sizeof *populationPMP2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationPMP2,tmpPopulation);

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
            sommeFitness = sommeFitness + (populationPMP1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationPMP2+i) = *(populationPMP1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationPMP2+i) = *(populationPMP1+j);
            }
        }
    }
    ///free(tmpPopulation);

}

///****************************************************************
/// Single point crossover function
///****************************************************************
void crossOverPMP(partitionPMP* populationPMP1, partitionPMP* populationPMP2)
{
    //printf("Cross Over \n");
    int i = 0,j;
    int choixInd1, choixInd2, choixLocus;
    ///***********************************************************************************
    for(i=0; i<elitismeRate; i++)
    {
        *(populationPMP2+i) = *(populationPMP1+i);
        (populationPMP2+i)->id = i;
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

        choixLocus = rnd(0,nbrNoeuds-1); /// -2 parce qu'on a commencer � partir du 0

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        {
            (populationPMP2+i)->genotype[j] = (populationPMP1+choixInd1)->genotype[j];
            (populationPMP2+i+1)->genotype[j] = (populationPMP1+choixInd2)->genotype[j];

        }
        ///*************************************************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        {
            (populationPMP2+i)->genotype[j] = (populationPMP1+choixInd2)->genotype[j];
            (populationPMP2+i+1)->genotype[j] = (populationPMP1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationPMP2+i)->id = i;
        (populationPMP2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }
    generatePopulationPMP(populationPMP2, taillePopulation-regeneration );
#else
    while (i < taillePopulation)
    {
        ///choixInd1 = rnd(elitismeRate,taillePopulation);
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            /// ce traitement �vite l'appariement d'un element avec lui m�me/
            /// choixInd2 = rnd(elitismeRate,taillePopulation);
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

        choixLocus = rnd(0,nbrNoeuds-1); /// -2 parce qu'on a commencer � partir du 0

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        {
            (populationPMP2+i)->genotype[j] = (populationPMP1+choixInd1)->genotype[j];
            (populationPMP2+i+1)->genotype[j] = (populationPMP1+choixInd2)->genotype[j];

        }
        ///*************************************************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        {
            (populationPMP2+i)->genotype[j] = (populationPMP1+choixInd2)->genotype[j];
            (populationPMP2+i+1)->genotype[j] = (populationPMP1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationPMP2+i)->id = i;
        (populationPMP2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }

#endif
}

///****************************************************************
/// find the best solution of the current population
///****************************************************************
int findTheBestSolutionPMP(partitionPMP *populationPMP)
{
    int maxFitness = populationPMP->coutCoupeNormalise;
    int indice=0;
    for(int i = 0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationPMP+i)->coutCoupeNormalise)
        {
            maxFitness = (populationPMP+i)->coutCoupeNormalise;
            indice = i;
        }
    }
#if mouchard
    printf("\n maxFitness = %0.4f \t indice = %d \n",maxFitness, indice);
#endif // mouchard
    return indice;
}

///****************************************************************
/// One single aleration mutation
///****************************************************************
void mutationPMP(partitionPMP* populationPMP)
{
    //printf("Mutation\n");
    int i,numeroGeneMute;
    double applicationOfMutation;
    for(i=elitismeRate; i<taillePopulation; i++)
    {
#if Windows
        applicationOfMutation= drand48ForWindows(0,1001);
        ///printf("applicationOfMutation = %0.4f \n",applicationOfMutation);
        ///mon_sleep(pause);
#else
        applicationOfMutation= drand48();
#endif // Windows
        if(applicationOfMutation <= tauxMutation)
        {
            /// choisir un g�ne al�atoirement pour lui appliquer la mutation
            numeroGeneMute = rnd(0,nbrNoeuds);
            ///((populationPMP+i)->genotype[numeroGeneMute] == 1)? 0 : 1;
            if((populationPMP+i)->genotype[numeroGeneMute] == 1)
            {
                (populationPMP+i)->genotype[numeroGeneMute] =0;
            }
            else
            {
                (populationPMP+i)->genotype[numeroGeneMute] =1;
            }

        }
    }
}

///****************************************************************
/// check the total flow sum function
///****************************************************************
float testerLaSommeDesFitnessPMP(partitionPMP* populationPMP)
{
    //printf("testerLaSommeDesFitness ...\n");
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationPMP+i)->fitness;
    }
    return sommeFitness;
}

///****************************************************************
/// Display the best solutions function
///****************************************************************
void displayTheBestSolutionPMP(partitionPMP* solutionDominante)
{
    int i;
    printf("id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
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
    printf("\n");
    for (int i = 0; i < nbrNoeuds; i++) {
        printf("%d\t",solutionDominante->phenotype[i]);
    }
    printf("\n");

    /*
    printf("\nMax clusters :%d\n",max_clusters);
    printf("\nMax clusters sizes :%d\n",max_sizeCluster);
    printf("\nClusters sizes :\n");
    for (int i = 0; i < solutionDominante->medians; i++) {
        printf("%d\t",solutionDominante->clusters[i].clusterSize);
    }
    printf("\nClusters formation :\n");
    for (int i = 0; i < solutionDominante->medians; i++) {
        for (int j = 0; j < solutionDominante->clusters[i].clusterSize; j++) {
            printf("%d\t",solutionDominante->clusters[i].clusterNoueds[j]);
        }
        printf("\n");
    }
    printf("\nConstraints : \n");
    for (int i = 0; i < 4; ++i) {
        printf("%d ", solutionDominante->constraintVector[i]);
    }
    */
}

///****************************************************************
/// check constraints and penalizing unfeasible solutions
///****************************************************************
void checkContrainstsPMP(partitionPMP *populationPMP)
{
    //printf("checkContrainstsPMP\n");
    int i,j,k;
    for(i=0; i<taillePopulation; i++)
    {
        ///initialisation de vecteur des contraintes et de la variable de contraintes violees
        (populationPMP+i)->contrainteViole=0;
        for(j=0; j < 4; j++)(populationPMP+i)->constraintVector[j]=0;
        ///***********************************************************************************

        if((populationPMP+i)->medians > max_clusters || (populationPMP+i)->medians < min_clusters)
        {
            (populationPMP+i)->constraintVector[0] += 1; /// le nombre de cluster n'est pas valide
        }
        if((populationPMP+i)->constraintVector[0] != 0)
        {
            (populationPMP+i)->contrainteViole++;
        }
        ///CHECK THE CLUSTER SIZE CONSTRAINTS
        for(j=0; j < (populationPMP+i)->medians; j++)
        {
            if((populationPMP+i)->clusters[j].clusterSize > max_sizeCluster ||
                (populationPMP+i)->clusters[j].clusterSize < min_sizeCluster)
            {
                (populationPMP+i)->constraintVector[1]++;
            }
        }
        if((populationPMP+i)->constraintVector[1] != 0)
        {
            (populationPMP+i)->contrainteViole++;
        }
        ///######################################
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstraintePMP(populationPMP);
    }

}

///****************************************************************
/// check the cohabitation and non cohabitation constraints
///****************************************************************
void checkCohabitationAndNonCohabitationConstraintePMP(partitionPMP *populationPMP)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        /// constraints vector is already initialize in the checkContrainstAndFitnessPenalizationPMP
        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationPMP+i)->phenotype[noeud1]!= (populationPMP+i)->phenotype[noeud2])
                {
                    (populationPMP+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes viol�es
            if((populationPMP+i)->constraintVector[2]!=0)
            {
                (populationPMP+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;

                if((populationPMP+i)->phenotype[noeud1]== (populationPMP+i)->phenotype[noeud2])
                {
                    (populationPMP+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes viol�es
            if((populationPMP+i)->constraintVector[3]!=0)
            {
                (populationPMP+i)->contrainteViole++;
            }
        }
    }

}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Assigning non medians vertices basing on the edges
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void pMedianEncodingAffectationParArretes(int numeroGraphe, metrics* metricsPMP)
{
    int bestSolution,iteration=1,bestSolutionIteration,ES=0;
    clock_t t1,t2;

    t1=clock();
if(migration==0) {
    generatePopulationPMP(populationPMP1,0);
    calculerMediansNonMediansVertices(populationPMP1);
    affectationDesNonMediansParArretes(populationPMP1);
    checkContrainstsPMP(populationPMP1);
    cutSizeAndFitnessPMP(populationPMP1);

    ///initialisation de la solution dominante
    bestSolution=findTheBestSolutionPMP(populationPMP1);
    *solutionDominantePMP= *(populationPMP1+bestSolution);
    bestSolutionIteration = 1;
    nbrApparition=1;
    writeBestSolutionInFilePMP(solutionDominantePMP,bestSolutionsOverIterationPMP,(10*migration+iteration));
}
    //displayTheBestSolutionPMP(solutionDominantePMP);
    for( ; iteration <= 10; iteration++) {
        printf(".");

        naturalSelectionPMP(populationPMP1,populationPMP2);
        crossOverPMP(populationPMP2, populationPMP1);
        mutationPMP(populationPMP1);
        calculerMediansNonMediansVertices(populationPMP1);
        affectationDesNonMediansParArretes(populationPMP1);
        checkContrainstsPMP(populationPMP1);
        cutSizeAndFitnessPMP(populationPMP1);
        ///=====================================================================================================
        /// dans ce qui suit, on va juste localiser la meilleur solution
        ///=====================================================================================================
        bestSolution=findTheBestSolutionPMP(populationPMP1);
        if((populationPMP1+bestSolution)->coutCoupeNormalise > solutionDominantePMP->coutCoupeNormalise )
        {
            *(solutionDominantePMP) = *(populationPMP1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else if((populationPMP1+bestSolution)->coutCoupeNormalise <= solutionDominantePMP->coutCoupeNormalise )
        {
            nbrApparition++;
        }
        writeBestSolutionInFilePMP(solutionDominantePMP,bestSolutionsOverIterationPMP,bestSolutionIteration);
        solutionsReproductionPMP(populationPMP1);
    }
    printf("|");
    t2=clock();
    //displayTheBestSolutionPMP(solutionDominantePMP);
    metricsPMP->runTime += (float)(t2-t1)/CLOCKS_PER_SEC;
    if(bestSolutionIteration >= 2){
        metricsPMP->ES += ((taillePopulation - elitismeRate) * (bestSolutionIteration - 2)) + (
                 taillePopulation + solutionDominantePMP->id - elitismeRate + 1);
    }
    else {
        metricsPMP->ES = solutionDominantePMP->id +1;
    }
}


///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Optimal solution is the best one over all the runs
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void writeOptimalSolutionInFilePMP(partitionPMP *solutionDominante,FILE* optimalSolutionFilePMP,
                                  int nbrRun, metrics metricsBE)
{
    if(optimalSolutionFilePMP == NULL) perror("the file was closed\n");
    int i,j,k;
    fprintf(optimalSolutionFilePMP,"%d,",nbrRun);
    fprintf(optimalSolutionFilePMP,"%d,%d,%d,%.3f,%.3f,",
        solutionDominante->coutCoupe,
        (sommeTotalFlux - solutionDominante->coutCoupe),
        sommeTotalFlux ,
        metricsBE.runTime,
        metricsBE.ES);
    fprintf(optimalSolutionFilePMP,"%d,%d,%d,",
        nbrNoeuds,
        nbrArretes,
        solutionDominante->contrainteViole);
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///RAJOUTER LES CONTRAINTES VIOLEES
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    fprintf(optimalSolutionFilePMP,"%d,%d,%d,%d,",
        solutionDominante->constraintVector[0],
        solutionDominante->constraintVector[1],
        solutionDominante->constraintVector[2],
        solutionDominante->constraintVector[3]);

    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///ECRIRE LA PARTITION CORRESPONDANTE A LA MEILLEUR SOLUTION
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK
    // int countNoued = 0;
    // for(j=0; j<solutionDominante->medians; j++){
    //     countNoued += solutionDominante->clusters[j].clusterSize;
    // }
    // if(countNoued != nbrNoeuds){
    //     printf("writing solutions  \n");
    //     printf("error : nbrNoued not conform ! %d != %d\n", nbrNoeuds, countNoued);
    //     exit(countNoued);
    // }
    ///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK
    for(j=0; j < solutionDominante->medians ; j++)
    {
        fprintf(optimalSolutionFilePMP,"{");
        for(k=0;k< solutionDominante->clusters[j].clusterSize; k++){
            fprintf(optimalSolutionFilePMP,"%d ",
                solutionDominante->clusters[j].clusterNoueds[k]);
        }
        if(j <= solutionDominante->medians - 1)fprintf(optimalSolutionFilePMP,"}\t");
    }
    fprintf(optimalSolutionFilePMP,"\n");
}


///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// This function is used by the qsort in a ascending way
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int compareCroissantFitnessPMP(void const *a, void const *b)
{
    partitionPMP const *pa = a;
    partitionPMP const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;;
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// This function is used by the qsort in descending way
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int compareDecroissantFitnessPMP(void const *a, void const *b)
{
    partitionPMP const *pa = a;
    partitionPMP const *pb = b;
    return  pa->coutCoupeNormalise - pb->coutCoupeNormalise;;
}


///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// write solution using the SVTC format
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int SVTC_Representation(int solution[])
{
    int i,j,k,m=nbrNoeuds, t1[m],t[m],sortIndex,tmp,exist,min;

    for(j=0; j<nbrNoeuds; j++)
    {
        t1[j] = -1;
    }

    t[0]=solution[0]; sortIndex=1;
    for(j=1;j<nbrNoeuds;j++){
        exist = 0;
        for(k=0;k<sortIndex;k++){
            if(solution[j]==t[k]){
                exist = 1;
                break;
            }
        }
        if(!exist){
            t[sortIndex]=solution[j];
            sortIndex++;
        }
    }

    for(j=0; j<sortIndex; j++)
    {
        for(k=0; k<nbrNoeuds; k++)
        {
            if (solution[k] == t[j] && t1[k] <0)
            {
                solution[k] = j;
                t1[k] = 1;
            }
        }
    }
    ///affichage de la solution reordonn�e
    for(i=0;i<nbrNoeuds;i++) printf("%d\t",solution[i]);
    printf("\n\n");


    return sortIndex; /// �a represente le nombre de cluster de la partition
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Calculate  solutionDominante clusters : I've to check this function
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void calculerLesClustersDeLaSolution(int C[],partitionBE *solutionDominante){

    int i,j,tmp, Ctmp[100],indice;
    for(i=0;i<nbrNoeuds;i++) Ctmp[i] = solutionDominante->phenotype[i];
    /// trie de phenotype
    for(i=0;i<nbrNoeuds-1;i++){
        for(j=i+1; j < nbrNoeuds;j++){
            if(Ctmp[i] > Ctmp[j]) {
                tmp = Ctmp[i];
                Ctmp[i] = Ctmp[j];
                Ctmp[j] = tmp;
            }
        }
    }


    C[0] = Ctmp[0];
    indice = 1;
    for(i=1; i<nbrNoeuds; i++)
    {
        if(Ctmp[i] != Ctmp[i-1])
        {
            C[indice] = Ctmp[i];
            indice++;
        }
    }
    solutionDominante->nbrCluster = indice;
}


///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
/// J'ai remarque que meme si les solutions sont realisables le nombre de clusters reste tout
/// de meme important du coup j'ai decide de contracter les clusters autant que possible
/// en respectant les contraintes de taille, et de non cohabitation.
/// iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
///contractionDesClusters
//////FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void contractionDesClusters(partitionBE *solutionDominante){

    printf("\ncontractionDesClusters ....\n");
    printf("cout de coupe avant contraction = %d\n",solutionDominante->coutCoupe);
    int i,j,k,nd,na,clustersNonContractable,C[1000];
    calculerLesClustersDeLaSolution(C,solutionDominante);
    ///**********************************************************************************

    for(i=0; i<solutionDominante->nbrCluster-1; i++){
        for(j=i+1; j<solutionDominante->nbrCluster; j++){
            if(solutionDominante->clusters[C[i]].clusterSize + solutionDominante->clusters[C[j]].clusterSize <= max_sizeCluster ){
                /// verification des contraintes de non cohabitation
                clustersNonContractable =0;
                for(k=0; k<nbrNonCohabitationConstraintes;k++){
                    nd = nonCohabitationConstraintes[k].nouedDepart;
                    na = nonCohabitationConstraintes[k].nouedArrive;
                    if( (solutionDominante->phenotype[nd] == C[i] && solutionDominante->phenotype[nd] == C[j]) ||
                       (solutionDominante->phenotype[nd] == C[j] && solutionDominante->phenotype[nd] == C[i])){
                            clustersNonContractable = 1;
                            printf("cluster %d et le cluster %d ne sont pas contracltable \n",C[i],C[j]);
                            break;
                    }
                }
                if(!clustersNonContractable){
                    printf("contraction des clusters %d et %d \n",C[i],C[j]);
                    for(k=0; k<nbrNoeuds;k++){
                        if(solutionDominante->phenotype[k] == C[j]) solutionDominante->phenotype[k] = C[i];
                    }
                    /// Mise � jour de la solution apr�s la contraction
                    solutionDominante->clusters[C[i]].clusterSize += solutionDominante->clusters[C[j]].clusterSize;
                    calculerLesClustersDeLaSolution(C,solutionDominante);
                }
            }
        }
    }
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// calculate the real value of the solutionDominantes cut size
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void realCutSizeForTheBestSolution(partitionPMP *solutionDominante){

    int i,j,sommeTotalCoutCoupe = 0;
    /// reinitialisation de la valriable apres chaque iteration (individus de la population)
    solutionDominante->coutCoupe = 0;
    for(i=0; i<nbrNoeuds-1; i++)
    {
        for(j=i+1; j<nbrNoeuds; j++)
        {
            if(solutionDominante->phenotype[i] == solutionDominante->phenotype[j] && fluxMatrix[i][j] > 0)
            {
                solutionDominante->coutCoupe = solutionDominante->coutCoupe +  fluxMatrix[i][j];
            }
        }
    }
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// calculateMaxDistance
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void calculateMaxDistance(int maxDistance[1000],int vertex, int phenotype[1000]){
    int i,j;
    for(i=0;i<nbrNoeuds;i++){
        maxDistance[i]=0;
    }
    /********************************************************************************/
    for(j=0;j<neighborMatrix[vertex][0];j++){
        if(phenotype[neighborMatrix[vertex][j]] != -1){
                //maxDistance[phenotype[neighborMatrix[vertex][j]]] += fluxMatrix[][];
        }
    }

}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Assignement of non medians using edges
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void affectationDesNonMediansParArretes(partitionPMP* populationPMP){
    //printf("affectationDesNonMediansParArretes\n");
    int i,l,maxDistance, indiceNotSeed, indiceSeed, assignedOpurtunity, nbrVoisins, voisin;
    int distanceArray[nbrNoeuds],nbrNotSeed,tmpNotSeed[nbrNoeuds],nbrFullCluster;

    // I'll see how to use it later
    // calculateNeighborsVector();

    for(i=0; i < taillePopulation; i++){
		for(int j=0;j<nbrNoeuds;j++){
			(populationPMP+i)->phenotype[j] = -1;
		}
        //*********************************************
        //if((populationPMP+i)->medians == 0) exit(5555);
        for(int l=0; l<(populationPMP+i)->medians; l++){
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // the cluster id is the median identifier
            (populationPMP+i)->clusters[l].clusterId = (populationPMP+i)->mediansVertices[l];
            (populationPMP+i)->clusters[l].clusterNoueds[0] = (populationPMP+i)->mediansVertices[l];
            (populationPMP+i)->clusters[l].clusterSize = 1;
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            (populationPMP+i)->phenotype[(populationPMP+i)->mediansVertices[l]] = l;
        }
        //*********************************************
        nbrNotSeed = (populationPMP+i)->nonMedians;
        nbrFullCluster = 0;
        for(l=0; l < nbrNotSeed; l++){
            /// the non affected vertices (non medians)
            tmpNotSeed[l] = (populationPMP+i)->nonMediansVertices[l];
        }
        //*********************************************
        int tmpNbrNotSeed;
        for(int j=0; j < (populationPMP+i)->nonMedians ;j++){
            if(tmpNotSeed[j] >= 0){ // the non median vertex is not assigned yet
               indiceNotSeed = tmpNotSeed[j];
               /// set to -1 the distanceArray vector
               for(int k=0;k<(populationPMP+i)->medians;k++)
                   distanceArray[k] = -1; // If needed I can get the median id from the mediansVertices array

               assignedOpurtunity = 0;
               nbrVoisins = neighborMatrix[indiceNotSeed][0];

               for(l=1; l<=nbrVoisins ; l++){
                 /// we must only check the neighborhood of the indiceNotSeed
                    voisin = neighborMatrix[indiceNotSeed][l];
                    if((populationPMP+i)->phenotype[voisin] != -1 /* this neighbor is already assigned */ &&
                           (populationPMP+i)->clusters[(populationPMP+i)->phenotype[voisin]].clusterSize < max_sizeCluster){
                        if(distanceArray[(populationPMP+i)->phenotype[voisin]] == -1)
                            distanceArray[(populationPMP+i)->phenotype[voisin]] = fluxMatrix[indiceNotSeed][voisin];
                        else
                            distanceArray[(populationPMP+i)->phenotype[voisin]] += fluxMatrix[indiceNotSeed][voisin];
                        assignedOpurtunity = 1;
                     }
               }
               if(assignedOpurtunity){
                   maxDistance = -1;
                   for(int k=0; k < (populationPMP+i)->medians ; k++){
                      if(distanceArray[k] != -1 && maxDistance < distanceArray[k]){
                          maxDistance = distanceArray[k];
                          indiceSeed = k;
                      }
                   }
                   if(maxDistance != -1){
                       (populationPMP+i)->phenotype[indiceNotSeed] = indiceSeed;
                       // Append the new node at the end of the array
                       (populationPMP+i)->clusters[indiceSeed].clusterNoueds[(populationPMP+i)->clusters[indiceSeed].clusterSize] = indiceNotSeed;
                       (populationPMP+i)->clusters[indiceSeed].clusterSize++;
                       if((populationPMP+i)->clusters[indiceSeed].clusterSize >= max_sizeCluster) nbrFullCluster++;
                            tmpNotSeed[j] = -1;
                            nbrNotSeed--;
                        }
                  }
                }
            }
        /// HEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHEREHERE
        //if the NbrNotSeed keep its value for at least two iteration
        //we assign the rest of the nodes to the non full clusters
        tmpNbrNotSeed = nbrNotSeed;
        while (nbrNotSeed > 0) {
                //iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
                // if all the neighbors of a given nodes are all assigned and their clusters are all full
                // we then should look for an unempty cluster to assign that vertex
                //iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
                int availableNeighbors[nbrNoeuds], nbrAvailableNeighbors, nbrAssignedNeighbors = 0, nbrNotAssignedNeighbors = 0;
                for (int l = 0; l < (populationPMP + i)->nonMedians; l++) {
                    nbrAvailableNeighbors = 0; nbrAssignedNeighbors = 0; nbrNotAssignedNeighbors=0;
                    if (tmpNotSeed[l] != -1) {// the node is not assigned yet
                        for (int k = 1; k <= neighborMatrix[tmpNotSeed[l]][0]; k++) {
                            // we check for each neighbor if it is already assigned and if its cluster host is full
                            if ((populationPMP + i)->phenotype[neighborMatrix[tmpNotSeed[l]][k]] != -1 &&
                                // -1 the neighbor has been assigned
                                (populationPMP + i)->clusters[(populationPMP + i)->phenotype[neighborMatrix[tmpNotSeed[
                                    l]][k]]].clusterSize >= max_sizeCluster) {
                                nbrAssignedNeighbors++;
                            } else if ((populationPMP + i)->phenotype[neighborMatrix[tmpNotSeed[l]][k]] != -1 &&
                                       // -1 the neighbor has been assigned
                                       (populationPMP + i)->clusters[(populationPMP + i)->phenotype[neighborMatrix[
                                           tmpNotSeed[l]][k]]].clusterSize < max_sizeCluster) {
                                availableNeighbors[nbrAvailableNeighbors] = neighborMatrix[tmpNotSeed[l]][k];
                                nbrAvailableNeighbors++;
                            } else {
                                nbrNotAssignedNeighbors++;
                            }
                        }

                        //If all the neighbors are assigned to full clusters
                        if (nbrAssignedNeighbors == neighborMatrix[tmpNotSeed[l]][0] || nbrNotAssignedNeighbors == neighborMatrix[tmpNotSeed[l]][0]) {
                            continue;
                        }
                        // we've neighbors assigned to not full clusters
                        for (int k = 0; k < (populationPMP + i)->medians; k++) {
                            distanceArray[k] = 0;
                            for (int m = 0; m < nbrAvailableNeighbors; m++) {
                                if ((populationPMP + i)->phenotype[availableNeighbors[m]] == k)
                                    distanceArray[k] += fluxMatrix[tmpNotSeed[l]][availableNeighbors[m]];
                            }
                        }
                        //here we've the distance array filled and we can use it
                        maxDistance = -1;
                        for (int k = 0; k < (populationPMP + i)->medians; k++) {
                            if (distanceArray[k] != 0 && maxDistance < distanceArray[k]) {
                                maxDistance = distanceArray[k];
                                indiceSeed = k;
                            }
                        }
                        if (maxDistance != -1) {
                            (populationPMP + i)->phenotype[tmpNotSeed[l]] = indiceSeed;
                            // Append the new node at the end of the array
                            (populationPMP + i)->clusters[indiceSeed].clusterNoueds[(populationPMP + i)->clusters[indiceSeed].clusterSize] = tmpNotSeed[l];
                            (populationPMP + i)->clusters[indiceSeed].clusterSize++;
                            if ((populationPMP + i)->clusters[indiceSeed].clusterSize >= max_sizeCluster) nbrFullCluster
                                    ++;
                            tmpNotSeed[l] = -1;
                            nbrNotSeed--;
                        }
                        //##############################################################################################
                    }
                }
                if (tmpNbrNotSeed == nbrNotSeed) {
                    break;
                } //the remaining nodes dont have an assinged neighbor
                tmpNbrNotSeed = nbrNotSeed;
        }
        if(nbrNotSeed > 0) {
            tmpNbrNotSeed = nbrNotSeed;
            for (int j = 0; j < (populationPMP+i)->nonMedians; j++) {
                if(tmpNotSeed[j] != -1) {
                    for (int k = 0; k < (populationPMP+i)->medians; k++) {
                        if((populationPMP+i)->clusters[k].clusterSize < max_sizeCluster) {
                            (populationPMP + i)->clusters[k].clusterNoueds[(populationPMP + i)->clusters[k].clusterSize]
                                    = tmpNotSeed[j];
                            (populationPMP+i)->phenotype[tmpNotSeed[j]] = k;
                            (populationPMP+i)->clusters[k].clusterSize++;
                            (populationPMP+i)->clusters[k].clusterId = (populationPMP+i)->mediansVertices[k];
                            tmpNotSeed[j] = -1;
                            nbrNotSeed--;
                            break;
                        }
                    }
                }
                else {
                    continue;
                }
                if(tmpNbrNotSeed == nbrNotSeed) {
                    break;
                }
                tmpNbrNotSeed = nbrNotSeed;
            }
            // medians = 0
            if(nbrNotSeed > 0) {
                (populationPMP+i)->medians = nbrNoeuds;
                (populationPMP+i)->nonMedians = 0;
                for (int j = 0; j < nbrNoeuds; j++) {
                    (populationPMP+i)->phenotype[j] = j;
                    (populationPMP+i)->clusters[j].clusterNoueds[0] = j;
                    (populationPMP+i)->clusters[j].clusterSize = 1;
                    (populationPMP+i)->clusters[j].clusterId = j;
                    tmpNotSeed[j] = -1;
                    nbrNotSeed--;
                }
            }
        }
    }

        //CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK
        /*
        i = 0;
        printf("max size cluster = %d\n", max_sizeCluster);
        printf("medians :\n");
        for (int k = 0; k < (populationPMP+i)->medians; k++) {
            printf("%d\t", (populationPMP+i)->mediansVertices[k]);
        }
        printf("\nnon medians :\n");
        for (int k = 0; k < (populationPMP+i)->nonMedians; k++) {
            printf("%d\t", (populationPMP+i)->nonMediansVertices[k]);
        }
        printf("\nphenotype :\n");

        for (int k = 0; k < nbrNoeuds ; k++) {
            printf("%d\t", (populationPMP+i)->phenotype[k]);
        }

        printf("\n");

        for (int k = 0; k < (populationPMP+i)->medians ; k++) {
            printf("clusterId = %d\t", (populationPMP+i)->clusters[k].clusterId);
            printf("cluster size = %d\n", (populationPMP+i)->clusters[k].clusterSize);
            printf("cluster vertices : \n");
            for (int l = 0; l < (populationPMP+i)->clusters[k].clusterSize; l++) {
                printf("%d\t",(populationPMP+i)->clusters[k].clusterNoueds[l]);
            }
            printf("\n");
        }
        exit(444);
        */
        //CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK
 	}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// This function  uses the same principale used by the two other encodings
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void writeBestSolutionInFilePMP(partitionPMP *solutionDominantePMP, FILE *bestSolutionsOverIterationPMP,int iteration)
{
    int i;
    fprintf(bestSolutionsOverIterationPMP,"%d,%d,%d,%d,%d,%d,%0.2f,%d,%d,",nbrRun,iteration,solutionDominantePMP->id,
            solutionDominantePMP->coutCoupe,(sommeTotalFlux - solutionDominantePMP->coutCoupe),solutionDominantePMP->coutCoupeNormalise,
            solutionDominantePMP->fitness, solutionDominantePMP->contrainteViole,solutionDominantePMP->medians);

    for(i=0; i<nbrNoeuds; i++)
    {
        fprintf(bestSolutionsOverIterationPMP,"%d ",solutionDominantePMP->genotype[i]);
    }
    fprintf(bestSolutionsOverIterationPMP,"\n");
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Reproduce a subset of 10 individuals of the population
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void solutionsReproductionPMP(partitionPMP* populationPMP) {
    int indice;
    qsort(populationPMP, taillePopulation, sizeof *populationPMP, compareDecroissantFitnessPMP);
    /// REPRODUCTION OF NEW INDIVIDUALS
    for (int i = taillePopulation - 10; i < taillePopulation; i++) {
        (populationPMP+i)->id =i;
        for(int j=0; j<nbrNoeuds; j++)
        {
            (populationPMP+i)->genotype[j] = 0;
        }
        (populationPMP+i)->medians =rnd(1,nbrNoeuds);
        for(int j=0; j<(populationPMP+i)->medians; j++) /// <= pour qu'on aura un �l�ment en plus pour le nombre des parties
        {
            do
            {
                indice = rnd(0,nbrNoeuds);
            }
            while((populationPMP+i)->genotype[indice] == 1);
            (populationPMP+i)->genotype[indice] = 1;
        }
        (populationPMP+i)->nonMedians = nbrNoeuds - (populationPMP+i)->medians;
    }
}

