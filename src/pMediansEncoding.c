#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "../headers/typeDeclaration.h"
#include "../headers/communesFunctions.h"
#include "../headers/pMediansEncoding.h"
#include "../headers/compilationConditionnelle.h"
#include "../headers/prog_params.h"
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
        neighborMatrix[1000][1000],
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

extern partitionPMP *populationPMP1 = NULL,*populationPMP2 = NULL, *solutionDominantePMP = NULL,*bestSolutionOverRunsPMP = NULL;

///***************************************************************************
///                     les corps des fonctions
///***************************************************************************

///****************************************************************************************************************************
///g�n�ration de la populationPMP initiale
void generatePopulationPMP(partitionPMP* populationPMP, int indiceFirstElt)
{

    int i,j,indice;
    /// g�n�ration des zero vertices
    for(i=indiceFirstElt; i< taillePopulation; i++)
    {
        (populationPMP+i)->id =i;
        /// affecter la valeur 0 � tous les noeuds de la solution
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationPMP+i)->genotype[j] = 0;
        }
        ///(populationPMP+i)->medians =2;
        (populationPMP+i)->medians =rnd(1,nbrNoeuds);
        ///(populationPMP+i)->medians =rnd(0,max_clusters); // laisser le nombre des medians libre donne des resumtats meilleurs que de le limit�
        ///printf("le nombre de medains dans cette solutions : %d \n",(populationPMP+i)->medians);
        /// choisir al�atoirement les noeuds medians
        for(j=0; j<(populationPMP+i)->medians; j++) /// <= pour qu'on aura un �l�ment en plus pour le nombre des parties
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
///****************************************************************************************************************************
/// affichage des population
void affichePopulationPMP(partitionPMP* populationPMP)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        if((populationPMP+i)->contrainteViole == 0)
        {
            printf("\nid = %d \t", (populationPMP+i)->id);
            printf("la genotypes est  = \t");
            for(j=0; j<nbrNoeuds; j++)
            {
                printf("%d ",(populationPMP+i)->genotype[j]);
            }
            printf("\n");
            /// affichage des medians de cette solution
            printf("  les medians sont  = \t");
            for(j=0; j<(populationPMP+i)->medians; j++)
            {
                printf("%d ",(populationPMP+i)->mediansVertices[j]);
            }
            printf("\n");

            printf("  les NOT  medians sont  = \t");
            for(j=0; j<(populationPMP+i)->nonMedians; j++)
            {
                printf("%d ",(populationPMP+i)->nonMediansVertices[j]);
            }
            printf("\n");

            printf("\n la phenotypes est  = \t");
            for(j=0; j<nbrNoeuds; j++)
            {
                printf("%d ",(populationPMP+i)->phenotype[j]);
            }
            printf("\n");
            printf("la fitness de la %d solution est : %0.2f",i,(populationPMP+i)->fitness );
            printf("\n**************************************************************************\n");
        }
        else if((populationPMP+i)->contrainteViole == 1)
        {
            printf("la %d solution a viole la contrainte de capacite \n",i);
        }
    }

}
///****************************************************************************************************************************
void calculerMediansNonMediansVertices(partitionPMP* populationPMP)
{

    int i,j,testValueGenotype;
#if mouchard
    printf("calculerMediansNonMediansVertices : ...\n");
#endif // mouchard
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
///****************************************************************************************************************************
void calculcoutCoupeeEtFitnessPMP(partitionPMP* populationPMP)
{
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

///****************************************************************************************************************************
/// Cette fonction assure la s�l�ction des individu pour la pprocedure d'appariement
/// elle �mite la roue de lottrie biais�
void naturalSelectionPMP(partitionPMP* populationPMP1,partitionPMP* populationPMP2)
{
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
///****************************************************************************************************************************
void crossOverPMP(partitionPMP* populationPMP1, partitionPMP* populationPMP2)
{
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
///****************************************************************************************************************************
int findTheBestSolutionPMP(partitionPMP *populationPMP)
{
    float maxFitness = populationPMP->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationPMP+i)->fitness)
        {
            maxFitness = (populationPMP+i)->fitness;
            indice = i;
        }
    }
#if mouchard
    printf("\n maxFitness = %0.4f \t indice = %d \n",maxFitness, indice);
#endif // mouchard
    return indice;
}
///****************************************************************************************************************************
void mutationPMP(partitionPMP* populationPMP)
{

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
///****************************************************************************************************************************
float testerLaSommeDesFitnessPMP(partitionPMP* populationPMP)
{
    printf("testerLaSommeDesFitness ...\n");
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationPMP+i)->fitness;
    }
    return sommeFitness;
}
///****************************************************************************************************************************
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
}

///************************************************************************************************************
void writeSolutionInFilePMP(partitionPMP *populationPMP, FILE *outputFilePop,int iteration)
{

    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///�criture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(populationPMP+i)->id,
                (populationPMP+i)->coutCoupe,(populationPMP+i)->fitness,(populationPMP+i)->coutCoupeNormalise,
                (populationPMP+i)->contrainteViole,(populationPMP+i)->medians);
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationPMP+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationPMP+i)->phenotype[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    /// fprintf(outputFilePop,"\n===============================================\n");
    fprintf(outputFilePop,"\n\n");

}
///************************************************************************************************************
void writeBestSolutionInFilePMP(partitionPMP *solutionDominante, FILE *outputFile,int iteration)
{
    int i;
    fprintf(outputFile,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominante->id,
            solutionDominante->coutCoupe,solutionDominante->fitness,solutionDominante->coutCoupeNormalise,
            solutionDominante->contrainteViole,solutionDominante->medians);
    for(i=0; i<nbrNoeuds; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->genotype[i]);
    }
    fprintf(outputFile,"\t\t\t\t");
    for(i=0; i<nbrNoeuds; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->phenotype[i]);
    }
    fprintf(outputFile,"\n");

}
///**************************************************************************************
void checkContrainstAndFitnessPenalizationPMP(partitionPMP *populationPMP)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        ///(populationPMP+i)->nbrCluster=0; /// le nombre de cluster doit etre initialis� � 1 pour prendre en cons�diration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes viol�es
        (populationPMP+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationPMP+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau r�sulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationPMP+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationPMP+i)->phenotype[k] == j)
                {
                    (populationPMP+i)->clustersSize[j]++;
                }
            }
            /**if ((population+i)->clustersSize[j] !=0)
            {
                (populationPMP+i)->nbrCluster++;
            }*/
        }

        ///if((populationPMP+i)->nbrCluster > max_clusters || (populationPMP+i)->nbrCluster < min_clusters)
        if((populationPMP+i)->medians > max_clusters || (populationPMP+i)->medians < min_clusters)
        {
            (populationPMP+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationPMP+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationPMP+i)->clustersSize[j]!=0)
            {
                if((populationPMP+i)->clustersSize[j]>max_sizeCluster || (populationPMP+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationPMP+i)->constraintVector[1]++;

                }
            }

        }

        if((populationPMP+i)->constraintVector[1] != 0)
        {
            (populationPMP+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstraintePMP(populationPMP);
    }

}
///************************************************************************************************************
void checkCohabitationAndNonCohabitationConstraintePMP(partitionPMP *populationPMP)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
	    /// Cette tache est d�j� r�alis� au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
        (populationPMP+i)->constraintVector[2]=0;
        (populationPMP+i)->constraintVector[3]=0;
		*/
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

///************************************************************************************************************
void pMedianEncodingAffectationParArretes(int nbrGeneration ,FILE *outputFilePMP,FILE *outputFilePopPMP,FILE *outputOptimalSolutionFilePMP,partitionPMP *populationPMP1,partitionPMP *populationPMP2, partitionPMP *solutionDominante)
{

    int bestSolution,iteration=1,bestSolutionIteration,ES=0;
    clock_t t1,t2;
    double temps;

    t1=clock();
    generatePopulationPMP(populationPMP1,0);
    calculerMediansNonMediansVertices(populationPMP1);
    affectationDesNonMediansParArretes(populationPMP1);
    checkContrainstAndFitnessPenalizationPMP(populationPMP1);
    calculcoutCoupeeEtFitnessPMP(populationPMP1);
    ///initialisation de la solution dominante
    bestSolution=findTheBestSolutionPMP(populationPMP1);
    *solutionDominante=*(populationPMP1+bestSolution);
    bestSolutionIteration = 1;
    nbrApparition=1;

    //for(iteration =2 ; iteration <= nbrGeneration; iteration++){
    do{

        naturalSelectionPMP(populationPMP1,populationPMP2);
        crossOverPMP(populationPMP2, populationPMP1);
        mutationPMP(populationPMP1);
        calculerMediansNonMediansVertices(populationPMP1);
        affectationDesNonMediansParArretes(populationPMP1);
        checkContrainstAndFitnessPenalizationPMP(populationPMP1);
        calculcoutCoupeeEtFitnessPMP(populationPMP1);

///=====================================================================================================
        /// dans ce qui suit, on va juste localiser la meilleur solution
///=====================================================================================================

        bestSolution=findTheBestSolutionPMP(populationPMP1);
        if((populationPMP1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise )
                /** le 23/12/2015 : pour avoir la meilleur solution avec le moins de contrainte viol�es*/
                /// && (populationPMP1+bestSolution)->contrainteViole   < solutionDominante->contrainteViole)
        {
            *(solutionDominante) = *(populationPMP1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else if((populationPMP1+bestSolution)->coutCoupeNormalise <= solutionDominante->coutCoupeNormalise )
        {
            nbrApparition++;
        }

        ///iteration++;
    }while(nbrApparition < max_steady_generation);

    t2=clock();
    displayTheBestSolutionPMP(solutionDominante);
    temps = (double)(t2-t1)/CLOCKS_PER_SEC;
    if(bestSolutionIteration>=2){
        ES = ((taillePopulation-elitismeRate)*(bestSolutionIteration-2))+(taillePopulation+ solutionDominante->id -elitismeRate+1);
    }
    else {
        ES = solutionDominante->id +1;
    }
    writeOptimalSolutionInFilePMP(solutionDominante,outputOptimalSolutionFilePMP,nbrRun, bestSolutionIteration, temps,ES);
    printf("le temps d execution est %lf \n", temps);

    fprintf(outputFilePMP,"\n");
    writeBestSolutionInFilePMP(solutionDominante,outputFilePMP,bestSolutionIteration);
}


///************************************************************************************************************
void affectationDesNonMediansParArretes_old(partitionPMP* populationPMP){

    int i,j,l,k,maxDistance,sommeDistances,indiceNotSeed,indiceSeed;
    clock_t t1,t2;
    double temps;
    int distanceArray[nbrNoeuds],nbrNotSeed,tmpNotSeed[nbrNoeuds],nbrFullCluster;



    for(i=0; i<taillePopulation; i++){

//***********************************************************************************************************
        for(l=0; l<(populationPMP+i)->medians; l++){
            (populationPMP+i)->clustersSize[(populationPMP+i)->mediansVertices[l]] = 1;
        }
//***********************************************************************************************************
        for(l=0; l<(populationPMP+i)->nonMedians; l++){
            tmpNotSeed[l] = 1;
        }
//***********************************************************************************************************
		for(j=0;j<nbrNoeuds;j++){
			(populationPMP+i)->phenotype[j] = j;
		}
//***********************************************************************************************************
        nbrNotSeed = (populationPMP+i)->nonMedians;
        nbrFullCluster =0;

        do{
            for(j=0; j < (populationPMP+i)->nonMedians;j++){
                if(tmpNotSeed[j]>0){
                    indiceNotSeed = (populationPMP+i)->nonMediansVertices[j];
                    for(k=0; k<(populationPMP+i)->medians; k++){
                        indiceSeed = (populationPMP+i)->mediansVertices[k];
                        distanceArray[k] = 0;
                        if( (populationPMP+i)->clustersSize[indiceSeed] < max_sizeCluster){
                            for(l=0; l<nbrNoeuds; l++){
                                if((populationPMP+i)->phenotype[l] == indiceSeed){
                                    if(fluxMatrix[indiceNotSeed][l] >= 0 ){
                                        distanceArray[k] += fluxMatrix[indiceNotSeed][l];
                                    }
                                }
                            }
                        }
                        else distanceArray[k] = -1;
                    }
                    maxDistance = distanceArray[0];
                    indiceSeed=(populationPMP+i)->mediansVertices[0];
                    for(k=1; k<(populationPMP+i)->medians; k++){
                        if(maxDistance < distanceArray[k] ) { maxDistance = distanceArray[k]; indiceSeed = (populationPMP+i)->mediansVertices[k]; }
                    }
                    if(maxDistance != -1){
                        (populationPMP+i)->phenotype[indiceNotSeed] = indiceSeed;
                        (populationPMP+i)->clustersSize[indiceSeed]++;
                        if((populationPMP+i)->clustersSize[indiceSeed] >= max_sizeCluster)nbrFullCluster++;
                        tmpNotSeed[j] = -1;
                        nbrNotSeed--;
                    }
                }
            }
            if(nbrFullCluster ==(populationPMP+i)->medians ) nbrNotSeed = 0;
        }while(nbrNotSeed > 0);
 	}
}


///************************************************************************************************************
void pMedianEncodingBaseSurClusterSize(int nbrGeneration ,FILE *outputFilePMP,FILE *outputFilePopPMP,FILE *outputOptimalSolutionFilePMP,partitionPMP *populationPMP1,partitionPMP *populationPMP2, partitionPMP *solutionDominante)
{

    int bestSolution,iteration=1,bestSolutionIteration,ES=0;
    clock_t t1,t2;
    double temps;

    t1=clock();

    generatePopulationPMP(populationPMP1,0);
    calculerMediansNonMediansVertices(populationPMP1);
    affectationDesNonMediansBaseSurClusterSize(populationPMP1);
    checkContrainstAndFitnessPenalizationPMP(populationPMP1);
    calculcoutCoupeeEtFitnessPMP(populationPMP1);
    ///initialisation de la solution dominante
    bestSolution=findTheBestSolutionPMP(populationPMP1);
    *solutionDominante=*(populationPMP1+bestSolution);
    bestSolutionIteration = 1;
    nbrApparition=1;

    //for(iteration =2 ; iteration <= nbrGeneration; iteration++){
   do{

        naturalSelectionPMP(populationPMP1,populationPMP2);
        crossOverPMP(populationPMP2, populationPMP1);
        mutationPMP(populationPMP1);
        calculerMediansNonMediansVertices(populationPMP1);
        affectationDesNonMediansBaseSurClusterSize(populationPMP1);
        checkContrainstAndFitnessPenalizationPMP(populationPMP1);
        calculcoutCoupeeEtFitnessPMP(populationPMP1);

        bestSolution=findTheBestSolutionPMP(populationPMP1);
        if((populationPMP1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise )
                /** le 23/12/2015 : pour avoir la meilleur solution avec le moins de contrainte viol�es*/
                /// && (populationPMP1+bestSolution)->contrainteViole   < solutionDominante->contrainteViole)
        {
            *(solutionDominante) = *(populationPMP1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else if((populationPMP1+bestSolution)->coutCoupeNormalise <= solutionDominante->coutCoupeNormalise )
        {
            nbrApparition++;
        }

        ///iteration++;
    }while(nbrApparition <max_steady_generation);


    t2=clock();
/**
    ///il faut modifier le cout de coupe de la solutionDominante (car j'ai rajouter maxFlow aux couple de coh )
    /// j'ai la somme r�elle des flux + les flux intra r�elle de la solution dominante ==> flux r�elle inter � calculer

    int i,nd, na,edgeNumber;
    for(i=0;i <nbrCohabitationConstraintes; i++){
        nd = cohabitationConstraintes[i].nouedDepart;
        na = cohabitationConstraintes[i].nouedArrive;
        edgeNumber = cohabitationConstraintes[i].numeroEdge;
        if(solutionDominante->phenotype[nd] == solutionDominante->phenotype[nd]){
            solutionDominante->coutCoupe = solutionDominante->coutCoupe - maxFlow;
        }
    }
*/
    displayTheBestSolutionPMP(solutionDominante);

    temps = (double)(t2-t1)/CLOCKS_PER_SEC;
    if(bestSolutionIteration>=2){
        ES = ((taillePopulation-elitismeRate)*(bestSolutionIteration-2))+(taillePopulation+ solutionDominante->id -elitismeRate+1);
    }
    else {
        ES = solutionDominante->id +1;
    }
    writeOptimalSolutionInFilePMP(solutionDominante,outputOptimalSolutionFilePMP,nbrRun, bestSolutionIteration, temps,ES);
    printf("le temps d execution est %lf \n", temps);

    fprintf(outputFilePMP,"\n\n");
    writeBestSolutionInFilePMP(solutionDominante,outputFilePMP,bestSolutionIteration);

}

///****************************************************************************************************************************
/// cette fonction affecte les noeuds non medians au noeuds medians, permettant la cr�ation des
/// PHENOTYPE DES SOLUTION CE QUI PERMET PAR LA SUITE LE CALCULE DES FITNESS
///void affectationDesNonMediansBaseSurLesArretes(partition* population,int taillePopulation, int nbrNoeuds, int fluxMatrix[nbrNoeuds][nbrNoeuds],int max_sizeCluster){
void affectationDesNonMediansBaseSurClusterSize(partitionPMP* populationPMP)
{

    int i,j,l,k,m,maxDistance,indiceNotSeed,indiceSeed,nbrSommets,nonCohabit;
    for(i=0; i<taillePopulation; i++)
    {
///***************************************************************************************************
        for(l=0; l<nbrNoeuds; l++)
        {
            (populationPMP+i)->phenotype[l] = l;
        }

        for(l=0; l<(populationPMP+i)->medians; l++)
        {
            (populationPMP+i)->clustersSize[(populationPMP+i)->mediansVertices[l]] = 1;
        }
///***************************************************************************************************

        if((populationPMP+i)->nonMedians > 0 && (populationPMP+i)->medians > 0) // sometimes we have some solutions that contains no medians or not medians vertices
        {
            for(j=0; j<(populationPMP+i)->nonMedians;j++)
            {
                maxDistance = -1;
                indiceNotSeed= (populationPMP+i)->nonMediansVertices[j];

                for(k=0; k<(populationPMP+i)->medians; k++)
                {
                    indiceSeed = (populationPMP+i)->mediansVertices[k];
                    nonCohabit = 0;
                    for(l=0;l<nbrNonCohabitationConstraintes;l++){
                        if((nonCohabitationConstraintes[l].nouedDepart == indiceSeed && nonCohabitationConstraintes[l].nouedArrive == indiceNotSeed)||
                           (nonCohabitationConstraintes[l].nouedDepart == indiceNotSeed && nonCohabitationConstraintes[l].nouedArrive == indiceSeed))
                        {nonCohabit = 1; break;}

                    }

                    if( fluxMatrix[indiceNotSeed][indiceSeed] >= 0
					&&
					(populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]] < max_sizeCluster
                    &&
                    !nonCohabit ) /// nonCohabit = 0 ==> les deux peuvent cohabiter
                    {
                        if(maxDistance == -1)
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[indiceSeed];

                        }
                        else  if(maxDistance < fluxMatrix[indiceNotSeed][indiceSeed])
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[indiceSeed];
                        }
                    }
                }
                if(maxDistance >= 0)
                {
                    indiceSeed = (populationPMP+i)->phenotype[indiceNotSeed] ;
                    (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]]++;
                }
                else
                {
                    for(k=0; k<(populationPMP+i)->medians; k++)
                    {
                        indiceSeed = (populationPMP+i)->mediansVertices[k];
                        if((populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceSeed]] < max_sizeCluster && !nonCohabit)
                        {
                            maxDistance = fluxMatrix[indiceNotSeed][indiceSeed];
                            (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[indiceSeed];
                            (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceNotSeed]]++;
                            break;
                        }
                    }
                    if(maxDistance == -1)
                    {
                        (populationPMP+i)->phenotype[indiceNotSeed] = (populationPMP+i)->phenotype[(populationPMP+i)->mediansVertices[(populationPMP+i)->medians-1]];
                        (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[indiceNotSeed]]++;
                    }
                }
            }
        }
    }

}


///**********************************************************************************************************
void writeOptimalSolutionInFilePMP(partitionPMP *solutionDominante,FILE* outputOptimalSolutionFilePMP,
                                   int nbrRun, int bestSolutionIteration, float runTime, int ES)
{
    int i;
    fprintf(outputOptimalSolutionFilePMP,"RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d |",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);


    fprintf(outputOptimalSolutionFilePMP," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",
            solutionDominante->constraintVector[0],solutionDominante->constraintVector[1],solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);



}
///**********************************************************************************************************

int compareCroissantFitnessPMP(void const *a, void const *b)
{

    partitionPMP const *pa = a;
    partitionPMP const *pb = b;

    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;;
}



///*********************************************************************************************************
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
///*************************************************************************************************************
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
///********************************************************************************************
/// J'ai remarqu� que m�me si les solutions sont r�alisables le nombre de clusters reste tout
/// de m�me important du coup j'ai d�cid� de contracter les clusters autant que possible
/// en respectant les contraintes de taille, et de non cohabitation.
///********************************************************************************************
void contractionDesClusters(partitionPMP *solutionDominante){


    printf("\ncontractionDesClusters ....\n");
    printf("cout de coupe avant contraction = %d\n",solutionDominante->coutCoupe);
    int i,j,k,nd,na,clustersNonContractable,C[1000];
    calculerLesClustersDeLaSolution(C,solutionDominante);
    ///**********************************************************************************

    for(i=0; i< solutionDominante->medians-1; i++){
        for(j=i+1; j<solutionDominante->medians; j++){
            if(solutionDominante->clustersSize[C[i]] + solutionDominante->clustersSize[C[j]] <= max_sizeCluster ){
                /// v�rification des contraintes de non cohabitation
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
                    solutionDominante->clustersSize[C[i]] += solutionDominante->clustersSize[C[j]];
                    calculerLesClustersDeLaSolution(C,solutionDominante);
                }
            }
        }
    }
}


///******************************************************************************************

int realCutSizeForTheBestSolution(partitionPMP *solutionDominante){

    int i,j,sommeTotalCoutCoupe = 0;

        solutionDominante->coutCoupe = 0; /// r�initialisation de la valriable apr�s chaque it�ration (individus de la population)
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
///*******************************************************************************************

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


///************************************************************************************************************
void affectationDesNonMediansParArretes(partitionPMP* populationPMP){
    int i,j,l,k,maxDistance,sommeDistances,indiceNotSeed,indiceSeed,assignedOpurtunity,nbrVoisins,voisin;
    clock_t t1,t2;
    double temps;
    int distanceArray[nbrNoeuds],nbrNotSeed,tmpNotSeed[nbrNoeuds],nbrFullCluster;

    for(i=0; i<taillePopulation; i++){
//***********************************************************************************************************
		for(j=0;j<nbrNoeuds;j++){
			(populationPMP+i)->phenotype[j] = -1;
		}
//***********************************************************************************************************
        for(l=0; l<(populationPMP+i)->medians; l++){
            (populationPMP+i)->clustersSize[(populationPMP+i)->mediansVertices[l]] = 1;
            (populationPMP+i)->phenotype[(populationPMP+i)->mediansVertices[l]] = (populationPMP+i)->mediansVertices[l];
        }
//***********************************************************************************************************
        for(l=0; l<(populationPMP+i)->nonMedians; l++){ /// the non affected vertices (non medians)
            tmpNotSeed[l] = 1;
        }
//***********************************************************************************************************
        nbrNotSeed = (populationPMP+i)->nonMedians;
        nbrFullCluster =0;

        do{
            for(j=0; j < (populationPMP+i)->nonMedians ;j++){
                if(tmpNotSeed[j]>0){
                    indiceNotSeed = (populationPMP+i)->nonMediansVertices[j];
                   /// set to -1 the distanceArray vector
                    for(k=0;k<(populationPMP+i)->medians;k++) distanceArray[(populationPMP+i)->mediansVertices[k]] = -1;
                    assignedOpurtunity = 0;
                    nbrVoisins = neighborMatrix[indiceNotSeed][0];

                    for(l=1; l<=nbrVoisins ; l++){ /// we must only check the neighborhood of the indiceNotSeed
                        voisin = neighborMatrix[indiceNotSeed][l];
                        if((populationPMP+i)->phenotype[voisin] != -1 /* le voisin est d�j� affect� */ &&
                           (populationPMP+i)->clustersSize[(populationPMP+i)->phenotype[voisin]] < max_sizeCluster){

                                if(distanceArray[(populationPMP+i)->phenotype[voisin]] == -1)
                                    distanceArray[(populationPMP+i)->phenotype[voisin]] = fluxMatrix[indiceNotSeed][voisin];
                                else
                                    distanceArray[(populationPMP+i)->phenotype[voisin]] += fluxMatrix[indiceNotSeed][voisin];

                                assignedOpurtunity = 1;

                        }
                    }
                    if(assignedOpurtunity){
                        maxDistance = -1;
                        for(k=0; k<(populationPMP+i)->medians ; k++){
                            if(distanceArray[(populationPMP+i)->mediansVertices[k]] != -1 &&
                                        maxDistance < distanceArray[(populationPMP+i)->mediansVertices[k]]){
                                maxDistance = distanceArray[(populationPMP+i)->mediansVertices[k]];
                                indiceSeed = (populationPMP+i)->mediansVertices[k];
                            }
                        }

                        if(maxDistance != -1){
                            (populationPMP+i)->phenotype[indiceNotSeed] = indiceSeed;
                            (populationPMP+i)->clustersSize[indiceSeed]++;
                            if((populationPMP+i)->clustersSize[indiceSeed] >= max_sizeCluster)nbrFullCluster++;
                            tmpNotSeed[j] = -1;
                            nbrNotSeed--;
                        }
                  }
                }
            }
            if(j==(populationPMP+i)->nonMedians && nbrNotSeed >0){
                for(k=0;k<nbrNoeuds;k++){
                    if((populationPMP+i)->phenotype[k]== -1){
                            for(l=0;l<(populationPMP+i)->medians;l++){
                                if((populationPMP+i)->clustersSize[(populationPMP+i)->mediansVertices[l]]< max_sizeCluster){
                                    (populationPMP+i)->phenotype[k] = (populationPMP+i)->mediansVertices[l];
                                    (populationPMP+i)->clustersSize[(populationPMP+i)->mediansVertices[l]]++;
                                    if((populationPMP+i)->clustersSize[(populationPMP+i)->mediansVertices[l]] >= max_sizeCluster)nbrFullCluster++;
                                    tmpNotSeed[j] = -1;
                                    nbrNotSeed--;
                                    break;
                                }
                            }
                    }
                }
            }
            if(nbrFullCluster ==(populationPMP+i)->medians ){
                for(j=0; j< nbrNoeuds; j++){
                    if((populationPMP+i)->phenotype[j] == -1 )
                        (populationPMP+i)->phenotype[j] = j;//(populationPMP+i)->mediansVertices[(populationPMP+i)->medians - 1];
                        //(populationPMP+i)->phenotype[j] = (populationPMP+i)->mediansVertices[(populationPMP+i)->medians - 1];
                }
                nbrNotSeed = 0;
            }
        }while(nbrNotSeed > 0);
 	}
}
///**********************************************************************************************




