#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../headers/typeDeclaration.h"
#include "../headers/communesFunctions.h"
#include "../headers/genetic_alg_SVTC.h"
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
        nbrRun;

extern float MBF, ART, AES, SDBF, BF[40];

extern edge *edgeVector,*fixedEdgeVector,cohabitationConstraintes[100],nonCohabitationConstraintes[100];
///***********************************************
extern partitionSVTC *populationSVTC1 = NULL,*populationSVTC2 = NULL, *solutionDominanteSVTC = NULL,*bestSolutionOverRunsSVTC = NULL;

///**************************************************************************************
///g�n�ration de la population initiale
void generatePopulationSVTC(partitionSVTC* populationSVTC, int indiceFirstElt)
{

    int i,j;
    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        for(j=0; j<nbrNoeuds; j++)
        {
            /// on va essayer de limiter le nombre de clusters pour voir ce que �a va donner
            ///(populationSVTC+i)->genotype[j] = rnd(0,nbrParties);
            (populationSVTC+i)->genotype[j] = rnd(0,nbrNoeuds); /// [0;nbrNoeuds-1]
        }
        (populationSVTC+i)->id =i;
    }
}
///**********************************************************************************
void agencementDesGenotype(partitionSVTC* populationSVTC)
{
    int i,j,k,m=nbrNoeuds, t1[m],t[m],sortIndex,tmp,exist,min;
    for(i=0; i<taillePopulation; i++){
        for(j=0; j<nbrNoeuds; j++)
        {
            t1[j] = -1;
        }

        t[0]=(populationSVTC+i)->genotype[0]; sortIndex=1;
        for(j=1;j<nbrNoeuds;j++){
            exist = 0;
            for(k=0;k<sortIndex;k++){
                if((populationSVTC+i)->genotype[j]==t[k]){
                    exist = 1;
                    break;
                }
            }
            if(!exist){
                t[sortIndex]=(populationSVTC+i)->genotype[j];
                sortIndex++;
            }
        }

        for(j=0; j<sortIndex; j++)
        {
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationSVTC+i)->genotype[k] == t[j] && t1[k] <0)
                {
                    (populationSVTC+i)->genotype[k] = j;
                    t1[k] = 1;
                }
            }
        }
    }
}


///**************************************************************************************
/// affichage des population
void affichePopulationSVTC(partitionSVTC* populationSVTC)
{
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        printf("\n id = %d \t cout de coupe = %d \t fitness = %0.4f \t",
               (populationSVTC+i)->id,(populationSVTC+i)->coutCoupe, (populationSVTC+i)->fitness);
        printf("le cout de coupe normalise = %d \n",(populationSVTC+i)->coutCoupeNormalise);
        printf("la partition est  = \t");
        for(j=0; j<nbrNoeuds; j++)
        {
            printf("%d ",(populationSVTC+i)->genotype[j]);
        }
        printf("\n");
        printf("le nombre des clusters est %d \n",(populationSVTC+i)->nbrCluster);
        printf("le nombre des contraintes violees est  = %d \n",(populationSVTC+i)->contrainteViole);
        printf("affichage des tailles des clusters \n");
        for(j=0; j<(populationSVTC+i)->nbrCluster; j++)
        {
            printf("%d ",(populationSVTC+i)->clustersSize[j]);
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
void calculCoutCoupeEtFitnessSVTC(partitionSVTC* populationSVTC)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationSVTC+i)->coutCoupe = 0; /// r�initialisation de la valriable apr�s chaque it�ration (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                /// on peut avoir deux noeuds du m�me cluster mais qui sont par relier par une arr�tes
                if((populationSVTC+i)->genotype[j] == (populationSVTC+i)->genotype[k] && fluxMatrix[j][k] >= 0)
                {
                    (populationSVTC+i)->coutCoupe = (populationSVTC+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationSVTC+i)->coutCoupeNormalise = (populationSVTC+i)->coutCoupe + ((nbr_constraint - (populationSVTC+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationSVTC+i)->coutCoupeNormalise;
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
            varianceCoutDeCoupeNormalise + pow(((populationSVTC+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajout�
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationSVTC+i)->expectedValue = (populationSVTC+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationSVTC+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationSVTC+i)->fitness = (float)((populationSVTC+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/SortedVTC/expectedValueSVTC.txt","w");
    for(i=0; i<taillePopulation ; i++)
    {
/**
        if(iteration == 1)
        {
            fprintf(file,"%d \t %d \t %0.2f \n",i,(populationSVTC+i)->coutCoupeNormalise,(populationSVTC+i)->expectedValue);
        }
*/
        (populationSVTC+i)->fitness = (float)((populationSVTC+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling
}
///**************************************************************************************
void naturalSelectionSVTC(partitionSVTC* populationSVTC1,partitionSVTC* populationSVTC2)
{
    int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<elitismeRate; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationSVTC1+maxFitness)->coutCoupeNormalise < (populationSVTC1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationSVTC2+i) = *(populationSVTC1+maxFitness);
        ///printf("(populationSVTC2+%d)->coutCoupeNormalise = %d \n",i,(populationSVTC2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationSVTC2, taillePopulation, sizeof *populationSVTC2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationSVTC2,tmpPopulation);

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
            sommeFitness = sommeFitness + (populationSVTC1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationSVTC2+i) = *(populationSVTC1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationSVTC2+i) = *(populationSVTC1+j);
            }
        }
    }
    ///free(tmpPopulation);

}
///**************************************************************************************
void crossOverSVTC(partitionSVTC* populationSVTC1, partitionSVTC* populationSVTC2)
{
    int i = 0,j =0;
    int choixInd1, choixInd2, choixLocus;
    ///***********************************************************************************
    for(i=0; i<elitismeRate; i++)
    {
        *(populationSVTC2+i) = *(populationSVTC1+i);
        (populationSVTC2+i)->id = i;
    }
    ///***********************************************************************************
    while (i < taillePopulation-regeneration)
    {
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );
        choixLocus = rnd(0,nbrNoeuds-1); /// l'intervalle sera entre 1 et nbrNoeuds-2 (la borne superieur -1)
        for(j=0; j<=choixLocus; j++)
        {
            (populationSVTC2+i)->genotype[j] = (populationSVTC1+choixInd1)->genotype[j];
            (populationSVTC2+i+1)->genotype[j] = (populationSVTC1+choixInd2)->genotype[j];
        }
        ///**********************************************************************
        for(j=choixLocus+1; j<nbrNoeuds; j++)
        {
            (populationSVTC2+i)->genotype[j] = (populationSVTC1+choixInd2)->genotype[j];
            (populationSVTC2+i+1)->genotype[j] = (populationSVTC1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationSVTC2+i)->id = i;
        (populationSVTC2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }
    generatePopulationSVTC(populationSVTC2, taillePopulation-regeneration);

}
///**************************************************************************************
int findTheBestSolutionSVTC(partitionSVTC *populationSVTC)
{
    /// la fitness est calcul� � partir des couts de coupe normalis�s
    float maxFitness = populationSVTC->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationSVTC+i)->fitness)
        {
            maxFitness = (populationSVTC+i)->fitness;
            indice = i;
        }
    }
    return indice;
}
///**************************************************************************************
void mutationSVTC(partitionSVTC* populationSVTC)
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
            numeroGeneMute = rnd(0,nbrParties);
            (populationSVTC+i)->genotype[numeroGeneMute] = rnd(0,nbrParties);

        }
    }

}
///**************************************************************************************
float testerLaSommeDesFitnessSVTC(partitionSVTC* populationSVTC)
{
#if mouchard
    printf("testerLaSommeDesFitness ...\n");
#endif // mouchard
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationSVTC+i)->fitness;
    }
    return sommeFitness;
}
///************************************************************************************************************
void displayTheBestSolutionSVTC(partitionSVTC* solutionDominante)
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

///**********************************************************************************
void writeSolutionInFileSVTC(partitionSVTC *populationSVTC, FILE *outputFilePop, int iteration)
{

    int i,j;
    for(i=0; i<taillePopulation; i++)
    {

        ///�criture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(populationSVTC+i)->id,
                (populationSVTC+i)->coutCoupe,(populationSVTC+i)->fitness,(populationSVTC+i)->coutCoupeNormalise,
                (populationSVTC+i)->contrainteViole,(populationSVTC+i)->nbrCluster);
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationSVTC+i)->genotype[j]);
        }
        fprintf(outputFilePop,"\t\t");

        for(j=0; j<nbr_constraint; j++)
        {
            fprintf(outputFilePop,"%d ",(populationSVTC+i)->constraintVector[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    fprintf(outputFilePop,"\n\n");

}
///**************************************************************************************
void writeBestSolutionInFileSVTC(partitionSVTC *solutionDominante, FILE *outputFile,int iteration)
{

    int i;
    fprintf(outputFile,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominante->id,
            solutionDominante->coutCoupe,solutionDominante->fitness,solutionDominante->coutCoupeNormalise,
            solutionDominante->contrainteViole,solutionDominante->nbrCluster);
    for(i=0; i<nbrNoeuds; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->genotype[i]);
    }
    fprintf(outputFile,"\n");
}
///======================================================================================
void sortedVertexToClusterEncoding(int nbrGeneration ,FILE *outputFileSVTC,FILE *outputFilePopSVTC,FILE *outputOptimalSolutionFileSVTC, partitionSVTC *populationSVTC1,partitionSVTC *populationSVTC2,partitionSVTC *solutionDominante)
{

    int bestSolution,iteration=1 ,bestSolutionIteration,ES = 0;
    clock_t t1,t2;
    double temps;

    t1=clock();


    generatePopulationSVTC(populationSVTC1,0);
    agencementDesGenotype(populationSVTC1);
    checkContrainstAndFitnessPenalizationSVTC(populationSVTC1);
    calculCoutCoupeEtFitnessSVTC(populationSVTC1);
    bestSolution=findTheBestSolutionSVTC(populationSVTC1);
    *solutionDominante=*(populationSVTC1+bestSolution);
    bestSolutionIteration= 1;
    nbrApparition=1;
    //for(iteration =2 ; iteration <= nbrGeneration; iteration++){
    do{
///=====================================================================================================
        /// application des trois op�rateurs de bases de g�n�tiques
        naturalSelectionSVTC(populationSVTC1,populationSVTC2);
        crossOverSVTC(populationSVTC2,populationSVTC1); /// nbrNoeuds repr�sente la taille des genotype
        mutationSVTC(populationSVTC1);
        agencementDesGenotype(populationSVTC1);
        checkContrainstAndFitnessPenalizationSVTC(populationSVTC1);
        calculCoutCoupeEtFitnessSVTC(populationSVTC1);
        bestSolution=findTheBestSolutionSVTC(populationSVTC1);
        if((populationSVTC1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise )
                /** le 23/12/2015 : pour avoir la meilleur solution avec le moins de contrainte viol�es*/
                /// && (populationSVTC1+bestSolution)->contrainteViole   < solutionDominante->contrainteViole)
        {
            *(solutionDominante) = *(populationSVTC1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else if((populationSVTC1+bestSolution)->coutCoupeNormalise <= solutionDominante->coutCoupeNormalise )
        {
            nbrApparition++;
        }
        ///iteration++;
    }while(nbrApparition <max_steady_generation);

    displayTheBestSolutionSVTC(solutionDominante);
    t2=clock();
    temps = (double)(t2-t1)/CLOCKS_PER_SEC;
    printf("le temps d execution est %lf \n", temps);

    if(bestSolutionIteration>=2){
        ES = ((taillePopulation-elitismeRate)*(bestSolutionIteration-2))+(taillePopulation+ solutionDominante->id -elitismeRate+1);
    }
    else {
        ES = solutionDominante->id +1;
    }

    writeOptimalSolutionInFileSVTC(solutionDominante,outputOptimalSolutionFileSVTC,nbrRun,bestSolutionIteration,temps, ES);

    fprintf(outputFileSVTC,"\n");
    writeBestSolutionInFileSVTC(solutionDominante,outputFileSVTC,bestSolutionIteration);

}

///**************************************************************************************
void checkContrainstAndFitnessPenalizationSVTC(partitionSVTC *populationSVTC)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        (populationSVTC+i)->nbrCluster=0; /// le nombre de cluster doit etre initialis� � 1 pour prendre en cons�diration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes viol�es
        (populationSVTC+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationSVTC+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau r�sulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationSVTC+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationSVTC+i)->genotype[k] == j)
                {
                    (populationSVTC+i)->clustersSize[j]++;
                }
            }
            if ((populationSVTC+i)->clustersSize[j] !=0)
            {
                (populationSVTC+i)->nbrCluster++;
            }
        }

        if((populationSVTC+i)->nbrCluster > max_clusters || (populationSVTC+i)->nbrCluster < min_clusters)
        {
            (populationSVTC+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationSVTC+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationSVTC+i)->clustersSize[j]!=0)
            {
                if((populationSVTC+i)->clustersSize[j]>max_sizeCluster || (populationSVTC+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationSVTC+i)->constraintVector[1]++;

                }
            }

        }

        if((populationSVTC+i)->constraintVector[1] != 0)
        {
            (populationSVTC+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstrainteSVTC(populationSVTC);
    }

}
///**************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteSVTC(partitionSVTC *populationSVTC)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
	    /// Cette tache est d�j� r�alis� au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
        (populationSVTC+i)->constraintVector[2]=0;
        (populationSVTC+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationSVTC+i)->genotype[noeud1]!= (populationSVTC+i)->genotype[noeud2])
                {
                    (populationSVTC+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes viol�es
            if((populationSVTC+i)->constraintVector[2]!=0)
            {
                (populationSVTC+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;
                if((populationSVTC+i)->genotype[noeud1]== (populationSVTC+i)->genotype[noeud2])
                {
                    (populationSVTC+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes viol�es
            if((populationSVTC+i)->constraintVector[3]!=0)
            {
                (populationSVTC+i)->contrainteViole++;
            }
        }


    }

}

///********************************************************************************************
void writeOptimalSolutionInFileSVTC(partitionSVTC *solutionDominante,FILE* outputOptimalSolutionFileSVTC,
                                  int nbrRun, int bestSolutionIteration, float runTime, int ES){

    fprintf(outputOptimalSolutionFileSVTC,"RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d | ",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);
fprintf(outputOptimalSolutionFileSVTC," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",
            solutionDominante->constraintVector[0],solutionDominante->constraintVector[1],solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);

}
///******************************************************************************************
int compareCroissantFitnessSVTC (void const *a, void const *b)
{

    partitionSVTC const *pa = a;
    partitionSVTC const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}
