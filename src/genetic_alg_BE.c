//
// Created by Ali on 01-Apr-24.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../headers/typeDeclaration.h"
#include "../headers/communesFunctions.h"
#include "../headers/genetic_alg_BE.h"
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
extern partitionBE *populationBE1,*populationBE2, *solutionDominanteBE,*bestSolutionOverRunsBE;
extern partitionDVTC *solutionDominanteDVTC;
///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
/// DECLARATION DES VARIABLES FILES POUR NOS FICHIERS
///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
extern char cheminBestSolutionBE[150],  cheminOptimalSolutionBE[150], cheminFichierDetailsProbleme[150];
extern FILE *bestSolutionsOverIterationBE, *optimalSolutionFileBE, *FichierDetailsProbleme;


///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
///FFFFFFFFFFF generatePopulationBE FFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
void generatePopulationBE(partitionBE* populationBE, int indexFirstIndividual)
{
    // TODO (Ali#): DONT FORGET TO AJUST THE P% OF INTRACLUSTER EDGES IN generatePopulationBE
    ///QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
    ///WHAT IS THE RELATIONSHIP BETWEEN P1 (nbrZerEdges and the whole length of the genotype) ?
    ///QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ

    /// IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    /// WE USE indexFirstIndividual TO SPECIFY FROM WHICH INDVIDUAL
    /// WE REPRODUCE NEW INDIVIDUALS TO PRESERVE DIVERSSIFICATION WIHTIN THE POP
    /// IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

    int allele,nbrZeroArretes, nbrOneArretes;

    nbrZeroArretes = nbrNoeuds / max_sizeCluster;
    if(nbrZeroArretes != (double)nbrNoeuds / max_sizeCluster)
        nbrZeroArretes++;

    /// FIX ALL THE EDGES OF CC TO 0 AND THOSE OF NCC TO 1
    for(int i=indexFirstIndividual; i<taillePopulation; i++)
    {
        /// THE SOLUTION ID
        (populationBE+i)->id =i;
        /// GENOTYPE INITIALIZATION
        for(int j=0; j<nbrArretes; j++) {
            (populationBE+i)->genotype[j] = 1;
        }
        /// FIX THE VALUES OF CC EDGES
        for(int j=0; j<nbrCohabitationConstraintes; j++)
        {
            for(int k=0; k<nbrArretes; k++)
            {
                if((cohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedDepart &&
                    cohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedArrive) ||
                   (cohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedArrive &&
                    cohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedDepart))
                {
                    (populationBE+i)->genotype[k] = 0;
                }
            }
        }

        /// FIX THE VALUES OF NCC EDGES
        for(int j=0; j<nbrNonCohabitationConstraintes; j++)
        {
            for(int k=0; k<nbrArretes; k++)
            {
                if((nonCohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedDepart &&
                    nonCohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedArrive) ||
                   (nonCohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedArrive &&
                    nonCohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedDepart) )
                {
                    (populationBE+i)->genotype[k] = 1;
                }
            }
        }

        for(int j=0; j<nbrZeroArretes; j++){
            do{
                allele = rnd(0,nbrArretes);
            }
            while((populationBE+i)->genotype[allele]==0);
            (populationBE+i)->genotype[allele] = 0;
        }

    }
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Display the whole population to track problems
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void affichePopulationBE(partitionBE* populationBE)
{
    FILE *pop;
    pop = fopen("./pop.txt","w");

    ///fprintf(pop,"affichePopulation ...\n");
    int i,j;
    for(i=0; i<taillePopulation; i++)
    {
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(pop,"%d ",(populationBE+i)->phenotype[j]);
        }
        fprintf(pop,"\n le nombre de contrainte viole = %d \n",(populationBE+i)->contrainteViole);

        fprintf(pop,"\n\n");
    }
}
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Natural selection process
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void naturalSelectionBE(partitionBE* populationBE1,partitionBE* populationBE2)
{
    float lotterie,sommeFitness;

    for(int i=0; i<taillePopulation; i++)
    {
        *(populationBE2+i) = *(populationBE1+i);
    }

    qsort (populationBE2, taillePopulation, sizeof *populationBE2, compareCroissantFitnessBE);
    for(int i=elitismeRate; i<taillePopulation-reproductionRate; i++)
    {
        lotterie= drand48ForWindows(0,101);
        sommeFitness = 0;
        for(int j=elitismeRate; j<taillePopulation-reproductionRate; j++)
        {
            sommeFitness = sommeFitness + (populationBE1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationBE2+i) = *(populationBE1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationBE2+i) = *(populationBE1+j);
            }
        }
    }
    generatePopulationBE(populationBE2, taillePopulation-reproductionRate-1);

}
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Crossover operator process
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void crossOverBE(partitionBE* populationBE1, partitionBE* populationBE2)
{

    int i,j;
    int choixInd1, choixInd2, choixLocus;
    ///***********************************************************************************
    for(i=0; i<elitismeRate; i++)
    {
        *(populationBE2+i) = *(populationBE1+i);
        (populationBE2+i)->id = i;
    }
    ///***********************************************************************************

    while (i < taillePopulation)
    {
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            /// ce traitement évite l'appariement d'un element avec lui même/
            choixInd2 = rnd(0,taillePopulation);
            ///choixInd2 = rnd(tauxReproduction,taillePopulation);
        }
        while(choixInd2 == choixInd1 );

        choixLocus = rnd(0,nbrArretes-2); /// -2 parce qu'on a commencer à partir du 0

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationBE2+i)->genotype[j] = (populationBE1+choixInd1)->genotype[j];
            (populationBE2+i+1)->genotype[j] = (populationBE1+choixInd2)->genotype[j];
        }
        for(j=choixLocus+1; j<nbrArretes; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la dernière valeure
        {
            (populationBE2+i)->genotype[j] = (populationBE1+choixInd2)->genotype[j];
            (populationBE2+i+1)->genotype[j] = (populationBE1+choixInd1)->genotype[j];
        }
        /// affectation de nouveau indice
        (populationBE2+i)->id = i;
        (populationBE2+i+1)->id = i+1;
        ///*****************************************************************
        i+=2;
    }

}
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Find the best solution among the population
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int findTheBestSolutionBE(partitionBE *populationBE)
{
    int i,indice=0,maxFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        /**
            should I add nbrCluster as a condition to select the bes solution
        */
        if(maxFitness < (populationBE+i)->coutCoupeNormalise)
        {
            maxFitness = (populationBE+i)->coutCoupeNormalise;
            indice = i;
        }
    }
    return indice;

}
///************************************************************************************************************
void mutationBE(partitionBE* populationBE)
{
    int i,numeroGeneMute;   ///int nbrAlleleMute,j;
    double applicationOfMutation;

    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    /// We dont apply mutation on elitists individuals
    /// 20/12/2022
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    for(i=elitismeRate; i<taillePopulation; i++)
    {
#if Windows
        applicationOfMutation = drand48ForWindows(0,101);
        ///nbrAlleleMute = rnd(1,nbrArretes/4);
#else
        applicationOfMutation = drand48();
#endif
        if(applicationOfMutation <= tauxMutation)
        {
            numeroGeneMute = rnd(0,nbrArretes);
            if((populationBE+i)->genotype[numeroGeneMute] == 1)
            {
                (populationBE+i)->genotype[numeroGeneMute] =0;
            }
            else
            {
                (populationBE+i)->genotype[numeroGeneMute] =1;
            }
        }
    }
}
///************************************************************************************************************
float testerLaSommeDesFitnessBE(partitionBE* populationBE)
{
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationBE+i)->fitness;
    }
    return sommeFitness;
}
///************************************************************************************************************
void displayTheBestSolutionBE(partitionBE* solutionDominante)
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
                printf("Cette contrainte n'est pas prise en charge par le système \n");

        }

    }
}
///************************************************************************************************************
void getPartitionFromSolutionBE_old(partitionBE *populationBE)
{
    ///CHECKCHECKCHECKCHECKCHECKCHECKCCHECKCHECKCHECKCHECKCHECKCHECKC
    FILE *phenotype;
    phenotype = fopen("E:/referencedGraphs/ThBE/phenotype.txt","a");
    ///CHECKCHECKCHECKCHECKCHECKCHECKCCHECKCHECKCHECKCHECKCHECKCHECKC
    int i,j,k,nd,na,numeroClustre; /// nd = noued de départ , na = noued d'arrivé
    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation des phénotype = partition par les numero des noeuds correspondant

        for(j=0; j<nbrNoeuds; j++)
        {
            (populationBE+i)->phenotype[j] = j;
        }
        /// détermination de la parition
        for(j=0; j<nbrArretes; j++)
        {
            ///printf("(populationBE+%d)->genotype[%d] = %d \n",i,j,(populationBE+i)->genotype[j]);
            /// ajouter le 23/12/2015 pour déminuer l'impact de l'arrête fictive sur les solution générer
            ///if((population+i)->genotype[j]==0 && fluxVector[j]>0){
            if((populationBE+i)->genotype[j]==0)
            {
                nd = (edgeVector+j)->nouedDepart;
                na = (edgeVector+j)->nouedArrive;
                if((populationBE+i)->phenotype[na]!= na) /// il a été modifier auparavant
                {

                    if((populationBE+i)->phenotype[na] < (populationBE+i)->phenotype[nd])
                    {
                        ///dans le cas ou le noued d'arrivé est inférieur au noued depart ==> toutes les occurence du noued
                        /// de départ prendront ceux du noeud d'arrivé
                        numeroClustre = (populationBE+i)->phenotype[nd];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationBE+i)->phenotype[k] == numeroClustre)
                            {
                                (populationBE+i)->phenotype[k] = (populationBE+i)->phenotype[na];
                            }
                        }
                    }
                    else if((populationBE+i)->phenotype[na] > (populationBE+i)->phenotype[nd])
                    {
                        numeroClustre = (populationBE+i)->phenotype[na];
                        for(k=0; k<nbrNoeuds; k++)
                        {
                            if((populationBE+i)->phenotype[k] == numeroClustre)
                            {
                                (populationBE+i)->phenotype[k] = (populationBE+i)->phenotype[nd];
                            }
                        }
                    }
                }
                else
                {
                    (populationBE+i)->phenotype[na] = (populationBE+i)->phenotype[nd];
                }
            }
        }
        ///CHECKCHECKCHECKCHECKCHECKCHECKCCHECKCHECKCHECKCHECKCHECKCHECKC
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(phenotype,"%d ",(populationBE+i)->phenotype[j]);
        }
        fprintf(phenotype,"\n ");
        ///CHECKCHECKCHECKCHECKCHECKCHECKCCHECKCHECKCHECKCHECKCHECKCHECKC
    }
    fclose(phenotype);
}
///************************************************************************************************************
///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK
/// This function is improved on Tuesday  December 20th, 2022 and it is under test
///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK
void getPartitionFromSolutionBE(partitionBE *populationBE)
{
    ///CHECKCHECKCHECKCHECKCHECKCHECKCCHECKCHECKCHECKCHECKCHECKCHECKC
    //FILE *phenotype;
    //phenotype = fopen("E:/referencedGraphs/ThBE/phenotype.txt","a");
    ///CHECKCHECKCHECKCHECKCHECKCHECKCCHECKCHECKCHECKCHECKCHECKCHECKC

    int i,j,k,nd,na,numeroClustre; /// nd = noued de départ , na = noued d'arrivé
    for(i=0; i<taillePopulation; i++)
    {
        /// initialisation des phénotype = partition par les numero des noeuds correspondant

        for(j=0; j<nbrNoeuds; j++)
        {
            (populationBE+i)->phenotype[j] = j;
        }
        /// détermination de la parition
        for(j=0; j<nbrArretes; j++)
        {
            if((populationBE+i)->genotype[j]==0)
            {
                ///IIIIIIIIIIIIIIIIIIIIIII
                /// TOUJOURS ON A nd < na
                ///IIIIIIIIIIIIIIIIIIIIIII
                nd = (edgeVector+j)->nouedDepart;
                na = (edgeVector+j)->nouedArrive;

                if((populationBE+i)->phenotype[nd]== nd  && (populationBE+i)->phenotype[na]== na) (populationBE+i)->phenotype[na] = (populationBE+i)->phenotype[nd];
                else if((populationBE+i)->phenotype[nd]!= nd  && (populationBE+i)->phenotype[na]== na) (populationBE+i)->phenotype[na] = (populationBE+i)->phenotype[nd];

                else if((populationBE+i)->phenotype[nd]== nd  && (populationBE+i)->phenotype[na]!= na)(populationBE+i)->phenotype[nd] = (populationBE+i)->phenotype[na];
                else if((populationBE+i)->phenotype[nd]!= nd  && (populationBE+i)->phenotype[na]!= na)
                {
                    if((populationBE+i)->phenotype[nd] < (populationBE+i)->phenotype[na])
                    {
                        ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                        /// AFFECTER TOUTES LES OCCURENCES DE NA A ND
                        ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                        for(j=0; j<nbrNoeuds; j++)
                        {
                            if((populationBE+i)->phenotype[j] == (populationBE+i)->phenotype[na])
                                (populationBE+i)->phenotype[j] = (populationBE+i)->phenotype[nd];
                        }
                    }
                    else if((populationBE+i)->phenotype[na] < (populationBE+i)->phenotype[nd])
                    {
                        ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                        /// AFFECTER TOUTES LES OCCURENCES DE ND A NA
                        ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                        for(j=0; j<nbrNoeuds; j++)
                        {
                            if((populationBE+i)->phenotype[j] == (populationBE+i)->phenotype[nd])
                                (populationBE+i)->phenotype[j] = (populationBE+i)->phenotype[na];
                        }
                    }
                }
            }
        }
        ///########## CALCULER LE NOMBRE DE COMPOSANTES CONNEXES #########
        /*
         int clusterIndex = 0;
         for(int j=0; j<nbrNoeuds; j++)
         {
             int clustersSize =0,cluster[200];
             for(int k=0; k<nbrNoeuds; k++)
             {
                 if ((populationBE+i)->phenotype[k] == j)
                 {
                     cluster[clustersSize] = k;
                     clustersSize++;
                 }
             }
             if (clustersSize !=0)
             {
                 (populationBE+i)->clusters[clusterIndex].clusterId = j;
                 (populationBE+i)->clusters[clusterIndex].clusterSize = clustersSize;
                 for(int l = 0; l < clustersSize; l++)
                     (populationBE+i)->clusters[clusterIndex].clusterNoueds[l] =  cluster[l];
                 clusterIndex++;
             }
         }
         */

        ///CHECKCHECKCHECKCHECKCHECKCHECKCCHECKCHECKCHECKCHECKCHECKCHECKC
        /*
        fprintf(phenotype,"\n ");
        for(j=0; j<nbrNoeuds; j++)
        {
           fprintf(phenotype,"%d ",(populationBE+i)->phenotype[j]);
        }
        */
        ///CHECKCHECKCHECKCHECKCHECKCHECKCCHECKCHECKCHECKCHECKCHECKCHECKC
    }
    //fclose(phenotype);
}
///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK

///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFF writeBestSolutionInFileBE FFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
void writeBestSolutionInFileBE(partitionBE *solutionDominante, FILE *outputFile,int iteration)
{
    int i;
    fprintf(outputFile,"%d,%d,%d,%d,%d,%d,%0.2f,%d,%d,",nbrRun,iteration,solutionDominante->id,
            solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),solutionDominante->coutCoupeNormalise,
            solutionDominante->fitness, solutionDominante->contrainteViole,solutionDominante->nbrCluster);
    for(i=0; i<nbrArretes; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->genotype[i]);
    }
    fprintf(outputFile,",");
    for(i=0; i<nbrNoeuds; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->phenotype[i]);
    }
    fprintf(outputFile,"\n");

}

///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
/// calculer le cout de coupe de chaque partition ainsi que son fitness
///____________________________________________________________________________________________
/// on va appliquer la formule suivante pour le calcule des cout de coupe
/// T'(C) = B - T(C)            <=> la minimisation devient une maximisation
/// Z(C)=T'(C)+(u-v(C)) ×B      <=>  calculer la nouvelle fitness par rapport aux contrainte
/// T'(C) : cout de coupe (maximiser les flux INTRA), u ; nombre de contrainte totale
/// v(C): nombre de contrainte violé,  B: la somme totale des flux
///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

void calculCoutCoupeEtFitnessBE(partitionBE* populationBE)
{
    int sommeTotalCoutCoupe = 0;

    for(int i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationBE+i)->coutCoupe = 0; /// réinitialisation de la valriable après chaque itération (individus de la population)

        for(int j=0; j<nbrNoeuds-1; j++)
        {
            for(int k=j+1; k<nbrNoeuds; k++)
            {
                if((populationBE+i)->phenotype[j] == (populationBE+i)->phenotype[k] && buckUpFluxMatrix[j][k] > 0)
                {
                    (populationBE+i)->coutCoupe = (populationBE+i)->coutCoupe +  buckUpFluxMatrix[j][k];
                }
            }
        }

        (populationBE+i)->coutCoupeNormalise = (populationBE+i)->coutCoupe + ((nbr_constraint - (populationBE+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationBE+i)->coutCoupeNormalise;
    }

    ///######################################################################################
    /// Scaling Fintnss : sigma scaling Melanie mitchelle référence Goldberg | le 26/11/2015
    /// pour régler le problème des valeurs négatives des expected values
    /// ==> utiliser sigma truncation de goldberg
    /// g(f) = f + (moyenne - c * sigma )
    /// 1 <= c <= 5, généralement initilisé à 2
    ///######################################################################################
    float varianceCoutDeCoupeNormalise = 0, moyenneCoutDeCoupeNormalise, sigma, sommeTotalOfExpectdValues = 0;
    int c= 2;

    moyenneCoutDeCoupeNormalise = sommeTotalCoutCoupe / taillePopulation;
    for(int i=0; i<taillePopulation; i++)
    {
        varianceCoutDeCoupeNormalise =
                varianceCoutDeCoupeNormalise + pow(((populationBE+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(int i=0; i<taillePopulation; i++)
    {
        /// calcul de sigma truncation
        (populationBE+i)->expectedValue = (populationBE+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationBE+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(int i=0; i<taillePopulation ; i++)
    {
        (populationBE+i)->fitness = (float)((populationBE+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
        ///(populationBE+i)->fitness *= 1000;
    }

}
///************************************************************************************************************
void calculCoutCoupeEtFitnessWithFlowVectorBE(partitionBE* populationBE)
{
#if mouchard
    printf("calculCoutCoupeEtFitness ...\n");
#endif // mouchard
    int i,j,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationBE+i)->coutCoupe = 0;
        for(j=0; j<nbrArretes; j++) ///parcourir le genotype =  solution codée en binaire
        {
            if((populationBE+i)->genotype[j]==0)  /// on cherche à maximiser les intraCluster donc les arrêtes 0 != 1
            {
                (populationBE+i)->coutCoupe = (populationBE+i)->coutCoupe + fluxVector[j];
            }
        }

        (populationBE+i)->coutCoupeNormalise = (populationBE+i)->coutCoupe + ((nbr_constraint - (populationBE+i)->contrainteViole)*sommeTotalFlux);

        if((populationBE+i)->coutCoupeNormalise < 0) {
            printf("Problem : nbrConstraints = %d violated constraints = %d\n", nbr_constraint, (populationBE+i)->contrainteViole);
            exit(-8989);
        }
        /// désormis on va utiliser le cout de coupe normaliser pour le calcule des fitness
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationBE+i)->coutCoupeNormalise;
    }
#if scaling

    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/binaryEncoding/expectedValueBE.txt","w");
    ///======================================================================================================
    /// Scaling Fintnss : sigma scaling Melanie mitchelle référence Goldberg | le 26/11/2015
    ///pour régler le problème des valeurs négatives des expected values ==> utiliser sigma truncation de goldberg
    /// g(f) = f + (moyenne - c * sigma )
    /// 1 <= c <= 5, généralement initilisé à 2
    ///======================================================================================================
    float varianceCoutDeCoupeNormalise = 0, moyenneCoutDeCoupeNormalise, sigma, sommeTotalOfExpectdValues = 0;
    int c= 2;

    moyenneCoutDeCoupeNormalise = sommeTotalCoutCoupe / taillePopulation;
    for(i=0; i<taillePopulation; i++)
    {
        varianceCoutDeCoupeNormalise =
            varianceCoutDeCoupeNormalise + pow(((populationBE+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajouté
    for(i=0; i<taillePopulation; i++)
    {
        /// calcul de sigma truncation
        (populationBE+i)->expectedValue = (populationBE+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        /*            if(iteration == 1){
                        fprintf(file,"%d \t %d\n",(populationBE+i)->contrainteViole,(populationBE+i)->coutCoupeNormalise);
                    }*/
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationBE+i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for(i=0; i<taillePopulation ; i++)
    {
        (populationBE+i)->fitness = (float)((populationBE+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }

#else
    for(i=0; i<taillePopulation ; i++)
    {
        (populationBE+i)->fitness = (float)((populationBE+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling

}
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFF binaryEncoding : THE MAIN FUNCTION FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
int  binaryEncoding(int numeroGraphe, metrics *metricsBE)
{

    int bestSolution,iteration=1, bestSolutionIteration,ES=0;
    clock_t t1,t2;
    double temps;
    float bestFitness;
    srand(time(NULL));

    if(nbrRun == 1){
        ///###################  OPTIMAL SOLUTION PER RUN ####################
        sprintf(cheminOptimalSolutionBE, "../data/results/BE/OptimalSolution/OptimalSolutionThBE%d.csv", numeroGraphe);
        optimalSolutionFileBE = fopen(cheminOptimalSolutionBE, "w");

        if (optimalSolutionFileBE == NULL) {
            perror("Impossible d'ouvrir le fichier OptimalSolutionBE en ecriture\n");
            exit(EXIT_FAILURE);
        }

        fprintf(optimalSolutionFileBE, "Run N,"
                                       "intra cut size , inter cut size, Total flows,"
                                       "Run time,ES,Nbr arêtes,"
                                       "Nbr vertices ,Average weights,Variances of weights ,SD of weights, Nbr Edgess,Nbr Violated constraintes, Cont T_M_C,Cont Cohab, Cont NON Coha, Partition \n");

    }
    if(migration == 0) {
        ///############ ALLOCATION DE LA MEMOIRE POUR LES VARIABLES POPULATION ##########
        if (populationBE1 == NULL) populationBE1 = (partitionBE *) malloc(taillePopulation * sizeof(partitionBE));
        if (populationBE1 == NULL) printf("memory allocation failed for the populationBE1 \n");

        if (populationBE2 == NULL) populationBE2 = (partitionBE *) malloc(taillePopulation * sizeof(partitionBE));
        if (populationBE2 == NULL) printf("memory allocation failed for the populationBE2 \n");

        if (solutionDominanteBE == NULL)solutionDominanteBE = (partitionBE *) malloc(sizeof(partitionBE));
        if (solutionDominanteBE == NULL) printf("memory allocation failed for the solutionDominanteBE \n");

        if (bestSolutionOverRunsBE == NULL)bestSolutionOverRunsBE = (partitionBE *) malloc(sizeof(partitionBE));
        if (bestSolutionOverRunsBE == NULL) printf("memory allocation failed for the bestSolutionOverRuns \n");

        ///#### OPEN OUTPUT FILES HERE TO ADD SOME ADDITIONAL INFO ####
        sprintf(cheminBestSolutionBE, "../data/results/BE/bestSolutions/bestSolutionThBE%d.csv", numeroGraphe);
        sprintf(cheminFichierDetailsProbleme, "../data/results/DetailsProblem/DetailsProblem%d.csv", numeroGraphe);

        bestSolutionsOverIterationBE = fopen(cheminBestSolutionBE, "w");
        FichierDetailsProbleme = fopen(cheminFichierDetailsProbleme, "w");

        if (bestSolutionsOverIterationBE == NULL) {
            perror("Impossible d'ouvrir le fichier bestSolutionBE.csv en ecriture\n");
            exit(EXIT_FAILURE);
        }

        if (FichierDetailsProbleme == NULL) {
            perror("Impossible d'ouvrir le fichier DetailsProblem.csv en ecriture 111\n");
            exit(EXIT_FAILURE);
        }


        ///###################  BEST SOLUTION PER ITERATION #################
        fprintf(bestSolutionsOverIterationBE,
                "Numéro de l'exécution,Numéro de l'itération,ID Solution Dominante,Cout de coupe intra,"
                "Cout de Coupe inter,Cout de coupe Normalisé,Fitness,Nbr contraintes violées,Nombre de cluster,Genotype,Phenotype \n");

        ///####################################################################
    /// The best solution of the previous run is kept in the new generation
    ///####################################################################
    t1=clock(); ///START THE CHRONOMETER

    generatePopulationBE(populationBE1,0);
    getPartitionFromSolutionBE(populationBE1);
    checkContrainstAndFitnessPenalizationBE(populationBE1);
    calculCoutCoupeEtFitnessBE(populationBE1);

    ///affichePopulationBE(populationBE1);
    ///initialisation de la solution dominante
    bestSolution=findTheBestSolutionBE(populationBE1);
    *solutionDominanteBE=*(populationBE1+bestSolution);

    bestSolutionIteration = 1;

    nbrApparition=1;
    writeBestSolutionInFileBE(solutionDominanteBE,bestSolutionsOverIterationBE,(10*migration + iteration));
} else{
        t1=clock(); ///START THE CHRONOMETER
    }
    ///############## GA  ITERATIONS #########################
     for(; iteration <= 10; iteration++)
   /// while(nbrApparition <= 30)
    {
        printf(".");
        naturalSelectionBE(populationBE1,populationBE2);
        crossOverBE(populationBE2,populationBE1);
        mutationBE(populationBE1);
        getPartitionFromSolutionBE(populationBE1);
        checkContrainstAndFitnessPenalizationBE(populationBE1);
        calculCoutCoupeEtFitnessBE(populationBE1);

        bestSolution=findTheBestSolutionBE(populationBE1);

        if((populationBE1+bestSolution)->coutCoupeNormalise > solutionDominanteBE->coutCoupeNormalise )
        {
            *(solutionDominanteBE) = *(populationBE1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else
        {
            nbrApparition++;
        }
        writeBestSolutionInFileBE(solutionDominanteBE,bestSolutionsOverIterationBE,bestSolutionIteration);
        //affichePopulationBE(populationBE1);

        ///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK
        //printf("CN = %d\n", (populationBE1+bestSolution)->coutCoupeNormalise);
        ///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK

        ///REPRODUCTION NEW SOLUTIONS FOR DIVERSIFICATION
        solutionsReproduction(populationBE1);

        //iteration++;
    }

    t2=clock();
    //printf("\n");
    metricsBE->runTime += (float)(t2-t1)/(float)CLOCKS_PER_SEC;


    if(bestSolutionIteration>=2)
    {
        metricsBE->ES += ((taillePopulation-elitismeRate)*(bestSolutionIteration-2))+(taillePopulation+ solutionDominanteBE->id -elitismeRate+1);
    }
    else
    {
        metricsBE->ES += solutionDominanteBE->id +1;
    }
}
///************************************************************************************************************
void checkContrainstAndFitnessPenalizationBE(partitionBE *populationBE)
{
    ///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK
    //int countNoued = 0;
    //FILE *clustersSize;
    //clustersSize = fopen("E:/referencedGraphs/ThBE/clustersSize.txt", "a");
    ///CHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECKCHECK
    int i,j,k;
    int tailleClusters;
    for(i=0; i<taillePopulation; i++)
    {
        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        ///initialisation de vecteur des contraintes et de la variable
        ///de contraintes violées
        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        (populationBE+i)->contrainteViole=0;
        for(j=0; j<3; j++)
        {
            (populationBE+i)->constraintVector[j]=0;
        }

        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
        /// CALCUL THE NUMBER OF CLUSTERS AND THEIR SIZES
        ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

        /// THIS TASK IS DONE IN THE : getPartitionFromSolution FUNCTION
        /** ------------------------------------------------------------*/
        (populationBE+i)->nbrCluster=0;
        for(j=0; j<nbrNoeuds; j++)
        {
            tailleClusters =0;
            for(k=0; k<nbrNoeuds; k++){
                if ((populationBE+i)->phenotype[k] == j){
                    (populationBE+i)->clusters[(populationBE+i)->nbrCluster].clusterNoueds[tailleClusters]=k;
                    tailleClusters++;
                }
            }
            if (tailleClusters !=0){
                (populationBE+i)->clusters[(populationBE+i)->nbrCluster].clusterId = (populationBE+i)->nbrCluster;
                (populationBE+i)->clusters[(populationBE+i)->nbrCluster].clusterSize = tailleClusters;
                (populationBE+i)->nbrCluster++;
            }
        }
        /*-----------------------------------------------------------------*/


        for(j=0; j < (populationBE+i)->nbrCluster; j++)
        {
            if((populationBE+i)->clusters[j].clusterSize  > max_sizeCluster)
            {
                (populationBE+i)->constraintVector[0]++;

            }
        }
        if((populationBE+i)->constraintVector[0] != 0)
        {
            (populationBE+i)->contrainteViole++;
        }
    }

    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstrainteBE(populationBE);
    }
}
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
///FFFFFFFFFFFFF checkCohabitationAndNonCohabitationConstraint FFFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
void checkCohabitationAndNonCohabitationConstrainteBE(partitionBE *populationBE)
{

    int i,j,noeud1,noeud2;
    for(i=0; i < taillePopulation; i++)
    {
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationBE+i)->phenotype[noeud1]!= (populationBE+i)->phenotype[noeud2])
                {
                    (populationBE+i)->constraintVector[1]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationBE+i)->constraintVector[1]!=0)
            {
                (populationBE+i)->contrainteViole++;
            }
        }

        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;

                if((populationBE+i)->phenotype[noeud1]== (populationBE+i)->phenotype[noeud2])
                {
                    (populationBE+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes violées
            if((populationBE+i)->constraintVector[2]!=0)
            {
                (populationBE+i)->contrainteViole++;
            }
        }
    }

}
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFF writeOptimalSolutionInFileBE FFFFFFFFFFFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
void writeOptimalSolutionInFileBE(partitionBE *solutionDominante,FILE* optimalSolutionFileBE,
                                  int nbrRun, metrics metricsBE,float moyenneDesPoids,
                                  float variancePoids, float standardDeviation)
{
    if(optimalSolutionFileBE == NULL) perror("the file was closed\n");
    int i,j,k;
    fprintf(optimalSolutionFileBE,"%d,",nbrRun);
    fprintf(optimalSolutionFileBE,"%d,%d,%d,%.3f,%.3f,%d,", solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,metricsBE.runTime, metricsBE.ES, fixedNbrArretes);
    fprintf(optimalSolutionFileBE,"%d,%0.2f,%0.2f,%0.2f,%d,%d,",nbrNoeuds,moyenneDesPoids,variancePoids,standardDeviation,nbrArretes,solutionDominante->contrainteViole);
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///RAJOUTER LES CONTRAINTES VIOLEES
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    fprintf(optimalSolutionFileBE,"%d,%d,%d,",solutionDominante->constraintVector[0],
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




///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// USED BY THE QSORT
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int compareCroissantFitnessBE (void const *a, void const *b)
{
    partitionBE const *pa = a;
    partitionBE const *pb = b;
    return  (pb->coutCoupeNormalise - pa->coutCoupeNormalise);
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// USED BY THE QSORT
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int compareDecroissantFitnessBE (void const *a, void const *b)
{
    partitionBE const *pa = a;
    partitionBE const *pb = b;
    return  (pa->coutCoupeNormalise - pb->coutCoupeNormalise);
}


///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
///FFFFFFFFFFF solutionsReproduction FFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
void solutionsReproduction(partitionBE* populationBE)
{
    int allele,nbrZeroArretes;
    nbrZeroArretes = nbrNoeuds / max_sizeCluster;
    if(nbrZeroArretes != (double)nbrNoeuds / max_sizeCluster)
        nbrZeroArretes++;

    qsort(populationBE, taillePopulation, sizeof *populationBE, compareCroissantFitnessBE);
    /// REPRODUCTION OF NEW INDIVIDUALS
    for(int i= taillePopulation - 10; i<taillePopulation; i++)
    {
        /// THE SOLUTION ID
        (populationBE+i)->id =i;
        /// GENOTYPE INITIALIZATION
        for(int j=0; j<nbrArretes; j++) {
            (populationBE+i)->genotype[j] = 1;
        }
        /// FIX THE VALUES OF CC EDGES
        for(int j=0; j<nbrCohabitationConstraintes; j++)
        {
            for(int k=0; k<nbrArretes; k++)
            {
                if((cohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedDepart &&
                    cohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedArrive) ||
                   (cohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedArrive &&
                    cohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedDepart))
                {
                    (populationBE+i)->genotype[k] = 0;
                }
            }
        }

        /// FIX THE VALUES OF NCC EDGES
        for(int j=0; j<nbrNonCohabitationConstraintes; j++)
        {
            for(int k=0; k<nbrArretes; k++)
            {
                if((nonCohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedDepart &&
                    nonCohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedArrive) ||
                   (nonCohabitationConstraintes[j].nouedDepart == edgeVector[k].nouedArrive &&
                    nonCohabitationConstraintes[j].nouedArrive == edgeVector[k].nouedDepart) )
                {
                    (populationBE+i)->genotype[k] = 1;
                }
            }
        }

        for(int j=0; j<nbrZeroArretes; j++){
            do{
                allele = rnd(0,nbrArretes);
            }
            while((populationBE+i)->genotype[allele]==0);
            (populationBE+i)->genotype[allele] = 0;
        }

    }
}


///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Migration from DVTC to BE
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
///FFFFFFFFFFF solutionsReproduction FFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
void migrationFromDVTCtoBE(){
    qsort(populationBE1, taillePopulation, sizeof *populationBE1, compareDecroissantFitnessBE);
    for(int i = 0; i <nbrNoeuds; i++)
        (populationBE1+taillePopulation-1)->phenotype[i] = solutionDominanteDVTC->genotype[i];

    (populationBE1+taillePopulation-1)->id = taillePopulation - 1;
    (populationBE1+taillePopulation-1)->coutCoupe = solutionDominanteBE->coutCoupe;
    (populationBE1+taillePopulation-1)->coutCoupeNormalise = solutionDominanteBE->coutCoupeNormalise;
    (populationBE1+taillePopulation-1)->contrainteViole = solutionDominanteBE->contrainteViole;
    (populationBE1+taillePopulation-1)->nbrCluster = solutionDominanteBE->nbrCluster;

    for(int i = 0; i <(populationBE1+taillePopulation-1)->nbrCluster; i++)
        (populationBE1+taillePopulation-1)->clusters[i] = solutionDominanteBE->clusters[i];

    for(int i = 0; i < nbr_constraint; i++)
        (populationBE1+taillePopulation-1)->constraintVector[i] = solutionDominanteBE->constraintVector[i];
    int sommeTotalCoutCoupe = 0;
    for(int i = 0; i < taillePopulation; i++)
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationBE1+i)->coutCoupeNormalise;

    float varianceCoutDeCoupeNormalise = 0, moyenneCoutDeCoupeNormalise, sigma, sommeTotalOfExpectdValues = 0;
    int c= 2;

    moyenneCoutDeCoupeNormalise = sommeTotalCoutCoupe / taillePopulation;
    for(int i=0; i<taillePopulation; i++)
    {
        varianceCoutDeCoupeNormalise =
                varianceCoutDeCoupeNormalise + pow(((populationBE1+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);

    /// calcul de sigma truncation
    (populationBE1+taillePopulation-1)->expectedValue = (populationBE1+taillePopulation-1)->coutCoupeNormalise +
                                                        (moyenneCoutDeCoupeNormalise - c*sigma);

    for(int i=0; i< taillePopulation;i++)
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationBE1+i)->expectedValue;

    for(int i=0; i<taillePopulation ; i++)
        (populationBE1+i)->fitness = (float)((populationBE1+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);

    //Encoding the solution using the binary encoding
    int edgeIndex;
    for(int i = 0; i < nbrArretes; i++ ) (populationBE1 + taillePopulation - 1)->genotype[i] = 1;
    for(int i = 0; i < (populationBE1+taillePopulation-1)->nbrCluster; i++) {
        for (int j = 0; j < (populationBE1 + taillePopulation - 1)->clusters[i].clusterSize - 1; j++) {
            for (int k = j + 1; k < (populationBE1 + taillePopulation - 1)->clusters[i].clusterSize; k++) {
                int nd = (populationBE1 + taillePopulation - 1)->clusters[i].clusterNoueds[j];
                int na = (populationBE1 + taillePopulation - 1)->clusters[i].clusterNoueds[k];
                if (fluxMatrix[nd][na] > 0) {
                    edgeIndex = searchForEdgeIndex(nd, na);
                }

                if (edgeIndex)
                    (populationBE1 + taillePopulation - 1)->genotype[edgeIndex] = 0;
            }
        }
    }
}


///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// searchForEdgeIndex
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int searchForEdgeIndex(int nd, int na){
    int found = 0;
    for(int i = 0; i <fixedNbrArretes; i++){
        if((edgeVector[i].nouedDepart == nd && edgeVector[i].nouedArrive == na) ||
                (edgeVector[i].nouedDepart == na && edgeVector[i].nouedArrive == nd)){
            return edgeVector[i].numeroEdge;
        }
    }
    if(!found)
        perror("This edge does not exist \n");
}