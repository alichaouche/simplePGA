#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../headers/typeDeclaration.h"
#include "../headers/communesFunctions.h"
#include "../headers/genetic_alg_FC.h"
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
extern partitionFC *populationFC1 = NULL,*populationFC2 = NULL, *solutionDominanteFC = NULL,*bestSolutionOverRunsFC = NULL;


///**************************************************************************************
///g�n�ration de la population initiale
void generatePopulationFC(partitionFC *populationFC, int indiceFirstElt)
{

    int i, j;
    for (i = indiceFirstElt; i < taillePopulation; i++)
    {
        ///srand(time(NULL));
        for (j = 0; j <= nbrNoeuds; j++) /// <= pour qu'on aura un �l�ment en plus pour le nombre des parties
        {
#if Windows
            (populationFC + i)->genotype[j] = drand48ForWindowsFC(0, 101);
#else
            (populationFC + i)->genotype[j] = drand48();
#endif // Windows
        }
        (populationFC + i)->id = i;
    }
}
///**************************************************************************************
/// affichage des population
void affichePopulationFC(partitionFC *populationFC)
{
    int i, j;
    for (i = 0; i < taillePopulation; i++)
    {
        printf("\nid = %d \t cout de coupe = %d \t fitness = %0.4f \n",
               (populationFC + i)->id, (populationFC + i)->coutCoupe, (populationFC + i)->fitness);
        printf("le cout de coupe normalise = %d \n", (populationFC + i)->coutCoupeNormalise);
        printf("le nombre des contraintes violees = %d \n", (populationFC + i)->contrainteViole);
        printf("le nombre des cluster de la partition = %d \n", (populationFC + i)->nbrCluster);
        printf("la genotype est  = \t");
        for (j = 0; j <= nbrNoeuds; j++)
        {
            printf("%0.2f ", (populationFC + i)->genotype[j]);
        }
        printf("\n");
        printf("la phenotype est  = \t");
        for (j = 0; j <= nbrNoeuds; j++)
        {
            printf("%d ", (populationFC + i)->phenotype[j]);
        }
        printf("\n_________________________________________________________________________\n");
    }
}
///************************************************************************************************************
/// calculer le cout de coupe de chaque partition ainsi que son fitness
///***************************************************************************************
/// on va appliquer la formule suivante pour le calcule des cout de coupe
/// T'(C) = B - T(C)            <=> la minimisation devient une maximisation
/// Z(C)=T'(C)+(u-v(C)) �B      <=>  calculer la nouvelle fitness par rapport aux contrainte
/// T'(C) : cout de coupe (maximiser les flux INTRA), u ; nombre de contrainte totale
/// v(C): nombre de contrainte viol�,  B: la somme totale des flux
///***************************************************************************************

void calculCoutCoupeEtFitnessFC(partitionFC *populationFC)
{

    int i, j, k, sommeTotalCoutCoupe = 0;
    for (i = 0; i < taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationFC + i)->coutCoupe = 0; /// r�initialisation de la valriable apr�s chaque it�ration (individus de la population)

        for (j = 0; j < nbrNoeuds - 1; j++)
        {
            for (k = j + 1; k < nbrNoeuds; k++)
            {
                if ((populationFC + i)->phenotype[j] == (populationFC + i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (populationFC + i)->coutCoupe = (populationFC + i)->coutCoupe + fluxMatrix[j][k];
                }
            }
        }

        (populationFC + i)->coutCoupeNormalise = (populationFC + i)->coutCoupe + ((nbr_constraint - (populationFC + i)->contrainteViole) * sommeTotalFlux);

        /// d�sormis on va utiliser le cout de coupe normaliser pour le calcule des fitness

        ///sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationFC+i)->coutCoupe;
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationFC + i)->coutCoupeNormalise;
    }

#if scaling

    ///======================================================================================================
    /// Scaling Fintnss : sigma scaling Melanie mitchelle r�f�rence Goldberg | le 26/11/2015
    ///pour r�gler le probl�me des valeurs n�gatives des expected values ==> utiliser sigma truncation de goldberg
    /// g(f) = f + (moyenne - c * sigma )
    /// 1 <= c <= 5, g�n�ralement initilis� � 2
    ///======================================================================================================
    float varianceCoutDeCoupeNormalise = 0, moyenneCoutDeCoupeNormalise, sigma, sommeTotalOfExpectdValues = 0;
    int c = 2;

    moyenneCoutDeCoupeNormalise = sommeTotalCoutCoupe / taillePopulation;
    for (i = 0; i < taillePopulation; i++)
    {
        varianceCoutDeCoupeNormalise =
            varianceCoutDeCoupeNormalise + pow(((populationFC + i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise), 2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajout�
    for (i = 0; i < taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationFC + i)->expectedValue = (populationFC + i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c * sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationFC + i)->expectedValue;
    }

    /// calcule du fitness
    /// le 26/11/2015 : calcul du fitness des solution avec leur expectedValues
    for (i = 0; i < taillePopulation; i++)
    {
        (populationFC + i)->fitness = (float)((populationFC + i)->expectedValue) / (float)(sommeTotalOfExpectdValues);
    }

#else
    FILE *file;
    file = fopen("D:/WorkingFolderAgProgramming/benchmark/fractionnalEncoding/expectedValueFE.txt", "w");

    for (i = 0; i < taillePopulation; i++)
    {
        (populationFC + i)->fitness = (float)((populationFC + i)->coutCoupeNormalise) / (float)(sommeTotalCoutCoupe);
    }
#endif // scaling
}
///**************************************************************************************
void naturalSelectionFC(partitionFC *populationFC1, partitionFC *populationFC2)
{
    int i, j, maxFitness;
    int tmpVectIndice[taillePopulation];
    for (i = 0; i < taillePopulation; i++)
        tmpVectIndice[i] = -1;
    float lotterie, sommeFitness;

    for (i = 0; i < elitismeRate; i++)
    {
        maxFitness = 0;
        for (j = 0; j < taillePopulation; j++)
        {
            if ((populationFC1 + maxFitness)->coutCoupeNormalise < (populationFC1 + j)->coutCoupeNormalise && tmpVectIndice[j] == -1)
            {
                maxFitness = j;
            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationFC2 + i) = *(populationFC1 + maxFitness);
        ///printf("(populationFC2+%d)->coutCoupeNormalise = %d \n",i,(populationFC2+i)->coutCoupeNormalise);
        ///mon_sleep(pause);
    }
    ///printf("\n\n\n");

    ///qsort (populationFC2, taillePopulation, sizeof *populationFC2, compareCroissantFitnessEE);
    ///sortingPopulationByFintness(populationFC2,tmpPopulation);

    for (i = elitismeRate; i < taillePopulation; i++)
    {
#if Windows
        lotterie = drand48ForWindows(0, 1001);
#else
        lotterie = drand48();
#endif // Windows
        sommeFitness = 0;

        for (j = 0; j < taillePopulation; j++)
        {
            sommeFitness = sommeFitness + (populationFC1 + j)->fitness;
            if (lotterie <= sommeFitness)
            {
                *(populationFC2 + i) = *(populationFC1 + j);
                break;
            }
            else if (j == taillePopulation - 1)
            {
                *(populationFC2 + i) = *(populationFC1 + j);
            }
        }
    }
    ///free(tmpPopulation);
}
///**************************************************************************************
void crossOverFC(partitionFC *populationFC1, partitionFC *populationFC2)
{
    int i, j;
    int choixInd1, choixInd2, choixLocus;
    ///***********************************************************************************
    for (i = 0; i < elitismeRate; i++)
    {
        *(populationFC2 + i) = *(populationFC1 + i);
        (populationFC2 + i)->id = i;
    }
    ///***********************************************************************************
#if NEW_CROSSOVER
    while (i < taillePopulation - regeneration)
    {
        ///choixInd1 = rnd(elitismeRate,taillePopulation); /** l'intervalle sera [0;taillePopulation-1] */
        choixInd1 = rnd(0, taillePopulation);
        do
        {
            /// ce traitement �vite l'appariement d'un element avec lui m�me/
            ///choixInd2 = rnd(elitismeRate,taillePopulation);
            choixInd2 = rnd(0, taillePopulation);
        } while (choixInd2 == choixInd1);

        choixLocus = rnd(0, nbrNoeuds - 1); /** l'intervalle sera [1;nbrNoeuds-2] */

        for (j = 0; j <= choixLocus; j++) /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        {
            (populationFC2 + i)->genotype[j] = (populationFC1 + choixInd1)->genotype[j];
            (populationFC2 + i + 1)->genotype[j] = (populationFC1 + choixInd2)->genotype[j];
        }
        ///********************************************************************************
        for (j = choixLocus + 1; j <= nbrNoeuds; j++) /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        {
            (populationFC2 + i)->genotype[j] = (populationFC1 + choixInd2)->genotype[j];
            (populationFC2 + i + 1)->genotype[j] = (populationFC1 + choixInd1)->genotype[j];
        }

        /// affectation de nouveau indice
        (populationFC2 + i)->id = i;
        (populationFC2 + i + 1)->id = i + 1;
        i += 2;
    }
    generatePopulationFC(populationFC2, taillePopulation - regeneration);
#else

    while (i < taillePopulation)
    {
        ///choixInd1 = rnd(elitismeRate,taillePopulation); /** l'intervalle sera [0;taillePopulation-1] */
        choixInd1 = rnd(0, taillePopulation); /** l'intervalle sera [0;taillePopulation-1] */
        do
        {
            /// ce traitement �vite l'appariement d'un element avec lui m�me/
            ///choixInd2 = rnd(elitismeRate,taillePopulation);
            choixInd2 = rnd(0, taillePopulation);
        } while (choixInd2 == choixInd1);

        choixLocus = rnd(0, nbrNoeuds - 1); /** l'intervalle sera [1;nbrNoeuds-2] */

        for (j = 0; j <= choixLocus; j++) /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        {
            (populationFC2 + i)->genotype[j] = (populationFC1 + choixInd1)->genotype[j];
            (populationFC2 + i + 1)->genotype[j] = (populationFC1 + choixInd2)->genotype[j];
        }
        ///********************************************************************************
        for (j = choixLocus + 1; j <= nbrNoeuds; j++) /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        {
            (populationFC2 + i)->genotype[j] = (populationFC1 + choixInd2)->genotype[j];
            (populationFC2 + i + 1)->genotype[j] = (populationFC1 + choixInd1)->genotype[j];
        }

        /// affectation de nouveau indice
        (populationFC2 + i)->id = i;
        (populationFC2 + i + 1)->id = i + 1;
        i += 2;
    }

#endif
}
///**************************************************************************************
void calculPhenotypeFC(partitionFC *populationFC)
{

    /**<
    � partir d'une solution g�n�r�e al�atoirement on d�duit la partition
    1- multiplier la derni�re valeur * max_clusters pour avoir le nombre de cluster en prennant la borne sup
    2- la valeur obtenu dans 1 est utilis�e pour d�terminer l'affectation des noeuds aux clusters
     */
    int i, j, nombreDesParties;
    for (i = 0; i < taillePopulation; i++) /// parcourir tous les individus de la population
    {
        nombreDesParties = ceil((populationFC + i)->genotype[nbrNoeuds] * (max_clusters - 1)); /// r�initialisation de la valriable apr�s chaque it�ration (individus de la population)
        (populationFC + i)->phenotype[nbrNoeuds] = nombreDesParties;
        for (j = 0; j < nbrNoeuds; j++)
        {
            (populationFC + i)->phenotype[j] = floor((populationFC + i)->genotype[j] * nombreDesParties) + 1;
        }
    }
}
///**************************************************************************************
int findTheBestSolutionFC(partitionFC *populationFC)
{
    float maxFitness = 0;
    int i, indice = 0;
    for (i = 0; i < taillePopulation; i++)
    {
        if (maxFitness < (populationFC + i)->fitness)
        {
            maxFitness = (populationFC + i)->fitness;
            indice = i;
        }
    }
    return indice;
}
///**************************************************************************************
void mutationFC(partitionFC *populationFC)
{

    int i, numeroGeneMute = -1;
    double applicationOfMutation;
    for (i = elitismeRate; i < taillePopulation; i++)
    {
#if Windows
        applicationOfMutation = drand48ForWindows(0, 1001);
#else
        applicationOfMutation = drand48();
#endif // Windows
        if (applicationOfMutation <= tauxMutation)
        {
            /// choisir un g�ne al�atoirement pour lui appliquer la mutation
            numeroGeneMute = rnd(0, nbrNoeuds + 1);
#if Windows
            (populationFC + i)->genotype[numeroGeneMute] = drand48ForWindows(0, 101);
#else
            (populationFC + i)->genotype[numeroGeneMute] = drand48();
#endif // Windows
        }
    }
}
///**************************************************************************************
float testerLaSommeDesFitnessFC(partitionFC *populationFC)
{

#if mouchard
    printf("testerLaSommeDesFitness ...\n");
#endif // mouchard
    int i;
    float sommeFitness = 0;
    for (i = 0; i < taillePopulation; i++)
    {
        sommeFitness = sommeFitness + (populationFC + i)->fitness;
    }
    return sommeFitness;
}
///**************************************************************************************
void displayTheBestSolutionFC(partitionFC *solutionDominante)
{
    int i;
    printf("id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
           solutionDominante->id, solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux - solutionDominante->coutCoupe));
    printf("le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for (i = 0; i < 4; i++)
    {
        switch (i)
        {

        case 0:
            if (solutionDominante->constraintVector[i] == 0)
                printf("contrainte sur le nombre de clusters : OK \n");
            else
                printf("contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if (solutionDominante->constraintVector[i] == 0)
                printf("contrainte sur la taille des clusters : OK \n");
            else
                printf("contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if (solutionDominante->constraintVector[i] == 0)
                printf("contrainte de cohabitation : OK \n");
            else
                printf("contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if (solutionDominante->constraintVector[i] == 0)
                printf("contrainte Non Cohabitation : OK \n");
            else
                printf("contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            printf("Cette contrainte n'est pas prise en charge par le syst�me \n");
        }
    }
}

///**************************************************************************************
void writeSolutionInFileFC(partitionFC *populationFC, FILE *outputFilePop, int iteration)
{

    int i, j;
    for (i = 0; i < taillePopulation; i++)
    {

        ///�criture de la meilleur solution dans le fichier output
        fprintf(outputFilePop, "%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t", iteration, (populationFC + i)->id,
                (populationFC + i)->coutCoupe, (populationFC + i)->fitness, (populationFC + i)->coutCoupeNormalise,
                (populationFC + i)->contrainteViole, (populationFC + i)->nbrCluster);
        for (j = 0; j <= nbrNoeuds; j++)
        {
            fprintf(outputFilePop, "%0.2f ", (populationFC + i)->genotype[j]);
        }
        fprintf(outputFilePop, "\t\t\t\t");
        for (j = 0; j <= nbrNoeuds; j++)
        {
            fprintf(outputFilePop, "%d ", (populationFC + i)->phenotype[j]);
        }
        fprintf(outputFilePop, "\n");
    }
    /// fprintf(outputFilePop,"\n===============================================\n");
    fprintf(outputFilePop, "\n\n");
}
///**************************************************************************************
void writeBestSolutionInFileFC(partitionFC *solutionDominante, FILE *outputFile, int iteration)
{

    int i;
    fprintf(outputFile, "%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t", iteration, solutionDominante->id,
            solutionDominante->coutCoupe, solutionDominante->fitness, solutionDominante->coutCoupeNormalise,
            solutionDominante->contrainteViole, solutionDominante->nbrCluster);
    for (i = 0; i <= nbrNoeuds; i++)
    {
        fprintf(outputFile, "%0.2f ", solutionDominante->genotype[i]);
    }
    fprintf(outputFile, "\t\t\t\t");
    for (i = 0; i <= nbrNoeuds; i++)
    {
        fprintf(outputFile, "%d ", solutionDominante->phenotype[i]);
    }
    fprintf(outputFile, "\n");
}

///=======================================================================================
void fractionalEncoding(int nbrGeneration, FILE *outputFileFC, FILE *outputFilePopFC, FILE *outputOptimalSolutionFileFC, partitionFC *populationFC1, partitionFC *populationFC2, partitionFC *solutionDominante)
{

    int bestSolution, iteration = 1, bestSolutionIteration, ES = 0;
    clock_t t1, t2;
    double temps;

    t1 = clock();

    /// g�n�ration de la population initiale
    generatePopulationFC(populationFC1, 0);
    calculPhenotypeFC(populationFC1);
    checkContrainstAndFitnessPenalizationFC(populationFC1);
    calculCoutCoupeEtFitnessFC(populationFC1);
    bestSolution = findTheBestSolutionFC(populationFC1);
    *solutionDominante = *(populationFC1 + bestSolution);
    bestSolutionIteration = 1;
    nbrApparition = 1;

    //for(iteration =2 ; iteration <= nbrGeneration; iteration++){
    do
    {
        ///=====================================================================================================
        naturalSelectionFC(populationFC1, populationFC2);
        crossOverFC(populationFC2, populationFC1); /// nbrNoeuds repr�sente la taille des genotype
        mutationFC(populationFC1);
        calculPhenotypeFC(populationFC1);
        checkContrainstAndFitnessPenalizationFC(populationFC1);
        calculCoutCoupeEtFitnessFC(populationFC1);
        ///=====================================================================================================
        /// dans ce qui suit, on va juste localiser la meilleur solution
        ///=====================================================================================================

        bestSolution = findTheBestSolutionFC(populationFC1);
        if ((populationFC1 + bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise)
        /** le 23/12/2015 : pour avoir la meilleur solution avec le moins de contrainte viol�es*/
        /// && (populationFC1+bestSolution)->contrainteViole   < solutionDominante->contrainteViole)
        {
            *(solutionDominante) = *(populationFC1 + bestSolution);
            nbrApparition = 1;
            bestSolutionIteration = iteration;
        }
        else if ((populationFC1 + bestSolution)->coutCoupeNormalise <= solutionDominante->coutCoupeNormalise)
        {
            nbrApparition++;
        }
        ///affichePopulation(populationFC1,taillePopulation,nbrNoeuds);

    } while (nbrApparition < max_steady_generation);

    t2 = clock();
    temps = (double)(t2 - t1) / CLOCKS_PER_SEC;
    displayTheBestSolutionFC(solutionDominante);
    if (bestSolutionIteration >= 2)
    {
        ES = ((taillePopulation - elitismeRate) * (bestSolutionIteration - 2)) + (taillePopulation + solutionDominante->id - elitismeRate + 1);
    }
    else
    {
        ES = solutionDominante->id + 1;
    }
    writeOptimalSolutionInFileFC(solutionDominante, outputOptimalSolutionFileFC, nbrRun, bestSolutionIteration, temps, ES);
    printf("le temps d execution est %lf \n", temps);

    fprintf(outputFileFC, "\n");
    writeBestSolutionInFileFC(solutionDominante, outputFileFC, bestSolutionIteration);
}

///**************************************************************************************
void checkContrainstAndFitnessPenalizationFC(partitionFC *populationFC)
{

    int i, j, k;

    for (i = 0; i < taillePopulation; i++)
    {
        (populationFC + i)->nbrCluster = 0; /// le nombre de cluster doit etre initialis� � 1 pour prendre en cons�diration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes viol�es
        (populationFC + i)->contrainteViole = 0;
        for (j = 0; j < 4; j++)
        {
            (populationFC + i)->constraintVector[j] = 0;
        }
        ///**********************************************************************************

        for (j = 0; j < nbrNoeuds; j++)
        {
            (populationFC + i)->clustersSize[j] = 0;
            for (k = 0; k < nbrNoeuds; k++)
            {
                if ((populationFC + i)->phenotype[k] == j)
                {
                    (populationFC + i)->clustersSize[j]++;
                }
            }
            if ((populationFC + i)->clustersSize[j] != 0)
            {
                (populationFC + i)->nbrCluster++;
            }
        }

        if ((populationFC + i)->nbrCluster > max_clusters || (populationFC + i)->nbrCluster < min_clusters)
        {
            (populationFC + i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationFC + i)->contrainteViole++;
        }

        for (j = 0; j < nbrNoeuds; j++)
        {
            if ((populationFC + i)->clustersSize[j] != 0)
            {
                if ((populationFC + i)->clustersSize[j] > max_sizeCluster || (populationFC + i)->clustersSize[j] < min_sizeCluster)
                {
                    (populationFC + i)->constraintVector[1]++;
                }
            }
        }

        if ((populationFC + i)->constraintVector[1] != 0)
        {
            (populationFC + i)->contrainteViole++;
        }
    }
    if (nbrCohabitationConstraintes != 0 || nbrNonCohabitationConstraintes != 0)
    {
        checkCohabitationAndNonCohabitationConstrainteFC(populationFC);
    }
}

///**************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteFC(partitionFC *populationFC)
{

    int i, j, noeud1, noeud2;
    for (i = 0; i < taillePopulation; i++)
    {
        ///24/11/2015 normalement cette taches est d�j� faite dans le fonction : checkContrainstAndFitnessPenalizationFC
        /**
        (populationFC+i)->constraintVector[2]=0;
        (populationFC+i)->constraintVector[3]=0;
        */
        if (nbrCohabitationConstraintes != 0)
        {
            for (j = 0; j < nbrCohabitationConstraintes; j++)
            {
                noeud1 = cohabitationConstraintes[j].nouedDepart;
                noeud2 = cohabitationConstraintes[j].nouedArrive;

                if ((populationFC + i)->phenotype[noeud1] != (populationFC + i)->phenotype[noeud2])
                {
                    (populationFC + i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes viol�es
            if ((populationFC + i)->constraintVector[2] != 0)
            {
                (populationFC + i)->contrainteViole++;
            }
        }

        if (nbrNonCohabitationConstraintes != 0)
        {
            for (j = 0; j < nbrNonCohabitationConstraintes; j++)
            {
                noeud1 = nonCohabitationConstraintes[j].nouedDepart;
                noeud2 = nonCohabitationConstraintes[j].nouedArrive;

                if ((populationFC + i)->phenotype[noeud1] == (populationFC + i)->phenotype[noeud2])
                {
                    (populationFC + i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes viol�es
            if ((populationFC + i)->constraintVector[3] != 0)
            {
                (populationFC + i)->contrainteViole++;
            }
        }
    }
}

///********************************************************************************************
void writeOptimalSolutionInFileFC(partitionFC *solutionDominante, FILE *outputOptimalSolutionFileFC,
                                  int nbrRun, int bestSolutionIteration, float runTime, int ES)
{

    fprintf(outputOptimalSolutionFileFC, "RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d | ",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe, (sommeTotalFlux - solutionDominante->coutCoupe), sommeTotalFlux, runTime, ES, solutionDominante->contrainteViole);

    fprintf(outputOptimalSolutionFileFC, " max_clusters = %d | \
            max_sizeCluster = %d  |  \
            cohabitation = %d | \
            non cohabitation = %d\n",

            solutionDominante->constraintVector[0],
            solutionDominante->constraintVector[1],
            solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);
}
///******************************************************************************************
int compareCroissantFitnessFC(void const *a, void const *b)
{

    partitionFC const *pa = a;
    partitionFC const *pb = b;
    return pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}

double drand48ForWindowsFC(int v1, int v2)
{
    ///return (double)(rand()) / (double)(RAND_MAX);
    return (double)(rnd(v1, v2)) / 100.0;
}

/**


void writeOptimalSolutionInFileFC(partitionFC *solutionDominante,FILE* outputOptimalSolutionFileFC)
{
    int i;
    fprintf(outputOptimalSolutionFileFC,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
            solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFileFC,"le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileFC,"contrainte sur le nombre de clusters : OK \n");
            else fprintf(outputOptimalSolutionFileFC,"contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileFC,"contrainte sur la taille des clusters : OK \n");
            else fprintf(outputOptimalSolutionFileFC,"contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileFC,"contrainte de cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileFC,"contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileFC,"contrainte Non Cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileFC,"contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            fprintf(outputOptimalSolutionFileFC,"Cette contrainte n'est pas prise en charge par le syst�me \n");

        }

    }

}
void checkContrainstAndFitnessPenalizationFC(partitionFC *populationFC)
{

    int i,j,k,m=nbrNoeuds,tabTmp[m],tmp;

    for(i=0; i<taillePopulation; i++)
    {
        (populationFC+i)->nbrCluster=1; /// le nombre de cluster doit etre initialis� � 1 pour prendre en cons�diration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes viol�es, pour chaque solution
        (populationFC+i)->contrainteViole=0;
        for(j=0;j<nbr_constraint;j++){(populationFC+i)->constraintVector[j]=0;}


        m=nbrNoeuds;
        ///copier le phenotype dans un deuxi�me tableau
        for(j=0; j<nbrNoeuds; j++)
        {
            tabTmp[j]=(populationFC+i)->phenotype[j];
        }
        ///trier le tableau r�sulant pour savoir le nombre de clusters
        /// par exemple : 2 2 3 4 1 1 2 3 ==> 1 1 2 2 2 3 3 4
        /// l'�tape suivante : 1 2 3 4 ==> 4 cluster, et combien de noeuds dans chaque cluster
        /// ceci pour v�rifier le respect des contraintes.

        for(j=0; j<m-1; j++)
        {
            for(k=j+1; k<m; k++)
            {
                if (tabTmp[k] < tabTmp[j])
                {
                    tmp = tabTmp[j];
                    tabTmp[j] = tabTmp[k];
                    tabTmp[k] = tmp;
                }
            }
        }
        ///calculer le nombre des cluser
        for(j=1; j<nbrNoeuds; j++)
        {
            if(tabTmp[j-1] != tabTmp[j])
            {
                (populationFC+i)->nbrCluster++;
            }
        }

        if((populationFC+i)->nbrCluster > max_clusters || (populationFC+i)->nbrCluster < min_clusters)
        {
            (populationFC+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationFC+i)->contrainteViole++;
        }

        /// v�rifier la taille des clusters

        for(j=0; j<m-1; j++)  /// je dois me limitter � m-1 pour que t arrive jusqu'� m
        {
            while(tabTmp[j] == tabTmp[j+1])
            {
                for(k=j; k<m; k++) /// la taille du vecteur change � chaque fois => m et non pas nbrNoeuds
                {
                    tabTmp[k] = tabTmp[k+1];
                }
                m--;
            }
        }
        for(j=0; j<m; j++)
        {
            (populationFC+i)->clustersSize[j]=0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if (tabTmp[j] == (populationFC+i)->phenotype[k])
                {
                    (populationFC+i)->clustersSize[j]=(populationFC+i)->clustersSize[j]+1;
                }
            }
            if((populationFC+i)->clustersSize[j]>max_sizeCluster || (populationFC+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a d�pass� la borne
                (populationFC+i)->constraintVector[1]++;

            }
        }

        if((populationFC+i)->constraintVector[1] != 0)
        {
            (populationFC+i)->contrainteViole++;
        }
   }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstrainteFC(populationFC);
    }

}

*/
