#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../headers/migration.h"
#include "../headers/prog_params.h"
#include "../headers/typeDeclaration.h"
#include "../headers/genetic_alg_PMP.h"
#include "../headers/genetic_alg_DVTC.h"
#include "../headers/genetic_alg_BE.h"
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
extern partitionDVTC *populationDVTC1,*populationDVTC2, *solutionDominanteDVTC,*bestSolutionOverRunsDVTC;
extern partitionBE *populationBE1,*populationBE2, *solutionDominanteBE,*bestSolutionOverRunsBE;

///***************************************************************************
///                     les corps des fonctions
///***************************************************************************

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Migration from DVTC to PMP
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void migrationFromDVTCtoPMP(){
    //qsort(populationPMP1, taillePopulation, sizeof *populationPMP1, compareCroissantFitnessPMP);
    for(int i = 0; i <nbrNoeuds; i++)
        (populationPMP1)->phenotype[i] = solutionDominanteDVTC->genotype[i];

    (populationPMP1)->id = 0;
    (populationPMP1)->coutCoupe = solutionDominanteDVTC->coutCoupe;
    (populationPMP1)->coutCoupeNormalise = solutionDominanteDVTC->coutCoupeNormalise;
    (populationPMP1)->contrainteViole = solutionDominanteDVTC->contrainteViole;
    (populationPMP1)->medians = solutionDominanteDVTC->nbrCluster;

    for(int i = 0; i < (populationPMP1)->medians; i++)
        (populationPMP1)->clusters[i] = solutionDominanteDVTC->clusters[i];

    for(int i = 0; i < nbr_constraint; i++)
        (populationPMP1)->constraintVector[i] = solutionDominanteDVTC->constraintVector[i];
    int sommeTotalCoutCoupe = 0;
    for(int i = 0; i < taillePopulation; i++)
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationPMP1+i)->coutCoupeNormalise;

    float varianceCoutDeCoupeNormalise = 0, moyenneCoutDeCoupeNormalise, sigma, sommeTotalOfExpectdValues = 0;
    int c= 2;

    moyenneCoutDeCoupeNormalise = sommeTotalCoutCoupe / taillePopulation;
    for(int i=0; i<taillePopulation; i++)
    {
        varianceCoutDeCoupeNormalise =
                varianceCoutDeCoupeNormalise + pow(((populationPMP1+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);

    /// calcul de sigma truncation
    (populationPMP1)->expectedValue = (populationPMP1+taillePopulation-1)->coutCoupeNormalise +
                                                        (moyenneCoutDeCoupeNormalise - c*sigma);

    for(int i=0; i< taillePopulation;i++)
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationPMP1+i)->expectedValue;

    for(int i=0; i < taillePopulation ; i++)
        (populationPMP1+i)->fitness = (float)((populationPMP1+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);

    //Encoding the solution using the PMEB encoding
    /*
     * Among the vertices we chose "medians" vertices those having the heigher degree
     * Within their cluster ==>
     */

    for (int i = 0; i < nbrNoeuds; i++) {
        (populationPMP1)->genotype[i] = 0;
    }

    for (int i = 0; i <  (populationPMP1)->medians; i++) {
        int maxDegreeIndex = (populationPMP1+i)->clusters[i].clusterNoueds[0],
            maxDegree = neighborMatrix[maxDegreeIndex][0];
        for (int j = 1; j <  (populationPMP1)->clusters[i].clusterSize; j++) {
            if(maxDegree < neighborMatrix[(populationPMP1)->clusters[i].clusterNoueds[j]][0]) {
                maxDegree = neighborMatrix[(populationPMP1)->clusters[i].clusterNoueds[j]][0];
                maxDegreeIndex = (populationPMP1)->clusters[i].clusterNoueds[j];
            }
        }
        (populationPMP1)->genotype[maxDegreeIndex] = 1;
    }
}


///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Migration from BE to PMP
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void migrationFromBEtoPMP(){
    //printf("\nMigration fron BE to PMP\n");
    qsort(populationPMP1, taillePopulation, sizeof *populationPMP1, compareCroissantFitnessPMP);
    for(int i = 0; i <nbrNoeuds; i++)
        (populationPMP1)->phenotype[i] = solutionDominanteBE->phenotype[i];

    (populationPMP1)->id = 0;
    (populationPMP1)->coutCoupe = solutionDominanteBE->coutCoupe;
    (populationPMP1)->coutCoupeNormalise = solutionDominanteBE->coutCoupeNormalise;
    (populationPMP1)->contrainteViole = solutionDominanteBE->contrainteViole;
    (populationPMP1)->medians = solutionDominanteBE->nbrCluster;

    for(int i = 0; i < (populationPMP1)->medians; i++)
        (populationPMP1)->clusters[i] = solutionDominanteBE->clusters[i];

    for(int i = 0; i < nbr_constraint; i++)
        (populationPMP1)->constraintVector[i] = solutionDominanteBE->constraintVector[i];
    int sommeTotalCoutCoupe = 0;
    for(int i = 0; i < taillePopulation; i++)
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationPMP1+i)->coutCoupeNormalise;

    float varianceCoutDeCoupeNormalise = 0, moyenneCoutDeCoupeNormalise, sigma, sommeTotalOfExpectdValues = 0;
    int c= 2;

    moyenneCoutDeCoupeNormalise = sommeTotalCoutCoupe / taillePopulation;
    for(int i=0; i<taillePopulation; i++)
    {
        varianceCoutDeCoupeNormalise =
                varianceCoutDeCoupeNormalise + pow(((populationPMP1+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);

    /// calcul de sigma truncation
    (populationPMP1)->expectedValue = (populationPMP1+taillePopulation-1)->coutCoupeNormalise +
                                                        (moyenneCoutDeCoupeNormalise - c*sigma);

    for(int i=0; i< taillePopulation;i++)
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationPMP1+i)->expectedValue;

    for(int i=0; i < taillePopulation ; i++)
        (populationPMP1+i)->fitness = (float)((populationPMP1+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);

    //Encoding the solution using the PMEB encoding
    /*
     * Among the vertices we chose "medians" vertices those having the heigher degree
     * Within their cluster ==>
     */

    for (int i = 0; i < nbrNoeuds; i++) {
        (populationPMP1)->genotype[i] = 0;
    }

    for (int i = 0; i <  (populationPMP1)->medians; i++) {
        int maxDegreeIndex = (populationPMP1)->clusters[i].clusterNoueds[0],
            maxDegree = neighborMatrix[maxDegreeIndex][0];
        for (int j = 1; j <  (populationPMP1)->clusters[i].clusterSize; j++) {
            if(maxDegree < neighborMatrix[(populationPMP1)->clusters[i].clusterNoueds[j]][0]) {
                maxDegree = neighborMatrix[(populationPMP1)->clusters[i].clusterNoueds[j]][0];
                maxDegreeIndex = (populationPMP1)->clusters[i].clusterNoueds[j];
            }
        }
        (populationPMP1)->genotype[maxDegreeIndex] = 1;
    }
}

///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
///FFFFFFFFFFF migrationFromBEtoDVTC FFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
void migrationFromBEtoDVTC(){
    //printf("\nMigration fron BE to DVTC\n");
    qsort(populationDVTC1, taillePopulation, sizeof *populationDVTC1, compareCroissantFitnessDVTC);
    for(int i = 0; i <nbrNoeuds; i++)
        (populationDVTC1)->genotype[i] = solutionDominanteBE->phenotype[i];

    (populationDVTC1)->id = taillePopulation - 1;
    (populationDVTC1)->coutCoupe = solutionDominanteBE->coutCoupe;
    (populationDVTC1)->coutCoupeNormalise = solutionDominanteBE->coutCoupeNormalise;
    (populationDVTC1)->contrainteViole = solutionDominanteBE->contrainteViole;
    (populationDVTC1)->nbrCluster = solutionDominanteBE->nbrCluster;

    for(int i = 0; i <(populationDVTC1)->nbrCluster; i++)
        (populationDVTC1)->clusters[i] = solutionDominanteBE->clusters[i];

    for(int i = 0; i < nbr_constraint; i++)
        (populationDVTC1)->constraintVector[i] = solutionDominanteBE->constraintVector[i];
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
    (populationDVTC1)->expectedValue = (populationDVTC1)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);

    for(int i=0; i< taillePopulation;i++)
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationDVTC1+i)->expectedValue;

    for(int i=0; i<taillePopulation ; i++)
    (populationDVTC1+i)->fitness = (float)((populationDVTC1+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);

}

///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
///FFFFFFFFFFF migrationFromBEtoDVTC FFFFFFFFFFFFFFF
///FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE
void migrationFromPMPtoDVTC(){
    //printf("\nMigration fron PMP to DVTC\n");
    qsort(populationDVTC1, taillePopulation, sizeof *populationDVTC1, compareCroissantFitnessDVTC);
    for(int i = 0; i <nbrNoeuds; i++)
        (populationDVTC1)->genotype[i] = solutionDominantePMP->phenotype[i];

    (populationDVTC1)->id = taillePopulation - 1;
    (populationDVTC1)->coutCoupe = solutionDominantePMP->coutCoupe;
    (populationDVTC1)->coutCoupeNormalise = solutionDominantePMP->coutCoupeNormalise;
    (populationDVTC1)->contrainteViole = solutionDominantePMP->contrainteViole;
    (populationDVTC1)->nbrCluster = solutionDominantePMP->medians;

    for(int i = 0; i < solutionDominantePMP->medians; i++)
        (populationDVTC1)->clusters[i] = solutionDominantePMP->clusters[i];

    for(int i = 0; i < nbr_constraint; i++)
        (populationDVTC1)->constraintVector[i] = solutionDominantePMP->constraintVector[i];
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
    (populationDVTC1)->expectedValue = (populationDVTC1)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);

    for(int i=0; i< taillePopulation;i++)
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationDVTC1+i)->expectedValue;

    for(int i=0; i<taillePopulation ; i++)
        (populationDVTC1+i)->fitness = (float)((populationDVTC1+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);

}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Migration from DVTC to BE
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void migrationFromDVTCtoBE(){
    //printf("\nMigration fron DVTC to BE\n");
    qsort(populationBE1, taillePopulation, sizeof *populationBE1, compareCroissantFitnessBE);
    for(int i = 0; i <nbrNoeuds; i++)
        (populationBE1)->phenotype[i] = solutionDominanteDVTC->genotype[i];

    (populationBE1)->id = 0;
    (populationBE1)->coutCoupe = solutionDominanteDVTC->coutCoupe;
    (populationBE1)->coutCoupeNormalise = solutionDominanteDVTC->coutCoupeNormalise;
    (populationBE1)->contrainteViole = solutionDominanteDVTC->contrainteViole;
    (populationBE1)->nbrCluster = solutionDominanteDVTC->nbrCluster;

    for(int i = 0; i <(populationBE1+taillePopulation-1)->nbrCluster; i++)
        (populationBE1)->clusters[i] = solutionDominanteDVTC->clusters[i];

    for(int i = 0; i < nbr_constraint; i++)
        (populationBE1)->constraintVector[i] = solutionDominanteDVTC->constraintVector[i];
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
    (populationBE1)->expectedValue = (populationBE1+taillePopulation-1)->coutCoupeNormalise +
                                                        (moyenneCoutDeCoupeNormalise - c*sigma);

    for(int i=0; i< taillePopulation;i++)
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationBE1+i)->expectedValue;

    for(int i=0; i<taillePopulation ; i++)
        (populationBE1+i)->fitness = (float)((populationBE1+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);

    //Encoding the solution using the binary encoding
    int edgeIndex;
    for(int i = 0; i < nbrArretes; i++ ) (populationBE1)->genotype[i] = 1;
    for(int i = 0; i < (populationBE1)->nbrCluster; i++) {
        for (int j = 0; j < (populationBE1)->clusters[i].clusterSize - 1; j++) {
            for (int k = j + 1; k < (populationBE1)->clusters[i].clusterSize; k++) {
                int nd = (populationBE1)->clusters[i].clusterNoueds[j];
                int na = (populationBE1)->clusters[i].clusterNoueds[k];
                if (fluxMatrix[nd][na] > 0) {
                    edgeIndex = searchForEdgeIndex(nd, na);
                }

                if (edgeIndex)
                    (populationBE1)->genotype[edgeIndex] = 0;
            }
        }
    }
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Migration from PMP to BE
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void migrationFromPMPtoBE(){
    //printf("\nMigration fron PMP to BE\n");
    //uTest_uTest_uTest_uTest_uTest_uTest_uTest_uTest_uTest_uTest_uTest

    //uTest_uTest_uTest_uTest_uTest_uTest_uTest_uTest_uTest_uTest_uTest
    qsort(populationBE1, taillePopulation, sizeof *populationBE1, compareCroissantFitnessBE);
    for(int i = 0; i < nbrNoeuds; i++)
        (populationBE1)->phenotype[i] = solutionDominantePMP->phenotype[i];

    (populationBE1)->id = 0;
    (populationBE1)->coutCoupe = solutionDominantePMP->coutCoupe;
    (populationBE1)->coutCoupeNormalise = solutionDominantePMP->coutCoupeNormalise;
    (populationBE1)->contrainteViole = solutionDominantePMP->contrainteViole;
    (populationBE1)->nbrCluster = solutionDominantePMP->medians;

    for(int i = 0; i < solutionDominantePMP->medians; i++)
        (populationBE1)->clusters[i] = solutionDominantePMP->clusters[i];

    for(int i = 0; i < nbr_constraint; i++)
        (populationBE1)->constraintVector[i] = solutionDominantePMP->constraintVector[i];
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
    (populationBE1)->expectedValue = (populationBE1+taillePopulation-1)->coutCoupeNormalise +
                                                        (moyenneCoutDeCoupeNormalise - c*sigma);

    for(int i=0; i< taillePopulation;i++)
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationBE1+i)->expectedValue;

    for(int i=0; i<taillePopulation ; i++)
        (populationBE1+i)->fitness = (float)((populationBE1+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);

    //Encoding the solution using the binary encoding
    int edgeIndex;
    for(int i = 0; i < nbrArretes; i++ ) (populationBE1)->genotype[i] = 1;
    for(int i = 0; i < (populationBE1)->nbrCluster; i++) {
        for (int j = 0; j < (populationBE1)->clusters[i].clusterSize - 1; j++) {
            for (int k = j + 1; k < (populationBE1)->clusters[i].clusterSize; k++) {
                int nd = (populationBE1)->clusters[i].clusterNoueds[j];
                int na = (populationBE1)->clusters[i].clusterNoueds[k];
                if (fluxMatrix[nd][na] > 0) {
                    edgeIndex = searchForEdgeIndex(nd, na);
                }

                if (edgeIndex)
                    (populationBE1)->genotype[edgeIndex] = 0;
            }
        }
    }
}