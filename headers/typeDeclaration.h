//
// Created by Ali on 01-Apr-24.
//

#ifndef SIMPLEPGA_TYPEDECLARATION_H
#define SIMPLEPGA_TYPEDECLARATION_H

    typedef unsigned long int ul;
    ///############ TYPE EDGE ################
    typedef struct
    {
        int numeroEdge;
        int nouedDepart;
        int nouedArrive;
        int poids;
    } edge;

    ///############ TYPE CLUSTER ################
    typedef struct
    {
        unsigned int clusterId;
        unsigned int clusterSize;
        unsigned int clusterNoueds[500];
    } clusterInfo;

    ///############# TYPE PARTITIONBE ############

    typedef struct
    {
        unsigned int id;
        //int *genotype;
        char genotype[100000];
        //int *phenotype;
        unsigned int phenotype[10000];
        ///************************************
        unsigned int contrainteViole;
        unsigned int nbrCluster;
        //clusterInfo *clusters;
        clusterInfo clusters[1000];
        unsigned int coutCoupeNormalise;
        char constraintVector[4];
        ///************************************
        unsigned int coutCoupe;
        float fitness;
        float expectedValue;
    } partitionBE;
///=============================================================================================


typedef struct
    {
        int id;
        float genotype[10000];
        int phenotype[10000];
        int contrainteViole;
        int nbrCluster;
        ul clustersSize[1000];
        int coutCoupeNormalise;
        int constraintVector[10];
        int coutCoupe;
        float fitness;
        float expectedValue;
    } partitionFC;
    ///=============================================================================================

    typedef struct
    {
        int id;
        int genotype[10000];
        int phenotype[10000];
        int mediansVertices[10000];
        int medians;
        int nonMediansVertices[10000];
        int nonMedians;
        int contrainteViole;
        ul clustersSize[1000];
        int coutCoupeNormalise;
        int constraintVector[10];
        int coutCoupe;
        float fitness;
        float expectedValue;
    } partitionPMP;
    ///=============================================================================================

    typedef struct
    {
        int id;
        int genotype[1000];
        int contrainteViole;
        int nbrCluster;
        //clusterInfo *clusters;
        clusterInfo clusters[1000];
        int clustersSize[1000];
        int coutCoupeNormalise;
        int constraintVector[10];
        int coutCoupe;
        float fitness;
        float expectedValue;
    } partitionDVTC;
    ///=============================================================================================

    typedef struct
    {
        int id;
        int genotypeEntier[100]; ///je vais utiliser ce vecteur pour générer des solution sans redandance.
        char genotype[100000];
        int phenotype[1000];
        int contrainteViole;
        int nbrCluster;
        ul clustersSize[1000];
        int coutCoupeNormalise;
        int constraintVector[10];
        int coutCoupe;
        float fitness;
        float expectedValue;
    } partitionVAE;

    ///==================================================================================================
    typedef struct
    {
        int id;
        int genotype[10000];
        int contrainteViole;
        int nbrCluster;
        ul clustersSize[1000];
        int coutCoupeNormalise;
        int constraintVector[10];
        int coutCoupe;
        float fitness;
        float expectedValue;
    } partitionSVTC;



    typedef struct {
        float runTime;
        float ES;
    }metrics;
#endif //SIMPLEPGA_TYPEDECLARATION_H
