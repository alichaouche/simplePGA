//
// Created by Ali on 01-Apr-24.
//

#ifndef SIMPLEPGA_GENETIC_ALG_BE_H
#define SIMPLEPGA_GENETIC_ALG_BE_H


void generatePopulationBE(partitionBE* populationBE, int indexFirstIndividual);
void generatePopulationRandomlyBE(partitionBE* populationBE);
void affichePopulationBE(partitionBE* populationBE);
void calculCoutCoupeEtFitnessWithFlowVectorBE(partitionBE* populationBE);

void naturalSelectionBE(partitionBE* populationBE1,partitionBE* populationBE2);
void crossOverBE(partitionBE* populationBE1, partitionBE* populationBE2);

int findTheBestSolutionBE(partitionBE *populationBE);
void mutationBE(partitionBE* populationBE);
float testerLaSommeDesFitnessBE(partitionBE* populationBE);
void displayTheBestSolutionBE(partitionBE* solutionDominante);
void getPartitionFromSolutionBE(partitionBE *populationBE);
void checkContrainstAndFitnessPenalizationBE(partitionBE *populationBE);
void writeSolutionInFileBE(partitionBE *populationBE, FILE *outputFilePop,int iteration);
void writeBestSolutionInFileBE(partitionBE *solutionDominante, FILE *outputFile,int iteration);
int  binaryEncoding(int numeroGraphe, metrics *metricsBE);
void checkCohabitationAndNonCohabitationConstrainteBE(partitionBE *populationBE);

void getPartitionFromSolutionWithoutRepetitionBE(partitionBE *populationBE);

void writeOptimalSolutionInFileBE(partitionBE *solutionDominante,FILE* optimalSolutionFileBE,
                                  int nbrRun, metrics metricsBE);
int compareCroissantFitnessBE (void const *a, void const *b);

int claculateThresholdsQuartiles(int *threshold, int thresholdQuartiles[3], int nbrThreshold);
void solutionsReproduction(partitionBE* populationBE);

int searchForEdgeIndex(int nd, int na);




#endif //SIMPLEPGA_GENETIC_ALG_BE_H
