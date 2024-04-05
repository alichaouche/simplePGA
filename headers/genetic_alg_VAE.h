#ifndef SORTEDCUTBASEDENCODING_H_INCLUDED
#define SORTEDCUTBASEDENCODING_H_INCLUDED



///=========================== Définitions des entêtes des fonction =======================
void vertexAssignmentEncoding(int nbrGeneration, FILE* outputFileSBGNE, FILE* outputFilePopSBGNE,FILE *outputOptimalSolutionFileSBGNE,partitionSBGNE *populationSBGNE1,partitionSBGNE *populationSBGNE2, partitionSBGNE *solutionDominante);
void genererPopulationInitialeSBGNE(partitionSBGNE *populationSBGNE , int indiceFirstElt);
void generateBinarySolutionSBGNE(partitionSBGNE* populationSBGNE, int indiceIndividu);
void homoginiserLesSolutionBinaire(partitionSBGNE* populaiton, int indiceIndividu);
void conversionSolutionBinaireToEntiereSBGNE(partitionSBGNE* populationSBGNE,int indiceIndividu);
void trieGenotypeEntieSBGNE(partitionSBGNE* populationSBGNE, int indiceIndividu);
int existanceDeSolutionSBGNE(partitionSBGNE* populationSBGNE, int indexNewSolution);
void getPartitionFromSolutionSBGNE(partitionSBGNE *populationSBGNE);
void calculCoutCoupeEtFitnessSBGNE(partitionSBGNE* populationSBGNE);
void checkContrainstAndFitnessPenalizationSBGNE(partitionSBGNE *populationSBGNE);
void checkCohabitationAndNonCohabitationConstrainteSBGNE(partitionSBGNE *populationSBGNE);
void naturalSelectionSBGNE(partitionSBGNE* populationSBGNE1,partitionSBGNE* populationSBGNE2);
void crossOverSBGNE(partitionSBGNE* populationSBGNE1, partitionSBGNE* populationSBGNE2);
int findTheBestSolutionSBGNE(partitionSBGNE *populationSBGNE);
void mutationSBGNE(partitionSBGNE* populationSBGNE);
float testerLaSommeDesFitnessSBGNE(partitionSBGNE* populationSBGNE);
void displayTheBestSolutionSBGNE(partitionSBGNE* solutionDominante);
void writeSolutionInFileSBGNE(partitionSBGNE *populationSBGNE, FILE *outputFilePop,int iteration);
void writeBestSolutionInFileSBGNE(partitionSBGNE *solutionDominante, FILE *outputFile,int iteration);
void affichePopulationSBGNE(partitionSBGNE* populationSBGNE);


void genererPopulationInitialeRandomlySBGNE(partitionSBGNE *populationSBGNE, int indiceFirstElt);
///void writeOptimalSolutionInFileSBGBE(partitionSBGNE *solutionDominante,FILE *outputOptimalSolutionFileSBGNE);
void writeOptimalSolutionInFileSBGNE(partitionSBGNE *solutionDominante,FILE* outputOptimalSolutionFileSBGNE,
                                     int nbrRun, int bestSolutionIteration, float runTime, int ES);

int compareCroissantFitnessSBGNE (void const *a, void const *b);

#endif // SORTEDCUTBASEDENCODING_H_INCLUDED
