#ifndef FRACTIONALENCODING_H_INCLUDED
#define FRACTIONALENCODING_H_INCLUDED


void generatePopulationFC(partitionFC* populationFC,int indiceFirstElt);
void affichePopulationFC(partitionFC* populationFC);
void calculCoutCoupeEtFitnessFC(partitionFC* populationFC);
void naturalSelectionFC(partitionFC* populationFC1,partitionFC* populationFC2);
void crossOverFC(partitionFC* populationFC1, partitionFC* populationFC2);
void calculPhenotypeFC(partitionFC* populationFC);
int findTheBestSolutionFC(partitionFC *populationFC);
void mutationFC(partitionFC* populationFC);
float testerLaSommeDesFitnessFC(partitionFC* populationFC);
void displayTheBestSolutionFC(partitionFC* solutionDominante);
void checkContrainstAndFitnessPenalizationFC(partitionFC *populationFC);
void writeSolutionInFileFC(partitionFC *populationFC,FILE *outputFilePop,int iteration);
void writeBestSolutionInFileFC(partitionFC *bestSolution, FILE *outputFile,int iteration);
void fractionalEncoding(int nbrGeneration,FILE *outputFileFC,FILE *outputFilePopFC,FILE *outputOptimalSolutionFileFC,partitionFC *populationFC1,partitionFC *populationFC2,partitionFC *solutionDominante );
void checkCohabitationAndNonCohabitationConstrainteFC(partitionFC *populationFC);
///void writeOptimalSolutionInFileFC(partitionFC *solutionDominante,FILE *outputOptimalSolutionFileFC);
void writeOptimalSolutionInFileFC(partitionFC *solutionDominante,FILE* outputOptimalSolutionFileFC,
                                  int nbrRun, int bestSolutionIteration, float runTime,int ES);
int compareCroissantFitnessFC (void const *a, void const *b);
double drand48ForWindowsFC(int v1,int v2);



#endif // FRACTIONALENCODING_H_INCLUDED
