#ifndef SORTEDVTC_H_INCLUDED
#define SORTEDVTC_H_INCLUDED


void generatePopulationSVTC(partitionSVTC* populationSVTC, int indiceFirstElt);
void affichePopulationSVTC(partitionSVTC* populationSVTC);
void calculCoutCoupeEtFitnessSVTC(partitionSVTC* populationSVTC);
void crossOverSVTC(partitionSVTC* populationSVTC1, partitionSVTC* populationSVTC2);
void naturalSelectionSVTC(partitionSVTC* populationSVTC1,partitionSVTC* populationSVTC2);
int findTheBestSolutionSVTC(partitionSVTC *populationSVTC);
void mutationSVTC(partitionSVTC* populationSVTC);
float testerLaSommeDesFitnessSVTC(partitionSVTC* populationSVTC);
void displayTheBestSolutionSVTC(partitionSVTC* solutionDominante);
void checkContrainstAndFitnessPenalizationSVTC(partitionSVTC *populationSVTC);

void agencementDesGenotype(partitionSVTC* populationSVTC);
void writeSolutionInFileSVTC(partitionSVTC *populationSVTC,FILE *outputFilePop,int iteration);
void writeBestSolutionInFileSVTC(partitionSVTC *bestSolution, FILE *outputFile,int iteration);
void sortedVertexToClusterEncoding(int nbrGeneration ,FILE *outputFileSVTC,FILE *outputFilePopSVTC,FILE *outputOptimalSolutionFileSVTC, partitionSVTC *populationSVTC1,partitionSVTC *populationSVTC2,partitionSVTC *solutionDominante);
void checkCohabitationAndNonCohabitationConstrainteSVTC(partitionSVTC *populationSVTC);
///void writeOptimalSolutionInFileSVTC(partitionSVTC *solutionDominante,FILE *outputOptimalSolutionFileSVTC);
int compareCroissantFitnessSVTC (void const *a, void const *b);
void writeOptimalSolutionInFileSVTC(partitionSVTC *solutionDominante,FILE* outputOptimalSolutionFileSVTC,
                                    int nbrRun, int bestSolutionIteration, float runTime, int ES);


#endif // SORTEDVTC_H_INCLUDED
