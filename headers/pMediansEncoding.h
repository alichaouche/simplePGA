#ifndef PMEDIANSPROBLEMBINARYENCODING_H_INCLUDED
#define PMEDIANSPROBLEMBINARYENCODING_H_INCLUDED


void generatePopulationPMP(partitionPMP* populationPMP, int indiceFirstElt);
void affichePopulationPMP(partitionPMP* populationPMP);
void calculerMediansNonMediansVertices(partitionPMP* populationPMP);
void affectationDesNonMedians(partitionPMP* populationPMP);
void calculeDeCutSizeEtFitnessPMP(partitionPMP* populationPMP);
void naturalSelectionPMP(partitionPMP* populationPMP1,partitionPMP* populationPMP2);
void crossOverPMP(partitionPMP* populationPMP1, partitionPMP* populationPMP2);
int findTheBestSolutionPMP(partitionPMP *populationPMP);
void mutationPMP(partitionPMP* populationPMP);
float testerLaSommeDesFitnessPMP(partitionPMP* populationPMP);
void displayTheBestSolutionPMP(partitionPMP* solutionDominante);
void checkContrainstAndFitnessPenalizationPMP(partitionPMP *populationPMP);
void writeSolutionInFilePMP(partitionPMP *populationPMP, FILE *outputFilePop,int iteration);
void writeBestSolutionInFilePMP(partitionPMP *solutionDominante, FILE *outputFile,int iteration);
void pMedianEncoding(int nbrGeneration,FILE *outputFilePMP,FILE *outputFilePopPMP,FILE *outputOptimalSolutionFilePMP);
void checkContrainstAndFitnessPenalizationPMP(partitionPMP *populationPMP);
///void writeOptimalSolutionInFilePMP(partitionPMP *solutionDominante,FILE *outputOptimalSolutionFilePMP);
int compareCroissantFitnessPMP (void const *a, void const *b);
void writeOptimalSolutionInFilePMP(partitionPMP *solutionDominante,FILE* outputOptimalSolutionFilePMP,
                                   int nbrRun, int bestSolutionIteration, float runTime,int ES);

void pMedianEncodingBaseSurClusterSize(int nbrGeneration ,FILE *outputFilePMP,FILE *outputFilePopPMP,FILE *outputOptimalSolutionFilePMP,partitionPMP *populationPMP1,partitionPMP *populationPMP2,partitionPMP *solutionDominante);
void affectationDesNonMediansBaseSurClusterSize(partitionPMP* populationPMP);

void affectationDesNonMediansParArretes(partitionPMP* populationPMP);
void pMedianEncodingAffectationParArretes(int nbrGeneration ,FILE *outputFilePMP,FILE *outputFilePopPMP,FILE *outputOptimalSolutionFilePMP,partitionPMP *populationPMP1,partitionPMP *populationPMP2,partitionPMP *solutionDominante);

void checkCohabitationAndNonCohabitationConstraintePMP(partitionPMP *populationPMP);
#endif // PMEDIANSPROBLEMBINARYENCODING_H_INCLUDED
