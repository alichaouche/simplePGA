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
void checkContrainstsPMP(partitionPMP *populationPMP);
void checkCohabitationAndNonCohabitationConstraintePMP(partitionPMP *populationPMP);
void writeSolutionInFilePMP(partitionPMP *populationPMP, FILE *outputFilePop,int iteration);
void writeBestSolutionInFilePMP_old(partitionPMP *solutionDominante, FILE *outputFile,int iteration);
void pMedianEncoding(int nbrGeneration,FILE *outputFilePMP,FILE *outputFilePopPMP,FILE *outputOptimalSolutionFilePMP);
///void writeOptimalSolutionInFilePMP(partitionPMP *solutionDominante,FILE *outputOptimalSolutionFilePMP);
int compareCroissantFitnessPMP (void const *a, void const *b);
int compareDecroissantFitnessPMP(void const *a, void const *b);
void writeOptimalSolutionInFilePMP(partitionPMP *solutionDominante,FILE* optimalSolutionFilePMP,
                                  int nbrRun, metrics metricsPMP);
void pMedianEncodingBaseSurClusterSize(int nbrGeneration, metrics *metricsPMP);
void affectationDesNonMediansBaseSurClusterSize(partitionPMP* populationPMP);

void affectationDesNonMediansParArretes(partitionPMP* populationPMP);
void pMedianEncodingAffectationParArretes(int numeroGraphe, metrics *metricsPMP);

void realCutSizeForTheBestSolution(partitionPMP *solutionDominante);
void writeBestSolutionInFilePMP(partitionPMP *solutionDominante, FILE *bestSolutionsOverIterationPMP,int iteration);
void solutionsReproductionPMP(partitionPMP* populationPMP);


#endif // PMEDIANSPROBLEMBINARYENCODING_H_INCLUDED
