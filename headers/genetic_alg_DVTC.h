#ifndef VERTEX_TO_CLUSTER_H_INCLUDED
#define VERTEX_TO_CLUSTER_H_INCLUDED





void generatePopulationDVTC(partitionDVTC* populationDVTC, int indiceFirstElt);

void affichePopulationDVTC(partitionDVTC* populationDVTC);
void calculCoutCoupeEtFitnessDVTC(partitionDVTC* populationDVTC);
void crossOverDVTC(partitionDVTC* populationDVTC1, partitionDVTC* populationDVTC2);
void naturalSelectionDVTC(partitionDVTC* populationDVTC1,partitionDVTC* populationDVTC2);
int findTheBestSolutionDVTC(partitionDVTC *populationDVTC);
void mutationDVTC(partitionDVTC* populationDVTC);
float testerLaSommeDesFitnessDVTC(partitionDVTC* populationDVTC);
void displayTheBestSolutionDVTC(partitionDVTC* solutionDominante);
void checkContrainstAndFitnessPenalizationDVTC(partitionDVTC *populationDVTC);


void writeSolutionInFileDVTC(partitionDVTC *populationDVTC,FILE *outputFilePop,int iteration);
void writeBestSolutionInFileDVTC(partitionDVTC *bestSolution, FILE *bestSolutionsOverIterationDVTC,int iteration);
void vertexToClusterEncoding(int numeroGraphe, metrics *metricsDVTC);
void checkCohabitationAndNonCohabitationConstrainteDVTC(partitionDVTC *populationDVTC);
///void writeOptimalSolutionInFileDVTC(partitionDVTC *solutionDominante,FILE *outputOptimalSolutionFileDVTC);
int compareCroissantFitnessDVTC (void const *a, void const *b);
void writeOptimalSolutionInFileDVTC(partitionDVTC *solutionDominante,FILE* outputOptimalSolutionFileDVTC,
                                  int nbrRun, metrics metricsDVTC);
void solutionsReproductionDVTC(partitionDVTC* populationDVTC);
void migrationFromBEtoDVTC();
int compareDecroissantFitnessDVTC (void const *a, void const *b);

#endif // VERTEX_TO_CLUSTER_H_INCLUDED
