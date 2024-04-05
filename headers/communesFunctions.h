//
// Created by Ali on 01-Apr-24.
//

#ifndef SIMPLEPGA_COMMUNESFUNCTIONS_H
#define SIMPLEPGA_COMMUNESFUNCTIONS_H
#include <stdarg.h>

typedef unsigned long ul;

double drand48ForWindows(int v1, int v2);
int rnd(int v1, int v2);
void mon_sleep( int nbr_seconds );
int compareTrieCroissant (void const *a, void const *b);
int compareTrieDecroissant(void const *a, void const *b);

void openingFile(int numeroGraphe);
void writeDetailsProblemInFile(FILE* detailsFile, int seuil, edge deletedEdges[],int nbrDeletedEdge);
void edgeListToAdjMatrix();
void getCohabitationCouples(int maxWeight);

void memoryAllocation(int, ...);
#endif //SIMPLEPGA_COMMUNESFUNCTIONS_H
