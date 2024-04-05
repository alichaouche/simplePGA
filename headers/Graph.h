//
// Created by Ali on 01-Apr-24.
//

#ifndef SIMPLEPGA_GRAPH_H
#define SIMPLEPGA_GRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include "typeDeclaration.h"
#include "prog_params.h"


void afficherMatriceFlux();
void connecterGraphe();
int isConnex();
int calculerNombreArrete();
void calculateEdgeVector();
void DisplayFlowVector();
int calculSommeTotalFlux();
void maxMinWeight(int *min, int *max);
int graphConnexion();
void graphConnectedness(int voisin[],int nonVoisin[],int v,int nv);
void deleteEdgesLowerThanThreshold(int seuil);
int getTheConnexComponentNumberBE(int vecteurBinaireDesAretes[], int vectorSize );
void newDeleteEdgesLowerThanThreshold(int edge[],int fixedNbrArretes);
void readFlowMatrixAndNbrNoeudsFromFlowMatrixFile(FILE *inputFile);
void readFlowMatrixAndNbrNoeudsFromEdgeListFile(FILE *edgeListFile );
void updateDeletedEdges(int *aretesSupprimees, int fixedNbrArretes, edge *deletedEdges, int *nbrAretesSupprimees );
int isFlowMatrixSymetric();
void displayEdgeVector();

#endif //SIMPLEPGA_GRAPH_H
