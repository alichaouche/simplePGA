//
// Created by Ali on 01-Apr-24.
//
#include "../headers/Graph.h"
///############################################
// Global Variables.
///############################################
extern int  nbrNoeuds,
        nbrArretes,
        fixedNbrArretes,
        fluxMatrix[tailleMax][tailleMax],
        neighborMatrix[tailleMax][tailleMax],
        buckUpFluxMatrix[tailleMax][tailleMax],
        *fluxVector,
        nbrCohabitationConstraintes,
        nbrNonCohabitationConstraintes;

extern edge *edgeVector,*fixedEdgeVector,cohabitationConstraintes[100],nonCohabitationConstraintes[100];
extern neighbors *neighborsVector;
///***********************************************

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Display the graph's flow matrix
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void afficherMatriceFlux()
{
    int i,j;
    for (i = 0; i < nbrNoeuds; i++)
    {
        ///for (j = i+1; j < n; j++) {
        for (j = 0; j < nbrNoeuds; j++)
        {
            if(i==j)
            {
                printf("%d|\t",fluxMatrix[i][j]);
            }
            else
            {
                printf("%d \t",fluxMatrix[i][j]);
            }
        }
        printf("\n");
    }
}
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Connect the graph if it isn't
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void connecterGraphe()
{
    int i;
    for (i = 0; i < nbrNoeuds-1; i++)
    {
        if(fluxMatrix[i][i+1] == -1)
        {
            fluxMatrix[i][i+1] = 0;
            fluxMatrix[i+1][i] = 0;
        }
    }
}
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Calculate the number of edges
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int calculerNombreArrete()
{
    int i,j;
    int nbrArretes = 0,nbrArretesTmp=0;

    for(i=0; i<nbrNoeuds-1; i++)
    {
        for(j=i+1; j<nbrNoeuds; j++)
        {
            if(fluxMatrix[i][j]!= -1 && i!=j)
            {
                nbrArretes++;
            }
        }
    }
    return nbrArretes;
}
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Calculate the EdgeVector
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void calculateEdgeVector(){
    int nbrArretesTmp = 0,i,j;
    for(i=0; i<nbrNoeuds-1; i++)
    {
        for(j=i+1; j<nbrNoeuds; j++)
        {
            if(fluxMatrix[i][j]!= -1 && i!=j)
            {
                edgeVector[nbrArretesTmp].numeroEdge = nbrArretesTmp;
                edgeVector[nbrArretesTmp].nouedDepart = i;
                edgeVector[nbrArretesTmp].nouedArrive = j;
                edgeVector[nbrArretesTmp].poids = fluxMatrix[i][j];
                nbrArretesTmp++;
            }
        }
    }
}


///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Display the flow vector
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void DisplayFlowVector()
{
    int i;
    for(int i=0; i<nbrArretes; i++)
    {
        printf("%d ",fluxVector[i]);
    }
    printf("\n");

}
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Calculate the total weight in the gr
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int calculSommeTotalFlux()
{
    int i,j,sommeFlux=0;
    for(i=0; i<nbrNoeuds; i++)
    {
        for(j=i+1; j<nbrNoeuds; j++)
        {
            if(buckUpFluxMatrix[i][j] > 0 )
            {
                sommeFlux+=buckUpFluxMatrix[i][j];
            }
        }
    }
    return sommeFlux;
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// min and max weight in the graph
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void maxMinWeight(int *min, int *max)
{
    int i,j;
    *max=buckUpFluxMatrix[0][0];
    for(i=0; i<nbrNoeuds-1; i++)
    {
        for(j=i+1; j<nbrNoeuds; j++)
        {
            if(*max < buckUpFluxMatrix[i][j]) *max = buckUpFluxMatrix[i][j];
        }
    }
    *min = *max;
    for(i=0; i<nbrNoeuds-1; i++)
    {
        for(j=i+1; j<nbrNoeuds; j++)
        {
            if(*min > buckUpFluxMatrix[i][j] && buckUpFluxMatrix[i][j] >= 0)
                *min = buckUpFluxMatrix[i][j];
        }
    }
}
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Check the connectedness of the graph
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int graphConnexion(){
    for(int i = 0; i< nbrNoeuds-1; i++){
        int connected=0;
        for(int j=i+1; j < nbrNoeuds; j++){
            if(fluxMatrix[i][j] >= 0){
                connected = 1;
                break;
            }
        }
        if(!connected )
            fluxMatrix[i][i+1] =0;
    }
}


///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Delete edges whose weight is lower to the threshold
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int isConnex(){
    int connected=0;
    for(int i = 0; i< nbrNoeuds-1; i++){
        connected=0;
        for(int j=i+1; j < nbrNoeuds; j++){
            if(fluxMatrix[i][j] >= 0){
                connected = 1;
                break;
            }
        }
        if(!connected )
            return 0;
    }
    return 1;
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Delete edges whose weight is lower to the threshold
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void deleteEdgesLowerThanThreshold(int seuil)
{
///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
/// Keep the edges between CC vertices and NCC vertices
///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    int i,j,k,buckUpWeight,cohabitaion = 0;
    for(i=0; i<nbrNoeuds; i++)
    {
        for(j=0; j<nbrNoeuds; j++)
        {
            if(fluxMatrix[i][j]<= seuil && i!=j && fluxMatrix[i][j] >= 0)
            {
                buckUpWeight = fluxMatrix[i][j];
                fluxMatrix[i][j]=-1;
                fluxMatrix[j][i]=-1;

                for(k=0; k<nbrCohabitationConstraintes; k++)
                {
                    if((cohabitationConstraintes[k].nouedDepart==i  &&
                        cohabitationConstraintes[k].nouedArrive == j)
                       ||
                       (cohabitationConstraintes[k].nouedDepart==j  &&
                        cohabitationConstraintes[k].nouedArrive == i))
                    {
                        cohabitaion = 1;
                        break;
                    }
                }

                if(!isConnex() || cohabitaion)
                {
                    fluxMatrix[i][j]=buckUpWeight;
                    fluxMatrix[j][i]=buckUpWeight;
                }
            }
        }
    }
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Calculate the nbr of connected componentes
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
int getTheConnexComponentNumberBE(int vecteurBinaireDesAretes[], int vectorSize )
{
    ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    ///vecteurBinaireDesAretes <=> aretesSupprimees
    ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    int i,j,k,nd,na,numeroClustre; /// nd = noued de d�part , na = noued d'arriv�
    int phenotype[nbrNoeuds],clustersSize[nbrNoeuds],nbrConnectedComponent;
    ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    /// initialisation des ph�notype = partition par les numero des noeuds correspondant
    ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

    ///########## INITIALISATION DE PHENOTYPE ###########
    for(i=0; i<nbrNoeuds; i++) phenotype[i] = i;


    ///########## DETERMINATION DE LA PARTITION #########
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    /// vectorSize = fixedNbrAretes
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    for(i=0; i<vectorSize; i++)
    {
        ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
        /// DANS CETTE SITUATION ON A 4 CAS DE FIGURES :
        /// [nd AFFECTE ET na AFFECTE] OU [nd AFFECTE ET na NOT AFFECTE]
        /// [na NOT AFFECTE ET nd AFFECTE] OU [na NOT AFFECTE ET nd NOT AFFECTE]
        ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
        if(vecteurBinaireDesAretes[i]==0)
        {
            ///IIIIIIIIIIIIIIIIIIIIIII
            /// TOUJOURS ON A nd < na
            ///IIIIIIIIIIIIIIIIIIIIIII
            nd = (fixedEdgeVector+i)->nouedDepart;
            na = (fixedEdgeVector+i)->nouedArrive;

            if(phenotype[nd]== nd  && phenotype[na]== na) phenotype[na] = phenotype[nd];
            else if(phenotype[nd]!= nd  && phenotype[na]== na) phenotype[na] = phenotype[nd];

            else if(phenotype[nd]== nd  && phenotype[na]!= na)phenotype[nd] = phenotype[na];
            else if(phenotype[nd]!= nd  && phenotype[na]!= na)
            {
                if(phenotype[nd] < phenotype[na])
                {
                    ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                    /// AFFECTER TOUTES LES OCCURENCES DE NA A ND
                    ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                    for(j=0; j<nbrNoeuds; j++)
                    {
                        if(phenotype[j] == phenotype[na])
                            phenotype[j] = phenotype[nd];
                    }
                }
                else if(phenotype[na] < phenotype[nd])
                {
                    ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                    /// AFFECTER TOUTES LES OCCURENCES DE ND A NA
                    ///IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                    for(j=0; j<nbrNoeuds; j++)
                    {
                        if(phenotype[j] == phenotype[nd])
                            phenotype[j] = phenotype[na];
                    }
                }
            }
        }
    }
    ///########## CALCULER LE NOMBRE DE COMPOSANTES CONNEXES #########
    nbrConnectedComponent = 0;
    for(i=0; i<nbrNoeuds; i++)
    {
        clustersSize[i] =0;
        for(j=0; j<nbrNoeuds; j++)
        {
            if (phenotype[j] == i)
            {
                clustersSize[i]++;
            }
        }
        if (clustersSize[i] !=0)
        {
            nbrConnectedComponent++;
        }
    }
    ///########## RETOURNER LE NOMBRE DE COMPOSANTES CONNEXES #########
    return nbrConnectedComponent;
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// delete edges new version
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void newDeleteEdgesLowerThanThreshold(int edge[],int fixedNbrArretes)
{

    int i,j,k,buckUpWeight;
    for(i=0; i<fixedNbrArretes; i++)
    {
        if(edge[i]==1)
        {
            buckUpWeight = fluxMatrix[fixedEdgeVector[i].nouedDepart][fixedEdgeVector[i].nouedArrive];
            fluxMatrix[fixedEdgeVector[i].nouedDepart][fixedEdgeVector[i].nouedArrive]=-1;
            fluxMatrix[fixedEdgeVector[i].nouedArrive][fixedEdgeVector[i].nouedDepart]=-1;

            for(k=0; k<nbrCohabitationConstraintes; k++)
            {
                if((cohabitationConstraintes[k].nouedDepart==fixedEdgeVector[i].nouedDepart  &&
                    cohabitationConstraintes[k].nouedArrive == fixedEdgeVector[i].nouedArrive)||
                   (cohabitationConstraintes[k].nouedDepart==fixedEdgeVector[i].nouedArrive  &&
                    cohabitationConstraintes[k].nouedArrive == fixedEdgeVector[i].nouedDepart))
                {
                    fluxMatrix[fixedEdgeVector[i].nouedDepart][fixedEdgeVector[i].nouedArrive]=-1;
                    fluxMatrix[fixedEdgeVector[i].nouedArrive][fixedEdgeVector[i].nouedDepart]=-1;

                    break;
                }
            }

        }

    }
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Read flow matrix and nbrNoueds from flow matrix file
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void readFlowMatrixAndNbrNoeudsFromFlowMatrixFile(FILE *inputFile){
    int i,j;
    fscanf(inputFile,"%d",&nbrNoeuds);
    for( i=0; i<nbrNoeuds; i++)
    {
        for( j=0; j < nbrNoeuds; j++)
        {
            fscanf(inputFile,"%d",&fluxMatrix[i][j]);
            /// cette partie est utile pour connecter le graphe
            if(i!=j && fluxMatrix[i][j]==0 )fluxMatrix[i][j] = -1;
            buckUpFluxMatrix[i][j] = fluxMatrix[i][j];
        }
    }
}

///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// Read flow matrix and nbrNoueds from edgeList file
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void readFlowMatrixAndNbrNoeudsFromEdgeListFile(FILE *edgeListFile ){

    int indFrom, indTo, weight;

    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///Put the cursor at the beginig  of the file
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii

    while(!feof(edgeListFile)){
        fscanf(edgeListFile,"%d %d %d",&indFrom,&indTo,&weight);
        //printf("[%d][%d] = [%d]\n",indFrom, indTo, weight);
        fluxMatrix[indFrom][indTo] = fluxMatrix[indTo][indFrom] = weight;
    }
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    ///check the matrix before write it down into the file
    ///iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
    int i,j,checkSymetric = 1;

    for(i = 0; i < indFrom; i++) {
        if (fluxMatrix[i][i] != 0){
            fluxMatrix[i][i] = 0;
        }
        for(j = 0; j < indFrom; j++) {
            if(fluxMatrix[i][j] != fluxMatrix[j][i]){
                printf("The matrix is not symetric  \n");
                checkSymetric = 0;
                break;
            }
        }
    }
    if(checkSymetric){
        printf("the matrix is symetric\n");
        exit(-5);
    }

    for(i = 0; i < indFrom; i++) {
        for(j = 0; j < indFrom; j++) {
            /// cette partie est utile pour connecter le graphe
            if(i != j && fluxMatrix[i][j]==0 )
                fluxMatrix[i][j] = fluxMatrix[j][i] = -1;
            buckUpFluxMatrix[i][j] = buckUpFluxMatrix[j][i] = fluxMatrix[i][j];
        }
    }

    nbrNoeuds = indFrom;

    ///###CLOSE THE FILES ###
    fclose(edgeListFile);
}

///FUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTION
/// update the deletedEdges vector
///FUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTION
void updateDeletedEdges(int *aretesSupprimees, int fixedNbrArretes, edge *deletedEdges, int *nbrAretesSupprimees ){
    for(int i = 0; i< fixedNbrArretes ; i++){
        if(aretesSupprimees[i] == 0){
            for(int j=0; j < *nbrAretesSupprimees; j++){
                if(deletedEdges[j].numeroEdge == i ){
                    int k = j;
                    while(k < *nbrAretesSupprimees-1){
                        deletedEdges[k] = deletedEdges[k+1];
                        k++;
                    }
                    *nbrAretesSupprimees--;
                    break;
                }
            }
        }
    }
}

///FUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTION
/// IS FLOW MATRIX SYMETRIC
///FUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTION
int isFlowMatrixSymetric(){

    for(int i=0; i<nbrNoeuds - 1; i++){
        for(int j = i+1; j < nbrNoeuds; j++){
            if(fluxMatrix[i][j] != fluxMatrix[j][i])
                return 0;
        }
    }
    return 1;
}


///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
/// display the edgeVector
///FCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCTFCT
void displayEdgeVector(){
    for(int i=0; i< fixedNbrArretes; i++){
        printf("edge n= %d = (%d,%d)\t",i,(edgeVector+i)->nouedDepart, (edgeVector+i)->nouedArrive);
    }
}


///FUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTION
/// CALCULATE NEIBHORS VECTORS
///FUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTION

void calculateNeighborsVector(){
    for(int i=0; i<nbrNoeuds; i++){
        neighborsVector[i].vertex = i;
        neighborsVector[i].nbrNeighbors = 0;
        neighborsVector[i].sumWeightsIncidentEdges = 0;
        for(int j = 0; j < nbrNoeuds; j++){
            if(i!=j && fluxMatrix[i][j] >= 0){
                neighborsVector[i].listOfNeighbors[neighborsVector[i].nbrNeighbors] = j;
                neighborsVector[i].nbrNeighbors++;
                neighborsVector[i].sumWeightsIncidentEdges +=  fluxMatrix[i][j];
            }
        }
    }
}


///FUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTION
/// CALCULATE NEIBHORHOOD MATRIX
///FUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTIONFUNCTION
void calculateNeighborhoodMatrix() {
    for (int i = 0; i < nbrNoeuds; i++) {
        int nbrNeighbors = 0;
        for (int j = 0; j < nbrNoeuds; ++j) {
            if(i!=j && fluxMatrix[i][j] > 0) {
                nbrNeighbors++;
                neighborMatrix[i][nbrNeighbors] = j;
            }
        }
        neighborMatrix[i][0] = nbrNeighbors;
    }
    //displayNeighborhoodMatrix();
}

void displayNeighborhoodMatrix() {
    for (int i = 0; i < nbrNoeuds; i++) {
        printf("Neighbors of node: %d\t",i);
        for (int j = 1; j <= neighborMatrix[i][0]; j++) {
            printf("%d\t",neighborMatrix[i][j]);
        }
        printf("\n");
    }
}
