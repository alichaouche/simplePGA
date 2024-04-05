
/**================= Sorted Binary Group Number Encoding ====================================================
* 04/12/2015 :
* la population initiale est maintenant pr�te
*============================================================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../headers/typeDeclaration.h"
#include "../headers/communesFunctions.h"
#include "../headers/genetic_alg_VAE.h"
#include "../headers/compilationConditionnelle.h"
#include "../headers/prog_params.h"

///############################################
// Global Variables.
///############################################
extern int sommeTotalFlux,
        nbrNoeuds,
        nbr_constraint,
        nbrArretes,
        fixedNbrArretes,
        nbrNoeudsInCluster,
        fluxMatrix[tailleMax][tailleMax],
        buckUpFluxMatrix[tailleMax][tailleMax],
        *fluxVector,
        max_clusters,
        min_clusters,
        min_sizeCluster,
        max_sizeCluster,
        bestSolution,
        nbrApparition,
        iteration,
        nbrCohabitationConstraintes,
        nbrNonCohabitationConstraintes,
        nbrRun;

extern float MBF, ART, AES, SDBF, BF[40];

extern edge *edgeVector,*fixedEdgeVector,cohabitationConstraintes[100],nonCohabitationConstraintes[100];
///***********************************************
extern partitionVAE *populationVAE1 = NULL,*populationVAE2 = NULL, *solutionDominanteVAE = NULL,*bestSolutionOverRunsVAE = NULL;


///======================================================================================
/// the main programm
void vertexAssignmentEncoding(int nbrGeneration, FILE* outputFileSBGNE, FILE* outputFilePopSBGNE,
     FILE *outputOptimalSolutionFileSBGNE,partitionSBGNE *populationSBGNE1,partitionSBGNE *populationSBGNE2, partitionSBGNE *solutionDominante)
{
    clock_t t1,t2;
    double temps;
    int   iteration=1,bestSolutionIteration,ES=0,i;

    t1=clock();
    genererPopulationInitialeSBGNE(populationSBGNE1,0);
    //genererPopulationInitialeRandomlySBGNE(populationSBGNE1,0);
    getPartitionFromSolutionSBGNE(populationSBGNE1);
    checkContrainstAndFitnessPenalizationSBGNE(populationSBGNE1);
    calculCoutCoupeEtFitnessSBGNE(populationSBGNE1);

    ///initialisation de la solution dominante
    bestSolution=findTheBestSolutionSBGNE(populationSBGNE1);
    *solutionDominante=*(populationSBGNE1+bestSolution);
    bestSolutionIteration = 1;
    nbrApparition=1;
    //iteration++;
    //for(iteration =2 ; iteration <= nbrGeneration; iteration++){
    do{

        naturalSelectionSBGNE(populationSBGNE1,populationSBGNE2);
        crossOverSBGNE(populationSBGNE2, populationSBGNE1);
        mutationSBGNE(populationSBGNE1);
        getPartitionFromSolutionSBGNE(populationSBGNE2);
        checkContrainstAndFitnessPenalizationSBGNE(populationSBGNE2);
        calculCoutCoupeEtFitnessSBGNE(populationSBGNE2);
        bestSolution=findTheBestSolutionSBGNE(populationSBGNE1);

        if((populationSBGNE1+bestSolution)->coutCoupeNormalise > solutionDominante->coutCoupeNormalise )
                /** le 23/12/2015 : pour avoir la meilleur solution avec le moins de contrainte viol�es*/
                /// && (populationSBGNE+bestSolution)->contrainteViole   < solutionDominante->contrainteViole)
        {
            *(solutionDominante) = *(populationSBGNE1+bestSolution);
            nbrApparition=1;
            bestSolutionIteration = iteration;
        }
        else if((populationSBGNE1+bestSolution)->coutCoupeNormalise <= solutionDominante->coutCoupeNormalise )
        {
            nbrApparition++;
        }

        ///iteration++;
    }while(nbrApparition <max_steady_generation);
    t2=clock();
    displayTheBestSolutionSBGNE(solutionDominante);

    temps = (double)(t2-t1)/CLOCKS_PER_SEC;
    if(bestSolutionIteration>=2){
        ES = ((taillePopulation-elitismeRate)*(bestSolutionIteration-2))+(taillePopulation+ solutionDominante->id -elitismeRate+1);
    }
    else {
        ES = solutionDominante->id +1;
    }
    writeOptimalSolutionInFileSBGNE(solutionDominante,outputOptimalSolutionFileSBGNE,nbrRun,bestSolutionIteration,temps,ES);
    printf("le temps d execution est %lf \n", temps);

    fprintf(outputFileSBGNE,"\n");
    writeBestSolutionInFileSBGNE(solutionDominante,outputFileSBGNE,bestSolutionIteration);
}
///**************************************************************************************
void genererPopulationInitialeSBGNE(partitionSBGNE *populationSBGNE, int indiceFirstElt)
{

    int i,nbrSolutionRepeter;

    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        (populationSBGNE+i)->id =i;
        generateBinarySolutionSBGNE(populationSBGNE, i);
    }
}

///**************************************************************************************
void genererPopulationInitialeRandomlySBGNE(partitionSBGNE *populationSBGNE, int indiceFirstElt)
{

    int i;
    ///printf("the genererPopulationInitialeSBGNE function ...\n");
    for(i=indiceFirstElt; i<taillePopulation; i++)
    {
        (populationSBGNE+i)->id =i;
        generateBinarySolutionSBGNE(populationSBGNE, i);
    }
}
/// la g�n�ration des solution binaires

void generateBinarySolutionSBGNE_Before_11012019(partitionSBGNE* populationSBGNE, int indiceIndividu)
{

    int numeroPartie,i,indiceNoeudAffecte,chromosomeSize = nbrParties * nbrNoeuds;
    for(i=0; i<chromosomeSize; i++)
    {
        //(populationSBGNE+indiceIndividu)->genotype[i] = 0;
        (populationSBGNE+indiceIndividu)->genotype[i] = rnd(0,2);
    }
/*
    for(i=0; i<nbrNoeuds; i++)
    {
        ///numeroPartie=rnd(0,max_clusters);
        numeroPartie=rnd(0,nbrParties);
        indiceNoeudAffecte = nbrNoeuds*numeroPartie + i;
        (populationSBGNE+indiceIndividu)->genotype[indiceNoeudAffecte] = 1;
    }
*/
    homoginiserLesSolutionBinaire(populationSBGNE,indiceIndividu);
}
///********************************************************************************************
/**
    apr�s la g�n�ration de la solution binaire, il va falloir proc�der � son homogination par le fait
    que deux noeuds n'ont pas le droit d'appartenir � deux cluster au m�me temps, voici la proc�dure :
    - si un noeud n'est gu�re affect� ==> il sera affect� au dernier cluster
    - si un noued est affect� � plusieurs cluster ==> il sera affecter au cluster du plus petit rang

    Important : Cette technique sera utilis�e lors de croiement
*/
void generateBinarySolutionSBGNE(partitionSBGNE* populationSBGNE, int indiceIndividu)
{

    int numeroPartie,i,j,indiceNoeudAffecte,chromosomeSize = nbrParties*nbrNoeuds,nd,na;
    for(i=0; i<chromosomeSize; i++){(populationSBGNE+indiceIndividu)->genotype[i] = 0;}
	for(i=0;i<nbrParties;i++)(populationSBGNE+indiceIndividu)->clustersSize[i]=0;
    for(i=0; i<nbrNoeuds; i++)
    {
        do{
			numeroPartie=rnd(0,nbrParties);
		}while((populationSBGNE+indiceIndividu)->clustersSize[numeroPartie]>=max_sizeCluster);
        indiceNoeudAffecte = nbrNoeuds*numeroPartie + i;

        (populationSBGNE+indiceIndividu)->genotype[indiceNoeudAffecte] = 1;
		(populationSBGNE+indiceIndividu)->clustersSize[numeroPartie]++;
    }

    for(i=0;i<nbrCohabitationConstraintes;i++){
            for(j=0; j<nbrParties; j++)
            {
                nd = cohabitationConstraintes[i].nouedDepart + (j*nbrNoeuds);
                na = cohabitationConstraintes[i].nouedArrive + (j*nbrNoeuds);
               if((populationSBGNE+indiceIndividu)->genotype[nd] == 1 || (populationSBGNE+indiceIndividu)->genotype[na] == 1){
                    (populationSBGNE+indiceIndividu)->genotype[na] = 1;
                    (populationSBGNE+indiceIndividu)->genotype[nd] = 1;
                }

            }
        }
    homoginiserLesSolutionBinaire(populationSBGNE,indiceIndividu);

}
void calculateClusterSize(char genotype[2000], int clustersSize[100]){
	int i,j,k;
	for(i=0;i<nbrParties;i++){
		clustersSize[i]=0;
		for(j=0;j<nbrNoeuds;j++){
			if(genotype[i*nbrNoeuds+j]==1)clustersSize[i]++;
		}
	}
}
/***************************************************************************************/
void homoginiserLesSolutionBinaire(partitionSBGNE* populationSBGNE, int indiceIndividu)
{

    int i,j,k,choixCluster,noeudAffecte ,numeroPartie,indiceNoeudAffecte,
    chromosomeSize=nbrParties*nbrNoeuds,nbrApparition,l;
    calculateClusterSize((populationSBGNE+indiceIndividu)->genotype, (populationSBGNE+indiceIndividu)->clustersSize);
    for(k=0; k<nbrNoeuds; k++){
        noeudAffecte = 0;
        ///for(j=0; j<max_clusters; j++)
        for(j=0; j<nbrParties; j++){
            if((populationSBGNE+indiceIndividu)->genotype[(j*nbrNoeuds) + k]  == 1 && (populationSBGNE+indiceIndividu)->clustersSize[j]<=max_sizeCluster){
                for(i=j+1; i<nbrParties; i++){
                    if((populationSBGNE+indiceIndividu)->genotype[(i*nbrNoeuds) + k] ==1){
                        (populationSBGNE+indiceIndividu)->genotype[(i*nbrNoeuds) + k] =0;
                        (populationSBGNE+indiceIndividu)->clustersSize[i]--;
                    }
                }
                noeudAffecte = 1;
                break;
            }
            else if((populationSBGNE+indiceIndividu)->genotype[(j*nbrNoeuds) + k]  == 1 && (populationSBGNE+indiceIndividu)->clustersSize[j]> max_sizeCluster){
                nbrApparition = 0;
                for(l=k;l<chromosomeSize;l+=nbrNoeuds)if((populationSBGNE+indiceIndividu)->genotype[l] == 1) nbrApparition++;
                if(nbrApparition > 1){
                        (populationSBGNE+indiceIndividu)->genotype[(j*nbrNoeuds) + k] = 0;
                        (populationSBGNE+indiceIndividu)->clustersSize[j]--;
                }
            }
        }
        if(!noeudAffecte)
        {
            //do{
                numeroPartie=rnd(0,nbrParties);
            //}while((populationSBGNE+indiceIndividu)->clustersSize[numeroPartie]>=max_sizeCluster);
            indiceNoeudAffecte = nbrNoeuds*numeroPartie + i;

            (populationSBGNE+indiceIndividu)->genotype[indiceNoeudAffecte] = 1;
            (populationSBGNE+indiceIndividu)->clustersSize[numeroPartie]++;
        }
    }
}


/// conversion des solution binaires g�n�r�es au solutions entiere pour illimin� la rendondance
///********************************************************************************************
void conversionSolutionBinaireToEntiereSBGNE(partitionSBGNE* populationSBGNE,int indiceIndividu)
{

    int i,j;
    ///for(i=0; i<max_clusters; i++)
    for(i=0; i<nbrParties; i++)
    {
        (populationSBGNE+indiceIndividu)->genotypeEntier[i] = 0;
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationSBGNE+indiceIndividu)->genotypeEntier[i] =
                (populationSBGNE+indiceIndividu)->genotypeEntier[i] +
                (pow(2,(nbrNoeuds-1-j))*(populationSBGNE+indiceIndividu)->genotype[(i*nbrNoeuds)+j]);
        }
    }
}

///**************************************************************************************
void trieGenotypeEntieSBGNE(partitionSBGNE* populationSBGNE, int indiceIndividu)
{

    int i,j,tmp;
    ///for(i=0; i<max_clusters-1; i++)
    for(i=0; i<nbrParties-1; i++)
    {
        ///for(j=i+1; j<max_clusters; j++)
        for(j=i+1; j<nbrParties; j++)
        {
            if((populationSBGNE+indiceIndividu)->genotypeEntier[i] > (populationSBGNE+indiceIndividu)->genotypeEntier[j])
            {
                tmp = (populationSBGNE+indiceIndividu)->genotypeEntier[i];
                (populationSBGNE+indiceIndividu)->genotypeEntier[i] = (populationSBGNE+indiceIndividu)->genotypeEntier[j];
                (populationSBGNE+indiceIndividu)->genotypeEntier[j] = tmp;
            }
        }
    }
}

///**************************************************************************************
int existanceDeSolutionSBGNE(partitionSBGNE* populationSBGNE, int indexNewSolution)
{

    int i=0,j,Exist = 0,theSameElement;

    while(i<indexNewSolution && Exist == 0)
    {
        theSameElement = 0;
        ///for(j=0; j<max_clusters; j++)
        for(j=0; j<nbrParties; j++)
        {
            if((populationSBGNE+i)->genotypeEntier[j] == (populationSBGNE+indexNewSolution)->genotypeEntier[j])
            {
                theSameElement++;
            }
        }
        ///if(theSameElement==max_clusters)
        if(theSameElement==nbrParties)
        {
            Exist = 1;
            break;
        }
        i++;
    }
    return Exist;
}


///**************************************************************************************
void getPartitionFromSolutionSBGNE(partitionSBGNE *populationSBGNE)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {

        /// d�termination de la parition
        ///for(j=0; j<max_clusters; j++)
        for(j=0; j<nbrParties; j++)
        {
            for(k=0; k<nbrNoeuds; k++)
            {
                if((populationSBGNE+i)->genotype[(j*nbrNoeuds)+k] == 1)
                {
                    (populationSBGNE+i)->phenotype[k] = j;
                }
            }
        }
    }
}


///***************************************************************************************
/// calculer le cout de coupe de chaque partition ainsi que son fitness
///***************************************************************************************
/// on va appliquer la formule suivante pour le calcule des cout de coupe
/// T'(C) = B - T(C)            <=> la minimisation devient une maximisation
/// Z(C)=T'(C)+(u-v(C)) �B      <=>  calculer la nouvelle fitness par rapport aux contrainte
/// T'(C) : cout de coupe (maximiser les flux INTRA), u ; nombre de contrainte totale
/// v(C): nombre de contrainte viol�,  B: la somme totale des flux
///***************************************************************************************

void calculCoutCoupeEtFitnessSBGNE(partitionSBGNE* populationSBGNE)
{
    int i,j,k,sommeTotalCoutCoupe = 0;

    for(i=0; i<taillePopulation; i++) /// parcourir tous les individus de la population
    {
        (populationSBGNE+i)->coutCoupe = 0; /// r�initialisation de la valriable apr�s chaque it�ration (individus de la population)
        for(j=0; j<nbrNoeuds-1; j++)
        {
            for(k=j+1; k<nbrNoeuds; k++)
            {
                if((populationSBGNE+i)->phenotype[j] == (populationSBGNE+i)->phenotype[k] && fluxMatrix[j][k] > 0)
                {
                    (populationSBGNE+i)->coutCoupe = (populationSBGNE+i)->coutCoupe +  fluxMatrix[j][k];
                }
            }
        }
        (populationSBGNE+i)->coutCoupeNormalise = (populationSBGNE+i)->coutCoupe + ((nbr_constraint - (populationSBGNE+i)->contrainteViole)*sommeTotalFlux);
        sommeTotalCoutCoupe = sommeTotalCoutCoupe + (populationSBGNE+i)->coutCoupeNormalise;
    }

#if scaling

    ///======================================================================================================
    /// Scaling Fintnss : sigma scaling Melanie mitchelle r�f�rence Goldberg | le 26/11/2015
    ///pour r�gler le probl�me des valeurs n�gatives des expected values ==> utiliser sigma truncation de goldberg
    /// g(f) = f + (moyenne - c * sigma )
    /// 1 <= c <= 5, g�n�ralement initilis� � 2
    ///======================================================================================================
    float varianceCoutDeCoupeNormalise = 0, moyenneCoutDeCoupeNormalise, sigma, sommeTotalOfExpectdValues = 0;
    int c= 2;

    moyenneCoutDeCoupeNormalise = sommeTotalCoutCoupe / taillePopulation;
    for(i=0; i<taillePopulation; i++)
    {
        varianceCoutDeCoupeNormalise =
            varianceCoutDeCoupeNormalise + pow(((populationSBGNE+i)->coutCoupeNormalise - moyenneCoutDeCoupeNormalise),2);
    }
    varianceCoutDeCoupeNormalise = varianceCoutDeCoupeNormalise / taillePopulation;
    sigma = sqrt(varianceCoutDeCoupeNormalise);
    /// tester la valeur ajout�
    for(i=0; i<taillePopulation; i++)
    {

        /// calcul de sigma truncation
        (populationSBGNE+i)->expectedValue = (populationSBGNE+i)->coutCoupeNormalise + (moyenneCoutDeCoupeNormalise - c*sigma);
        sommeTotalOfExpectdValues = sommeTotalOfExpectdValues + (populationSBGNE+i)->expectedValue;
    }

    for(i=0; i<taillePopulation ; i++)
    {
        (populationSBGNE+i)->fitness = (float)((populationSBGNE+i)->expectedValue)/(float)(sommeTotalOfExpectdValues);
    }
#else
    for(i=0; i<taillePopulation ; i++)
    {
        (populationSBGNE+i)->fitness = (float)((populationSBGNE+i)->coutCoupeNormalise)/(float)(sommeTotalCoutCoupe);
    }
#endif// scaling



}
///**************************************************************************************
void checkContrainstAndFitnessPenalizationSBGNE(partitionSBGNE *populationSBGNE)
{

    int i,j,k;

    for(i=0; i<taillePopulation; i++)
    {
        (populationSBGNE+i)->nbrCluster=0; /// le nombre de cluster doit etre initialis� � 1 pour prendre en cons�diration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes viol�es
        (populationSBGNE+i)->contrainteViole=0;
        for(j=0; j<4; j++)
        {
            (populationSBGNE+i)->constraintVector[j]=0;
        }
        ///**********************************************************************************


        ///trier le tableau r�sulant
        for(j=0; j<nbrNoeuds; j++)
        {
            (populationSBGNE+i)->clustersSize[j] =0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if ((populationSBGNE+i)->phenotype[k] == j)
                {
                    (populationSBGNE+i)->clustersSize[j]++;
                }
            }
            if ((populationSBGNE+i)->clustersSize[j] !=0)
            {
                (populationSBGNE+i)->nbrCluster++;
            }
        }

        if((populationSBGNE+i)->nbrCluster > max_clusters || (populationSBGNE+i)->nbrCluster < min_clusters)
        {
            (populationSBGNE+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationSBGNE+i)->contrainteViole++;
        }

        for(j=0; j<nbrNoeuds; j++)
        {
            if((populationSBGNE+i)->clustersSize[j]!=0)
            {
                if((populationSBGNE+i)->clustersSize[j]>max_sizeCluster || (populationSBGNE+i)->clustersSize[j]< min_sizeCluster)
                {
                    (populationSBGNE+i)->constraintVector[1]++;

                }
            }

        }

        if((populationSBGNE+i)->constraintVector[1] != 0)
        {
            (populationSBGNE+i)->contrainteViole++;
        }
    }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0)
    {
        checkCohabitationAndNonCohabitationConstrainteSBGNE(populationSBGNE);
    }

}

///**************************************************************************************
void checkCohabitationAndNonCohabitationConstrainteSBGNE(partitionSBGNE *populationSBGNE)
{

    int i,j,noeud1,noeud2;
    for(i=0; i<taillePopulation; i++)
    {
	    /// Cette tache est d�j� r�alis� au niveau de la fonction checkCohabitationAndNonCohabitationConstrainteBE
		/**
        (populationSBGNE+i)->constraintVector[2]=0;
        (populationSBGNE+i)->constraintVector[3]=0;
		*/
        if(nbrCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrCohabitationConstraintes; j++)
            {
                noeud1=cohabitationConstraintes[j].nouedDepart;
                noeud2=cohabitationConstraintes[j].nouedArrive;

                if((populationSBGNE+i)->phenotype[noeud1]!= (populationSBGNE+i)->phenotype[noeud2])
                {
                    (populationSBGNE+i)->constraintVector[2]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes viol�es
            if((populationSBGNE+i)->constraintVector[2]!=0)
            {
                (populationSBGNE+i)->contrainteViole++;
            }
        }


        if(nbrNonCohabitationConstraintes!=0)
        {
            for(j=0; j<nbrNonCohabitationConstraintes; j++)
            {
                noeud1=nonCohabitationConstraintes[j].nouedDepart;
                noeud2=nonCohabitationConstraintes[j].nouedArrive;

                if((populationSBGNE+i)->phenotype[noeud1]== (populationSBGNE+i)->phenotype[noeud2])
                {
                    (populationSBGNE+i)->constraintVector[3]++;
                }
            }
            /// cette instruction permet d'avoir le nombre de contraintes viol�es
            if((populationSBGNE+i)->constraintVector[3]!=0)
            {
                (populationSBGNE+i)->contrainteViole++;
            }
        }


    }

}

///**************************************************************************************
/// la s�l�ction naturelle des individus pour la nouvelle g�n�ration
void naturalSelectionSBGNE(partitionSBGNE* populationSBGNE1,partitionSBGNE* populationSBGNE2)
{
   int i,j,maxFitness;
    int tmpVectIndice[taillePopulation];
    for(i=0;i<taillePopulation;i++) tmpVectIndice[i]=-1;
    float lotterie,sommeFitness;

    for(i=0; i<elitismeRate; i++){
        maxFitness = 0;
        for(j=0;j<taillePopulation;j++){
            if((populationSBGNE1+maxFitness)->coutCoupeNormalise < (populationSBGNE1+j)->coutCoupeNormalise
               && tmpVectIndice[j] == -1 ){
                maxFitness = j;

            }
        }
        tmpVectIndice[maxFitness] = 0;
        *(populationSBGNE2+i) = *(populationSBGNE1+maxFitness);
    }

    for(i=elitismeRate; i<taillePopulation; i++)
    {
#if Windows
        lotterie= drand48ForWindows(0,1001);
#else
        lotterie= drand48();
#endif // Windows
        sommeFitness = 0;

        for(j=0; j<taillePopulation; j++)
        {
            sommeFitness = sommeFitness + (populationSBGNE1+j)->fitness;
            if (lotterie<=sommeFitness)
            {
                *(populationSBGNE2+i) = *(populationSBGNE1+j);
                break;
            }
            else if( j==taillePopulation-1)
            {
                *(populationSBGNE2+i) = *(populationSBGNE1+j);
            }
        }
    }
}
///**************************************************************************************
///le croisement des solution s�l�ctionner
void crossOverSBGNE(partitionSBGNE* populationSBGNE1, partitionSBGNE* populationSBGNE2)
{

/// le croisement se fait au niveau des points d'articulation du g�notype

    int i = 0,j,k;
    int choixInd1, choixInd2, choixLocus=0, chromosomeSize = nbrNoeuds*nbrParties;
    ///***********************************************************************************
    for(i=0; i<elitismeRate; i++)
    {
        *(populationSBGNE2+i) = *(populationSBGNE1+i);
        (populationSBGNE2+i)->id = i;
    }
    ///***********************************************************************************
    while (i < taillePopulation-regeneration)
    {
        choixInd1 = rnd(0,taillePopulation);
        do
        {
            choixInd2 = rnd(0,taillePopulation);
        }
        while(choixInd2 == choixInd1 );
        choixLocus = rnd(0,chromosomeSize-1); /// -1 pour s'arr�ter � moins 2
        //choixLocus = rnd(0,nbrParties-1); /// -1 pour s'arr�ter � moins 2

        for(j=0; j<=choixLocus; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        {
            (populationSBGNE2+i)->genotype[j] = (populationSBGNE1+choixInd1)->genotype[j];
            (populationSBGNE2+i+1)->genotype[j] = (populationSBGNE1+choixInd2)->genotype[j];
        }

        ///************************************************************************************

        //for(j=choixLocus+1; j<nbrParties; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        for(j=choixLocus+1; j<chromosomeSize; j++)  /// l'appel est fait avec nbrNoeuds+1 pour inclure la derni�re valeure
        {
            (populationSBGNE2+i)->genotype[j] = (populationSBGNE1+choixInd2)->genotype[j];
            (populationSBGNE2+i+1)->genotype[j] = (populationSBGNE1+choixInd1)->genotype[j];
        }

        /// affectation de nouveau indice
        (populationSBGNE2+i)->id = i;
        (populationSBGNE2+i+1)->id = i+1;

        homoginiserLesSolutionBinaire(populationSBGNE2, i);
        homoginiserLesSolutionBinaire(populationSBGNE2, i+1);
        ///*****************************************************************
        i+=2;
    }
    generateBinarySolutionSBGNE(populationSBGNE2, taillePopulation-regeneration);
}
///**************************************************************************************
int findTheBestSolutionSBGNE(partitionSBGNE *populationSBGNE)
{
    float maxFitness = (populationSBGNE+0)->fitness;
    int i,indice=0;
    for(i=0; i<taillePopulation ; i++)
    {
        if(maxFitness < (populationSBGNE+i)->fitness)
        {
            maxFitness = (populationSBGNE+i)->fitness;
            indice = i;
        }
    }
    return indice;

}
///**************************************************************************************
void mutationSBGNE(partitionSBGNE* populationSBGNE)
{

    int i,j,numeroGeneMute, choixCluster, choixNoeud,chromosomeSize = nbrNoeuds*nbrParties;
    double applicationOfMutation;

    for(i=elitismeRate; i<taillePopulation; i++)
    {
#if Windows
        applicationOfMutation = drand48ForWindows(0,1001);
#else
        applicationOfMutation = drand48();
#endif // Windows
        if(applicationOfMutation <= tauxMutation)
        {
            /// choisir un g�ne al�atoirement pour lui appliquer la mutation
            ///choixCluster = rnd(0,max_clusters);
            choixCluster = rnd(0,nbrParties);
            choixNoeud = rnd(0,nbrNoeuds);
            numeroGeneMute = (choixCluster*nbrNoeuds)+choixNoeud;
            if((populationSBGNE+i)->genotype[numeroGeneMute] == 1)
            {
                (populationSBGNE+i)->genotype[numeroGeneMute] =0;
            }
            else
            {
                (populationSBGNE+i)->genotype[numeroGeneMute] =1;
            }
            ///apr�s l'application de la mutation, il va falloire remettre le genotype en ordre
            ///for(j=0; j<max_clusters; j++)
            for(j=0; j<nbrParties; j++)
            {
                if(j*nbrNoeuds+choixNoeud != numeroGeneMute && (populationSBGNE+i)->genotype[numeroGeneMute] == 1)
                {
                    (populationSBGNE+i)->genotype[j*nbrNoeuds+choixNoeud]=0;
                }
                else if(j*nbrNoeuds+choixNoeud != numeroGeneMute && (populationSBGNE+i)->genotype[numeroGeneMute] == 0)
                {
                    (populationSBGNE+i)->genotype[j*nbrNoeuds+choixNoeud]=1;
                }

            }

        }
    }
}
///**************************************************************************************
float testerLaSommeDesFitnessSBGNE(partitionSBGNE* populationSBGNE)
{
    int i;
    float sommeFitness = 0;
    for(i=0; i<taillePopulation ; i++)
    {
        sommeFitness = sommeFitness + (populationSBGNE+i)->fitness;
    }
    return sommeFitness;
}
///**************************************************************************************
void displayTheBestSolutionSBGNE(partitionSBGNE* solutionDominante)
{
    int i;
    printf("id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
           solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    printf("le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte sur le nombre de clusters : OK \n");
            else printf("contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte sur la taille des clusters : OK \n");
            else printf("contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte de cohabitation : OK \n");
            else printf("contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) printf("contrainte Non Cohabitation : OK \n");
            else printf("contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            printf("Cette contrainte n'est pas prise en charge par le syst�me \n");

        }

    }
}

///**************************************************************************************
void writeSolutionInFileSBGNE(partitionSBGNE *populationSBGNE, FILE *outputFilePop,int iteration)
{

    int i,j,k;
    for(i=0; i<taillePopulation; i++)
    {

        ///�criture de la meilleur solution dans le fichier output
        fprintf(outputFilePop,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,(populationSBGNE+i)->id,
                (populationSBGNE+i)->coutCoupe,(populationSBGNE+i)->fitness,(populationSBGNE+i)->coutCoupeNormalise,
                (populationSBGNE+i)->contrainteViole,(populationSBGNE+i)->nbrCluster);
        ///for(j=0; j<max_clusters; j++)
        for(j=0; j<nbrParties; j++)
        {
            for(k=0; k<nbrNoeuds; k++)
            {
                fprintf(outputFilePop,"%d ",(populationSBGNE+i)->genotype[j*nbrNoeuds+k]);
            }
        }
        fprintf(outputFilePop,"\t\t\t\t");
        for(j=0; j<nbrNoeuds; j++)
        {
            fprintf(outputFilePop,"%d ",(populationSBGNE+i)->phenotype[j]);
        }
        fprintf(outputFilePop,"\n");
    }
    fprintf(outputFilePop,"\n\n");
}
///**************************************************************************************
void writeBestSolutionInFileSBGNE(partitionSBGNE *solutionDominante, FILE *outputFile,int iteration)
{

    int i;
    fprintf(outputFile,"%d\t%d\t%d\t%0.2f\t%d\t%d\t%d\t",iteration,solutionDominante->id,
            solutionDominante->coutCoupe,solutionDominante->fitness,solutionDominante->coutCoupeNormalise,
            solutionDominante->contrainteViole,solutionDominante->nbrCluster);
    ///for(i=0; i<max_clusters*nbrNoeuds; i++)
    for(i=0; i<nbrParties*nbrNoeuds; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->genotype[i]);
    }
    fprintf(outputFile,"\t\t\t\t");
    for(i=0; i<nbrNoeuds; i++)
    {
        fprintf(outputFile,"%d ",solutionDominante->phenotype[i]);
    }
    fprintf(outputFile,"\n");

}

///**************************************************************************************
void affichePopulationSBGNE(partitionSBGNE* populationSBGNE)
{
    printf("affichePopulation...\n");
    int i,j,k,l;
    for(i=0; i<taillePopulation; i++)
    {

        printf("\nid = %d | cout de coupe = %d | fitness = %0.4f",
               (populationSBGNE+i)->id,(populationSBGNE+i)->coutCoupe, (populationSBGNE+i)->fitness);

        /// affichage de vecteur de la partition
        /**       printf("\nla partition est  = \n");
               for(l=0;l<nbrPartie;l++){
                   printf("%d ",(populationSBGNE+i)->genotypeEntier[l]);
               }
        */

        printf("\nla solution generer est  = \n");
        ///for(j=0; j<max_clusters; j++)
        for(j=0; j<nbrParties; j++)
        {
            for(k=0; k<nbrNoeuds; k++)
            {
                ///printf("%d ",(populationSBGNE+i)->genotype[k+(j*max_clusters)]);
                printf("%d ",(populationSBGNE+i)->genotype[k+(j*nbrParties)]);
            }
            printf("\n");
        }
        printf("\n");
        /// affichage de vecteur de la partition
        printf("\nla partition est  = \n");
        for(l=0; l<nbrNoeuds; l++)
        {
            printf("%d ",(populationSBGNE+i)->phenotype[l]);
        }
        printf("\n===============================================================================\n");

    }
}

///********************************************************************************************
void writeOptimalSolutionInFileSBGNE(partitionSBGNE *solutionDominante,FILE* outputOptimalSolutionFileSBGNE,
                                     int nbrRun, int bestSolutionIteration, float runTime, int ES)
{
    int i;
    fprintf(outputOptimalSolutionFileSBGNE,"RUN NUMBER = %d | ITERATION N = %d | SOLUTION IDENTIFIER = %d | CUT SIZE INTRA = %d | CUT SIZE INTER = %d | SOMME TOTAL DES FLUX = %d | RUN TIME = %f | ES = %d | NBR VIOLATED CONSTRAINTE = %d | ",
            nbrRun, bestSolutionIteration, solutionDominante->id, solutionDominante->coutCoupe,(sommeTotalFlux - solutionDominante->coutCoupe),sommeTotalFlux ,runTime, ES, solutionDominante->contrainteViole);
 fprintf(outputOptimalSolutionFileSBGNE," max_clusters = %d | max_sizeCluster = %d  |  cohabitation = %d | non cohabitation = %d\n",
            solutionDominante->constraintVector[0],solutionDominante->constraintVector[1],solutionDominante->constraintVector[2],
            solutionDominante->constraintVector[3]);


}
///******************************************************************************************
int compareCroissantFitnessSBGNE (void const *a, void const *b)
{

    partitionSBGNE const *pa = a;
    partitionSBGNE const *pb = b;
    return  pb->coutCoupeNormalise - pa->coutCoupeNormalise;
}
/**************************************************************************************

///**************************************************************************************
void writeOptimalSolutionInFileSBGNE(partitionSBGNE *solutionDominante,FILE* outputOptimalSolutionFileSBGNE)
{
    int i;
    fprintf(outputOptimalSolutionFileSBGNE,"id = %d | normalised cut size = %d | cut size intrA = %d cut size intEr = %d ",
            solutionDominante->id,solutionDominante->coutCoupeNormalise, solutionDominante->coutCoupe, (sommeTotalFlux-solutionDominante->coutCoupe));
    fprintf(outputOptimalSolutionFileSBGNE,"le contraintes violees = %d \n", solutionDominante->contrainteViole);
    for(i=0; i<4; i++)
    {
        switch (i)
        {

        case 0:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileSBGNE,"contrainte sur le nombre de clusters : OK \n");
            else fprintf(outputOptimalSolutionFileSBGNE,"contrainte sur le nombre de clusters : NOT OK \n");
            break;
        case 1:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileSBGNE,"contrainte sur la taille des clusters : OK \n");
            else fprintf(outputOptimalSolutionFileSBGNE,"contrainte sur la taille des clusters : NOT OK \n");
            break;
        case 2:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileSBGNE,"contrainte de cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileSBGNE,"contrainte de cohabitation : NOT OK \n");
            break;
        case 3:
            if(solutionDominante->constraintVector[i] ==0) fprintf(outputOptimalSolutionFileSBGNE,"contrainte Non Cohabitation : OK \n");
            else fprintf(outputOptimalSolutionFileSBGNE,"contrainte Non Cohabitation : NOT OK \n");
            break;
        default:
            fprintf(outputOptimalSolutionFileSBGNE,"Cette contrainte n'est pas prise en charge par le syst�me \n");

        }

    }

}

void checkContrainstAndFitnessPenalizationSBGNE(partitionSBGNE *populationSBGNE)
{

    int i,j,k,m=nbrNoeuds,tabTmp[m],tmp;

    for(i=0; i<taillePopulation; i++)
    {
        (populationSBGNE+i)->nbrCluster=1; /// le nombre de cluster doit etre initialis� � 1 pour prendre en cons�diration le premier cluster

        ///initialisation de vecteur des contraintes et de la variable de contraintes viol�es
        (populationSBGNE+i)->contrainteViole=0;
        for(j=0;j<nbr_constraint;j++){(populationSBGNE+i)->constraintVector[j]=0;}


        m=nbrNoeuds;
        ///copier le phenotype dans un deuxi�me tableau
        for(j=0; j<nbrNoeuds; j++)
        {
            tabTmp[j]=(populationSBGNE+i)->phenotype[j];
        }
        ///trier le tableau r�sulant
        for(j=0; j<m-1; j++)
        {
            for(k=j+1; k<m; k++)
            {
                if (tabTmp[k] < tabTmp[j])
                {
                    tmp = tabTmp[j];
                    tabTmp[j] = tabTmp[k];
                    tabTmp[k] = tmp;
                }
            }
        }
        ///calculer le nombre des cluser
        for(j=1; j<nbrNoeuds; j++)
        {
            if(tabTmp[j-1] != tabTmp[j])
            {
                (populationSBGNE+i)->nbrCluster++;
            }
        }
        if((populationSBGNE+i)->nbrCluster > nbrPartie || (populationSBGNE+i)->nbrCluster < min_clusters)
        {
            (populationSBGNE+i)->constraintVector[0] = 1; /// le nombre de cluster n'est pas valide
            (populationSBGNE+i)->contrainteViole++;
        }

        /// v�rifier la taille des clusters

        for(j=0; j<m-1; j++)  /// je dois me limitter � m-1 pour que t arrive jusqu'� m
        {
            while(tabTmp[j] == tabTmp[j+1])
            {
                for(k=j; k<m; k++) /// la taille du vecteur change � chaque fois => m et non pas nbrNoeuds
                {
                    tabTmp[k] = tabTmp[k+1];
                }
                m--;
            }
        }
        for(j=0; j<m; j++)
        {
            (populationSBGNE+i)->clustersSize[j]=0;
            for(k=0; k<nbrNoeuds; k++)
            {
                if (tabTmp[j] == (populationSBGNE+i)->phenotype[k])
                {
                    (populationSBGNE+i)->clustersSize[j]=(populationSBGNE+i)->clustersSize[j]+1;
                }
            }
            if((populationSBGNE+i)->clustersSize[j]>max_sizeCluster || (populationSBGNE+i)->clustersSize[j]< min_sizeCluster)
            {
                ///24.07.2015 ---> l'ajout de vecteur de contraintes = le nombre de cluster dont la taille a d�pass� la borne
                (populationSBGNE+i)->constraintVector[1]++;

            }
        }

        if((populationSBGNE+i)->constraintVector[1] != 0)
        {
            (populationSBGNE+i)->contrainteViole++;
        }
   }
    if(nbrCohabitationConstraintes!=0 || nbrNonCohabitationConstraintes !=0){
        checkCohabitationAndNonCohabitationConstrainteSBGNE(populationSBGNE);
    }

}

*/
/**
void genererPopulationInitialeSBGNE(partitionSBGNE *populationSBGNE)
{

    int i,nbrSolutionRepeter;
    ///printf("the genererPopulationInitialeSBGNE function ...\n");
    for(i=0; i<taillePopulation; i++)
    {
        (populationSBGNE+i)->id =i;
        ///printf("(populationSBGNE+i)->id = %d",(populationSBGNE+i)->id);
        /// affecter la valeur 0 � tous les noeuds de la solution

        nbrSolutionRepeter =0;

        do
        {
            generateBinarySolutionSBGNE(populationSBGNE, i);
            conversionSolutionBinaireToEntiereSBGNE(populationSBGNE,i);
            trieGenotypeEntieSBGNE(populationSBGNE, i);
            nbrSolutionRepeter++;
        }
        while(existanceDeSolutionSBGNE(populationSBGNE, i));
        ///printf("cr�ation du %d individus apr�s %d tentatives...!\n",i,nbrSolutionRepeter);

    }
}*/
