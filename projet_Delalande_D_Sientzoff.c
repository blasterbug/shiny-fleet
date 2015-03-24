/**
 ** Théo DELALANDE-DELARBRE
 ** Benjamin SIENTZOFF
 ** 601A
 ** Shiny fleet
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glpk.h>
#include "shiny_parser.h"

int main( int argc, char *argv[] )
{
	/*Initialisation des données du problème*/
	donnees donprob;

	/*Nombre de variables x_ij*/
	int nbvar;
	/*Nombre de contraintes de base (seulement (1) et (2)*/
	int nbcont;

	/*Initialisation des données et de la matrice creuse pour GLPK */
	glp_prob *prob;
	int *ia;
	int *ja;
	double *ar;

	/*Nombre de coeff non-nuls dans la matrice de contraintes*/
	int nbcreux;

	nbcreux = 2*donprob.n*donprob.n;
}

