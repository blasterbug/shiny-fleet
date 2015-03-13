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

typedef struct {
	int n; /*Nombre de lieux à visiter n*/
	int *d; /*Matrice des distances c_ij*/
	int nbvar; /*Nombre de variables x_ij*/
	int nbcont; /*Nombre de contraintes de base (seulement (1) et (2)*/
} donnees;


int main(int argc, char *argv[])
{
	
}
