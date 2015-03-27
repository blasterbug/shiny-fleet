/**
 ** Théo DELALANDE-DELARBRE
 ** Benjamin SIENTZOFF
 ** 601A
 ** Shiny fleet
 ** Lecture des données pour le projet
 **/
#ifndef parser
#define parser

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	/* Nombre de lieux à visiter n */
	int n;
	/* Matrice des distances c_ij */
	int *d;
} donnees;

/**
 ** Parcer pour les données
 **/
void shiny_reader( char* filename, donnees *dat )
{
	/* fichier à lire */
	FILE *fichier;
	/* indice */
	int i;
	/* entiers lus et convertis */
	int val, res;
	/* ouverture du fichier en lecture */
	fichier = fopen( filename,"r" );
	
	/* lecture du nombre de lieux à visiter */
	res = fscanf( fichier, "%d", &val );
	//printf( "%d\n", res );
	dat->n = val;
	
	/* allocation mémoire pour les données */
	dat->d = (int*) malloc( dat->n * dat->n * sizeof(int) );
	
	/* remplissage du tableau des distances */
	for( i=0; i < dat->n*dat->n; i++ )
	{
		res = fscanf( fichier, "%d", &val );
		dat->d[i] = val;
		//printf("%d : %d\n", i, val);
	}
	
	/* fermeture du fichier */
	fclose( fichier );
}

#endif
/** fin du fichier **/
