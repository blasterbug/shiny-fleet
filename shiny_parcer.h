/**
 ** Théo DELALANDE-DELARBRE
 ** Benjamin SIENTZOFF
 ** 601A
 ** Shiny fleet
 ** Lecture des données pour le projet
 **/
#ifndef parcer
#define parcer

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glpk.h>


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
	/* entiers lus er convertis */
	int val, res;
	
	/* ouverture du fichier en lecture */
	fichier = fopen( filename,"r" );
	
	/* lecture du nombre de lieux à visiter */
	res = fscanf( fichier, "%d", &val );
	dat->n = val;
	
	/* remplissage du tableau des distances */
	for( i=0; i < dat->n*dat->n; i++ )
	{
		res = fscanf( fichier, "%d", &val );
		dat->d[i] = val;
	}
	
	/* fermeture du fichier */
	fclose( fichier );
}
#endif
