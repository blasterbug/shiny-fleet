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
	shiny_reader( argv[1], &donprob);
	
	
	
	/*Nombre de variables x_ij (var binaire : passage de i à j)*/
	int nbvar;
	/*Nombre de contraintes de base (seulement (1) et (2))*/
	int nbcont;

	/*Initialisation des données et de la matrice creuse pour GLPK */
	glp_prob *prob;
	int *ia;
	int *ja;
	double *ar;
	
	/*Resultats*/
	double z;
	double *x;

	/*Nombre de coeff non-nuls dans la matrice de contraintes*/
	int nbcreux;
	
	/*Indices de boucles utilisés ici et là*/
	int i,j=0;

	nbvar = donprob.n*donprob.n;
	nbcont = 2*donprob.n;
	nbcreux = 2*donprob.n*donprob.n;


	prob = glp_create_prob(); /*Allocation mémoire pour le problème*/
	glp_set_prob_name(prob, "trajet"); /* affectation d'un nom */
	glp_set_obj_dir(prob, GLP_MIN); /* Il s'agit d'un problème de minimisation */
	
	/* Déclaration du nombre de contraintes (nombre de lignes de la matrice des contraintes) */
	glp_add_rows(prob, nbcont); 


	/* Bornes sur les contraintes */
	for(i=1;i<=nbcont;i++)
	{
		glp_set_row_bnds(prob, i, GLP_FX, 1.0, 1.0); 
	}
	
	/* Déclaration du nombre de variables */
 	glp_add_cols(prob, nbvar);

	/* On précise le type des variables */
	for(i=1;i<=nbvar;i++)
	{
		/* bornes éventuelles sur les variables, et type */
		glp_set_col_bnds(prob, i, GLP_DB, 0.0, 1.0); 
		/* bornes sur les variables, comme sur les contraintes */
		glp_set_col_kind(prob, i, GLP_BV);
	} 
	
	/* définition des coefficients des variables dans la fonction objectif */
	/*Buggue certainement*/
	for (i=1; i<=nbvar; i++) glp_set_obj_coef(prob,i,donprob.d[i]);

	/* allocation */
	ia = (int *) malloc ((1 + nbcreux) * sizeof(int));
	ja = (int *) malloc ((1 + nbcreux) * sizeof(int));
	ar = (double *) malloc ((1 + nbcreux) * sizeof(double));

	/*remplissage matrice contrainte*/
	int pos = 1;
	/*Contraintes (1)*/
	for(i=1; i<=donprob.n; i++){
		for (j=0; j<donprob.n; j++){
			ia[pos] = i;
			ja[pos] = i+donprob.n;
			ar[pos] = 1.0;
			pos++;
		}
	}
	
	/*Contraintes (2)*/
	for(i=1; i<=donprob.n; i++){
		for (j=0; j<donprob.n; j++){
			ia[pos] = i+donprob.n;
			ja[pos] = i+j*donprob.n;
			ar[pos] = 1.0;
			pos++;
		}
	}
	
	/*Check matrice de contrainte*/
	for(i=1; i<pos; i++)
	{
		printf("ia[%d] = %d; ja[%d] = %d; ar[%d] = %f;\n", i, ia[i], i, ja[i], i, ar[i]);
	}
	
	/* chargement de la matrice dans le problème */
	glp_load_matrix(prob,nbcreux,ia,ja,ar); 
	
	/* écriture de la modélisation dans un fichier*/
	glp_write_lp(prob,NULL,"trajet.lp");

	/* Résolution, puis lecture des résultats */
	glp_simplex(prob,NULL);	glp_intopt(prob,NULL); /* Résolution */
	
	z = glp_mip_obj_val(prob); /* Récupération de la valeur optimale.*/
	
	/* Récupération de la valeur des variables */
	x = (double *) malloc (nbvar * sizeof(double));
	for(i = 0;i < nbvar; i++) x[i] = glp_mip_col_val(prob,i+1); 

	printf("z = %lf\n",z);
	/*for(i = 0;i < nbvar;i++) printf("x%c = %d, ",'B'+i,(int)(x[i] + 0.5)); /* un cast est ajouté, x[i] pourrait être égal à 0.99999... */ 
	puts("");

	/* libération mémoire */
	glp_delete_prob(prob);
	free(ia);
	free(ja);
	free(ar);
	free(x);

	return 0;
}

