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
	if(argc < 2){
		perror("Nombre d'arguments incorrect\n");
		exit(EXIT_FAILURE);
	}
	
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
	double z; //Res fonction objectif
	double *x; //Variables

	/*Nombre de coeff non-nuls dans la matrice de contraintes*/
	int nbcreux;
	
	/*Indices de boucles utilisés ici et là*/
	int i,j=0;
	int pos;

	nbvar = donprob.n*donprob.n;
	nbcont = 2*donprob.n;
	nbcreux = 2*donprob.n*donprob.n;
	printf("nbvar=%d, nbcont=%d, nbcreux=%d \n", nbvar, nbcont, nbcreux);

	prob = glp_create_prob(); /*Allocation mémoire pour le problème*/
	glp_set_prob_name(prob, "trajet"); /* affectation d'un nom */
	glp_set_obj_dir(prob, GLP_MIN); /* Il s'agit d'un problème de minimisation */
	/*
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_OFF; *//* Paramètre de GLPK dans la résolution d'un PL en variables continues afin de couper les affichages (dans lesquels on se noierait) */
	/*glp_iocp parmip;
	glp_init_iocp(&parmip);
	parmip.msg_lev = GLP_MSG_OFF;*/ /* Paramètre de GLPK dans la résolution d'un PL en variables entières afin de couper les affichages (dans lesquels on se noierait) */

	
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
	for (i=1; i<=nbvar; i++) glp_set_obj_coef(prob,i,donprob.d[i-1]);

	/* allocation */
	ia = (int *) malloc ((1 + nbcreux) * sizeof(int));
	ja = (int *) malloc ((1 + nbcreux) * sizeof(int));
	ar = (double *) malloc ((1 + nbcreux) * sizeof(double));

	/*remplissage matrice contrainte*/
	pos = 1;

	for(i=1; i<=donprob.n; i++){
		for (j=1; j<=donprob.n; j++){
			/*Contraintes (1)*/
			ia[pos] = i;
			ja[pos] = j+(i-1)*donprob.n;
			ar[pos] = 1.0;
			pos++;

			/*Contraintes (2)*/
			ia[pos] = j+donprob.n;
			ja[pos] = j+(i-1)*donprob.n;
			ar[pos] = 1.0;
			pos++;
		}
	}
	
	/*Check matrice de contrainte*/
	/*
	for(i=1; i<pos; i++)
	{
		printf("ia[%d] = %d; ja[%d] = %d; ar[%d] = %f;\n", i, ia[i], i, ja[i], i, ar[i]);
	}*/
	
	/* chargement de la matrice dans le problème */
	glp_load_matrix(prob,nbcreux,ia,ja,ar); 
	
	/* écriture de la modélisation dans un fichier*/
	///ERREUR ICI
	//glp_write_lp(prob,NULL,"trajet.lp");

	/* Résolution, puis lecture des résultats */
	glp_simplex(prob,NULL);
	glp_intopt(prob,NULL); /* Résolution */
	
	z = glp_mip_obj_val(prob); /* Récupération de la valeur optimale.*/
	
	/* Récupération des valeurs de variables */
	x = (double *) malloc (nbvar * sizeof(double));
	for(i = 0;i < nbvar; i++) x[i] = glp_mip_col_val(prob,i+1); 

	printf("z = %lf\n",z);
	for(i = 0;i < donprob.n;i++){
		for(j = 0;j < donprob.n;j++){
		printf("x%d_%d = %d, ",i,j,(int)(x[i*donprob.n+j] + 0.5)); /* un cast est ajouté, x[i] pourrait être égal à 0.99999... */
		}
		puts("");
	}
	

	/* libération mémoire */
	glp_delete_prob(prob);
	free(ia);
	free(ja);
	free(ar);
	free(x);

	return 0;
}

