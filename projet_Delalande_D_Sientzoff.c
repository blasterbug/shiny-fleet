/**
 ** Théo DELALANDE-DELARBRE
 ** Benjamin SIENTZOFF
 ** 601A
 ** Shiny fleet
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <glpk.h>
#include "shiny_parser.h"
#include "sous_cycles.h"

/* Déclarations pour le compteur de temps CPU */
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>


struct timeval start_utime, stop_utime;

void crono_start()
{
	struct rusage rusage;
	
	getrusage(RUSAGE_SELF, &rusage);
	start_utime = rusage.ru_utime;
}

void crono_stop()
{
	struct rusage rusage;
	
	getrusage(RUSAGE_SELF, &rusage);
	stop_utime = rusage.ru_utime;
}

double crono_ms()
{
	return (stop_utime.tv_sec - start_utime.tv_sec) * 1000 +
    (stop_utime.tv_usec - start_utime.tv_usec) / 1000 ;
}


int main( int argc, char **argv )
{
	if(argc < 2){
		perror("Nombre d'arguments incorrect\n");
		exit(EXIT_FAILURE);
	}
	
	/*Initialisation des données du problème*/
	donnees donprob;
	
	
	/*Nombre de variables x_ij (var binaire : passage de i à j)*/
	int nbvar;
	/*Nombre de contraintes ajoutées pour obtenir une solution composée d'un unique cycle*/
	int nbcontr;

	/*Initialisation des données et de la matrice creuse pour GLPK */
	glp_prob *prob;
	int *ia;
	int *ja;
	double *ar;
	
	/*Resultats*/
	double z; // Res fonction objectif
	double *x; // Variables
	int *xcast; // Variables castées en entiers (pour le traitement)

	/*Nombre de coeff non-nuls dans la matrice de contraintes*/
	int nbcreux;
	
	/*Indices de boucles utilisés ici et là*/
	int i,j=0;
	int pos;
	
	/* Tableau contenant le sous-cycle à casser */
	int *tab_cycle; 
	
	double temps;
	int nbsol = 0; /* Compteur du nombre d'appels au solveur GLPK */ 
	
	/*Chargement des données d'apreès le fichier*/
	shiny_reader( argv[1], &donprob);
	
	crono_start(); /* Lancement du compteur*/

	nbvar = donprob.n*donprob.n;
	nbcontr = 2*donprob.n;
	nbcreux = (2*donprob.n*donprob.n)-2*donprob.n;
	printf("nbvar=%d, nbcontr=%d, nbcreux=%d \n", nbvar, nbcontr, nbcreux);

	prob = glp_create_prob(); /*Allocation mémoire pour le problème*/
	glp_set_prob_name(prob, "trajet"); /* affectation d'un nom */
	glp_set_obj_dir(prob, GLP_MIN); /* Il s'agit d'un problème de minimisation */
	
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_OFF; /*Paramètre de GLPK dans la résolution d'un PL en variables continues afin de couper les affichages (dans lesquels on se noierait) */
	glp_iocp parmip;
	glp_init_iocp(&parmip);
	parmip.msg_lev = GLP_MSG_OFF;/* Paramètre de GLPK dans la résolution d'un PL en variables entières afin de couper les affichages (dans lesquels on se noierait) */

	
	/* Déclaration du nombre de contraintes (nombre de lignes de la matrice des contraintes) */
	glp_add_rows(prob, nbcontr); 


	/* Bornes sur les contraintes */
	for(i=1;i<=nbcontr;i++)
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
		/* bornes sur les variables, ///me sur les contraintes */
		glp_set_col_kind(prob, i, GLP_BV);
	} 
	
	/* définition des coefficients des variables dans la fonction objectif */
	for (i=1; i<=nbvar; i++) glp_set_obj_coef(prob,i,donprob.d[i-1]);

	/* allocation */
	ia = (int *) malloc((1 + nbcreux) * sizeof(int));
	ja = (int *) malloc((1 + nbcreux) * sizeof(int));
	ar = (double *) malloc((1 + nbcreux) * sizeof(double));

	/*remplissage matrice contrainte*/
	pos = 1;

	for(i=1; i<=donprob.n; i++){
		for (j=1; j<=donprob.n; j++){
			if(i!=j){
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
	}
	
	/*Check matrice de contrainte*/
	/*
	for(i=1; i<=nbcreux; i++)
	{
		printf("ia[%d] = %d; ja[%d] = %d; ar[%d] = %f;\n", i, ia[i], i, ja[i], i, ar[i]);
	}*/
	
	/* chargement de la matrice dans le problème */
	glp_load_matrix(prob,nbcreux,ia,ja,ar); 
	
	/* écriture de la modélisation dans un fichier*/
	///glp_write_lp(prob,NULL,"trajet.lp");

	/* Résolution, puis lecture des résultats */
	nbsol++;
	glp_simplex(prob,&parm);
	glp_intopt(prob,&parmip); /* Résolution */
	
	z = glp_mip_obj_val(prob); /* Récupération de la valeur optimale.*/
	
	/* Récupération des valeurs de variables */
	x = (double *) malloc(nbvar * sizeof(double));
	for(i = 0;i < nbvar; i++) x[i] = glp_mip_col_val(prob,i+1); 

	/// printf("z = %lf\n",z);
	
	xcast = (int *) malloc(nbvar * sizeof(int));
	int xnew;
	for(i = 0;i < donprob.n;i++){
		for(j = 0;j < donprob.n;j++){
		xnew=(int)(x[i*donprob.n+j] + 0.5);
		/* un cast est ajouté, x[i] pourrait être égal à 0.99999... */
		/// printf("x%d_%d = %d, ",i,j,xnew);
		xcast[i*donprob.n+j]=xnew;
		}
		/// puts("");
	}
	
	tab_cycle = (int *) calloc( sizeof(int), donprob.n * sizeof(int)); // Allocation dans fonction

	// Récupérer le cycle dans tab_cycle
	int long_cycle = plus_petit_cycle(xcast, donprob.n, tab_cycle);
	/// printf("Cycle de taille %d\n",long_cycle);
	
	while(long_cycle < donprob.n)
	{
		
		nbcreux += long_cycle;		
		
		// ****ajout contraintes
		
		/* Ajout d'une ligne de contrainte */
		glp_add_rows(prob, 1); 
		nbcontr++; 

		/* Bornes sur la nouvelle contrainte pour casser le cycle */
		glp_set_row_bnds(prob, nbcontr, GLP_UP, long_cycle-1.0, long_cycle-1.0); 
		
		// ****realloc
		/* allocation */
		///printf("Nouveau nbcreux : %d\n",nbcreux);
		ia = (int *) realloc(ia, (100 + nbcreux)* sizeof(int));
		ja = (int *) realloc(ja, (100 + nbcreux)* sizeof(int));
		ar = (double *) realloc(ar, (100 + nbcreux)* sizeof(double));
			
		/* Remplissage matrice contrainte */
		/* Transition fin/début du cycle */
		ia[pos] = nbcontr;
		ja[pos] = ((tab_cycle[long_cycle-1])*donprob.n)+tab_cycle[0]+1;
		ar[pos] = 1.0;
		pos++;
			
		/* Autres transitions */
		for(i=0;i<long_cycle-1;i++){
			ia[pos] = nbcontr;
			ja[pos] = (tab_cycle[i]*donprob.n)+tab_cycle[i+1]+1;
			ar[pos] = 1.0;
			pos++;
		}
		
		// ****repeter la resolution
		/* A mettre dans une fonction */
			
		/* Check matrice de contrainte */
		for(i=1; i<=nbcreux; i++)
		{
			/// printf("ia[%d] = %d; ja[%d] = %d; ar[%d] = %f;\n", i, ia[i], i, ja[i], i, ar[i]);
		}
		/* chargement de la matrice dans le problème */
		glp_load_matrix(prob,nbcreux,ia,ja,ar); 
	
		/* écriture de la modélisation dans un fichier*/
		///glp_write_lp(prob,NULL,"trajet.lp");

		/* Résolution, puis lecture des résultats */
		nbsol++;
		glp_simplex(prob,&parm);
		glp_intopt(prob,&parmip); /* Résolution */
	
		z = glp_mip_obj_val(prob); /* Récupération de la valeur optimale.*/
	
		/* Récupération des valeurs de variables */
		x = (double *) malloc(nbvar * sizeof(double));
		for(i = 0;i < nbvar; i++) x[i] = glp_mip_col_val(prob,i+1); 

		/// printf("z = %lf\n",z);
	
		xcast = (int *) malloc(nbvar * sizeof(int));

		for(i = 0;i < donprob.n;i++){
			for(j = 0;j < donprob.n;j++){
			xnew=(int)(x[i*donprob.n+j] + 0.5);
			/* un cast est ajouté, x[i] pourrait être égal à 0.99999... */
			/// printf("x%d_%d = %d, ",i,j,xnew);
			xcast[i*donprob.n+j]=xnew;
			}
			/// puts("");
		}
		
		// recalcul pluspetitcycle
		long_cycle = plus_petit_cycle(xcast, donprob.n, tab_cycle );
	}
	
	/* Résolution achevée, arrêt du compteur de temps et affichage des résultats */
	crono_stop();
	temps = crono_ms()/1000,0;

	// Affichage de la solution
	printf("z = %lf\n",z);
	
	printf( "\n Itinéraire calculé :\n" );
	int xi = 0 ;
	int depart, destination;
	depart=0;
	for (xi = 0 ; xi < donprob.n ; xi++ )
	{
		destination = depart;
		//On cherche la valeur 1 sur la "ligne"
		while( xcast[ destination ] != 1 ) destination++;
		printf( "%d -> ", depart/donprob.n);
		depart = (destination % donprob.n) * donprob.n;
	}
	puts(" ");
	
	printf("Temps : %f\n",temps);	
	printf("Nombre d'appels à GPLK : %d\n",nbsol);
	printf("Nombre de contraintes ajoutées : %d (%d contraintes de base)\n",nbcontr, 2*donprob.n);

	/* libération mémoire */
	glp_delete_prob(prob);
	free(ia);
	free(ja);
	free(ar);
	free(x);
	free(xcast);
	free(tab_cycle);

	return EXIT_SUCCESS;
}

