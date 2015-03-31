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

/**
 * Trouve le plus petit cycle après résolution
 *	@param[in] valeursvar tableau des variables après résolution
 *	@param[in] nbdest nombre de destinations dans le problème
 *	@param[out] boucle_min tableau avec les cycles
 * 	@param[out] n Taille du plus petit cycle trouvé
 */
int plus_petit_cycle( const int* valeursvar, const int nbdest, int* boucle_min )
{
	int* boucle_courante = (int*) malloc((nbdest) * sizeof(int));
	int boucle_courante_long = 0;
	//boucle_min = (int*) malloc((nbdest) * sizeof(int)); // Bug ?! WTF!
	int boucle_min_long = nbdest;
	bool* desti_visit = (bool*) malloc((nbdest) * sizeof(bool));
	
	int i, j, suiv;
	int * desti_succ = (int *) malloc((nbdest) * sizeof(int));
	
	// initialisation des destinations à vérifier et des tableaux
	for (i=0; i< nbdest; i++){
		boucle_courante[i] = 0;
		boucle_min[i] = 0;
		desti_visit[i]=false;
	}
	
	// Récolte le successeur de chaque destination
	for (i=0; i< nbdest; i++){
		j=nbdest*i;
		//On cherche la valeur 1 sur la "ligne"
		while (valeursvar[j]!=1.)
		{ 
			j++;
		}
		desti_succ[i] = j % nbdest;
	}
	
	// Affichage des successeurs pour debug
	// for(i=0; i<nbdest; i++) printf("%d->%d\n",i,desti_succ[i]); 

	// recherche des cycles
	// chaque sommet
	for (i=0; i<nbdest; i++){
		// éliminer les sommets déjà visités
		if(!desti_visit[i]){
			suiv = i;
			// pour chaque somet visité
			// on en visite au moins un
			do{
				desti_visit[ suiv ] = true;
				// TODO: modulo taille du tableau
				boucle_courante[ boucle_courante_long ] = suiv;
				boucle_courante_long++;
				suiv = desti_succ[ suiv ];
			} while (suiv != i);
			
			// Affichage du cycle détecté pour debug
			// for (j=0; j< boucle_courante_long; j++) printf("%d,",boucle_courante[j]);
			// puts(" ");
			
			// Si le cycle détecté est plus petit que l'ancien cycle
			if (boucle_courante_long<boucle_min_long){
				// alors on le conserve
				for (j=0; j< nbdest	; j++){
					boucle_min[j] = boucle_courante[j];
				}
				// MàJ taille du cycle
				boucle_min_long = boucle_courante_long;
			}
			
			// Nettoie boucle_courante
			for (j=0; j< boucle_courante_long; j++){
				boucle_courante[j] = 0;
			}
			boucle_courante_long=0;
		}
	}
	
	// Affichage du plus petit cycle détecté pour debug
	/// printf("Plus petit cycle détecté :");
	/// for (j=0; j< boucle_min_long; j++) printf("%d, ", boucle_min[j]);
	/// puts(" ");
	
	// libération de la mémoire
	free(boucle_courante);
	free(desti_visit);
	free(desti_succ);
	
	return boucle_min_long;
}

int main( int argc, char **argv )
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
	double z; // Res fonction objectif
	double *x; // Variables
	int *xcast; // Variables castées en entiers (pour le traitement)

	/*Nombre de coeff non-nuls dans la matrice de contraintes*/
	int nbcreux;
	
	/*Indices de boucles utilisés ici et là*/
	int i,j=0;
	int pos;
	
	int *tab_cycle; // Tableau contenant le sous-cycle à casser 

	nbvar = donprob.n*donprob.n;
	nbcont = 2*donprob.n;
	nbcreux = 2*donprob.n*donprob.n;
	printf("nbvar=%d, nbcont=%d, nbcreux=%d \n", nbvar, nbcont, nbcreux);

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
	for(i=1; i<=nbcreux; i++)
	{
		printf("ia[%d] = %d; ja[%d] = %d; ar[%d] = %f;\n", i, ia[i], i, ja[i], i, ar[i]);
	}*/
	
	/* chargement de la matrice dans le problème */
	glp_load_matrix(prob,nbcreux,ia,ja,ar); 
	
	/* écriture de la modélisation dans un fichier*/
	///glp_write_lp(prob,NULL,"trajet.lp");

	/* Résolution, puis lecture des résultats */
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
	int long_plus_petit_cycle = plus_petit_cycle(xcast, donprob.n, tab_cycle);
	/// printf("Cycle de taille %d\n",long_plus_petit_cycle);
	
	while(long_plus_petit_cycle < donprob.n)
	{
		nbcont++; 
		nbcreux += long_plus_petit_cycle;
		
		
		
		// ****ajout contraintes
		
		/* Ajout d'une ligne de contrainte */
		glp_add_rows(prob, 1); 


		/* Bornes sur la nouvelle contrainte pour casser le cycle */
		glp_set_row_bnds(prob, nbcont, GLP_UP, long_plus_petit_cycle-1.0, long_plus_petit_cycle-1.0); 
		
		// ****realloc
		/* allocation */
		///printf("Nouveau nbcreux : %d\n",nbcreux);
		ia = (int *) realloc(ia, (100 + nbcreux)* sizeof(int));
		ja = (int *) realloc(ja, (100 + nbcreux)* sizeof(int));
		ar = (double *) realloc(ar, (100 + nbcreux)* sizeof(double));
			
		/* Remplissage matrice contrainte */
		/* Transition fin/début du cycle */
		ia[pos] = nbcont;
		ja[pos] = ((tab_cycle[long_plus_petit_cycle-1])*donprob.n)+tab_cycle[0]+1;
		ar[pos] = 1.0;
		pos++;
			
		/* Autres transitions */
		for(i=0;i<long_plus_petit_cycle-1;i++){
			ia[pos] = nbcont;
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
		long_plus_petit_cycle = plus_petit_cycle(xcast, donprob.n, tab_cycle );
	}

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

