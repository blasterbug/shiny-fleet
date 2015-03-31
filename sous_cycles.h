/**
 ** Théo DELALANDE-DELARBRE
 ** Benjamin SIENTZOFF
 ** 601A
 ** Shiny fleet
 ** Détection du plus-petit sous-cycle d'un ensemble de variables
 **/

#ifndef SOUSCYCLES
#define SOUSCYCLES

#include <stdio.h>
#include <stdlib.h>

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
#endif
/** fin du fichier **/
