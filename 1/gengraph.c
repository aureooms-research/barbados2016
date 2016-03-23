// -*- coding: utf-8 -*-
/*


   Graph Generator                                   © Cyril Gavoille


     Use: gengraph [-option] graph_name [parameter_list]


     A free command-line program to generate graphs in many formats:
     plain text, .dot, .pdf, .fig, .svg ... Formats like .pdf or .fig
     are produced thanks to GraphViz and allow visualization.  Easy to
     install, there is a single .c file to compile. There is a on-line
     manual available in the command (for the moment the manual is in
     French only and included at the end of the source).

                                -----

     Petit programme libre en ligne de commande permettant de générer
     des graphes dans plein de formats différents: texte, .dot, .pdf,
     .fig, .svg ... Certains formats (comme .pdf ou .jpg), produits
     grâce à GraphViz, permettent de visualiser les graphes.  Très
     simple à installer, il n'y a qu'un seul fichier .c à compiler.
     Pour le reste il faut voir l'aide en ligne qui est aussi incluse
     à la fin du source.


   Comment l'installer / le compiler ?

     (MacOs)  gcc gengraph.c -o gengraph
     (Linux)  gcc gengraph.c -lm -o gengraph

*/

#include <stdio.h>  /* pour printf(), sprintf() ... */
#include <stdlib.h> /* pour system(), strtod(), RAND_MAX ... */
#include <string.h> /* pour strcomp() ... */
#include <unistd.h> /* pour getpid() ... */
#include <math.h>   /* pour sqrt(), cos(), acos() ... */
#include <float.h>  /* pour DBL_MAX, DBL_DIG ... */
#include <limits.h> /* pour INT_MAX, LONG_MAX ... */
#include <time.h>   /* pour time() */
#include <assert.h> /* pour assert() */


#define DEBUG(I) //if(1){ I }

/* types */

/* graphe ou famille de graphes */
typedef struct _graph {
  int id;    // numéro du graphe, utilisé pour les familles de graphes 
  int n;     // nb de sommets, <0 si non défini
  int m;     // nb d'arêtes, <0 si non défini
  int sort;  // vrai ssi les listes d'adjacences sont triées
  int *d;    // d[u]=degré du sommet u
  int **L;   // L[u][i]=i-ème voisin du sommet u, i=0...d[u]-1
  int sym;   // vrai ssi les listes d'adjacence du graphe sont symétriques
  double **W;// W[u][i]=poids du i-ème voisin du sommet u, i=0...d[u]-1
  double *xpos,*ypos; // tableau de positions des sommets (graphe géométrique)
  int f;     // nombre de graphes de la famille, =0 si graphe normal
  struct _graph **G; // G[i]=i-ème graphe de la famille, i=0..f-1

  // les champs suivants sont utilisés pour communiquer des valeurs
  // ou des résultats à travers les appels de fonctions

  int int1;  // paramètre: entier
  int* pint1;// paramètre: tableau d'entiers
} graph;

typedef int test(graph*); /* fonction de test sur un graphe */

/* chemin simple d'un graphe G */
typedef struct{
  int n;  /* nombre de sommets du chemin */
  int *P; /* liste ordonnée des sommets du chemin */
  int *V; /* V[u]=i si le sommet u de G est le i-ème (i dans [0,G->n[)
	     sommet du chemin (V[P[i]]=i), V[u]=-1 si u n'est dans le
	     chemin  */
  int **aux; /* tableau auxiliaire utilisé (et géré) par NextPath() */
  int nG;    /* nombre de sommets du graphe G, sert pour free_path() */
} path;

/* code pour le type "list" */

enum{
  T_NODE,  // item est un sommet u
  T_EDGE,  // item est un sommet v connecté au précédent de la liste par une arête: u-v
  T_ARC,   // item est un sommet v connecté au précédent de la liste par un arc: u->v
  T_ID,    // item est l'identifiant d'un nouveau graphe d'une famille
  T_NB,    // item est le nombre de graphes de la famille, item toujours en tête de liste
  T_OPENE, // item est le premier sommet v d'un groupe arête: u-(v ...
  T_OPENA, // item est le premier sommet v d'un groupe arc: u->(v ...
  /* token non encore géré */
  T_UNIV,  // u-* = item est un sommet universel u-(0 ... n)
  T_UNOUT, // u->* item est un sommet universel u->(0 ... n-1)
  T_UNIN,  // *->u item est un sommet universel u<-(0 ... n-1)
};

/* code pour les fonctions de hashage */

enum{
  H_PRIME,
  H_SHUFFLE,
  H_MOD,
};

typedef int adjacence(int,int); /* fonction d'adjacence */
typedef struct _list { int item,type; struct _list *next; } list; /* liste chaînée */
typedef struct{ double x,y; } point; /* un point du plan */
typedef struct{ int u,v; } edge; /* une arête */
typedef struct{ int x,y; } couple; /* une paire d'entiers */
typedef struct{ int x,y; long z; } triplet; /* un triplet d'entiers (avec lon long) */
typedef struct{ unsigned int r:8,g:8,b:8; } color;
/* une couleur est un int accessible par les champs r,g,b de 8 bits */

/* constantes pour la ligne de commande */

#define DIMAX   256  /* nb maximum de paramètres pour un graphe */
#define CMDMAX  1024 /* nb maximum de caractères sur la ligne de commande */
#define NAMEMAX 20   /* nb maximum de caractères d'un nom de sommet */
#define PARAMSIZE 1024 /* taille mémoire des buffers FPARAM et CPARAM (en octets) */

/* constantes pour le format dot (-visu) */

double VSIZESTD=0.05; /* taille standard des sommets */
double VSIZEXY=0.12;  /* taille des sommets pour les graphes géométriques */
double LEN=1.0; /* longueur des arêtes avec dot neato */
#define GRAPHPDF g.pdf /* nom du graphe de sortie par défaut */
#define VSIZEK 5.5 /* rapport entre taille max et min des sommets (-vsize) */
#define C32 32 /* Coefficient par lequel sont multipliées les
	          coordonnées des points dans le cas des graphes
	          géométriques et le format dot. Plus la valeur est
	          élevée, plus les sommets paraissent petits et les
	          arêtes fines. */
#define STR_DISTRIB "▪" /* pour affichage des distributions */
#define LEN_DISTRIB (60.0) /* longueur max pour l'affiche des distributions */
#define LEN_BARRE 80 /* longueur d'une barre */
#define CHRONOMAX 5 /* nombre maximum de chronomètres */
#define RS_NI_FP "name-independent, fixed-port model"
#define RS_L_FP  "labeled, fixed-port model"

/* codes pour les formats de sorties possibles */

enum{
  F_standard,
  F_list,
  F_matrix,
  F_smatrix,
  F_dot,
  F_xy,
  F_no
};

/* codes pour les algorithmes possibles. La variable CHECK vaut l'une
   de ces constantes. Par défaut, CHECK=CHECK_OFF(=0). Si CHECK>0,
   alors le graphe sera stocké en mémoire. Si CHECK>1, alors en plus
   l'algorithme spécifique sera appliqué. */

enum{
  CHECK_OFF=0, // valeur par défaut
  CHECK_ON,    // graphe chargé en mémoire
  CHECK_BFS,
  CHECK_DFS,
  CHECK_NCC,
  CHECK_BELLMAN,
  CHECK_DEGENERATE,
  CHECK_GCOLOR,
  CHECK_DEG,
  CHECK_ISO,
  CHECK_SUB,
  CHECK_ISUB,
  CHECK_MINOR,
  CHECK_TWDEG,
  CHECK_TW,
  CHECK_PS1,
  CHECK_PS1b,
  CHECK_PS1c,
  CHECK_PS1x,
  CHECK_GIRTH,
  CHECK_RADIUS,
  CHECK_DIAMETER,
  CHECK_PATHS,
  CHECK_PATHS2,
  CHECK_MAINCC,
  CHECK_KCOLOR,
  CHECK_KCOLORSAT,
  CHECK_KINDEPSAT,
  CHECK_INFO,
  CHECK_RS_CLUSTER,
  CHECK_RS_DCR,
  CHECK_RS_TZRPLG,
  CHECK_RS_BC,
};

/* codes pour la variable XYtype indiquant le type de génération des
   points (XPOS,YPOS) des graphes géométriques. */

enum{
  XY_FILE, // points à partir d'un fichier
  XY_UNIF, // points selon une loi uniforme
  XY_PLAW, // points selon une loi en puissance
  XY_PERM, // points selon une permutation π, (i,π(i))
  XY_MESH, // points sur une grille Xmesh x Ymesh
  XY_USER, // points définis par la fonction d'adjacence
};


/* fonctions inline */

/* Attention!: utiliser fmin() pour des doubles, car
   min(double,double) ne provoque par de warning */

static inline int min(int i,int j){ return (i<j)?i:j; }
static inline int max(int i,int j){ return (i>j)?i:j; }

/* macros & pseudo-fonction utiles */

#define EQUAL(s)   (strcmp(ARGV[i],s)==0)
#define PREFIX(s)  (strncmp(ARGV[i],s,strlen(s))==0)
#define PLURIEL(n) (((n)>1)?"s":"")
#define SWAP(x,y,z)  (z)=(x),(x)=(y),(y)=(z) /* échange x et y, via z */
#define XORSWAP(x,y) ((x)^=(y),(y)^=(x),(x)^=(y)) /* comme SWAP() mais sans z */
#define VIDE(s)    *s='\0'    /* vide le tableau de caractères s */ 
#define ESTVIDE(s) (*s=='\0') /* vrai ssi s est vide */
#define NONVIDE(s) (*s!='\0') /* vrai ssi s est non vide */
#define MEM(mem,pos,type) (*(type*)(mem+pos)) /* écrit/lit la mémoire mem[pos] */
#define RAND01  ((double)random()/(double)RAND_MAX) /* réel aléatoire dans [0,1[ */

/* permet l'expansion d'une macro. Ex: scanf("%"xstr(DMAX)"s",buffer); */
#define xstr(s) str(s)
#define str(s) #s

/* Alloue à la variable T un tableau de n valeurs du type de T[] */
#define _ALLOC(T,n) if(!((T)=malloc((n)*sizeof(*(T))))) Erreur(3)
#define ALLOC(T,n) do{ _ALLOC(T,n); } while(0)

/* Comme ALLOC mais initialise en plus le tableau avec le term z */
#define ALLOCZ(T,n,z)						\
  do{								\
    _ALLOC(T,n);						\
    int _i;for(_i=0;_i<(n);_i++) (T)[_i]=(z);			\
  }while(0)

/* Un realloc() qui évite de faire un free() non souhaité à un
   pointeur (si n=0), et donc qui évite les erreurs "double free". */
#define REALLOC(P,n) P=realloc(P,max(n,1)*sizeof(*P))

/* Alloue à la variable T une matrice de n x s valeurs, c'est-à-dire
   un tableau de n tableaux de s valeurs du type de T[][]. */
#define ALLOCMAT(T,n,s)			\
  do{					\
    int _k=0,_n=n;			\
    _ALLOC(T,_n);			\
    for(;_k<_n;_k++) _ALLOC(T[_k],s);	\
  }while(0)

#define ALLOCMATZ(T,n,s)		\
  do{					\
    int _k=0,_n=n;			\
    _ALLOC(T,_n);			\
    for(;_k<_n;_k++) ALLOCZ(T[_k],s,z);	\
  }while(0)

/* Les variantes NALLOC permettent, en plus de ALLOC, de déclarer le
   pointeur T sur des valeurs de type t.

   Ex: int *T;      ->   NALLOC(int,T,n);
       ALLOC(T,n);

   Attention de ne pas mettre ces macros dans un bloc { ... }, car
   sinon la portée de T se limitera à ce block. Les variantes NALLOC
   ne peuvent remplacer des ALLOC qui sont dans les configurations
   suivantes:

   - la déclaration et le malloc() ne sont pas dans le même bloc, comme:
     int *T;
     if(...)
       ALLOC(T,...);

   - la déclaration est après l'étiquette d'un goto ou d'un case, comme:
     label:
       int *T;
       ALLOC(T,...);
*/
#define NALLOC(t,T,n) t *T;ALLOC(T,n);
#define NALLOCZ(t,T,n,z) t *T;ALLOCZ(T,n,z);
#define NALLOCMAT(t,T,n,s) t **T;ALLOCMAT(T,n,s);
#define NALLOCMATZ(t,T,n,s,z) t **T;ALLOCMAT(T,n,s,z);

/* Libère un tableau de pointeurs. On ne peut pas faire de fonction,
   le type de X n'étant pas pré-déterminé on ne peut écrire le
   prototype de la fonction. */
#define FREE2(X,n)			\
  do{					\
    if(X){				\
      int _j=0;				\
      for(;_j<(n);_j++) free(X[_j]);	\
      free(X);				\
    }					\
  }while(0)

/* Comme FREE2() mais pour un tableau de dimension 3.
   Ex: FREE3(X,n,m[_i]); */
#define FREE3(X,n,m)			\
  do{					\
    if(X){				\
      int _i=0;				\
      for(;_i<(n);_i++) FREE2(X[_i],m);	\
      free(X);				\
    }					\
  }while(0)

/* ajoute une arête uv au graphe G en supposant qu'il y a assez de place */
#define ADD_EDGE(G,u,v)	   	\
  do{				\
    G->L[u][G->d[u]++]=v;	\
    G->L[v][G->d[v]++]=u;	\
  }while(0)

/* comme ADD_EDGE */
#define ADD_ARC(G,u)	   	\
  do{				\
    G->L[u][G->d[u]++]=v;	\
  }while(0)

/* macros de débuggage */

/* affiche une variable v de type int */
#define PRINT(v)					\
  do{							\
    printf(#v " = %i\n",v);				\
  }while(0)

/* affiche une variable v de type pointeur */
#define PRINTP(v)					\
  do{							\
    printf(#v " = %p\n",v);				\
  }while(0)

/* affiche un tableau d'entiers */
#define PRINTT(T,n)					\
  do{							\
    int _i;						\
    printf(#T " =");					\
    if((T)==NULL) printf(" NULL"); else			\
      for(_i=0;_i<(n);_i++) printf(" %i",(T)[_i]);	\
    printf("\n");					\
  }while(0)

/* affiche une liste chaînée */
#define PRINTLIST(L)		        	\
  do{						\
    list *_L=L;                                 \
    printf(#L " = {");				\
    while(_L!=NULL){				\
      printf(" (%i,%i)",_L->item,_L->type);     \
      _L=_L->next;                              \
      } printf(" }\n");				\
  }while(0)

// affiche un aperçu d'une liste de n valeurs
// ex: APERCU(_i*_i,10,3,2) affichera
// "[ 0 1 2 ... 64 81 ]"
#define APERCU(F,n,K1,K2)						\
  do{									\
    int _i;								\
    printf("[ ");							\
    for(_i=0;_i<n;_i++)							\
      if((_i<K1)||(_i>n-K2-1)||(K1==n-K2-1))				\
	printf("%i ",F);						\
      else{ if((_i>=K1)&&(_i<n-K2)){ printf("... "); _i=n-K2-1; }}	\
    printf("]\n");							\
  }while(0)

#define PAUSE scanf("%*c") /* appuyer sur [RETURN] pour continuer */
#define STRTOI(s) ((int)strtol(s,NULL,10)) 
#define STRTOD(s) strtod(s,NULL)

/* variables globales */

adjacence *adj;  /* nom de la fonction d'adjacence */
int N;           /* N=nb de sommets du graphes avant suppression */
int NF;          /* nb de sommets final du graphes (donc après suppression) */
int *V;          /* étiquette des sommets, en principe V[i]=i */
int *VF;         /* VF[j] est l'indice i du j-ème sommet non supprimé */
int *INC;        /* INC[i]=deg(i). Si =0, alors i est un sommet isolé */
int ARGC;        /* variable globale pour argc */
char **ARGV;     /* variable globale pour argv */
int PARAM[DIMAX];/* liste des paramètres du graphe */
double DPARAM[DIMAX]; /* idem PARAM, mais si paramètres de type double */
char SPARAM[64]; /* pour un paramètre de graphe de type chaîne de caractères */
char PARAM_PAL[64]; /* mot définissant le dégradé de couleur pour -vcolor pal */
void* CPARAM=NULL; /* liste de paramètres (pointeur tout type, en octets) pour -check */
void* FPARAM=NULL; /* liste de paramètres (pointeur tout type, en octets) pour -filter */
int CVALUE;      /* sert pour la valeur dans -filter */
int PVALUE;      /* =1 ssi on affiche la valeur du paramètre dans -filter */
test *FTEST;     /* pour l'option -filter */
double DELE=0.0; /* proba de supprimer une arêtes */
double DELV=0.0; /* proba de supprimer un sommet */
double REDIRECT=0.0; /* proba de rediriger une arête */
int SHIFT=0;     /* début de la numérotation des sommets */
int DIRECTED=0;  /* vrai ssi graphe orienté for(i=0 ... for(j=...) ) */
int LOOP=0;      /* vrai ssi ont n'affiche pas les boucles */
int PERMUTE=0;   /* vrai ssi -permute */
int NOT=0;       /* vrai ssi -not */
int POS=0;       /* vrai ssi -pos */
int FAST=0;      /* vrai ssi -fast */
int CHECK=0;     /* vrai ssi option -check */
int VARIANT=0;   /* variante de l'option -variant */
int STAR=0;      /* paramètre de -star */
int APEX=0;      /* paramètre de -apex */
int LABEL=0;     /* vrai ssi affiche les labels des sommets (-label) */
int NORM=2;      /* norme pour les graphes géométriques: L_2 par défaut */
int FORMAT=F_standard; /* type de la sortie, par défaut le format standard */
int HEADER=0;    /* par défaut pas de préambule, sinon HEADER=1 */
int WIDTH=12;    /* nb maximum d'arêtes ou de sommets isolés affichés par ligne */
unsigned SEED;   /* graîne du générateur aléatoire */
char NAME[NAMEMAX]; /* nom du sommet retourné par adj(i,ADJ_NAME) */
char *DOTFILTER="neato"; /* nom du filtre "dot" par défaut */
char *CAPTION=NULL; /* légende du graphe */
char *STRsep1;   /* séparateur d'arête: "" (pour standard) ou ";" (dot) */
char *STRsep2;   /* autre séparateur */
char *STRedge;   /* caractère pour une arête: "-" (pour standard) ou "--" (dot) */

char *FILEXY;    /* nom de fichier des coordonnées */
double *XPOS=NULL,*YPOS=NULL; /* tableaux (X,Y) pour les graphes géométriques */
double BOXX=-1.0,BOXY; /* pour le redimensionement: <0 signifie pas de redim. */
double NOISEr=-1.0,NOISEp; /* paramètre pour -xy noise: <0 signifie pas de "noise" */
double XMIN=0,YMIN=0,XMAX=1,YMAX=1; /* Bounding Box par défaut */
int ROUND=DBL_DIG; /* arrondi de XPOS/YPOS à 10^-ROUND près */
/* DBL_DIG (=15) est considérée comme valeur impossible pour ROUND et signifie aucun arrondi */
int XYtype=XY_UNIF; /* type de génération des points, uniforme par défaut */
int XYunique=0; /* =1 ssi on élimine les points doubles */
int XYgrid=0; /* =1 ssi on affiche une grille grisée */
int Xmesh=0,Ymesh=0; /* dimension de la grille pour l'option -xy mesh */
double XYvsize=1.0; /* facteur d'échelle pour la taille des sommets dans F_dot */
int SEEDk;  /* nombre de graînes pour génération des points */
double SEEDp; /* loi puissance pour les graines pour génération des points */
double *XSEED=NULL,*YSEED=NULL; /* tableaux de doubles pour les graînes */

int WRAP[DIMAX]; /* WRAP[k]=1 ssi la dimension k est "torique" */
int LOADC=0; /* vrai ssi graphe "loadc file" */
graph *LOAD=NULL;   /* graphe auxilière pour "load file", "bdrg ..." */
graph *GF=NULL;     /* graphe pour l'option -check */
graph *FAMILY=NULL; /* graphe pour l'option -filter */
int **REP; /* Associe à chaque sommet i une représentation sous forme
              de tableau d'entiers (REP[i] est donc un tabeau). Sert
              pour les représentations implicites de graphes, les
              arbres, les graphes d'intersections, etc. */
double **DREP; /* comme REP mais avec des doubles */
int HASH=H_PRIME; /* fonction de hashage par défaut */
int VSIZE=0; /* code pour la taille des sommets */
int VCOLOR=0; /* code pour la couleur les sommets */
char COLORCHAR[]="redjykugfocatbhsqvmpinzxlw";
color COLORBASE[]={ /* HTML Color Names */
  /* palette des couleurs de bases, l'ordre doit être celui de COLORCHAR */
  {255,  0,  0}, // r=red
  {210,105, 30}, // e=chocolate
  {255,140,  0}, // d=darkorange
  {255,165,  0}, // j=orange
  {255,255,  0}, // y=yellow
  {240,230,140}, // k=khaki
  {154,205, 50}, // u=yellowgreen
  {  0,255,  0}, // g=green (lime)
  { 34,139, 34}, // f=forestgreen
  {128,128,  0}, // o=olive
  {  0,255,255}, // c=cyan
  {127,255,212}, // a=aquamarine
  {  0,128,128}, // t=teal
  {  0,  0,255}, // b=blue
  {255,105,180}, // h=hotpink
  {250,128,114}, // s=salmon
  {255,192,203}, // q=pink
  {238,130,238}, // v=violet
  {255,  0,255}, // m=magenta
  {128,  0,128}, // p=purple
  { 75,  0,130}, // i=indigo
  {  0,  0,128}, // n=navy
  {  0,  0,  0}, // z=black
  {128,128,128}, // x=gray
  {230,230,250}, // l=lavender
  {255,255,255}  // w=white
};
color *PALETTE=COLORBASE; /* palette par défaut */
#define COLORNB    ((int)sizeof(COLORCHAR)-1) /* nb de couleurs de la palette de base */
int NPAL=COLORNB; /* taille de la palette par défaut */
struct{ int mode,u,v; } SCENARIO; /* scenario pour l'option "-check routing" */

/***********************************

         ROUTINES EN VRAC

***********************************/


void Erreur(int erreur){
  char *s;

  switch(erreur){
  case  1: s="Erreur : option -xy non reconnue."; break;
  case  2: s="Erreur : option non reconnue."; break;
  case  3: s="Erreur : espace mémoire insuffisant."; break;
  case  4: s="Erreur : nombre trop grand de paramètres."; break;
  case  5: s="Erreur : format de sortie inconnu."; break;
  case  6: s="Erreur : paramètre incorrect."; break;
  case  7: s="Erreur : ouverture du fichier impossible."; break;
  case  8: s="Erreur : tableau de coordonnées inexistant."; break;
  case  9: s="Erreur : option -vcolor non reconnue."; break;
  case 10: s="Erreur : nom de graphe inconnu ou erreur dans les paramètres."; break;
  case 11: s="Erreur : le graphe doit être connexe."; break;
  case 12: s="Erreur : option -check non reconnue."; break;
  case 13: s="Erreur : format de famille de graphes invalide."; break;
  case 14: s="Erreur : option -filter non reconnue."; break;
  case 15: s="Erreur : graphe(s) non trouvé(s)."; break;
  case 16: s="Erreur : plage de valeurs incorrecte."; break;
  case 17: s="Erreur : nom de sommets trop grand."; break;
  case 18: s="Erreur : opérateurs -star et -apex incompatibles."; break;
  case 19: s="Erreur : nom de fichier trop long."; break;
  case 20: s="Erreur : nombre de couleurs dans la palette trop important."; break;
  case 21: s="Erreur : code inconnue dans la fonction SortInt()."; break;
  case 22: s="Erreur : sommet de valeur négative dans le format standard."; break;
  case 23: s="Erreur : code inconnu dans la fonction routing_test()."; break;
  case 24: s="Erreur : option -visu incompatible avec -format no."; break;
  case 25: s="Erreur : la variante -fast n'est pas implémentée pour ce graphe."; break;
  case 26: s="Erreur : numéro de chronomètre incorrect."; break;
  case 27: s="Erreur : plusieurs options -check sur la ligne de commande."; break;
  case 28: s="Erreur : mauvais format du fichier d'entrée."; break;
  case 29: s="Erreur : loadc et -not/-permute sont incompatibles."; break;
  case 30: s="Erreur : loadc devrait être suivi d'un option -check."; break;
  case 31: s="Erreur : le graphe doit être simple et non-orienté, essayez -check info."; break;
  case 32: s="Erreur : problème dans la génération du graphe."; break;
  case 33: s="Erreur : dépassement arithmétique."; break;
  case 34: s="Erreur : la séquence de degré n'est pas graphique."; break;
  case 35: s="Erreur : -caption de devrait contenir qu'une occurence de format %XXX."; break;
  case 36: s="Erreur: le graphe doit comporter au moins une arête.";
  default: s="Erreur : code d'erreur inconnue."; /* ne devrait jamais arriver */
  }

  printf("%s\n",s);
  exit(EXIT_FAILURE);
}


/* codes pour GraphFromArray() */
enum{
  GFA_END   =-1, // fin du graphe
  GFA_CUT   =-2, // fin d'une séquence
  GFA_PATH  =-3, // chemin de sommets consécutifs
  GFA_HAM   =-4, // cycle Hamiltonien
  GFA_STAR  =-5, // étoile
  GFA_WHEEL =-6, // camembert
};


int GraphFromArray(int i,int j,int *T){
/*
  Fonction générique permttant de générer des graphes sans paramètre
  de petite taille. Plutôt que d'utiliser une matrice, c'est un
  tableau T qui définit les adjacences du graphe G.  En gros,
  GraphFromArray(i,j,T)=1 (soit "i-j") ssi i et j se suivent dans le
  tableau T, c'est-à-dire s'il existe un indice k tel que T[k]=i et
  T[k+1]=j ou le contraire.  Les arêtes du graphe sont ainsi couvertes
  par des chemins éventuellement non élémentaires.

  Les valeurs négatives du tableau ont des significations
  particulières (voir ci-dessous la signification des codes
  GFA_xxx). Le principe de l'algorithme de décodage est le suivant: on
  lit le tableau progressivement valeur après valeur, et soit x le
  dernier sommet mémorisé, c'est-à-dire une valeur entière >=0. Si on
  lit une valeur y>=0 alors y est un sommet et le graphe possède
  l'arête x-y. On remplace alors x par y et on recommence. Si y<0,
  c'est un code GFA_xxx et le comportement est alors spécifique
  comme décrit ci-dessous:

  GFA_END: fin du tableau et donc du graphe.

  GFA_CUT: fin d'une séquence. Le prochain sommet lu ne sera pas
     connecté au dernier sommet mémorisé du tableau. Sans ce code, la
     prochaine valeur >=0 lue est toujours connectée au dernier sommet
     mémorisé. En fait, toute valeur <0 qui n'est pas un des codes
     reconnue à le même effet que GFA_CUT.

  GFA_PATH: chemin de sommets consécutifs, "a,GFA_PATH,n" représente
      un chemin de n arêtes suivant le sommet a. C'est équivalent à
      "a,a+1,...,a+n,GFA_CUT".

  GFA_HAM: cycle Hamiltonien, équivalent à "0,GFA_PATH,N-1,N-1,0,GFA_CUT".

  GFA_STAR: étoile. La séquence "a,GFA_STAR,a_1,...,a_n,GFA_CUT"
      représente une étoile de centre a et de feuilles
      a_1,...,a_n. C'est équivalent à
      "a,a_1,a,a_2,a,...,a,a_n,GFA_CUT". On peut remplacer le GFA_CUT
      par n'importe qu'elle valeur négative qui peut ainsi être
      enchaînée.

  GFA_WHEEL: comme GFA_STAR, sauf qu'en plus les sommets a_i,a_{i+1}
      sont adjacents.

  Une bonne stratégie pour trouver un code assez court pour un graphe
  (c'est-à-dire un tableau T de petite taille) est de l'éplucher par
  degré maximum jusqu'à obtenir des chemins ou des chemins de
  cycles. Les sommets isolés n'ont pas à être codés. Les sommets ainsi
  épluchés peuvent être codés par GFA_STAR ou GFA_WHEEL, les cycles,
  les chemins, et les chemins de cycles par GFA_PATH.

  Ex:          5-6-7
               |  /      -> code = 0,GFA_PATH,10,GFA_CUT,4,8,GFA_END
       0-1-2-3-4-8-9-10

*/

  int k=0,a,n;

  while(T[k]!=GFA_END){
    if(T[k]==GFA_PATH){
      a=T[k-1];
      n=T[++k];
      if((a<=i)&&(i<a+n)&&(j==(i+1))) return 1; /* suppose que i<j */
      goto next;
    }
    if(T[k]==GFA_HAM){
      if(j==(i+1)) return 1; /* suppose que i<j */
      if((i==0)&&(j==(N-1))) return 1;/* utilise N ici */
      goto next;
    }
    if(T[k]==GFA_STAR){
      a=T[k-1];
      while(T[++k]>=0){
	if((i==a)&&(j==T[k])) return 1;
	if((j==a)&&(i==T[k])) return 1;
      }
      continue;
    }
    if(T[k]==GFA_WHEEL){
      a=n=T[k-1];
      while(T[++k]>=0){
	if(((i==a)||(i==n))&&(j==T[k])) return 1;
	if(((j==a)||(j==n))&&(i==T[k])) return 1;
	n=T[k];
      }
      continue;
    }
    if((i==T[k])&&(j==T[k+1])) return 1;
    if((j==T[k])&&(i==T[k+1])) return 1;
  next:
    k++;
  }
  return 0;
}


void Permute(int *T,int n){
/*
  Permute aléatoirement les n premiers éléments de T.
*/
  int i,j,k;
  for(i=0;i<n;i++){
    j=i+(random()%(n-i));
    SWAP(T[i],T[j],k);
  }
}


void NAME_Vector(int *R,int n,char *sep,char *par,int d,char *f){
/*
  Ecrit dans la variable globale NAME le nom représenté par les n
  premiers entiers (lettres) du tableau R. Les lettres (chiffres) sont
  séparées par la chaîne "sep" (qui peut être vide). Le mot est
  parenthésé par la chaîne "par" qui peut soit vide soit doit être
  composé de deux caractères: par[0]=parenthèse gauche,
  par[1]=parenthèse droite. Si d>0, les lettres sont écrites dans le
  sens croissant des indices (vers la droite), sinon c'est dans le
  sens des indices décroissant. La chaîne f indique le format
  d'affichage (avec printf) de l'entier REP[i]. En général, il s'agit
  de "%i".

  Ex: R={3,6,1}, n=3, f="%i".
  d=1, sep="," et par="{}": alors NAME="{3,6,1}"
  d=1, sep="-" et par=""  : alors NAME="3-6-1"
  d=1, sep=""  et par=""  : alors NAME="361"
  d=0, sep=""  et par=""  : alors NAME="163"
*/
  int i,b,c;
  char s[NAMEMAX];
  char vide[2]="\0\0"; /* vide=00 */
  
  VIDE(NAME); /* vide la chaîne */
  if(d>0){ c=1;i=0;b=n; }else{ c=-1;i=n-1;b=-1; }

  /* parcoure R dans un sens ou l'autre */
  while(i!=b){
    sprintf(s,f,R[i]);
    strcat(NAME,s);
    i += c;
    if(i!=b) strcat(NAME,sep);
  }

  /* met les parenthèses */
  if(ESTVIDE(par)) par=vide;
  s[0]=par[0]; s[1]='\0';
  strcat(s,strcat(NAME,par+1));
  strcpy(NAME,s);
  return;
}


void NAME_Base(int u,int b,int n,char *sep,char *par,int d)
/*
  Comme NAME_Vector(...,n,sep,par,d,"%i") sauf que NAME est l'écriture
  de la valeur de u en base b.

  Ex: R=361, b=10, n=3 et d=1.
  sep="," et par="{}": alors NAME="{3,6,1}"
  sep="-" et par="": alors NAME="3-6-1"
  sep="" et par="": alors NAME="361"
*/
{
  int i;
  int R[NAMEMAX];

  if(b<2) return;
  for(i=0;i<n;i++){
    R[i]=u%b;
    u /= b;
  }

  NAME_Vector(R,n,sep,par,d,"%i");
  return;
}


int LoadXY(char *file){
/*
  Remplit les tableaux XPOS et YPOS à partir d'un fichier (file), et
  renvoie le nombre de points lues. Le programme peut être amené à
  lire un point de trop (et des coordonnées vides) si le "n" du
  fichier est > aux nombres de points qu'il contient réellement et que
  le fichier contient un retour de ligne finale.
*/
  FILE *f=stdin;
  int n,i;

  if(strcmp(file,"-")) f=fopen(file,"r"); /* ouvre en lecture */
  if(f==NULL) Erreur(7);
  i=fscanf(f,"%i",&n); /* lit n */
  if((n<0)||(i<0)) n=0;
  ALLOC(XPOS,n);
  ALLOC(YPOS,n);
  for(i=0;(i<n)&&(!feof(f));i++){
    if(fscanf(f,"//%*[^\n]\n")>0) continue;
    fscanf(f,"%lf %lf",XPOS+i,YPOS+i);
  }
  fclose(f);
  return i; /* i=minimum entre n et le nb de valeurs lues */
}


list *new_list(void){
/*
  Crée une nouvelle liste (qui est renvoyée) comprenant une célulle
  (sentinelle) avec le champs ->next=NULL.
*/
  NALLOC(list,L,1);
  L->next=NULL;
  return L;
}


list *Insert(list *p,int v,int t){
/*
  Écrit (v,t) dans l'élément de la liste chaînée pointée par p qui ne
  doit pas être NULL, crée un nouvel élément qui est chaîné à la suite
  de p. On renvoit le nouvel élément crée.
*/
  p->item=v;
  p->type=t;
  return p->next=new_list(); /* nouvel élément vide */
}


graph *new_graph(int n){
/*
  Renvoie un objet de type "graph". Les champs sont initialisés à
  leurs par défaut. Si n>0, alors les tableaux ->d et ->L de taille n
  sont alloués, mais pas les n tableaux ->L[u].
*/

  NALLOC(graph,G,1);
  G->n=0;
  G->m=-1;
  G->sort=0;
  G->id=-1;
  G->d=NULL;
  G->L=NULL;
  G->W=NULL;
  G->xpos=NULL;
  G->ypos=NULL;
  G->f=0;
  G->sym=1;
  G->G=NULL;

  G->pint1=NULL;
  G->int1=-1;

  if(n>0){
    G->n=n;
    ALLOC(G->d,n);
    ALLOC(G->L,n);
  }

  return G;
}


void free_graph(graph *G)
/*
  Libère G et tous ses tableaux.  Dans le cas d'une famille G, chaque
  graphe est aussi libéré (de manière récursive).
*/
{
  if(G==NULL) return;

  /* Remarque: ce n'est pas grave de faire free() sur un ptr NULL */

  int i;
  free(G->d);
  free(G->pint1);
  free(G->xpos);
  free(G->ypos);
  FREE2(G->L,G->n);
  FREE2(G->W,G->n);
  for(i=0;i<G->f;i++) free_graph(G->G[i]);
  free(G->G);
  free(G);
}


path *new_path(graph *G,int *P,int n){
/*
  Créer un chemin d'un graphe G, défini par une liste P de n sommets
  de G. Attention, le pointeur P est utilisé par le chemin (struct
  path) renvoyé. Il ne faut pas détruire P après cet appel. P sera
  libéré par free_path(). Si P=NULL, alors le champs ->P de taille n
  est alloué, et le champs ->n=0. C'est une façon de créer un chemin
  simple et vide. Dans ce cas n représente la longueur maximum d'un
  chemin de G. Le champs ->V (de taille G->n) est initialisé selon P,
  si P<>NULL, ou à -1 sinon.
*/
  if(G==NULL) return NULL;

  NALLOC(path,Q,1);
  ALLOCZ(Q->V,G->n,-1);
  Q->aux=NULL;
  Q->nG=G->n;

  if(P==NULL){
    ALLOC(Q->P,n);
    Q->n=0;
  }
  else{
    int i;
    Q->P=P;
    for(i=0;i<n;i++) Q->V[P[i]]=i;
    Q->n=n;
  }

  return Q;
}


void free_path(path *P){
  if(P==NULL) return;
  free(P->P);
  free(P->V);
  FREE2(P->aux,P->nG);
  free(P);
}


/* structure de données pour une forêt enracinée */
typedef struct{
  int n;       /* nombre de sommets */
  int nroot;   /* nombre de racines, c'est-à-dire d'arbres de la forêt */
  int *lroot;  /* lroot[i]=i-ème racine de la forêt, i=0..nroot-1 */
  int *height; /* height[i]=hauteur du i-ème arbre, i=0..nroot-1 */
  int *root;   /* root[u]=racine de u (root[u]=u si u racine) */
  int *parent; /* parent[u]=parent de u, -1 si u racine */
  int *nchild; /* nchild[u]=nombre de fils de u */
  int **child; /* child[u][i]=i-ème fils de u, i=0..nchild[u]-1 */
  int *dfs;    /* dfs[u]=ordre dfs de u (dfs[u]=i <=> order[i]=u) */
  int *order;  /* order[i]=u si u est le i-ème sommet dans le parcours pré-fixe */
  int *post;   /* post[i]=u si u est le i-ème sommet dans le parcours post-fixe */
  int *weight; /* weight[u]=nombre de descendents de u, u compris */
  int *depth;  /* depth[u]=profondeur à la racine de u */
  int *light;  /* light[u]=plus proche ancêtre (propre) léger de u, =-1 si u racine */ 
  int *apex;   /* apex[u]=premier sommet de la branche lourde de u */
  int *heavy;  /* heavy[u]=1 ssi l'arête entre u et son parent est lourde, =0 si u racine */
} tree;


tree *new_tree(int n){
/*
  Renvoie un objet de type "tree", un arbre enraciné à n sommets. Les
  champs sont initialisés à leurs par défaut. Si n>0, alors les
  tableaux simple de taille n sont alloués, mais pas les doubles
  tableaux (child).
*/

  NALLOC(tree,T,1);
  T->n=max(n,0);
  T->nroot=-1;
  T->lroot=NULL;
  T->height=NULL;
  T->root=NULL;
  T->parent=NULL;
  T->nchild=NULL;
  T->child=NULL;
  T->dfs=NULL;
  T->order=NULL;
  T->post=NULL;
  T->weight=NULL;
  T->depth=NULL;
  T->light=NULL;
  T->apex=NULL;
  T->heavy=NULL;

  if(n>0){
    ALLOC(T->lroot,n);
    ALLOC(T->height,n);
    ALLOC(T->root,n);
    ALLOC(T->parent,n);
    ALLOC(T->nchild,n);
    ALLOC(T->child,n);
    ALLOC(T->dfs,n);
    ALLOC(T->order,n);
    ALLOC(T->post,n);
    ALLOC(T->weight,n);
    ALLOC(T->depth,n);
    ALLOC(T->light,n);
    ALLOC(T->apex,n);
    ALLOC(T->heavy,n);
  }

  return T;
}


void free_tree(tree *T)
/*
  Libère un arbre T et tous ses tableaux.
*/
{
  if(T==NULL) return;

  free(T->lroot);
  free(T->height);
  free(T->root);
  free(T->parent);
  free(T->nchild);
  FREE2(T->child,T->n);
  free(T->dfs);
  free(T->order);
  free(T->post);
  free(T->weight);
  free(T->depth);
  free(T->light);
  free(T->apex);
  free(T->heavy);
  free(T);
}


tree *MakeTree(int *P,int n,int code){
/*

  NON UTILISEE, NON TESTEE

  Construit une forêt enracinée (structure tree) depuis une relation
  de parentée P à n sommets (tableau d'entiers P[u]=père(u) ou -1 s'il
  n'en a pas). Certains champs de la structure sont initialisés
  suivant la valeur binaire de code:

  code&1=0: la relation de parentée P est dupliquée
        =1: le pointeur P est simplement copié

  code&2=0: le champs height est calculé
        =1: le champs height est calculé puis supprimé

  code&4=0: le champs dfs est calculé
        =1: le champs dfs est calculé puis supprimé

  code&8=0: le champs order est calculé
        =1: le champs order est calculé puis supprimé

  code&16=0: le champs post est calculé
         =1: le champs post est calculé puis supprimé

  code&32=0: le champs depth est calculé
         =1: le champs depth est calculé puis supprimé

	 -----

  code&64=0: le champs weight est calculé
         =1: le champs weight n'est pas calculé

  code&128=0: le champs lroot est calculé
          =1: le champs lroot n'est pas calculé

  code&256=0: le champs light est calculé
          =1: le champs light n'est pas calculé

  code&512=0: le champs apex est calculé
          =1: le champs apex n'est pas calculé

  code&1024=0: le champs heavy est calculé
           =1: le champs heavy n'est pas calculé

  Les champs height,dfs,order,post,depth sont calculés pendant le
  parcours (qui n'est pas effectué si tous les bits de ces champs sont
  à 1).  Donc mettre le bit correspondant à 1 ne fait pas gagner de
  temps, mais de l'espace. Les complexités en temps et en espace sont
  en O(n).
*/
  if(n<=0) return NULL; /* problème */

  tree *T=new_tree(n); /* T->n=n>0 */
  int u,p,i;

  /* copie le tableau A ou pas */
  if(code&1) T->parent=P;
  else ALLOCZ(T->parent,n,P[_i]);

  /* calcule le nombre de racines et le nombre de fils pour chaque
     sommets */

  ALLOCZ(T->nchild,n,0);
  for(u=0;u<n;u++){
    p=T->parent[u]; /* p=père(u) */
    if(p<0){
      T->root[u]=u;
      T->nroot++;
    }else T->nchild[u]++;
  }

  /* alloue les tableaux T->child et remet à zéro T->nchild */

  for(u=0;u<n;u++)
    if(T->nchild[u]){
      ALLOC(T->child[u],T->nchild[u]);
      T->nchild[u]=0;
    }
  
  /* remplit les tableaux T->child et rétablit T->nchild */
  /* remplit aussi la liste des racines */

  ALLOC(T->lroot,T->nroot);
  T->nroot=0; /* on va le recalculer */
  for(u=0;u<n;u++){
    p=T->parent[u]; /* p=père(u) */
    if(p<0) T->lroot[T->nroot++]=u; /* ajoute u aux racines */
    else T->child[p][T->nchild[p]++]=p; /* ajoute u aux fils de p */
  }

  /* parcoure la forêt (si nécessaire) */

  if((code&(4+8+16+32+64))==0) goto maketree_suite; /* saut si les bits 2 à 6 sont nuls */

  NALLOC(int,pile,n); /* pile */
  NALLOC(int,next,n); /* next[u]=compteur courant du prochain voisin de u à visiter */
  int sp=-1; /* sp=sommet de la pile, pile[sp]=dernier élément empilé */
  int v,dfs=0; /* dfs=date de première visite */
  p=0; /* ordre post-fixe */

  for(i=0;i<T->nroot;i++){ /* pour chaque racine faire ... */

    u=T->root[i];
    pile[++sp]=u; /* empile la racine i */
    next[u]=0; /* 1er voisin de u à visiter */
    T->dfs[u]=dfs; /* un sommet visité */
    T->order[dfs++]=u; /* ordre pré-fixe des sommets */
    T->depth[u]=0; /* profondeur du sommet u */
    T->height[i]=0; /* hauteur de la i-ème racine */

    while(sp>=0){ /* tant que la pile n'est pas vide */
      u=pile[sp]; /* u=sommet courant sur la pile */
      if(next[u]<T->nchild[u]){ /* on visite le voisin de u d'indice next[u] */
	v=T->child[u][next[u]++]; /* v=voisin de u à empiler */
	pile[++sp]=v; /* on empile le voisin */
	T->dfs[v]=dfs; /* date de première visite de v */
	T->order[dfs++]=v; /* ordre du sommet */
	T->depth[v]=T->depth[u]+1; /* hauteur de v */
	T->height[i]=max(T->height[i],T->depth[u]);
      }else{ /* on a visité tous les voisins de u */
	T->post[p++]=pile[sp--]; /* on dépile u, ordre post-fixe des sommets */
      }
    }
    
  }
  
  free(pile);
  free(next);

  /* calcule le poids des sommets */

  if((!(code&64))||(code&256)||(code&512)||(code&1024)){
    ALLOCZ(T->weight,n,1); /* tout le monde a poids 1 au départ */
    for(i=0;i<n;i++){
      u=T->post[i]; /* parcours post-fixe */
      p=T->parent[u]; /* p=père de u */
      if(p>=0) T->weight[p] += T->weight[u]; /* le père reçoit le poids de son fils */
    }
  }

  /* calcule l'ancêtre léger de chaque sommets (il faut les poids) */

  if(!(code&256)){
    ALLOC(T->light,n);
    for(i=0;i<n;i++){
      u=T->order[i]; /* parcours préfixe */
      p=T->parent[u]; /* p=père de u */
    }
  }

 maketree_suite:

  if(code&2)   { free(T->height); T->height=NULL; }
  if(code&4)   { free(T->dfs);    T->dfs=   NULL; }
  if(code&8)   { free(T->order);  T->order= NULL; }
  if(code&16)  { free(T->post);   T->post=  NULL; }
  if(code&32)  { free(T->depth);  T->depth= NULL; }
  if(code&128) { free(T->lroot);  T->lroot= NULL; }
  if(code&256) { free(T->light);  T->light= NULL; }
  if(code&512) { free(T->apex);   T->apex=  NULL; }
  if(code&1024){ free(T->heavy);  T->heavy= NULL; }

  return T;
}


int fcmp_int(const void *P,const void *Q)
/* Compare deux entiers, pour qsort(). */
{  
  return *(int*)P - *(int*)Q;
}


int fcmp_point(const void *P,const void *Q)
/* Compare deux points, pour qsort(). */
{
  if(((point*)P)->x < ((point*)Q)->x) return -1;
  if(((point*)P)->x > ((point*)Q)->x) return 1;
  if(((point*)P)->y < ((point*)Q)->y) return -1;
  if(((point*)P)->y > ((point*)Q)->y) return 1;
  return 0;
}


int fcmp_profile(const void *P,const void *Q)
/*
  Compare deux profiles, pour qsort(). Les profiles de plus grande
  longueur sont classés avant les plus courts, ceux-ci étant plus
  discriminant.
*/
{
  int* A=*(int**)P;
  int* B=*(int**)Q;

  if(*A>*B) return -1; /* si long(A)>long(B), alors A<B */
  if(*A<*B) return 1; /* si long(A)<long(B), alors A>B */
  /* ici, profiles de même longueur n=A[0] */

  int u=2,n=*A; /* surtout ne pas utiliser A[1] */

  for(;u<n;u++){
    if(A[u]<B[u]) return -1;
    if(A[u]>B[u]) return 1;
  }

  return 0;
}


int fcmp_graphid(const void *P,const void *Q)
/* Compare les identifiants de deux graphes. Sert pour qsort() et
   bsearch(). Ici, P et Q sont des pointeurs de (graph*). */
{
  return (*(graph**)P)->id - (*(graph**)Q)->id;
}


enum{
  R_EQ,   // code =x
  R_INF,  // code <x
  R_SUP,  // code >x
  R_INT,  // code x-y
  R_TRUE  // code t
};

int ReadRange(char *s,int *R)
/*
  Lit une chaîne de caractères décrivant un intervalle de valeurs
  entières, et renvoie dans le tableau d'entiers R les valeurs et les
  codes correspondant pour que la fonction InRange(x,R) puisse
  fonctionner. En quelque sorte cette fonction prépare la fonction
  InRange(). On ignore les caractères non reconnus (pas d'erreur). On
  renvoie le nombre d'opérations décodées, c'est-à-dire le nombre de
  valeurs écrites dans le tableau R, nombre qui est aussi écrit dans
  la première entrée de R.

  Ex: s="1,-3,5-7,>30,<50" (on interprète les "," comme des "ou")
  => R={12,R_EQ,1,R_EQ,-3,R_INT,5,7,R_SUP,30,R_INF,50} (12=taille(R))

  La valeur des codes d'opérations (R_EQ, R_INF, ...) sont données par
  l'enum ci-dessus.
*/
{
  if(R==NULL) return 0;
  if(s==NULL){ R[0]=R_EQ; return 0; }

  int i,r,p,x,start,c;
  i=x=c=0;
  r=start=p=1;

  /* r=indice de R[] */
  /* i=indice de s[] */
  /* x=valeur entière lue */
  /* c=1 ssi le code d'opération a été détecté */
  /* start=1 ssi on va commencer à lire un entier */
  /* p=1 ou -1, signe de x */

  while(s[i]!='\0'){
    if(s[i]=='='){ R[r++]=R_EQ; c=start=p=1; }
    if(s[i]=='<'){ R[r++]=R_INF; c=start=p=1; }
    if(s[i]=='>'){ R[r++]=R_SUP; c=start=p=1; }
    if(s[i]=='-'){
      if(start) p=-p;
      else{ R[r++]=R_INT; R[r++]=x; c=start=p=1; }
    }
    if(s[i]=='t'){ x=R_TRUE; c=r=1; break; } /* t=true, pour avoir false faire "not" et "t" */
    if(s[i]=='p') PVALUE=1; /* pas de code pour "p" */
    if(s[i]==','){
      if(c==0) R[r++]=R_EQ; /* code '=' par défaut */
      R[r++]=x; c=0; start=p=1;
    }
    if(('0'<=s[i])&&(s[i]<='9')){
      if(start) start=0;
      x=x*10+p*(s[i]-'0'); /* on lit x en base 10 en tenant compte du signe p */
    }
    if(start) x=0;
    i++;
  }

  if(PVALUE==i){ x=R_TRUE;c=1; } /* si s="p", alors comme "t" */
  if(c==0) R[r++]=R_EQ;
  R[r++]=x; /* on ajoute le dernier opérande */
  R[0]=r;
  return r;
}


int InRange(int x,int* R)
/*
  Détermine si x appartient aux valeurs décrites par le "range" R.
  R[0] est la taille de R, R[0] compris.
*/
{
  int i,n,t;
  n=R[t=0]; /* n=taille(R) */
  i=1; /* commence à lire R[1] */
  CVALUE=x;

  while(i<n){
    switch(R[i++]){ /* lit le code d'opération */
    case R_EQ  : t=(x==R[i++]); break;
    case R_INF : t=(x<R[i++]); break;
    case R_SUP : t=(x>R[i++]); break;
    case R_INT : t=((R[i]<=x)&&(x<=R[i+1])); i+=2; break;
    case R_TRUE: return 1;
    default: Erreur(16); /* ne devrait jamais se produire */
    }
    if(t) break;
  }
  return t;
}


void PrintGraphList(graph *G)
/*
  Affiche le graphe G sous la forme d'une liste d'adjacence. Tient
  compte de SHIFT.
*/
{
  if(G==NULL){ printf("NULL\n"); return; }
  int u,d,i,n=G->n;

  for(u=0;u<n;u++){
    printf("%i:",u+SHIFT);
    for(i=0,d=G->d[u];i<d;i++){
      printf(" %i",G->L[u][i]+SHIFT);
    }
    printf("\n");
  }
  return;
}


void PrintGraphMatrix(graph *G)
/*
  Affiche le graphe G sous la forme d'une matrice d'adjacence complète
  ou triangulaire supérieure (en tennant compte du FORMAT, smatrix ou
  matrix). La complexité en espace est seulement de O(n).
*/
{
  int u,d,i,z,t,n=G->n;

  NALLOCZ(int,M,n,0);
  t=(FORMAT==F_smatrix);

  for(u=z=0;u<n;u++){
    if(t) z=u;
    for(i=0,d=G->d[u];i<d;M[G->L[u][i++]]=1);
    for(i=0;i<n;i++)
      if(i<z) printf(" ");
      else printf("%c",'0'+M[i]);
    for(i=0;i<d;M[G->L[u][i++]]=0); /* remet rapidement M[] tout à 0 */
    printf("\n");
  }

  free(M);
  return;
}


void PrintPath(graph *G,path *P)
/*
  Affiche le chemin P d'un graphe G.
  Sert pour le débugage.
*/
{
  if((G==NULL)||(P==NULL))
    printf("NULL\n");
  else{
    int i,j,u,d;
    for(i=0;i<P->n;i++)
      if(P->V[P->P[i]]!=i) break;
    if(i<P->n) goto error;
    for(u=0;u<G->n;u++)
      if((P->V[u]>=0)&&(P->P[P->V[u]]!=u)) break;
    if(u<G->n) goto error;
    printf("P->aux:");
    if(P->aux==NULL) printf(" NULL\n");
    else{
      printf("\n");
      for(i=0;i<P->n;i++){
	u=P->P[i];
	d=P->aux[u][0];
	printf("  %i:",u);
	for(j=1;j<=d;j++){
	  printf(" %i",P->aux[u][j]);
	}
	printf("\n");
      }
    }
  }
  return;
  
 error:
  printf("Chemin incohérent.\n");
  return;
}


int *SortGraph(graph *G,int code)
/*
  Force le tri des listes d'adjacence d'un graphe G, c'est-à-dire pour
  chaque sommet u, G->L[u] est une liste d'entiers triés par ordre
  croissant. Le champs G->sort est mis à jour.  L'algorithme effectue
  un simple appel à qsort(). Sa complexité est à peu près en
  O(n+m*log(m/n)).

  Si code=0, on s'arrête après l'étape du tri. Sinon (code<>1) on
  lance une étape de vérification du graphe: présence de multi-arêtes,
  de boucles, etc. Le temps de la vérification est comparable à celui
  du tri. Le résultat de la vérification est un tableau de
  statistiques S de taille fixe (déclaré en static qui ne doit pas
  être libéré par l'appelant), ayant la signification suivante:

    S[0]=nombre de boucles
    S[1]=nombre de multi-arcs
    S[2]=nombre d'arcs (avec multi-arcs et boucles)
    S[3]=nombre d'adjacence non-symétriques
    S[4]=nombre de voisins d'ID < 0
    S[5]=nombre de voisins d'ID >= n
    S[6]=1 ssi G est simple et non-orienté
    S[7]=degré maximum
    S[8]=degré minimum
    S[9]=nombre de sommets isolés

  A l'issue de la vérification, G->sym est mise à jour.
*/
{
  if(G==NULL) return NULL;
  int u,n=G->n;

  /* trie G */
  for(u=0;u<n;u++)
    qsort(G->L[u],G->d[u],sizeof(int),fcmp_int);

  G->sort=1;
  if((code==0)||(G->n==0)) return NULL;
  
  /* satistiques sur G */
  static int S[10]; /* static est important, car on fait return S */
  int v,i,d,w;
  
  memset(S,0,sizeof(S)); /* initialise les stats à 0 */
  S[7]=S[8]=G->d[0]; /* il faut G->d<>NULL */
  
  for(u=0;u<n;u++){ /* parcoure le graphe */
    d=G->d[u]; S[7]=max(S[7],d); S[8]=min(S[8],d);
    S[9] += (d==0); /* un sommet isolé */
    S[2]+=d; /* ajoute le nombre de voisins */
    w=-1; /* w=voisin précédant le voisin courant v */
    for(i=0;i<d;i++){ /* pour chaque voisin */
      v=G->L[u][i];
      if(u==v) S[0]++; /* une boucle */
      if(v==w) S[1]++; /* une multi-arête */
      w=v; /* mémorise le dernier voisin rencontré */
      if(v<0){ S[4]++; continue; } /* un voisin négatif */
      if(v>=n){ S[5]++; continue; } /* un voisin trop grand */
      if(bsearch(&u,G->L[v],G->d[v],sizeof(int),fcmp_int)==NULL) S[3]++; /* un arc asymétrique */
    }
  }

  S[6]=((S[0]+S[1]+S[3]+S[4]+S[5])==0); /* vrai ssi G simple et non-orienté */
  G->sym=(S[3]==0);
  return S;
}


void GraphRealloc(graph *G,int *D){
/*
  Redimensionne le graphe G à G->n sommets suivant le tableau de degré
  D. On réajuste en premier les tableaux G->d et G->L pour qu'ils
  aient une taille G->n, puis on réajuste les listes d'adjacences des
  sommets de G suivant le tableau des degrés D (qui doit être de
  taille au moins G->n). Si D[u] est plus petit que G->d[u], alors la
  liste G->L[u] est tronquée. Si D[u] est plus grand que G->d[u],
  alors G->L[u] est réajusté. Le degré G->d[u] est initialisé au
  minimum G->d[u] et de D[u]. NB: le nombre d'arêtes G->m, qui a pu
  changer, est réinitialisé à -1. G->sort n'est pas changé car l'ordre
  des listes G->L n'est pas modifié.

  Pour plonger G dans un graphe complet faire:
    NALLOCZ(int,D,n,n-1);
    GraphRealloc(G,D);
    free(D);
*/
  int u,d,n=G->n;
  for(u=0;u<n;u++){
    d=D[u];
    REALLOC(G->L[u],d);
    G->d[u]=min(G->d[u],d);
  }

  /* Il ne faut pas réajuster G->d et G->L avant la boucle for(u=...)
     car potentiellement on libère G->d et G->L. Or il est possible
     d'avoir D=G->d. */

  REALLOC(G->d,n);
  REALLOC(G->L,n);
  G->m=-1; /* le nombre d'arêtes n'est plus à jour */
  return;
}


graph *new_fullgraph(int n){
/*
  Retourne un graphe G comme new_graph(n), mais en plus alloue G->L[u]
  de taille max(n-1,1), et initialise G->d[u]=0 pour tous les sommets
  u. Une fois le graphe construit, on peut rédimensionner le graphe
  grâce à GraphRealloc, comme dans l'exemple:

    graph *G=new_fullgraph(n);
      ...
      ADD_EDGE(G,u1,v1);
      ADD_EDGE(G,u2,v2);
      ...
    GraphRealloc(G,G->d);
    ...
    free_graph(G);
*/

  if(n<1) return NULL;
  graph *G=new_graph(n);
  int u,n1=max(n-1,1);

  for(u=0;u<n;u++){
    G->d[u]=0;
    ALLOC(G->L[u],n1);
  }
  
  return G;
}


graph *ExtractSubgraph(graph *G,int *T,int k,int code){
/*
  Construit, à partir d'un graphe G et d'une liste T de k sommets, un
  nouveau graphe S correspondant au sous-graphe de G induit par les
  sommets de T (code=1) ou de V(G)\T (si code=0). Les sommets de S
  sont dans [0,k[ (ou [0,n-k[ si code=0).

  On peut ainsi faire une copie C du graphe G simplement en faisant:
  graph *C=ExtractSubgraph(G,NULL,0,0);

  Effet de bord: S->pint1 est alloué si T<>NULL. Dans ce cas on
  renvoie dans S->pint1 un tableau X de taille G->n indiquant la
  renumérotation de G: pour tout sommet u de G (u dans [0,G->n[)
  S->pint1[u]=0 si u est un sommet abscent de S et S->pint1[u]=d>0 si
  u est numéroté d-1>=0 dans S. Le nombre d'arêtes S->m du graphe S
  renvoyé est à jour. L'ordre relatif des listes de G est préservé. En
  particulier, si G->sort=1, alors le sous-graphe sera aussi
  trié. G->sym est aussi copié.
*/
  
  if(G==NULL) return NULL;
  int n=G->n;
  int u,v,d,i,s,ns,m;
  graph *S;

  NALLOC(int,X,n);
  for(u=1-code,i=0;i<n;i++) X[i]=u;
  if(T!=NULL) for(i=0;i<k;i++) X[T[i]] ^= 1;
  for(i=d=0;i<n;i++) if(X[i]) X[i]=++d; 
  /* ici X[i]=0 si i est un sommet à supprimer */
  /* ici X[i]=d>0 si i doit être renuméroté en d-1>=0 */

  ns=(code)?k:n-k;
  S=new_fullgraph(ns);

  for(s=u=m=0;u<n;u++)
    if(X[u]){ /* si u existe, s=X[u]-1 */
      d=G->d[u];
      for(i=0;i<d;i++){
	v=G->L[u][i];
	if(X[v]){ m++; S->L[s][S->d[s]++]=X[v]-1; } /* si v existe */
      }
      s++;
    }

  /* réduit la taille des listes */
  GraphRealloc(S,S->d);

  S->pint1=X;
  S->sort=G->sort;
  S->sym=G->sym;
  S->m=(m>>1);
  return S;
}


graph *List2Graph(list *L,int code){
/*
  Retourne un graphe G simple à partir d'un graphe défini par une
  liste L de codes (voir File2List() pour le codage précis du type
  "list"). Certaines opérations sont effectuées sur L en fonction de
  la valeur binaire de code:

  - code&1 =1: optimisation des listes du graphe (tri par ordre croissant)
           =0: sans optimisation
  - code&2 =1: auto-détection du shift dans L (pour "load file")
           =0: pas d'auto-détection du shift
  - code&4 =1: gestion d'un sous-graphe (V,NF) => code&2=0
           =0: pas de sous-graphe

  Les codes suivants servent à List2Family():

  - code&8 =1: tri de la famille suivant les identifiants (sert pour List2Family)
           =0: pas de tri de la famille
  - code&16=1: ne libère pas la liste L (sert pour List2Family)
           =0: libère la liste L
  - code&32=1: renvoie toujours un graphe, le 1er si c'est une famille
           =0: renvoie une famille si c'est une famille

  Pour calculer le graphe (et sa liste d'adjacence) on effectue
  plusieurs passes sur L: une passe pour déterminer n; une autre pour
  calculer les degrés des sommets; et une 3e pour remplir G et pour
  éventuellement libérer la liste L.
*/
  if(L==NULL) return NULL; /* si pas de cellule, ne rien faire */

  int u,v,x,n,*D;
  graph *G;
  list *p;

  u=INT_MAX;
  if(code&4){ /* si sous-graphe définit par (V,NF) */
    p=L;
    while(p){
      p->item=V[p->item]-SHIFT;
      p=p->next;
    }
    n=NF; /* on connaît n */
  }
  else{ /* sinon, on calcule n, et on lit les valeurs min (=u) et
	   valeur max (=v) de L */
    p=L; v=0;
    while(p){
      x=p->item;
      if(x<u) u=x;
      if(x>v) v=x;
      p=p->next;
    }
    if(code&2){ /* on décale les valeurs dans [0,n[ */
      p=L;
      while(p){
	p->item -= u;
	p=p->next;
      }
      n=v+1-u;
    }else{
      if((u<0)||(v<0)) Erreur(22); /* il ne devrait pas avoir ici de valeur < 0 */
      n=v+1;
    }
  }

  ALLOCZ(D,n,0);

  /* on lit les degrés (sortant) des sommets, et les met dans le
     tableau D. NB: la variable u n'est pas initialisé, car on passe
     toujours d'abord par un item de type T_NODE */
  
  p=L; x=1; /* x=1 ssi on n'est PAS dans un groupe */
  while(p){
    v=p->item;
    if(p->type==T_NODE) x=1;
    else if(p->type==T_EDGE) { D[u]++; D[v]++; }      /* u-v */
    else if(p->type==T_OPENE){ D[u]++; D[v]++; x=0; } /* u-(v */
    else if(p->type==T_ARC)    D[u]++;                /* u->v */
    else if(p->type==T_OPENA){ D[u]++; x=0; }         /* u->(v */
    if(x) u=v;
    p=p->next;
  }

  /* initialise la liste d'adjacence G. On se sert plus tard de D[u]
     pour indiquer la prochaine position libre dans G[u][]. */

  G=new_graph(n); /* G->n=n, alloue G->d et G->L */
  for(u=0;u<n;u++){
    ALLOC(G->L[u],D[u]); /* alloue une liste pour chaque sommet */
    G->d[u]=D[u]; /* G->d[u]=deg(u) */
    D[u]=0; /* prochaine position libre dans G[u] */
  }

  /* Remplit G. On met aussi à jour G->sym (orienté ou pas). On
     pourrait tester à la volée si les listes sont triées et mettre à
     jour G->sort. */
  
  p=L; x=1; /* x=1 ssi on n'est PAS dans un groupe */
  while(p){
    v=p->item;
    if(p->type==T_NODE) x=1;
    else if(p->type==T_EDGE) { G->L[u][D[u]++]=v; G->L[v][D[v]++]=u; }      /* u-v */
    else if(p->type==T_OPENE){ G->L[u][D[u]++]=v; G->L[v][D[v]++]=u; x=0; } /* u-(v */
    else if(p->type==T_ARC)  { G->L[u][D[u]++]=v; G->sym=0; }               /* u->v */
    else if(p->type==T_OPENA){ G->L[u][D[u]++]=v; G->sym=0; x=0; }          /* u->(v */
    if(x) u=v;
    p=p->next;
  }

  /* libère L si bit-4 à 1 */
  if(!(code&16)){
    p=L;
    while(p){
      L=p;
      p=p->next;
      free(L);
    }
  }
  
  free(D); /* plus besoin de D */
  if(code&01) SortGraph(G,0);
  return G;
}


graph *List2Family(list *L,int code){
/*
  Transforme une liste en famille de graphes.  Si L représente un
  graphe simple (pas de type T_NB ou T_ID), alors un graphe simple est
  retournée (plutôt qu'une famille à un seul élément). Donc,
  List2Family() généralise List2Graph(). On utilise List2Graph() comme
  sous-routine. Pour "code" voir List2Graph().
    
  Effet de bord:
  - la famille est triée par ID croissant si code&8=1
  - la liste L est libérée si code&16=0
  - on retourne un graphe si code&32=1 (plutôt qu'une famille)
*/
  if(L==NULL) return NULL; /* liste vide */
  if(L->type!=T_NB) return List2Graph(L,code); /* si graphe */

  /* ici on a donc une famille */
  int n=L->item; /* nb de graphes dans la liste */
  list *T;
  
  if(n<=0){ /* famille vide */
    if(code&16) /* libère éventuellement L */
      while(L){
	T=L->next; /* ici L<>NULL */
	free(L);
	L=T;
      }
    return NULL;
  }

  int i,id;
  graph *F=new_graph(0);
  list *P;

  F->f=n;
  ALLOC(F->G,n); /* F->G[0..n[: tableau de n pointeurs de graphes */
  T=L; L=L->next;
  if(!(code&16)) free(T); /* on libère l'élément (n,T_NB) */

  /* ici L=début du 1er graphe de la famille */
  for(i=0;i<n;i++){ /* pour chaque graphe */
    /* ici L pointe sur un élément (id,T_ID) */
    if((L==NULL)||(L->type!=T_ID)) Erreur(13);
    id=L->item; /* identifiant du graph i */
    T=L->next;
    if((code&16)==0) free(L); /* on libère l'élément (id,T_ID) */
    P=L=T; /* P=L=T=tête courante du graphe i */
    while((L)&&(L->type!=T_ID)){ P=L; L=L->next; } /* cherche la fin du graphe i */
    /* ici le graphe i va de T à P */
    P->next=NULL; /* on coupe la liste */
    F->G[i]=List2Graph(T,code); /* T=liste du graphe i */
    F->G[i]->id=id; /* Attention! F->G[i] n'existe qu'après l'appel à List2Graph() */
    if(code&16) P->next=L; /* recolle la liste si on ne souhaite pas la libérer */
    }

  /* éventuellement trie la famille suivant les IDs */
  if(code&8) qsort(F->G,F->f,sizeof(graph*),fcmp_graphid);

  /* extrait le premier graphe */
  if(code&32){
    graph *G=ExtractSubgraph(F->G[0],NULL,0,0); /* copie le premier graphe */
    free_graph(F); /* libère complètement la famille F */
    F=G; /* F=premier graphe */
  }

  return F;
}

list *File2List(char *file){
/*
  Lit le fichier "file" contenant un graphe (ou une famille) au format
  standard, orientés ou non, et retourne le contenu dans une liste.
  Tient compte de -shift mais pas des noms originaux (-label 1). Dans
  le cas d'une famille de graphes, il est possible de spécifier un
  "range" pour "file" avec la forme: "file:range" où "range" est une
  liste de valeurs ayant la même signification que pour "-filter F id
  value". Par exemple, "file:5" spécifie le graphe d'identifiant 5, et
  "file:5-8" est la famille contenant les graphes d'identifiant
  5,6,7,8. Notez que "-:5" est le graphe d'identifiant 5 de la famille
  lue depuis l'entrée standard.

  Chaque élément de la liste est une paire d'entiers (item,type) où
  "type" précise le rôle joué par l'entier "item". Voir l'enum pour
  une description des types.

  Si, par exemple, le graphe est "0-1 2 1-2-3 5->4" alors la liste
  retournée sera { (0,T_NODE), (1,T_EDGE), (2,T_NODE), (1,T_NODE),
  (2,T_EDGE), (3,T_EDGE), (5,T_NODE), (4,T_ARC) }.

  Si le graphe est "0-(2 3) 4 5->(6 7)" alors la liste retournée sera
  { (0,T_NODE), (2,T_OPENE), (3,T_OPENE), (4,T_NODE), (5,T_NODE),
  (6,T_OPENA), (7,T_ARC) }. Autrement dit, T_OPENE ou T_OPENA décrive
  un groupe d'arêtes ou d'arcs.

  NB: "i-j-(k ...)" est correct mais pas "i-(j ...)-k".

  La fonction est généralisée à la lecture d'une famille de graphes.
  Si le fichier contient "[5] 0-1 [8] 0->2-1" alors la liste
  contiendra { (2,T_NB), (5,T_ID), (0,T_NODE), (1,T_EDGE), (8,T_ID),
  (0,T_NODE), (2,T_ARC), (1,T_EDGE) }, où le premier élément (n,T_NB)
  signifie qu'il s'agit d'une famille de n graphes, et où (u,T_ID)
  signifie que u est l'identifiant du graphe à venir.
*/
  FILE *f;
  list *T; /* tête de la liste */
  list *L; /* élément courant */
  list *P; /* sauvegarde le dernier élément (sentinelle qui faudra supprimer) */
  int read=1; /* pour InRange(), par défaut on lit tout */
  char *r=NULL,c[2],*s;
  int range[CMDMAX]={2,R_TRUE}; /* par défaut: range toujours vrai */
  unsigned v; /* valeur lue */
  long p; /* position dans le fichier f */
  int n=0; /* nb de graphes dans la famille */
  int t=-1; /* t<0 si on est pas dans un groupe */

  T=P=L=new_list(); /* crée la liste */

  /* TO DO: si file commence par " ", alors file représente le graphe
     lui-même (pour future option -add/-del). Par exemple, file="
     5-6,7-8-0". Dans ce cas on écrit un fichier temporaire avec "5-6
     7-8-0" et on continue normalement. On détruit ensuite ce
     fichier. NB: le contenu de file est modifié. */
  /*
    char *s;
    if(file[0]==' '){
    for(v=0; v<strlen(file); v++) if(file[v]==',') file[v]=' ';
    s=strdup("/tmp");
    f=fopen(s,"rw");
    fputs(f,file);
    rewind(f); // pas la peine de fermer le fichier
    file=s;
  }
  */
  
  /* ouverture du fichier: file ou file:range */

  f=strcmp(file,"-")? fopen(file,"r"):stdin;
  if(f==NULL){ /* on a pas réussit à ouvrir file */
    fclose(f); /* il faut fermer le fichier, même si c'est NULL ! */
    if((r=strchr(file,':'))==NULL) Erreur(7); /* est-ce file:range ? */
    *r='\0'; /* coupe file en (préfixe,r=range) */
    f=strcmp(file,"-")? fopen(file,"r"):stdin;
    if(f==NULL) Erreur(7); /* on a pas réussit à ouvrir file */
    *r++ = ':'; /* rétablit le nom de fichier original de file */
    ReadRange(r,range); /* lecture du range */
  }

  /* lecture du fichier */

  while(!feof(f)){
    p=ftell(f); /* ftell() vaut toujours -1 si f=stdin (non seekable) */
    fscanf(f,"//%*[^\n]\n"); /* essaye de lire "//" */
    if(ftell(f)>=p+2) continue; /* on a lu au moins 2 caractères -> commentaire */
    fseek(f,p,SEEK_SET);
    if(fscanf(f," [%u]",&v)>0){
      read=InRange(v,range);
      if(read){ L=Insert(P=L,v,T_ID); n++; }
      continue;
    }
    if(read){
      fseek(f,p,SEEK_SET);
      if(fscanf(f,"-%u",&v)>0){ L=Insert(P=L,v,T_EDGE); continue; }
      fseek(f,p,SEEK_SET);
      if(fscanf(f,">%u",&v)>0){ L=Insert(P=L,v,T_ARC); continue; }
      fseek(f,p,SEEK_SET);
      if(fscanf(f,"-( %u",&v)>0){ L=Insert(P=L,v,T_OPENE);
	p=ftell(f);
	if(fscanf(f," %1[)]c",c)>0) t=-1;
	else{ t=T_EDGE; fseek(f,p,SEEK_SET); }
	continue;
      }
      fseek(f,p,SEEK_SET);
      if(fscanf(f,">( %u",&v)>0){ L=Insert(P=L,v,T_OPENA);
	p=ftell(f);
	if(fscanf(f," %1[)]c",c)>0) t=-1;
	else{ t=T_ARC; fseek(f,p,SEEK_SET); }
	continue;
      }
      fseek(f,p,SEEK_SET);
      if(fscanf(f," %u",&v)>0){
	p=ftell(f);
	if(fscanf(f," %1[)]c",c)>0){ L=Insert(P=L,v,t); t=-1; }
	else{ L=Insert(P=L,v,(t<0)?T_NODE:t); fseek(f,p,SEEK_SET); }
	continue;
      }
      
      /* ici on a rien trouvé: est-ce une erreur de format ? */
      fseek(f,p,SEEK_SET);
      s=fgets(c,2,f); /* lit au plus un caractère */
      if((s==NULL)||(c[0]==' ')||(c[0]=='\n')) continue; /* ok si ' ' ou '\n' */
      Erreur(28); /* mauvais format sinon */
    }

    /* ici on a rien trouvé, mais read est faux */
    fseek(f,p,SEEK_SET); /* on a rien trouvé */
    fscanf(f," %*c"); /* lit au moins un caractère */
  }
  
  fclose(f); /* on ferme le fichier */
  free(L); /* supprime le dernier élément (la sentinelle) */
  if(L==T) return NULL; /* si on a lu aucun élément */
  P->next=NULL;

  if(n>0){ /* il s'agit d'une famille */
    /* on ajoute un nouvel élément en tête de la liste indiquant le
       nombre de graphes de la famille */
    L=new_list();
    L->item=n; /* nombre de graphes de la famille */
    L->type=T_NB;
    L->next=T;
    T=L; /* nouvelle tête de liste */
  }

  return T; /* on retourne la tête */
}


graph *File2Graph(char *file,int code){
/*
  Renvoie un graphe (ou une famille) à partir d'un fichier. Pour
  "code" voir List2Graph() & List2Family(). La liste intermédiaire
  calculée par File2List() est toujours libérée.
*/
  graph *G=List2Family(File2List(file),code&(63-16)); /* annule le bit-4 */
  if(G==NULL) Erreur(15);
  return G;
}


double PosAspect(void){
/*
  Donne le coefficient par lequel les positions XPOS/YPOS vont être
  multipliées pour le format dot pour permettre une taille de sommets
  raisonable par rapport à la longueur des arêtes. On tient compte de
  N et de BOXX et BOXY.
*/
  double w=C32*sqrt((double)N); /* la largeur est en sqrt(N) */
  if((BOXX>0.0)&&(BOXY>0.0)) w /= min(BOXX,BOXY);
  if(LABEL>0) w *= 3.0; /* augmente l'aspect si besoin des LABELs (et POS) */
  return w;
}

void BoundingBox(){
/*
  Calcule XMIN,YMIN,XMAX,YMAX des tableaux XPOS/YPOS.  Il faut que
  N>0.
*/
  int i;
  XMIN=XMAX=XPOS[0];
  YMIN=YMAX=YPOS[0];
  for(i=1;i<N;i++){
    XMIN=fmin(XMIN,XPOS[i]);
    YMIN=fmin(YMIN,YPOS[i]);
    XMAX=fmax(XMAX,XPOS[i]);
    YMAX=fmax(YMAX,YPOS[i]);
  }
}


void InitXY(void){
/*
  Initialise les tableaux XPOS et YPOS suivants les options
  éventuelles de -xy. Les variables suivantes peuvent être mise à
  jour: N, XMIN,XMAX,YMIN,YMAX (BoundingBox), VSIZESTD, VSIZEXY.
*/

  int i,k;
  double sx,sy,tx,ty;

  for(;;){ /* pour pouvoir faire un break; */

    if(XYtype==XY_USER) break; /* coordonnées définies par l'utilisateur */

    if(XYtype==XY_FILE){ /* charge à partir d'un fichier et met à jour N */
      N=LoadXY(FILEXY);
      break;
    }

    if(N<0) N=0;
    ALLOC(XPOS,N);
    ALLOC(YPOS,N);

    if(XYtype==XY_UNIF){ /* uniforme dans [0,1[ */
      for(i=0;i<N;i++){
	XPOS[i]=RAND01;
	YPOS[i]=RAND01;
      }
      break;
    }
    
    if(XYtype==XY_PLAW){ /* loi puissance autour des graines choisies dans [0,1[ */
      ALLOC(XSEED,SEEDk);
      ALLOC(YSEED,SEEDk);
      sx=sy=0.0; /* calcule (sx,sy), le barycentre des SEEDk graines */
      for(i=0;i<SEEDk;i++){
	XSEED[i]=RAND01;sx+=XSEED[i];
	YSEED[i]=RAND01;sy+=YSEED[i];
      }
      sx /= (double)SEEDk; sy /= (double)SEEDk;
      sx -= 0.5; sy -= 0.5;
      for(i=0;i<SEEDk;i++){ /* centre par rapport au barycentre */
	XSEED[i] -= sx; /* enlève le barycentre puis décale de 0.5 */
	YSEED[i] -= sy;
      }
      /* on génère les points autour des graines */
      tx=sqrt(log((double)(SEEDk+1.0))/(double)(SEEDk+1.0)); /* rayon r */
      for(i=0;i<N;i++){
	k=random()%SEEDk;    /* on choisit la graine numéro k au hasard */
	sx=2.0*M_PI*RAND01;  /* angle aléatoire */
	sy=tx*pow(RAND01,SEEDp); /* longueur aléatoire */
	XPOS[i]=XSEED[k]+sy*cos(sx);
	YPOS[i]=YSEED[k]+sy*sin(sx);
      }
      break;
    }

    if(XYtype==XY_PERM){ /* permutation de [0,N[ */
      NALLOC(int,P,N);
      for(i=0;i<N;i++) XPOS[i]=(double)(P[i]=i);
      Permute(P,N);
      for(i=0;i<N;i++) YPOS[i]=(double)P[i];
      free(P);
      break;
    }
    
    if(XYtype==XY_MESH){ /* mesh de paramètre Xmesh x Ymesh*/
      for(i=0;i<N;i++) XPOS[i]=(double)((int)(i%Xmesh));
      for(i=0;i<N;i++) YPOS[i]=(double)((int)(i/Ymesh));
      break;
    }
  }

  /* ici on devrait normalement avoir XPOS,YPOS <> NULL */
  if((XPOS==NULL)||(YPOS==NULL)) Erreur(8);

  if(NOISEr>0.0) /* "noise" doit être avant "box" */
    if((XYtype!=XY_PERM)&&(XYtype!=XY_MESH))
      for(i=0;i<N;i++){
	sx=2.0*M_PI*RAND01; /* angle aléatoire */
	sy=NOISEr*pow(RAND01,NOISEp); /* longueur aléatoire */
	XPOS[i] += sy*cos(sx); /* décale XPOS */
	YPOS[i] += sy*sin(sx); /* décale YPOS */
      }

  if((BOXX>0.0)&&(BOXY>0.0)){ /* "box" doit être après "noise" */
    BoundingBox(); /* calcule les BB */
    tx=(XMAX-XMIN); if(tx==0.0) tx=0.5;
    ty=(YMAX-YMIN); if(ty==0.0) ty=0.5;
    sx=2.0*sqrt(N+1); /* le +1 est pour être sûr d'avoir sx>0 */
    tx /= sx; /* tx=largeur de la bande vide */
    ty /= sx; /* ty=hauteur de la bande vide */
    sx=BOXX/(XMAX-XMIN+2.0*tx); /* sx ne peut pas être nul */
    sy=BOXY/(YMAX-YMIN+2.0*ty); /* sy ne peut pas être nul */

    for(i=0;i<N;i++){
      XPOS[i]=sx*(XPOS[i]-XMIN+tx);
      YPOS[i]=sy*(YPOS[i]-YMIN+ty);
    }
  }

  if(ROUND<DBL_DIG){ /* arrondit éventuellement les coordonnées */
    sx=pow(10.0,ROUND);
    for(i=0;i<N;i++){
      XPOS[i]=rint(XPOS[i]*sx)/sx;
      YPOS[i]=rint(YPOS[i]*sx)/sx;
    }
  }

  if(XYunique){ /* élimine les doubles, en triant les points */
    NALLOC(point,P,N);
    for(i=0;i<N;i++) P[i].x=XPOS[i],P[i].y=YPOS[i];
    qsort(P,N,sizeof(point),fcmp_point); /* tri les points */
    point p=P[0]; p.x -= 1.0; /* ici le point p <> du 1er élément */
    for(i=k=0;i<N;i++)
      if(fcmp_point(P+i,&p)){ /* copie que si différent de l'élément p */
	p=P[i];
	XPOS[k]=p.x;
	YPOS[k++]=p.y;
      }
    free(P);
    if(k<N){
      N=k;
      REALLOC(XPOS,N);
      REALLOC(YPOS,N);
    }
  }

  /* on calcule les BB */
  BoundingBox();

  /* mise à jour de la taille des sommets */
  VSIZESTD *= XYvsize;
  VSIZEXY  *= XYvsize;

  return;
}


color *GradColor(color *T,int n,int m){
/*
  Retourne un tableau de m couleurs formant un dégradé obtenu à partir
  d'un tableau T de n couleurs. Pour avoir un dégradé simple d'une
  couleur T[0] à T[1] il faut initialiser T[0],T[1] et poser n=2. Pour
  avoir un dégradé cyclique, il suffit de répéter la couleur T[0] en
  dernière position de T (et ajouter 1 à n, donc d'avoir
  T[n-1]=T[0]). Il faut dans tous les cas n>1 et m>0. On peut avoir
  m<n. Dans ce cas on prend la première couleur de T, puis la i-ème
  couleur est (i*(n-1))/(m-1).
*/
  color c1,c2;
  int i,j,k,r,q,n1,m1,dr,dg,db;

  if(T==NULL) return NULL; /* normalement ne sert à rien */

  NALLOC(color,P,m);
  c2=P[0]=T[0];
  if(m==1) return P;
  /* maintenant m >= 2 */

  m1=m-1; n1=n-1; /* valeurs utilisées souvent */

  if(m<=n){ /* cas où il y a moins de couleurs demandées que dans T */
    for(i=1;i<m;i++) /* m-1 fois */
      P[i]=T[(i*n1+m1-1)/m1]; /* le "+m-2" est pour arrondir à l'entier sup. */
    return P;
  }

  /*
    Cas m>n.  Soient B_1,B_2,...B_(n-1) les n-1>0 blocs de couleurs,
    B_i commençant juste après la couleurs T[i-1] et se terminant avec
    la couleur T[i]. On doit répartir m-1 couleurs dans ces n-1 blocs,
    la couleurs T[0] étant déjà dans P. On met alors
    floor((m-1)/(n-1)) couleurs par blocs, et une de plus pour B_i si
    i<=(m-1)%(n-1).
   */
  r=m1%n1; /* il reste r couleurs */
  q=(m1/n1)+(r>0); /* nombre de couleurs par blocs (+1 au départ si r>0) */
  for(i=j=k=1;i<n;i++){ /* on traite le bloc B_i, P[k]=prochaine couleur libre */
    c1=c2;   /* c1=T[i-1] */
    c2=T[i]; /* c2=T[i] */
    dr=c2.r-c1.r;
    dg=c2.g-c1.g;
    db=c2.b-c1.b;
    for(j=1;j<=q;j++){ /* extrapolation linéaire de q points dans ]c1,c2] */
      P[k].r=c1.r+(j*dr)/q;
      P[k].g=c1.g+(j*dg)/q;
      P[k].b=c1.b+(j*db)/q;
      k++;
    }
    if(i==r) q--; /* une couleur de moins pour B_{r+1}...B_{n-1} */
  }
  return P;
}


int graphical(int *S,int k){
/*
  Vérifie si la suite (n_1,d_1,n_2,d_2,...,n_k,d_k) est graphique ou
  pas, c'est-à-dire s'il existe au moins un graphe simple ayant
  exactement n_i sommets de degré d_i. On renvoie une valeur <0 si la
  séquence n'est pas graphique, et sinon renvoit le nombre n de
  sommets du graphe (n=sum_i n_i). Il faut que n_i,d_i >=0.

  L'algorithme est basé sur le test d'Erdős and Gallai (1960) qui ont
  prouvé qu'une séquence de degrés (t_1,...,t_n) est graphique ssi la
  somme des degrés des sommets est paire et la suite vérifie la
  propriété suivante (ici t_i est le degré du sommet i=1..n):

  (1)   sum_{i=1}^r t_i <= r*(r-1) + sum_{i=r+1}^n min{r,t_i}

  pour chaque entier r=1..n-1 (Skiena 1990, p. 157). Cette propriété
  se généralise aux graphes orientés.  Tripathi et Vijay (2003) on
  montré qu'il suffit de le vérifier pour les valeurs de r où les t_i
  changent.

  Si on note r_j=n[0]+...+n[j], il faut donc vérifier (1) seulement
  pour r=r_0, r=r_1,..., r=r_j, ..., r=r_{k-1}. A l'étape j=0...k-1,
  l'équation (1) se réécrit donc en:
  
  (2) sum_{i=0}^j n[j]*d[j] <= r_j*(r_j-1) + sum_{i=0}^{k-1} n[i]*min{r_j,d[i]}

*/

  if((S==NULL)||(k<=0)) return -1;
  
  /* range les valeurs de S dans n[i] et d[i] */
  NALLOCZ(int,n,k,S[2*_i]);   /* n[0] ... n[k-1]: les n_i */
  NALLOCZ(int,d,k,S[2*_i+1]); /* d[0] ... d[k-1]: les d_i */

  int i;
  int j=0; /* j=somme des degrés des sommets */
  int nb=0; /* nb=nombre de sommets=somme des n[i] */

  for(i=0;i<k;i++){
    nb += n[i];
    j += n[i]*d[i];
  }
  
  if(j&01) nb=-1; /* somme de degré impaire => pas graphique */
  else{ /* somme des degrés paire */
    int rj=0; /* rj=sum_{i=0}^j n[i] */
    int s1=0; /* s1=sum_{i=0}^{j-1} n_i*d_i */
    int s2;   /* s2=sum_{i=j+1}^{k-1} n[i]*min{r_j,d[i]} */
    
    for(j=0;j<k;j++){ /* étape j=0..k-1 */
      s1 += n[j]*d[j];
      rj += n[j];
      for(i=j+1,s2=0;i<k;i++) s2 += n[i]*min(n[j],d[i]);
      if(s1>rj*(rj-1)+s2){ nb=-1; break; } /* équation fausse */
    }
  }
  
  free(n);
  free(d);
  return nb;
}


/***********************************

         ROUTINES POUR LES
          FONCTIONS adj()

***********************************/


double Norme(int i,int j){
/*
  Calcule la distance entre les points (XPOS[i],YPOS[i]) et
  (XPOS[j],YPOS[j]) selon la norme définie par la constante NORM
  (=1,2,...).  Pour L_2 (NORM=2) c'est le carré de la norme qui est
  retourné pour éviter l'extraction d'une racine carrée.

  NORM=1 -> L_1
  NORM=2 -> L_2
  NORM=3 -> L_max
  NORM=4 -> L_min
*/
  double x=XPOS[i]-XPOS[j];
  double y=YPOS[i]-YPOS[j];

  if(NORM==2) return (x*x + y*y); /* norme L_2 */
  
  if(x<0.0) x=-x; /* x=fabs(x) */
  if(y<0.0) y=-y; /* y=fabs(y) */

  if(NORM==1) return (x+y);    /* norme L_1 */
  if(NORM==3) return fmax(x,y); /* norme L_max */
  if(NORM==4) return fmin(x,y); /* norme L_min */
  return 0.0; /* si aucune NORM trouvée */
}


double distgone(int u,int v,int i,int p,int k,double w){
/*
  Calcule la distance P_i(u,v). Il s'agit de la "distance p-gone" (un
  polygone régulier à p cotés) relative à la direction i (axe d'angle
  i*2pi/k) entre les points u et v, restreint au cône de visibilité
  d'angle w*(p-2)*pi/p (d'angle w*pi si p est infini, c'est-à-dire
  p<3), w=0...1 (voir aussi la définition du thetagone dans l'aide en
  ligne). Ici, k est le nombre de directions.  L'algorithme est en
  O(1).

  Soient a_j (j=0...p-1) les sommets du p-gone P_i de rayon unité avec
  a_0=u et numérotés consécutivement en tournant dans le sense
  positif. Soit c le centre de P_i. Donc dist(u,c)=1, les a_j sont sur
  un cercle de rayon unité. On remarque que l'angle (u,a_j) est
  indépendant du rayon du p-gone. En fait, l'angle (a_j,u,a_{j+1})
  vaut la moitié de l'angle (a_j,c,a_{j+1}), soit pi/p. Et donc
  l'angle entre deux cotés consécutifs d'un p-gone vaut
  (p-2)*pi/p. Pour le calcul de distgone(u,v) on fait:

  1. Trouver la direction j tq (u,v) est dans le cône de visibilité et
     dans la région [(u,a_j),(u,a_{j+1})[.  Si j n'existe pas, alors
     on renvoit une distance infinie. Si p est infini, a_j est
     simplement sur la droite (u,v).

  2. On calcule l'intersection v' entre les droites (u,v) et
     (a_j,a_{j+1}). Si p est infini, v'=a_j. L'intersection existe
     forcément. Eventuellement v est sur la droite (u,a_j).

  3. distgone(u,v)=dist(u,v)/dist(u,v').

*/
  int j;
  double xu,xv,dxc,dxa,dxb,dxv;
  double yu,yv,dyc,dya,dyb,dyv;
  double hv,A,Ac,Ap,Aw;

  xu=XPOS[u];yu=YPOS[u]; /* coordonnées de u */
  xv=XPOS[v];yv=YPOS[v]; /* coordonnées de v */
  Ac=(double)i*2.0*M_PI/(double)k; /* angle (u,c), c=centre de P_i */
  dxc=cos(Ac);dyc=sin(Ac); /* coordonnées du centre c dans le repère u */

  dxv=xv-xu;dyv=yv-yu; /* coordonnées de v dans repère u */
  hv=hypot(dxv,dyv); /* |v-u|=dist(u,v) */
  if(hv==0.0) return 0.0; /* si u,v ont les mêmes coordonnées */

  /*
    Rappel: Si a et b sont deux vecteurs, alors le produit scalaire
    (dot product) est le réel a.b = xa*xb + ya*yb = |a|*|b|*cos(a,b),
    où |a|=sqrt(xa^2 + ya^2). Donc le signe de cos(a,b) est celui de
    xa*xb + ya*yb. Pour calculer sin(a,b) il faut faire une rotation
    de +pi/2 au vecteur a et calculer cos(a',b) où a'=(-ya,xa). Donc
    sin(a,b) = cos(a',b) = (-ya*xb + xa*yb) / (|a|*|b|). Et le signe
    de sin(a,b) est celui de xa*yb - ya*xb.
  */

  /* Aw=demi-angle du cône de visibilité */
  Aw=0.5*w*M_PI; /* si p infini */
  if(p>2) Aw *= (double)(p-2)/(double)p; /* si p est fini */

  /*
    Il faut bien sûr que (u,v) soit dans le cône de visibilité. La
    bissectrice de ce cône est l'axe (u,c) et son angle est
    w*(p-2)pi/p (w*pi si p infini). On note (w1,u,w2) le cône en
    question. Il faut que (u,v) soit entre (u,w1) (compris) et (u,w2)
    (non compris). Donc si sin(w1,u,v) < 0 ou si sin(w2,u,v) > 0 alors
    il n'existe pas de j (et donc on retourne une distance infinie).
  */

  A=Ac-Aw; /* A=angle (c,w1) */
  dxa=cos(A);dya=sin(A); /* coordonnées de w1 relatif à u */
  if(dxa*dyv<dya*dxv) return DBL_MAX; /* v avant w1 */

  A=Ac+Aw; /* A=angle (c,w2) */
  dxa=cos(A);dya=sin(A); /* coordonnées de w2 relatif à u */
  if(dxa*dyv>=dya*dxv) return DBL_MAX; /* v après ou sur w2 */

  /*
    Ici v est dans le cône de visibilité et donc la droite (uv)
    intersecte P_i en un point v'.
  */

  Ac -= M_PI; /* Ac=Ac-pi */

  /* Cas p infini */
  if(p<3){
    /*
      On raisone dans le repère u.  On pose c'=(dyc,-dxc),
      c'est-à-dire une rotation de -pi/2 de (u,v). On a |uc'|=1. On
      calcule l'angle A=(uv,uc'), en fait cos(A). On obtient v' en
      tournant autour de c d'un angle Ac-M_PI+2A ...
     */
    A=acos((dxv*dyc-dyv*dxc)/hv);
    A=Ac+2.0*A;
    dxa=dxc+cos(A);dya=dyc+sin(A);
    return hv/hypot(dxa,dya);
  }

  /*
    Cas p fini.  On cherche j de sorte qu'en tournant dans le sens
    positif, le vecteur (u,v) soit compris entre (u,a_j) (compris) et
    (u,a_{j+1}) (non compris). La droite (u,v) intersecte le segment
    [a_j, a_{j+1}[. L'indice j recherché est l'entier tq: (j-1)*pi/p
    <= angle(a_1,u,a_j) < j*pi/p. Et donc, j = 1 +
    floor{arcos(a_1,u,v)/(pi/p)}.
  */

  Ap=2.0*M_PI/(double)p; /* valeur souvent utilisée */
  A=Ac + Ap; /* angle (c,a_1) */
  dxa=dxc+cos(A);dya=dyc+sin(A); /* coordonnées de a_1 relatif à u */

  /* Aw=cos(a_1,u,v) = (a_1-u).(v-u) / dist(a_1,u)*dist(u,v) */
  Aw=(dxa*dxv+dya*dyv)/(hypot(dxa,dya)*hv);
  j=(int)((acos(Aw)*(double)p)/M_PI); /* en fait, la variable j vaut "j-1" */
  A += Ap*(double)j; /* angle (c,a_j): on part de a_1, donc on décale de j-1 cônes */
  dxa=dxc+cos(A);dya=dyc+sin(A); /* coordonnées de a_j relatif à u */
  A += Ap; /* angle (c,a_{j+1}) */
  dxb=dxc+cos(A)-dxa;dyb=dyc+sin(A)-dya; /* vecteur (a_j,a_{j+1}) */

  /*
    Calcule l'unique intersection entre la droite Dv définie par le
    vecteur (u,v) et la droite Dj définie par le vecteur
    (a_j,a_{j+1}). Dans le repère u, les équations sont:

    Dv: dyv*X - dxv*Y = 0
    Dj: dyb*X - dxb*Y = B avec B=dyb*dxa-dxb*dya,

    car a_j=(dxa,dya) appartient à Dj.
    L'intersection (x0,y0) est, dans le repère u:
    en faisant dxv*Dj-dxb*Dv, on a: x0=dxv*B/(dxv*dyb-dxb*dyv)
    en faisant dyb*Dv-dyv*Dj, on a: y0=dyv*B/(dxv*dyb-dxb*dyv)
  */

  A=(dyb*dxa-dxb*dya)/(double)(dxv*dyb-dxb*dyv); /* A=B/(...) */
  return hv/hypot(dxv*A,dyv*A);
}


enum{
  DYCK_WALK,
  DYCK_WORD,
  DYCK_TREE,
  DYCK_KTREE,
};


int* Dyck(int *R,int n,int k,int code)
/*
  Construit de manière aléatoirement uniforme un mot de k-Dyck composé
  de n segments montant (codé par k), chacun de longueur k, et de kn
  pas descendant (codé par 0). On doit avoir n>0 et k>0. Le résultat
  est donc en général un mot de n+kn lettres composés de n valeurs k
  et de kn valeurs 0. La propriété de ce mot est que tout préfixe
  possède au moins autant de pas montant que descendant. Le mot
  classique de Dyck et les arbres binaires sont obtenus pour k=1.

  Le résultat est renvoyé et mis dans R, qui est alloué si R=NULL.
  Suivant la valeur de code on renvoit des variantes du mot de
  Dyck. Attention, la taille de R, fonction de n et k, varie avec code
  (voir ci-dessous).

    code    taille  résultat
  -------------------------------------------------------------------
  DYCK_WALK  n+nk   mot ayant n valeurs k et kn valeurs 0
  DYCK_WORD  n+nk   même mot mais décalé sur un minimum
  DYCK_TREE  nk+1   arbre DFS à nk+1 noeuds
  DYCK_KTREE n+nk+1 arbre (k+1)-ary à n+nk+1 noeuds dont n de degré k+1
  -------------------------------------------------------------------

  On peut obtenir un mot de longueur p ayant q valeurs non nul
  aléatoire uniforme avec DYCK_WALK en choisissant n=q et k=p/q-1.

  Exemple pour n=5 et k=1:
  code=DYCK_WALK -> 1101001010
                                            3   4           3   4
                       /\/\                / \ / \           \ /
                    /\/    \/\        1   2   2   2   5     1 2 5
  code=DYCK_WORD -> 1011010010       / \ /         \ / \     \|/
  code=DYCK_TREE -> [-1,0,0,2,2,0]  0   0           0   0     0

  L'algorithme principal, en O(kn), consiste à tirer aléatoirement n
  valeurs k et k*n valeurs 0. On tire k avec une probabilité
  proportionelle au nombre de valeurs k restant à tirer sur le nombre
  total de valeurs restant à tirer. C'est un tirage aléatoire
  uniforme.

  L'arbre DFS (DYCK_TREE) est construit de sorte que le mot de Dyck
  forme le parcours DFS de l'arbre: 1 on descend le long d'une arête,
  et 0 on en revient. Les propriétés de cet arbre sont nombreuses. NB:
  Si k>1, on fait comme si on avait k pas de 1. Par exemple:
  - R[0] = -1, car 0 est la racine.
  - nk = dernière feuille de l'arbre, R[nk] est donc sont père
  - si u-R[u]>1 alors u démarre une nouvelle branche

  Exemple pour n=2 et k=2: arbre ternaire à 2 sommets internes

                         /\           4 5 6
                      /\/  \           \|/
                     /      \         1 2 3
  code=DYCK_WORD  -> 2 02 000          \|/
  code=DYCK_KTREE -> [-1,0,0,0,2,2,2]   0

  Exemple pour n=3 et k=1: arbre binaire à 3 sommets internes
  code=DYCK_WORD  -> 110100
  code=DYCK_KTREE -> [-1,0,0,1,1,4,4]

        5   6
         \ /
      3   4
       \ /
        1   2
         \ /
          0

  Exemple pour n=3 et k=2: arbre ternaire à 3 sommets internes
  code=DYCK_WORD  -> 202200000
  code=DYCK_KTREE -> [-1,0,0,0,2,2,2,4,4,4]

      7 8 9
       \|/
        4 5 6
         \|/
        1 2 3
         \|/
          0

  Pour construire l'arbre (k+1)-ary (DYCK_KTREE) on procède selon un
  DFS modifié (on pose les fils avant la récursion), en lisant
  succéssivement les n+kn valeurs du mot de k-Dyck. Plus précisément,
  depuis le sommet courant u: si on lit k, on ajoute k+1 fils à u, le
  nouveau sommet courant devenant le 1er fils de u. Si on lit 0, on
  remonte les ancêtres de u jusqu'à trouver un ancêtre v dont le fils
  f menant à u ne soit pas le dernier fils (le (k+1)-ième). Le nouveau
  sommet courant est alors f+1. Pour déterminer si f est le dernier
  fils (ou pas) il suffit de tester si f%(k+1)=0 ou pas.
*/
{
  const int m=n+k*n; /* longueur du mot de Dyck */
  int *B; /* mot de Dyck */
  int i;

  if(R==NULL){
    if(code==DYCK_WALK)  ALLOC(R,m);
    if(code==DYCK_WORD)  ALLOC(R,m);
    if(code==DYCK_TREE)  ALLOC(R,m-n+1);
    if(code==DYCK_KTREE) ALLOC(R,m+1);
  }
  if(code==DYCK_WALK) B=R; else ALLOC(B,m);

  /* DYCK_WALK: construit B */
  /* commun à DYCK_WORD et DYCK_TREE */
  
  int t=n; /* t=nombre de k à tirer */
  int s=m; /* s=nombre total de valeur à tirer */

  for(i=0;s>0;s--,i++){
    B[i]=((random()%s)<t); /* k avec proba proportionnel à t */
    if(B[i]) B[i]=k,t--; /* enlève un k */
  }
  if(code==DYCK_WALK) return R; /* NB: ici R=B */
  
  /* cherche la position r dans B de la hauteur minimum */
  /* commun à DYCK_WORD et DYCK_TREE */

  int r=-1;
  int h=0;
  
  for(i=s=0;i<m;i++){
    s += B[i]? k : -1; /* s=hauteur courante = +k ou -1 */
    if(s<h) h=s,r=i;   /* h=hauteur minimum */
  }
  r=(r+1)%m; /* la position r recherchée */
  
  /* DYCK_WORD: décale le mot de B vers R */
  
  if(code==DYCK_WORD){
    for(i=0;i<m;i++) R[i]=B[(i+r)%m];
    goto fin_dyck;
  }

  /* DYCK_TREE: décale et construit un arbre DFS, R[v]=parent du sommet v */
  /* Si k>1, on fait comme si on avait k pas de 1 */

  if(code==DYCK_TREE){
    int u=0; /* u=dernier sommet visité (=racine) */
    int v=1; /* v=prochain sommet à visiter (=sommet 1) */
    R[0]=-1; /* père de la racine */
  
    for(i=0;i<m;i++) /* parcoure toutes les valeurs de B */
      if(B[t=(i+r)%m])
	while(B[t]){ /* tantque c'est > 0, on monte */
	  R[v]=u; /* père de v=dernier sommet visité */
	  u=v++; /* met à jour le dernier sommet visité */
	  B[t]--;
	}
      else u=R[u]; /* si c'est un 0, on descend */
    goto fin_dyck;
  }

  /* DYCK_KTREE: on a n noeuds internes de degré k+1 */

  if(code==DYCK_KTREE){
    int u=0,s=1; /* u=sommet courant, indexe dans R */
    R[0]=-1;     /* père de la racine */
    for(i=0;i<m;i++) /* parcoure toutes les valeurs de B */
      if(B[(i+r)%m]){
	for(t=0;t<=k;t++) R[s++]=u; /* si on lit k, on pose k+1 fils */
	u=s-t; /* le nouveau sommet courant est le 1er fils = s-k-1=s-t */
      }
      else{
	while(u%(k+1)==0) u=R[u]; /* si on lit 0, on remonte tous les derniers fils */
	u++; /* le nouveau sommet courant est u+1 */
      }
    goto fin_dyck;
  }

 fin_dyck:
  free(B);
  return R;
}


int NextPermutation(int *P,int n,int *C)
/*
  Génère, à partir d'une permutation P, la prochaine dans l'ordre
  lexicographique suivant les contraintes définies par le tableau C
  (voir ci-après). Mettre C=NULL s'il n'y a pas de contrainte
  particulière. On renvoie 1 ssi la dernière permutation a été
  atteinte. Dans ce cas la permutation la plus petite selon l'ordre
  lexicographique est renvoyée. On permute les éléments de P que si
  leurs positions sont entre C[j] et C[j+1] (exclu) pour un certain
  indice j. On peut initialiser P avec ALLOCZ(P,k,_i) ou si le tableau
  P existe déjà avec NextSet(P,-1,k).

  Ex: C={2,3,5}. Les permutations sous la contrainte C sont:
  (on peut permuter les indices {0,1}{2}{3,4})

                 0 1 2 3 4 (positions dans P)
	      P={a,b,c,d,e}
	        {b,a,c,d,e}
		{a,b,c,e,d}
		{b,a,c,e,d}
  
  Evidemment, il y a beaucoup moins de permutations dès que les
  contraintes augmentent. Par exemple, si C contient k intervalles de
  même longueur, alors le nombre de permutations sera de (n/k)!^k au
  lieu de n!. Le rapport des deux est d'environ k^n.

  Concrêtement, pour:
  - n=9 et k=3, on a 216 permutations au lieu de 362.880 (k^n=19.683)
  - n=12 et k=3, on a 13.824 permutations au lieu de 479.001.600 (k^n=531.441)

  Le dernier élément de C doit être égale à n-1 (sentinelle), le
  premier étant omis car toujours = 0. Donc C est un tableau à au plus
  n éléments. Si C=NULL, alors il n'y a pas de contrainte
  particulière, ce qui est identique à poser C[0]=n.

  On se base sur l'algorithme classique (lorsque C=NULL, sinon on
  l'applique sur l'intervalle de positions [C[j],C[j+1][):

  1. Trouver le plus grand index i tel que P[i] < P[i+1].
     S'il n'existe pas, la dernière permutation est atteinte.
  2. Trouver le plus grand indice j tel que P[i] < P[j].
  3. Echanger P[i] avec P[j].
  4. Renverser la suite de P[i+1] jusqu'au dernier élément.

*/
{
  int i,j,a,b,c,T[1];

  if(C==NULL){
    T[0]=n;
    C=T;
  }

  b=C[i=j=0]; /* j=indice de la prochaine valeur à lire dans C */
  c=-1;

  /* étape 1: on cherche l'intervalle [a,b[ avec i tq P[i]<P[i+1] */
 etape1:
  for(a=i;i<b-1;i++) if(P[i]<P[i+1]) c=i; /* on a trouvé un i tq P[i]<P[i+1] */
  if(c<0){ /* le plus grand i tq P[i]<[i+1] n'existe pas */
    for(i=a;i<b;i++) P[i]=i; /* on réinitialise P[a]...P[b-1] */
    if(b==n) return 1; /* alors on a fini d'examiner C */
    b=C[++j]; /* b=nouvelle borne */
    goto etape1;
  }
  i=c; /* i=le plus grand tq P[i]<P[i+1] avec a<=i,i+1<b */

  /* étape 2: cherche j=le plus grand tq P[i]<P[j] */
  for(j=i+1;j<b;j++) if(P[i]<P[j]) c=j;
  j=c;

  /* étape 3: échange P[i] et P[j] */
  SWAP(P[i],P[j],c);

  /* étape 4: renverse P[i+1]...P[b-1] */
  for(++i,--b;i<b;i++,b--) SWAP(P[i],P[b],c);

  return 0;
}


int NextSet(int *S,int n,int k){
/*
  Calcule les sous-ensembles de k entiers de [0,n[. Si n<0, alors S
  est initialisé au plus petit ensemble possible, soit S={0,1,2,
  ...,k-1}. L'idée est de maintenir une sorte de compteur S qui
  représente le sous-ensemble courant et qui va passer par tous les
  sous-ensembles possibles. Les éléments sont rangés dans l'ordre
  croissant. On renvoie 1 ssi S était le dernier sous-ensemble. Dans
  ce cas l'ensemble le plus petit est écrit dans S.

  On peut facilement en déduire un algorithme pour générer des
  multi-ensembles, comme par exemple tous les multi-ensembles du type
  [2,0,0,2,1,0] comprennant 3 fois 0, 1 fois 1 et 2 fois 2:
  NextMultiSet(S,C,k) avec C=[3,1,2] (voir la fonction ggosset()).

  La stratégie pour "incrémenter" S est la suivante : on essaye
  d'incrémenter le plus petit élément de S[i] tout en restant
  strictement inférieur à l'élément suivant S[i+1]. Si c'est possible
  on a trouvé le nouveau sous-ensemble. Sinon, on réinitialise
  S[0],...,S[i] à leur plus petites valeurs: S[0]=0,...,S[i]=i. Si on
  a pas pu incrémenter S[k-1] sans dépasser n on a atteint le dernier
  sous-ensemble.

  L'algorithme est en O(k) dans le pire des cas, mais de manière
  amortie c'est beaucoup moins car on incrémente moins souvent S[j]
  que S[i] si j>i.
*/
  int i=0,j,s;

  if(n<0){
    for(;i<k;i++) S[i]=i;
    return 0;
  }
  
  while(i<k){
    s=++S[j=i++];
    if(i==k){
      if(s<n) return 0;
      goto fin;
    }
    if(s<S[i]) return 0;
  fin:
    S[j]=j;
  }
  return 1;
}


int NextArrangement(int *S,int *P,int n,int k){
/*
  Permet de générer tous les arrangements de k entiers de
  [0,n[. L'arrangement est représenté par les tableaux S et P de k
  entiers. S représente un ensemble de k entiers et P une permutation
  de [0,k[. Ainsi, l'arrangement A=(4,2,7,3) est représenté par
  S=(2,3,4,7) et P=(2,0,3,1). Autrement dit A[i]=S[P[i]] pour tout
  i=0...k-1.

  L'idée est d'essayer d'incrémenter le double compteur S,P. On essaye
  d'incrémenter P en premier avec NextPermutation(). Si on est arrivé
  à la fin de P, on incrémente S avec NextSet(). Si n<0, alors S et P
  sont initialisés au plus petit arrangement possible, soit
  S=P=(0,1,2, ...,k-1). On renvoie 1 ssi S,P était le dernier
  arrangement. Dans ce cas l'arrangement le plus petit est écrit dans
  S,P.
*/
  if(n<0){
    int i;
    for(i=0;i<k;i++) S[i]=P[i]=i;
    return 0;
  }
  if(NextPermutation(P,k,NULL)) return NextSet(S,n,k);
  return 0;
}


int *NextPart(int *S,int n,int s,int *C)
/*
  Permet de générer toutes les suites S de n>0 entiers >=0 dont la
  somme fait s et dont la i-ème part S[i] ne dépasse pas C[i]. Il faut
  que s <= sum_{i=0}^(n-1) C[i], n=1 est possible.

  S est la suite courante de somme s et on renvoie dans S la prochaine
  suite (la fonction renvoie aussi S). On renvoie NULL si on a atteint
  la dernière suite, et on remplit S avec le première suite. Si
  S=NULL, alors S est allouée et initialisée à la première suite. La
  première suite de somme s est obtenue en remplissant autant que
  possible les parts S[n-1],S[n-2],...

  L'algorithme est le suivant:
   1. on détermine le plus grand indice j tq S[j]>0
   2. si j=0, alors on a finit: on fait 5. avec x=s+1 et i=-1
   3. on détermine le plus grand indice i<j tq S[i] peut être augmenté
   4. on calcule x = sum_{j=i+1}^(n-1) S[i]
   5. on incrémente S[i]
   6. on remplit S[i+1]...S[n-1] avec la première suite de somme x-1

  Exemple: s=n=5

  C=1 2 2 1 1
  S=0 1 2 1 1
    0 2 1 1 1
    0 2 2 0 1
    0 2 2 1 0
    1 0 2 1 1
    1 1 1 1 1
    1 1 2 0 1
    1 1 2 1 0
    1 2 0 1 1
    1 2 1 0 1
    1 2 1 1 0
    1 2 2 0 0
*/
{
  int x,i,j,r;

  i=0;
  r=(S==NULL);
  if(r) ALLOC(S,n);
  else i=n-1;
  
  /* calcule le plus grand indice i tq S[i]>0 */
  while((i>0)&&(S[i]==0)) i--;
  x=S[i--]; /* rem: si i=0, alors i -> -1 */

  /* calcule le plus grand indice j<i tq S[j] peut être augmenté */ 
  while((i>=0)&&(S[i]==C[i])) x += S[i--];

  if(i>=0){ S[i]++; s=x-1; } /* si cet indice n'existe pas i<0 => FIN */
  
  /* écrit la première suite de somme s dans S[i+1]...S[n-1] */
  for(j=n-1;j>i;j--){
    x=max(s,0);
    x=min(C[j],x);
    S[j]=x;
    s -= x;
  }

  /* on retourne S sauf si i<0 et r=0 (<=> FIN ) */
  return ((i<0)&&(!r))? NULL : S;
}


int SetCmp(int *T1,int *T2,int n1,int n2)
/*
  Compare deux tableaux d'entiers T1 et T2 de taille n1 et n2 triés
  par ordre croissant. Les tableaux peuvent être de taille nulle. La
  valeur renvoyée est un entier interprété en binaire comme suit:

  bit-0: 1 ssi T1 intersecte T2 (possède au moins 1 élément commun)
  bit-1: 1 ssi T1 égale T2
  bit-2: 1 ssi T1 est strictement inclu dans T2
  bit-3: 1 ssi T2 est strictement inclu dans T1

  Les valeurs possibles sont donc: 0,1,2,3,4,5,8,9 (=[0,9]\{6,7})
  La complexité est O(n1+n2).
*/
{
  if(n1==0) return (n2==0)? 2:4;
  if(n2==0) return 8;
  /* ici T1 et T2 contiennent au moins 1 élément */

  if((T1[n1-1]<T2[0])||(T2[n2-1]<T1[0])) return 0; /* cas trivial de disjonction */

  int i1,i2,r;
  i1=i2=0;
  r=14; /* tous les bit à 1 sauf b0 */

  while((i1<n1)&&(i2<n2)){
    if(T1[i1]==T2[i2]){
      i1++; i2++; /* T1 et T2 ont un élément en commun */
      r |= 1; continue; /* met b0 */
    }
    r &= 13; /* annule b1 (15-2) car T1<>T2 */
    if(T1[i1]<T2[i2]){
      i1++; /* T1 contient des éléments qui ne sont pas dans T2 */
      r &= 11; /* annule b2 (15-4) car T1 ne peux pas contenir T1 */
    }else{
      i2++; /* T2 contient des éléments qui ne sont pas dans T1 */
      r &= 7; /* annule b3 (15-8) car T2 ne peux pas contenir T1 */
    }
  }

  if(i1<n1) r &= 9; /* annule b2 et b1 (15-4-2) */
  if(i2<n2) r &= 5; /* annule b3 et b1 (15-8-2) */
  if(r&2)   r &= 3; /* annule b3 et b2 (15-8-4) */

  return r;
}


int SetSearch(int u,int *T,int n,int sort){
/*
  Cherche l'entier u dans le tableau T de taille n. Si u est dans T,
  on renvoie son indice sinon on renvoie -1. Si sort=1, alors on
  suppose T trié par ordre croissant et la recherche est dichotomique,
  sinon la recherche est linéaire.
*/
  if(sort){ /* recherche dichotomique */
    int *t=bsearch(&u,T,n,sizeof(int),fcmp_int);
    return t? (int)(t-T) : -1;
  }
  /* recherche linéaire */
  int i;
  for(i=0;i<n;i++) if(u==T[i]) return i;
  return -1;
}


int Binom(int n,int k)
/*
  Calcule l'entier B={n choose k}.  L'algorithme utilisé ici est en
  O(k). Il n'utilise que des multiplications et divisions entières sur
  des entiers en O(B), sans aucun tableau.

  L'algorithme classique issu du Triangle de Pascal, qui lui est en
  O(n*k), utilise un espace en O(k) (tableaux d'entiers en O(B)). Par
  contre il n'utilise que des additions sur des entiers en O(B).

  Principe: B = n x (n-1) x ... x (n-k+1) / k x (k-1) x ... x 1

  On réecrit B en (((((n/1) x (n-1)) / 2) x (n-2)) / 3) ...
  c'est-à-dire en multipliant le plus grand numérateur et en divisant
  par le plus petit numérateur. Chaque terme est ainsi un certain
  binomial, et donc toujours un entier.

  Catalan(n) = Binom(2n,n)/(n+1). On peut aussi calculer ce nombre
  avec Catalan(0)=1 et Catalan(n) = (2*(2n-1)*Catalan(n-1)/(n+1).
*/
{
  int B,s,i;
  for(B=s=n,i=2;i<=k;i++) B=(B*(--s))/i;
  return B;
}

/* code pour la fonction SortInt() */
enum{
  SORT_INC,     /* tri croissant selon les valeurs */
  SORT_DEC,     /* tri décroissant selon les valeurs */
  SORT_INC_RANK,/* donne le rang des valeurs de T dans l'ordre croissant */
  SORT_DEC_RANK,/* donne le rang des valeurs de T dans l'ordre décroissant */
  SORT_FREQv,   /* donne la fréquence des valeurs de T */
  SORT_FREQe,   /* donne la fréquence des éléments de T */
  SORT_INDEXi,  /* donne l'indice des éléments de T dans l'ordre croissant */
  SORT_INDEXe   /* donne les éléments de T dans l'ordre croissant */
};

int *SortInt(int *T,int *R,int n,int a,int *m,const int code)
/*
  Trie par ordre croissant un tableau T non NULL de taille n>0 dont
  les valeurs sont des entiers de [a,a+*m[ si m<>NULL. Sinon, les
  valeurs de T sont supposées être dans [a,a+n[. Pour simplifier le
  texte qui suit, je note m pour dire *m.

  La complexité en temps est O(n+m), ce qui est mieux que qsort(). Le
  tableau T n'est pas modifié. Le résultat est rangé dans le tableau R
  qui doit être de taille au moins n, et est aussi renvoyé par la
  fonction. Il R=NULL, alors R est alloué et retourné par la
  fonction. Pour être compétitif avec qsort() pour m=n il faut,
  théoriquement, n>32.

  L'opération de tri dépend de "code" qui est interprétée comme suit
  (la complexité est le nombre d'étapes dans les boucles "for"):

    v = valeur d'un élément de T, v dans [a,a+m[
    d = nombre de valeurs v distinctes de T, d dans [1,min(n,m)]
    e = indice du tableau T, e dans [0..n[
    i = indice du tableau R, i dans [0..n[
    r = rang dans un tableau, r dans [0..d[

  - code=SORT_FREQv: renvoie un tableau F[0..m[ de fréquence des
    valeurs de T, c'est-à-dire que F[i] est le nombre de fois où la
    valeur i+a apparaît dans T. Le tableau R et la variable m ne sont
    pas modifiés. Complexité: m+n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        F = [0 1 3 2 0 1 1] de taille m

  - code=SORT_FREQe: renvoie dans R[0..n[ un tableau de fréquence des
    éléments de T, c'est-à-dire où R[e] est le nombre de fois où T[e]
    apparaît dans T. La variable m n'est pas modifiée. Complexité:
    m+2n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        R = [ 3  1  3  2  3  1  1  2] de taille n

  - code=SORT_INC ou SORT_DEC: renvoie dans R[0..n[ le tableau T trié
    par ordre croissant (ou décroissant). La variable m n'est pas
    modifiée. Complexité: 2m+2n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        R = [11 12 12 12 13 13 15 16] de taille n

  - code=SORT_INC_RANK ou SORT_DEC_RANK: renvoie dans R[0..n[ un
    tableau de rangs où R[e] est le rang r dans [0..d[ de l'élément
    T[e] dans la version triée dans l'ordre croissant (ou décroissant)
    de T. Le tableau R est modifié et on renvoie d dans m. Complexité:
    3m+2n.

    Ex: T = [12 11 12 13 12 15 16 13]  avec n=8, a=10, m=7
        R = [ 1  0  1  2  1  3  5  2] et m=5

  - code=SORT_INDEXi: renvoie dans R[0..n[ un tableau d'indices où
    R[i]=e est l'élément de T en position i dans la version triée par
    ordre croissant de T. Pour obtenir un tri de T il suffit de lister
    T[R[i]] pour i=0..n-1. La variable m n'est pas modifiée.
    Complexité: 2m+2n.

    Ex: T = [12 11 12 13 12 18 15 11]  avec n=8, a=10, m=9
        R = [ 1  7  0  2  4  3  6  5]

  - code=SORT_INDEXe: renvoie dans R[0..n[ un tableau d'indices où
    R[e] est la position de T[e] dans la version triée par ordre
    croissant de T. La variable m n'est pas modifiée. Complexité:
    2m+2n.

    Ex: T = [12 11 12 13 12 18 15 11]  avec n=8, a=10, m=9
        R = [ 2  0  3  5  4  7  6  1]
*/
{
  int i,r,t;
  int k=(m==NULL)? n:*m;

  /* initialise F[0..m[ */
  NALLOCZ(int,F,k,0); /* coût: m */
  
  /* calcule F[i]=fréquence de la valeur v=i+a dans T */
  for(i=0;i<n;i++) F[T[i]-a]++; /* coût: n */

  if(code==SORT_FREQv) return F;

  /* alloue R, si nécessaire */
  if(R==NULL) ALLOC(R,n);

  if(code==SORT_FREQe){
    for(i=0;i<n;i++) R[i]=F[T[i]-a]; /* coût: n */
    free(F); return R;
  }

  if(code==SORT_INC){ /* R=tableau T trié, ordre croissant */
    for(i=r=0;i<k;i++) for(t=F[i];t>0;t--) R[r++]=i+a; /* coût: m+n */
    free(F); return R;
  }

  if(code==SORT_DEC){ /* R=tableau T trié, ordre décroissant */
    for(i=0,r=n;i<k;i++) for(t=F[i];t>0;t--) R[--r]=i+a; /* coût: m+n */
    free(F); return R;
  }

  /* calcule F[i]=nombre de valeurs de T qui sont < i+a.     
     Ex:  T = [2 1 2 3 2 8 5]  avec n=7, m=9, a=0
	  F = [0 0 1 4 5 5 6 6 6]
  */
  for(i=r=0;i<k;i++){ t=F[i]; F[i]=r; r += t; } /* coût: m */

  if(code==SORT_INDEXi){
  for(i=0;i<n;i++) R[F[T[i]-a]++]=i; /* coût: n */
    free(F); return R;
  }

  if(code==SORT_INDEXe){
    for(i=0;i<n;i++) R[i]=F[T[i]-a]++; /* coût: n */
    free(F); return R;
  }

  /* pour SORT_INC_RANK ou SORT_DEC_RANK */
  /* calcule dans F[i] le rang dans [0,d[ du nb de valeurs < i+a */
  
  for(t=r=-1,i=0;i<k;i++){ /* coût: m */
    if(F[i]!=t) { t=F[i]; r++; }
    F[i]=r;
  }
  if(m!=NULL) *m=r+1; /* m=nb de valeurs différentes de T */

  if(code==SORT_INC_RANK){
    for(i=0;i<n;i++) R[i]=F[T[i]-a]; /* coût: n */
    free(F); return R;
  }
  
  if(code==SORT_DEC_RANK){
    for(i=0;i<n;i++) R[i]=r-F[T[i]-a]; /* coût: n */
    free(F); return R;
  }

  free(F);
  free(R);
  Erreur(21);
  return NULL;
}


/***********************************

       ROUTINES SUR
       LES GRAPHES

***********************************/


void PrintGraph(graph *G)
/*
  Affiche un graphe ou une famille de graphes au format standard sous
  forme compacte. Utilise WIDTH. Effet de bord: le (ou les) graphes
  sont triés, et donc G->sort=1 en sortie.
*/
{
  if(G==NULL){ printf("NULL\n"); return; }

  int u,v,i,k,n,ligne,nk=(G->f>0)?G->f:1;
  graph *H;
  int *P;

  for(k=0;k<nk;k++){

    if(G->f>0){
      H=G->G[k];
      printf("[%i]",H->id);
    } else H=G;

    SortGraph(H,0);
    n=H->n;
    ALLOCZ(P,n,0);
    i=u=ligne=0;
    v=-1;

    while(i<n){
      /* si u==i, alors u=tête d'un nouveau chemin */
      while((v<u)&&(P[u]<H->d[u])) v=H->L[u][P[u]++];
      if(v<u){ /* on a pas trouvé v => fin d'un chemin */
	if(H->d[u]==0){ /* cas des sommets isolés */
	  printf(" %i",u);
	  if(++ligne==WIDTH){ printf("\n"); ligne=0; }
	}
	u=(i+=(u==i));
	v=-1;
      }
      else{ /* u a un voisin v>u */
	if((u==i)||(ligne==0)) printf(" %i",u); /* on affiche la tête */
	printf("-%s%i",(G->sym)?"":">",v); /* on affiche -v ou ->v */
	if(++ligne==WIDTH){ printf("\n"); ligne=0; }
	u=v; /* on continu avec v */
	v=-1;
      }
    } /* fin du while */

    if(ligne>0) printf("\n"); /* newline si fini avant la fin de ligne */
    free(P);
  }

  G->sort=1; /* effet de bord */
  return;
}


int NbEdges(graph *G)
/*
  Retourne le nombre d'arêtes d'un graphe symétrique G ou bien le
  champs G->m s'il est positif. Si G->m<0, alors G->m est mis à jour à
  partir de la somme des G->d[i].
*/
{ int m=G->m;
  if(m<0){
    int i,n=G->n;
    for(i=m=0;i<n;i++)
      m += G->d[i];
    G->m=(m>>=1);
  }
  return m;
}


int Degree(graph *G,int maxi)
/*
  Renvoie le degré maximum (si maxi=1) ou minimum (si maxi=0) d'un
  graphe G. On renvoie -1 si G est nul ou est une famille de graphes.
*/
{ if((G==NULL)||(G->f>0)) return -1;
  int i=1,n=G->n,d=G->d[0];
  for(;i<n;i++)
    if(maxi) d=max(d,G->d[i]);
    else d=min(d,G->d[i]);
  return d;
}


/***********************************

       BFS, DFS, ...
       (DIJKSTRA)

***********************************/


typedef struct{
  int root;  /* racine ou source du BFS. */
  int radius;/* eccentricité du sommet source. */
  int *D;    /* D[u]=distance entre u et root. Les sommets u avec
	        D[u]=-1 sont à une distance infinie de root (situation
	        initiale par défaut). En entrée, les sommets avec
	        D[u]=-2 doivent être considérés comme inexsitant dans
	        G. Si D=NULL, D est alloué puis initialisé à -1. */
  int *P;    /* P[u]=père de u dans un arbre BFS de racine root, avec
	        P[root]=-1. Seuls les sommets u avec D[u]>=0 ont un
	        père défini (sauf si u=root). Si P=NULL, alors P est
	        alloué. Il n'est pas nécessaire de l'initialiser. */
  int n;     /* nombre de sommets parcourus, nombre d'éléments dans la file */
  int *file; /* contenu de la file = liste des n sommets parcourus. La
		taille de ce tableau est toujours le nombre de sommets
		du graphe. */
  int cycle; /* longueur du plus petit cycle rencontré lors du
		parcours. Cette valeur (>2) est indépendente du
		parcours spécifique des voisins (ie de tout BFS). Elle
		ne dépend que de la structure non-étiquetée de G et de
		la source. Si cycle<0, alors cette longueur est
		infinie.  On peut déterminer la maille du graphe en
		prennant le minimum obtenu pour chacun des sommets du
		graphe comme source. Cette valeur n'est pas définie si
		G est orienté (-1). */
  int clean; /* permet d'initialiser à -1 les sommets parcourus dans
		le précédent bfs(). Par défaut, clean=0, et le tableau
		D n'est pas initialisé (sauf si D=NULL, dans ce cas il
		est initialisé complètement à -1). Si clean=2, alors
		on remet D[u] à -1 seulement pour les sommets u de
		file. C'est plus rapide si bfs() est utilisé plusieurs
		fois avec le même paramètre et tous les sommets ne
		sont pas visités (vol ou hmax >0). Si clean=1, alors
		on force l'initialisation complète de D à -1, puis on
		passe clean=2 (pour les appels suivants). */
  int tf;    /* tête de file, ie nombre d'éléments parcourus dont tous
		les voisins ont été enfilés (tf<=n). */
  int cont;  /* cont=1 si on poursuit le bfs() précédant là où on
		s'était arrêté dans le cas d'un arrêt par hmax. Par
		défaut cont=0. Après le 1er bfs() où cont=1, cont
		passe à 2. Cela sert à augmenter progressivement la
		hauteur jusqu'à obtenir une certaine condition. */
  /* le calcul des deux champs suivants rajoutent au plus n tests à la complexité */
  int vmax;  /* arrête le parcours lorsque vmax sommets on été parcourus */
  int hmax;  /* arrête le parcours lorsqu'un sommet de hauteur > hmax est atteint */
} param_bfs;


param_bfs *new_param_bfs(void)
/*
  Crée et initialise une structure pour la fontion bfs(). C'est bfs()
  qui se charge, éventuellement, d'allouer les tableaux P et D. Le
  tableau D peut être utilisé pour effacer des sommets.
*/
{
  NALLOC(param_bfs,X,1);
  X->D=X->P=X->file=NULL;
  X->radius=X->n=X->clean=X->cont=X->tf=0;
  X->cycle=X->root=X->vmax=X->hmax=-1;
  return X;
}


void free_param_bfs(param_bfs *X){
  if(X==NULL) return;
  free(X->D);
  free(X->P);
  free(X->file);
  free(X);
}


param_bfs *bfs(graph *G,int source,param_bfs *X)
/*
  Effectue un parcours en largeur (BFS) d'un graphe G (orienté ou non)
  depuis le sommet source. Les résultats du parcours (comme les
  distances, le père, le nombre de sommets parcourus, etc.) sont
  stockées dans la variable X qui est renvoyée. Si X=NULL, X est
  d'abord allouée.

  On peut également réaliser un BFS seulement sur un sous-ensemble de
  sommets (sous-graphe), limiter le parcours à une profondeur donnée
  (en fixant X->hmax), ou à un nombre de sommets parcourus (en fixant
  X->vmax). Il est important que la table X->D soient correctement
  initialisée pour les sommets à parcourir. Mêmes NULL, les tables
  X->P et X->file ne sont pas initialisés. Seules les valeurs des X->n
  sommets parcourus (ceux dans X->file) sont garanties. A noter que le
  parcours est lancé depuis la source, même si source est dans T.

  On l'utilise comme ceci:

    param_bfs *p=bfs(G,s,NULL); // BFS depuis le sommet s
    ... // p->D[u]=distance entre s et u
    ... // p->P[u]=père de u
    ... // p->cycle=longueur du plus petit cycle passant par s
    ... //
    free_param_bfs(p); // libère la variable crée p

  Pour réaliser un BFS d'un sous-graphe G\T de G (évitant les k
  sommets du tableau T):

    param_bfs *p=new_param_bfs();   // par défaut p->clean=0
    ALLOCZ(p->D,G->n,-1);           // alloue p->D et met tout à -1
    for(i=0;i<k;i++) p->D[T[i]]=-2; // les sommets à enlever doivent être à -2
    bfs(G,s,p);
    ... // p->D[u]=distance entre s et u dans G\T
    ... // p->P[u]=père de u, ou -1 s'il n'existe pas
    ... // p->cycle=longueur du plus petit cycle dans G\T passant par s
    ... //
    free_param_bfs(p);

  Pour faire des appels multiples à bfs() et optimiser la complexité:

    param_bfs *p=new_param_bfs();
    p->clean=1; // initialisation complète puis partielle de p->D
    ...
    p->vmax=100; // pour parcourir au plus 100 sommets
    bfs(G,u,p);  // initialise complètement p->D à -1 avant le bfs
    ...
    p->hmax=3;   // pour faire un bfs à distance au plus 3 de v
    bfs(G,v,p);  // p->D sera initialisé en initialisant seulement les sommets
    ...          // parcourus du bfs() précédant (donc au plus 100 sommets)
    ...
    bfs(G,w,p);  // initialisation partielle de p->D
    ...
    free_param_bfs(p);

  La complexité est en le nombre d'arêtes dans la boule des sommets
  parcourus à condition que X->D ne soit pas NULL, puisque sinon une
  initialisation de X->D en O(n) sera effectuée. C'est un point
  important si on lance plusieurs BFS partiels à partir de sommets
  d'un même graphe. Pour être efficace, à partir du 2e bfs() il est
  faut rétablir X->D[u]=-1 pour tous les sommets u de la X->file. On
  peut le réaliser en mettant X->clean=1. La complexité pour chaque
  appel (à part le premier qui initialise complètement X->D à -1)
  reste en le nombre d'arêtes de la boule des sommets parcourus.

  ALGORITHME:
  - prendre le sommet u en tête de file
  - pour tout ses voisins v non marqués:
    - enfiler v
    - marquer v
  
  Si D<>NULL, alors D n'est pas initialisé (sauf si clean>0). Dans ce
  cas, il faut que D[u]=-1 pour tout sommet u, sauf pour les sommets
  que l'on souhaite ne pas visiter où il faut D[u]<>-1 (typiquement,
  D[u]=-2). On peut initialiser partiellement ou complètement D avec
  clean=1 ou clean=2.

  Pour déterminer X->cycle, on calcule la longueur L du tout premier
  cycle crée. On remarque que X->cycle peut être peut-être L ou L-1.
  Il est L-1 si plus tard, sur le même niveau, on rencontre deux
  sommets voisins sur ce même niveau.
*/
{
  int i,u,v,d,ff,tf,h,n=G->n;
  int vmax,hmax;

  if(X==NULL) X=new_param_bfs(); /* NB: X->n=X->clean=0 */
  if(X->P==NULL) ALLOC(X->P,n); /* alloue tableau si P==NULL */
  if(X->file==NULL) ALLOC(X->file,n); /* alloue la file */
  
  if(X->D==NULL) ALLOCZ(X->D,n,-1); /* alloue et initialise D en O(n) */
  else{ /* initialisation partielle ou pas de D */
    if(X->cont<2){
      if(X->clean==1) for(u=0;u<n;u++) X->D[u]=-1; /* initialisation complète */
      if(X->clean==2) for(i=0;i<X->n;i++) X->D[X->file[i]]=-1; /* initialisation partielle */
    }
  }
  if(X->clean==1) X->clean=2; /* la prochaine fois, initialisation partielle de X->D */

  if(X->cont==2){
    tf=X->tf;
    ff=X->n;
  }else{
    tf=0; /* tf=tête de file, pointe sur la tête */
    ff=0; /* ff=fin de file, pointe sur prochain élément libre */
    if(X->cont==1) X->cont=2; /* la prochaine fois on continue */
    X->root=source; /* la racine est la source */
    X->P[source]=-1; /* pas de père pour la source */
    X->D[source]=0; /* distance=0 pour la source */
    X->file[ff++]=source; /* enfile le 1er sommet (=source), même s'il est supprimé */
    X->cycle=(G->sym)? n+1 : 0; /* X->cycle non défini (=0) si G orienté */
  }
  
  hmax=(X->hmax==-1)? n : X->hmax; /* si X->hmax non défini */
  vmax=(X->vmax==-1)? n : X->vmax; /* si X->vmax non défini */
  h=(G->sym)? 1+(X->cycle/2) : 0;   /* h=hauteur à partir de laquelle le plus court
				       cycle ne peut plus apparaître */

  while(tf<ff){
    u=X->file[tf]; /* défile la tête de file */
    if(X->D[u]>=hmax) break; /* fin si on a atteint une hauteur >=
				hmax. Dans ce cas, tous ceux de
				hauteur <= hmax on été enfilés. */
    tf++;
    for(i=0,d=G->d[u];i<d;i++){ /* pour tout voisin v de u */
      v=G->L[u][i];
      if(X->D[v]==-1){ /* si v voisin non marqué, si =-2 on saute le sommet */
	X->P[v]=u; /* le père de v est u */
	X->D[v]=X->D[u]+1; /* hauteur(v)=1+hauteur(père(v)) */
	X->file[ff++]=v; /* enfile v */
	if(ff>vmax){ tf=ff; i=d; } /* fin si parcouru vmax sommets */
      }else /* si v a déjà été visité (ou ne doit pas être visité) */
	if((X->D[u]<h)&&(v!=X->P[u])){ /* sinon X->cycle ne peut plus être améliorée */
	  h=X->D[u]+1; /* pas au delà de X->D[u] */
	  X->cycle=min(X->cycle,h+X->D[v]);
	}
    }
  }

  if(X->cycle>n) X->cycle=-1; /* si > n, alors pas trouvé de cycle -> -1 */
  X->n=ff; /* nb de sommets parcourus */
  X->tf=tf; /* nb de sommets parcourus dont les voisins ont été enfilés */
  X->radius=X->D[u]; /* hauteur du dernier sommet défilé (le plus éloigné) */
  
  /* c'est une mauvaise idée de faire ici un REALLOC(X->file,ff)
     car lors d'appels suivants avec la même structure, X->file n'aura
     plus la bonne taille ! */

  return X;
}


typedef struct{
  int nc; // nc=nombre de composantes connexes du graphe
  int na; // na=nombre de sommets d'articulation (ou cut-vertex)
  int *C; // C[u]=couleur de la composante de u, entier de [0,nc[ ou sommet supprimé
  int *P; // P[u]=parent de u dans une forêt couvrante, P[racine]=-1
  int *R; // R[i]=i-ème sommet racine dans la forêt couvrante, i dans [0,nc[
  int *d; // d[u]=date de début de visite du sommet u, entier de [0,n[
  int *A; // A[u]=vrai ssi u est un sommet d'articulation, u dans [0,n[
} param_dfs;


param_dfs *new_param_dfs(int n){
/*
  Attention ! le tableau X->C n'est pas alloué par cette fonction,
  même si n>0. Cela doit être fait par l'utilisateur.
*/

  NALLOC(param_dfs,X,1);

  X->C=X->P=X->R=X->d=X->A=NULL;
  X->nc=X->na=0;

  if(n>0){
    // X->C est alloué par l'utilisateur ou par dfs()
    ALLOC(X->P,n);
    ALLOC(X->R,n);
    ALLOC(X->A,n);
    ALLOC(X->d,n);
  }
  return X;
}


void free_param_dfs(param_dfs *X){
  if(X==NULL) return;
  free(X->C);
  free(X->P);
  free(X->R);
  free(X->A);
  free(X->d);
  free(X);
}


param_dfs *dfs(graph *G,int source,param_dfs *X)
/*
  Effectue un parcours en profondeur du graphe G depuis le sommet
  source. Version non récursive. On détermine également tous les
  sommets d'articulations (voir la définition de param_dfs pour lire
  le résultat). On l'utilise comme suit:

  param_dfs *p=dfs(G,s,NULL); // DFS dans G depuis s
  ...
  free_param_dfs(p);

  ou alors, pour un DFS dans G évitant les sommets de T:

  param_dfs *p=new_param_dfs(G->n);
  ALLOCZ(p->C,G->n,-1);
  for(i=0;i<G->n;i++) p->C[T[i]]=-2;
  dfs(G,s,p);
  ...
  free_param_dfs(p);

  Si p->C[u]=-2, alors le sommet u n'est pas parcouru (il est
  virtuellement supprimé de G). Les autres sommets v (non supprimés)
  doivent avoir p->C[v]=-1. Si p->C==NULL (par défaut), alors ce
  tableau est alloué et initialisé à -1. Il sera libéré par
  free_param_dfs(p).
*/
{
  if(G==NULL) return NULL;
  int u,i,d,v,k,t,n=G->n,r=0;
  int tete,nc,na,b;

  if(X==NULL){ r=1; X=new_param_dfs(n); }
  if(X->C==NULL) ALLOCZ(X->C,n,-1);
  for(i=0;i<n;i++) X->A[i]=0;

  nc=na=0;
  NALLOC(int,pile,n);  /* pile */
  NALLOC(int,next,n);  /* next[u]=prochain voisin de u à visiter */
  NALLOC(int,level,n); /* level[u]=... */
  t=tete=-1;

  for(k=0;k<n;k++,source=(source+1)%n)
    /* on essaye tous les sommets à partir de source */
    if(X->C[source]==-1){ /* si ==-2 ou >=0 alors on saute le sommet */
      pile[++tete]=source;
      next[source]=0; /* premier voisin à visiter */
      X->P[source]=-1;
      X->R[nc]=source;

      while(tete>=0){ /* tant que la pile n'est pas vide */
	u=pile[tete]; /* u=sommet courant */
	i=next[u]; /* indice du prochain voisin de u à visiter */
	if(i==0){
	  X->C[u]=nc; /* couleur de la composante courante */
	  level[u]=X->d[u]=++t; /* date de début de visite */
	}
	d=G->d[u]; /* degré de u */
	b=1; /* sentiennelle pour savoir si on a trouvé un v */
	while(i<d){ /* on cherche le prochain voisin v de u non visité */
	  v=G->L[u][i++]; /* pour tous les voisins v de u */
	  if(X->C[v]==-1){ /* si v n'a jamais été visité */
	    if((u==source)&&(t>X->d[u])&&(!X->A[u])) na++,X->A[u]=1; /* u=cut-vertex */
	    X->P[v]=u; /* père(v)=u */
	    pile[++tete]=v; /* on empile v */
	    next[v]=b=0; /* le prochain voisin de v à visiter */
	    next[u]=i; /* le prochain voisin de u à visiter */
	    break;
	  } else /* v existe et a déjà été visité */
	    if((X->C[v]>=0)&&(v!=X->P[u])) /* si (u,v) est un arc de retour */
	      level[u]=min(level[u],X->d[v]);
	}
	if(b){ --tete; /* il n'y a plus de voisin v: on dépile u pour toujours */
	  if((v=(X->P[u]))>=0){ /* si u n'est pas une racine, il a un père v */
	    level[v]=min(level[v],level[u]); /* met à jour level du père de u */
	    if((v!=source)&&(level[u]>=X->d[v])&&(!X->A[v])) na++,X->A[v]=1; /* v=cut-vertex */
	  }
	}
      } /* fin du while(i<d) */

      nc++; /* ici nc=nb de composantes visitées */
    }

  X->nc=nc;
  X->na=na;
  free(pile);
  free(next);
  free(level);

  /* on réduit le tableau ->R au minimum que si dfs() l'a
     alloué. Sinon, il ne faut pas toucher aux pointeurs de X */
  if(r) REALLOC(X->R,nc);

  return X;
}


typedef struct{
  int n;        // nb de sommets du graphe
  int source;   // source
  int *pere;    // tableau des parents: pere[u]
  double *dist; // tableau de distance: dist[u]
} param_bellman;


param_bellman *new_param_bellman(int n){
  NALLOC(param_bellman,p,1);
  p->n=n;
  p->source=0;
  if(n>0){
    ALLOC(p->pere,n);
    ALLOC(p->dist,n);
  }
  return p;
}


void free_param_bellman(param_bellman *p){
  if(p==NULL) return;
  free(p->pere);
  free(p->dist);
  free(p);
  return;
}


param_bellman *Bellman_Ford(graph *G,int source){
  /*
    Calcule les distances du sommet source à tous les autres de G
    selon l'algorithme de Bellman-Ford. On retourne un tableau ->pere
    et ->dist. Il permet d'avoir des poids négatifs sur les arcs. Il
    ne faut pas qu'il y ait de cycle absorbant.
    
    On utilise en interne 4 tableaux de taille n.

    A faire: gérer les sommets éteints.
  */
  
  if(G==NULL) return NULL;
  if(G->W==NULL) return NULL;

  int u,v,n=G->n;
  int i,i1,i2,d;
  int *tmp;
  param_bellman *p=new_param_bellman(n);
  double w;

  NALLOC(int,pile1,n);
  NALLOC(int,pile2,n);
  NALLOCZ(int,vec1,n,1); /* vec1[u]=0 ssi u est dans pile1 */
  NALLOCZ(int,vec2,n,1); /* idem mais pour vec2 */

  for(u=0;u<n;u++){
    p->dist[u]=DBL_MAX;
    p->pere[u]=-1;
  }

  p->dist[source]=0.0;
  p->pere[source]=p->source=source;

  i1=i2=0;
  pile1[i1++]=source; /* empile la source */
  vec1[source]=0; /* source est dans la pile1 */

  while(i1>0){
    while(i1>0){
      u=pile1[--i1]; /* dépile u */
      vec1[u]=1; /* u n'est plus dans la pile1 */
      d=G->d[u]; /* d=degré de u */
      for(i=0;i<d;i++){
	v=G->L[u][i];
	w=p->dist[u]+G->W[u][i];
	if(w<p->dist[v]){
	  p->pere[v]=u;
	  p->dist[v]=w;
	  if(vec2[v]){ /* si v pas dans P2, empiler v dans P2 */
	    pile2[i2++]=v;
	    vec2[v]=0;
	  }
	}
      }
    }
    SWAP(i1,i2,i);
    SWAP(pile1,pile2,tmp);
    SWAP(vec1,vec2,tmp);
  }

  free(pile1);
  free(pile2);
  free(vec1);
  free(vec2);
  return p;
}


/*
  Dijkstra(G,W,s):

  1. Init(G,s)
  2. S={}
  3. Q=V(G)
  4. while Q<>{}
  5.  do u=Extract-Min(Q)
  6.     S=S u {u}
  7.     for each v in Adj[u]: Relax(u,v,W)

  Init(G,s):
  1. for each vertex v in V(G):
  2.   do d[v]=infinity
  3.      pi[v]=NULL
  4. d[s]=0

  Relax(u,v,W):
  1. if d[v]>d[u]+W[u,v]
  2.   then d[v]=d[u]+W[u,v]
  3.        pi[v]=u

  Extract-Min(Q):
  extract vertex u with smallest d[v]
*/

int WeightGraph(graph *G){
  /*
    Calcule les longueurs (poids) des arêtes dans le cas d'un graphe
    géométrique (tout en créant le tableau les tableaux G->W). Il faut
    que les tableaux XPOS,YPOS existent. Retourne 1 ssi tout c'est
    bien passé.
  */
  if((G==NULL)||(G->xpos==NULL)) return 0;

  int u,v,i,d,n=G->n;
  double dx,dy;
  ALLOC(G->W,n);

  for(u=0;u<n;u++){
    d=G->d[u];
    ALLOC(G->W[u],d);
    for(i=0;i<d;i++){
      v=G->L[u][i];
      dx=G->xpos[u]-G->xpos[v];
      dy=G->ypos[u]-G->ypos[v];
      G->W[u][i]=sqrt(dx*dx+dy*dy);
    }
  }
  return 1;
}

/***********************************

       ISOMORPHISM, SUBGRAPH,
       MINOR, INDUCEDSUBGRAPH,
       PATHS, PS1, TREEWIDTH ...

***********************************/


int *Isomorphism(graph *G,graph *H){
/*
  Renvoie un tableau P <> NULL ssi G est isomorphe à H. Si tel est le
  cas, le tableau P indique le morphisme de H vers G. Après l'appel,
  le graphe H est modifié: ses listes sont triées si H->sort est faux
  (sauf si G=H - même pointeur). Le graphe G n'est par contre pas
  modifié. Dans H->int1 est retourné le nombre de tests effectués.
  Moins le graphe possède de symétrie, plus faible est le nombre de
  tests (et rapide est la fonction).

  On applique l'algorithme suivant. Pour chacun des deux graphes et
  chaque sommet u, on calcule son "profile" un tableau noté
  profile[u]: profile[u][i+2] = le nombre de sommets à distance i de
  u, en commençant à partir de i=1. Donc, profile[u][3] est le degré
  de u. Ceci est calculé par un simple BFS, les indices 0,1,2 étant
  réservés. On utilise profile[u][0] pour coder la taille du tableau
  profile[u], profile[u][1] pour stocker le nom de u (plutôt que le
  nombre de sommet à distance 0 qui est toujours 1) et profile[u][2]
  pour stocker la longueur du plus petit cycle passant par u. Cette
  dernière information étant calculée "gratuitement" par le BFS.

  On ordonne ensuite les sommets des graphes suivant les profiles des
  sommets avec qsort(). On ne renumérote pas les sommets dans le
  graphe, mais plutôt on donne un ordre: c'est possible avec qsort()
  car le nom original u est dans profile[u][1] le sommet u. Si bien
  que profile[i][1] sera le sommet u de rang i. On priviligie les
  profiles le grande taille (que l'on classe en premier) ce qui est
  plus discriminant. Les isomorphismes (=permutations) à tester ne
  concernent que les sommets de même profile. On construit les
  contraintes dans un tableau C, C[j] indiquant que les sommets de
  rang C[j-1] à C[j] (exclu) ont même profile, et sur lequel se base
  NextPermutation().

  Sur le même principe, on pourrait imaginer un profile plus complexe,
  comme suit: à chaque distance i et sommet u, on calcule le graphe
  G[u][i] induit par les sommets à distance i de u. On peut alors
  calculer le profile de chaque sommet de G[u][i] et ordonner les
  sommets selon celui-ci.
*/
  
  H->int1=0; /* par défaut, 0 tests */
  if((G==NULL)||(H==NULL)) return NULL;
  if(G->n!=H->n) return NULL;

  int *P; /* isomorphisme final */
  int n=G->n;

  if(G==H){ /* isomorphisme trivial si même emplacement mémoire */
    ALLOCZ(P,n,_i);
    return P;
  }

  param_bfs *param=new_param_bfs(); /* pour le BFS */
  int **profile,**profileG=NULL,**profileH;
  int *R,*C; /* permutation et contraintes (sur les rangs) */
  int u,v,t,r,i;
  graph *M;

  for(M=G;;){ /* on fait la même chose pour M=G puis M=H */
    ALLOC(profile,n); /* profile[u] */
    for(u=0;u<n;u++){ /* faire un BFS pour tous les sommets u */
      bfs(M,u,param); /* le premier BFS va allouer param->D et param->P */
      t=3+param->radius; /* taille du tableau profile[u] */
      ALLOCZ(profile[u],t,0); /* initialise le profile */
      for(v=0;v<n;v++){
	i=param->D[v];
	if(i>0) profile[u][i+2]++; /* compte les sommets à distance i>0 de v */
	param->D[v]=-1; /* réinitialise les distances pour les BFS suivants */
      }
      profile[u][0]=t; /* taille du profile */
      profile[u][1]=u; /* nom du sommet, pour qsort() */
      profile[u][2]=param->cycle; /* maille */
    }
    qsort(profile,n,sizeof(int*),fcmp_profile); /* trie les profiles */

    if(M==H){ profileH=profile; break; } /* on s'arête si M=H */
    profileG=profile; /* on refait la boucle pour H */
    M=H;
  }
  free_param_bfs(param);

  /* on verifie que profileG "=" profileH */
  for(u=0;u<n;u++)
    if(fcmp_profile(profileG+u,profileH+u)){
      P=NULL;
      goto fin_noniso;
    }

  /* calcule les contraintes */
  /* ici les profiles de G et H sont identiques */
  ALLOC(C,n);
  R=profile[0]; /* R=profile du premier sommet. NB: profile=profileH. */
  for(u=t=0;u<n;u++){
    if(fcmp_profile(&R,profile+u)){ /* si profiles différent */
      R=profile[u];
      C[t++]=u;
    }
  }
  C[t]=n;

  ALLOC(P,n);
  ALLOCZ(R,n,_i); /* initialise l'isomorphisme sur les rangs des profiles */
  if(!H->sort) SortGraph(H,0); /* on trie H pour la recherche dichotomique */
  H->int1=0; /* compte le nb de tests */

  /* vérifie, pour chaque permutation P, que P(G)=H */

  do{
    H->int1++;

    /* calcule P en fonction de R */
    for(r=0;r<n;r++) P[profileG[r][1]]=profileH[R[r]][1];
    for(r=0;r<n;r++){ /* on commence par les profiles longs: r=0,1,2, ... */
      u=profileG[r][1]; /* v=P[u] est le sommet correspondant à u dans H */
      for(i=0,t=G->d[u];i<t;i++) /* on regarde si chaque voisin de u est dans H->L[v] */
	if(bsearch(P+(G->L[u][i]),H->L[P[u]],t,sizeof(int),fcmp_int)==NULL){
	  /* alors élément non trouvé */
	  i=t;r=n; /* prochaine permutation à tester */
	}
    } /* ici r=n (trouvé) ou n+1 (pas encore trouvé) */
    if(r==n) goto fin_iso; /* on a trouvé un isomorphisme P */
  }
  while(!NextPermutation(R,n,C));

  /* si on arrive ici, c'est qu'on a pas trouvé l'isomorphisme */
  free(P);
  P=NULL;

 fin_iso:
  free(R);
  free(C);

 fin_noniso:
  FREE2(profileG,n);
  FREE2(profileH,n);

  return P;
}


edge *ListEdges(graph *G){
/*
  Construit la liste des arêtes de G, chaque arête uv ne figure qu'une
  seule fois. On a, pour tout i, E[i].u < E[i].v. Le champs G->m est
  mis à jour.
*/
  int u,v,i,j,d;
  int m=NbEdges(G);
  int n=G->n;

  NALLOC(edge,E,m);
  for(u=j=0;u<n;u++){
    for(i=0,d=G->d[u];i<d;i++){
      v=G->L[u][i];
      if(u<v){ E[j].u=u; E[j++].v=v; }
    }
  }
  return E;
}


graph *Subgraph(graph *G,graph *H){
  /*
    Détermine si G est un sous-graphe de H s'ils ont même nombre de
    sommets. On renvoie un sous-graphe S de H isomorphe à G, et on
    renvoie NULL si H ne possède pas G comme sous-graphe.

    Effets de bord:
    - les listes de G sont triées et donc G->sort=1,
    - H->int1 contient le nombre total de tests effectués,
    - S->pint1 contient l'isomorphisme de S vers G, et donc de H vers G.

    L'algorithme est le suivant: on teste d'abord si la séquence des
    degrés de H est bien supérieure à celle de G (ceci prend un temps
    O(n)). Si c'est le cas, on effectue, pour tous les sous-graphes S
    de H qui ont autant d'arêtes que G, un test d'isomorphisme entre S
    et G grâce à isomorphisme(S,G).
  */
  int n=H->n;
  H->int1=0; /* nb de tests = 0 par défaut */
  if(n!=G->n) return NULL; /* pas le même nombre de sommets */

  /* on trie en O(n) les deux listes */
  int *Eh=SortInt(H->d,NULL,n,0,NULL,SORT_INC);
  int *Eg=SortInt(G->d,NULL,n,0,NULL,SORT_INC);

  int i;
  /* on s'arrête si, pour un rang i donné, degH(i)<degG(i) */
  for(i=0;i<n;i++) if(Eh[i]<Eg[i]) break;
  free(Eh);
  free(Eg);
  if(i<n) return NULL; /* G ne peut pas être un sous-graphe de H */

  int mg=NbEdges(G);
  int mh=NbEdges(H);
  graph *S=new_graph(n); /* sous-graphe de H, alloue S->d et S->L */
  edge *E=ListEdges(H); /* liste des arêtes de H: e_j=(u,v) -> E[j].u=u et E[j].v=v */
  int *B; /* les arêtes de S, sous-ensemble d'indices d'arêtes de H */
  int *P; /* isomorphisme S->G */
  int u,v,j,d;

  /* réserve l'espace pour S, sous-graphe de H */
  for(i=0;i<n;i++) ALLOC(S->L[i],H->d[i]);

  /* initialise le sous-ensemble d'arêtes B de H avec mg arêtes */
  ALLOC(B,mg);
  NextSet(B,-1,mg);
  d=0; /* d=compteur pour le nombre total de tests */

  do{

    /* remplit S avec les arêtes de B */
    for(u=0;u<n;u++) S->d[u]=0; /* position libre pour le sommet u de S */
    for(i=0;i<mg;i++){ j=B[i]; /* j=numéro de la i-ème arête de B */
      u=E[j].u; v=E[j].v; /* l'arête j de E est (u,v) */
      S->L[u][S->d[u]++]=v; /* ajoute v comme voisin de u et incrémente deg(u) */
      S->L[v][S->d[v]++]=u; /* ajoute u comme voisin de v et incrémente deg(v) */
    }
    
    /* il vaut mieux que G soit le 2e paramètre, car il va être trié
       la première fois par Isomorphism(), puis plus jamais grâce au
       test de G->sort, alors que S serait trié à chaque appel */
    P=Isomorphism(S,G);
    d += 1+G->int1; /* on ajoute 1 pour donner le nombre d'ensembles testés */
  }while((P==NULL)&&(!NextSet(B,mh,mg)));

  H->int1=d; /* nombre total de tests */

  if(P==NULL){ /* on a pas trouvé de sous-graphe de H isomorphe à G */
    free_graph(S);
    S=NULL;
  }
  else S->pint1=P; /* S isomorphe à G, sous-graphe de H */

  free(E);
  free(B);
  return S;
}


graph *MatrixToGraph(int **M,int n){
  /*
    Renvoie le graphe correspondant à la matrice d'adjacence n x n
    symétrique où seule la partie inférieure est utilisée. Les listes
    du graphe de retour sont triées et les champs ->m et ->sort sont
    mise à jour.
  */
  if((M==NULL)||(n<=0)) return NULL;
  int u,v,m;
  graph *G=new_fullgraph(n);

  for(u=m=0;u<n;u++)
    for(v=0;v<u;v++)
      if(M[u][v]){
	ADD_EDGE(G,u,v); /* ajoute u-v et v-u */
	m++; /* une arête de plus */
      }
  
  /* réduit les listes */
  GraphRealloc(G,G->d);

  G->m=m;
  G->sort=1;
  return G;
}


graph *GraphOfColor(graph *G,int *col,int k){
  /*
    Renvoie un graphe C dont les sommets sont les entiers de [0,k[
    (les valeurs du tableau col[] qui doit avoir une taille G->n) et
    les arêtes les paires uv telle qu'il existe une arête xy de G avec
    col[x]=u et col[y]=v. La valeur C->m est mise à jour, et les
    listes de C sont triées (C->sort=1).
  */
  if((k<0)||(col==NULL)||(G==NULL)) return NULL;

  int u,v,cu,cv,i,d,n=G->n;
  graph *C; /* le graphe des couleurs renvoyé */

  /* matrice d'adjacence inférieure à 0 */
  NALLOCMAT(int,M,k,k-1); /* matrice d'adjacence du graphe des couleurs */
  for(u=0;u<k;u++)
    for(v=0;v<u;M[u][v++]=0);

  for(u=0;u<n;u++) /* parcourt G et remplit la partie inférieure de M */
    for(i=0,d=G->d[u];i<d;i++){
      v=G->L[u][i];
      if(u<v){ /* si cu=cv, on ne fait rien */
	cu=col[u];
	cv=col[v];
	if(cu>cv) M[cu][cv]=1;
	if(cv>cu) M[cv][cu]=1;
      }
    }
  
  C=MatrixToGraph(M,k);
  FREE2(M,k);
  return C;
}


int FindSet(int x,int *p)
/* routine pour UNION-FIND pour Minor() */
{
  if(x!=p[x]) p[x]=FindSet(p[x],p);
  return p[x];
}


int *Minor(graph *H,graph *G)
/*
  Détermine si H est un mineur de G. Si c'est le cas, un tableau T est
  renvoyé, sinon NULL est renvoyé. Le tableau T code un modèle du
  mineur de H dans G. Plus précisément, T[u]=c, où u est un sommet de
  G, est le nom du sommet c de H si bien que l'ensemble des sommets u
  de G tel que T[u]=c forme le super-noeud c.

  L'algorithme est le suivant: on effectue, pour tous les ensembles de
  contractions d'arêtes de G produisant un mineur avec autant de
  sommets que H, un test de sous-graphe grâce à Subgraph().
*/
{
  graph *C; /* graphe des couleurs = graphe contracté */
  graph *S=NULL; /* sous-graphe éventuellement isomorphe à C */
  int *B; /* sous-ensemble (d'indices) d'arêtes de G à contracter */
  int *couleur; /* les sommets de même couleurs seront contractés */
  int *rang; /* pour union-find rapide */
  edge e;
  
  int nh=H->n;
  int ng=G->n;
  int c=ng-nh; /* c=nb de contractions à effectuer */
  int t;
  H->int1=t=0; /* initialise le nb de tests */
  if(c<0) return NULL; /* pas de mineur */

  edge *E=ListEdges(G); /* E=liste des arêtes de G, met à jour G->m */
  int mg=G->m;
  if(c>mg){ /* fini s'il faut contracter plus d'arêtes qu'il y a dans G */
    free(E);
    return NULL;
  }

  int i,u,v,x,y;
  int test=((c<<1)<ng); /* vrai ssi on liste suivant les arêtes ou suivant les sommets */
  ALLOC(B,c); NextSet(B,-1,c); /* B=premier sous-ensemble de c arêtes de E */
  ALLOC(couleur,ng); /* couleur des sommets */
  ALLOC(rang,ng); /* rang des sommets, sert pour UNION-FIND */

  /*
    On pourrait générer des sous-ensembles acycliques d'arêtes
    directement, en combinant NextSet() et le test d'acyclicité.
   */

  do{

    t++; /* on compte 1 test par ensemble B testés */

    /* initialise la couleur et le rang des sommets */
    for(u=0;u<ng;u++){ couleur[u]=u; rang[u]=0; }

    /* on teste rapidement (avec UNION-FIND) si B contient un cycle */
    for(i=0;i<c;i++){ e=E[B[i]]; /* e=i-ème arête de B */
      u=e.u; couleur[u]=x=FindSet(u,couleur);
      v=e.v; couleur[v]=y=FindSet(v,couleur);
      if(x==y) break; /* y'a un cycle, sinon on fait UNION */
      if(rang[x]>rang[y]) couleur[y]=x;
      else{ couleur[x]=y; if(rang[x]==rang[y]) rang[y]++; }
    }

    if(i==c){ /* si B est acyclique, on fait un test de sous-graphe */
      if(test)
	/* on met à jour la couleur de chaque sommet. Suivant les
	   valeurs respectives de c et ng (test) on met à jour soit
	   suivant les arêtes ou suivant les sommets. */
	for(i=0;i<c;i++){
	  e=E[B[i]]; /* e=i-ème arête de B */
	  u=e.u; couleur[u]=FindSet(u,couleur);
	  v=e.v; couleur[v]=FindSet(v,couleur);
	}
      else
	for(u=0;u<ng;u++) couleur[u]=FindSet(u,couleur);

      /* on recadre les couleurs dans [0,c[. Complexité: 4ng */
      for(i=0;i<ng;rang[i++]=0);
      for(i=0;i<ng;i++) rang[couleur[i]]++; /* rang=fréquence */
      for(i=u=0;i<ng;i++) /* repère les couleurs manquantes */ 
	if(rang[i]==0) u++; else rang[i]=u;
      /* ici rang[i]=nb de zéros (=valeurs manquantes) dans rang[0..i[ */
      for(i=0;i<ng;i++) couleur[i] -= rang[couleur[i]];

      C=GraphOfColor(G,couleur,nh);
      S=Subgraph(H,C); /* avant l'appel, S=NULL nécessairement */
      t += C->int1;
      free_graph(C); /* on a plus besoin du graphe des couleurs */
    }

  }while((S==NULL)&&(!NextSet(B,mg,c)));

  H->int1=t;
  free(B);
  free(E);

  /* on a rien trouvé */
  if(S==NULL){
    free(couleur);
    free(rang);
    return NULL;
  }
  
  /* on a trouvé un mineur, on construit le modèle dans rang[] */
  for(u=0;u<ng;u++) rang[u]=S->pint1[couleur[u]];
  free_graph(S);
  free(couleur);
  return rang;
}


int *InducedSubgraph(graph *H,graph *G)
/*
  Indique si H est un sous-graphe induit de G.  La fonction renvoie un
  ensemble X de sommets tel que G[X] est ismomorphe à H. Evidemment
  |X|=H->n. On renvoie dans G->int1 le nombre de tests effectués, et
  dans G->pint1 l'isomorphisme entre H et G[X].

  L'algorithme consiste à générer tous les ensembles X possibles de
  |V(H)| sommets et à tester l'isomorphisme entre G[X] et H.
 */
{
  if((G==NULL)||(H==NULL)) return NULL;
  int ng=G->n,nh=H->n;
  if(nh>ng) return NULL;

  graph *S;
  int *P,t=0;
  NALLOC(int,X,nh);
  NextSet(X,-1,nh); /* premier sous-ensemble */

  do{
    t++;
    S=ExtractSubgraph(G,X,nh,1);
    P=Isomorphism(S,H);
    t += H->int1;
    free_graph(S);
  }while((P==NULL)&&(!NextSet(X,ng,nh)));

  G->int1=t;
  free(G->pint1); /* pour éviter les fuites mémoires */
  G->pint1=P;
  if(P==NULL){ free(X); X=NULL; }
  return X;
}


int NextPath(graph *G,path *P,int j)
/*
  Cette fonction (récursive) permet de trouver tous les chemins
  simples entre deux sommets. Si 1 est renvoyé, c'est que tout c'est
  bien passé, sinon on renvoie 0. On l'utilise comme ceci:

   path *P=new_path(G,NULL,G->n); // crée un chemin vide P mais avec G->n sommets possibles
   P->P[0]=x; // origine du chemin
   P->P[1]=y; // destination du chemin, y=x possible
   r=NextPath(G,P,-1); // crée le premier chemin, r=0 s'il n'existe pas, sinon r=1
   r=NextPath(G,P,0);  // calcule le chemin suivant, r=0 s'il n'existe pas, sinon r=1
   r=NextPath(G,P,0);  // calcule le chemin suivant, r=0 s'il n'existe pas, sinon r=1
   ...
   free_path(P); // libère le chemin P

  Plus précisément, étant donnés un chemin P=x-...-v_j-...-y du graphe
  G et un sommet v_j du chemin (v_j=j-ème sommet du chemin), la
  fonction tente de compléter P par un autre chemin de v_j à y évitant
  x-...-v_(j-1). Si ce chemin à été trouvé, alors la partie de v_j à y
  de P est mise à jour et on renvoie 1. Sinon, le chemin est coupé
  après v_j et on renvoie 0. Dans tous les cas P est un chemin à jour
  de G. Si j<0, alors on initialise P par un chemin allant de
  x=P->P[0] à y=P->P[1].

  Algorithme: on essaye d'améliorer en premier le sous-chemin de
  v_(j+1) à y (récursivement). Si cela n'est pas possible, on calcule
  un nouveau chemin de v_j à y passant par un autre voisin v de v_j
  (autre que v_(j+1)) et évitant le chemin x-...-v_j. On passe en
  revue ainsi tous les voisins v de v_j qui ne sont pas dans
  P[0]...P[j]. Si aucun des voisins n'a un tel chemin, on retourne 0.

  Comme il faut tester les voisins v de v_j qu'une seule fois (pour un
  chemin P[0]...P[j] donné), on utilise le champs aux[v_j][i] de P qui
  donne le i-ème voisin libre de v_j avec la convention que
  aux[v_j][0] est le nombre de voisins encore possibles.

  Effet de bord: P est mis à jour.
*/
{
  if((P==NULL)||(G==NULL)) Erreur(-1); /* ne devrait jamais arriver */
  param_bfs *p;
  int i,x,y,u,v,n;

  if(j<0){ /* initialisation du premier chemin */
    n=G->n;
    if(P->aux==NULL) ALLOCMAT(P->aux,n,n);
    for(u=0;u<n;P->V[u++]=-1); /* vide le chemin */
    x=P->P[0];
    y=P->P[1];
    if((x<0)||(y<0)||(x>=n)||(y>=n)){ P->n=0; p=NULL; goto fin_0; } /* sommets inexistant */
    p=bfs(G,x,NULL); /* calcule le chemin */
    i=p->D[y];
    if(i<0){
      P->V[x]=0; /* x est en première position dans P */
      P->n=1;
      goto fin_0;
    }
    /* on initialise aux[x] et aux[y] qu'une seule fois */
    P->aux[y][0]=0;
    j=-1; /* pour que la longueur soit correcte */
    goto fin_1;
  }

  n=P->n;

  if(j+1==n) return 0; /* si x=y, alors pas de prochain chemin */
  if(NextPath(G,P,j+1)) return 1; /* c'est le NextPath à partir du voisin de v_j */

  /* Ici on ne peut pas augmenter v_(j+1)...y. Donc, il faut calculer
     un premier chemin depuis le prochain voisin v disponible de u=v_j
     et y qui évite x=P[0]-...-P[j]. Si cela ne marche pas avec le
     voisin v, il faut essayer le suivant, etc. */

  /* efface depuis P la fin du chemin P[j+1]...P[n-1] */
  /* pour ne garder que P[0]...P[j] */
  for(i=j+1;i<n;i++) P->V[P->P[i]]=-1;

  /* rem: ici t<d */
  p=new_param_bfs();
  ALLOC(p->D,G->n);
  y=P->P[n-1];
  u=P->P[j];
  i=-1;

  while((P->aux[u][0])&&(i<0)){ /* tant qu'il y a encore des voisins de u non testés */
    v=P->aux[u][(P->aux[u][0])--]; /* lit et enlève le dernier voisin dispo */
    if(P->V[v]>=0) continue; /* si le voisin est dans P[0] ... P[j] on le saute */
  
    /* initialise p->D: on doit le faire à chaque appel */
    for(i=0;i<G->n;i++) p->D[i]=-1;
    /* on enlève P[0]...P[j] du graphe pour le BFS */
    for(i=0;i<=j;i++) p->D[P->P[i]]=-2;

    /* calcule un chemin de v à y dans G\(P[0]-...-P[j]) */
    bfs(G,v,p);

    /* a-t-on trouvé un chemin ? */
    i=p->D[y]; /* si i>=0, c'est oui */
  }

  /* ici i=distance de u=v_j à y */
  if(i<0){ /* on a pas trouvé de chemin */
  fin_0:
    free_param_bfs(p);
    return 0;
  }

  /* on a trouvé un chemin, on met à jour P */
 fin_1:
  P->n=i+j+2; /* nouvelle longueur */
  /* ajoute à P[0]...P[j] le chemin de P[j+1] à y en partant de y */
  while(i>=0){
    P->P[i+j+1]=y;
    P->V[y]=i+j+1;
    y=p->P[y];
    i--;
  }

  /* initialise aux[u] pour tous les sommets u de P[j+1] à y non
     compris. Pour P[j] on le fait au moment de la lecture de v dans
     aux[u], et pour y on le fait une seule fois à la création
     (avec aux[y][0]=0) */

  for(i=j+1,n=P->n-1;i<n;i++){ /* si j=-1, alors on fait pour x aussi */
    u=P->P[i];
    P->aux[u][0]=0;
    for(v=0;v<G->d[u];v++){
      j=G->L[u][v]; /* j=voisin de u */
      /* si pas dans P on l'ajoute comme voisin possible */
      if(P->V[j]<0) P->aux[u][++P->aux[u][0]]=j;
    }
  }

  free_param_bfs(p);
  return 1;
}


int NextPath2(graph *G,path *P,int j)
/*
  Comme NextPath(G,P,j) sauf que si le nouveau chemin P' renvoyé a le
  même ensemble de sommets sur les positions allant de j à n-1, on
  passe au chemin suivant. Donc on renvoie le prochain chemin de P
  (s'il existe) qui a un ensemble de sommets différents de P->P[j] à
  P->P[n-1].
*/
{
  if(j<0) return NextPath(G,P,j);
  if(P==NULL) Erreur(-1);

  int i,n=P->n,m=n-j;
  NALLOC(int,T,m);
  for(i=j;i<n;i++) T[i-j]=P->P[i]; /* copie les n-j derniers sommets de P */

  if(!NextPath(G,P,j)){
    free(T);
    return 0;
  }

  if(n==P->n){
    for(i=0;i<m;i++)
      if(P->V[T[i]]==-1){ /* ici on a un chemin différent */
	free(T);
	return 1;
      } /* ici on a le même chemin qu'au départ */
    free(T);
    return NextPath2(G,P,j);
  }

  free(T); /* ici on a un chemin différent */
  return 1;
}


#define LCONF 9  /* nb de bits pour coder un sommet dans le graphe des conflits */
#define CONFMAX (1<<LCONF) /* = 2^LCONF = nb max de sommets du graphe des conflits */
#define CONFMASK (CONFMAX-1) /* = CONFMAX-1 = mask pour un sommet */
#define CONFC 2 /* constante pour modifier le code des noeuds à 0. Il faut CONFC>1 */
#define NCMAX 512 /* taille maximum de la liste de composantes maximales. NCMAX<CONFMAX */

/* ensemble de variables pour gérer le graphe des conflits */
typedef struct{
  int n;   /* nb de sommets du graphe d'origine, sert pour la règle 5 */
  int nbi; /* nb de noeuds avec code=-1 (indéterminée) */
  int nbzi;/* nb de noeuds de code 0 indépendants */
  int outmem; /* 1 ssi il y a eut un dépassement de CONFMAX ou NCMAX */
  int paire[CONFMAX]; /* 1er noeud de la paire xy */
  int path[CONFMAX];  /* 1er noeud du chemin P, ne sert que pour PrintConflit() */
  int code[CONFMAX];  /* code du noeud i. Valeurs possibles: -1,0,1 */
  int nbc[CONFMAX];   /* nb de noeuds de la paire i dont le code est <1 */
  int *comp[CONFMAX];  /* liste des noeuds de la composante de i */
  int *cmax[NCMAX+1]; /* liste des composantes maximales */
  int tmax[NCMAX+1]; /* taille des composantes maximales */
  int ncmax; /* taille de la liste cmax */
  int nodemax; /* noeud correspondant à la composante maximale de la liste cmax */
  int valmax; /* valeur (VRAI/FAUX) des composantes maximales */
  int tcomp[CONFMAX]; /* taille de comp[i] */
  int x[CONFMAX],y[CONFMAX]; /* paire du noeud i, ne sert que pour PrintConflit() */
  graph *G; /* graphe des conflits */
  /* Attention! les noeuds des listes d'adjacence de G peuvent être >
     CONFMAX et donc aussi > G->n. On s'en sert pour coder un type
     d'arête. Donc si u est un noeud de G et v=G->L[u][i], alors le
     i-ème voisin de u est (v&CONFMASK) et le type de l'arête uv est
     (v>>LCONF). */
} conflit;


void PrintConflit(conflit *c)
/*
  Affiche le graphe des conflits, et rien si c->G->n=0.
*/
{
  if(c->G->n==0) return;

  int u,v,i,t,d;
  char T1[]="|=><";
  char T2[]="-01*";
  char T3[]="_/.";

  // UTF-8:
  // char T3[]="━┛o";
  // \xe2\x94\x80: ─
  // \xe2\x94\x81: ━
  // \xe2\x94\x98: ┘
  // \xe2\x94\x99: ┙
  // \xe2\x94\x9a: ┚
  // \xe2\x94\x9b: ┛

  printf("code (x,y) [comp] node-path-paire: node[type]\n");
  for(u=0;u<c->G->n;u++){
    t=c->code[u];
    if(c->nbzi>0){ /* si on a des zéros indépendants */
      if(t==0) t=2; /* c'est un zéro indépendant */
      else if(t==CONFC) t=0; /* c'est un zéro marqué -> remet sa valeur d'origine */
    }
    printf(" %c ",T2[t+1]);
    if(c->paire[u]==u)
      printf("(%i,%i) [",c->x[u],c->y[u]);
    else printf("\t [");
    for(v=0;v<c->tcomp[u];v++){
      printf("%i",c->comp[u][v]);
      if((c->n>9)&&(v<c->tcomp[u]-1)) printf(",");
    }
    printf("]\t%s%i",(v<5)?"\t":"",u);
    if(u>c->path[u]) printf("%c \t:",T3[1]); else{
      printf("%c",T3[0]);
      if(u>c->paire[u]) printf("%c\t:",T3[1]); else{
	printf("%c%c\t:",T3[0],T3[2]);
      }
    }
    for(i=0,d=min(c->G->d[u],WIDTH);i<d;i++){
      v=c->G->L[u][i];
      t=(v>>LCONF);
      printf(" %i%c",v&CONFMASK,T1[t]);
    }
    if(c->G->d[u]>WIDTH) printf("...");
    printf("\n");
  }
  printf("#nodes in conflict graph: %i\n",c->G->n);
  printf("#heavy indep. components: %i\n",c->nbzi);
  printf("#unspecified values: %i\n",c->nbi);
  if(c->outmem){
    printf("!!! Out of Memory !!!\n");
    if(c->ncmax>=NCMAX)
      printf("#nodes in conflit graph exceeded (%i)\n",CONFMAX);
    else
      printf("#maximal components exceeded (%i)\n",NCMAX);
  }
  return;
}


void ps1_delxy(conflit *c,int w)
/*
  Supprime du graphe des conflits c la dernière paire créee, et met
  ainsi à jour c. Ici w est l'indice du premier noeud de la dernière
  paire dans le graphe des conflits. Il faut supprimer tous les noeuds
  u >= w et tous les arcs entre les noeuds u et v<w.

  La liste des composantes maximales de la dernière paire sera effacée
  à la fin de la paire courante (voir le code après nextxy:).
*/
{
  int u,i,v,j;

  for(u=w;u<c->G->n;u++){
    for(i=0;i<c->G->d[u];i++){ /* pour tous les arcs v->u */
      v=c->G->L[u][i]&CONFMASK; /* v est le i-ème voisin de u */
      if(v<w) /* seulement pour les noeuds v avant w */
	for(j=0;j<c->G->d[v];j++)
	  /* on cherche et supprime u de la liste de v */
	  if((c->G->L[v][j]&CONFMASK)==u)
	    c->G->L[v][j]=c->G->L[v][--(c->G->d[v])]; /* supprime u */
    }
    free(c->G->L[u]); /* supprime la liste de u */
    free(c->comp[u]); /* supprime la composante de u */
  }

  c->G->n=w; /* met à jour le nombre de noeuds du graphe des conflits */
  return;
}


int ps1_addmax(int *C,int t,conflit *c)
/*
  Ajoute une composante C (de taille t) à la liste des composantes
  maximales déjà rencontrées, et maintient cette propriété. Il faut
  que la liste (c->cmax) soit de taille suffisante pour accueillir une
  nouvelle composante. La fonction renvoie VRAI ssi la composante C a
  été ajoutée à la liste.

  Algorithme:
    1. Pour toute composante X de L (=la liste):
       1.1. Si C est inclue ou égale à X alors FIN
       1.2. Si C contient X alors supprimer X de L
    2. Ajouter C à L

  On remarque aussi qu'il n'est pas possible de rencontrer plusieurs
  fois la même composante C avec des valeurs différentes. Si cela
  arrive, c'est pour deux chemins différents, disons Q1 et Q2. Soit X
  les voisins de C qui sont à la fois dans Q1 et Q2.  X se décompose
  en segments, chacun étant un sous-chemin de Q1 inter Q2. Tout sommet
  voisin de C doit être dans X sinon C aurait un sommet de trop. Donc
  tout chemin de G contenant X tel que C est une composante de G\P,
  doit être parallèle à Q1 et Q2. En particulier Q1 et Q2 sont
  parallèles ... et alors [A FINIR] ? Bon, bah en fait c'est
  possible. Q1 peut se réduire à une arête xy et Q2 à xzy (un sommet
  de plus). Et avec Q1 c'est vrai, mais avec Q2 cela devient faux. On
  a des contre-exemple avec des sommets (z) de degré deux.
*/
{
  int **L=c->cmax,n=c->ncmax; /* L=liste et n=|L|*/
  int *T=c->tmax; /* T=taille de chaque composante */
  int i,r;

  /* on passe en revue chaque composante la liste L */
  for(i=0;i<n;i++){
    r=SetCmp(C,L[i],t,T[i]); /* compare C et L[i] */
    if(r&6) return 0; /* C est strictement contenu (&4) ou égale (&2) à L[i] */
    if((r&8)&&(n>0)){ /* L[i] est contenu (strictement) dans C et L non vide */
      free(L[i]); /* libère ce tableau de sommets précédemment alloué par un ps1_addmax() */
      /* si L[i] dernier élément de L (i=n-1), alors il n'y a plus rien à faire */
      n--;
      if(i<n){ /* si c->cmax[i] pas dernier élément, il faut déplacer le dernier */
	L[i]=L[n]; L[n]=NULL; /* pour éviter plus tard un double free. NB: i<>n */ 
	T[i]=T[n]; /* nombre de sommets dans L[i] */
	i--; /* pour recommencer avec ce nouveau L[i] */
      }
    }
  }

  /*ajoute C à la fin de L */
  ALLOCZ(L[n],t,C[_i]); /* alloue et copie C dans L[i] */
  T[n]=t; /* taille de C */
  c->ncmax=n+1; /* nouvelle taille de L */
  return 1; /* on a ajouté C à L */
}


int ps1_push(int x,int v,conflit *c)
/*
  Affecte la valeur v (=0 ou 1) dans le noeud x et propage, en
  appliquant les règles décrites ci-après, à tous ses voisins et
  récursivement à tout le graphe. Renvoie 1 si une contradiction
  apparaît, ou bien si tous les noeuds de la paire de x ont pour code
  1. Sinon, on renvoie 0 (terminaison standard).

  Attention! Il faut appliquer cette fonction que si la paire de x est
  complète. Donc il ne faut pas l'appliquer si on vient de créer le
  noeud x, car on n'est pas sûr de ne pas faire un "Out of Memory" un
  peu plus tard sur la même paire.

  Effet de bord: met à jour plusieurs champs de la variable c, mais ne
  modifie pas le graphe c->G.

  Règles pour l'arc x->y:

  x=noeud x
  y=noeud y
  v=valeur 0 ou 1 qui vient d'être écrite dans x
  t=type de l'arête: 0=(Tx|Ty), 1=(|Tx|=|Ty|), 2=(Tx<Ty), 3=(Tx>Ty)
  c=graphe des conflits

  (Attention! "Tx|Ty" signifie que les composantes Tx et Ty sont
  disjointes.)

  règle 1: si Tx|Ty (disjoint) et v=0, alors écrire 1 dans y
  règle 2: si Tx=Ty, alors écrire v dans y
  règle 3: si Tx<Ty et v=0, alors écrire 0 dans y
  règle 4: si Tx>Ty et v=1, alors écrire 1 dans y
  règle 5: si Tx|Ty, v=1 et |Tx|+|Ty|=n, alors écrire 0 dans y

  On applique également la "règle du dernier -1", à savoir si la paire
  de x, après avoir écrit v, possède exactement 1 seul noeud de valeur
  -1, alors on écrit 0 dans ce noeud.

  La "règle du max" et la "règle de l'influence des voisins" sont
  appliquées plus directement par la fonction ps1() lors de la
  création d'une paire. Elles n'ont pas lieu d'être lors de la
  propagation.
*/
{
  int i,d,y,t;

  if(c->code[x]==v) return 0; /* rien à faire */
  if(c->code[x]==1-v) return 1; /* contradiction */
  /* ici on a code[x]==-1 */
  c->code[x]=v; /* écrit la valeur v */
  c->nbi--; /* et une valeur indéterminée en moins ! */
  t=(c->nbc[c->paire[x]] -= v); /* diminue seulement si on a écrit 1 */
  if(t==0) return 1; /* la paire de x est bonne, elle contient que des 1 ! */

  /* applique les règles 1-5 à tous les arcs sortant de x */
  for(i=0,d=c->G->d[x];i<d;i++){ /* pour tous les voisins de x */
    y=c->G->L[x][i]; /* y=i-ème voisin de x */
    t=(y>>LCONF); /* t=type de l'arc x->y: 0,1,2,3 */ 
    y &= CONFMASK; /* y=numéro du noeud voisin de x */

    /* applique les règles */
    switch(t){
    case 0:
      if((v==0)&&(ps1_push(y,1,c))) return 1; /* règle 1 */
      if((v==1)&&(c->tcomp[x]+c->tcomp[y]==c->n)&&(ps1_push(y,0,c))) return 1; /* règle 5 */
      break;
    case 1: if(ps1_push(y,v,c)) return 1; break; /* règle 2 */
    case 2: if((v==0)&&(ps1_push(y,0,c))) return 1; break; /* règle 3 */
    case 3: if((v==1)&&(ps1_push(y,1,c))) return 1; break; /* règle 4 */
    }
  }

  /* règle du dernier -1 ? */
  if(c->nbc[c->paire[x]]==1){
    /* on cherche alors l'unique noeud x de la paire courante qui est
       de code < 1.  NB: ce noeud existe forcément */
    x=c->paire[x]; /* x=premier noeud de la paire de x */
    while(c->code[x]==1) x++; /* passe au suivant si le code est 1 */
    return ps1_push(x,0,c); /* écrit 0 dans le noeud x */
  }

  return 0;
}

/* pour le débugage de PS1() */
/* N=niveau de récursion, POS=numéro de ligne */
#define PRINTS do{				\
    int _i;					\
    printf("%03i:%02i  ",++POS,N);		\
    for(_i=0;_i<3*N;_i++) printf(" ");		\
  }while(0)

int PS1(graph *G,path *P,int version)
/*
  P est un chemin d'un graphe G, et G\P est formé d'une seule
  composante connexe. Il faut G et P <> NULL (initialisés avec
  new_graph() et new_path()). Renvoie 1 si P peut "séparer" G (voir
  explications ci-dessous). Il y a une optimisation avec le graphe des
  conflits si version>0. On utilise cette fonction comme ceci:

  path *P=new_path(G,NULL,G->n); // P=chemin vide, sans sommet
  int r=PS1(G,P); // r=0 ou 1
  free_path(P);
  
  Effet de bord: G->int1 retourne le nombre de tests (nombre de
  paires, nombre de chemins testés, et nombre de passes dans le graphe
  des conflits). Les autres champs de G et de P ne sont pas modifiés.

  Le paramètre "version" indique la variante du test:
  - version=0: sans le graphe des conflits
  - version=1: avec le graphe des conflits
  - version=2: comme version=1 mais sans le graphe des conflits lors de la récursivité
  - version=3: comme version=1 mais avec l'écriture de valeurs dans le graphe des conflits

  Améliorations possibles:

  - Si G n'est pas biconnexe, alors on pourrait tester si toutes ses
    composantes biconnexes sont bien évaluée à vraie. Si P
    n'intersecte pas une composante biconnexe B de G, alors il faut
    évaluer PS1(B,{}).

  - G\P étant connexe, on pourrait déjà supprimer les sommets de P qui
    forment un segment de P contenant une extrémité de P et qui n'ont
    pas de voisin dans G\P. En particulier, si une des extrémités de P
    est de degré 1, on peut la supprimer.

  - Privilégier les paires de sommets qui ne sont pas adjacents (pour
    diminuer la taille les composantes et avoir des chemins plus
    longs). Plus généralement, on pourrait privilégier les paires de
    sommets les plus distants. Ceci dit, on ne gagne probablement pas
    grand chose, et une renumérotation aléatoire des sommets devrait
    suffir pour ne pas traité les paires de sommets voisins en
    priorité.

  - On pourrait tester des cas simples pour G: arbre (tester si m=n-1,
    on sait que G est connexe), clique (tester si m=n(n-1)/2: si n<=4
    sommets alors vraie, sinon faux). (Ces tests sont déjà
    implémentés). Plus dur: G est outerplanar. En fait, si G est un
    arbre de cycles (chaque composante connexe est un cycle ou un K4),
    alors le test est vrai. C'est donc vrai en particulier si
    m=n. Pour tester si G est un arbre de cycle, il suffit de faire un
    DFS, puis de vérifier que si (u,x) et (u,y) sont deux arêtes qui
    ne sont pas dans l'arbre du DFS, alors x et y ne sont pas ancêtres
    l'un de l'autre (??? Pourquoi ???).

  - Pour tous les chemins possibles testés récursivement pour G, ne
    tester effectivement que ceux qui correspondent à des
    sous-ensembles de sommets différents puisque le résultat sera le
    même (mêmes composantes connexes). Pour cela, il faut gérer une
    table pour mémoriser les sous-ensembles testés. Notons que si un
    chemin est induit alors cela ne sert à rien de le mémoriser, il ne
    pourra jamais être rencontré de nouveau. On pourrait coder un
    sous-ensemble de n sommets par un entier sur n bits (n < 32 ou 64
    donc). La recherche/insertion pourrait être une recherche dans un
    arbre binaire de recherche.

 EXPLICATIONS:

  Soit P et Q deux chemins d'un graphe G. On dit que Q est parallèle à
  P, noté Q//P, s'il n'y a pas d'arête entre P\Q et G\(Q u P).

  Quelques propriétés:

  - La relation // n'est pas symétrique.
  - Tout chemin est parallèle au chemin vide.
  - Tout chemin est parallèle à lui-même.
  - La relation // n'est pas transitive.

  Pour la dernière proposition, on peut vérifier en fait que si R//Q
  et Q//P, alors on a R//P sauf s'il existe une arête uv telle que u
  in (P inter Q)\R et v in Q\(R u P).

  Soit P et Q deux chemins de G. On dit que Q dérive de P, noté Q///P,
  s'il existe une suite P_0,P_1,...,P_n de chemins de G avec P_0=Q et
  P_n=P tels que P_{i-1}//P_i pour chaque i=1..n.

  Quelques propriétés:

  - Si Q//P, alors Q///P. En particulier, Q///{} puisque Q//{}.
  - On peut avoir R//Q//P, soit R///P, sans avoir R//P (cf. ci-dessus).
  - Si R///Q et Q///P, alors R///P.

  On dit que Q///P dans un graphe valué (G,w) si tous les chemins
  P_0,...,P_n (en particulier Q et P) sont des plus courts chemins
  selon w.

  Soit P un chemin de G et w une valuation de G. On définit le
  potentiel pour P selon w comme score(P,w) := max_C { w(C)*|V(G)|+|C|
  } où le maximum est pris sur toute composante connexe C de G\P.

  Lemme 1. Supposons que G\P a une seule composante connexe, et soit Q
  un chemin de G parallèle à P différent de P. Alors, pour chaque
  valuation w de G, soit P est un demi-séparateur de G ou bien
  score(Q,w) < score(P,w).

  Preuve. Supposons que P ne soit pas un demi-séparateur de G pour la
  valuation w. Soit C la composante de G\P, et posons n = |V(G)|. Par
  définition, score(P,w) = w(C)*n + |C|. On a w(C) > w(G)/2, et donc
  w(P) < w(G)/2. Il suit que w(P) < w(C). Soit C' une composante de
  G\Q telle que w(C')*n+|C'| = score(Q,w). Comme Q est parallèle à P,
  soit C' est contenue dans C soit C' est contenue dans P. En effet,
  si C' intersecte C et P, alors C' contient une arête uv avec u in C
  et v in P. Bien sûr uv not in Q. Cela contredit le fait qu'il existe
  pas d'arête entre P\Q et G\(QuP).

  Si C' est contenue dans P, alors w(C') <= w(P) < w(C). Il suit que
  w(C') <= w(C)-1, soit w(C')*n <= w(C)*n - n. Clairement |C'| < n +
  |C|. D'où w(C')*n+|C'| < w(C)*n+|C|, soit score(Q,w) < score(P,w).

  Si C' est contenue dans C, alors w(C') <= w(C) et |C'| < |C| car
  Q<>P. Il suit que w(C')*n + |C'| < w(C)*n + |C|, soit score(Q,w) <
  score(P,w).

  Dans les deux cas nous avons prouvé que score(Q,w) < score(P,w).
  QED

  Soit G un graphe et P un chemin de G tel que G\P est composé d'une
  seule composante connexe (en particulier G est connexe). On définit
  PS1(G,P) le prédicat qui est VRAI ssi pour toute pondération w de G
  telle que P est un plus court chemin il existe dans (G,w) un chemin
  demi-séparateur qui dérive de P.

  Lemme 2. G est dans PS1 ssi PS1(G,{}) = VRAI.

  Preuve. En effet, en réécrivant la définition de PS1(G,{}) on déduit
  que PS1(G,{}) est VRAI ssi pour toute pondération w de G il existe
  dans (G,w) un chemin demi-séparateur qui dérive du chemin vide (tout
  chemin dérive du chemin vide). Notons que c'est nécessairement un
  plus court chemin de (G,w). C'est précisemment la définition de la
  classe PS1. QED

  L'objectif est d'avoir un test noté ps1(G,P) qui implémente
  PS1(G,P), disons qu'il s'en approche. La propriété souhaitée est que
  si ps1(G,P) est VRAI, alors PS1(G,P) aussi. En particulier, si
  ps1(G,{}) est VRAI, alors G est dans PS1.

  Algorithme pour le test ps1(G,P):

  On renvoie VRAI s'il existe une paire x,y de sommets où y n'est pas
  dans P telle que tout chemin Q de x à y:
  1. Q est parallèle à P, et
  2. pour toute composante C de G\(QuP), ps1(CuQ,Q)=VRAI.

  Lemme 3. Si ps1(G,P)=VRAI, alors PS(G,P)=VRAI.

  Preuve [A FINIR]. Par induction sur le nombre de sommets hors de
  P. Soit C est l'unique composante G\P. Si |C|=0, alors P est un
  demi-séparateur de G et donc PS(G,P) est VRAI. Supposons le lemme
  vrai pour tout entier < |C|.

  On suppose donc que ps1(G,P) est VRAI. Soit x,y la paire de sommets
  telle que y n'est pas dans P et où tous les chemins de x à y sont
  parallèles à P. En particulier, pour chaque valuation w, où P est un
  plus court chemin, tout plus court chemin Q selon w entre x et y est
  parallèle à P. Comme Q est différent de P (à cause du choix de y),
  on peut appliquer le lemme 1, et donc soit P est un demi-séparateur,
  soit score(Q,w) < score(P,w). On peut supposer qu'on est pas dans le
  premier cas, c'est-à-dire que P n'est pas un demi-séparateur de G,
  puisque sinon PS(G,P)=VRAI et le lemme est prouvé.

  Si Q est un demi-séparateur pour G, alors PS(G,P) est VRAI puisque Q
  est parallèle à P. Supposons que Q n'est pas un demi-séparateur pour
  G, et soit C' la composante de G\Q telle que w(C')>w(G)/2.

  On peut appliquer l'induction sur (C'uQ,Q) car comme Q<>P,
  |C'|<|C|. Posons G'=C'uQ. D'après le test ps1(G,P), ps1(G',Q) est
  VRAI. Donc par induction PS(G',Q)=VRAI et G' contient un chemin
  demi-séparateur pour la valuation w, disons P', parallèle à
  Q. Montrons d'abord que dans G, P' est parallèle à P.

  ...

  w(C')<w(G)/2 ...

  QED

  Lemme 4. Si ps1(G,P)=VRAI, alors soit il existe une paire x,y de
  sommets de G telle que tout chemin de x à y contient P, ou bien il
  n'existe pas de sommet de P ayant trois voisins inclus dans un cycle
  de G\P. [PAS SÛR, A VÉRIFIER]

  Preuve. [A FINIR] Soit x,y une paire de sommets de G avec y pas dans
  P telle que tous les chemins de x à y soient parallèles à P. Soient
  Q un tel chemin. Supposons que Q ne contient pas P. ...  QED

  Dit autrement, si on n'a pas cette propriété, alors ps1(G,P)=FAUX et
  il est inutile de tester tous les chemins Q possibles.  Remarquons
  que si G est 2-connexe, alors il ne peut avoir de paire x,y de
  sommets (avec y pas dans P) où tout chemin de x et y contient P. A
  montrer: ps1(G,P)=VRAI ssi tout les composantes 2-connexes G' de G
  on a ps1(G',P inter G)=VRAI ...

  [A VOIR]

  - u := ps1(CuQ,Q)
  - ajoute (C,u) à la liste L des composantes maximales (voir ps1_addmax)
  - si u = FAUX, ajouter un nouveau noeud au graphe des conflits (GC)
  - recommencer avec le chemin Q suivant de mêmes extrémités xy (s'il existe)
  - s'il n'y a plus de tel chemin:

    On essaye d'appliquer les trois règles (max, influence, dernier):
    - la règle du max en tenant compte de la liste L.
      si |L|<>1, alors on ne peut pas appliquer la règle du max
      sinon, L={(C,u)}
        si u = VRAI, on peut supprimer la paire xy de GC
        sinon, alors on peut appliquer la règle du max
    - la règle d'influence des voisins
    - la règle du dernier -1

    Si l'application d'une de ces règles (avec ps1_push) produit une
    contradiction dans GC, alors on a trouvé une bonne paire xy, et on
    renvoie VRAI.  Sinon, s'il existe une autre paire xy on recommence
    avec cette pnouvelle paire.

  A la fin du traitement de toutes les paires:

  - soit il reste des indéterminées (-1), et il faut les éliminer. On
    les force d'abord à 0. S'il y a une contradiction (en faisant des
    ps1_push), on en déduit que la valeur doit être 1 (et on fait
    ps1_push). Sinon, on force la valeur à 1. Si y a une
    contradiction, on déduit que la valeur doit être 0 (et on fait
    ps1_push). Sinon, s'il y a une valeur initialement indéterminée
    qui passe à la même valeur pour le forcage à 0 et à 1, on élimine
    cette indéterminée (et on fait ps1_push). On fait ainsi pour
    chaque indéterminée. Chaque fois qu'une indéterminée est éliminée,
    on recommence la recherche sur l'ensemble des indéterminées (et
    pas seulement sur celles qui restent). Si on arrive ainsi à
    éliminer toutes les indéterminées on peut passer à la
    suite. Sinon, on ne peut pas faire grand chose à part essayer tous
    les système possibles ... Voir MINISAT+ (système pseudo booléens)
    et Sugar qui transforme du système linéaire ou CSP en SAT.

  - soit il n'y a plus d'interminées. Dans ce cas on peut déduire un
    système d'équations linéaire indépendantes où les inconnues sont
    les poids des sommets. Si le système n'a pas de solution, alors
    renvoyer VRAI. Sinon, la solution trouvée peut renseigner sur une
    valuation possible pour prouver que G est éventuellement pas dans
    PS1(G,P).

 */

{ 
  G->int1=1; /* compte les tests */

  DEBUG(
	N++;
	int u;int v;
	if(G==GF) N=POS=0;
	printf("\n");
	PRINTS;printf("version=%i G=%p n=%i ",version,G,G->n);PRINTT(P->P,P->n);
	PRINTS;printf("G=");
	for(u=0;u<G->n;u++){
	  printf("%i:[",u);
	  for(v=0;v<G->d[u];v++){
	    printf("%i",G->L[u][v]);
	    if(v<(G->d[u])-1) printf(" ");
	  } printf("] ");
	} printf("\n");
	);
  
  int n=G->n;

  /* Ici on élimine un certain nombre de cas faciles à tester.
     Attention ! vérifier que G est dans PS1 dans ces cas là ne suffit
     pas (ce qui revient à vérifier que PS1(G,{})=VRAI). Il faut être
     certain qu'on a en fait PS1(G,P)=VRAI. */

  if(n-(P->n)<3){
    DEBUG(PRINTS;printf("-> n-|P|<3, moins de 3 sommets hors P\n"););
    DEBUG(N--;);
    return 1;
  }
  /* Par hypothèse P est léger, donc la composante hors de P est
     lourde. Il suffit de prendre un chemin entre ces au plus deux
     sommets hors de P. */

  /* lors d'un appel récursif, le nombre d'arêtes G->m est déjà mis à
     jour car G provient alors du résultat de ExtractSubgraph(). Donc
     NbEdges(G) est calculé en temps constant dès le 2e appel. */

  if((!P->n)&&(n<6)&&(NbEdges(G)<10)){
    DEBUG(PRINTS;printf("-> |P|=0 et n<6 et pas K_5\n"););
    DEBUG(N--;);
    return 1;
  }

  /* ici on a au moins 3 sommets de G qui sont hors de P, et P est non
     vide. Alors, comme P contient au moins deux somemts, cela fait
     que G possède au moins 5 sommets. */

  if(NbEdges(G)==((n*(n-1))>>1)){
    DEBUG(PRINTS;printf("-> clique avec > 2 sommets hors P\n"););
    DEBUG(N--;);
    return 0;
  }
  /* Ici G est clique avec n>4 et n-|P|>2. Dans le cas où tous les
     poids sont à 1 (sommets et arêtes) alors, il n'existe aucun
     chemin Q parallèle P permettant de progresser. Notons qu'une
     clique G avec n sommets donne PS1(G,P)=VRAI s'il y a 0, 1 ou 2
     sommets dans G\P. */

  if(NbEdges(G)<=n){
    DEBUG(PRINTS;printf("-> m<=n\n"););
    DEBUG(N--;);
    return 1;
  }
  /* Dans ce cas, il s'agit d'un arbre avec un cycle. Soit C la
     composante lourde de G\P, et x une des extrémités de P. La
     composante C est connectée à P par au plus deux sommets, disons u
     et v (u=v possible) et u le plus près de x. Soit y le voisin de v
     dans C. On définit alors P' comme le plus court chemin de G
     allant de x à y. Il est facile de voir que P' longe P depuis x
     puis entre dans C soit par u soit par v. Dans les deux cas les
     sommets de C ne peuvent être connectés à aucun sommet de P\P',
     les deux seules arêtes connectant C à P étant détruite par P'. */

  /* ici G possède au moins 5 sommets et 6 arêtes. */

  int x,y,i,u,v,w,d;
  path *Q=new_path(G,NULL,n);
  path *R=new_path(G,NULL,n);
  param_dfs *p=new_param_dfs(n); /* p->C n'est pas alloué */
  graph *C;

  ALLOC(p->C,n);   /* pour le DFS avec sommets supprimés */
  NALLOC(int,T,n); /* pour la composante C de G */
  NALLOC(int,M,n); /* pour la compatiblité des chemins P et Q */

  conflit c; /* ensembles de variables pour gérér le graphe des conflits */
  int npaire,npath,goodxy;
  /* npaire = 1er noeud dans le graphe des conflits de la paire courante */
  /* npath = 1er noeud dans GC du chemin courant pour la paire courante */
  /* goodxy = 1 ssi la paire xy est bonne, ie. tous les chemins et comp. sont ok */

  if(version>0){
    c.n=n; /* nombre de sommets du graphe G */
    c.G=new_graph(CONFMAX); /* graphe des conflits, alloue G->d et G->L */
    c.G->n=c.nbi=c.nbzi=c.ncmax=c.outmem=0; /* c.G->n=nb de noeuds déjà crées */
  }

  /* pour toutes les paires de sommets x,y de G avec y pas dans P */

  for(x=0;x<n;x++)
    for(y=0;y<n;y++){

      if(P->V[y]>=0) continue; /* y ne doit pas être dans P */
      if((P->V[x]<0)&&(y<=x)) continue; /* si x et y pas dans P, ne tester que le cas x<y */

      /* ici on a soit:
	 1) x dans P et y pas dans P, ou bien
	 2) x et y pas dans P et x<y
      */

      (G->int1)++; /* +1 pour chaque paire testée */
      goodxy=1; /* par défaut on suppose la paire xy comme bonne */

      /* calcule le 1er chemin entre x et y */
      Q->P[0]=x; Q->P[1]=y; /* initialise les extrémités du chemin Q */
      if(!NextPath(G,Q,-1)) goto fin_ps1; /* fin si pas de chemin x->y, impossible si G connexe */
      if((version>0)&&(!c.outmem)){
	npaire=c.G->n; /* initialise le 1er noeud de la paire courante */
	c.nbc[npaire]=0; /* nombre de valeurs < 1 pour la paire courante */
	/* on efface la liste des composantes maximales, s'il y en avait */
	for(u=0;u<c.ncmax;u++) free(c.cmax[u]);
	c.ncmax=0;
      }

      DEBUG(PRINTS;printf("Paire (%i,%i) G=%p n=%i\n",x,y,G,n););

      do{ /* pour tous les chemins Q entre x et y */
	(G->int1)++; /* +1 pour chaque sommet testé */

	DEBUG(PRINTS;printf("Essai du chemin ");PRINTT(Q->P,Q->n););

	/* On vérifie que Q est parallèle avec P. Il faut qu'aucun
	   sommet de P\Q n'ait de voisin en dehors de P ou de Q. */

	for(i=0;i<P->n;i++){
	  if(Q->V[u=P->P[i]]>=0) continue; /* on est aussi dans Q */
	  /* ici on est dans P\Q */
	  for(v=0;v<G->d[u];v++){ /* vérifie chaque voisin v de u */
	    if(P->V[G->L[u][v]]>=0) continue; /* si v est dans P */
	    if(Q->V[G->L[u][v]]>=0) continue; /* si v est dans Q */
	    DEBUG(PRINTS;printf("-> chemin Q non parallèle à P\n"););
	    goodxy=0; /* cette paire n'est pas bonne */
	    goto nextxy; /* aller à la prochaine paire */
	  }
	}

	/* ici Q est parallèle à P */

	DEBUG(PRINTS;printf("-> chemin Q parallèle à ");PRINTT(P->P,P->n););

	/* on vérifie que pour chaque composante C de G\Q,
	   PS1(CuQ,Q)=VRAI. Notons que les sommets de P\Q sont
	   léger. Il ne faut donc pas les considérer. Il faut donc
	   prendre les composantes de G\(QuP).*/

	/* on enlève de G les sommets de QuP pour calculer les
	   composantes de G\(QuP) */
	for(u=0;u<n;u++) /* initialise le dfs */
	  if((P->V[u]<0)&&(Q->V[u]<0)) p->C[u]=-1; else p->C[u]=-2;
	dfs(G,0,p); /* calcule les composantes de G\(QuP) */
	DEBUG(PRINTS;printf("#composantes dans G\\(QuP) à tester: %i\n",p->nc););
	/* si aucune composante inutile de tester récursivement ce
	   chemin, G\(QuP) est vide. On peut passer au prochain chemin */
	if(p->nc==0) goto nextQ;
	/* ici, il y a au moins une composante dans G\(QuP) */

	d=R->n=Q->n; /* d=nombre de sommets du chemin Q */
	for(u=0;u<d;u++) T[u]=Q->P[u]; /* T=liste des sommets de Q */
	if((version>0)&&(!c.outmem)) npath=c.G->n; /* initialise le chemin Q au noeud courant */

	DEBUG(PRINTS;printf("Q = ");PRINTT(T,d););

	/* pour chaque composante de G\Q */

	for(i=0;i<p->nc;i++){  /* T[d..v[=i-ème composante de G\(QuP) u Q */
	  for(u=0,v=d;u<n;u++) if(p->C[u]==i) T[v++]=u;
	  /* T[0..v[=sommets de Q u C */
	  /* T[0..d[=sommets de Q, T[d..v[=sommets de la composante */
	  /* NB: les sommets de T[d..v[ sont dans l'ordre croissant */

	  DEBUG(PRINTS;printf("C\\Q=CC(%i) = ",i);PRINTT(T+d,v-d););

	  C=ExtractSubgraph(G,T,v,1); /* crée C=G[T] avec v=|T| sommets */
	  /* Attention! Q n'est pas un chemin de C, mais de G. On crée
	     donc un nouveau chemin R dans C qui est équivalent à Q */
	  for(u=0;u<v;u++) R->V[u]=-1; /* Attention! boucle avec v=C->n sommets */
	  for(u=0;u<d;u++) R->P[u]=C->pint1[Q->P[u]]-1;
	  for(u=0;u<d;u++) R->V[R->P[u]]=u;

	  DEBUG(PRINTS;printf("Q = ");PRINTT(R->P,R->n););

	  /* appel récursif avec une nouvelle version 0 ou 1:
	     - si version=0: sans le graphe des conflits -> 0
	     - si version=1: avec le graphe des conflits -> 1
	     - si version=2: comme version=1 mais sans le graphe
	       des conflits lors de la récursivité -> 0
	     - si version=3: comme version=1 mais avec l'écriture
	       de valeurs dans le graphe des conflits -> 1
	  */
	  u=PS1(C,R,version%2);

	  DEBUG(PRINTS;printf("PS1(CuQ,Q)=%i\n",u););
	  DEBUG(if(c.outmem){PRINTS;printf("PROBLEME MEMOIRE !\n");});

	  G->int1 += C->int1; /* met à jour le nombre de tests (ajoute au moins 1) */
	  free_graph(C); /* libère C qui ne sert plus à rien */

	  /* à faire que si graphe des conflits et pas eut de problème mémoire */
	  if((version>0)&&(!c.outmem)){
	    
	    /* on vérifie si la même composante n'existe pas déjà dans
	       la paire de c.G->n. Si c'est le cas, on passe à la
	       prochaine composante. NB: la composante de c.G->n est
	       dans T[d...v[ et sa taille vaut v-d. */

	    for(w=npaire;w<c.G->n;w++)
	      if(SetCmp(c.comp[w],T+d,c.tcomp[w],v-d)&2)
		goto nextC; /* si même composante, alors prochaine composante */
	    /* ici w=c.G->n */
	    /* si u=0, w sera le nouveau noeud du graphe des conflits */

	    /* ajoute T[d..v[ à la liste des composantes maximales (c.cmax) */
	    if(ps1_addmax(T+d,v-d,&c)){ c.nodemax=w; c.valmax=u; }
	    /* On stocke dans c->nodemax le noeud de la dernière
            composante ajoutée à la liste des composantes maximales.
            Si à la fin du traitement de la paire xy la liste a
            exactement une seule composante, alors c->nodemax sera
            forcément le noeud de la dernière composante ajoutée et
            c->valmax sa valeur (VRAIE/FAUSSE). */
	    if(c.ncmax>=NCMAX){ /* dépassement du nb de composantes maximales */
	      c.outmem=1; /* Ouf of Memory */
	      /* la paire xy n'est pas complète. On l'enlève, sinon
	         l'analyse des conflits pourrait ne va pas être correcte */
	      ps1_delxy(&c,npaire); /* supprime la dernière paire créee */
	      goodxy=0; /* cette paire n'est pas bonne */
	      goto nextxy; /* change de paire */
	    }
	  }
	  /* ici outmem=0 */

	if(u) goto nextC; /* si u=VRAI, on a finit le traitement de la composante */
	goodxy=0; /* une composante n'est pas bonne pour la paire xy */
	if(version==0) goto nextxy; /* si pas graphe de conflits, alors changer de paire */

	/* ici u=FAUX et version>0 */
	/* on va donc essayer d'ajouter le noeud w au graphe des conflits */
	/* ici pas de problème de mémoire, donc on ajoute le noeud w */

	(c.G->n)++; /* un noeud de plus */
	if(c.G->n==CONFMAX){ /* Out of Memory */
	  c.outmem=1;
	  ps1_delxy(&c,npaire); /* supprime la dernière paire créee, qui n'est pas bonne */
	  goto nextxy; /* ici goodxy=0 */
	}
	/* ici v=|T|, d=|Q|, w=nouveau noeud à créer */
	
	/* Attention! ne pas utiliser ps1_push(w,...) juste après la
	   création de w, car on est pas certain que la paire de w
	   sera complète ... On ne peut faire ps1_push(w,...)
	   seulement après le while(NextPath...). Idem pour
	   l'application de la règle du max ! */

	/* initialise le nouveau noeud w */

	DEBUG(PRINTS;printf("Ajoût du noeud %i dans le graphe des conflits\n",w););

	c.x[w]=x; c.y[w]=y; /* mémorise la paire x,y (pour l'affichage) */
	c.path[w]=npath; c.paire[w]=npaire; /* w rattaché à npath et à npaire */
	c.code[w]=-1; c.nbi++; /* un noeud de plus à valeur = -1 */
	(c.nbc[npaire])++; /* un noeud de plus à valeur < 1 pour cette paire */
	c.G->d[w]=0; /* initialise le degré de w */
	ALLOC(c.G->L[w],CONFMAX); /* crée la liste des voisins de w */
	ALLOCZ(c.comp[w],v-d,T[d+(_i)]); /* crée & copie la composante de w */
	c.tcomp[w]=v-d; /* taille de la composante de w */
	
	/* y'a-t'il des arcs vers les noeuds précédant w ? */
	/* calcule les arêtes sortantes de w et de leurs types */
	
	for(u=0;u<w;u++){ /* pour chaque noeud < w, v=type de l'arc */
	  if(u>=npath) v=0; /* test rapide: si u>=napth, alors clique */
	  else{ /* sinon calcule l'intersection des composantes */
	    v=SetCmp(c.comp[u],c.comp[w],c.tcomp[u],c.tcomp[w]);
	    /* ici les valeurs possibles pour SetCmp: 0,1,3,5,9 car T1 et T2 != {} */
	    if(v==1) v=-2; /* si v=1, T1 intersecte T2, et donc pas d'arête u-w */
	    /* avant: v: 0=disjoint, 3=égalité, 5=T1 sub T2, 9=T2 sub T1 */
	    v >>= 1; v -= (v==4); /* rem: si avant v=-2, alors après v=-1 */
	    /* après: v: 0=disjoint, 1=égalité, 2=T1 sub T2, 3=T2 sub T1 */
	  } /* si v<0, alors pas d'arête */
	  if(v>=0){ /* ajoute les arcs u->w et w->u */
	    c.G->L[u][c.G->d[u]++]=(v<<LCONF)|w; /* u->w */
	    if(v>1) v ^= 1; /* si v=2, alors v->3, et si v=3, alors v->2 */
	    c.G->L[w][c.G->d[w]++]=(v<<LCONF)|u; /* w->u (asymétrie pour v=2,3) */
	  }
	} /* fin du for(u=...) */
	
	nextC:; /* prochaine composante (next i) */
	} /* fin du for(i=0...) */

      nextQ:; /* prochain chemin Q */
      }while(NextPath(G,Q,0));
      
      /* si ici goodxy=1, c'est que pour la paire xy, tous les chemins
	 entre xy sont parallèles à P et qu'on a jamais vue de
	 composante à FAUX. Dans ce cas on a trouvé une bonne paire et
	 on a fini. */

      if(goodxy) goto fin_ps1; /* termine l'algo. */
      if(version==0) goto nextxy; /* si pas graphe des conflits */
      /* ici version>0 et goodxy=0 */
      
      /* Ici on a examiné tous les chemins Q pour la paire xy et on
      crée au moins un noeud dans le graphe des conflits. Ici tous les
      noeuds crées ont comme valeur -1. Il reste à essayer d'appliquer
      trois règles:

	 1) la règle du max: on écrit 0 s'il existe une composante de
	 xy contenant toutes les autres (y compris celle évaluées
	 récursivement à VRAI). Pour cela il faut avoir exactement une
	 seule composante maximale dans c.cmax qui doit être à
	 FAUX. Si elle est à VRAI il faut supprimer la paire xy (car
	 la plus lourde est VRAI, donc toutes les autres sont légères
	 !) et passer à la suivante.

	 2) la règle d'influence des voisins: vérifier si des voisins
	 v de chaque noeud u (de la paire xy) ne peuvent pas
	 influencer u. Par exemple, si la composante de v est la même
	 que celle de u il faut alors écrire cette valeur dans le code
	 de u.

	 3) la règle du dernier -1: si le nb de valeurs < 1 est
	 exactement 1, alors il faut écrire 0 dans ce noeud là. Il
	 faut, comme pour la règle d'influence des voisins, balayer
	 tous les sommets u de la paire xy.

	 NB: on peut traiter dans un même parcours les deux dernières
	 règles.
      */

      v=c.paire[c.G->n-1]; /* v=1er noeud de la dernière paire */
      /* NB: c.outmem=0 et c.G->n > 0 */

      /* règle du max */
      if(c.ncmax==1){ /* sinon pas de règle du max */
	if(c.valmax){ /* supprime la paire xy */
	  ps1_delxy(&c,v);
	  goto nextxy;
	}
	/* applique la règle du max: on écrit 0 dans la composante maximale */
	if(ps1_push(c.nodemax,0,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
      }

      /* règle d'influence des voisins + règle du dernier -1 */
      for(u=v;u<c.G->n;u++){ /* pour tout noeud u de la dernière paire */
	if(c.code[u]>=0) continue; /* plus rien à faire */
	/* NB: si code[u]=0 ou 1, c'est qu'on a forcément déjà fait
	   un ps1_push() sur u. Et donc la règle d'influence des
	   voisins n'a pas a être testée sur u */
	else /* code=-1 */
	  if(c.nbc[c.paire[v]]==1){ /* règle du dernier -1 */
	    /* tous les sommets sont à 1, sauf u qui est ici a -1 */
	    if(ps1_push(u,0,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
	    else break; /* on peut sortir du for(u...), on a tout testé */
	  }
	
	for(i=0,d=c.G->d[u];i<d;i++){ /* scan tous les voisins v de u */
	  v=c.G->L[u][i]&CONFMASK; /* v=i-ème voisin de u */
	  w=c.code[v]; if(w<0) continue; /* on ne déduit rien */
	  /* efface le code de v pour pouvoir le réécrire avec
	     ps1_push().  Il faut prendre des précaution pour que c
	     soit cohérent (il faut mettre à jour .nbi et .nbc). NB:
	     l'ancien code de v est w=0 ou 1, celui de u est aussi 0 ou 1 */
	  c.code[v]=-1;
	  c.nbi++; /* met à jour le nombre d'indéterminées */
	  c.nbc[c.paire[v]] += w; /* augmente s'il y avait 1 */
	  /* réécrit le code pour v */
	  if(ps1_push(v,w,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
	  /* NB: la propagation de cette écriture va tester l'arc
	     v->u et donc modifier éventuellement le code de u */
	}
      }

    nextxy:;
    } /* fin du for(x,y=...) */
  
  /* Ici on a testé toutes paires xy et aucune n'est bonne */
  /* NB: goodxy=0 */

  if((version>0)&&(!c.outmem)){

    /* Traitement final du graphe des conflits. On pourrait ici
       essayer d'éliminer les derniers -1. Seulement, il semble
       qu'aucune des règles ne peut plus être appliquées. */

  loop_ps1:
    do{
      if(c.nbi==0) break; /* fini: il ne reste plus de -1 à modifier */
      x=c.nbi; /* mémorise le nb total de valeurs à -1 */
      //
      // que faire ? essayer toutes les valeurs possibles pour les -1 ?
      //
      // G->int1++;
    }
    while(c.nbi<x); /* tant qu'on a enlevé au moins un -1, on recommence */
    
    /* ps1x */
    if(version==3){ /* NB: si appel récursif, alors version!=3 */
      i=MEM(CPARAM,0,int);
      if(i>0){ /* lit les valeurs à l'envers */
	u=MEM(CPARAM,(2*i-1)*sizeof(int),int);
	v=MEM(CPARAM,(2*i)*sizeof(int),int);
	MEM(CPARAM,0,int)--;
	if(ps1_push(u,v,&c)){ goodxy=1; goto fin_ps1; } /* fin si contradiction */
      goto loop_ps1;
      }
    }
    
    /*
      on cherche les noeuds "0" indépendants, correspondant à des
      composantes lourdes (donc à des inégalités).
      
      Algorithme: on ne sélectionne que les noeuds de code=0 et qui
      n'ont pas de voisins v de code aussi 0 avec v<u (inclus). Ils
      doivent donc être de code 0 et minimal par rapport aux autres
      voisins de code 0. Si un tel noeud existe, on marque tous ces
      voisins de code 0 qui ne pourront plus être sélectionnés.
    */

    c.nbzi=0; /* =nb de zéros indépendants */
    for(u=0;u<c.G->n;u++){ /* pour chaque noeud du graphe des conflits */
      if(c.code[u]) continue; /* saute le noeud si pas 0 */
      for(i=0,d=c.G->d[u];i<d;i++){ /* scan tous les voisins u (qui est forcément à 0) */
	v=c.G->L[u][i]; w=c.code[v&CONFMASK];
	if(!((w==0)||(w==CONFC))) continue; /* on ne regarde que les voisins de code=0 */
	if((v>>LCONF)==3) break; /* ici i<d */
      }
      if(i==d){ /* on a trouvé un noeud u de code 0 et sans aucun arc u->v avec v<u */
	c.nbzi++; /* un de plus */
	for(i=0,d=c.G->d[u];i<d;i++){ /* modifie le code de tous les voisins à 0 de u */
	  v=c.G->L[u][i]&CONFMASK; /* v=voisin de u */
	  if(c.code[v]==0) c.code[v]=CONFC; /* on ne modifie pas deux fois un voisin */
	}
      }else c.code[u]=CONFC; /* u n'est pas bon, on ne laisse pas la valeur 0 */
    }

  }
  /* ici on a fini le traitement du graphe des conflits */

 fin_ps1: /* termine l'algo avec goodxy=0 ou 1 */

  if(version>0){ /* affichage et libération du graphe des conflits */
    if((!goodxy)&&(G==GF)) PrintConflit(&c); /* G=GF => affiche que le 1er niveau de récurrence */
    /* efface la liste des composantes maximales */
    for(u=0;u<c.ncmax;u++) free(c.cmax[u]);
    /* efface le graphe des conflits */
    for(u=0;u<c.G->n;u++) free(c.comp[u]);
    free_graph(c.G); /* c.G->L après c.G->n n'ont pas été alloués ou sont déjà libérés */
  }

  /* efface les tableaux alloués */
  free_path(Q);
  free_path(R);
  free_param_dfs(p);
  free(T);
  free(M);

  DEBUG(N--;);
  return goodxy;
}


/***********************************

       ROUTINES POUR LES
        ROUTING SCHEMES

***********************************/


char *TopChrono(int i)
/*
  Met à jour le chronomètre interne numéro i (i=0..CHRNONMAX-1) et
  renvoie sous forme de char* le temps écoulé depuis le dernier appel
  à la fonction pour le même chronomètre. La précision dépend du temps
  mesuré. Il varie entre la seconde et le 1/1000 de seconde. Plus
  précisément le format est le suivant:

  1h00'00"  si le temps est > 60' (seconde)
  1'00"0    si le temps est > 1'  (1/10e de seconde)
  1"00      si le temps est > 1"  (1/100e de seconde)
  0"000     si le temps est < 1"  (1/1000e de seconde)

  Pour initialiser et mettre à jour tous les chronomètres (dont le
  nombre vaut CHRONOMAX), il suffit d'appeler une fois la fonction,
  par exemple avec TopChrono(0). Si i<0, alors les pointeurs alloués
  par l'initialisation sont désalloués.
*/
{
  if(i>=CHRONOMAX) Erreur(26);

  /* variables globales, locale à la fonction */
  static int first=1; /* =1 ssi c'est la 1ère fois qu'on exécute la fonction */
  static char *str[CHRONOMAX];
  static clock_t last[CHRONOMAX];
  int j;

  if(i<0){ /* libère les pointeurs */
    if(!first) /* on a déjà alloué */
      for(j=0;j<CHRONOMAX;j++)
	free(str[j]);
    first=1;
    return NULL;
  }

  if(first){ /* première fois, on alloue et on renvoie TopChrono(i) */
    first=0;
    for(j=0;j<CHRONOMAX;j++){
      str[j]=malloc(sizeof("00h00'00\""));
      last[j]=clock();
    }
  }

  /* t=temps en 1/1000" écoulé depuis le dernier appel à TopChrono(i) */
  long t=(long)(1000*(clock()-last[i])/CLOCKS_PER_SEC);
  last[i]=clock(); /* met à jour le chrono interne */
  if(t<0L) t=0L;
  
  for(;;){ /* pour faire un break */
    if(t<1000L){ /* t<1" */
      sprintf(str[i],"0\"%03li",t);
      break;
    }
    t /= 10L; /* t en centième de seconde */
    if(t<6000L){ /* t<60" */
      sprintf(str[i],"%li\"%02li",t/100L,t%100L);
      break;
    }
    t /= 10L; /* t en dixième de seconde */
    if(t<36000L){ /* t<1h */
      sprintf(str[i],"%li'%02li\"%li",t/360L,(t/10L)%60L,t%10L);
      break;
    }
    t /= 10L; /* t en seconde */
    if(t>=360000L){ /* t>99h59'59" */
      sprintf(str[i],"99h59'59\"");
      break;
    }
    sprintf(str[i],"%lih%02li'%02li\"",t/3600L,(t/60L)%60L,t%60L);
  }
  
  return str[i];
}


/* type pour une fonction f(u,v,T) renvoyant la longueur d'une route
   de u vers v étant données la table de routage globale T */
typedef int(*rt_length)(int,int,void*);


/* structure de données pour la table de routage d'un sommet */
typedef struct{
  int n;      /* taille des listes */
  int *node;  /* liste de sommets */
  int *dist;  /* distances */
  int radius; /* rayon d'une boule */
  int vpd;    /* voisin par défaut */
} table;


/*

pour avoir plusieurs noms différents pour une même structre

typedef union{
  struct{
    int n;
    int *dist;
  };
  struct{
    int n;
    int *color;
  };
} tableau;

  tableau *T=malloc(sizeof(*T));
  T->n=0;
  T->dist=NULL;
  printf("->y=%p\n",T->color);

*/

table *new_table(int n){
/*
  Crée un objet de type table où les champs sont initialisés à leur
  valeurs par défaut, les pointeurs étant initialisés à NULL. Si n>0,
  alors les pointeurs (->node et ->dist) de taille n sont alloués, et
  le champs ->n est initialisé à n.
*/
  NALLOC(table,X,1);
  X->n=0;
  X->node=NULL;
  X->dist=NULL;
  X->radius=-1;
  X->vpd=-1;

  if(n>0){
    X->n=n;
    ALLOC(X->node,n);
    ALLOC(X->dist,n);
  }
  
  return X;
}


void free_table(table *X){
  if(X==NULL) return;
  free(X->node);
  free(X->dist);
  free(X);
  return;
}


/* tables de routage globale pour rs_cluster() */
typedef struct{
  int n; // pour libérer les n tables
  table **B; // les boules
  table **S; // les spanning trees depuis les landmarks
  table **R; // les sommets en charge des landmarks
  table **W; // tables des voisins dans le cluster
  int *H; // hahs des sommets
  int *C; // pour routage vers voisins de couleur i
  int center; // centre du cluster (pour VARIANT=1)
} rs_cluster_tables;


void free_rs_cluster_tables(rs_cluster_tables *X){
/*
  Libère les tables créées par rs_cluster().
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->S[u]);
    free_table(X->R[u]);
    free_table(X->W[u]);
  }
  free(X->B);
  free(X->S);
  free(X->R);
  free(X->W);
  free(X->H);
  free(X->C);
  free(X);
  return;
}


/* tables de routage globale pour rs_dcr() */
typedef struct{
  table **B; // B[u]=tables de voisinage indexées par les sommets
  table **W; // W[u]=tables des représentants indexées par les couleurs 
  param_bfs **S; // S[u]=bfs(u) pour u landmark (sinon =NULL)
  int *H; // H[u]=hash du sommet u
  int *C; // C[u]=couleur du sommet u
  int **dist; // distance partielle (issues des landmarks)
  int n; // pour libérer les n tables
} rs_dcr_tables;


void free_rs_dcr_tables(rs_dcr_tables *X){
/*
  Libère les tables créées par rs_dcr().  On a pas besoin de libérer
  les pointeurs de X->dist, car ce sont ceux de X->S->D.
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->W[u]);
    free_param_bfs(X->S[u]); // libère aussi X->dist
  }
  free(X->B);
  free(X->W);
  free(X->S);
  free(X->H);
  free(X->C);
  free(X->dist);
  free(X);
  return;
}


/* tables de routage pour rs_tzrplg() */
typedef struct{
  table **B; // boules
  table **L; // landmarks
  int* label; // les étiquettes des sommets
  int n; // pour libérer les n tables
} rs_tzrplg_tables;


void free_rs_tzrplg_tables(rs_tzrplg_tables *X){
/*
  Libère les tables créees par rs_cluster().
*/
  if(X==NULL) return;
  int u;
  for(u=0;u<X->n;u++){
    free_table(X->B[u]);
    free_table(X->L[u]);
  }
  free(X->B);
  free(X->L);
  free(X);
  return;
}


/* tables de routage pour rs_bc() */
typedef struct{
  int n; // pour libérer les n tables
} rs_bc_tables;


void free_rs_bc_tables(rs_bc_tables *X){
/*
  Libère les tables créees par rs_bc().
*/
  if(X==NULL) return;
  free(X);
  return;
}


int fcmp_couple_x(const void *P,const void *Q)
/* Compare le champs x de deux couples d'entiers, pour qsort(). */
{
  return ((couple*)P)->x - ((couple*)Q)->x;
}


int fcmp_stretch(const void *P,const void *Q)
/*
  Compare les ratios P.x/P.y et Q.x/Q.y, dans le cas de ratios
  irréductibles. Sert pour routing_test().
*/
{
  return ((triplet*)P)->x*((triplet*)Q)->y - ((triplet*)P)->y*((triplet*)Q)->x;
}


void SortTable(table *X)
/*
  Trie la table X par ordre croissant du champs X->node[], les autres
  champs (->dist[]) étant rangés selon ce nouvel ordre s'ils ne sont
  pas NULL.
*/
{
  if((X==NULL)||(X->n<2)) return; /* pas de tri à faire */
  if(X->dist==NULL){
    qsort(X->node,X->n,sizeof(int),fcmp_int);
    return;
  }
  
  NALLOC(couple,T,X->n);

  int i;
  for(i=0;i<X->n;i++){
    T[i].x=X->node[i];
    T[i].y=X->dist[i];
  }

  qsort(T,X->n,sizeof(couple),fcmp_couple_x);
  
  for(i=0;i<X->n;i++){
    X->node[i]=T[i].x;
    X->dist[i]=T[i].y;
  }

  free(T);
  return;
}

/* affiche le min/max et la moyenne d'un tableau d'entiers */
// Ex:
//     MINMAXMOY(T[_i],n,T[_i]>0,"size of T");
//     Calcule la moyenne des termes T[_i] tq T[_i]>0
//     pour i=0..n-1. Le texte sert pour l'affiche.
//
#define MINMAXMOY(term,n,condition,string)				\
  do{									\
    int _i,m0,m1,cpt;long sum=0L;					\
    for(_i=cpt=0;_i<(n);_i++)						\
      if(condition){							\
        sum+=(long)term;						\
	if(cpt){ m0=min(m0,term); m1=max(m1,term); }			\
	else m0=m1=term;						\
	cpt++;								\
      }									\
    printf("- minimum %s: %i\n",string,m0);				\
    printf("- maximum %s: %i\n",string,m1);				\
    printf("- average %s: %.2lf (%li/%i)\n",string,			\
	   cpt?(double)sum/(double)cpt:0,sum,cpt);			\
  }while(0)


/* affiche la fréquence, le min et le max d'un tableau d'entiers */
// Ex:
//     FREQMINMAX(F,k,T,n,"size of T");
//       Calcule dans F la fréquence des éléments d'un tableau T
//       à n éléments, les valeurs de T étant des entiers de [0,k[.
//       Attention ! le tableau F est alloué par et initialisé par
//       la macro. Le texte sert pour l'affiche.
//
#define FREQMINMAX(F,k,T,n,string)					\
  do{									\
    int u,v,r,i;							\
    F=SortInt(T,NULL,n,0,&k,SORT_FREQv);				\
    for(u=v=F[i=0];i<k;i++) u=min(u,F[i]),v=max(v,F[i]);		\
    r=n/k+((n%k)>0); /* ceil(n/k) */					\
    printf("- balance ratio of %s: %.0lf%% (ceil{n/k}=%i max=%i min=%i)\n", \
	   string,100.0*u/(double)v,r,v,u);				\
  }while(0)								\


void ruling(char *s,int n){
/*
  Affiche n fois la chaîne de caractère s
*/
  int i;
  for(i=0;i<n;i++) printf("%s",s);
}
#define BARRE do{ruling("―",LEN_BARRE);printf("\n");}while(0)
#define RULING(x) ruling(STR_DISTRIB,min(LEN_DISTRIB*x,LEN_DISTRIB-5))


void PrintDistrib(int *Z,int n,int r,char *message){
/*
  Affiche la distribution des valeurs entières contenues dans le
  tableau Z (de taille n>0) selon r>0 intervalles et avec le message
  m. Les intervalles consécutifs nuls sont regroupés. Donc le nombre
  d'intervalles affichés peut être < r.

  Exemple avec m="routing table size":

  - routing table size distribution: 10 ranges
    [339,381[ 	01% ▪ [168] 
    [381,423[ 	08% ▪▪▪▪ [698] 
    [423,465[ 	15% ▪▪▪▪▪▪▪▪▪ [1282] 
    [465,507[ 	57% ▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪▪ [4848] 
    [507,549[ 	11% ▪▪▪▪▪▪ [975] 
    [549,591[ 	02% ▪ [230] 
    [591,633[ 	02% ▪ [182] 
    [633,675[ 	00%  [18] 
    [675,717[ 	00%  [4] 
    [717,752[ 	00%  [1] 
  - minimum routing table size: 339
  - maximum routing table size: 751
  - average routing table size: 481.48 (4047360/8406)
*/

  int i,j;
  int m0=Z[0]; /* calcule m0=min_i{Z[i]} */
  int m1=Z[0]; /* calcule m1=max_i{Z[i]} */
  for(i=0;i<n;i++){ m0=min(m0,Z[i]); m1=max(m1,Z[i]); }
  m1++; /* m1=max{Z}+1 */
  int d=m1-m0; /* d=plage des valeurs, d>=1 */
  r=min(r,d); /* r=le nombre d'intervalles maximum ne peut dépasser d */
  d=((d-1)/r)+1; /* d=ceil(d/r)=taille d'un intervalle, d>=1 car d>=r>=1 */  
  NALLOCZ(int,F,r,0); /* initialise F */
  for(i=0;i<n;i++) F[(Z[i]-m0)/d]++; /* compte les valeurs de chaque intervalle */
  int k=r; /* k=nombre d'intervalles qui vont être affichés */
  for(i=1;i<r;i++) if((F[i]==0)&&(F[i-1]==0)) k--;
  printf("- %s distribution: %i range%s\n",message,k,PLURIEL(k));

  /* affichage des itnervalles, éventuellement en les fusionnant */
  double x;
  for(i=0;i<r;){
    j=i+1; while((j<r)&&(F[j]==0)) j++; /* j=prochain intervalle */
    x=(double)F[i]/(double)n;
    printf("  [%i,%i[ \t%02i%% ",m0+i*d,min(m0+j*d,m1),(int)(100.0*x));
    RULING(x);
    printf(" [%i] \n",F[i]);
    i=j; /* j=i+1 ou plus */
  }
  free(F); /* ne sert plus à rien */
  MINMAXMOY(Z[_i],n,1,message);
  
  return;
}


char *millier(long i){
/*
  Renvoie une chaîne de caractères correspondant à l'entier long i
  écrit avec le séparateur des milliers ','. Ne pas faire de free()
  sur le pointeur renvoyé.

  Ex: > printf("n=%s\n",millier(-123456789));
      > n=-123,456,789
*/  
  static char r[64];
  char *p,*s=r; /* r=résultat final, s=r ou r+1 si i<0 */
  int n;
  
  snprintf(r,sizeof(r),"%li",i); /* écrit i dans r */
  if(i<0) s++; /* si i<0, r="-123..." */
  n=strlen(s); /* n=longueur(r) sans le signe */
  p=s+n-3;
  
  while(p>s){
    memmove(p+1,p,n);
    *p=',';
    n+=4; /* 4 chiffres avec la virgule */
    p-=3; /* paquet de 3 chiffres */
  }
  
  return r;
}


int lg(int n){
/*
  Renvoie le plus petit entier k tel que n < 2^k. Autrement dit k est
  tel que 2^{k-1} <= n < 2^k. D'où:

            / 1+floor{log_2(n)} si n>0     n  | 0 1 2 3 4 5 6 7 8 9 ...
   lg(n) := |                            -----o------------------------
            \ 0 sinon                    lg(n)| 0 1 2 2 3 3 3 3 4 4 ...

  Lorsque n>0, c'est aussi le nombre de bits dans l'écriture binaire
  de n, et c'est également un de plus que la position du bit de poids
  fort.
*/
  int k=0;
  while(n>0) n>>=1,k++;
  return k;
}


int hash_prime(int u){
/*
  Renvoie un hash pour u, un entier de [0,p[ où p=2^31-1=0x7FFFFFFF
  est un nombre premier. Si u<0, alors les constantes A,B (qui
  accélèrent les calculs) sont initialisées et on renvoie -1.  On
  prend la fonction de hashage de Carter & Wegman (1978) avec
  h(x)=((A*x+B) mod p) mod k) où 0<A,B<p et p>=n est premier. La
  probabilité de collision est alors floor(p/k)*(p-1). Il est
  important de faire les calculs en "long unsigned" car la valeur
  A*x+B dépasse 32 bits en général. Sinon, le calcul de (A*x+B)%p est
  erroné (s'il est réalisé en 32 bits).
*/
  static const long unsigned p=0x7FFFFFFF;
  static long unsigned A,B;

  if(u<0){
    A=(1+random())%p;
    B=random()%p;
    return -1;
  }

  return (A*(long unsigned)u+B)%p;
}


int hash_shuffle(int u,int n){
/*
  Calcule un entier p(u) de [0,n[ qui est une permutation
  p:[0,n[->[0,n[ basée sur deux constantes aléatoires R1,R2 de
  [0,n[. Le temps de calcul est d'environ (logn)/2. Si u<0, alors les
  constantes R1,R2,K,N0,N1 (qui accélèrent les calculs) sont
  initialisées et on renvoie -1. Cette permutation est basée sur deux
  "shuffles" qui sont les permutations P0 et P1 suivantes:

    P0(u)=(u>>1)+((u&01)==0)*N0; // avec N0=floor(n/2)
    P1(u)=(u>>1)+((u&01)==1)*N1; // avec N1=ceil(n/2)

  Pour n=10, la permutation P0(u) donne: 0-5-1-6-2-7-3-8-4-9. Donc on
  écrit 0,1,2,... aux positions paires 0,2,4,6,...

  Pour n=10, la permutation P1(u) donne: 5-0-6-1-7-2-8-3-9-4. Donc on
  écrit 0,1,2,... aux positions impaires 1,3,5,6,...

  On applique la permutation P0 ou P1 suivant les bits de poids faible
  de R1, en répétant cette opération K fois, où K=(logn)/2 est environ
  la moitié des bits nécessaires pour écrire n en binaire. En effet,
  on retrouve la valeur de départ après avoir effectué logn fois la
  permutation P0 ou logn fois la permutation P1.

  Enfin, on inverse le résultat (u->n-u) et on le décale de R2
  position modulo n (u->(u+R2)%n). Ainsi la permutation finale de u,
  p(u), peut valoir n'importe quelle valeur entre [0,n[.
*/

  /* variables globales / locales */
  static int R1=0;
  static int R2=0;
  static int N0=0;
  static int N1=0;
  static int K=0;

  if(u<0){
    R1=random()%n; /* R1 contient au moins k bits aléatoires */
    R2=random()%n; /* pour le décalage final */
    N0=N1=(n>>1);  /* N0=floor(n/2) */
    N1 += (n&01);  /* N1=ceil(n/2) */
    K=lg(n)/2;     /* K=la moitié de floor(logn) */
    return -1;
  }
  
  /* calcule u=p(u) */
  int i,b,r=R1;
  for(i=0;i<K;i++){
    b=(r&01); /* lit un bit b de R1 et applique P0 ou P1 */
    u=(u>>1)+((u&01)==b)*(b?N1:N0); /* u/2, u/2+N0 ou u/2+N1 */
    r>>=1; /* le supprime le bit lu */
  }

  /* inverse et décalage aléatoire */
  return (n+R2-u)%n;
}


int *MakeHash(int *H,int n,int k,int M){
/*
  Calcule un hash h:[0,n[->[0,k[ avec k<=n selon la méthode M (voir
  HASH). Remplit et renvoie le tableau H de taille n (ou bien un
  nouveau tableau alloué si H=NULL). Si M n'est pas une valeur
  reconnue, on fait comme si M=H_MOD. On renvoie NULL si n<=0.
*/

  if(n<=0) return NULL;
  if(H==NULL) ALLOC(H,n);
  
  int u;
  
  switch(M){
    
  case H_PRIME:
    hash_prime(-1); /* initialisation */
    for(u=0;u<n;u++) H[u]=hash_prime(u)%k;
    break;
    
  case H_SHUFFLE:
    hash_shuffle(-1,n); /* initialisation */
    for(u=0;u<n;u++) H[u]=hash_shuffle(u,n)%k;
    break;
    
  case H_MOD:
  default:
    for(u=0;u<n;u++) H[u]=u%k;
    break;
  }
  
  return H;
}


rs_cluster_tables *rs_cluster(graph *G,int k)
/*
  Calcule pour le graphe G le routing scheme cluster de paramètre
  k>0. On renvoie les tables calculées. Le stretch est toujours <=
  5. Cependant lorsque k=1 le stretch devient <= 3.

  Lorsque VARIANT=1, le cluster est fixé au voisinage du center. Puis
  on utilise le routage dans le cluster comme dans une étoile, sans
  tables W.

  Lorsque VARIANT=2, les tables B sont vidées, le stretch max n'est
  alors plus forcément borné. Si de plus k=1, le routage s'effectue
  donc selon un arbre de racine le center.
*/
{
  param_bfs *X0,*X;
  int u,v,i,j,r,t,center;
  table **S,**B,**W,**R;
  int *C,*Z,*H,*F;

  printf("\nCLUSTER\n");
  BARRE;
  
  TopChrono(1); /* reset du chrono */
  
  /* trouve un sommet center de degré max, v=deg(center) */
  for(u=v=center=0;u<G->n;u++)
    if(G->d[u]>v) v=G->d[center=u];

  if(VARIANT==1) k=v+1; /* tous les voisins */
  k=min(k,v+1); /* k=taille effective du cluster=#colors */
  printf("- cluster size: %i\n",k);
  printf("- degree of the center: %i (id:%i)\n",v,center);

  /* construit le cluster C de taille k à partir des k-1 sommets de
     plus haut degré voisins du center */

  /* construit une table H des degrés des voisins de center */
  /* on en profite pour calculer le degré min pour optimiser SortInt() */
  u=v; /* u=futur degré min, v=deg(center) */
  ALLOC(H,v); /* H[i]=degré du i-ème voisin du center, taille(H)=deg(center) */
  for(i=0;i<v;i++){ H[i]=G->d[G->L[center][i]]; u=min(u,H[i]); }

  /* trie selon les listes des degrés H: Z[i]=indice dans H du degré de rank i */
  r=v-u+1; /* H[i] valeur dans [u,v]=[u,u+r[ */ 
  Z=SortInt(H,NULL,v,u,&r,SORT_INDEXi); /* H[Z[v-1]]=+grand degré parmi les v voisins du center */
  free(H); /* H ne sert plus à rien, seul Z sert encore */

  /* remplit C selon les degrés décroissant */
  ALLOC(C,k);
  C[0]=center; /* met le center dans C, forcément de plus haut degré */
  for(i=1;i<k;i++) C[i]=G->L[center][Z[--v]]; /* les k-1 voisins de + haut degré */
  free(Z); /* Z ne sert plus à rien */
  
  /* affiche un aperçu des degrés du cluster */
  printf("- cluster degree: ");
  APERCU(G->d[C[_i]],k,10,2);

  /* on trie C pour déterminer la couleur de u dans C rapidement. La
     couleur de u dans C est alors i=SetSort(C,k,1). Attention ! la
     couleur du center n'est pas forcément 0. */
  qsort(C,k,sizeof(int),fcmp_int);
  
  X0=bfs(G,center,NULL); /* calcule un BFS depuis center */
  printf("- eccentricity of the center: %i\n",X0->radius);
  printf("- time to construct the cluster with its BFS: %s\n",TopChrono(1));

  /*
    Calcule, pour chaque sommet u de C, le nombre de sommets qui ne
    sont pas dans C et qui ont u comme plus proche ancêtre. On calcule
    ensuite le maximum de ces nombres.

    On utilise un simple tableau H indiquant si un sommet a déjà été
    visité ou pas. De plus, si c'est un sommet de C, H[u] indique
    combien de sommets lui sont directement descendant. Au départ
    H[u]=1 si u est dans C, et H[u]=-1 sinon. Puis, pour chaque sommet
    u qui n'est pas déjà visité ou dans C, on marque u et on remonte
    dans l'arbre u=parent(u) jusqu'à trouver soit un sommet de l'arbre
    déjà visité ou bien un sommet de C. Pour marqué un sommet u que
    l'on visite, on pose H[u]=0. Si on tombe sur un sommet du cluster,
    on incrémente H[u]. Donc H[u]>0 si u est dans C, H[u]<0 si u n'a
    pas été visité, et H[u]==0 sinon (pas dans C mais déjà
    visité). Cet algo prend un temps O(n) car chaque arête de l'arbre
    n'est visitée qu'une seule fois.
   */
  ALLOCZ(H,G->n,-1);
  for(i=0;i<k;H[C[i++]]=1);
  for(u=0;u<G->n;u++){
    v=u;
    while(H[v]<0){ /* v n'est pas marqué */
      H[v]=0; /* on marque v */
      v=X0->P[v]; /* v=parent(v) */
    }
    if(H[v]>0) H[v]++;
    }
  j=0; /* cherche u de C avec H[u] maximum */
  for(i=1;i<k;i++) if(H[C[i]]>H[C[j]]) j=i;
  printf("- maximum number of cluster direct descendants: %i (id: %i)\n",H[C[j]],C[j]);
  free(H);

  /* efface les pères des sommets de C pour test d'appartenance rapide */
  for(i=0;i<k;i++) X0->P[C[i]]=-2; /* si -2 alors dans C */

  DEBUG(PRINTT(C,k););

  /****************************/
  /* calcule des stats / hash */
  /****************************/

  /* calcule le nombre d'arêtes intra et extra C, de voisins de C */
  TopChrono(1);
  ALLOCZ(H,G->n,0);
  t=r=0;
  for(i=0;i<k;i++){
    u=C[i]; // u=i-ème sommet de C
    for(j=0;j<G->d[u];j++){ /* pour chaque voisin de u */
      v=G->L[u][j]; /* v=j-ème voisin de u */
      if(X0->P[v]==-2) t++; // arête dans C
      else{ r++; H[v]=1; } // arête hors C; marque les sommets
    }
  }
  printf("- #edges inside the cluster: %i\n",t/2);
  printf("- #edges outgoing the cluster: %i\n",r);
  for(u=t=0;u<G->n;u++) t+=H[u];
  printf("- #neighbors of the cluster: %i\n",t);
  printf("- their average #neighbors in the cluster: %.2lf (%i/%i)\n",(double)r/(double)t,r,t);
  
  /* construit le hash H[u]=0..k-1 de chaque sommet u=0..n-1 */
  /* le sommet C[i] est en charge des sommets de hash i=0..k-1 */

  MakeHash(H,G->n,k,HASH);
  FREQMINMAX(F,k,H,G->n,"the hash");
  DEBUG(PRINTT(H,G->n);PRINTT(F,k););
  printf("- time to compute stats and hash: %s\n",TopChrono(1));
  BARRE;
  
  /************************/
  /* calcule les tables B */
  /************************/
  
  /*
    Tables seulement définies pour les sommets u qui ne sont pas dans
    C. Il s'agit dans un premier temps des sommets à distance < r de u
    où r est la distance du plus proche sommet de C. Cependant, si le
    plus proche sommet est center (cela ne peut pas se produire si C
    contient le centre et tous ces voisins, en particulier si
    VARIANT=1), il convient alors de vérifier si à distance r-1 de u
    il n'y a pas un sommet de C. Si oui, il faut alors prendre r-1 et
    non r. Eventuellement tree(center) peut-être modifié. Le sommet u
    est supprimé de B[u]. Dans un second temps on ajoute tous les
    voisins de u à B[u] si |B[u]| était vide, c'est-à-dire si u était
    voisin de C. Le choix par défaut (->vpd) est le voisin de u allant
    de vers le sommet de B[u] le plus proche du center (dans
    tree(center)).
  */

  TopChrono(1);
  ALLOCZ(B,G->n,NULL); /* tableau de n tables B (vides au départ) */
  X=new_param_bfs(); /* pour ne pas le faire à chaque sommet u */
  X->clean=1; /* initialisation complète de X->D, puis partielle */

  int bmax,bsum; /* stat pour les "vraie" boules */
  bmax=bsum=0;
  
  for(u=0;u<G->n;u++){ /* pour chaque sommet u du graphe */

    /* Cherche d'abord le rayon r de B[u] pour calculer ensuite plus
       rapidement B[u] avec un BFS de profondeur r. Pour cela on
       remonte dans tree(center) jusqu'à trouvé un sommet de
       C. Cependant, ce n'est pas forcément le plus proche de u. Après
       le BFS depuis u, il faut alors vérifier si dans la dernière
       couche du BFS on a pas un sommet de C. Si oui, on décrémente le
       rayon et on modifie tree(center) de sorte à trouver un sommet
       de C plutôt et éviter de futurs problèmes. */

    t=v=u; r=-1;
    while(X0->P[v]!=-2){ /* si v pas dans C */
      v=X0->P[t=v]; /* t=v et v=parent(t) dans l'arbre X0->P */
      r++;
    }
    /* Ici:
       v=1er sommet ancêtre dans tree(center) de u qui est dans C
       t=dernier sommet de B[u] avant C (=fils de v dans tree(center))
       r=dist(u,t)=dist(u,v)-1
    */
    if(r<0) continue; /* rien à faire si on était parti d'un sommet de C */
    
    B[u]=new_table(0); /* on crée la table B[u] vide, NB: B[u]->n=0 */
    if(VARIANT==2){ B[u]->vpd=v; continue; }
      
    if(r>0){ /* ne rien faire ici si r=0 (=si u voisin de C) */
      X->hmax=r; bfs(G,u,X); /* calcule un BFS depuis u de profondeur r */
      i=B[u]->n=(X->n)-1; /* taille de B[u] sans u, si r est correct */
      if(v==center){ /* il faut vérifier si u n'a pas de sommet v dans C à distance r */
	/* on part des sommets i les plus éloignés de u, donc à distance r */
	/* NB: il y a toujours dans X->D au moins un sommet à distance < r, c'est u */
	while((X->D[X->file[i]]==r)&&(X0->P[X->file[i]]!=-2)) i--;
	if(X0->P[X->file[i]]==-2){ /* on a trouvé un sommet dans C à distance r */
	  v=X->file[i]; /* v=sommet dans C le plus proche de u */
	  t=X->P[v]; /* t=dernier sommet de B[u] avant C (=père de v dans le BFS de u) */
	  X0->P[t]=v; /* modifie le père de t dans tree(center) pour accélérér les tests suivants */
	  while(X->D[X->file[i]]==r) i--; /* on peut effacer les sommets de i+1 à X->n-1 */
	  B[u]->n=i; /* la taille de B[u] est i+1 avec u, donc i sans */
	  r--; /* r=dist(u,v)-1=dist(u,t) */
	}
      }
    }
    
    /* Ici:
       v=le plus proche sommet de u qui est dans C
       t=dernier sommet de B[u] avant C (=fils de v dans tree(center))
       r=dist(u,t)=dist(u,v)-1>=0
       B[u]->n=(taille de boule de rayon r sans u)-1
     */

    /* pour les stats avant ajoût des voisins et optimisation */
    bmax=max(bmax,1+B[u]->n);
    bsum += (1+B[u]->n);

    if(r==0){ /* si u voisin de C, alors B[u]=tous les voisins de u (sans u donc) */
      B[u]->vpd=v; /* v=voisin de u (dans C) */
      B[u]->n=G->d[u];
      B[u]->radius=1; /* rayon de B[u] */
      ALLOCZ(B[u]->node,B[u]->n,G->L[u][_i]); /* copie les voisins dans B[u]->node */
      ALLOCZ(B[u]->dist,B[u]->n,1); /* tous à distance 1 */
      continue; /* on passe au sommet u suivant */
    }

    /* ici r>0, et donc B[u]->n>0, et t<>u */
    /* on cherche le voisin par défaut de u, cad le fils de u dans
       BFS(u) qui mène à t. Là, on pourrait tester le père dans X0
       pour éventuellement le corriger et éviter les problèmes de
       correction de r en r-1 comme plus haut. */
    while(X->P[t]!=u) t=X->P[t];
    B[u]->vpd=t; /* voisin de u par défaut, ne peut pas être dans C */
    B[u]->radius=r; /* rayon de la boule */
    
    /* copie les sommets (et leur distance) du BFS couvrant B[u] en supprimant u */
    ALLOCZ(B[u]->node,B[u]->n,X->file[_i+1]); /* B[u]->node[i]=i-ème sommet de B[u] */
    ALLOCZ(B[u]->dist,B[u]->n,X->D[B[u]->node[_i]]); /* B[u]->dist[i]=distance du i-ème à u */
  }
  printf("- time to construct tables B: %s\n",TopChrono(1));
  printf("- maximum non-cluster ball size: %i\n",bmax);
  printf("- average non-cluster ball size: %.2lf (%i/%i)\n",
	 (double)bsum/(double)G->n,bsum,G->n);
  MINMAXMOY(B[_i]->n,G->n,B[_i],"table B size");

  /*************************/
  /* optimise les tables B */
  /*************************/
  
  /*
    Dans un premier temps on met des -1 aux sommets à supprimer.
    Ensuite on réorganise les tables B en gardant l'ordre des sommets
    restant. NB: on peut avoir B[u]<>NULL et B[u]->n=0. B[u]->radius
    n'est pas mis à jour et représente toujours la distance entre u et
    le sommet de B[u] le plus loin avant "clean". Il est important
    qu'avant l'appel, les sommets de B[u] soient rangés dans l'ordre
    du bfs(), donc non triées.
  */

  X->clean=1; /* initialisation complète de X->D, puis partielle */
  for(u=0;u<G->n;u++){
    if(B[u]){ /* il faut que la table B[u] existe */
      /* si deg[u]<=1 */
      if(G->d[u]<2) t=0;
      else{
	/* si deg(u)>1 */
	v=B[u]->vpd; /* v=voisin par défaut */
	X->hmax=(B[u]->radius)-1;
	bfs(G,v,X); /* calcule un BFS depuis v de profondeur rayon de B[u] -1 */
	for(i=0;i<X->n;i++){ /* passe en revue les sommets du bfs(v...) */
	  t=X->file[i]; /* t=sommet du bfs(v,...) */
	  r=SetSearch(t,B[u]->node,B[u]->n,0); /* r=indice tq B[u]->node[r]=t */
	  if((r>=0)&&(B[u]->dist[r]==1+X->D[t])) B[u]->node[r]=-1; /* supprime t */
	}
	/* supprime les -1 de B[u] en décalant B[u]->node et B[u]->dist */
	for(i=t=0;i<B[u]->n;i++)
	  if(B[u]->node[i]>=0){
	    B[u]->node[t]=B[u]->node[i];
	    B[u]->dist[t]=B[u]->dist[i];
	    t++;
	  }
      }
      /* redimensionne les tables à t = nouvelle taille de B[u] */
      B[u]->n=t;
      REALLOC(B[u]->node,t);
      REALLOC(B[u]->dist,t);
    }
  }
  printf("- time to clean tables B: %s\n",TopChrono(1));
  MINMAXMOY(B[_i]->n,G->n,B[_i],"table B size");

  DEBUG(for(u=0;u<G->n;u++)
	  if(B[u]){
	    printf("u=%i B[u]->vpd=%i B[u]->radius=%i ",u,B[u]->vpd,B[u]->radius);
	    if(B[u]) PRINTT(B[u]->node,B[u]->n);
	  });
  
  /************************/
  /* calcule les tables S */
  /************************/

  /*
    Tables seulement définies pour les sommets u de C sauf center.
    S[u] est la liste des sommets de hash i où u=C[i] et qui sont à
    distance au plus r=min(2,logn/loglogn) de u. On ne met ni u ni
    center dans S[u]. NB: Il est possible d'avoir dans S[u] des
    sommets de C (et même tous sauf u et center).
  */
  
  r=max(2,lg(G->n));
  r=ceil((double)r/(double)lg(r));
  const int depthS=max(2,r); /* peut-être > hauteur(tree(center)) */

  ALLOCZ(S,G->n,NULL); /* tableau de n tables S (vides au départ) */
  X->hmax=depthS; /* profondeur max du BFS */
  X->clean=1; /* initialisation complète de X->D, puis partielle */

  for(i=0;i<k;i++){ /* pour chaque sommet u du cluster sauf center */
    u=C[i]; if(u==center) continue;
    bfs(G,u,X); /* calcule un BFS depuis u<>center de profondeur depthS */
    t=min(X->n-1,F[i]-(H[u]==i)-(H[center]==i)); /* t=nombre max de sommets dans S[u] */
    S[u]=new_table(t); /* table S pour u */
    S[u]->n=0; /* on se sert de S[u]->n comme indice */
    /* on va parcourir les sommets du bfs() de profondeur <= depthS et
       ne garder que ceux de hash = i (ici i = hash(u) = hash(C[i]) */
    S[u]->vpd=-1;
    /* S[u]->vpd est l'indice (dans S[u]->node) du 1er sommet v de
       S[u] qui est à une distance >= depthS-1 de u.  Dans S[u], les
       sommets sont rangés selon le bfs(), donc en fonction de leur
       distance à u. On s'en sert plus tard pour accélérer la
       suppression dans R[u] des sommets v déjà dans S[u] (lire
       commentaires sur le remplissage de R[u]). NB: si aucun sommet
       n'est de hash = i (mauvaise fonction de hashage !), alors
       S[u]->vpd<0. */
    for(j=1;(j<X->n)&&(S[u]->n<t);j++){ /* pour chaque v du bfs(), v<>u */
      /* lorsque S[u] contient t sommets, on peut s'arrêter. */
      v=X->file[j]; /* v=j-ème sommet de la file, v<>u */
      if((H[v]==i)&&(v!=center)){ /* v a le bon hash et v<>center */
	S[u]->node[S[u]->n]=v;
	S[u]->dist[S[u]->n]=X->D[v];
	if((X->D[v]>=depthS-1)&&(S[u]->vpd<0))
	  S[u]->vpd=S[u]->n; /* 1er sommet de hash i à distance >= depthS-1 */
	S[u]->n++; /* un sommet de plus dans S[u] */
      }
    }
    if(S[u]->n){ /* redimensionne les tables de S[u] pour pas gaspiller */
      REALLOC(S[u]->node,S[u]->n);
      REALLOC(S[u]->dist,S[u]->n);
    }else{ /* si aucun sommet dans S[u], on efface S[u] */
      free_table(S[u]);
      S[u]=NULL;
    }
  }
  free_param_bfs(X); /* plus de BFS à faire, X ne sert plus à rien */
  printf("- time to construct tables S: %s (radius %i)\n",TopChrono(1),r);
  
  DEBUG(
	for(i=0;i<k;i++){
	  u=C[i];
	  printf("u=%i hash=%i |S[u]|=%i depth=%i\n",u,i,S[u]?S[u]->n:0,depthS);
	  if(S[u]){
	    PRINTT(S[u]->node,S[u]->n);   
	    PRINTT(S[u]->dist,S[u]->n);
	  }
	}
	);
  
  /************************/
  /* calcule les tables R */
  /************************/
  
  /*
    Tables seulement définies pour les sommets u de C avec u=C[i].
    R[u] est la liste des sommets dont le hash est i et qui ne sont
    pas déjà dans S.  On ne met ni u ni center dans R[u]. C'est donc
    comme S[u] mais sans la restriction de distance. Pour les sommets
    de R[u], on appliquera le routage dans tree(center).
  
    Les sommets de hash=couleur(center) sont tous dans R[center] (sauf
    center lui-même si hash(center)=couleur(center)). Si u<>center,
    alors pour qu'un sommet v (de hash i) soit dans R[u] il faut qu'il
    soit à une distance >= depthS de center. S'il est à distance =
    depthS-1 (ou inférieure), alors en passant par center on obtient
    une route de u à v de longueur <= (depthS-1)+1=depthS ce qui n'est
    pas mieux que d'utiliser S[u]. Le sommet v sera donc déjà dans
    S[u]. Maintenant, les sommets candidats v à distance >= depthS de
    center ne peuvent être qu'à distance au moins depthS-1 de u. S'ils
    étaient à distance < depthS-1, alors ils seraient à distance <
    depthS de center. Or, ils sont déjà à distance >= depthS. Donc
    pour chercher si v est déjà dans S[u], on peut se contenter de
    vérifier les sommets de S[u] à distance >= depthS-1 de u. Ces
    sommets candidats sont stockés à partir de l'indice S[u]->vpd dans
    S[u]->node (car S n'est pas encore retriée). Dans R[u]->dist on
    stocke la distance réalisée dans l'arbre pour faire u->v.
  */
  
  ALLOCZ(R,G->n,NULL); /* tableau de tables R (vides au départ) */
  for(i=0;i<k;i++){ /* seulement pour les k sommets u de C */
    u=C[i]; /* u=sommet de C */
    t=F[i]-(H[u]==i); /* t=#sommets à mettre dans R[u] (moins éventuellement u) */
    if((u!=center)&&(H[center]==i)) t--; /* enlève u et aussi center si H[center]=i */
    if(S[u]) t-=S[u]->n; /* on enlève ceux déjà dans S[u] */
    if(t>0){ /* si table pas vide */
      R[u]=new_table(t); /* table R pour u */
      R[u]->n=0; /* on se sert de R[u]->n comme indice */
    }
  }
  free(F); /* ne sert plus à rien */
  /* on remplit les tables R[u]->node et R[u]->dist */
  /* R[u]->dist[i]=dist. dans tree(center) entre u et le i-ème sommet de R[u] */
  for(v=0;v<G->n;v++){ /* écrit chaque sommet v dans la bonne table R[u] */
    u=C[H[v]]; /* u=sommet de C en charge du hash de v (donc de couleur H[v]) */
    if((u==v)||(v==center)) continue; /* on ne met ni u ni center dans R[u] */
    if((X0->D[v]<depthS)&&(u!=center)) continue; /* seulement pour v à dist > depthS-1 de center */
    /* On cherche si v n'est pas déjà dans S[u]. Il suffit de chercher
       des sommets v à distance >= depthS-1 de u */ 
    if(S[u]){ /* cherche v dans S[u], si elle n'est pas vide bien sûr */
      j=S[u]->vpd; /* 1er sommet à distance >= depthS-1 de u */
      if(j<0) j=S[u]->n; /* ne rien faire si pas de 1er sommet (est-ce nécessaire ?) */
      for(;j<S[u]->n;j++) if(S[u]->node[j]==v) break; /* v est dans S[u] */
    }
    if((S[u]==NULL)||(j==S[u]->n)){ /* v n'a pas été trouvé dans S[u] */
      R[u]->node[R[u]->n]=v; /* ajoute v à R[u] */
      // Calcule R[u]->dist: on remonte de v jusqu'au cluster et on
      // voit si l'on passe par u ou pas. La distance d de v à u dans
      // tree(center) est la suivante:
      //   si u=center, alors          d=X0->D[v]
      //   si u<>center et t=u, alors  d=X0->D[v]-1
      //   si u<>center et t<>u, alors d=X0->D[v]+1
      if(u==center) /* remontée inutile si u est le center */
	R[u]->dist[R[u]->n]=X0->D[v];
      else{ /* on remonte à partir de v jusqu'à t=sommet dans C */
	t=v; while(X0->P[t]!=-2) t=X0->P[t]; /* si t pas dans C, t=parent(t) */
	R[u]->dist[R[u]->n]=X0->D[v]+((t==u)?-1:1);
      }
      R[u]->n++; /* un sommet de plus dans R[u] */
    }
  }
  printf("- time to construct tables R: %s\n",TopChrono(1));
  free_param_bfs(X0); /* ne sert plus à rien */

  DEBUG(
	for(i=0;i<k;i++){
	  u=C[i];
	  printf("u=%i hash=%i |R[u]|=%i\n",u,i,R[u]?R[u]->n:0);
	  if(R[u]) PRINTT(R[u]->node,R[u]->n);   
	  if(R[u]) PRINTT(R[u]->dist,R[u]->n);
	}
	);
    
  /************************/
  /* calcule les tables W */
  /************************/

  /*
    Tables seulement définies pour tous les sommets u de C, W[u] donne
    la liste des couleurs (=indice dans C) des voisins de u dans C. La
    couleur du center n'est jamais ajouté à W[u], car c'est le port
    par défaut. NB: Le nom du sommet de couleur i est C[i]. On ne se
    sert pas de W[u]->dist. Si W[u] est vide, alors on ne la supprime
    pas de façon à avoir le port par défaut. Dans le cas VARIANT=1,
    toutes les tables sont vides (routage dans l'étoile sans table
    pour le center et les feuilles).
  */

  ALLOCZ(W,G->n,NULL); /* tableau des tables W (vides au départ) */
  for(i=0;i<k;i++){ /* seulement pour les sommets u de C, y compris center */
    u=C[i]; /* u=sommet de C de couleur i */
    W[u]=new_table(0); /* table W pour u, NB: W[u]->n=0 */
    ALLOC(W[u]->node,k-1); /* au plus k-1 voisins */
    W[u]->vpd=center; /* voisin par défaut = center */
    if(VARIANT==1) continue; /* toutes les tables seront vides */
    for(j=0;j<G->d[u];j++){ /* pour chaque voisin de u */
      v=G->L[u][j]; /* v=j-ème voisin de u */
      if(B[v]||(v==center)) continue; /* on veut v dans C et v<>center */
      W[u]->node[W[u]->n++]=SetSearch(v,C,k,1); /* ajoute la couleur de v à W[u]->node */
    }
    if(W[u]->n) /* réajuste la taille */
      REALLOC(W[u]->node,W[u]->n);
  }
  printf("- time to construct tables W: %s\n",TopChrono(1));

  DEBUG(
	for(i=0;i<k;i++){
	  u=C[i];
	  printf("u=%i hash=%i |W[u]|=%i\n",u,i,W[u]?W[u]->n:0);
	  if(W[u]) PRINTT(W[u]->node,W[u]->n);   
	}
	);

  /*******************/
  /* trie les tables */
  /*******************/
  for(u=0;u<G->n;u++){
    SortTable(B[u]);
    SortTable(R[u]);
    SortTable(S[u]);
    SortTable(W[u]);
  }
  printf("- time to sort tables B,R,S,W: %s\n",TopChrono(1));
  BARRE;
  
  /* calcule de la taille Z[u] des tables de chaque sommet u */
  ALLOCZ(Z,G->n,1); /* taille=1 au moins pour chaque sommet u */
  for(u=0;u<G->n;u++){ /* pour chaque sommet u */
    if(B[u]) Z[u] += B[u]->n;
    if(S[u]) Z[u] += S[u]->n;
    if(R[u]) Z[u] += R[u]->n;
    if(W[u]) Z[u] += W[u]->n;
  }

  /* affiche taille min/max et moyenne des différentes tables */
  MINMAXMOY(B[_i]->n,G->n,B[_i],"table B size");
  MINMAXMOY(S[C[_i]]->n,k,S[C[_i]],"table S size");
  MINMAXMOY(R[C[_i]]->n,k,R[C[_i]],"table R size");
  MINMAXMOY(W[C[_i]]->n,k,W[C[_i]],"table W size");

  /* affiche la distribution des tailles de table */
  PrintDistrib(Z,G->n,10,"routing table size");
  free(Z); /* ne sert plus à rien */
  BARRE;

  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_cluster_tables,RT,1);
  RT->B=B;
  RT->S=S;
  RT->R=R;
  RT->W=W;
  RT->C=C; // besoin pour le routage avec W
  RT->H=H; // besoin du hash des sommets
  RT->n=G->n; // besoin pour libérer les n tables
  RT->center=center; // besoin pour VARIANT=1

  return RT;
}


int rs_cluster_length(int u,int v,rs_cluster_tables *X){
/*
  Renvoie le nombre de sauts du routage selon les tables générées par
  rs_cluster() pour router un message de u à v, ou bien -1 si la route
  n'a pu être déterminée. Dans X on a toutes les tables nécessaire au
  schéma, notamment les tables B,S,R,W. Si u<0 alors on teste la
  validité des tables (et on renvoie une valeur non-nulle en cas
  d'erreur).

  Amélioration possible: Lorsque u a des voisins dans C, alors, plutôt
  que d'aller vers ->vdp, choisir un landmark de C parmi un ensemble
  prédéterminé de p = ceil(2m/n) = O(1) sommets (le degré moyen)
  d'éventuellement la bonne couleur h(v). Aller vers ->vdp seulement
  si cette recherche à échoué. Pour le schéma théorique, on peut
  découper la table B en deux, B1 et B2, où B2 serait les sommets
  voisins de u dans C, et doubler B2 (double accès par sommet et par
  couleur) de façon à garantir un temps constant dans tous les cas, et
  pas seulement un temps égale au degré moyen (garanti seulement avec
  grande proba dans les RPLG.
*/

  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->S==NULL) return 1;
    if(X->R==NULL) return 1;
    if(X->W==NULL) return 1;
    if(X->H==NULL) return 1;
    if(X->C==NULL) return 1;
    return 0;
  }
  
  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivé
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }
  int i;

  // si u n'est pas dans C
  
  if(X->B[u]){
    // v dans B[u] ?
    i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
    // si v dans B[u]:u -> v
    if(i>=0){ DEBUG(printf("in B[u]\n");); return X->B[u]->dist[i];}
    // si v pas dans B[u]: u -> B[u]->vpd -> v
    DEBUG(printf("via B[u]->vpd=%i\n",X->B[u]->vpd););
    return 1 + rs_cluster_length(X->B[u]->vpd,v,X);
  }
  
  // si u est dans C

  if(X->S[u]){ // il faut S non vide, donc u<>center
    // v dans S[u] ?
    i=SetSearch(v,X->S[u]->node,X->S[u]->n,1);
    // si v est dans S
    if(i>=0){ DEBUG(printf("in S[u]\n");); return X->S[u]->dist[i];}
  }

  // si v n'est pas dans S
  // v dans R[u] ?

  if(X->R[u]){ // il faut R non vide
    i=SetSearch(v,X->R[u]->node,X->R[u]->n,1);
    // si v est dans R, on cherche si v descendant de u dans tree(center)
    if(i>=0){ DEBUG(printf("in R[u]\n");); return X->R[u]->dist[i]; }
  }
  
  // si v n'est pas ni dans R ni dans S
  // H[v] dans W[u] ?

  if(X->W[u]){ /* ne peut pas être vide */
    i=SetSearch(X->H[v],X->W[u]->node,X->W[u]->n,1);
    // si H[v] est dans W[u]: u -> C[H[v]] -> v
    // si H[v] pas dans W[u]: u -> center -> v
    if((i>=0)||((VARIANT==1)&&(u==X->center))) u=X->C[X->H[v]];
    else u=X->W[u]->vpd;
    DEBUG(
	  if(u==X->C[X->H[v]])
	    printf("%s: u -> C[H[v]]=%i\n",(VARIANT)?"variant 1":"H[v] in W[u]",u);
	  else
	    printf("%s: u -> W[u]->vpd=%i\n",(VARIANT)?"variant 1":"H[v] not in W[u]",u);
	  );
    return 1 + rs_cluster_length(u,v,X);
  }

  // y'a un problème
  DEBUG(printf("fail: W[u] does not exist\n"););
  return -1;
}


rs_dcr_tables *rs_dcr(graph *G,int k){
/*
  Calcule pour le graphe G le routing scheme dcr de paramètre k>0
  correspondant au nombre de couleurs. On renvoie les tables
  calculées. Le stretch est toujours <= 5.

  Si une couleur n'existe pas, alors cela marche quand même: les
  boules couvrent tour le graphe et le stretch est alors 1.
*/

  int u,v,w,i,j,t,c,nc,nl;

  printf("\nDCR\n");
  BARRE;
  printf("- wanted number of colors: %i\n",k);
  TopChrono(1); /* reset du chrono */

  /**********************************/
  /* calcule C & H, couleur et hash */
  /**********************************/

  /* C[u]=0..k-1, couleur du sommet u */
  NALLOCZ(int,C,G->n,random()%k);

  /* H[u]=0..k-1, hash du sommet u */
  int *H=MakeHash(NULL,G->n,k,HASH);

  DEBUG(PRINTT(C,G->n);PRINTT(H,G->n););

  /* affichage de stats sur C & H */
  int *F;
  FREQMINMAX(F,k,C,G->n,"the colors");
  nl=F[0]; /* nl=nombre de landmarks, ceux de couleurs 0 */
  for(i=nc=0;i<k;i++) nc += (F[i]>0); /* nc=nombre de couleurs différentes, nc<=k */
  free(F); /* nécessaire car réalloué par le prochain FREQMINMAX */
  FREQMINMAX(F,k,H,G->n,"the hash"); /* fréquence des hashs (sert pour taille des tables) */

  printf("- real number of colors: %i\n",nc);
  printf("- number of landmarks: %i\n",nl);
  printf("- time to compute hash, color & stats: %s\n",TopChrono(1));
  BARRE;

  /*********************************************/
  /* calcule S, la liste des bfs des landmarks */
  /*********************************************/

  /* S[u]=bfs(u,G), pour u landmark, NULL sinon */
  
  NALLOC(param_bfs*,S,G->n);
  for(u=0;u<G->n;u++)
    if(C[u]) S[u]=NULL;
    else{
      S[u]=bfs(G,u,NULL);
      free(S[u]->file); /* libère les pointeurs inutilisés */
      S[u]->file=NULL; /* important pour le free_param_bfs() plus tard */
    }
  printf("- time to construct landmark bfs (array S): %s\n",TopChrono(1));

  DEBUG(
	for(u=0;u<G->n;u++)
	  if(C[u]==0){
	    PRINT(u);
	    PRINTT(S[u]->P,G->n);
	    PRINTT(S[u]->D,G->n);
	    printf("\n");
	  }
	);
  
  /*****************************/
  /* calcule les tables B et W */
  /*****************************/

  /*
    B[u]=boule de voisinage de u, la plus "petite" contenant toutes
         les couleurs (liste des sommets et leur distance à u), la
         dernière couronne étant ordonnée selon les identifiants des
         sommets.

    B[u]->node[i]=i-ème sommet de la boule de u
    B[u]->dist[i]=distance entre u et B[u]->node[i]
    B[u]->vpd=plus proche landmark de u (plus proche sommet de couleur 0 dans B[u])
    W[u]->node[c]=next-hop vers le 1er sommet de couleur c de B[u]->node
    (la table W, indexée par les couleurs, sert pour avoir le temps constant)

    Algorithme:

      On fait un bfs(G,u) par couche (avec cont=1). On trie chaque
      couche par identifiant (si on a visité au moins assez de
      sommets), et on regarde à chaque fois le premier moment où l'on
      voit une couleur donnée (remplissage de W). Au passage on
      détermine, lorsque la couleur est 0, le plus proche landmark t
      de u. On s'arrête quand on a visité toutes les couleurs ou bien
      que tous les sommets du graphe ont été visités (cela peut
      arriver si une couleur n'est portée par aucun sommet).
   */

  NALLOCZ(table*,B,G->n,NULL);  // B=tableau de n tables B (vides au départ)
  NALLOCZ(table*,W,G->n,NULL);  // W=tableau de n tables W (vides au départ)
  param_bfs *X=new_param_bfs(); // X=le résultat du bfs() depuis u
  NALLOC(int,T,G->n); /* T=tableau pour la dernière couche de sommets */
  X->clean=1; /* initialisation complète des distances, puis partielle */

  for(u=0;u<G->n;u++){ /* calcule B[u] & W[u] pour chaque u */
    c=0; /* c=nombre de couleurs (différentes) déjà rencontrées */
    X->cont=1; /* évite de recalculer le début de l'arbre */
    X->hmax=0; /* au départ on visite seulement u */
    W[u]=new_table(0); /* W[u]->node[c]=-1 si la couleur c n'a pas été visitée */
    ALLOCZ(W[u]->node,k,-1); /* on n'utilise pas ->dist */
    W[u]->n=k; /* taille de W[u] */
    do{
      
      bfs(G,u,X); /* on parcoure jusqu'à distance hmax de u */
      for(j=0,i=X->tf;i<X->n;i++) T[j++]=X->file[i]; /* copie la dernière couche dans T */
      if(c+j>=nc) qsort(T,j,sizeof(int),fcmp_int); /* trie la dernière couche (=T) */
      /* on ne trie T que si on a le potentiel pour avoir toutes les
	 couleurs, c'est-à-dire si le nombre de sommets de la dernière
	 couche (=j) + le nombre de couleurs déjà rencontrées (=c) est
	 au moins nc */

      DEBUG(
	    PRINT(u);PRINT(X->n);PRINT(X->tf);
	    PRINTT(T,j);PRINT(c);
	    printf("\n");
	    );
      
      /* remplit W en parcourant les sommets de la dernière couche */
      /* le plus proche landmark de u, lorsque rencontré, est stocké dans t */
      for(i=0;(i<j)&&(c<k);i++){
	v=T[i]; /* v=i-ème sommet de la dernière couche */
	if(W[u]->node[C[v]]>=0) continue; /* couleur déjà rencontrée */
	/* ici la couleur C[v] n'a jamais été rencontrée */
	w=v; /* si v=u, alors il faut effacer le -1 dans W[u] */
	if(v!=u) while(X->P[w]!=u) w=X->P[w]; /* cherche le next-hop w depuis u vers v */
	W[u]->node[C[v]]=w; /* met le next-hop dans W[u] */
	if(C[v]==0) t=v; /* si v est un landmark, alors c'est le +proche de u */
	c++; /* une couleur de plus */
      }
      X->hmax++; /* couche suivante */

      /* On sort de la boucle si: soit on a toutes les couleurs dans
	 W[u]->node, ou bien on a visité tout le graphe, ce qui est
	 possible si toutes les couleurs n'étaient pas représentées.
	 NB: ici, dans tous les cas, i est nombre de sommets de la
	 dernière couche à recopier partiellement. */
      
    }while((c<k)&&(X->n<G->n));
    
    /* on construit B[u] à partir de X->file et de T */
    /* attention ! on ne met pas u dans B[u] */
    B[u]=new_table((X->tf)+i-1); /* table de B[u] (sans u) */
    B[u]->vpd=t; /* t=plus proche landmark de u */

    j=0; /* indice pour B[u]->node[] */
    t=1; /* indice pour X->file[], on saute u */
    for(;t<X->tf;j++){ /* copie X->file sauf la dernière couche */
      v=X->file[t++];
      B[u]->node[j]=v;
      B[u]->dist[j]=X->D[v];
    }
    for(t=0;t<i;j++){ /* copie la partie traitée de la dernière couche */
      v=T[t++];
      B[u]->node[j]=v;
      B[u]->dist[j]=X->D[v];
    }
   }
  printf("- time to construct tables B & W: %s\n",TopChrono(1));
  free_param_bfs(X);
  free(T); /* ne sert plus à rien */

  /* tri des tables B, important pour le routage */
  for(u=0;u<G->n;u++) SortTable(B[u]);
  printf("- time to sort tables B: %s\n",TopChrono(1));
  BARRE;
  
  DEBUG(
	for(u=0;u<G->n;u++){
	  PRINT(u);
	  PRINT(B[u]->vpd);
	  PRINTT(B[u]->node,B[u]->n);
	  PRINTT(W[u]->node,W[u]->n);
	  printf("\n");
	}
	);
  
  /* taille totale des tables */
  NALLOCZ(int,Z,G->n,1); /* taille=1 au moins pour chaque sommet u */
  for(u=0;u<G->n;u++){ /* pour chaque sommet u */
    Z[u] += B[u]->n; /* taille de B */
    Z[u] += nl-(C[u]==0); /* table des landmarks */
    Z[u] += k-1; /* taille de W (rien pour la couleur 0) */
    Z[u] += F[C[u]]-(C[u]==H[u]); /* table des sommets dont le hash est la couleur de u */
  } /* taille théorique: k*H(k) + n/k + k + n/k, voir aussi func1() */

  /* affiche taille min/max et moyenne des différentes tables */
  MINMAXMOY(B[_i]->n,G->n,1,"table B size");
  MINMAXMOY(nl-(C[_i]==0),G->n,1,"landmark table size");
  MINMAXMOY(k-1,G->n,1,"table W size");
  MINMAXMOY(F[C[_i]]-(C[_i]==H[_i]),G->n,1,"own color/hash table size");
  free(F); /* ne sert plus à rien */

  /* affiche la distribution des tailles de table */
  PrintDistrib(Z,G->n,10,"routing table size");
  free(Z); /* ne sert plus à rien */
  printf("- theoretical average: 2√(n*ln(n*ln(n))) = %i\n",
	 (int)(2.0*sqrt(G->n*log(G->n*log(G->n)))));
  BARRE;
  
  /* assemble les tables en une seule pour le retour */  
  NALLOC(rs_dcr_tables,RT,1);
  RT->B=B;
  RT->W=W;
  RT->S=S;
  RT->H=H;
  RT->C=C;
  RT->n=G->n;
  ALLOCZ(RT->dist,G->n,S[_i]?S[_i]->D:NULL); /* distance partielle */

  return RT;
}


rs_tzrplg_tables *rs_tzrplg(graph *G,double t){
/*
  Calcule les tables de routage selon le schéma tz_rplg pour le graphe
  G, une adaptation du schéma de Thorup et Zwick avec les landmarks
  sur les sommets de plus haut degré.

  Cet algorithme est spécialisé pour les graphes de type RPLG(n,t). Le
  paramètre t attendu (power-law exponent) est celui du graphe G
  fournit. Les performances de ce schéma sont meilleures si le
  paramètre t est le bon, mais marche quel que soit t>1.5. Les valeurs
  de t et de VARIANT permettent de calculer le nombre de landmarks.
  L'ensemble des landmarks est aussi nommé "core" du graphe par Chung
  & Lu.

  On calcule deux tables: B et L. La table B[u], pour chaque sommet u,
  contient tous les sommets strictement plus proche que son plus
  proche landmark. La table L[u], pour chaque sommet u, contient le
  next-hop vers de u vers le landmark le plus proche de u selon
  l'arbre de +cc enraciné dans le landmark.

  Rem: dans [CSTW12], pour les BC graphs, pour les moyennes ils
  prennent probablement n=10K et non n=+grande composante connexe
  (7K), ce qui change les choses.
*/

  printf("\nTZ RPLG\n");
  BARRE;

  TopChrono(1); /* reset du chrono */

  double gamma; /* gamma = paramètre fonction de t */
  int core_size; /* core_size = nombre de landmarks */
  int i,u;

  /*** Calcul des landmarks ***/
  /* C[i]=liste des landmarks, i=0..core_size-1 */

  // calcule D = liste triée (ordre croissant) des degrés des sommets
  // ou permutation aléatoire des sommets (si VARIANT=2)
  int *D;
  if(VARIANT<2) D=SortInt(G->d,NULL,G->n,1,NULL,SORT_INDEXi);
  if(VARIANT==2){ // permutation aléatoire
    ALLOCZ(D,G->n,_i);
    Permute(D,G->n);
  }

  // calcule core_size, dépend de VARIANT et de t
  if((VARIANT==0)&&(t>1.5)){
    gamma=(double)(t-2.0)/(double)(2.0*t-3.0);
    core_size=ceil(pow(G->n,gamma));
  }
  if((VARIANT==1)&&(t>1.5)){
    // C = { u : deg(u)>(n^gamma')/4 }, gamma'=(1-gamma)/(t-1)=1/(2t-3) 
    gamma=1.0/(double)(2.0*t-3.0);
    u=ceil(pow(G->n,gamma));
    // cherche le 1er sommet de degré <= u
    i=G->n-1; // part de la fin (haut degré)
    while((i>=0)&&(G->d[D[i--]]>u));
    core_size=i/4;
  }
  if((VARIANT==2)&&(t>0.0)) core_size=(int)t;
  if(t<0.0) core_size=(int)(-t);
  if(t==0.0) core_size=ceil(sqrt(G->n));

  core_size=max(core_size,0);    // core_size >= 0
  core_size=min(core_size,G->n); // core_size <= G->n

  // C=liste des landmarks
  NALLOCZ(int,C,core_size,D[G->n-_i-1]); // de +haut degré (ou sommets aléatoires)
  
  // affiche C et sa densité
  for(i=u=0;i<core_size;i++) u += G->d[C[i]];
  printf("- core degree: ");
  APERCU(G->d[D[G->n-_i-1]],core_size,10,2);
  free(D); /* ne sert plus à rien */
  printf("- core size: %i",core_size);
  if((VARIANT==0)&&(t>1.5)) printf(" (n^%g)",gamma);
  if((VARIANT==1)&&(t>1.5)) printf(" (deg>n^%g)",gamma);
  if(t==0.0) printf(" (sqrt(n))");
  printf("\n");
  printf("- sum of cluster's degrees: %i",u);
  if(u*100>G->n) printf(" (%.2lfn)",(double)u/(double)G->n);
  printf("\n- time to construct core: %s\n",TopChrono(1));
  
  /*** Construction des tables L ***/

  // L[u]->node[i] = identifiant du i-ème landmark
  // L[u]->dist[i] = distance entre u et L[u]->node[i]
  // lmin[u] = distance entre u et son plus proche landmark
  // label[u] = indice i du landmark le plus proche de u

  NALLOCZ(table*,L,G->n,new_table(core_size)); // n tables L de taille core_size
  NALLOCZ(int,lmin,G->n,G->n); // par défaut lmin[u]=n
  NALLOC(int,label,G->n);
  
  int l;
  param_bfs *X;
  
  for(i=0;i<core_size;i++){ // on fait un bfs() pour chaque landmark
    l=C[i]; // l=i-ème landmark
    X=bfs(G,l,NULL); // bfs() depuis l
    for(u=0;u<G->n;u++){
      L[u]->node[i]=l;
      L[u]->dist[i]=X->D[u];  // dist(u,l)
      if(X->D[u]<lmin[u]){ // landmark plus proche ?
	lmin[u]=X->D[u]; // NB: le rayon de B[u] sera lmin[u]-1
	label[u]=i; /* NB: L[u]->node[label[u]]=landmark le plus proche de u */ 
      }
    }
    free_param_bfs(X);
  }
  free(C);
  printf("- time to construct tables L: %s\n",TopChrono(1));
  
  /*** Construction des tables B ***/
  
  NALLOCZ(table*,B,G->n,NULL); // tableau de n tables B (vides au départ)
  X=new_param_bfs(); // pour bfs() depuis u
  X->clean=1; // initialisation complète des distances, puis partielle
  
  /* construit la boule B[u] de rayon lmin[u]-1 avec bfs(u,.) */

  for(u=0;u<G->n;u++){
    X->hmax=max(0,lmin[u]-1); // la boule contient les sommets
			      // strictement plus proche que le coeur
    B[u]=new_table(0);
    // si B n'est pas vide, alors on fait un BFS pour calculer B
    if(X->hmax){
      bfs(G,u,X);
      B[u]->n=(X->n)-1; // taille de B[u] (sans u)
      /* copie les sommets (et leur distance) du bfs() en supprimant u */
      ALLOCZ(B[u]->node,B[u]->n,X->file[_i+1]); /* =i-ème sommet de B[u] */
      ALLOCZ(B[u]->dist,B[u]->n,X->D[B[u]->node[_i]]); /* =distance du i-ème à u */
    }
    // si B est vide alors:
    else B[u]->n=0;
  }
  free_param_bfs(X);
  free(lmin); /* ne sert plus à rien */
  
  printf("- time to construct tables B: %s\n",TopChrono(1));
  
  /* tri des tables B */
  for(u=0;u<G->n;u++) SortTable(B[u]);
  printf("- time to sort tables B: %s\n",TopChrono(1));
  BARRE;

  /* taille totale des tables */
  NALLOCZ(int,Z,G->n,1); /* taille=1 au moins pour chaque sommet u */
  for(u=0;u<G->n;u++){ /* pour chaque sommet u */
    if(B[u]) Z[u] += B[u]->n;
    if(L[u]) Z[u] += L[u]->n;
  }
  
  /* affiche taille min/max et moyenne des différentes tables */
  MINMAXMOY(B[_i]->n,G->n,B[_i],"table B size");
  MINMAXMOY(L[_i]->n,G->n,L[_i],"table L size");

  /* affiche la distribution des tailles de table */
  PrintDistrib(Z,G->n,10,"routing table size");
  free(Z); /* ne sert plus à rien */
  BARRE;
  
  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_tzrplg_tables,RT,1);
  RT->B=B;
  RT->L=L;
  RT->label=label;
  RT->n=G->n; // besoin pour libérér les n tables
  
  return RT;
}


int rs_tzrplg_length(int u,int v,rs_tzrplg_tables *X){
/*
  Renvoie le nombre de sauts du routage selon les tables générées par
  rs_tzrplg() pour router un message de u à v, ou bien -1 si la route
  n'a pu être déterminée. Dans X on a toutes les tables nécessaire au
  schéma, notamment les tables B et L. Si u<0 alors on teste la
  validité des tables (et on renvoie une valeur non-nulle en cas
  d'erreur).
*/

  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->L==NULL) return 1;
    return 0;
  }

  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivé
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }

  // routage dans la boule de u
  // v dans B[u] ?

  if(X->B[u]){
    int i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
    // si v dans B[u]: u -> v
    if(i>=0){
      DEBUG(
	    printf("in B[u], distance %i: (", X->B[u]->dist[i]);
	    int _v;
	    for(_v=0;_v<X->B[u]->n;_v++)
	      printf("%i,", X->B[u]->node[_v]);
	    printf(")\n");
	    );
      return X->B[u]->dist[i];
    }
  }

  // routage via le plus proche landmark de v
  // route: u -> lv -> v (aucun raccourci, voir Algo. 1 dans CSTW12)

  if(X->L[u]){
    // On doit récupérer 2 entrées: L[u][lv] et L[v][lv]. Pour cela on
    // doit d'abord trouver l'identité de lv qui se trouve dans
    // l'étiquette de v:
    int lv = X->label[v];
    if(X->L[v])
        return X->L[u]->dist[lv] + X->L[v]->dist[lv];
  }

  // y'a un problème, pas d'entrée pour v
  DEBUG(printf("fail: no route found when routing from %i to %i \n",u,v););
  return -1;
}


rs_bc_tables *rs_bc(graph *G, int k){
/*
  Schéma de routage selon Brady-Cowen 2006. Ici k est le rayon du Core
  (=C), la zone "dense" du graphe G. Le stretch est additif est <=
  2k. Dans l'article d'origine, k=d/2, mais prendre un rayon d/2 quand
  d est impair n'est pas facile.

  Principe: On construit un arbre BFS (=T) enraciné dans la partie
  dense de G, un peu comme rs_cluster() à partir d'un sommet de grand
  degré. Le Core (=C) est la boule de rayon k depuis la racine de
  T. Le Fringe de G est F=G\C, le graphe induit par les sommets à
  distance >= k de r, ceux qui ne sont donc pas dans C.

  On construit une liste (=L) de BFS couvrant G ainsi qu'une forêt
  (=H) de BFS de G comme suit. Au départ, L={T}, et H est la forêt T
  restreinte à F, c'est-à-dire H=T\C. Puis, pour chaque arête {u,v} de
  F\H, on vérifie si l'ajoût de {u,v} à H crée un cycle ou pas (se
  servir de FindSet() et de la détection de cycle comme dans
  Minor()). Si c'est non, on met à jour la forêt H en lui ajoutant
  {u,v}. Si c'est oui, on calcule un BFS de racine u (ou v) qu'on
  ajoute à L. En pratique il sera sans doute intéressant de rendre
  asymétrique le graphe de façon à ne considérer qu'une seule fois
  chaque des arêtes (voir HalfGraph()).

  Une fois toutes les arêtes {u,v} de F ainsi balayées, on cacule pour
  chaque composantes connexes de H un BFS qu'on ajoute à la liste
  L. Un sommet de F est donc traversé par l'arbre T, un arbre de H, et
  autant arbres que d'arêtes qui ont crées un cycle dans l'étape
  précédante. Pour un sommet de C, c'est un de moins, H ne couvrant
  pas C.

  L'algorithme de routage de u à v consiste simplement à router dans
  l'arbre A de L contenant u et v et qui minimise dist_A(u,v),
  distance qu'on peut calculer avec nca_bfs().

  En fait, tous les arbres de L\H couvrent G, et donc contiennent u et
  v. Dans la pratique, on peut donc vérifier dans un premier temps si
  u et v sont dans la même composante de H et le cas échéant calculer
  une première estimation de la distance entre u et v, puis de ne plus
  vérifier aucun arbre de H. Le plus pratique est donc de gérer
  séparément T, H et la liste des autres BFS.
*/

  int u,v,center;

  printf("\nBRADY-COWEN\n");
  BARRE;
  
  TopChrono(1); /* reset du chrono */
  
  /* trouve un sommet center de degré max, v=deg(center) */
  for(u=v=center=0;u<G->n;u++)
    if(G->d[u]>v) v=G->d[center=u];
  printf("- degree of the center: %i (id:%i)\n",v,center);


  /* tri des tables B */
  ;;;
  BARRE;

  /* taille totale des tables */
  NALLOCZ(int,Z,G->n,1); /* taille=1 au moins pour chaque sommet u */
  ;;;
  
  /* affiche taille min/max et moyenne des différentes tables */
  ;;;

  /* affiche la distribution des tailles de table */
  ;;;
  free(Z); /* ne sert plus à rien */
  BARRE;
  
  /* assemble les tables en une seule pour le retour */
  NALLOC(rs_bc_tables,RT,1);
  RT->n=G->n; // besoin pour libérér les n tables
  
  return RT;
}


int rs_bc_length(int u,int v,rs_bc_tables *X){
/*
  Renvoie le nombre de sauts du routage selon les tables générées par
  rs_bc() pour router un message de u à v, ou bien -1 si la route n'a
  pu être déterminée. Dans X on a toutes les tables nécessaire au
  schéma, notamment les tables ... Si u<0 alors on teste la validité
  des tables (et on renvoie une valeur non-nulle en cas d'erreur).
*/
  
  return -1;
}


int nca_bfs(int u,int v,param_bfs *X){
/*
  Calcule le plus petit ancêtre commun entre u et v dans l'arbre donné
  par le bfs X. Algorithme (de complexité dist(u,v)): on remonte
  d'abord le sommet le plus profond. Lorsque les deux sommets sont à
  la même profondeur on teste si c'est les mêmes et si oui on a
  trouver le nca, ou alors on les remontes tous les deux.
*/
  if(X->D[u]<X->D[v]){ int w; SWAP(u,v,w); }
  
  /* ici u est plus profond que v (ou à même profondeur) */
  while(X->D[u]!=X->D[v]) u=X->P[u];

  /* ici u et v sont à la même profondeur */
  while(u!=v) u=X->P[u],v=X->P[v];
	
  return u;
}


int rs_dcr_length_rec(int u,int v,int w,int a,rs_dcr_tables *X){
/*
  Fonction récursive donnant la longueur de la route de u à v avec les
  tables X initialisées par dcr, et l'en-tête (w,a):
  - w=landmark intermédiaire (s ou t, w=s par défaut)
  - a=nca(u,v) dans l'arbre du landmark w (a<0 par défaut)
*/
  
  DEBUG(printf("  \nu=%i v=%i: ",u,v););

  // on est arrivé
  if(u==v){ DEBUG(printf("u=v\n");); return 0; }
  int i;

  // v dans B[u] ? (NB: B[u] existe toujours)
  
  i=SetSearch(v,X->B[u]->node,X->B[u]->n,1);
  // si v dans B[u]:u -> v
  if(i>=0){ DEBUG(printf("in B[u]\n");); return X->B[u]->dist[i];}
  
  // v n'est pas dans B[u]
  // est-ce que v est dans la table des landmarks ?

  if(X->C[v]==0){ // 0=couleur des landmarks
    DEBUG(printf("%i is a landmark\n",v););
    param_bfs *Y=X->S[v];
    return Y->D[u]; // dist(u,v), v=landmark
  }

  // v n'est pas landmark
  // est-on arrivé au responsable de v ?
  
  if((a<0)&&(X->C[u]!=X->H[v])){ // si pas encore au responsable
    u=X->W[u]->node[X->H[v]]; // next-hop de u vers le responsable de v
    DEBUG(printf("go to %i the closest node of color %i, the hash of %i (a=%i and w=%i)",u,X->H[v],v,a,w););
    return 1 + rs_dcr_length_rec(u,v,w,a,X);
  }

  // on est passé par (ou on est sur) le responsable de v ?
  
  if(a<0){ // on est arrivé au responsable de v
    int w1=X->B[w]->vpd; // landmark le +proche de la source (=w si a<0)
    int w2=X->B[v]->vpd; // landmark le +proche de cible (=v)
    param_bfs *Y1=X->S[w1]; // l'arbre de w1
    param_bfs *Y2=X->S[w2]; // l'arbre de w2
    int a1=nca_bfs(u,v,Y1); // ancêtre commun entre u et v dans Y1
    int a2=nca_bfs(u,v,Y2); // ancêtre commun entre u et v dans Y2
    int d1=Y1->D[u]+Y1->D[v]-2*Y1->D[a1]; // d1=dist(u,v) dans Y1
    int d2=Y2->D[u]+Y2->D[v]-2*Y2->D[a2]; // d2=dist(u,v) dans Y2
    if(d1<d2) w=w1,a=a1; else w=w2,a=a2;
    
    /* la ligne suivante est une optimisation: si d1=d2, alors on a
       intérêt de choisir l'ancêtre a1 ou a2 le plus loin de u cela
       laisse plus de chances à l'algo de court-circuiter le routage
       dans l'arbre. */
    if((d1==d2)&&(Y1->D[u]-Y1->D[a1]>=Y2->D[u]-Y2->D[a2])) w=w1,a=a1;
    
    DEBUG(
	  printf("node in charge of %i reached\n",v);
	  printf("go to a=%i, the nca of %i and %i in the tree of landmark w=%i\n",a,u,v,w);
	  );
  }

  // ici on a:
  //  w=landmark intermédiaire (landmark de s ou de t)
  //  a=ancêtre commun entre u et v dans l'arbre de w
  if(u==a){ // si u est arrivé au bon ancêtre
    DEBUG(printf("ancestor a=%i reached\n",a););
    return X->S[w]->D[v] - X->S[w]->D[a];
  }

  // ici il faudrait ajouter les racourcis, si un des ancêtres de v
  // dans l'arbre de w (stockés dans l'étiquette de v) est dans B[u]
  
  u=X->S[w]->P[u]; // on remonte dans l'arbre de w
  DEBUG(printf("go up to a=%i in the tree of landmark w=%i",a,w););
  return 1 + rs_dcr_length_rec(u,v,w,a,X);
}


int rs_dcr_length(int u,int v,rs_dcr_tables *X){
/*
  Renvoie le nombre de sauts du routage selon les tables générées par
  rs_dcr() pour router un message de u à v, ou bien -1 si la route n'a
  pu être déterminée. Dans X on a toutes les tables nécessaires au
  schéma. Si u<0, on réalise quelques tests de bases sur les tables
  X. La fonction fait essentiellement appel à une fonction recursive
  où la gestion d'un en-tête est nécessaire.
*/

  if(u<0){
    if(X==NULL) return 1;
    if(X->B==NULL) return 1;
    if(X->W==NULL) return 1;
    if(X->S==NULL) return 1;
    if(X->H==NULL) return 1;
    if(X->C==NULL) return 1;
    return 0;
  }

  return rs_dcr_length_rec(u,v,u,-1,X);
}


enum{
  SC_NONE,
  SC_ALL,
  SC_ONE,
  SC_NPAIRS,
  SC_PAIR,
  SC_EDGES,
  SC_UV
};


static inline long route_uv(graph *G,
			    int u,int v,
			    int h,int hmax,
			    int **dist,long **stat,int *L,int *M,param_bfs *X){
/*
  Remplit la table de statistiques stat[][] et aussi de distance
  dist[][] des sommets testés avec h=longueur du routage de u à
  v. Dans le scenario SC_EDGES les distances ne sont pas mise à jour
  (la distance étant toujours 1). Renvoie une valeur non-nulle si une
  erreur est survenue, 0 sinon.
*/
  DEBUG(printf("  #hop=%i\n",h););
  if((h<0)||(h>hmax)) return 1; /* y'a un problème, sommet suivant */ 

  int j,k,l,t;

  /* calcule k=dist(u,v) */
  if(SCENARIO.mode==SC_EDGES) k=1; /* pas de bfs() dans ce cas */
  else{ /* est-ce que dist(u,v) est déjà connue ? */
    if((dist[u]==NULL)&&(dist[v]==NULL)){
      bfs(G,u,X);
      ALLOCZ(dist[u],G->n,X->D[_i]); /* alloue dist[u] et y copie X->D[] */
    }
    k=dist[u]? dist[u][v] : dist[v][u]; /* k=dist(u,v) */
  }
  
  if(h<k) return 1; /* il y'a un problème */

  /* ici il faut faire +1 dans stat[k][j] où j=h-k */
  j=h-k;
  if(j>=L[k]){ /* tableau stat[k] n'est pas assez grand */
    l=L[k]; /* sauvegarde L[k] */
    L[k]=max(j+1,l)<<1; /* on double la taille courante ou de h-k+1 */
    if(stat[k]==NULL) ALLOCZ(stat[k],L[k],0L); /* première fois */
    else{
      REALLOC(stat[k],L[k]); /* aggrandit stat[k] */
      for(t=l;t<L[k];t++) stat[k][t]=0L; /* initialise la fin du tableau */
    }
  }
  stat[k][j]++; /* une route de longueur h=k+j de plus */
  M[k]=max(M[k],j); /* détour max pour la distance k */

  return 0; /* tout c'est bien passé */
}


int pgcd(int a,int b){
/*
  Renvoie le plus grand commun diviseur de a et b.
*/
  int c;
  while(a){ c=a; a=b%a; b=c; }
  return b;
}


int routing_test(graph *G,void *T,rt_length length,int hmax,int **distp){
/*
  Teste un scenario de routage pour le graphe G avec les tables de
  routage T et la fonction de longueur length(), puis affiche les
  distributions des longueurs de route, de distance et de stretch. Le
  scenario est décrit par la variable globale SCENARIO. La matrice
  distp[][] est une matrice partielle de distance du graphe
  (éventuellement calculées lors de la construction de T), ce qui peut
  accélérer les tests.

  Plus précisément, pour tout sommet u de G, distp[u], si non NULL,
  doit être la distance de u vers tous les sommets de G. Attention !
  un vecteur partiel n'est pas autorisé. Il ne faut donc pas utilisé
  de bfs() partiel par exemple. Il est possible d'avoir distp[u]=NULL
  (matrice partielle <> vecteur partiel). Si distp=NULL, alors aucune
  distance n'est pré-calculées. Le tableau distp n'est pas
  modifiée. Seules les vecteurs de distance calculés par
  routing_test() sont libérés, charge à l'appelant de libérer les
  vecteurs de distp.

  La fonction renvoie une valeur non-nulle si une erreur s'est
  produite, et 0 sinon. La valeur de hmax indique la longueur maximum
  de routage. Cela permet de détecter des routes trop longues pour
  être correctes (et éviter une boucle infinie). Si hmax<0, alors on
  prend hmax=2n comme valeur par défaut.
*/
  if(SCENARIO.mode==SC_NONE) return 0;
  
  TopChrono(1);
  printf("\nROUTING TEST\n");
  BARRE;

  /* vérification de la table T */
  if(length(-1,0,T)) return 1;

  int n=G->n;
  int u,v,i,j,k,t;
  int dmax; /* plus longue distance */
  int lmax; /* plus longue route */
  double x; /* pour les affichage de distributions */

  /* type "long" pour le comptage */
  long p,s;
  long err=0L; /* nb de routage erronés */

  param_bfs *X=new_param_bfs(); /* pour le calcul de la distance */
  X->clean=1; /* pour appels multiples */

  /* dist[u][v]=distance entre u et v=0..n-1, n valeur différentes au plus */
  /* on a dist[u]=NULL si on a pas encore calculé de bfs(u,...) */
  NALLOCZ(int*,dist,n,distp?distp[_i]:NULL); /* dist[u]=distp[u] ou NULL */

  /* stat[k][j]=#routage entre sommets à distance k et de longueur de k+j (j=détour) */
  NALLOCZ(long*,stat,n,NULL); /* stat[k]=NULL si pas encore eut de stat pour dist k */
  /* L[k]=taille du tableau stat[k] */
  NALLOCZ(int,L,n,0); /* taille nulle par défaut */
  /* M[k]=détour maximum (=longueur-k) rencontrée pour la distance k, M[k]<= L[k] */
  NALLOCZ(int,M,n,-1); /* par défaut M[k]<0 */
  if(hmax<0) hmax=n<<1;

  switch(SCENARIO.mode){
    
  case SC_ALL:
    p=(long)n*(long)(n-1); /* cast "(long)" important */
    if(2*lg(n)>(int)(8*sizeof(p)-1)) Erreur(33); /* dépassement arithmétique */
    printf("- all-to-all pairs: %s pairs ...",millier(p));
    fflush(stdout);
    for(u=0;u<n;u++)
      for(v=0;v<n;v++)
	if(u!=v) err += route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X);
    break;

  case SC_NPAIRS:
    printf("- n random pairs: %s pair%s ...",millier(n),PLURIEL(n));
    fflush(stdout);
    for(i=0;i<n;i++){
      u=random()%n;
      v=random()%n;
      if(u!=v) err += route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X);
    }
    break;
  
  case SC_PAIR:
    p=SCENARIO.u;
    printf("- random pairs: %s pair%s ...",millier(p),PLURIEL(p));
    fflush(stdout);
    for(i=0;i<SCENARIO.u;i++){
      u=random()%n;
      v=random()%n;
      if(u!=v) err += route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X);
    }
    break;
  
  case SC_EDGES:
    p=NbEdges(G)<<1;
    printf("- all neighbor pairs: %s pair%s ...",millier(p),PLURIEL(p));
    fflush(stdout);
    for(u=0;u<n;u++){
      p=G->d[u];
      for(i=0;i<p;i++){
	v=G->L[u][i];
        err += route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X);
      }
    }
    break;
  
  case SC_ONE:
    u=SCENARIO.u;
    p=n-1;
    if(u<0) u=random()%n;
    printf("- one-to-all pairs from %i: %s pair%s ...",u,millier(p),PLURIEL(p));
    fflush(stdout);
    for(v=0;v<n;v++)
      if(u!=v) err += route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X);
    break;
  
  case SC_UV:
    u=SCENARIO.u;
    v=SCENARIO.v;
    if(u<0) u=random()%n;
    if(v<0) v=random()%n;
    printf("- routing from %i to %i: 1 pair ...",u,v);
    fflush(stdout);
    err += route_uv(G,u,v,length(u,v,T),hmax,dist,stat,L,M,X);
    break;

  default:
    Erreur(23);
  }
  printf(" Ok (%s)\n",TopChrono(1));

  /* libère les tableaux devenus inutiles: X et L */
  free_param_bfs(X);
  free(L);
  
  /* on ne libère que les distances calculées par routing_test() */
  for(u=0;u<G->n;u++) if((distp==NULL)||(distp[u]==NULL)) free(dist[u]);
  free(dist);

  /*
    Affiche la distribution:
    - des longueurs de route,
    - des distances,
    - des stretch.
  */

  /* détermine dmax (=distance max) et lmax (=longueur routage max) */
  dmax=lmax=0;
  for(k=0;k<n;k++)
    if(stat[k]){
      dmax=max(dmax,k);
      lmax=max(lmax,k+M[k]);
  }

  /* pour un petit gain mémoire: réajuste les tableaux dépendant de la distance k */
  /* NB: le tableau L n'existe plus */
  REALLOC(stat,dmax+1); /* stat[0..dmax] */
  REALLOC(M,dmax+1);    /*    M[0..dmax] */

  /* calcule le nombre p de routages non erronés */
  for(k=p=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	p += stat[k][j];

  printf("- #failed routings: %s",millier(err));
  if(err) printf(" (%i%%)",p?(int)ceil((100.0*err)/(double)(p+err)):100);
  printf("\n");
  if(p==0L){
    printf("  all routings failed!\n");
    goto fin_routing_test;
  }

  BARRE;
  printf("- route length distribution:\n");

  /* calcule F[i]=nombre de routages de longueur i */
  NALLOCZ(long,F,lmax+1,0L);
  for(k=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	F[k+j] += stat[k][j];
  
  /* s=somme totale des longueurs */
  for(i=0,s=0L;i<=lmax;i++) s += (long)i*F[i];

  for(i=0;i<=lmax;i++){
    if(F[i]==0L) continue;
    x=(double)F[i]/(double)p;
    printf("    %i \t%02i%% ",i,(int)(100.0*x));
    RULING(x);
    printf(" [%li] \n",F[i]);
  }
  printf("- average route length: %.02lf (%li/%li)\n",(double)s/(double)p,s,p);
  printf("- maximum route length: %i\n",lmax);
  free(F);

  BARRE;
  printf("- distance distribution:\n");
  ALLOCZ(F,dmax+1,0L); /* F[k]=nombre de distances de longueur k */
  for(k=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	F[k] += stat[k][j];

  /* s=somme totale des distances */
  for(i=0,s=0L;i<=dmax;i++) s += (long)i*F[i];

  for(i=0;i<=dmax;i++){
    if(F[i]==0) continue;
    x=(double)F[i]/(double)p;
    printf("    %i \t%02i%% ",i,(int)(100.0*x));
    RULING(x);
    printf(" [%li] \n",F[i]);
  }
  printf("- average distance: %.02lf (%li/%li)\n",(double)s/(double)p,s,p);
  printf("- maximum distance: %i\n",dmax);
  free(F);

  BARRE;
  printf("- stretch distribution:\n");

  /* t=borne sup sur le nombre de stretch différents */
  for(k=t=0;k<=dmax;k++)
    if(stat[k]) t += M[k]+1;

  triplet *ptr;
  triplet e={0,0,0L}; /* triplet nul */
  NALLOCZ(triplet,P,t,e); /* P=tableau de triplets (j,k,stat[k][j]) */

  /* on construit P */
  /* t=nombre de triplets ajoutés à P */
  for(k=t=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++)
	if(stat[k][j]){
	  if(k){ u=pgcd(k,j); e.x=j/u; e.y=k/u; }
	  else{ e.x=0; e.y=1; } /* pour dist=0, le stretch est 1+0 */
	  e.z=stat[k][j];
	  ptr=bsearch(&e,P,t,sizeof(triplet),fcmp_stretch);
	  if(ptr) P[(int)(ptr-P)].z += e.z;  /* e est déjà dans P, ajoute e.z */
	  else{
	    P[t++]=e; /* on ajoute une nouvelle entrée e à P */
	    qsort(P,t,sizeof(triplet),fcmp_stretch); /* trie le nouveau tableau */
	  }
	}

  /* affiche les stretchs contenus dans P */
  for(i=0;i<t;i++){
    x=(double)(P[i].z)/(double)p;
    printf("  %0.3lf (%i/%i)",1.0+(double)P[i].x/(double)P[i].y,P[i].x+P[i].y,P[i].y);
    printf("\t%02i%% ",(int)(100.0*x));
    RULING(x);
    printf(" [%li] \n",P[i].z);
  }
  free(P);

  /* stretch moyen et stretch max */
  /* on pourrait l'avoir directement avec P, mais avec stat[][] on
     peut aussi avoir la distance max qui atteint le stretch max (on a
     perdu l'info à cause du pgcd). */
  double r=0.0; /* r=somme total des stretch */
  double smax=1.0; /* stretch=max */

  u=v=1;
  for(k=0;k<=dmax;k++)
    if(stat[k])
      for(j=0;j<=M[k];j++){
	x=1.0; if(k) x += (double)j/(double)k; /* stretch=1.0 si k=0 */
	r += (double)stat[k][j]*x;
	if(x>=smax){ smax=x; u=j; v=k; }
      }
  printf("- average stretch: %.03lf (%.02lf/%li)\n",(double)r/(double)p,r,p);
  printf("- maximum stretch: %.02lf (%i/%i)\n",1.0+(double)u/(double)v,u+v,v);

 fin_routing_test:
  BARRE;
  free(M);
  FREE2(stat,dmax+1);
  return 0;
}


/***********************************

       FIN ROUTINES POUR LES
          ROUTING SCHEMES

***********************************/


void RemoveVertex(graph *G,int u)
/*
  Supprime toutes les occurences du sommet u dans le graphe G, sans
  renuméroter les sommets et sans réallouées les listes. Après cet
  appel, les listes de G ne sont plus forcément triées, mais la liste
  de G->L[u] existe toujours. Simplement G->d[u] est à zéro.

  Effet de bord: modifie G->sort, G->d[.], G->L[.]. Attention !  G->n
  n'est pas modifié.
*/
{
  int i,j,v,du=G->d[u],dv;
  G->d[u]=0;
  for(i=0;i<du;i++){
    v=G->L[u][i];
    dv=G->d[v];
    for(j=0;j<dv;j++)
      if(G->L[v][j]==u)	G->L[v][j]=G->L[v][--dv];
    G->d[v]=dv;
  }
  G->sort=0;
  return;
}


int AdjGraph(graph *G,int u,int v)
/*
  Renvoie 1 ssi dans le graphe G le sommet v apparaît dans la liste
  d'adjacence de u. Si les listes du graphe G sont est triées (G->sort
  est vrai), une dichotomie est effectuée (donc en log(deg(u))) sinon
  c'est un parcours séquentiel (en deg(u)).
*/
{
  /* recherche dichotomique si le graphe est trié */
  if(G->sort)
    return (bsearch(&v,G->L[u],G->d[u],sizeof(int),fcmp_int)!=NULL);

  /* sinon, recherche séquentielle */
  int i;
  for(i=G->d[u];i>0;)
    if(G->L[u][--i]==v) return 1; /* on a trouvé v dans la liste de u */
  return 0;
}


int Treewidth(graph *H,int code)
/*
  Calcul approché ou excate de la treewidth de H. Si code=0, on
  calcule une borne supérieure (heuristique degmin). Si code=1 on
  calcule la valeur exacte. H n'est pas modifié. Dans H->int1 on
  renvoie le nb de tests effectués.
*/
{
  if(H==NULL) return 0;
  int n,j,t,tw,l,v,n1;
  int d,r,i,u,k;
  int d0,r0,i0,u0;
  graph *G;

  H->int1=0;
  n=H->n;
  n1=n-1;
  tw=(code==1)? Treewidth(H,0) : min(n1,NbEdges(H));

  /* tw=upper bound sur tw. On a intérêt d'avoir la valeur la plus
     faible possible, car on va alors éliminer les mauvais ordres
     d'éliminations plus rapidement.

    tw max en fonction du nb m d'arêtes:
    m=0 => tw=0
    m=1 => tw=1
    m=2 => tw=1
    m=3 => tw<=2
    m=4 => tw<=2
    m=5 => tw<=2?
    m=6 => tw<=3
    ...
   */

  NALLOCZ(int,P,n,_i); /* permutation initiales des sommets */
  NALLOCZ(int,D,n,n1);

  do{
    G=ExtractSubgraph(H,NULL,0,0); /* on copie H */
    GraphRealloc(G,D); /* plonge G dans le graphe complet */
    k=0; /* tw par défaut si G sans d'arêtes */

    for(j=0;j<n1;j++){ /* n-1 fois */

      H->int1++;
      u0=P[i0=j];
      d0=r0=G->d[u0];
      if(code==0){
	for(i=j;i<n1;i++){
	  u=P[i]; d=G->d[u];
	  if(d<=d0){
	    /* calcule r=nb d'arêtes dans N(u) */
	    for(r=l=0;l<d;l++)
	      for(t=l+1;t<d;t++)
		r += AdjGraph(G,G->L[u][l],G->L[u][t]);
	    if((d<d0)||(r<r0)){ d0=d;r0=r;i0=i;u0=u; }
	  }
	} /* fin for(i=...) */
	P[i0]=P[j]; /* décale P[i], que si code=0 */
      }
      k=max(k,d0); /* met à jour tw */
      if(k>=tw) goto nextP; /* on ne fera pas moins */
      RemoveVertex(G,u0); /* supprime u */
      /* remplace N(u) par une clique */
      for(i=0;i<d0;i++){
	u=G->L[u0][i];
	for(t=i+1;t<d0;t++){
	  v=G->L[u0][t];
	  if(!AdjGraph(G,u,v)) ADD_EDGE(G,u,v);
	}
      }
      
    } /* fin for(j=...) */

    tw=min(tw,k);
  nextP:
    free_graph(G);
    if((code==0)||(tw<2)) break; /* si tw(G)=0 ou 1, alors on ne
				    trouvera pas plus petit */
  }while(!NextPermutation(P,n,NULL));

  free(D);
  free(P);
  return tw;
}


int *Prune(graph *G,int *z)
/*
  Algorithme permettant de supprimer les sommets du graphe G par
  dégrés minimum. Un ordre d'élémination des sommets est renvoyé sous
  la forme d'un tableau de taille n. Plus précisément, on renvoie un
  tableau de sommets T de sorte que u=T[i] est un sommet de degré
  minimum du graphe G\(T[0]...T[i-1]). La complexité est linéaire,
  i.e. O(n+m). Si *z<>NULL, on retourne dans z le degré maximum
  rencontré donc la dégénérescence du graphe.
*/
{
  int i,j,u,v,d,k,p,r,t,c,n1,n=G->n;

  /*
    1. On construit les tableaux suivants:

    T[0..n[ = tableau de sommets triés par degré croissant (le résultat)
    R[0..n[ = tableau de positions dans T des sommets
    D[0..n[ = tableau des degrés des sommets
    P[0..n[ = tableau des positions dans T du premier sommet de degré donné

    2. On traite les sommets dans l'ordre T[0],T[1],... Supposons que
    u=T[i]. Alors, pour chaque voisin v dans G, encore existant (donc
    situé après i dans le tableau T), on met à jour la structure
    T,R,D,P.
   */

  /* initialise T,R,D,P */
  if(n<1) return NULL;
  NALLOCZ(int,D,n,G->d[_i]);
  int *R=SortInt(D,NULL,n,0,NULL,SORT_INDEXe);
  NALLOC(int,T,n); for(u=0;u<n;u++) T[R[u]]=u;
  NALLOCZ(int,P,n,-1);
  for(i=n1=n-1;i>=0;i--) P[D[T[i]]]=i; /* i=n-1,...,0 */

  for(c=i=0;i<n1;i++){ /* pour chaque sommet, pas besoin de faire le dernier */
    u=T[i]; /* u = sommet de degré (D[u]) minimum */
    if(D[u]>c) c=D[u]; /* mémorise le plus grand degré rencontré */
    d=G->d[u]; /* d = deg(u) */
    for(j=0;j<d;j++){ /* pour chaque voisin */
      v=G->L[u][j]; /* v=voisin de u */
      r=R[v]; /* v=T[r]; */
      if(r>i){ /* mettre à jour que si v existe encore */
	k=D[v]--; /* le degré de v diminue, k=D[v] juste avant */
	p=P[k--]++; /* on va échanger v et T[p]. Les sommets de degré k commencent juste après */
	if(P[k]<0) P[k]=p; /* y'a un sommet de degré k en plus */
	if(p>i){ /* si p est avant i, ne rien faire. Cel arrive que si k=0 */
	  t=T[p]; /* t = 1er sommet de de degré k */
	  T[p]=v; /* dans T, on avance v en position p (à la place de t donc) */
	  T[r]=t; /* on met t à la place de v */
	  R[v]=p; /* on met à jour la nouvelle position de v dans T */
	  R[t]=r; /* on met à jour la nouvelle position de t dans T */
	}
      }
    }
  }

  if(z!=NULL) *z=c; /* retourne la dégénérescence */
  free(R);
  free(D);
  free(P);
  return T;
}


int *GreedyColor(graph *G,int *R)
/*
  Algorithme permettant de colorier de manière gloutonne un graphe G à
  n sommets. La complexité en temps est linéaire, i.e. O(n+m). Il faut
  G non NULL. On renvoie un tableau C[0..n[ où C[u] est la couleur du
  sommet u, un entier entre 0 et n-1. Le tableau R[0..n[ donne un
  ordre (inverse) de visite des sommets: On colorie le sommet u avec
  la plus petite couleur possible dans le graphe induit par les
  sommets situé après u dans R (on commence donc avec R[n-1]). Si
  R=NULL, alors l'ordre R[i]=i est utilisé. On renvoie dans G->int1 la
  couleur maximum utilisée (donc le nb de couleurs-1). Cette valeur
  vaut -1 si n<1.  */
{
  int i,j,u,d,l,c,n=G->n;

  /*
    On utilise 3 tableaux:

    C[0..n[ = tableau final des couleurs, au départ =-1
    L[0..n-1[ = liste de couleurs utilisées par le sommet courant
    M[0..n-1[, M[i]=1 ssi la couleur i est utilisée par un voisin du
    sommet courant, au départ =0. On met une sentiennelle à M si bien
    que toujours on M[n-1]=0.

  */

  G->int1=-1;
  if(n<1) return NULL;
  NALLOCZ(int,C,n,-1);
  NALLOCZ(int,M,n,0);
  i=n-1;
  NALLOC(int,L,i);

  for(;i>=0;i--){ /* en partant de la fin */
    u=(R==NULL)? i : R[i];
    d=G->d[u]; /* d = deg(u) */

    /* on liste dans L et M les couleurs rencontrées */
    for(j=l=0;j<d;j++){ /* pour chaque voisin v de u */
      c=C[G->L[u][j]]; /* c=couleur du voisin v */
      if(c>=0){ /* voisin déjà colorié ? */
	L[l++]=c; /* on ajoute la couleur c à la liste */
	M[c]=1; /* la couleur c n'est pas à prendre */
      }
    }
    
    /* on cherche la 1ère couleur à 1 (=non-rencontrée) */
    j=0; while(M[j]) j++; /* s'arrête toujours à cause de la sentiennelle */
    C[u]=j; /* j est la couleur recherchée pour u */
    G->int1=max(G->int1,j); /* couleur maximum rencontrée */

    /* il faut ré-initialiser rapidement M */
    for(j=0;j<l;j++) M[L[j]]=0; /* on efface les 1 qu'on à mis dans M */
  }
  
  free(L);
  free(M);
  return C;
}


void HalfGraph(graph *G){
/*
  Transforme G en un graphe asymétrique en ne gardant que les arcs
  u->v tels que v<u. En particulier les boucles sont supprimés. Les
  listes d'adjacence sont aussi triées.  Cela permet de faire des
  parcours de graphe deux fois plus rapidement, par exemple, pour
  vérifier qu'une coloration est propre.
*/
  if(G==NULL) return;
  if(!G->sym) return;

  int n=G->n,i,u,d;
  NALLOC(int,D,n); /* tableau des nouveaux degrés */

  SortGraph(G,0);
  for(u=0;u<n;u++){
    d=G->d[u];
    for(i=0;i<d;i++)
      if(G->L[u][i]>=u) break;
    D[u]=i; /* coupe la liste à i */
  }
  GraphRealloc(G,D); /* modifie toutes les listes d'adjacence */
  free(D);
  return;
}


void kColorSat(graph *G,int k){
/*
  Affiche les contraintes SAT pour la k-coloration de G au format
  Dimacs. Attention ! G est modifié (manque la moitié de ses arcs).

  Le nombre de variables est n*k.
  Le nombre de clause est n+m*k.
*/

  if(G==NULL) return;
  int n,c,u,i,d,m;
  
  n=G->n;
  m=NbEdges(G);
  printf("p cnf %i %i\n",n*k,n+m*k);

  /*
    Variables: chaque sommet i a k variables numérotés x(i,c) avec
    i=0..n-1 et c=0..k-1. Il faut x(i,c)>0. On pose x(i,c)=1+k*i+c.

    Contraintes sommets: il faut x(i,0) v ... v x(i,k-1), et la
    conjonction pour tous les sommets.

    Containtes arêtes: pour chaque arête {i,j} et couleur c il ne faut
    pas que x(i,c)=x(j,c). Autrement dit il faut -x(i,c) v -x(j,c). Il
    faut la conjonction de toutes les couleurs c et toutes les arêtes
    {i,j}.
  */

  /* liste chaque sommet */
  for(u=0;u<n;u++){
    for(c=0;c<k;c++) printf("%i ",1+u*k+c);
    printf("0\n");
  }

  HalfGraph(G); /* enlève la moitié des arcs */

  /* liste chaque arc */
  for(u=0;u<n;u++){
    d=G->d[u];
    for(i=0;i<d;i++)
      for(c=0;c<k;c++)
	printf("-%i -%i 0\n",1+u*k+c,1+G->L[u][i]*k+c);
  }
  
  return;
}


void kIndepSat(graph *G,int k){
/*
  Affiche les contraintes SAT pour ensemble indépendant de taille k
  pour le graphe G au format Dimacs. Attention ! G est modifié (manque
  la moitié de ses arcs).

  Variables:

    Pour chaque i=0..n-1:
    X(i)=1 ssi le sommet i est dans l'ensemble indépendant.
    Si i-j est une arête alors -X(i) v -X(j).
    NB: pour Vertex Cover, il faudrait suffit d'ajouter X(i) v X(j).

    Pour chaque t=0..n et b=0..k:
    Y(t,b)=1 ssi la somme des t variables X_0+...+X(t-1) est au
    moins b. On a:

       Y(n,k) = 1
       Y(t,0)=1 pour tout t>=0
       Y(t,b)=0 pour tout 0<=t<b et b>0
       Y(t,b) => Y(t-1,b) v (Y(t-1,b-1) ^ X(t-1))

       <=> -Y(t,b) v Y(t-1,b) v Y(t-1,b-1) ET
           -Y(t,b) v Y(t-1,b) v X(t-1)

    #variables X(i): n
    #variables Y(t,b): (n+1)*(k+1)
    #variables totales: n+(n+1)*(k+1)

*/

  if(G==NULL) return;
  int n,i,j,t,d,b,m;
  
  n=G->n;
  m=NbEdges(G);

  // il faut entrer le nombre exactes de variables et de clauses
  printf("p cnf %i %i\n",n+(n+1)*(k+1),m+n+2+k*(k+1)/2+k*(2*n+1-k));

  /* numéros des variables: attention de ne surtout pas utiliser la
     variable numéro 0 qui signifie "fin de ligne" */

#define X(i)   ((i)+1)             // numéro de la variable X(i)
#define Y(t,b) (n+1+(t)*(k+1)+(b)) // numéro de la variable Y(t,b)

  HalfGraph(G); /* rend le graphe simple et asymétrique */

  // pour chaque arêtes i-j: -X(i) v -X(j)
  // #clauses: m
  for(i=0;i<n;i++){
    d=G->d[i];
    for(j=0;j<d;j++)
      printf("-%i -%i 0\n",X(i),X(G->L[i][j]));
  }

  // Y(n,k)=1
  // #clause: 1
  printf("%i 0\n",Y(n,k));

  // cas b=0 et t>=0: Y(t,0)=1
  // #clauses: n+1
  for(t=0;t<=n;t++) printf("%i 0\n",Y(t,0));

  // cas b>=1 et 0<=t<b: Y(t,b)=0
  // #clauses: 1+2+3+...+k = k*(k+1)/2
  for(b=1;b<=k;b++)
    for(t=0;t<b;t++)
      printf("-%i 0\n",Y(t,b));

  // cas b>=1 et t>=b: récurrence
  // #clauses: 2*sum_{b=1}^k (n-b+1) = 2*(n+1)*k - k(k+1) = k*(2*n+1-k)
  for(b=1;b<=k;b++)
    for(t=b;t<=n;t++){
      printf("-%i %i %i 0\n",Y(t,b),Y(t-1,b),Y(t-1,b-1));
      printf("-%i %i %i 0\n",Y(t,b),Y(t-1,b),X(t-1));
    }

  return;
}


int *kColor(graph *G,int k){
/*
  Algorithme permettant de colorier en au plus k couleurs un graphe G,
  si c'est possible. La complexité est (n+m)*k^{n-1} dans le pire des
  cas.  On renvoie un tableau C[0..n[ où C[u] est la couleur du sommet
  u, un entier entre 0 et k-1. On renvoie NULL s'il n'est pas possible
  de colorier G en k couleurs, si G=NULL ou k<1. On renvoie dans
  G->int1 la couleur maximum utilisée (qui peut être < k-1). On
  utilise toutes les couleurs de [0,G->int1]. On symétrise le graphe,
  qui est donc modifié, afin d'enlever la moitié des arêtes à
  vérifier.

  La stratégie est la suivante. On part d'une coloration initiale des
  sommets C=[0,...,0,k-1] où la couleur du dernier sommet u=n-1 est
  fixée par C[n-1]=k-1. Puis on vérifie si C est propre ou non. Si
  c'est non, on incrémente C comme un compteur, et on recommence avec
  la coloration suivante. Ainsi toutes les colorations possibles de G
  sont passées en revue. On teste toujours en priorité la dernière
  arête qui a posé problème avant de vérifier tout le graphe.

  Optimisations possibles à faire:
  
  1. Réduction de données. On peut supprimer récursivement les sommets
     de degré < k. On pourra toujours les ajouter à la fin si la
     coloration a réussie.

  2. On peut décomposer le graphe en composantes connexes, il y a
     ainsi moins de colorations possibles à tester.

  3. On peut renuméroter les sommets selon un parcours BFS. De cette
     façon la vérification pourrait être accélérée lorsqu'on change
     une seule couleur.

*/

  int n,c,i,d,u,v,b,*T;
  if((G==NULL)||(k<1)) return NULL;
  HalfGraph(G); /* enlève la moitié des arêtes */
  n=G->n;

  NALLOCZ(int,C,n,0); /* C[u]=couleur du sommet u */
  C[n-1]=k-1; /* on peut supposer que le sommet n-1 a une couleur fixée */
  if(n<2) goto fin_kcolor; /* s'il y a un seul sommet */
  b=1; /* b=vrai ssi la dernière arête coloriée est propre */

  do{
    /* vérifie si C est une coloration propre */
    if(b){ /* on le fait que si b est vrai */
      for(u=0;u<n;u++){ /* pour chaque sommet */
	c=C[u]; d=G->d[u]; /* degré et couleur de u */
	for(i=0;i<d;i++) /* pour chaque voisin de i */
	  if(c==C[G->L[u][i]]) break; /* coloration pas propre */
	if(i<d) break; /* coloration pas propre */
      }
      if(u==n) goto fin_kcolor; /* la coloration est propre */
    }
    /* ici l'arête (u,i) n'est pas correctement coloriée */
    
    /* on change C en incrémentant le compteur */
    v=0;
  loop_kcolor:
    C[v]++;
    if(C[v]==k){
      C[v++]=0;
      if(v==n) break;
      goto loop_kcolor;
    }
    
    /* est-ce que l'arête (u,i) est toujours mal coloriée ?
       si oui, pas la peine de vérifier tout le graphe */
    b=(C[u]!=C[G->L[u][i]]); /* b=vrai ssi (u,i) est propre */

  }while(v<n);

  /* aucune coloration trouvée */
  free(C);
  return NULL;

 fin_kcolor:
  /* on a trouvé une coloration propre */
  /* on réduit l'espace des couleurs utilisées.  NB: ici on a encore
     C[n-1]=k-1 */

  ALLOCZ(T,k,-1);
  /* si T[c]=i, alors la couleur c se renumérote en i. Si c n'a jamais
     été rencontrée, alors i=-1 */

  for(i=u=0;u<n;u++){
    if(T[C[u]]<0) T[C[u]]=i++; /* la couleur C[u] n'a pas jamais été vue */ 
    C[u]=T[C[u]]; /* ici la couleur C[u] a déjà été vue */
  }

  free(T); /* plus besoin de T */
  G->int1=i-1; /* couleur max utilisée */
  return C;
}


int *power_law_seq(int n,double t,int *T){
/*
  Ecrit dans le tableau T une distribution de degré pour un graphe à
  n>0 sommets selon une lois en puissance d'exposant t>0. La
  distribution est codée par une suite (n_1,d_1,...,n_k,d_k,-1) de
  paires (n_i,d_i) qui signifie n_i sommets de degré d_i. Si T=NULL,
  alors T est alloué et renvoyé, sinon il être assez grand pour
  recevoir la distribution. Cette taille ne dépasse jamais 2n+1.

  La distribution est la suivante:
  - d_1=1, n_1 = ⌊exp(a)⌋ + n-s(a)
  - d_i=i, n_i = ⌊exp(a)/i^t⌋ pour i dans [2,p(a)]
  - où a est un réel minimisant |n-s(a)| avec
    s(a) := sum_{i=1}^p(a) ⌊exp(a)/i^t⌋
    p(a) := ⌊exp(a/t)⌋

  Attention ! Dans l'article original de [Lu01], il y a une erreur
  dans le signe du terme correctif r=n-s(a) pour les sommets de degré
  1. Il faut faire +r et non -r comme c'est écrit.

  Le nombre de paires (n_i,d_i) dans T est exactement p(a). Si n>0,
  alors T contient au moins 1 paire (et au plus n), car les d_i sont
  différents. On a s(0)=1 car p(0)=1. Aussi si a>=ln(n)+1, alors s(a)
  >= floor{exp(a)} >= n.

  En choisissant a0=0 et a1=ln(n)+1, on a alors s(a0) <= n <=
  s(a1). La fonction s(a) étant croissante, pour minimiser n-s(a) on
  réalise un recherche binaire pour a dans [a0,a1]. Le nombre
  d'itérations est au plus 64 (car sizeof(double)*8=64).
*/
  if((t<=0.0)||(n<=0)) return NULL;

  /* calcule la valeur de 'a' optimale */
  
  double a0=0.0,a1=log(n)+1.0; // intervalle pour a
  double a,b,e; // a=milieu de [a0,a1], b=meilleur a
  int p,i,j,s,cont=1; // cont=1 ssi on continue le calcul
  int r=0; // r=valeur minimum (=0 au départ)
  
  do{
    a=(a0<a1)? (a0+a1)/2.0 : a0; // si a0>=a1, on calcule s puis on s'arrête
    if((a==a0)||(a==a1)) cont=0; // intervalle trop faible, on calcule s puis on s'arrête
    e=exp(a); p=(int)exp(a/t); // p=nombre de paires, NB: p>=1
    for(s=0,i=1;i<=p;i++) s += (int)(e/pow(i,t)); // calcule s=s(a)
    if(s==n) cont=0; // valeur optimale, on va avoir b=a
    if((r==0)||(abs(n-s)<abs(r))) b=a,r=n-s; // NB: si s=n, alors b=a
    if(s<n) a0=a; // on est avant n
    if(s>n) a1=a; // on est après n
  }while(cont);

  /* ici on a calculé b, le meilleur a */

  p=(int)exp(b/t); // nombre de paires
  if(T==NULL) ALLOC(T,2*p+1);
  e=exp(b);

  /* écrit la distribution dans T */
  
  for(j=0,i=1;i<=p;i++){
    T[j++]=(int)(e/pow(i,t)); // NB: si i=1, T[0]=floor(exp(b))
    T[j++]=i;
  }
  T[j]=-1; // -1 terminal
  T[0]+=r; // correction pour les sommets de degré 1

  return T;
}


/***********************************

           GRAPHES DE BASE

***********************************/


#define ADJ_END  -1
#define ADJ_NAME -2
#define ADJ_INIT -3
#define ADJ_ERR  -4

/* 
   Les fonctions d'adjacences adj(i,j) doivent avoir les fonctionalités
   suivantes :

   1. adj(ADJ_INIT,j): initialise l'adjacence (appelée avant la génération)
   2. adj(ADJ_END,j) : finalise l'adjacence (appelée après la génération)
   3. adj(ADJ_NAME,j): construit l'étiquette du sommet j (via la variable NAME).
   4. adj(i,j) pour i,j>=0, retourne 1 ssi i est adjacent à j

   Rem1: adj(ADJ_INIT,j) doit toujours retourner le nombre de sommets
   du graphe et fixer la variable globale N à cette valeur. Certaines
   fonctions d'adjacence suppose que i<j. Si une valeur incorrecte des
   paramètres est détectée, adj(ADJ_INIT,j) devrait renvoyer ADJ_ERR.

   Rem2: C'est une mauvaise idée que d'utiliser des variables
   statiques dans les fonctions d'ajacence car elles peuvent être
   appelées entre elles par d'autres fonctions. Par exemple
   arboricity(i,j) est utilisée par plusieurs autres fonctions, sans
   forcément initialiser des variables statiques qui seraient
   initialisées par arboricity(ADJ_INIT,j).

   On devrait systématiquement, lors de l'initialisation
   adj(ADJ_INIT,j), avoir le test "if(N<=0) return N=0;" (le "=0" dans
   le "if(N<=0)" est important) permettant de générer le graphe vide
   sans rien faire d'autre surtout pour les graphes allouant de la
   mémoire.

   Les graphes utilisant le tableau REP[u] (représentation implicite)
   pour chaque sommet u devraient se terminer par la libération des
   tableaux, par exemple:

   if(i<0){
     if(i==ADJ_END) return kautz(i,j);
     ...;
   }

   Les graphes géométriques qui utilisent les tableaux XPOS et YPOS
   devraient se terminer par leur libération, par exemple:

   if(i<0){
     if(i==ADJ_END) return gabriel(i,j);
     ...;
   }

   A FAIRE: regrouper tous les "if(i<0)" et "if(j<0)" comme ceci:

   int adj(int i,int j){

     if(i<0){
       if(i==ADJ_NAME){
       ...;
       return i;
       }
       if(i==ADJ_END){
       ...;
       return i;
       }
       if(i==ADJ_INIT){
       ...;
       return N;
       }
       return ADJ_ERR;
     }

     ...;
   }
  
*/


int load(int i,int j)
/*
  Graphe défini par un fichier (ou l'entrée standard). A
  l'initialisation, il est chargé en mémoire dans la variable LOAD de
  type "graph". Suivant la valeur de LOAD->sort, le test d'adjacence
  est linéaire ou logarithmique en min{deg(i),deg(j)}.
*/
{
  if(j<0){
    if(j==ADJ_END){
      free_graph(LOAD);
      LOAD=NULL;
    }
    return 0;
  }

  if(i<0){
    LOAD=File2Graph(SPARAM,2); /* remplit LOAD à partir du fichier */
    if(LOAD->f>0){ /* si c'est une famille, on sélectionne le premier graphe */
      graph *G=ExtractSubgraph(LOAD->G[0],NULL,0,0); /* copie le premier graphe */
      free_graph(LOAD); /* libère complètement la famille LOAD */
      LOAD=G; /* LOAD=premier graphe */
    }
    if(!LOAD->sym) DIRECTED=1; /* si le graphe est asymétrique, affichage DIRECTED */
    return N=LOAD->n;
  }

  /* pour avoir du min{deg(i),deg(j)} en non-orienté */
  if((!DIRECTED)&&(LOAD->d[i]>LOAD->d[j])){ int k; SWAP(i,j,k); }

  return AdjGraph(LOAD,i,j);
}


int prime(int i,int j){
  if(j<0) return 0;
  if(i<0) return N=max(PARAM[0],0);
  return (i>1)? ((j%i)==0) : 0;
}


int paley(int i,int j){
/*
  Le résidu est r=|i-j|. Pour savoir si r est un carré, on teste s'il
  existe un entier k<=(n-1)/2 tel que (k*k)%n=r.
*/

  if(j<0) return 0;
  const int n=PARAM[0];
  if(i<0) return N=max(n,0);

  const int q=n/2;
  const int r=abs(i-j);
  int k;

  for(k=1;k<=q;k++)
    if(r==((k*k)%n)) return 1;
  return 0;
}


int mycielski(int i,int j){   // suppose i<j
  int ki,kj,b,k;

  if(j<0) return 0;
  if(i<0){
    k=PARAM[0];
    if(k<2) return N=0;
    return N=3*(1<<(k-2))-1;
  }

  ki=ceil(log2((double)(i+2)/3.));
  kj=ceil(log2((double)(j+2)/3.));
  k=3*(1<<kj)-2; /* rem: k est pair */
  b=(j==k);
  if(ki==kj) return b;
  if(b) return 0;
  j -= (k>>1);
  if(j==i) return 0;
  if(i<j) return mycielski(i,j);
  return mycielski(j,i);
}


int windmill(int i,int j){
  if(j<0) return 0;
  if(i<0) return N=(PARAM[0]<<1)+1;
  return (i==0)||((i&01)&&(j==i+1));
}


int matching(int i,int j){
  if(j<0) return 0;
  if(i<0) return N=max(PARAM[0]<<1,0);
  return (i==j-1)&&(j&1); /* utilise i<j */
}


int ring(int i,int j){
  if(j<0) return 0;

  const int n=PARAM[0];

  if(i<0) return N=max(n,0);

  const int t=PARAM[1]+2;
  int k=2;

  for(;k<t;k++)
    if((j==((i+PARAM[k])%n))||(i==((j+PARAM[k])%n))) return 1;

  return 0;
}


int cage(int i,int j){
  if(j<0) return 0;

  const int n=PARAM[0];
  const int k=PARAM[1];

  if(i<0){
    if(k<1) Erreur(6);
    return N=max(n,0);
  }

  if( (j==(i+1)%n)||(i==(j+1)%n) ) return 1;
  if( (j==(i+PARAM[(i%k)+2])%n)||(i==(j+PARAM[(j%k)+2])%n) ) return 1;

  return 0;
}


int crown(int i,int j){
  if(j<0) return 0;
  int k=PARAM[0];
  if(i<0) return N=(k<<1);
  return ((i<k)&&(j>=k)&&(i!=(j-k))); /* utilise i<j */
}


int fan(int i,int j){
  if(j<0) return 0;
  const int p=PARAM[0];
  const int q=PARAM[1];
  if(i<0) return N=p+q;
  return ((j==i+1)&&(j<p))||((i<p)&&(j>=p)); /* utilise i<j */
}


int chess(int i,int j){
  if(j<0){
    if(j==ADJ_NAME) NAME_Base(i,PARAM[0],2,",","()",1);
    return 0;
  }
  int p=PARAM[0];
  int q=PARAM[1];
  if(i<0) return N=p*q;
  int x=PARAM[2];
  int y=PARAM[3];
  int xi=i%p;
  int yi=i/p;
  int xj=j%p;
  int yj=j/p;  
  return ((abs(xi-xj)==x)&&(abs(yi-yj)==y))||((abs(xi-xj)==y)&&(abs(yi-yj)==x));
}


int grid(int i,int j){
  int x,y,k,z,p,h,b,d=PARAM[0];

  if(j<0){
    if(j==ADJ_NAME){
      int R[NAMEMAX],k;
      z=1; /* z=vrai ssi toutes les dimensions sont < 10 */
      for(k=0;k<d;k++){
	b=PARAM[k+1];
	R[k]=i%b;
	i /= b;
	z &= (b<11);
      }
      NAME_Vector(R,d,(z?"":","),(z?"":"()"),0,"%i");
    }
    return 0;
  }

  if(i<0){
    N=1;
    for(k=0;k<d;k++){
      p=PARAM[k+1];
      WRAP[k]=(p<0);
      p=abs(p);
      PARAM[k+1]=p;
      N *= p;
    }
    return N;
    }

  z=h=k=b=0;

  while((k<d)&&(b<2)&&(h<2)){
    p=PARAM[k+1];
    x=i%p;
    y=j%p;
    h=abs(x-y);
    if(h==0) z++;
    if((h>1)&&(WRAP[k])&&(((x==0)&&(y==p-1))||((y==0)&&(x==p-1)))) h=1;
    if(h==1) b++;
    i /= p;
    j /= p;
    k++;
  }

  return (b==1)&&(z==(d-1));
}


int rplg(int i,int j){
/*
  PARAM[0]=n
  DPARAM[0]=t
  DPARAM[1]=sum_i w_i
  DREP[i][0]=degré moyen du sommet i = (N/(i+1))^(1/(t-1))
*/

  if(j<0){
    if(j==ADJ_END) if(N>0) FREE2(DREP,N);
    return 0;
  }

  if(i<0){
    int k;
    double c,s,n;

    N=PARAM[0]; if(N<=0) return N=0;
    n=(double)N;
    ALLOCMAT(DREP,N,1);

    s=0.0;
    c=1.0/(DPARAM[0]-1.0);
    for(k=0;k<N;k++) s += (DREP[k][0]=pow(n/((double)k+1.0),c));
    DPARAM[1]=s;

    return N;
  }

  return (RAND01 < (DREP[i][0]*DREP[j][0]/DPARAM[1]));
}


int butterfly(int i,int j){
  int d=PARAM[0]+1; /* d=dim+1 */

  if(j<0) return 0;
  if(i<0){
    d--;          /* d=dim */
    N=1; N <<= d; /* N=2^dim */
    N *= d+1;     /* N=(dim+1)*2^dim */
    return N;
  }

  int x=i/d;i%=d; /* i -> (x,i) = (mot binaire,niveau) */
  int y=j/d;j%=d; /* j -> (y,j) = (mot binaire,niveau) */
  
  if(j==i+1) return (x==y) || ((x^y)==(1<<i));
  if(i==j+1) return (x==y) || ((x^y)==(1<<j));
  return 0;
}


int debruijn(int i,int j){
  int x,y,b=PARAM[1];

  if(j<0) return 0;

  if(i<0){
    int k;
    int d=PARAM[0];
    for(N=1,k=0;k<d;k++) N *= b;
    return N;
  }

  x=j-(i*b)%N;
  y=i-(j*b)%N;
  return ((0<=x)&&(x<b))||((0<=y)&&(y<b));
}


int barbell(int i,int j){

  if(j<0) return 0;

  const int n1=abs(PARAM[0]);
  const int n2=abs(PARAM[1]);
  const int p=PARAM[2];

  if(i<0) return N=n1+n2+p-1;

  /* utilise le fait que i<j */

  if(j<n1){
    if(PARAM[0]<0) return (j==i+1)||((i==0)&&(j==n1-1)); /* cycle 1 */
    return 1; /* clique 1 */
  }
  if(i>=n1-1+p){
    if(PARAM[1]<0) return (j==i+1)||((i==n1-1+p)&&(j==n1+n2+p-2)); /* cycle 2 */
    return 1; /* clique 2 */
  }
  if((n1-1<=i)&&(j<n1+p)) return (j-i==1); /* chemin */
  return 0;
}


int shuffle(int i,int j){
  int n=PARAM[0];
  if(j<0){
    if(j==ADJ_NAME) NAME_Base(i,2,n,"","",-1);
    return 0;
  }
  n=(1<<n);
  if(i<0) return N=n;
  
  if((i>>1)==(j>>1)) return 1;
  n>>=1; // n=N/2
  if(j==((i<<1)%N+(i>=n))) return 1;
  if(j==((i>>1)+((i&01)?n:0))) return 1;
  return 0;
}


int kautz(int i,int j){
  /*
    A chaque sommet i qui est entier de [0,b*(b-1)^(d-1)[, on associe
    une représentation (s_1,...,s_d) codée sous la forme d'un entier
    de [0,b^d[. Alors i adjacent à j ssi i et j sont adjacents dans le
    De Bruijn.
   */

  if(j<0){
    if(j==ADJ_END) if(N>0) FREE2(REP,N);
    return 0;
  }

  int d,b,k,l,q,r,s,t,x;
  if(i<0){
    d=PARAM[0];
    N=x=b=PARAM[1];
    t=b-1;
    for(k=1;k<d;k++) {
      N *= t;
      x *= b;
    }
    PARAM[2]=x; /* nb de sommets de De Bruijn */
    PARAM[3]=N; /* nb sommets de Kautz */
    if(N<=0) return N=0;
    ALLOCMAT(REP,N,1);

    for(l=0;l<N;l++){ /* pour tous les sommets faire .. */
      /* On voit un sommet l comme (r,s2...sd) avec r dans [0,b[ et
	 s_i de [0,b-1[. On le converti en (x1,...,xd) avec xd=r et
	 x_(d-1)=s2 ou s2+1 suivant si s2<r ou pas, etc. */
      r=x=l%b;
      q=l/b;
      for(k=1;k<d;k++){
	s=q%t;
	s+=(s>=r);
	r=s;
	q=q/t;
	x=x*b+r;
      }
      REP[l][0]=x;
    }
    return N;
  }

  N=PARAM[2]; /* modifie le nb de sommets */
  x=debruijn(REP[i][0],REP[j][0]);
  N=PARAM[3]; /* rétablit le nb de sommets */
  return x;
}


int ggosset(int i,int j){
  /*
    ggosset p k d_1 v_1 ... d_k v_k
    mais PARAM = p k d d_1 v_1 ... d_k v_k

    Sommet: permutations du vecteur (v_1...v_1, ..., v_k...v_k) (ou de
    son opposé) telles que le nombre de valeurs entières v_t est
    d_t. On pose REP[i] = vecteur du sommet i qui est de taille
    d=d_1+...+d_k. On a l'arête i-j ssi le produit scalaire entre
    REP[i] et REP[j] vaut p.

    Pour calculer tous les vecteurs de taille d (les sommets) à partir
    de (v_1...v_1,...,v_k...v_k) on procède comme suit (codage d'un
    multi-ensemble à k valeurs):

    On choisit les positions de v_1 (il y en a Binom(d,d_1) possibles,
    puis on choisit les positions de v_2 parmi les d-d_1 restantes (il
    y en a Binom(d-d_1,d_2)), puis les positions de v_3 parmi les
    d-d_1-d_2 retantes, etc. Pour chaque t, les positions des d_t
    valeurs v_t sont codées par un sous-ensemble S[t] de [0,d-1] de
    taille d_t.
   */

  int t,p,d=PARAM[2];

  if(j<0){
    if(j==ADJ_END) return kautz(i,j);
    if(j==ADJ_NAME) NAME_Vector(REP[i],d,",","[]",1,"%i");
    return 0;
  }

  if(i<0){
    int m,u,v,c,**S,k=PARAM[1];

    /* Calcule N */
    m=d; N=2; if(d<=1) Erreur(6); /* paramètre incorrect */
    for(t=3;m>0;t += 2){
      p = PARAM[t]; /* p=d_i */
      if(p<0) Erreur(6); /* paramètre incorrect */
      N *= Binom(m,p);
      m -= p;
    }
    if(N<=0) return N=0;
    ALLOCMAT(REP,N,d); /* vecteur de taille d représentant les sommets */
    NALLOC(int,P,d); /* tableau intermédiaire */
    ALLOCMAT(S,k,d); /* on réserve k tableaux (sous-ensembles) de taille <= d */
    for(t=0;t<k;t++) NextSet(S[t],-1,d); /* initialise les sous-ensembles */
    /* Note: taille |S[t]|=PARAM[(t<<1)+3] et v_t=PARAM[(t<<1)+4] */

    for(u=0;u<N;u+=2){
      /* Pour chaque sommet u on fait:

	 1. on remplit REP[u] et REP[u+1] à partir des sous-ensembles S[0]...S[k-1]
	 2. on incrémente les sous-ensembles S[0]...S[k-1]
	 
	 Pour l'étape 1, il faut passer par un tableau intermédiaire P
	 puis remplir REP[u] et REP[u+1] à partir de P. Supposons d=5,
	 k=3, et S[0]={1,3} S[1]={1}, S[2]={0,1}.  On met dans P les
	 indices t des S[t] aux bons emplacements comme suit:

           - au départ P est vide: P={-,-,-,-,-}
	   - puis on ajoute S[0]:  P={-,0,-,0,-}
	   - puis on ajoute S[1]:  P={-,0,1,0,-}
	   - puis on ajoute S[2]:  P={2,0,1,0,2}
      */

      /* Calcule P */
      for(t=0;t<d;t++) P[t]=-1; /* tableau vide au départ */
      for(t=0;t<k;t++){ /* pour chaque sous-ensemble S[t] */
	m=-1;
	for(p=v=0;p<PARAM[(t<<1)+3];p++){ /* on parcoure S[t] */
	  /* mettre t dans la S[t][p]-ème case vide de P à partir de l'indice v */
	  /* utilise le fait que dans S[t] les nombres sont croissant */
	  /* v=position courante de P où l'on va essayer d'écrire t */
	  /* c=combien de cases vides de P à partir de v faut-il encore sauter ? */
	  /* si P[v]=-1 et c=0 alors on écrit P[v]=t */
	  c=S[t][p]-m;
	  m=S[t][p]; /* mémorise l'ancienne valeur de S[t][p] */
	  while(c>0){
	    if(P[v]<0){
	      c--;
	      if(c==0) P[v]=t; /* écrit t et passe à la casse suivante */
	    }
	    v++;
	  }
	}
      }

      /* Remplit REP[u] et REP[u+1] grâce au tableau P (et incrémenter u) */
      for(t=0;t<d;t++){
	v=PARAM[(P[t]<<1)+4]; /* valeur v_t à écrire */
	REP[u][t]=v;
	REP[u+1][t]=-v;
      }

      /* Incrémente S[0]...S[k-1] grâce à NextSet() */
      t=0; /* on commence avec S[0] */
      m=d; /* S[t] dans {0...d-2} */
      v=3; /* PARAM[v] = taille de S[t] */
      while((NextSet(S[t],m,PARAM[v]))&(t<k)){
	t++; /* si S[t] fini on passe à S[t+1] */
	m -= PARAM[v];
	v += 2;
      }
      /* si t==k alors on a atteint le dernier sommet */
    }

    free(P);
    FREE2(S,k);
    return N;
  }

  /* Calcule le produit scalaire de REP[i] et REP[j] */
  for(p=t=0;t<d;t++) p += REP[i][t]*REP[j][t]; 
  return (p==PARAM[0]);
}


int gpstar(int i,int j)
/*
  REP[i][0...n-1] = représentation de la permutation du sommet i.
*/
{
  int k,c;
  const int n=PARAM[0];

  if(j<0){
    if(j==ADJ_END) return kautz(i,j);
    if(j==ADJ_NAME){
      if(n<10) NAME_Vector(REP[i],n,"","",1,"%i");
      else NAME_Vector(REP[i],n,",","()",1,"%i");
    }
    return 0;
  }
  if(i<0){
    for(N=1,k=2;k<=n;k++) N *= k;
    ALLOCMAT(REP,N,n); /* ici N>0 */
    NALLOCZ(int,P,n,_i); /* initialise P */

    /* génère toutes les permutations */
    for(c=0;c<N;c++){
      for(k=0;k<n;k++) REP[c][k]=P[k]+1; /* copie P dans REP */
      NextPermutation(P,n,NULL);
    }

    free(P);
    return N;
  }

  /* Distance de Hamming: on compte les différences entre les tableaux
     REP[i] et REP[j] */
  for(k=c=0;k<n;k++) c += (REP[i][k]!=REP[j][k]);
  return (c==PARAM[1]);
}


int linial(int i,int j){
/*
  REP[i][0...t-1] = représentation le nom du sommet i (sa vue).
*/

  const int t=PARAM[1];

  if(j<0){
    if(j==ADJ_END) return kautz(i,j);
    if(j==ADJ_NAME){
      if(PARAM[0]<10) NAME_Vector(REP[i],t,"","",1,"%i");
      else NAME_Vector(REP[i],t,",","()",1,"%i");
    }
    return 0;
  }

  int k,u;

  if(i<0){
    const int m=PARAM[0];
    if((m<t)||(t<1)) return N=0; /* graphe vide si n<t */
    for(N=1,k=m-t+1;k<=m;k++) N *= k; /* calcul N */
    ALLOCMAT(REP,N,t);
    NALLOC(int,S,t);
    NALLOC(int,P,t);
    NextArrangement(S,P,-1,t); /* initialisation de S et P */
    for(u=0;u<N;u++){ /* génère tous les arrangements */
      for(k=0;k<t;k++) REP[u][k]=S[P[k]];
      NextArrangement(S,P,m,t);
    }
    free(P);
    free(S);
    return N;
  }

  if((REP[i][0]!=REP[j][t-1])||(PARAM[0]==t)){
    for(u=1,k=0;u<t;u++,k++) if(REP[i][u]!=REP[j][k]) break;
    if(u==t) return 1;
  }
  if(i<j) return linial(j,i);
  return 0;
}


int linialc(int i,int j){
  if(j<0) return linial(i,j);
  if(i<0){
    int m,t,u,k,v,m1,x,y;
    m=PARAM[0];
    t=PARAM[1];
    N=m; m1=m-1;
    for(N=m,u=m1,k=1;k<t;k++) N *= u; /* calcule N=m*(m-1)^t */
    ALLOCMAT(REP,N,t);
    /* on transforme u en (x0,x1,...x_t) avec x0 in [0,m[ et x_i in [0,m-1[ */
    for(u=0;u<N;u++){
      x=REP[u][0]=(u%m);
      for(v=u/m,k=1;k<t;k++){
	y = v%m1; /* y in [0,m-1[ */
	v /= m1;
	x=REP[u][k]=y+(y>=x); /* si x=y, on incrémente y */
      }
    }
    return N;
  }

  return linial(i,j);
}


int pancake(int i,int j){
/*
  REP[i][0...n-1] = représentation de la permutation du sommet i.
*/
  int k,l;

  if((j<0)||(i<0)) return gpstar(i,j);

  /* Test d'adjacence à partir de REP[i] et REP[j] */
  k=PARAM[l=0];
  do k--; while(REP[i][k]==REP[j][k]); /* i<>j, donc on s'arrête toujours */
  while((k>=0)&&(REP[i][k]==REP[j][l])) k--,l++; /* teste le "reversal" */
  return (k<0); /* adjacent si préfixe=reversal */
}


int bpancake(int i,int j){
/*
  REP[i][0...n-1] = représentation de la permutation signée du sommet
  i. Il s'agit d'une valeur de {1,2,...,n,-1,-2,...,-n} (on évite
  soigneusement 0, car -0=+0.
*/
  int k,l;

  if(j<0){
    if(j==ADJ_END) return kautz(i,j);
    if(j==ADJ_NAME) NAME_Vector(REP[i],PARAM[0],"","",1,"%+i");
    return 0;
  }
  if(i<0){
    int p,q,t,c,u,n=PARAM[0];
    for(t=1,k=p=2;k<=n;k++){ t *= k; p <<= 1; } /* calcule t=n! et p=2^n */
    N=p*t; /* N=nb de sommets, forcément >0 */
    ALLOCMAT(REP,N,n); /* permutations signées représentant les sommets */
    NALLOCZ(int,P,n,_i); /* initialise une permutation P (non signée) */

    /* Génère toutes les permutations signées. On représente les
       signes par les bits de l'entier q=0...2^n-1. Si bit à 1 -> +1,
       et si bit à 0 -> -1 */
    for(c=u=0;c<t;c++){ /* répète n! fois */
      for(q=0;q<p;q++,u++) /* répète 2^n fois */
	for(k=0,l=1;k<n;k++,l<<=1){ /* l=mask=bit-0,bit-1,bit-2...,bit-(n-1) */
	  REP[u][k]=P[k]+1; /* copie P dans REP avec le signe +/-1 */
	  if(q&l) REP[u][k]=-REP[u][k]; /* copie P dans REP avec le signe +/-1 */
	}
      NextPermutation(P,n,NULL);
    }

    free(P);
    return N;
  }

  /* Test d'adjacence à partir de REP[i] et REP[j] */
  k=PARAM[l=0];
  do k--; while(REP[i][k]==REP[j][k]); /* i<>j, donc on s'arrête toujours */
  while((k>=0)&&(REP[i][k]==-REP[j][l])) k--,l++; /* teste le "reversal" */
  return (k<0); /* adjacent si préfixe=-reversal */
}


int pstar(int i,int j){
/*
  REP[i][0...n-1] = représentation de la permutation du sommet i.
*/

  if((j<0)||(i<0)) return gpstar(i,j); /* c'est comme gpstar() */

  /* Il faut deux différences dont le premier chiffre */
  return (REP[i][0]!=REP[j][0]) && (gpstar(i,j));
}


int gabriel(int i,int j){
/*
  L'initialisation et la terminaison sont communes à beaucoup de
  graphe géométriques. L'adjacence est en O(N). Il est importnat de
  tester aussi les sommets supprimés (par -delv), sinon le résultat
  n'est pas un sous-graphe.
*/
  if(j<0){
    free(XPOS),free(YPOS);
    free(XSEED),free(YSEED);
    XPOS=YPOS=XSEED=YSEED=NULL;
    return 0;
  }
  
  if(i<0){
    N=PARAM[0]; if(N<=0) return N=0;
    InitXY();
    return N;
  }

  const int n=PARAM[0];
  double xi,yi;

  xi=XPOS[i],yi=YPOS[i]; /* on sauvegarde les coordonnées de i */
  XPOS[i]=(xi+XPOS[j])/2.0; /* on change les coordonnées de i */
  YPOS[i]=(yi+YPOS[j])/2.0; /* i <- milieu du segment i-j */
  const double r=Norme(i,j); /* rayon du disque centré au milieu de i-j */

  int z;
  for(z=0;z<n;z++) /* teste même les sommets supprimés */
    if((Norme(i,z)<r)&&(z!=i)) z=n; /* ici i=milieu i-j */
  /* si on sort par le if, alors z=n+1, sinon z=n */

  XPOS[i]=xi,YPOS[i]=yi; /* on remet les coordonnées de i */
  return (z==n);
}


int pat(int i,int j){
/*
  Graphe issu du jeu de Pat Morin. Un sommet correspond à un sommet
  (x,y) d'une des k grilles. On suppose que i<j, ce qui revient à dire
  que la grille de i est placée avant ou est égale à celle de j.

  Exemple: p=q=3 et k=4

  06 07 08  15 16 17  24 25 26  33 34 35
  03 04 05  12 13 14  21 22 23  30 31 32
  00 01 02  09 10 11  18 19 20  27 28 29

  Grille 0  Grille 1  Grille 2  Grille 3

*/
  if(j<0) return gabriel(i,j);

  const int p=PARAM[0]; // p=colonne
  const int q=PARAM[1]; // q=ligne
  const int r=PARAM[2]; // r=round
  const int pq=p*q;

  if(i<0){
    N=pq*r;
    if(N<=0) return N=0;
    if((p<=0)||(q<=0)||(r<=0)) Erreur(6);

    /* Détermine les coordonnées des points un controler le
       dessin. Les grilles sont placées à des hauteurs de plus en plus
       grandes pour "voir" les arêtes. Au départ (z=0), la première
       grille est mise à une hauteur 0, puis à une hauteur 1, puis à
       une hauteur 3, puis à une hauteur 6, etc. La hauteur de la
       grille pour z quelconque est z(z+1)/2. En fait, on utilise
       z(z+1)/2.5 pour une meilleure lisibilité. */

    ALLOC(XPOS,N);
    ALLOC(YPOS,N);
    int u,x,y,z;
    double h;

    for(u=0;u<N;u++){// pour tous les sommets, faire:
      z=u/(pq);      // z=grille
      x=u%p;         // x=colonne
      y=(u%(pq))/p;  // y=ligne
      h=z*(z+1)/2.5; // h=décalage vers le haut de la grille z
      XPOS[u]=(double)(x+z*p)/(double)max(p,q);
      YPOS[u]=(double)(y+h*q)/(double)max(p,q);
    }

    XYtype=XY_USER; /* coordonnées fixées par l'utilisateur */
    InitXY(); // pour les options -xy noise/scale ... */
    return N;
  }

  // z=grille
  const int zi = i/pq;
  const int zj = j/pq;

  /* ici: zi<=zj car i<j */

  // x=colonne
  const int xi = i%p;
  const int xj = j%p;

  // y=ligne
  const int yi = (i%pq)/p;
  const int yj = (j%pq)/p;

  if(zj==zi)
    return ((xj>=xi)&&(yj<=yi)) || ((xj<=xi)&&(yj>=yi));
  
  /* ici: zj>zi */
  return ((xj>=xi)&&(yj==yi)) || ((xj==xi)&&(yj>=yi));
}


int thetagone(int i,int j){
/*
  L'adjacence est en O(k*N).
*/
  if(j<0) return gabriel(i,j);

  const int p=PARAM[1];
  const int k=PARAM[2];

  if(i<0){
    if(p<3) PARAM[1]=-1; /* p inifini */
    if(k<1) return N=0;
    return gabriel(i,j);
  }

  /*
    Adjacence: pour tous les axes t, on cacule P_t(i,j)=distgone(),
    puis on détermine s'il existe un autre sommet z avec P_t(i,z) plus
    petit. Si c'est non (et que P_t(i,j) est finie), alors i est
    adjacent à j, sinon ils ne le sont pas (mais j peut être adjacent
    à i !).
   */

  int t,z;
  double d;

  const int n=PARAM[0];
  const double w=DPARAM[0];

  for(t=0;t<k;t++){ /* pour tous les axes t */
    d=distgone(i,j,t,p,k,w); /* calcule P_t(i,j) */
    if(d<DBL_MAX){ /* distance infinie */
      for(z=0;z<n;z++) /* pour tous les autres sommets z, même supprimés ! */
	if((z!=i)&&(distgone(i,z,t,p,k,w)<d)) z=n; /* z plus proche ? */
      /* si oui, on arrête: P_t(i,j) contient z */
      if(z==n) return 1; /* on a pas trouvé de sommet z plus proche que j */
    }
  }

  /*
    A priori ici il n'y a pas d'arête entre i et j. Il faut cependant
    tester aussi adj(j,i) car la distance P_t(i,j) n'est pas
    symétrique.
  */
  if(i>j) return 0;
  return thetagone(j,i);
}


int udg(int i,int j){
/*
  NB: deux points peuvent avoir les même coordonnées.
*/
  if(j<0) return gabriel(i,j);

  if(i<0){
    if(NORM==2) DPARAM[0] *= DPARAM[0]; /* si norme L_2, alors r=r^2 */
    return gabriel(i,j);
  }

  return (Norme(i,j)<=DPARAM[0]);
}


int rng(int i,int j){
/*
  Adjacence en O(N).
*/
  if((j<0)||(i<0)) return gabriel(i,j);
  
  int z;
  const int n=PARAM[0];
  const double r=Norme(i,j);

  for(z=0;z<n;z++) /* teste même les sommets supprimés */
      if(fmax(Norme(i,z),Norme(j,z))<r) return 0;

  return 1; /* fonction de distance symétrique */
}


int nng(int i,int j){
/*
  Adjacence en O(N).
*/
  if((j<0)||(i<0)) return gabriel(i,j);

  int z;
  const int n=PARAM[0];
  const double r=Norme(i,j);

  for(z=0;z<n;z++) /* teste même les sommets supprimés */
      if((Norme(i,z)<r)&&(z!=i)) z=n; /* alors d(z,i)<d(i,j) */

  if(z==n) return 1;

  /* Avant de dire que adj(i,j)=0, il faut tester adj(j,i), car le
     test n'est pas symétrique. */

  if(i>j) return 0;
  return nng(j,i);
}


int hexagon(int i,int j){
/*
  Utilise le fait que i<j.
    
  On voit le graphe comme une grille de p+1 lignes de 2q+2 colonnes
  allant du coin en bas à droite (0,0) au coin (p,2q+1) (en notation
  (ligne,colonne)), et dans laquelle deux sommets ont été supprimés: le
  coin (0,2q+1) et le coin (p,2q+1) si p est impair, le coin (p,0) si
  p est pair. Les numéros sont consécutifs sur une ligne, de haut en
  bas.

 Ex:

 hexagon 3 2   hexagon 2 3

 o-o-o-o-o x   o-o-o-o-o-o-o x
 |   |   |     |   |   |   |
 o-o-o-o-o-o   o-o-o-o-o-o-o-o
   |   |   |     |   |   |   |
 o-o-o-o-o-o   x o-o-o-o-o-o-o
 |   |   |
 o-o-o-o-o x

*/
  int li,lj,ci,cj;

  if(j<0) return 0;

  const int p=PARAM[0];
  const int q=PARAM[1];
  const int t=(q<<1)+2; /* longueur d'une ligne */

  if(i<0) return N=(p+1)*t-2;

  if(i>=t-1) i++; /* on insère virtuellement le coin (0,2q+1) */
  if(j>=t-1) j++;
  if((p&01)==0){ /* si p est impair, on a rien à faire */
    if(i>=p*t) i++; /* on insère virtuellement le coin (p,0) si p est pair */
    if(j>=p*t) j++;
  }

  /* on calcule les coordonnées de i et j placés sur cette
     grille (avec les coins manquant) */

  li=i/t;ci=i%t;
  lj=j/t;cj=j%t;
  
  /* utilise le fait que i<j: dans le dernier cas lj=li+1 */
  return ((li==lj)&&(abs(ci-cj)==1)) || ((ci==cj)&&(lj==(li+1))&&((ci&01)==(li&01)));
}


int whexagon(int i,int j){
/* Utilise le fait que i<j et hexagon(i,j) */
  int li,ci,lj,cj;

  if(j<0) return 0;

  const int p=PARAM[0];
  const int q=PARAM[1];

  if(i<0) return N=hexagon(i,j)+p*q;

  int t=N-p*q; /* nb de sommets de l'hexagone */

  /* teste si i et j sont dans l'hexagone */
  if((i<t)&&(j<t)) return hexagon(i,j);

  /* teste si i et j sont hors l'hexagone */
  if((i>=t)&&(j>=t)) return 0;

  /* on a i dans l'hexagone et j en dehors car i<j */
  lj=(j-t)/q;cj=(j-t)%q; /* j est le centre de l'hexagone (lj,cj) */
  t=(q<<1)+2; /* longueur d'une ligne de l'hexagone */

  /* on calcule les coordonnées de i dans l'hexagone */
  if(i>=t-1) i++; /* on corrige */
  if(((p&01)==0)&&(i>=p*t)) i++; /* on recorrige */
  li=i/t;ci=i%t; /* (li,ci) coordonnées de i */

  return ( ((li==lj)||(li==(lj+1))) && (abs(2*cj+1+(lj&01)-ci)<2) );
}


int hanoi(int i,int j){
/*
  Adjacence: on écrit i (et j) en base b, mot de n lettres. i et j
  sont adjacents ssi i=Puv...v et j=Pvu...u où P est un préfixe
  commun, et u,v des lettres qui se suivent (modulo b).
*/
  int ri,rj,u,v,k;

  const int n=PARAM[0];
  const int b=PARAM[1];

  if(j<0){
    if(j==ADJ_NAME){
      if(b<10) NAME_Base(i,b,n,"","",1);
      else NAME_Base(i,b,n,",","",1);
    }
    return 0;
  }

  if(i<0){
    if(b<2) return N=0;
    for(N=k=1;k<=n;k++) N *= b; /* calcule le nombre de sommets N=b^n */
    return N;
  }

  /* on égraine les chiffres en base b */
  
  for(ri=rj=k=0;k<n;k++){
    if(ri!=rj) break; /* on s'arrête dès qu'on diffère, on garde k */
    ri=i%b; i/=b; /* ri=dernier chiffre, i=i sans dernier chiffre */
    rj=j%b; j/=b; /* rj=dernier chiffre, j=j sans dernier chiffre */
  }
  if((((ri+1)%b)!=rj)&&(((rj+1)%b)!=ri)) return 0; /* alors pas voisin */

  u=ri; v=rj; /* ici u et v sont consécutifs (mod b) */
  for(;k<n;k++){
    ri=i%b; i/=b;
    rj=j%b; j/=b;
    if(ri!=v) return 0; /* pas bon */
    if(rj!=u) return 0; /* pas bon */
  }

  return 1;
}


int sierpinski(int i,int j){
/*
  On utilise REP[i][k] pour représenter le sommet i. C'est un mot d'au
  plus n lettres de [0,b[. On pose REP[i][n]=L où L est la longueur du
  mot.

  Le mot du sommet i représente la suite des cycles auxquels il
  appartient, sauf s'il est l'un des b sommets initiaux. Prennons
  b=3. Les 3 sommets du triangle sont 0,1,2. Si n>1, les autres
  sommets commenceront tous par 0. On a alors 3 autres triangles,
  numérotés 0,1,2, le triangle i à pour sommet le sommet i. Les 3
  sommets internes sont 00,01,02. Le sommet 00 partage les triangles 0
  et 1, 01 les triangles 1 et 2, 02 les triangles 2 et 0. Si n>2, tous
  les autres sommets commenceront par 00, 01 ou 02 suivant le
  sous-triangle auquel ils appartiennent. Dans le sous-triangle 0, les
  3 sommets internes seront 000, 001, 002. Etc.

  Ex: n=2 b=3         0 
                     /  \
                   00 -- 02
                  /  \  /  \
                 1 -- 01 -- 2

  Adjacence:

  CAS 1: extrémité (|i|=1 et |j|=n)
  Si i=x et n>1, alors il est voisin de:
  - 0 x^{n-2} x
  - 0 x^{n-2} (x-1)
  Si i=x et n=1, alors c'est le CAS 2.
  
  Soit P=le plus préfixe commun entre i et j, et k=|P|

  CAS 2: sommet du triangle le plus interne (|i|=|j|=n)
  Si i=Px, k=n-1, et n>0, alors il est voisin de:
  - P (x+1)
  - P (x-1)

  CAS 3: sommet entre deux triangles (1<|i|<n et |j|=n)
  Si i=Px, alors:

  CAS 3.1: i=P=Qx est préfixe de j (k=|i|).
  Alors il est voisin de:
  - P (x+1)^{n-k-1} x
  - P (x+1)^{n-k-1} (x+1)

  CAS 3.2: i=Px.
  Alors il est voisin de:
  - P (x+1) x^{n-p-2} x
  - P (x+1) x^{n-p-2} (x-1)

*/
  int k,t,r,x,c,li,lj;

  const int n=PARAM[0];
  const int b=PARAM[1];

  if(j<0){
    if(j==ADJ_END) return kautz(i,j);
    if(j==ADJ_NAME){
      if(b<10) NAME_Vector(REP[i],REP[i][n],"","",1,"%i");
      else NAME_Vector(REP[i],REP[i][n],",","()",1,"%i");
    }
    return 0;
  }

  if(i<0){
    if((b<=2)||(n==0)) return N=0; /* graphe vide, graphe non défini */
    for(N=b,k=2;k<=n;k++) N=b*N-b;
    ALLOCMAT(REP,N,n+1); /* REP[t][k]=k-ème lettre du sommet t */

    for(t=0;t<N;t++){ /* calcule les noms */
      x=t;
      k=r=0;
      if(x<b) REP[t][0]=x; /* un des b sommets du 1er cycle? */
      else{
	x -= b;	
	while(1){
	  REP[t][k++]=r;
	  if(x<b) { REP[t][k]=x; break; }
	  x -= b; r=x%b; x /= b;
	}
      }
      REP[t][n]=k+1-(n==0); /* longueur du mot, corrige si n=0 */
    }

    return N;
  }

  li=REP[i][n]; /* longueur de i */
  lj=REP[j][n]; /* longueur de j */
  /* Propriété: si i<j, alors li<=lj */
  if(lj<n) return 0;

  /* CAS 1 */
  if((li==1)&&(n>1)&&(REP[j][0]==0)){
    x=REP[i][0];
    for(c=t=1;t<=n-2;t++) if(REP[j][t]!=x) { c=0; break; }
    if(c){ /* c=vrai ssi j=0x^(n-2) */
      if(REP[j][n-1]==x) return 1;
      if(REP[j][n-1]==((x-1+b)%b)) return 1;
    }
  }

  /* calcule k=longueur du préfixe commun */
  for(k=0;(k<li)&&(k<lj)&&(k<n);k++)
    if(REP[j][k]!=REP[i][k]) break;

  /* CAS 2 */
  if((li==n)&&(k==n-1)){
    x=REP[i][k];
    if(REP[j][k]==((x+1)%b)) return 1;
    if(REP[j][k]==((x-1+b)%b)) return 1;
  }

  /* CAS 3 */
  if((li==1)||(li==n)) return 0;
  x=REP[i][li-1];
  /* ici on a 1<|i|<n, |j|=n, et x=dernière lettre de i */

  /* CAS 3.1 */
  if(k==li){
    for(c=1,t=k;t<=n-2;t++) if(REP[j][t]!=((x+1)%b)) { c=0; break; }
    if(c){
      if(REP[j][n-1]==x) return 1;
      if(REP[j][n-1]==((x+1)%b)) return 1;
    }
  }
  
  /* CAS 3.2 */
  if((k==li-1)&&(REP[j][k]==((x+1)%b))){
    for(c=1,t=k+1;t<=n-2;t++) if(REP[j][t]!=x) { c=0; break; }
    if(c){
      if(REP[j][n-1]==x) return 1;
      if(REP[j][n-1]==((x-1+b)%b)) return 1;
    }
  }

  return 0; /* si i>j, alors on ne fait rien */
}


int rpartite(int i,int j)
/*
  On se sert du fait que i<j pour le test d'adjacence. Les sommets
  sont numérotés consécutivement dans chacune des parts. Utilise le
  tableau WRAP en interne: WRAP[k]=a_1+...+a_k est la somme partielle,
  WRAP[0]=0. Donc les sommets de la part i sont numérotés de WRAP[i-1]
  à WRAP[i] (exclu).
*/
{
  int k,r,s;

  if(j<0) return 0;

  r=PARAM[0];

  if(i<0){
    N=WRAP[0]=0;
    for(k=1;k<=r;k++){
      N += PARAM[k];
      WRAP[k]=N;
    }
    return N;
  }

  /*
    Pour calculer adj(i,j), avec i<j, on calcule les numéros de la
    part de i et de j: adj(i,j)=1 ssi part(i)<>part(j). Pour celà, on
    fait une recherche dichotomique de la part de i (c'est-à-dire d'un
    k tq WRAP[k]<=i<WRAP[k+1]). La complexité est O(logr), r=nb de
    parts.
  */

  /* Cherche la part de i dans [s,r[: au départ c'est dans [0,r[ */
  s=0;
  k=r/2;
  while(s<r){
    if(i<WRAP[k]){
      if(j>=WRAP[k]) return 1; /* i et j sont dans des parts différentes */
      r=k; /* ici i<j<WRAP[k]: on cherche dans [s,k[ */
      k=(s+k)>>1;
    }else{ /* ici WRAP[k]<=i<j */
      if(j<WRAP[k+1]) return 0; /* i et j sont dans la part k */
      s=k; /* ici WRAP[k]<=i<j: on cherche dans [k,r] */
      k=(k+r)>>1; 
    }
  }
  return (j>=WRAP[k+1]); /* ici i est dans la part k. On vérifie si j aussi */
}


int aqua(int i,int j)
/*
  On se sert de REP[i][0..n], n=PARAM[0].
*/
{
  int n=PARAM[0]; /* n=nb de paramètres */
  if(j<0){
    if(j==ADJ_END) return kautz(i,j);
    if(j==ADJ_NAME) NAME_Vector(REP[i],n,",","",1,"%i");
    return 0;
  }

  int k,x,y;
  int *C=PARAM+1; /* C=tableau de contraintes */

  if(i<0){
    x=*C; /* s=PARAM[1]=premier terme=la somme */
    if(n<0) n=0;
    int *S=NextPart(NULL,n,x,C);
    N=0; /* calcule d'abord N pour ALLOCREP() */
    do N++; while(NextPart(S,n,x,C)!=NULL);
    ALLOCMAT(REP,N,n); /* ici N>0 */
    N=0; /* calcule les sommets */
    do{
      for(k=0;k<n;k++) REP[N][k]=S[k];
      N++;
    }while(NextPart(S,n,x,C)!=NULL);

    return N;
  }

  if(i==j) return 0;

  /* compte et mémorise les différences entre REP[i] et REP[j]: il y en a au moins 2 */
  for(k=0,x=y=-1;k<n;k++)
    if(REP[i][k]!=REP[j][k]){
      if((x>=0)&&(y>=0)) return 0; /* si plus de 2 différences, alors pas d'arc */
      if(x<0) x=k; else if(y<0) y=k;
    }

  /* soit on a versé x vers y */
  /* k=quantité que l'on peut verser de x à y */
  k=min(C[y],REP[i][x]+REP[i][y])-REP[i][y];
  if((REP[j][y]==REP[i][y]+k)&&(REP[j][x]==REP[i][x]-k)) return 1;

  /* soit on a versé y vers x */
  /* k=quantité que l'on peut verser de y à x */
  k=min(C[x],REP[i][y]+REP[i][x])-REP[i][x];
  if((REP[j][x]==REP[i][x]+k)&&(REP[j][y]==REP[i][y]-k)) return 1;

  return 0;
}


int icosahedron(int i,int j){
  int T[]={
    5,0,1,2,3,4,5,6,7,8,
    9,10,11,6,4,7,3,8,2,
    9,1,10,5,1,2,0,3,4,0,
    5,6,10,9,11,9,11,8,11,7,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=12;
  return GraphFromArray(i,j,T);
}


int rdodecahedron(int i,int j){
  int T[]={
    0,1,2,3,4,5,0,6,
    7,1,7,8,2,8,9,3,9,
    10,4,10,11,5,11,6,
    7,12,9,12,11,GFA_CUT,
    0,13,2,13,4,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=14;
  return GraphFromArray(i,j,T);
}


int cuboctahedron(int i,int j){
  int T[]={
    0,1,2,3,4,5,0,6,7,2,
    8,9,5,6,10,7,1,11,0,GFA_CUT,
    3,10,4,9,11,8,3,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=12;
  return GraphFromArray(i,j,T);
}


int tutte(int i,int j){
  int T[]={
    0,GFA_PATH,8,
    4,8,9,3,9,10,11,2,
    11,12,1,12,13,14,10,
    14,7,6,15,13,GFA_CUT,
    0,16,GFA_PATH,7,
    23,19,23,24,18,24,25,26,17,
    26,27,16,27,28,29,25,
    29,22,21,30,28,GFA_CUT,
    0,31,GFA_PATH,7,
    38,34,38,39,33,39,40,41,32,
    41,42,31,42,43,44,40,
    44,37,36,45,43,GFA_CUT,
    15,20,GFA_CUT,
    30,35,GFA_CUT,45,5,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=46;
  return GraphFromArray(i,j,T);
}


int herschel(int i,int j){
  int T[]={
    0,GFA_PATH,10,
    10,3,2,9,0,7,6,1,4,GFA_CUT,
    8,5,10,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=11;
  return GraphFromArray(i,j,T);
}


int goldner_harary(int i,int j){
  int T[]={
    0,1,2,0,3,4,1,5,6,7,2,8,3,
    1,6,2,3,9,6,10,1,9,2,10,GFA_CUT,
    4,9,5,GFA_CUT,
    8,9,7,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=11;
  return GraphFromArray(i,j,T);
}


int fritsch(int i,int j){
  int T[]={GFA_HAM,0,4,6,3,7,2,4,GFA_CUT,1,8,1,7,6,8,5,0,2,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=9;
  return GraphFromArray(i,j,T);
}


int soifer(int i,int j){
  int T[]={GFA_HAM,0,6,4,2,6,0,7,5,8,1,3,8,5,3,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=9;
  return GraphFromArray(i,j,T);
}


int poussin(int i,int j){
  int T[]={
    GFA_HAM,0,2,4,6,8,10,12,14,2,
    0,13,11,9,5,8,12,7,3,14,7,GFA_CUT,
    4,1,9,0,11,GFA_CUT,1,5,GFA_CUT,3,6,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=15;
  return GraphFromArray(i,j,T);
}


int errara(int i,int j){
  int T[]={
    0,GFA_STAR,10,9,5,15,11,16,GFA_CUT,
    1,GFA_STAR,6,16,12,13,4,7,GFA_CUT,
    2,GFA_WHEEL,6,7,8,9,10,GFA_CUT,
    3,GFA_WHEEL,11,12,13,14,15,11,GFA_CUT,
    4,GFA_STAR,7,8,5,14,13,GFA_CUT,
    5,GFA_STAR,8,9,14,15,GFA_CUT,
    16,GFA_STAR,6,11,12,10,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=17;
  return GraphFromArray(i,j,T);
}


int kittell(int i,int j){
  int T[]={
    0,GFA_PATH,22,
    10,12,14,5,13,11,3,10,2,9,1,6,0,
    5,16,6,17,7,1,8,22,20,18,7,GFA_CUT,
    5,15,21,14,22,9,12,GFA_CUT,
    9,14,GFA_CUT,18,8,20,GFA_CUT,
    19,GFA_STAR,15,16,17,21,GFA_CUT,
    2,0,3,0,4,13,4,11,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=23;
  return GraphFromArray(i,j,T);
}


int frucht(int i,int j){
  int T[]={GFA_HAM,0,4,5,3,2,7,6,8,9,11,10,1,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=12;
  return GraphFromArray(i,j,T);
}


int moser(int i,int j){
  int T[]={0,1,2,3,0,4,5,6,0,GFA_CUT,1,3,2,5,4,6,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=7;
  return GraphFromArray(i,j,T);
}


int markstrom(int i,int j){
  int T[]={
    8,0,GFA_PATH,12,5,13,GFA_PATH,7,
    3,21,22,23,0,23,1,2,21,22,18,GFA_CUT,
    15,20,14,13,4,GFA_CUT,
    6,12,7,GFA_CUT,17,19,GFA_CUT,
    9,11,10,16,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=24;
  return GraphFromArray(i,j,T);
}


int robertson(int i,int j){
  int T[]={GFA_HAM,0,4,8,13,17,2,6,10,14,0,
    1,9,16,5,12,1,GFA_CUT,3,11,18,7,15,3,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=19;
  return GraphFromArray(i,j,T);
}


int headwood4(int i,int j){
  int T[]={
    0,GFA_WHEEL,1,5,4,3,13,12,2,GFA_CUT,
    1,GFA_WHEEL,5,6,7,8,2,GFA_CUT,
    2,GFA_WHEEL,8,9,10,11,GFA_CUT,
    16,GFA_WHEEL,4,15,17,19,5,GFA_CUT,
    12,GFA_WHEEL,11,20,19,18,GFA_CUT,
    14,15,3,14,13,18,17,14,18,GFA_CUT,
    22,GFA_WHEEL,6,21,24,23,7,GFA_CUT,
    5,20,6,11,21,10,24,9,23,8,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=25;
  return GraphFromArray(i,j,T);
}


int wiener_araya(int i,int j){
  int T[]={
    0,GFA_STAR,1,4,12,15,GFA_CUT,
    1,GFA_PATH,39,
    41,GFA_STAR,20,23,36,GFA_CUT,
    18,36,35,16,17,1,2,19,GFA_CUT,
    3,21,22,5,4,8,7,26,27,9,10,29,28,39,25,24,6,GFA_CUT,
    8,12,11,31,32,13,14,34,33,37,38,23,GFA_CUT,
    30,40,33,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=42;
  return GraphFromArray(i,j,T);
}


int zamfirescu(int i,int j){
  int T[]={
    0,GFA_PATH,47,
    0,9,5,1,17,1,2,19,18,22,21,20,39,3,4,41,
    40,47,43,42,6,7,44,31,32,28,27,34,33,
    45,46,38,37,21,GFA_CUT,
    8,30,29,10,11,27,26,13,14,24,25,35,
    36,23,16,15,0,12,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=48;
  return GraphFromArray(i,j,T);
}


int hatzel(int i,int j){
  int T[]={
    0,GFA_PATH,55,
    46,30,31,11,12,8,7,26,27,9,10,29,28,
    47,48,44,45,38,39,43,56,52,51,23,22,
    5,6,24,25,50,49,56,GFA_CUT,
    8,4,0,12,13,32,33,37,36,40,41,42,53,
    54,20,21,3,2,19,18,55,41,GFA_CUT,
    14,34,35,16,17,1,0,15,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=57;
  return GraphFromArray(i,j,T);
}


int harborth(int i,int j){
  int T[]={
    0,GFA_PATH,19,GFA_CUT,19,
    0,20,1,21,2,22,3,23,4,24,5,
    25,6,26,7,27,8,28,9,29,10,
    30,11,31,12,32,13,33,14,34,15,
    35,16,36,17,37,18,38,19,39,0,
    39,38,40,41,37,36,35,42,41,43,
    44,45,21,20,45,43,40,39,GFA_CUT,
    46,44,22,23,24,46,25,26,27,
    47,48,28,29,48,49,47,46,GFA_CUT,
    50,49,51,30,31,51,50,32,33,34,
    42,50,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=52;
  return GraphFromArray(i,j,T);
}


int bidiakis(int i,int j){
  int T[]={
    0,GFA_HAM,
    0,6,GFA_CUT,1,5,GFA_CUT,11,7,GFA_CUT,
    10,2,GFA_CUT,9,3,GFA_CUT,8,4,GFA_CUT,
    GFA_END};
  if(j<0) return 0;
  if(i<0) return N=12;
  return GraphFromArray(i,j,T);
}


int cricket(int i,int j){
  int T[]={0,1,2,3,1,4,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=5;
  return GraphFromArray(i,j,T);
}


int moth(int i,int j){
  int T[]={0,1,2,3,4,1,5,1,3,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=6;
  return GraphFromArray(i,j,T);
}


int cross(int i,int j){
  int T[]={0,1,2,3,GFA_CUT,4,1,5,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=6;
  return GraphFromArray(i,j,T);
}


int tgraph(int i,int j){
  int T[]={0,1,2,3,2,4,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=5;
  return GraphFromArray(i,j,T);
}


int bull(int i,int j){
  int T[]={0,1,2,3,4,3,1,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=5;
  return GraphFromArray(i,j,T);
}


int hgraph(int i,int j){
  int T[]={0,1,2,3,GFA_CUT,4,1,2,5,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=6;
  return GraphFromArray(i,j,T);
}


int rgraph(int i,int j){
  int T[]={0,1,2,3,4,1,5,GFA_END};
  if(j<0) return 0;
  if(i<0) return N=6;
  return GraphFromArray(i,j,T);
}


int flower_snark(int i,int j){
  if(j<0){
    if(j==ADJ_NAME) sprintf(NAME,"%c%i",((i&3)==0)? 'c' : 't'+(i&3),i>>2);
    return 0;
  }
  const int k=PARAM[0];
  if(i<0) return N=(k<<2);

  const int u=(j>>2)-(i>>2); /* i<j, donc u>=0 */
  i &= 3;
  j &= 3;
  if(u==0) return (i==0);
  if((u==1)&&(i==j)) return (i>0);
  if(u!=(k-1)) return 0;
  i*=j;
  return ((i==1)||(i==6));
}


int gear(int i,int j){
  if(j<0) return 0;
  int n=PARAM[0];
  if(i<0){
    n<<=1;  /* double le paramètre n */
    PARAM[0]=n; PARAM[1]=PARAM[2]=2; PARAM[3]=0;
    return N=n+1;
  }
  return ((j<n)&&cage(i,j))||((j==n)&&((i&1)==0)); /* suppose i<j */
}


int clebsch(int i,int j){
  if((i<0)||(j<0)) return grid(i,j);
  if(((i|j)==N-1)&&((i&j)==0)) return 1; /* sommets opposés */
  return grid(i,j);
}


int arboricity(int i,int j)
/*
  Utilise REP pour la représentation implicite. REP[i][0..k-1] sont
  les k pères du sommet i. Si REP[i][j]<0, c'est que le père j de i
  n'existe pas (i est une racine de la forêt j par exemple).
*/
{
  if(j<0) return kautz(i,j);

  const int k=PARAM[1];
  int t;

  if(i<0){ /* calcule N et REP[i][0..k-1] */
    int v;
    N=PARAM[0]; if(N<=0) return N=0;
    if(k<1) Erreur(6);   /* il faut k>0 */
    ALLOCMAT(REP,N,k);   /* REP=représentation finale */
    NALLOCZ(int,P,N,_i); /* P=permutation aléatoire */
    NALLOC(int,T,N);     /* T=arbre aléatoire */

    for(t=0;t<k;){
      Dyck(T,N-1,1,DYCK_TREE); /* calcule un arbre aléatoire T */
      T[0]=0; /* racine = 0 plutôt que -1, car P[-1] n'existe pas ! */
      for(v=0;v<N;v++) REP[P[v]][t]=P[T[v]]; /* copie et permute le père */
      REP[P[0]][t]=-1; /* père de la racine = -1 */
      if(++t<k) Permute(P,N); /* si pas fini, on permute aléatoirement X */
    }

    free(T);
    free(P);
    return N;
    }

  for(t=0;t<k;t++) /* j est un père de i ou le contraire */
    if((REP[i][t]==j)||(REP[j][t]==i)) return 1;
  return 0;
}


int kpage(int i,int j)
/*
  Utilise REP pour la représentation implicite. Chaque sommet i
  possède 2k pères: REP[i][2p] et REP[i][2p+1] sont les 2 pères du
  sommet i de la page p (p-ème outerplanar), p=0...k-1. Pour le
  générer, on fait l'union de k outerplanars connexes enracinés et
  plan ("outer-plan" en fait).  Chacun est numérotés selon un parcours
  de la face extérieure avec une permutation circulaire aléatoire,
  sauf le premier. Pour l'adjacence on fait comme pour un graphe
  d'arboricité 2k. Le pagenumber d'un ne peut dépasser ceil(n/2), mais
  des valeurs de k plus grandes sont quand même autorisées.

  Le tirage de chaque outerplanar est construit à partir d'un arbre
  bicolorié (couleur 0 ou 1). Le 2e parent de u existe si u est
  colorié (disons à 1). Il vaut alors le prochain sommet non-related
  de u. On ne tient pas compte de la dernière branche, ces sommets
  n'ayant jamais de 2e parent. Ce tirage est uniforme sur les
  "outer-plan" connexes.
*/
{
  if(j<0){ PARAM[1]>>=1; return kautz(i,j); }
  if(i<0){
    int u,p,q,t,z,k,c;
    N=PARAM[0]; if(N<=0) return N=0;
    k=PARAM[1]; /* k=nombre de pages <= ceil{n/2} */
    if(k<1) Erreur(6);
    PARAM[1]<<=1; /* double le paramètre pour arboricity n 2k */
    ALLOCMAT(REP,N,PARAM[1]);
    NALLOC(int,T,N); /* T=arbre aléatoire */

    /*
      On veut, pour chaque page p:
      REP[u][2p+0]=père1 de u
      REP[u][2p+1]=père2 de u
    */

    for(p=q=0;p<k;p++,q++){ /* pour chaque page p=0...k-1. Attention! q=2*p */
      Dyck(T,N-1,1,DYCK_TREE); /* calcule un arbre DFS aléatoire T */
      
      /* c=permutation circulaire aléatoire pour les noms de sommets */
      if(p) c=random()%N; else c=0; /* aucune permutation pour p=0 */

      /* calcule le père1 des sommets, le père dans T */
      REP[c][q]=-1; /* la racine n'a pas de père */ 
      for(t=1;t<N;t++) /* parcoure les sommets t de T selon le DFS, sauf la racine */
	REP[(t+c)%N][q]=(T[t]+c)%N;

      /* calcule le père2 des sommets */

      /* Principe: le sommet courant est t. Chaque fois qu'on démarre
	 une nouvelle branche (t-T[t]>1), alors on parcoure les
	 sommets u allant de t-1 à T[t] non compris (la branche donc)
	 et décide pour chaque u de le connecter ou pas vers t (t qui
	 est donc le prochain non-related de u). NB: On ne parcoure
	 qu'au plus deux fois chacun des sommets. */

      q++; /* pour le père2 */
      REP[c][q]=-1; /* la racine n'a pas de père */
      for(t=1;t<N;t++){ /* parcoure les sommets t de T selon le DFS, sauf la racine */
	u=(t+c)%N; /* u=sommet du graphe correspondant à t */
	REP[u][q]=-1; /* par défaut, pas de père2 */
	if(t-T[t]>1){ /* ici t démarre une nouvelle branche */
	  z=t-1; /* dernier sommet de la branche précédante */
	  while(z!=T[t]){ /* parcoure la branche précédante */
	    if(random()&01) REP[(z+c)%N][q]=u; /* ajoute un voisin vers u si coloré */
	    z=T[z]; /* descend le long de la branche */
	  }
	}
      }
    }

    DEBUG(
	  for(p=0;p<k;p++){
	    ruling("―",10);printf("\n");
	    printf("page %i:\n",p);
	    for(u=0;u<N;u++)
	      printf("%i: %i %i\n",u,REP[u][2*p],REP[u][2*p+1]);
	  }
	  ruling("―",10);printf("\n");
	  );

    free(T);
    return N;
  }

  return arboricity(i,j);
}


int planar(int i,int j){
/*
  REP[u][0..1] = représentation implicite de u, ces deux pères.
  Attention ! PARAM[1] est modifié (et perdu) à l'initialisation à
  cause de l'utilisation de arboricity().
*/

  if(j<0) return kautz(i,j);
  if(i<0){
    int n=PARAM[0]; // n=nombre de faces
    int f=PARAM[1]; // f=taille maximum des faces internes
    int d=PARAM[2]; // d=degré des sommets internes
    int w=(f<0);    // w=vrai ssi face de taille au plus r

    f=abs(f);
    if((f<3)||(n<=0)) return N=0;

    /* le nombre maximum de sommets du graphe est a priori
       f+(n-1)*(f-2) = n(f-2)+2: 1 cycle de taille f au départ et on
       crée, pour chaque autre face, au plus f-2 nouveaux sommets */

    int k=n*(f-2)+2; /* nb max de sommets */
    if(d<0) d=k; /* pas de contraintes de degré */
    NALLOCZ(int,A,k,-1); /* A[u]=1er voisin de u, aucun par défaut */
    NALLOCZ(int,B,k,-1); /* B[u]=2e voisin de u, aucun par défaut */
    NALLOC(int,C,k); /* C[i]=u=i-ème sommets de la face extérieure */
    NALLOC(int,T,k); /* pour mettre à jour C */
    NALLOCZ(int,D,k,2); /* D[u]=degré du sommet u, 2 par défaut */

    /* initialisation de la 1ère face, la face extérieure */
    int fr=(w)? 3+(random()%(f-2)) : f; /* fr=taille de la 1ère face */
    for(N=0;N<fr;N++){
      A[N]=(N+1)%fr; // père vers le prochain du cycle
      C[N]=N; // face extérieure
    }
    int c=N; // c=|C|=nombre de sommets de la face extérieure */
    int a,b,t,s,p,u,*Z;
    
    /* ajoute toutes n-1 autres faces */
    /* ici N est le nouveau sommet courant */

    while((--n)>0){ // tant qu'il reste une face à faire
      fr=(w)? 3+(random()%(f-2)) : f; // fr=taille de la face
      k=random()%c; // C[k]=un sommet de C au hasard
      A[N]=C[k]; // 1er père de C[k]
      D[C[k]]++; // un voisin de plus pour C[k]

      /* on va tester les voisins valides (degré interne >= d) autour
	 de C[k], puis en choisir un au hasard. Enfin on le connectera
	 par un chemin jusqu'au nouveau sommet courant, N, et mettra à
	 jour la face extérieure. */

      p=min(fr-2,c/2); // p=nombre maximum de sommets de part et
                       // d'autre de C[k] dont il faut tester le degré
      
      /* teste les successeurs de C[k] sur C: C[k+1]...C[k+p] */
      for(a=1;a<=p;a++)	if(D[C[(k+a)%c]]<d) break;
      if(a>p) a--;
      // ici on peut se connecter à n'importe quel sommet entre C[k+1]..C[k+a] 

      /* teste les prédécesseurs de C[k] sur C: C[k-1]...C[k-p] */
      for(b=1;b<=p;b++)	if(D[C[(k+c-b)%c]]<d) break;
      if(b>p) b--;
      // ici on peut se connecter à n'importe quel sommet entre C[k-1]..C[k-b] 

      t=random()%(a+b); // choisit 1 indice parmi [0,a+b[
      if(t<a) s=1; else{ s=-1; t-=a; } 
      p=fr-t-3; // p=nombre de sommets du chemin, N non compris, p=0 possible

      // t random dans [0,a[ si s>0
      // t random dans [0,b[ si s<0
      // s>0 => chemin = C[k]-N-(N+1)- ... - (N+p)-C[k+t+1]
      // s<0 => chemin = C[k-t-1]-(N+p)- ... - (N+1)-N-C[k]
      
      b=(k+s*(t+1)+c)%c; // C[b]=C[k+t+1] ou C[k-t-1]=sommet à connecter
      D[C[b]]++; // un voisin de plus pour C[b]
      for(u=N;u<N+p;u++) B[u]=u+1; // adjacence du chemin
      B[u]=C[b]; // le dernier sommet du chemin pointe sur C[b], u=N possible

      /* mise à jour de C: on supprime d'abord les sommets de C qui
	 vont passer à l'intérieur, puis on calcule dans T la nouvelle
	 face extérieur (en sautant les sommets effacés de C):
	 - si s>0, il faut supprimer C[k+1]...C[k+t]
	 - si s<0, il faut supprimer C[k-t]...C[k-1]
      */
      for(u=1,k+=c;u<=t;u++) C[(k+u*s)%c]=-1; // supprime des sommets de C

      /* créer la nouvelle face extérieure T, a=|T| */
      for(u=a=0;u<=p;u++) T[a++]=N++; // commence par les p+1 nouveaux sommets
      for(u=0,b+=c;u<c;u++) // le reste de l'ancienne face extérieure
	if(C[(b+u*s)%c]>=0) T[a++]=C[(b+u*s)%c];

      SWAP(C,T,Z); // échange C et T
      c=a; // nouvelle taille de C
    }

    free(T);
    free(C);
    free(D);

    /* recopie A,B dans REP */
    ALLOCMAT(REP,N,2);
    for(u=0;u<N;u++){
      REP[u][0]=A[u];
      REP[u][1]=B[u];
    }
    free(A);
    free(B);
    PARAM[1]=2; /* pour utiliser l'arboricité 2 */
    return N;
  }

  return arboricity(i,j);
}


int NextDyck(int *X,int n){
/*
  Calcule le prochain mot de Dyck de longueur 2n contenu dans X (qui
  doit donc être un tableau de taille au moins 2n). Renvoie 1 ssi le
  dernier mot de Dyck a été atteint. Si n<0, alors X est initialisé au
  mot X=(10)^n.  Les mots sont énumérés selon le nombre croissant de 1
  consécutifs à gauche.  L'algorithme est celui de
  https://github.com/cassioneri/Dyck

  Pour n=4:

  10101010 10101100 10110010 10110100 10111000 11001010 11001100
  11010010 11010100 11011000 11100010 11100100 11101000 11110000

*/
  int const m=2*n-1;
  int y=0;
  int x=0;
  int i;

  if(n<0){
    x=1;
    y=-(n<<1);
    for(i=0;i<y;i++){ X[i]=x; x=1-x; }
    return (y==2);
  }

  for(i=m;i>0;i--)
    if(X[i]){
      if(X[i-1]) x++;
      else{
	X[i-1]=1;
	X[i]=0;
	for(y=y-x;y!=0;y--) X[++i]=0;
	while(i<m){
	  X[++i]=1;
	  X[++i]=0;
	}
	return 0;
      }
    }else y++;

  return 1; /* dernier mot atteint */
}


int flip(int i,int j)
/*
  REP[i] = mot de Dyck = [ 1 0 1 1 0 0 ]. Ici n est le nombre de un du
  mot = PARAM[0]-2 = 3.  Le mot représente un arbre binaire complet
  (chaque noeud interne à deux fils exactement). On obtient le codage
  de l'arbre en mot de Dyck par un parcours DFS et en écrivant 1 si
  l'arête parcouru mène à un fils gauche et 0 si elle mène à un fils
  droit. On écrit rien lorsqu'on parcoure les arêtes vers le père.

  Rotation sur le noeud x:

               B   C                      A   B
                \ /                        \ /
             A   o           ->             o   C
              \ /                            \ /
               x                              x

  U=[... 1 A 0 1 B 0 C ...]  ->  V=[... 1 1 A 0 B 0 C ...]

  Algorithme d'adjacence entre les mots U et V:

  1. On cherche le plus grand suffixe commun. Soit i la position tel
     que U[i]<>V[i] et U[i+1...2n-1] = V[i+1...2n-1]. Pour celà on
     remonte de 2n-1 jusqu'à la première différence.

  2. On échange éventuellement U et V de sorte que U[i]=1 et V[i]=0.

  3. On essaye de lire [... 1 A 0 ...] dans V à partir de i et dans U
     à partir de i-1, toujours selon les indices décroissant. En même
     temps qu'on vérifie que V[i]=U[i-1], on calcule la hauteur h avec
     +1 si V[i]=0 et -1 si V[i]=1 (car on lit le mot à l'envers). On
     s'arête dès que h<0. Soit i l'indice dans V où h<0 pour la
     première fois.

  4. On vérifie alors que V[j-1]=1, puis que U[0..j-2]=V[0..j-2]. On
     conclut que U et V sont adjacents.

  https://fr.wikipedia.org/wiki/Nombre_de_Catalan#Chemins_sous-diagonaux_dans_le_carr.C3.A9
  https://en.wikipedia.org/wiki/Tree_rotation
*/
{
  if(j<0){
    if(j==ADJ_END) kautz(i,j);
    if(j==ADJ_NAME) NAME_Vector(REP[i],2*(PARAM[0]-2),"","",1,"%i");
    return 0;
  }
  
  if(i<0){
    int n,t,u;
    /* calcule N =  */
    n=PARAM[0]-2;
    if(n<=0) return N=0;
    t=(n<<1); /* t=2n */
    N=Binom(t,n)/(n+1); /* N=Catalan(n) */
    ALLOCMAT(REP,N,t);
    NALLOC(int,X,t);
    NextDyck(X,-n); /* initialise le 1er mot */
    t *= sizeof(int);
    for(u=0;u<N;u++){
      bcopy(X,REP[u],t); /* copie le mot de Dyck courrant vers REP[u] */
      NextDyck(X,n); /* calcule le mot suivant */
    }
    free(X);
    return N;
  }
  
  int* U=REP[i];
  int* V=REP[j];
  int* W;
  int h=1; /* hauteur de [ ... 1 A 0 ... ] */
  int k=2*(PARAM[0]-2)-1; /* k=dernière position de U ou V */

  /* calcule suffixe */
  while(U[k]==V[k]) k--;
  if(V[k]) SWAP(U,V,W); /* échange U et V */

  /* k=position du 0 dans V */
  while((V[k]==U[k-1])&&(h>0)) h += 1-2*V[--k];

  /* problème ? */
  if((V[--k]==0)||(h>0)) return 0;

  /* préfixe */
  while((k>=0)&&(U[k]==V[k])) k--;
  return (k<0);
}


int linegraph(int i,int j)
/*
  Chaque sommet i possède 2 couleurs prises dans 1...k. Les sommets i
  et j sont adjacents si une couleur de l'un est une couleur de
  l'autre.  Utilise REP pour la représentation des couleurs.
*/
{
  if(j<0) return kautz(i,j);
  if(i<0){
    int u,k;
    N=PARAM[0]; if(N<=0) return N=0;
    k=PARAM[1]; if(k<=0) return k=1;
    ALLOCMAT(REP,N,2);
    for(u=0;u<N;u++){
      REP[u][0]=random()%k;
      REP[u][1]=random()%k;
    }
    return N;
  }

  return ((REP[i][0]==REP[j][0])||(REP[i][0]==REP[j][1])||
	  (REP[i][1]==REP[j][0])||(REP[i][1]==REP[j][1]));
}


int ringarytree(int i,int j)
/*
  Le sommet i est un chemin P(i)=x_1,x_2,... allant de la racine (=0)
  à i, chaque lettre x_t est le numéro du fils, numéro dans [0,r[ pour
  la racine et dans [0,k[ pour les noeuds internes. P(0)={} est vide.

  Pour que i et j soient voisins, avec i<j, il faut:

  - si p=0 (seulement connexion dans l'arbre): P(j)=P(i),x contient
    une seule lettre supplémentaire

  - si p=1 (p=0 et chemin entre noeuds de même niveau):
    P(i)=C,x_1,...,x_k
    P(j)=C,y_1,...,y_k
    y_1=1+x_1
    y_t=0 et x_t=k-1 pour tout t>1

  - si p=2 (p=1 et cycle entre noeuds de même niveau):
    x_1=0 et y_1=r-1 ou k-1 (suivant si C={} ou pas)
    x_t=0 et y_t=k-1 pour tout t>1

    Principe: on calcule P(i) et P(j) en parallèle, lettre par lettre.
*/
{
  if(j<0){
    if(j==ADJ_NAME){
      if(i<=0) sprintf(NAME,"ε");
      else{
	int d=PARAM[2]; /* d=r puis d=k */
	int b=1,t=N,f; /* f=fils */
	char s[NAMEMAX];VIDE(s);
	VIDE(NAME);
	while(i>0){
	  i--,t--,t/=d,f=i/t;
	  sprintf(s,"%i",f); /* écrit le fils */
	  if((t>1)&&(min(PARAM[1],PARAM[2])>9)) strcat(s,","); /* ajoute une "," */
	  strcat(NAME,s); /* ajoute au nom courant */
	  i%=t;
	  if(b) d=PARAM[1],b=0;
	}
      }
    }
    return 0;
  }

  int h=PARAM[0];
  int k=PARAM[1];
  int r=PARAM[2];
  
  if(i<0){
    if((h<0)||(k<0)||(r<0)) return N=0;
    if((h==0)||(r==0)) return N=1;
    if(k==0) PARAM[0]=h=1; /* si k=0 et h>0, alors la hauteur est 1 */
    if(k<=1) N=1+r*h;
    else{ /* ici k>1, h>0 et r>0 */
      int t; /* taille sous-arbre = 1+r+r^2+...+r^{h-1} = (r^h-1)/(k-1) */
      for(t=0,N=1;t<h;t++) N *= k; /* après cette boucle, N=r^h */
      N=1+r*(N-1)/(k-1); /* N=racine + r x (taille sous-arbre) */
    }
    return N;
  }
  
  /* calcule le préfixe commun de P(i) et P(j) */
  
  int fx,fy;
  int x=i,y=j; /* copies de i et j, NB: x<y */
  int t=N; /* t=taille de T, l'arbre où x et y sont */
  int d=r; /* d=nombre de fils de la racine de T, d=r puis d=k */
  int b=1; /* b=1 la 1ère fois, puis b=0 */
  int p=PARAM[3]; /* p=0,1,2 */
  
  /* ici x et y sont dans un arbre T de taille t ayant d fils. La
     taille des fils de T est t'=(t-1)/d. Pour trouver le fils fx de T
     contenant x il suffit de faire (x-1)/t'. Le suffixe de x dans ce
     nouveau sous-arbre est alors (x-1)%t'. */

  fx=fy=0;

  /* tant que x et y sont tout deux dans T (ils viennent du même fils
     fx=fy, et aucun d'eux n'est la racine de T: on calcule alors
     leurs fils fx et fy et met à jour la nouvelle taille de T ainsi
     que le suffixe de x et y dans ce nouveau T. */
  
  while((fx==fy)&&(x>0)&&(y>0)){
    x--,y--,t--; /* on enlève la racine de T */
    t/=d; /* t=taille des fils de T */
    fx=x/t,fy=y/t; /* fils des sous-arbres de x et de y */
    x%=t,y%=t; /* suppression du préfixe de x et de y */
    if(b) d=k,b=0; /* d=#fils des fils de T */
  }

  /*
    Ici la situation est la suivante: x et y sont dans des sous-arbres
    isomorphes à T de taille t, dans les sous-arbres des fils fx et
    fy. Chacun de ces sous-arbres est à d fils. Si fx=fy c'est que x
    est la racine (x=0) car y=0 impossible puisque x<y.

                       o
                   fx / \ fy
                     o   o
                    / \ / \
                     x   y

  */

  /* fx=fy: x et y sont dans le même arbre T */
  if(fx==fy) return ((y-1)%((t-1)/d)==0); /* y fils de x ? */
  if(p==0) return 0;

  /* fx<>fy: x et y sont dans des sous-arbres différents */

  /* si 1er et dernier sous-arbres (voisins dans le cycle), on échange x et y */
  if((p==2)&&(fx==0)){
    b=(t==(N-1)/r)? r:k; /* fx,fy fils de niveau 1 ou > 1 ? */
    if(fy==b-1){ SWAP(x,y,b); fy=1; }
  }
  
  if(fy-fx>1) return 0; /* sous-arbres qui ne sont pas voisins */

  /* x et y sont dans des sous-arbres voisins, chacun de taille t */

  for(;;){
    if((x==0)&&(y==0)) return 1; /* racines voisines */
    if((x==0)||(y==0)) return 0; /* pas même niveau */
    x--,y--,t--,t/=k,fx=x/t,fy=y/t,x%=t,y%=t; /* met à jour x,y,t,fx,fy */
    if((fx!=0)&&(fy!=k-1)) return 0; /* il faut x fils 0 et y fils d-1 */
  }
}


int rarytree(int i,int j)
/*
  Utilise REP. Les sommets sont numérotés selon un DFS modifié: on
  pose les fils avant la récursivité (voir Dyck()).
*/
{

  if(j<0) return kautz(i,j);
  if(i<0){
    int n=PARAM[0]; if(n<=0) return N=0;
    int b=PARAM[1]; if(b<2) return N=0;
    int z=PARAM[2]; if((z!=0)&&(z!=1)) return N=0; /* z=0 ou 1 */
    int *B=Dyck(NULL,n,b-1,DYCK_KTREE); /* B=arbre b-aire aléatoire avec n noeuds internes */
    ALLOCMAT(REP,N=b*n+1+z,1); /* représentation implicite */
    for(b=0;b<N-z;b++) REP[b][0]=B[b]; /* copie l'arbre B dans REP, |B|=N-z */
    free(B);
    if(z) REP[N-1][0]=0; /* le dernier sommet pointe vers la racine 0 */
    PARAM[1]=1; /* pour test d'aboricité k=1 */
    return N;
    }

  return arboricity(i,j);
}


int ktree(int i,int j){
/*
  Utilise REP pour la représentation implicite. REP[i][0...k-1] sont
  les k pères du sommet i. Cette fonction utilise un 3e paramètre
  caché, PARAM[2] qui vaut 0 s'il faut générer un arbre aléatoire et 1
  s'il faut générer un chemin.
*/
  
  if(j<0) return kautz(i,j);
  if(i<0){ /* calcule N et REP[i][0..k-1] */
    int *T;
    int k,t,p,w,x,y;
    k=PARAM[1];
    N=PARAM[0]; if(N<=0) return N=0;
    if(PARAM[2]) ALLOCZ(T,N-k,_i-1); /* chemin de N-k noeuds */
    else T=Dyck(NULL,N-k-1,1,DYCK_TREE); /* arbre de N-k noeuds */
    ALLOCMAT(REP,N,k); /* représentation implicite */

    /* Chacun des k+1 sommets de la racine (numéros de 0 à k) ont pour
       pères tous les autres sommets (=> clique) */

    for(t=0;t<=k;t++) /* pour les k+1 sommets */
      for(p=0;p<k;p++) /* pour les k pères 0..k-1 */
	REP[t][p]=(t+p+1)%(k+1); /* il faut sauter t */

    /* On utilise le fait que les noeuds de T forment un DFS. En
       traitant les sommets dans l'ordre on est sûr que le père est
       déjà traité */
 
    for(t=k+1;t<N;t++){ /* on démarre à k+1, les k+1 sommets de la
			   racine sont déjà traités, tout sommet a donc une racine */
      p=T[t-k]; /* p=noeud père du sommet t du graphe */
      w=random()%(k+1); /* indice d'un des sommets du noeud père qui ne sera pas choisi */
      p += k; /* p=nom réel du père dans le graphe */
      REP[t][0]=p; /* remplit avec le père lui-même */
      x=0; /* x=indice des pères pour REP[p], x=0..k-1 */
      y=(w>0); /* y=prochain indice des pères pour REP[t]. Si w=0 on
		  saute le père */
      while(x<k){
	REP[t][y++]=REP[p][x++];
	if(w==x) y--;
      }
    }

    free(T);  /* libère l'arbre */
    return N;
    }

  return arboricity(i,j); /* même fonction */
}


int apollonian(int i,int j){
/*
  Utilise REP pour la représentation implicite, REP[i][0,1,2] sont les
  3 pères du sommet i. On pourrait factoriser avec polygon().
*/
  if(j<0) return kautz(i,j);
  if(i<0){ /* calcule N et REP[i][0..2] */
    N=PARAM[0]; if(N<4) return N=0; /* il faut au moins 4 sommets */
    PARAM[1]=3; /* pour test d'arboricité, 3 pères */
    const int n=N-3; /* n=nombre de sommets internes */
    const int m=3*n+1; /* nombre de sommets de l'arbre ternaire */
    int *P=Dyck(NULL,n,2,DYCK_KTREE); /* arbre ternaire à n sommets internes */

    /*
      Principe de la construction. On part d'un arbre ternaire à n
      sommets internes (dont la racine), comme dans rarytree(). La
      racine correspond à un K_4 dont le centre est le sommet 3. Puis
      le k-ième noeud interne de l'arbre (donc qui n'est pas une
      feuille) correspond à un nouveau K_4 dont le centre est un
      nouveau sommet numéroté k et qui est connecté à un triangle
      parent. Il y en a trois possibles suivant que le numéro du fils
      où l'on est.

                             0          3    (=triangle:012, centre:3)
                            /|\        /|\
                           1 2 3      4 . .  (=triangle:301, centre:4)
                          /|\        /|\
                         4 5 6      . 5 .    (=triangle:401, centre:5)
                          /|\        /|\
                         7 8 9      . . .

	P = [-,0,0,0,1,1,1,5,5,5] (= père dans l'arbre)
	C = [3,4,-,-,-,5,-,-,-,-] (= centre des triangles des sommets internes de l'arbre)
	
    */

    ALLOCMAT(REP,N,3); /* représentation implicite */
    int u,p,c;

    /* calcule C[u]=numéro du centre du triangle ou 0 si feuille, pour
       tout noeud u de l'arbre */

    NALLOC(int,C,m);
    for(u=1;u<m;u++) C[u]=0,C[P[u]]=1; /* par défaut u est une feuille, et son père non */
    C[0]=c=3; /* racine = centre 3 */
    for(u=1;u<m;u++) if(C[u]) C[u]=++c; /* met le numéro de centre */
    
    for(u=0;u<4;u++)   /* pour le premier tirangle */
      for(c=0;c<3;c++) /* pour les 3 fils de chaque sommet du K_4 */
	REP[u][c]=(u+c+1)%4;

    /* on calcule le triangle REP[c], pour chaque centre c=3...N.
       REP[c][0..2] représente le triangle dont le centre est c */

    for(u=1;u<m;u++){ /* on parcoure les noeuds de l'arbre */
      c=C[u];
      if(c){ /* si u est un noeud interne */
	p=C[P[u]]; /* p=centre du père de u dans l'arbre */
	REP[c][0]=p; /* un sommet du triangle est le centre */
	REP[c][1]=REP[p][u%3]; /* on en prend deux autres parmis le triangle du père */
	REP[c][2]=REP[p][(u+1)%3]; /* en fonction du numéro du fils -> u%3 */
      }
    }

    free(C);
    free(P);
    return N;
    }

  return arboricity(i,j); /* même fonction */
}


int polygon(int i,int j){
/*
  Utilise REP pour la représentation implicite. Chaque sommet u
  possède 2 pères: REP[u][0] et REP[u][1]. Utilise PARAM[1]=2 pour
  l'adjacence: arboricity(i,j).

  La construction est similaire à apollonian(), avec un arbre binaire
  au lieu de ternaire, les K_3 remplaçant les K_4. On pourrait
  factoriser polygon() et apollonian() en passant en paramètre un
  booléen, via PARAM[2] par exemple: 0 pour polygon() et 1 pour
  apollonian(). Pas clair qu'on puisse encore généraliser cette
  construction où apollonian() et polygon() seraient des instances.
*/

  if(j<0) return kautz(i,j);
  if(i<0){
    int u,c,p;
    N=PARAM[0]; if(N<3) return N=0;
    ALLOCMAT(REP,N,2); /* chaque sommet a deux pères */
    PARAM[1]=2; /* pour arboricity(i,j) */
    const int n=N-2; /* n=nombre de sommets internes de l'arbre binaire */
    const int m=2*n+1; /* m=nombre total de sommets de l'arbre */
    int *T=Dyck(NULL,n,1,DYCK_KTREE); /* calcule un arbre binaire aléatoire T */

    /* Principe: on parcours T selon un parcours en profondeur
       modifié, on pose les deux fils avant la récursion.

       Dans l'exemple ci-dessous, la racine correspond au triangle
       (0,a,b). Puis, le fils gauche interne (=1) correspond au
       triangle (1,0,a). Le fils gauche interne 4 correspond au
       triangle (4,1,0). Les autres sommets feuilles ne correspondent
       à aucun triangle. Dans le graphe, les numéros des sommets sont
       décalés de +2.

                             0           
			    / \            a---b      0---1
                           1   2          / \ /      / \ /
                          / \            1---0      3---2
                         3   4            \ /        \ /
                            / \            4          4
                           5   6

       Plus précisément, à chaque nouveau fils u>0 de T qui est un
       fils interne on associe un nouveau sommet c du graphe, un coin
       du triangle. Si u est un fils gauche de v=T[u] alors on
       connecte c aux sommets v et REP[v][0] de G. Si u est un fils
       droit, on connecte c aux sommets v et REP[v][1]. La parité d'un
       noeud détermine s'il s'agit d'un fils droit ou gauche.

	T = [-,0,0,1,1,4,4] (= père dans l'arbre)
	C = [2,3,-,-,-,4,-] (= derniers sommets des triangles des noeuds internes de l'arbre)
    */

    NALLOC(int,C,m);
    for(u=1;u<m;u++) C[u]=0,C[T[u]]=1; /* par défaut u est une feuille, et son père non */
    C[0]=c=2; /* dernier sommet du triangle de la racine */
    for(u=1;u<m;u++) if(C[u]) C[u]=++c; /* met le numéro de centre */
    
    REP[0][0]=REP[0][1]=0;    /* le sommet 0 n'a aucun père */
    REP[1][0]=-1,REP[1][1]=0; /* le sommet 1 a un seul père, 0 */
    REP[2][0]=0,REP[2][1]=1;  /* le sommet 2 a pour pères 0 et 1 */

    for(u=1;u<m;u++){ /* on parcoure les noeuds de l'arbre */
      c=C[u];
      if(c){ /* si u est un noeud interne */
	p=C[T[u]]; /* p=centre du père de u dans l'arbre */
	REP[c][0]=p;
	REP[c][1]=REP[p][u%2];
      }
    }

    free(C);
    free(T);
    return N;
  }

  return arboricity(i,j);
}


int treep(int i,int j){
/*
  Attention ! REP[0][i] = père de i dans l'arbre, et non REP[i][u].
  Les feuilles de cet arbre sans sommet de degré deux sont les sommets
  de numéro < p.
*/

  if(j<0){
    if(j==ADJ_END) if(N>0) free(REP[0]);
    return 0;
  }
  if(i<0){
    const int p=PARAM[0]; /* nb de feuilles */
    if(p<3) return N=0;
    ALLOCMAT(REP,1,(p<<1)-1); /* REP[0] = 1 tableau d'au plus 2p-2 entiers */
    int *P=REP[0]; /* raccourci pour REP[0] */

    /*
      Principe du calcul d'un arbre aléatoire à p feuilles et sans
      sommets de degré deux.

      L'idée est de construire l'arbre à partir des feuilles en
      déterminant le père de ces sommets, tout en garantissant qu'un
      sommet interne soit bien de degré au moins deux (sans son père).

      Soit A la liste des sommets qui n'ont encore pas de père. Cette
      liste de taille au plus p contient initialement toutes les
      feuilles: A={0,1,...,p-1}. Tant que |A|>2 on répète la procédure
      suivante:

      1. On tire un tableau aléatoire R de taille |A| entiers >=0 dont
         la somme fait |A| et possédant au moins une valeur > 1. Pour
         cela on fixe une case aléatoire à 2 puis on répète |A|-2 fois
         R[random()%|A|]++.

      2. Puis on détermine les pères de certains sommets de A en
         fonction de R. L'idée est que si R[k]>1, alors les sommets
         A[u]...A[u+R[k]-1] vont recevoir le même père, un nouveau
         sommets crée et réinjecté à la liste A. Plus précisément, on
         parcoure séquentiellement R. Si R[k]=0, alors on passe
         simplement à k+1 sans rien faire d'autre. Si R[k]=1, on
         maintient A[u] dans la liste A on passe au prochain sommet
         u+1 de A. Le père du sommet A[u] n'est alors toujours pas
         fixé. Par contre, si R[k]>1, alors on crée un nouveau sommet
         v que l'on ajoute à la fin de la liste A, et les R[k] sommets
         A[u]...A[u+R[k]-1] ont alors tous pour père le sommet v.

      Lorsque |A|<=2 on ne rajoute plus aucun sommet, et le nombre de
      sommets N est alors déterminé. On fixe alors que le sommet A[0]
      est racine. Si |A|=2 alors le père de A[1] est fixé à A[0]. Il
      n'est pas possible d'avoir |A|=0 car tout regroupement de sommet
      crée au moins un nouveau sommet dans père (et donc ajouté à A).

      Une complexité de O(p^2) est possible car |A| diminue seulement
      d'une unité si à chaque étape R ne contient qu'une seule valeur
      > 1. Cependant, en moyenne O(log(p)) étapes suffisent, soit
      O(p*log(p)) en tout, car il est facile de voir que R contient
      une fraction de |A| valeurs > 2. (Dans la représentation unaire
      de R il y a n/2 blocks et la moitié sont de longueur > 2.) Et
      pour chacun de tels blocks tous sauf 1 seront enlevé de A.
    */

    int t;
    NALLOCZ(int,A,p,_i); /* tableau des sommets actifs */
    NALLOC(int,R,p); /* tableau des valeurs random */
    int u; /* u=indice dans A du prochain sommet actif */
    int q; /* q=indice dans A du prochain nouveau sommet actif, q<=u */
    int k; /* k=indice dans R, k>=u */
    int a; /* a=taille de A, a<=p */
    N=a=p; /* N=nb courant de sommets dans le graphe, N>=p */

    while(a>2){
      for(k=0;k<a;R[k++]=0); /* tableau R à zéro */
      R[random()%a]=2; /* met un "2" quelque part */
      for(k=2;k<a;k++) R[random()%a]++; /* incrémente a-2 fois des positions de R */
      for(k=u=q=0;k<a;k++){ /* parcoure les valeurs de R */
	if(R[k]==0) continue;
	if(R[k]==1) { A[q++]=A[u++]; continue; }
	t=u+R[k]; /* ici t>=2 */
	for(;u<t;u++) P[A[u]]=N; /* P[A[u]]=père de A[u]=nouveau sommet */
	A[q++]=N++; /* un sommet de plus, et un nouveau actif de plus */
      }
      a=q; /* nouvelle taille de A=nb de nouveaux sommets actifs */
    }

    P[A[0]]=-1;
    if(a==2) P[A[1]]=A[0];

    free(A);
    free(R);
    REP[0]=REALLOC(P,N); /* recalibrage du tableau REP[0] */
    
    return N;
  }

  return ((REP[0][i]==j)||(REP[0][j]==i)); /* arbre */
}


int halin(int i,int j){
/*
  Utilise treep().
*/

  if(j<0) return treep(i,j);
  if(i<0) return N=treep(i,j);
  
  const int p=PARAM[0]; /* nb de feuilles */
  if((i<p)&&(j<p)) return ((j==((i+1)%p))||(i==((j+1)%p))); /* cycle */
  return treep(i,j); /* arbre */
}


int permutation(int i,int j){
/*
  REP[0][i] = permutation du sommet i.
*/
  if(j<0){
    if(j==ADJ_END) if(N>0) free(REP[0]);
    if(j==ADJ_NAME) sprintf(NAME,"(%i,%i)",i,REP[0][i]);
    return 0;
  }
  if(i<0){
    int u;
    N=PARAM[0]; if(N<=0) return N=0;
    ALLOCMAT(REP,1,N); /* REP[0] = 1 tableau de N entiers */
    /* Génère dans une permutation aléatoire de [0,N[ */
    for(u=0;u<N;u++) REP[0][u]=u;
    Permute(REP[0],N);
    return N;
  }

  return ((i-REP[0][i])*(j-REP[0][j])<0);
}


int interval(int i,int j){
/*
  A chaque sommet i correspond un intervalle [a,b] de [0,2N[, avec
  a=REP[i][0] et b=REP[i][1].
*/

  if(j<0){
    if(j==ADJ_END) return kautz(i,j);
    if(j==ADJ_NAME) sprintf(NAME,"[%i,%i]",REP[i][0],REP[i][1]);
    return 0;
  }
  if(i<0){
    N=PARAM[0]; if(N<=0) return N=0;
    const int m=(N<<1);
    int k,x;

    /* génère un intervalle REP[k] pour k, [a,b] dans [0,2N[ avec a<=b */
    ALLOCMAT(REP,N,2);
    for(k=0;k<N;k++){
      x=random()%m;
      REP[k][0]=x;
      REP[k][1]=x+random()%(m-x);
    }
    return N;
  }

  return ( ((REP[i][0]<=REP[j][0])&&(REP[j][1]<=REP[i][1])) ||
	   ((REP[j][0]<=REP[i][0])&&(REP[i][1]<=REP[j][1])) );
}


int sat(int i,int j){
  /*
    [0,2n[: les variables positives et négatives
    [2n+i*k,2n+(i+1)*k[: clause numéro i (i=0 ... PARAM[1])
   */

  if(j<0) return 0;

  const int n=PARAM[0]<<1;
  int k;

  if(i<0){
    if(n<=0) Erreur(6);
    N=n+PARAM[1]*PARAM[2];
    ALLOCMAT(REP,N,1);

    /* chaque sommet-clause est connecté à une variable (positive ou négative) */
    for(k=n;k<N;k++) REP[k][0]=random()%n;

    return N;
  }

  if(i>j) SWAP(i,j,k);
  /* maintenant i<j */

  k=PARAM[2];
  if(j<n) return ((j==i+1)&&(j&01)); /* i-j et j impaire */
  if(i>=n) return (j-i<=k); /* dans la même clique ? */
  return (REP[j][0]==i);
}


int gpetersen(int i,int j){
/*
  u_i dans [0,n[ et v_i dans [n,2n[, N=2n.
  Utilise le fait que i<j.
*/

  if(j<0){
    if(j==ADJ_NAME) sprintf(NAME,"%c%i",(i<N/2)?'u':'v',i%(N/2));
    return 0;
  }
  const int n=PARAM[0];
  if(i<0) return N=max(2*n,0);

  /* u_i-v_i, j>i */
  if(j==(i+n)) return 1;

  /* sinon, par d'arête entre u_i et v_j */
  if((i<n)&&(n<=j)) return 0;

  /* u_i-u_{i+1 mod n}, ici i<j<n */
  if(i<n) return (j==(i+1)%n)||(i==(j+1)%n);

  /* v_i-v_{i+r mod n} */
  const int r=PARAM[1];
  i -= n;
  j -= n;
  /* ici i,j<n mais j<i possible*/
  return (j==(i+r)%n)||(i==(j+r)%n);
}


int antiprism(int i,int j){
  if((i<0)||(j<0)) return gpetersen(i,j);
  return gpetersen(i,j)||(j==PARAM[0]+(i+1)%PARAM[0]);
}


int deltohedron(int i,int j){
/*
  Les sommets de [0,n[ forme le cycle avec n=PARAM[0].
  Utilise i<j.
*/
  if(j<0) return 0;
  const int n=PARAM[0];
  if(i<0){
    if(n<=0) Erreur(6);
    return N=n+2;
  }
  return (j==(i+1)%n) || ((j==n-1)&&(i==0)) || ((j==n)&&(i%2==0)) || ((j==n+1)&&(i%2==1));
}


int kneser(int i,int j){
/*
  REP[i][0...k-1] sont les ensembles représentant les sommets.
*/
  int v,x,y;
  const int k=PARAM[1]; /* ici k>=0 */

  if(j<0){
    if(j==ADJ_END) return kautz(i,j);
    if(j==ADJ_NAME) NAME_Vector(REP[i],k,",","{}",1,"%i");
    return 0;
  }

  if(i<0){
    const int n=PARAM[0];
    N=Binom(n,k); /* on a besoin de N pour allouer REP */
    if(N<=0) return N=0; /* si N=0, fin: graphe vide */
    ALLOCMAT(REP,N,k); /* on connaît N */
    if(N==1) return N=1; /* si N=1, fin: graphe à 1 sommet */
    NALLOC(int,S,k);
    NextSet(S,-1,k); /* premier sous-ensemble */
    for(x=0;x<N;x++){ /* pour tous les sommets x du graphe */
      for(y=0;y<k;y++) REP[x][y]=S[y]; /* copie dans REP[x] */
      NextSet(S,n,k); /* sous-ensemble suivant */
    }
    free(S);
    return N;
  }

  /*
    Calcule si l'intersection possède au plus r éléments. L'algorithme
    ici est en O(k) en utilisant le fait que les éléments de REP sont
    rangés dans l'ordre croissant.
   */

  /* v=nb d'élements commun */
  const int r=PARAM[2];
  v=x=y=0; /* indices pour i et j */

  while((x<k)&&(y<k)&&(v<=r))
    if(REP[i][x]==REP[j][y]) v++,x++,y++;
    else if(REP[i][x]<REP[j][y]) x++; else y++;

  return (v<=r);
}


int rig(int i,int j){
/*
  REP[i][0]=taille t_i de l'ensemble associé au sommet i
  REP[i][1...t_i]=ensemble associé au sommet i
*/

  if(j<0){
    if(j==ADJ_END) return kautz(i,j);
    if(j==ADJ_NAME) NAME_Vector(REP[i]+1,REP[i][0],",","{}",1,"%i");
    return 0;
  }

  int x,y,k,t;

  if(i<0){
    N=PARAM[0]; if(N<=0) return N=0;
    const double p=DPARAM[0];
    k=PARAM[1];
    NALLOC(int,S,k+1); /* ensemble S[1...k] temporaire pour un sommet */
    ALLOC(REP,N);
    for(x=0;x<N;x++){ /* pour chaque sommet x */
      t=0; for(y=1;y<=k;y++) if(RAND01<p) S[++t]=y; /* t=taille de S */
      ALLOC(REP[x],t+1); /* espace pour le sommet x */
      REP[x][0]=t; /* écrit S dans REP[x][1...t] */
      for(y=1;y<=t;y++) REP[x][y]=S[y];
    }
    free(S);
    return N;
  }

  /*
    Détermine si l'intersection de REP[i][1...] et REP[j][1...] est
    vide ou pas.  L'algorithme utilise le fait que les éléments de REP
    sont rangés dans un ordre croissant.
   */

  x=y=1; /* indices pour les ensemble de i et j */
  k=REP[i][0]; /* taille de l'ensemble de i */
  t=REP[i][0]; /* taille de l'ensemble de j */

  while((x<=k)&&(y<=t)){
    if(REP[i][x]==REP[j][y]) return 1;
    if(REP[i][x]<REP[j][y]) x++; else y++;
  }

  return 0;
}

int bdrg(int i,int j){
/*
  PARAM[]={n_1,d_1, ... n_t,d_t,-1}.  Crée dans LOAD un graphe dont la
  distribution des degrés est donnés par PARAM[]. La construction est
  basée sur une permutation aléatoire du tableau T des demi-arêtes
  (matching). Les arêtes du graphe sont alors les arêtes simples entre
  T[2*i] et T[2*i+1]. La taille de T est sum d_i. Si elle n'est pas
  paire, on supprime une demi-arête à un sommet ayant un d_i>0. C'est
  toujours correct, car si d_i=0 pour tous les i, c'est que sum_i d_i
  est paire.

  Considérons l'exemple suivant:

    - 1 sommet  de degré 1 (le sommet 0)
    - 2 sommets de degré 2 (les sommets 1 et 2)
    - 1 sommets de degré 3 (le sommet 3)

    => PARAM[]={1,1,2,2,1,3,-1}
    => T[]=[ 0 1 1 2 2 3 3 3 ] (on répète le sommet u deg(u) fois)
    => permute(T) = [ 2 3 2 0 3 1 3 3 ]
    => arêtes:        2-3,2-0,3-1,3-3
    => arêtes simples: 2-3, 2-0, 3-1
*/

  if(i<0){
    int c=0,u,u0,v,k,p;
    while(PARAM[c++]>=0); /* c=2*t=nombre de valeurs dans PARAM[] sans le -1 */
    c--; if(c&01) Erreur(6); /* il faut un nombre paires de valeurs */
    for(k=N=0;k<c;k+=2) N+=PARAM[k]; /* N=nombre de sommets */
    if(N==0) return 0; /* graphe vide */
    LOAD=new_graph(N); /* crée le graphe (non vide) */
    for(k=u=0;k<c;k+=2){ /* alloue les listes */
      if(PARAM[k+1]>0) u0=u; /* mémorise un sommet avec un degré d_i>0 */
      for(p=0;p<PARAM[k];p++,u++){
	LOAD->d[u]=PARAM[k+1]; /* degré du sommet u */
	ALLOC(LOAD->L[u],LOAD->d[u]); /* liste d'adjacence de u */
      }
    }
    for(v=k=0;k<c;k+=2) v+=PARAM[k]*PARAM[k+1]; /* v=sum d_i*n_i=2*|E| */
    if(v&01){ /* si la somme est impaire, on enlève une demi-arête au sommet u0 */
      LOAD->d[u0]--; /* u0 est nécessairement défini si la somme était impaire */
      REALLOC(LOAD->L[u0],LOAD->d[u0]); /* raccourcie la liste de u0 */
      v--; /* une demi-arête de moins */
    }
    c=v; /* c=nombre de demi-arêtes */
    NALLOC(int,T,c); /* T=tableau des demi-arêtes */
    for(u=k=0;u<N;u++) /* parcoure tous les arcs et remplit T */
      for(p=0;p<LOAD->d[u];p++) T[k++]=u;
    Permute(T,c); /* mélange aléatoire les demi-arêtes */
    memset(LOAD->d,0,LOAD->n*sizeof(int)); /* met à zéro tous les degrés */
    for(p=0;p<c;){ /* parcoure T et ajoute les arêtes simples */
      u=T[p++]; v=T[p++]; if(u==v) continue; /* pas de boucle */
      if(LOAD->d[u]>LOAD->d[v]) SWAP(u,v,k); /* pour une recherche plus rapide */
      if(SetSearch(v,LOAD->L[u],LOAD->d[u],0)<0) /* v dans LOAD->L[u] ? */
	ADD_EDGE(LOAD,u,v); /* non: ajoute u-v et met à jour la taille des listes */
    }
    free(T);
    GraphRealloc(LOAD,LOAD->d);
    return N;
  }

  return load(i,j);
}


int fdrg(int i,int j){
/*
  PARAM[]={n_1,d_1, ... n_t,d_t,-1}.  Crée dans LOAD un graphe dont la
  distribution des degrés est donnés par PARAM[]. La construction est
  basée sur ...
*/

  if(i<0){
    int c=0,u,v,k,p;
    while(PARAM[c++]>=0); /* c=2*t=nombre de valeurs dans PARAM[] sans le -1 */
    c--; if(c&01) Erreur(6); /* il faut un nombre paires de valeurs */
    N=graphical(PARAM,c<<1); /* N=nombre de sommets */
    if(N<0) Erreur(34); /* séquence non graphique */
    N=0; // A FINIR ...
    
    if(N==0) return 0; /* graphe vide */
    LOAD=new_graph(N); /* crée le graphe (non vide) */
    for(k=u=0;k<c;k+=2) /* alloue les listes */
      for(p=0;p<PARAM[k];p++,u++){
	LOAD->d[u]=PARAM[k+1]; /* degré du sommet u */
	ALLOC(LOAD->L[u],LOAD->d[u]); /* liste d'adjacence de u */
      }
    for(v=k=0;k<c;k+=2) v+=PARAM[k]*PARAM[k+1]; /* v=sum d_i */
    return N;
  }

  return load(i,j);
}


int kout(int i,int j){
/*
  REP[i]=tableau des voisins de i. Si REP[i][j]<0, alors c'est la fin
  de la liste.  Attention ! Si i<j, alors REP[i] ne peut pas contenir
  le sommet j (qui a été crée après i). Pour le test d'adjacence, il
  faut donc tester si i est dans REP[j], et pas le contraire !
  On a forcément k>0.
*/
  
  if(j<0) return kautz(i,j);

  int x,y,k=PARAM[1];

  if(i<0){
    int r,d,z;

    N=PARAM[0]; if(N<=0) return N=0;
    ALLOCMAT(REP,N,k);
    NALLOC(int,T,N);

    REP[0][0]=-1;     /* le sommet 0 est seul !*/
    for(x=1;x<N;x++){ /* x=prochain sommet à rajouter */
      r=min(x,k);     /* le degré de x sera au plus r=min{x,k}>0 */
      d=1+random()%r; /* choisir un degré d pour x: d=[1,...,r] */
      for(y=0;y<x;y++) T[y]=y; /* tableau des voisins possibles */
      r=x;              /* r-1=index du dernier élément de T */
      for(y=0;y<d;y++){ /* choisir d voisins dans T */
	z=random()%r;   /* on choisit un voisin parmi ceux qui restent */
        REP[x][y]=T[z]; /* le y-ème voisin de x est T[z] */
        T[z]=T[--r];    /* enlève T[z], met le dernier élément de T à sa place */
      }
      if(d<k) REP[x][d]=-1; /* arrête la liste des voisins de x */
    }

    free(T);
    return N;
  }

  if(i>j) SWAP(i,j,x);
  /* maintenant i<j, donc j a été crée après i */
  
  for(y=0;y<k;y++){
    if(REP[j][y]==i) return 1;
    if(REP[j][y]<0) return 0;
  }
  return 0;
}


int expander(int i,int j){
/*
  REP[i]=tableau des k>0 successeurs de i. Il est possible que le même
  voisin apparaisse plusieurs fois dans REP[i].

  Algorithme:
   1. on part du tableau T[i]=i, pour i=0..n-1
   2. on forme le cycle T[0]-T[1]-...-T[n-1]-T[0]
   3. on permute seulement les n-1 premières valeurs de T
   4. on recommence en 2.

  Il faut répéter k fois la ligne 2, cependant on a pas besoin
  d'effectuer la dernière permutation de la ligne 3. Rem: il est
  inutile de permuter les n cases de T, seule les n-1 première cases
  suffisent car il s'agit de permutations circulaires.
*/

  if(j<0) return kautz(i,j);
  if(i<0){
    int t,c,k;

    N=PARAM[0]; if(N<=0) return N=0;
    k=PARAM[1]; if(k<1) Erreur(6); /* k>0 */
    ALLOCMAT(REP,N,k);
    NALLOCZ(int,T,N,_i); /* tableau pour une permutation aléatoire */

    for(c=0;c<k;){ /* répète k fois, pour chaque cycle */
      for(t=0;t<N;t++) REP[T[t]][c]=T[(t+1)%N]; /* suit et copie le cycle numéro c */
      if(++c<k) Permute(T,N-1); /* permutation aléatoire (inutile le dernier coup) */
      /* seul N-1 éléments ont besoin d'être permutés */
    }

    free(T);
    return N;
  }

  return arboricity(i,j);
}


/********************

  OPERATEURS UNAIRES

  C'est un graphe comme les autres (défini par la fonction d'adjacence
  adj), sauf qu'il prend comme paramètre un autre graphe (qui a pour
  fonction d'adjacence la variable globale ADJ0). L'adjacence adj() se
  calcule donc en fonction de ADJ0(). En général, les opérateurs
  doivent mettre POS=0, sinon les positions XPOS/YPOS pourraient ne
  pas être définies pour certains indices.

  - apex(i,j): ajoute des sommets universels
  - star(i,j): ajoute des sommets de degré 1 à un graphe.

*********************/

adjacence *ADJ0; /* graphe d'origine pour l'opération star() */
int N0; /* nb de sommets du graphe initial */
int *TAB0; /* tableau auxiliaire */

int simule(adjacence *f,int i,int j,int n){
/*
  Appelle la fonction d'adjacence f(i,j) sans modifier la valeur courante N.
    N=nombre de sommets du graphe courant
    f=fonction d'adjacence du graphe d'origine
    n=nombre de sommets du graphe d'origine
*/
  int t,r;
  t=N; /* sauve le N courant */
  N=n; /* met le N du graphe initial, car f peut dépendre de cette valeur */
  r=f(i,j); /* calcul adj(i,j) */
  N=t; /* remet le N du graphe courant */
  return r; /* retourne la valeur */
}

int apex(int i,int j){
/*
  Les sommets de 0..N0-1 sont ceux du graphe initial. Ceux de numéro
  >= N0 sont les apices (donc de degré N0). On se sert du fait que
  i<j. Désactive POS.
*/
  if(i<0){
    N0=simule(ADJ0,i,j,0); /* N0=N du graphe initial */
    N=N0+APEX; /* N=nb de sommets du graphe final, k=nb de sommets ajoutés */
    POS=0; /* il ne faut pas afficher les positions, cela n'a pas de sens */
    return N; /* valeur de retour pour star() */
  } 

  if(j<0){
    simule(ADJ0,i,j,N0);
    return 0;
  }

  if((i<N0)&&(j<N0)) return simule(ADJ0,i,j,N0);
  if((i>=N0)&&(j>=N0)) return 0; /* les sommets ajoutés ne sont pas voisins */
  return 1; /* ici i<N0<=j, j est un apex */
}


int star(int i,int j){
/*
  Les sommets de 0...N0-1 sont ceux du graphe initial. Ceux de numéro
  >= N0 sont ceux de degré 1.  Se sert du fait que i<j. Désactive POS.
*/
  if(i<0){
    int k,t;
    N0=simule(ADJ0,i,j,0); /* N0=N du graphe initial */
    k=STAR; if(k<0) k *= (-N0);
    N=N0+k; /* N=nb de sommets du graphe final, k=nb de sommets ajoutés */
    ALLOC(TAB0,k); /* k=nb de sommets ajoutés */
    for(t=0;t<k;t++){
      if(STAR>0) TAB0[t]=random()%N0;
      else TAB0[t]=t/(-STAR);
    }
    POS=0; /* il ne faut pas afficher les positions, cela n'a pas de sens */
    return N; /* valeur de retour pour star() */
  } 

  if(j<0){
    simule(ADJ0,i,j,N0);
    if(j==ADJ_END) free(TAB0);
    return 0; /* valeur de retour pour star() */
  }

  if((i<N0)&&(j<N0)) return simule(ADJ0,i,j,N0);

  if((i>=N0)&&(j>=N0)) return 0; /* les sommets ajoutés ne sont pas voisins */

  /* ici i<N0<=j */
  return (TAB0[j-N0]==i);
}


/***********************************

        FONCTIONS DU PROGRAMME
              PRINCIPAL

***********************************/

int InitVertex(int n){
/*
  Remplit le tableau V[i] donnant l'étiquette finale du sommet i et
  supprime les sommets suivant la valeur DELV. Utilise aussi SHIFT.
  Si PERMUTE est vrai V[] est remplit d'une permutation aléatoire de
  SHIFT+[0,n[. Si V[i]=-1 cela signifie que i a été supprimé (DELV).
  La fonction retourne le nombre de sommets final du graphe,
  c'est-à-dire le nombre d'étiquettes >=0. Si k sommets ont été
  supprimés, alors les valeurs de V[] sont dans SHIFT+[0,n-k[.

  Initialise aussi le tableau VF[j], avec j=0...n-k, de sorte que
  VF[j]=i si VF[i] est le j-ème sommet non supprimé. Il est important
  que VF[] ait une taille de N au départ. Un realloc() le
  redimensionne plus tard dans la fonction.
*/

  int i,j,k,r; /* r=n-(nb de sommets supprimés) */
  long seuil;

  /* supprime les sommets */
  if(DELV<0.0){ /* ici on en supprime exactement |DELV| sommets */
    if(DELV<-(double)N) DELV=-(double)N;
    for(i=0;i<n;i++) VF[i]=i; /* on va se servir temporairement de VF */
    r=-(int)DELV; /* les r premières valeurs de VF seront les sommets à supprimer */
    for(i=0;i<r;i++){
      j=i+(random()%(n-i));
      SWAP(VF[i],VF[j],k);
    }
    for(i=0;i<r;i++) V[VF[i]]=-1; /* on supprime ces r sommets */
    for(i=r=0;i<n;i++) /* on remplit V et VF */
      if(V[i]>=0) { VF[r]=i; V[i]=r++; }
  }
  else{ /* ici on supprime chaque sommet avec proba DELV */
    seuil=(double)DELV*(double)RAND_MAX;
    for(i=r=0;i<n;i++)
      if(random()<seuil) V[i]=-1;
      else { VF[r]=i; V[i]=r++; }
  } /* dans les deux cas, r=nb de sommets restant */

  /* réajuste le tableau VF à la taille minimum */
  REALLOC(VF,r);

  if(PERMUTE) Permute(V,n);

  /* ne rien faire si SHIFT=0 */
  if(SHIFT>0)
    for(i=0;i<r;i++)
      V[VF[i]] += SHIFT;

  return r;
}


void ScanINC(int *dmax,int *dmin,int *m){
/*
  Calcule, en fonction des tableaux INC[] et VF[], le degré max, le
  degré min et le nombre d'arêtes du graphe (final).
*/
  int d,i;
  *m=*dmax=0;
  *dmin=INT_MAX;
  if(NF<=0) *dmin=0;
  for(i=0;i<NF;i++){
    d=INC[VF[i]]; /* d=degré du i-ème sommet existant */
    *m += d;
    if(d>*dmax) *dmax=d;
    if(d<*dmin) *dmin=d;
  }
  *m >>= 1;
  return;
}


char *MakeCMD(char *s,int deb,int fin){
/*
  Routine permettant de recomposer la ligne de commande. On ajoute à
  la fin de la chaîne s les arguments ARGV[i] pour i=deb à i<fin. Si
  s=NULL alors un pointeur statique sur la chaîne est renvoyée,
  pointeur qui n'a donc pas besoin d'être libéré.

  Si un argument comprend un espace, il est alors parenthésé par '...'
  de façon à être interprété comme un seul argument. Les arguments
  sont séparés par un espace. Le dernier argument est toujours suivi
  d'un espace, éventuellement inutile.
*/
  static char r[CMDMAX];
  if(s==NULL){ s=r; VIDE(s); }
  
  int i;
  for(i=deb;i<fin;i++)
    if(index(ARGV[i],' ')) /* si l'argument est en plusieurs mots */
      strcat(strcat(strcat(s,"'"),ARGV[i]),"' ");
    else strcat(strcat(s,ARGV[i])," ");

  return s;
}


/* pour DateHeure() */
#define DATE_FORMAT "%d/%m/%Y - %Hh%M'%S"
#define SIZE_DATE 22 // avec le 0 final


char *DateHeure(void){
/*
  Renvoie la date et l'heure courante. Il n'est pas (et il ne faut pas
  !) faire de free() sur le pointeur retourné.
*/
  static char date[SIZE_DATE];
  time_t t=time(NULL);
  struct tm *tm=localtime(&t);
  strftime(date,SIZE_DATE*sizeof(char),DATE_FORMAT,tm);

  date[SIZE_DATE-1]=0;
  return date;
}


void Header(int c){
/*
  Affiche et calcule le préambule (date, nom du graphe, nb de
  sommets). Suivant la valeur de c, le nombre d'arêtes est donnée.  Si
  bit-0 de c=1, alors l'entête de base est affichée.  Si bit-1 de c=1,
  alors on affiche le nombre d'arêtes, le degré min et max.

*/

  /* affichage de la date, nom du graphe, ligne de commande et de n */
  if(c&01){
    printf("//\n");
    printf("// %s - seed=%u\n",DateHeure(),SEED);
    printf("// %s\n",MakeCMD(NULL,0,ARGC));
    printf("// n=%i",NF);
  }

  /* affichage du nombre d'arêtes, maxdeg et mindeg */
  if(c&2){
    int maxdeg,mindeg,nbedges;
    ScanINC(&maxdeg,&mindeg,&nbedges);
    if(!(c&01)) printf("//\n//");
    printf(" m=%i",nbedges);
    printf(" maxdeg=%i",maxdeg);
    printf(" mindeg=%i",mindeg);
  }

  if(c) printf("\n//\n");
  return;
}


char *ComputeNAME(int i){
/*
  Fonction servant deux fois dans Out(i,j).  On retourne le pointeur
  NAME après avoir modifier son contenu. Dans le cas du format dot, il
  faut laisser les indices, les noms sont affichés à la fin.
*/
  sprintf(NAME,"%i",(V==NULL)?i:V[i]); /* par défaut NAME="V[i]" */
  if((FORMAT==F_standard)&&(LABEL==1)) adj(i,ADJ_NAME); /* NAME = nom de i */
  if(strlen(NAME)>NAMEMAX) Erreur(17); /* ici NAME<>NULL */
  return NAME;
}


void Out(int i,int j){
/*
  Affiche l'arête i-j suivant le format FORMAT.

  Si i<0 et j<0, alors la fonction Out() est initialisée.
  Si i<0 (et j>=0), alors c'est la terminaison de la fonction Out().
  Si j<0 (et i>=0), alors affiche i seul (sommet isolé).
  Sinon, affiche l'arête i-j

  Si HEADER=1, alors Out() doit se charger d'afficher l'en-tête.

  Si CHECK>0 alors Out() doit se charger de créer et de remplir la
  liste d'adjacence du graphe GF et de déterminer son nombre de
  sommets NF. Pour cela une liste (de type "list") est créee et
  progressivement rempli. A la fin, on calcule GF avec List2Graph().

  Variables globales modifiées:
  - N, GF, NF, VF, NAME
  - CAPTION, NPAL, PALETTE

  Autres variables globales utilisées:
  - CHECK, FORMAT, ROUND, WIDTH
  - HEADER, DIRECTED, VCOLOR
  - XMIN, XMAX, YMIN, YMAX
  - VSIZEK, VSIZESTD, VSIZEXY
  - POS, LABEL, LEN, XPOS, YPOS
  - PARAM_PAL, CPARAM
  - COLORCHAR, COLORBASE
*/
  int x,y,z;
  double w;
  
  static list *L0,*L1,*L2; /* pour -check: tête, dernier, avant-dernier */
  static int cpt;  /* compteur d'arêtes ou de sommets isolés par ligne */
  static int last; /* extrémité de la dernière arête affichée */
  
  /* format de précision par défaut pour l'affichage de XPOS/YPOS dans
     le format F_xy, soient 6 digits par défaut */
  char fmt[17]="%lf %lf\n";

  /*-----------------------------------
    Initialise la fonction (i<0 et j<0)
  -----------------------------------*/
  if((i<0)&&(j<0)){

    cpt=0;
    last=-1;
    if(CHECK) L0=L1=L2=new_list(); /* initialise la liste */

    switch(FORMAT){

    case F_standard:
      if(HEADER) Header(1);
      STRsep1="";
      STRsep2=" ";
      STRedge=(DIRECTED)?"->":"-";
      return;

    case F_dot:
      if(HEADER) Header(1);
      printf("%sgraph {\n",(DIRECTED)?"di":"");
      if(CAPTION){
	printf("\tgraph [label=\"%s\"];\n",CAPTION);
	free(CAPTION);
	CAPTION=NULL;
      }
      if(VCOLOR&0x10){ /* "list" */
	printf("\tgraph [layout=nop, splines=line];\n");
	printf("\tnode [height=1.0, width=0.4, style=filled, shape=rectangle];\n");
	return;
      }
      if(POS){
	w=PosAspect();
	/* layout=nop: pour dire à dot de ne pas re-calculer les positions */
	printf("\tgraph [layout=nop, splines=line, bb=\"%.2lf,%.2lf,%.2lf,%.2lf\"];\n",
	       w*XMIN-2*VSIZEK*VSIZEXY,w*YMIN-2*VSIZEK*VSIZEXY,
	       w*XMAX+2*VSIZEK*VSIZEXY,w*YMAX+2*VSIZEK*VSIZEXY);
	/* si on ne met pas le "2*" on a des sommets tronqués ... */
	if(XYgrid){
	  /* affiche éventuellement une sous-grille, avant le graphe
	     pour qu'elle apparaisse en arrière-plan */
	  if(XYgrid<0){
	    if(XYtype==XY_PERM) XYgrid=N; /* si -xy permutation */
	    else XYgrid=1+(int)(sqrt((double)N));
	  }
	  printf("\n\tsubgraph {\n");
	  printf("\t node [label=\"\", height=0, width=0, color=gray];\n");
	  printf("\t edge [color=gray];");
	  z=N; /* premier sommet de la grille, numéroté après ceux de G */
	  double rx=(double)(XMAX-XMIN)/(double)(XYgrid-1); /* pas de la grille en X */
	  double ry=(double)(YMAX-YMIN)/(double)(XYgrid-1); /* pas de la grille en Y */
	  for(y=0;y<XYgrid;y++){
	      for(x=0;x<XYgrid;x++){
		/* affiche le sommet courant (x,y) ainsi que deux
		   arêtes incidentes vers les voisins (x+1,y) et
		   (x,y+1), s'ils existent */
		printf("\n\t %i [pos=\"%lf,%lf\"];",z,
		       w*(XMIN+(double)x*rx),w*(YMIN+(double)y*ry));
		if(x+1<XYgrid) printf("\t%i--%i;",z,z+1); /* arête vers (x+1,y) */
		if(y+1<XYgrid) printf("\t%i--%i;",z,z+XYgrid); /* arête vers (x,y+1) */
		z++; /* prochain sommet */
	      }
	    }
	  printf("\n}\n\n");
	}
      }
      printf("\tnode [");
      if(LABEL==0) printf("label=\"\", shape=point, "); /* sans label */
      w=POS?VSIZESTD:VSIZEXY; /* taille des sommets */
      printf("height=%.2lf, width=%.2lf];\n",w,w);
      if(strcmp(DOTFILTER,"neato")==0) printf("\tedge [len=%.2lf];\n",LEN);
      STRsep1=";";
      STRsep2="; ";
      STRedge=(DIRECTED)?"->":"--";
    }
    return;
  }

  /*--------------------
    Termine la fonction
  ----------------------*/
  if(i<0){
    if(CHECK){ /* on crée GF en fonction de la liste L0 */
      free(L1); /* supprime la sentienelle (dernier élément) de L0 */
      if(L0==L1) GF=new_graph(0); /* si premier = dernier alors graphe vide */
      else{
	L2->next=NULL; /* coupe la liste à l'avant dernier élément qui a été supprimer */
	GF=List2Graph(L0,4); /* initialise GF et NF */
      }
      NF=GF->n;
    }

    switch(FORMAT){

    case F_standard:
    case F_dot:
      if((VCOLOR&0x10)==0){ /* court-circuite l'affichage du graphe si "-vcolor list" */
	if(cpt>0) printf("%s\n",STRsep1); /* newline si fini avant la fin de ligne */
	if(FORMAT==F_standard){ /* fin si format standard */
	  if(HEADER) Header(2);
	  return;
	}

	if(POS||(LABEL>0)){
	  w=PosAspect();
	  printf("\n");
	  for(y=0;y<NF;y++){
	    x=VF[y]; /* le sommet x existe */
	    printf("%i [",V[x]);
	    if(POS) /* Note: XPOS/YPOS existent forcément si POS=1 */
	      printf("pos=\"%lf,%lf\"",w*XPOS[x],w*YPOS[x]);
	    if(LABEL>0){
	      VIDE(NAME); /* on vide NAME avant de le calculer */
	      if(LABEL==1){
		adj(x,ADJ_NAME); /* calcule le nom dans NAME */
		if(strlen(NAME)>NAMEMAX) Erreur(17);
	      }
	      printf("%slabel=\"%s\"",(POS?", ":""),(NONVIDE(NAME)? NAME : "\\N"));
	    }
	    printf("];\n");
	  }
	}
	
	if(VSIZE&&(NF>0)){ /* taille en fonction du degré */
	  double alpha,smin;
	  ScanINC(&x,&z,&y); /* x=degmax, z=degmin */
	  smin=POS?VSIZESTD:VSIZEXY;
	  alpha=(x==z)? 0.0 : smin*(VSIZEK-1)/((double)(x-z));
	  printf("\n");
	  for(y=0;y<NF;y++){
	    x=VF[y]; /* le sommet x existe */
	    w=smin + alpha*(INC[x]-z);
	    printf("%i [height=%.2lf, width=%.2lf];\n",V[x],w,w);
	  }
	}
      } /* fin du if((VCOLOR&0x10)==0) ... */

      if(VCOLOR&&(NF>0)){ /* couleur des sommets */
	color c={0,0,0},*P; /* couleur noire par défaut */
	int *D;

	if(VCOLOR&0x8){ /* option "pal" on initialise la PALETTE */
	  NPAL=(int)strlen(PARAM_PAL);
	  if(NPAL==1) { /* PARAM_PAL="x" alors PARAM_PAL="xx" */
	    PARAM_PAL[1]=PARAM_PAL[0];
	    PARAM_PAL[NPAL=2]='\0';
	  }
	  ALLOC(PALETTE,NPAL); /* PALETTE = tableau de NPAL "color" */
	  for(y=z=0;y<NPAL;y++){ /* z=prochaine entrée libre dans PALETTE */
	    x=(int)(index(COLORCHAR,PARAM_PAL[y])-COLORCHAR); /* x=index de la couleur */
	    x/=sizeof(*COLORCHAR); /* normalement inutile */
	    if((0<=x)&&(x<COLORNB)) PALETTE[z++]=COLORBASE[x]; /* on a trouvé la couleur */
	  }
	  if(z<2){ PALETTE[0]=PALETTE[1]=c; z=2; } /* si pas trouvé au moins 2 couleurs */
	  NPAL=z; /* NPAL=nb de couleurs trouvées */
	}
	
	if(VCOLOR&0x17){ /* fonction de couleur: 1,2,3,4,5 ou "-vcolor list" */
	  if((VCOLOR&0x7)>2){ /* si 3="degm", 4="randg"  ou 5="kcolor" */
	    if((VCOLOR&0x7)==3){ /* si "degm" */
	      int *T=Prune(GF,NULL);
	      D=GreedyColor(GF,T); /* calcule les couleurs selon l'ordre T */
	      y=1+GF->int1; /* y=nb de couleurs */
	      free(T); /* on a plus besoin de T */
	    }
	    if((VCOLOR&0x7)==4){ /* si "randg" */
	      NALLOCZ(int,T,NF,_i);
	      Permute(T,NF); /* T=permutation aléatoire */
	      D=GreedyColor(GF,T); /* calcule les couleurs selon l'ordre T */
	      y=1+GF->int1; /* y=nb de couleurs */
	      free(T); /* on a plus besoin de T */
	    }
	    if((VCOLOR&0x7)==5){ /* si "kcolor" */
              y=MEM(CPARAM,0,int); /* y=nb de couleur */
	      D=kColor(GF,y);
              if(D==NULL){ ALLOCZ(D,NF,0); y=1; } /* pas trouvé -> une seule couleur */
	    }
	  }
	  else{ /* si 1="deg" ou 2="degr" */
	    ScanINC(&x,&z,&y); /* calcule x=degmax, z=degmin */
	    y=x-z+1; /* y=nb a priori de couleurs nécessaires */
	    ALLOCZ(D,NF,INC[VF[_i]]-z);
	    if((VCOLOR&0x7)==2){ /* si "degr" */
	      int *R=SortInt(D,NULL,N,0,&y,SORT_INC_RANK);
	      /* après SortInt(), y=nb de degré différents */
	      free(D);
	      D=R; /* on substitue R à D */
	    }
	  }
	  /* ici D[x]=indice de la couleur du sommet x, et y=nb de couleurs */
	  P=GradColor(PALETTE,NPAL,y); /* calcule une palette P de y couleurs. NB: ici NPAL>1 */
	  printf("\n");
	  if(VCOLOR&0x10){ /* si "-vcolor list" */
	    for(x=0;x<y;x++){
	      c=P[x];
	      for(z=0;z<COLORNB;z++) /* on cherche c dans COLORBASE */
		if((COLORBASE[z].r==c.r)&&(COLORBASE[z].g==c.g)&&(COLORBASE[z].b==c.b)) break;
	      printf("\t%i [pos=\"%i,0\", color=\"#%.02x%.02x%.02x\", label=\"%c\", fontcolor=%s];\n",x,10+28*x,c.r,c.g,c.b,(z<COLORNB)?COLORCHAR[z]:' ',(c.r+c.g+c.b<150)?"white":"black");
	    }
	  }else{ /* si pas "-vcolor list" */
	    for(x=0;x<NF;x++){
	      c=P[D[x]]; /* c=couleur du degré (ou du rang) de x */
	      printf("%i [style=filled, fillcolor=\"#%02x%02x%02x\"];\n",V[VF[x]],c.r,c.g,c.b);
	    }
	  }

	  free(D);
	  free(P);
	}
	
	if(VCOLOR&0x8){
	  free(PALETTE); /* la PALETTE ne sert plus à rien */
	  PALETTE=COLORBASE; /* remet à la valeur par défaut, au cas où */
	  NPAL=COLORNB; /* taille de la palette par défaut */
	}
      }

      printf("}\n"); /* si F_dot on ferme le "}" */
      if(HEADER) Header(2); /* affiche les arêtes */
      return;
      
    case F_no:
      if(HEADER) Header(3);
      return;

    case F_list:
      if(HEADER) Header(3);
      PrintGraphList(GF);
      return;
      
    case F_matrix:
    case F_smatrix:
      if(HEADER) Header(3);
      PrintGraphMatrix(GF);
      return;

    case F_xy:
      if(HEADER) Header(3);
      if((XPOS==NULL)||(YPOS==NULL)) Erreur(8);
      printf("%i\n",NF); /* nombre de sommets final */
      if(ROUND<DBL_DIG){ /* construit le nouveau format */
	int r=max(ROUND,0); /* si ROUND<0 -> 0 */
	if(ROUND<10){
	  strcpy(fmt,"%.0*lf %.0*lf\n"); /* remplace '*' par r */
	  fmt[3]=fmt[10]=(char)('0'+r);
	}else{
	  strcpy(fmt,"%.0**lf %.0**lf\n"); /* remplace '**' par r */
	  fmt[3]=fmt[11]=(char)('0'+(r/10));
	  fmt[4]=fmt[12]=(char)('0'+(r%10));
	}
      }
      for(y=0;y<NF;y++){
	x=VF[y]; /* le sommet x existe */
	printf(fmt,XPOS[x],YPOS[x]);
      }
      return;

    default: Erreur(5); /* normalement sert à rien */
    }
  }

  /*-----------------------------------------
    Affichage à la volée: "i-j", "i" ou "-j"
  -------------------------------------------*/

  if(CHECK){
    L1=Insert(L2=L1,i,T_NODE); /* on ajoute i */
    if(j>=0) L1=Insert(L2=L1,j,(DIRECTED)?T_ARC:T_EDGE); /* on ajoute ->j ou -j */
    if(VCOLOR&0x10) return; /* ne rien faire d'autre si "-vcolor list".
			       NB: CHECK>0 dans ce cas */
  }

  if((FORMAT==F_standard)||(FORMAT==F_dot)){
    if(j<0) last=-1; /* sommet isolé, donc last sera différent de i */
    if(i!=last) printf("%s%s",(cpt==0)? "" : STRsep2,ComputeNAME(i));
    last=j; /* si i=last, alors affiche -j ou ->j. Si j<0 alors last<0 */
    if(j>=0) printf("%s%s",STRedge,ComputeNAME(j));
    if(++cpt==WIDTH){
      cpt=0; last=-1;
      printf("%s\n",STRsep1);
    }
  } /* si format matrix, smatrix etc., ne rien faire (c'est fait à la fin: i<0) */

  return;
}


double CheckProba(double v){
/*
  Renvoie min(1.0,max(v,0.0)), une probabilité entre [0,1].
*/
  if(v<0.0) return 0.0;
  if(v>1.0) return 1.0;
  return v;
}


/***************************************

           FONCTIONS LIEES
  A L'ANALYSE DE LA LIGNE DE COMMANDE

***************************************/


void Grep(int i){
/*
  Cherche le mot ARGV[i] dans l'aide contenu dans le source du
  programme, puis termine par exit().

  Plusieurs cas peuvent se présenter:
  
  Cas 1: ARGV[i] est une option, ou aucun ARGV[j] avec j<i n'est une
  option.  Alors on affiche l'aide allant de "^....ARGV[i]($| )" à
  "^$".

  Cas 2: ARGV[j] est une option mais pas ARGV[i] avec j<i. Dans ce
  cas, soit mot=ARGV[j]" "ARGV[j+1]" "..." "ARGV[i]. Alors on affiche
  l'aide allant de "[ ]{7}mot($| )" à "^$" ou "^[ ]{7}-fin" avec
  fin défini de sorte que mot ne soit pas un préfixe.

  En fait on bride la recherche de l'option à ARGV[i-1] ou ARGV[i-2]
  seulement, si bien que j=i, i-1 ou i-2.
*/

  int j=i,t,k;

  do{ // calcule j=i, i-1 ou i-2
    if((i>0)&&(*ARGV[i-0]=='-')){ j=i-0; break; }
    if((i>1)&&(*ARGV[i-1]=='-')){ j=i-1; break; }
    if((i>2)&&(*ARGV[i-2]=='-')){ j=i-2; break; }
  }while(0);

  // construit mot

  NALLOC(char,mot,CMDMAX); VIDE(mot);
  for(t=j;t<=i;t++){
    strcat(mot,ARGV[t]);
    if(t<i) strcat(mot," ");
  }

  // construit fin (à faire seulement si j<i)
  // si ARGV[j]="-option abc xy"
  // alors fin="-option ([^a]|a[^b]|ab[^c]abc) ([^x]|x[^y])"

  NALLOC(char,fin,CMDMAX); VIDE(fin);
  if(j<i){
    strcat(fin,ARGV[j]);
    strcat(fin," (");
    t=j+1; k=0;
    while(ARGV[t][k]){
      strncat(fin,ARGV[t],k);
      strcat(fin,"[^");
      strncat(fin,ARGV[t]+k,1);
      strcat(fin,"]|");
      k++;
      if(ARGV[t][k]==0){
	if(t==i) break;
	strcat(fin,ARGV[t]);
	strcat(fin,") (");
	t++; k=0;
      }
    }
    fin[strlen(fin)-1]=')';
  }

  // construit la commande s

  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' ");
  strcat(strcat(s,*ARGV),".c|");
  /* rem: sed -E n'est pas standard */

  if(j==i)
    strcat(strcat(strcat(s,"sed -nE '/^[.]{4}"),mot),"($|[ ])/,/^$/p'|");
  else{
    strcat(strcat(strcat(s,"sed -nE '/^[ ]{7}"),mot),"($|[ ])/,/(^$)|(^[ ]{7}");
    strcat(s,strcat(fin,")/p'|tail -r|sed -n '2,$p'|"));
    // les tail -r permettent de supprimer la dernière ligne
    // le awk est pour enlever éventuellement l'avant dernière ligne "...."
    strcat(s,"awk '{if(NR>1)print $0;else if(!match($0,/^[.]{4}/))print $0;}'|");
    strcat(s,"tail -r|");
  }

  strcat(s,"sed -E 's/^[.]{4}/    /g'");
  strcat(s,"| awk '{n++;print $0;}END{if(!n) ");
  strcat(s,"print\"Erreur : argument incorrect.\\nAide sur ");
  strcat(s,mot);
  strcat(s," non trouvée.\\n\";}'");
  printf("\n");
  system(s);

  if(j<i) printf("\n");
  free(s);
  free(mot);
  free(fin);
  exit(EXIT_SUCCESS);  
}


char *GetArgInc(int *i){
/*
  Retourne ARGV[i], s'il existe, puis incrémente i.  Si l'argument
  n'existe pas ou si ARGV[i]="?", on affiche l'aide en ligne sur le
  dernier argument demandé et l'on s'arrête.
*/
  
  if(((*i)==ARGC)||(strcmp(ARGV[*i],"?")==0)) Grep((*i)-1);
  return ARGV[(*i)++];
}


void CheckHelp(int *i){
/*
  Incrémente i puis vérifie si l'argument est "?". Si tel est le cas,
  une aide est affichée sur ARGV[i] (avant incrémentation).  Cette
  fonction devrait être typiquement appellée lorsqu'on vérifie l'aide
  pour un graphe sans paramètre. Sinon c'est fait par GetArgInc().
*/
  
  (*i)++;
  if(((*i)!=ARGC)&&(strcmp(ARGV[*i],"?")==0)) Grep((*i)-1);
  return;
}


void Help(int i){
/*
  Affiche:
  - l'aide complète si ARGV[i] est "-help" ou "?", ou
  - les paragraphes contenant ARGV[i].
*/
  
  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' "); /* filtre l'aide */
  strcat(strcat(s,*ARGV),".c | "); /* enlève 1ère et dernière ligne */
  strcat(s,"sed -e 's/\\/\\*[#] ###//g' -e 's/### [#]\\*\\///g' ");
  i++;
  if((i==ARGC)||(ARGV[i][0]=='?'))
    strcat(s,"-e 's/^[.][.][.][.][.]*/    /g'|more"); /* aide complète */
  else{
    strcat(s,"|awk 'BEGIN{x=\".....\"}/^[.]{4}./{x=$0} /");
    strcat(strcat(s,ARGV[i]),"/{if(x!=\".....\")print substr(x,5)}'|sort -u");
  }
  system(s);
  free(s);
  exit(EXIT_SUCCESS);
}


void ListGraph(void){
/*
  Affiche les graphes possibles et puis termine.
*/

  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' ");
  strcat(strcat(s,*ARGV),".c| ");
  strcat(s,"sed -e 's/\\/\\*[#] ###//g' -e 's/### [#]\\*\\///g'| ");
  strcat(s,"grep '^[.][.][.][.][^-.]'| sed 's/^[.][.][.][.]//g'");
  system(s);
  free(s);
  exit(EXIT_SUCCESS);
}


void Version(void){
/*
  Affiche la version du programme.
*/

  NALLOC(char,s,CMDMAX); VIDE(s);
  strcat(s,"sed -n '/*[#] ###/,/### #/p' ");
  strcat(strcat(s,*ARGV),".c| ");
  strcat(s,"sed -n 's/.*[-] v\\([0-9]*[.][0-9]*\\) [-].*/\\1/p'");
  system(s);
  free(s);
  exit(EXIT_SUCCESS);
}


void PipeDot(int j){
/*
  Gère l'option "-format dot<xxx>".

  Rem: on pourrait utiliser popen() au lieu de réécrire la ligne de
  commande et de lancer system().
*/
  char type[5];

  CheckHelp(&j);j--; /* vérifie l'aide, ici ARGV[j]="dot<xxx>" */
  strcpy(type,ARGV[j]+3); /* type=<xxx> */
  ARGV[j][3]='\0'; /* ARGV[j]="dot" plutôt que "dot<xxx>" */

  /*
    On réécrit la ligne de commande:
    1. en remplaçant "-format dot<xxx>" par "-format dot"
    2. puis en ajoutant à la fin: "| dot -T<xxx> -K <filter>"
  */
  
  char *s=MakeCMD(NULL,0,ARGC);
  strcat(strcat(strcat(strcat(s,"| dot -T"),type)," -K "),DOTFILTER);
  system(s);
  exit(EXIT_SUCCESS);
}


void Visu(int j){
/*
  Gère l'option "-visu".

  On réécrit la ligne de commande:
  1. en remplaçant "-visu" par "-format dot"
  2. puis en ajoutant à la fin: "| dot -Tpdf -K <filter> -o g.pdf"

  Si le FORMAT est F_no (c'est le cas si l'on a fait "-check maincc"
  ou "loadc" par exemple), alors il y a un problème puisqu'il faut
  qu'un graphe soit généré.
*/
  
  CheckHelp(&j);j--; /* vérifie l'aide, ici ARGV[j]="-visu" */
  if(FORMAT==F_no) Erreur(24);

  /* on reconstruit dans s la ligne de commande */
  char *s=MakeCMD(NULL,0,j);
  strcat(s,"-format dot ");
  MakeCMD(s,j+1,ARGC);
  strcat(strcat(strcat(s,"| dot -Tpdf -K "),DOTFILTER)," -o "xstr(GRAPHPDF));
  system(s);
  exit(EXIT_SUCCESS);
}


void MainCC(int j){
/*
  Gère l'option "-maincc".

  On réécrit la ligne de commande en remplaçant "-maincc" par "-check
  maincc | ./gengraph load - -fast".
*/

  CheckHelp(&j);j--; /* vérifie l'aide, ici ARGV[j]="-maincc" */

  /* on reconstruit dans s la nouvelle ligne de commande */
  char *s=MakeCMD(NULL,0,j);
  strcat(s,"-check maincc | ./gengraph load - -fast ");
  MakeCMD(s,j+1,ARGC);
  system(s);
  exit(EXIT_SUCCESS);
}


int gsub(char *s,const char *t,const char *r){
/*
  Remplace, dans la chaîne s, toutes les occurences de t par r et
  renvoie le nombre de remplacements. La chaîne s est modifiée
  directement. Il est donc nécessaire que s ait suffisament de place.
*/

  const int lr=strlen(r); // longueur de t
  const int lt=strlen(t); // longueur de r
  const int d=lt-lr;

  char *p=s+strlen(s)+1; // pointeur sur la fin de s avec son '\0'
  int n=0; // nombre de remplacements effectués

  while((s=strstr(s,t))){
    n++; // une occurrence de plus
    p -= d; // met à jour le pointeur de fin de s
    memmove(s+lr,s+lt,p-s); // déplace la fin
    memcpy(s,r,lr); // copie le remplacement
  }

  return n;
}


/***********************************

           OPTIONS -ALGO

***********************************/


void PrintMorphism(char *s,int *P,int n)
/*
  Affiche le tableau P de n éléments sous la forme:
  i0->j0   i1->j1 ...
  i8->j8   i9->j8 ...
*/
{
  const int k=8; /* nb de "->" affichés par ligne */
  int i;
  printf("%s",s); /* normalement printf(s) est ok, mais Warning sur certains systèmes */
  for(i=0;i<n;i++)
    printf("%i->%i%s",i,P[i],(((i%k)<(k-1))&&(i<n-1))?"\t":"\n");
  return;
}


/***********************************

           FONCTIONS DE TESTS
              POUR -FILTER

***********************************/


int ftest_minus(graph *G)
/*
  Retourne VRAI ssi G n'est pas dans la famille F, c'est-à-dire G
  n'est isomorphe à aucun graphe de F.
*/
{
  graph *F=MEM(FPARAM,0,graph*);
  if(F==NULL) return (G==NULL);

  int i,*P;
  for(i=0;i<F->f;i++){
    P=Isomorphism(G,F->G[i]);
    free(P);
    if(P!=NULL) return 0;
  }
  
  return 1;
}


int ftest_minus_id(graph *G)
/*
  Retourne VRAI ssi F2 ne contient aucun graphe d'identifiant égale à
  celui de G. La complexité est en O(log|F2|). Il est important que F2
  soit triée par ordre croissant des ID.
*/
{
  graph *F=MEM(FPARAM,0,graph*);
  if(F==NULL) return 0;
  return (bsearch(&G,F->G,F->f,sizeof(graph*),fcmp_graphid)==NULL);
}


int ftest_unique(graph *G)
/*
  Retourne VRAI ssi la sous-famille F allant des indices i+1 à F->f ne contient
  pas G, où i=MEM(FPARAM,0,int).

  Effet de bord: MEM(FPARAM,0,int) est incrémenté.
*/
{
  int i = (MEM(FPARAM,0,int) += 1);
  int *P;

  for(;i<FAMILY->f;i++){
    P=Isomorphism(G,FAMILY->G[i]);
    free(P);
    if(P!=NULL) return 0;
  }

  return 1;
}


int ftest_minor(graph *G)
/*
  Retourne VRAI ssi H est mineur de G.
*/
{
  int *T=Minor(MEM(FPARAM,0,graph*),G);
  if(T==NULL) return 0;
  free(T);
  return 1;
}


int ftest_minor_inv(graph *G)
{
  int *T=Minor(G,MEM(FPARAM,0,graph*));
  if(T==NULL) return 0;
  free(T);
  return 1;
}


int ftest_sub(graph *G)
/*
  Retourne VRAI ssi H est sous-graphe de G avec même nb de sommets.
*/
{
  graph *C=Subgraph(MEM(FPARAM,0,graph*),G);
  if(C==NULL) return 0;
  free_graph(C);
  return 1;
}


int ftest_sub_inv(graph *G)
{
  graph *C=Subgraph(G,MEM(FPARAM,0,graph*));
  if(C==NULL) return 0;
  free_graph(C);
  return 1;
}


int ftest_isub(graph *G)
/*
  Retourne VRAI ssi H est sous-graphe induit de G.
*/
{
  int *C=InducedSubgraph(MEM(FPARAM,0,graph*),G);
  if(C==NULL) return 0;
  free(C);
  return 1;
}


int ftest_isub_inv(graph *G)
{
  int *C=InducedSubgraph(G,MEM(FPARAM,0,graph*));
  if(C==NULL) return 0;
  free(C);
  return 1;
}


int ftest_iso(graph *G)
/*
  Retourne VRAI ssi H est isomorphe à G. Aucun intérêt de faire
  programmer iso-inv.
*/
{
  int *T=Isomorphism(MEM(FPARAM,0,graph*),G);
  if(T==NULL) return 0;
  free(T);
  return 1;
}


int ftest_id(graph *G)
{ return InRange(G->id,FPARAM); }


int ftest_vertex(graph *G)
{ return InRange(G->n,FPARAM); }


int ftest_edge(graph *G)
{ return InRange(NbEdges(G),FPARAM); }


int ftest_degmax(graph *G)
{ return InRange(Degree(G,1),FPARAM); }


int ftest_degmin(graph *G)
{ return InRange(Degree(G,0),FPARAM); }


int ftest_deg(graph *G)
{ int u,b=1,n=G->n;
  for(u=0;(u<n)&&(b);u++) b=InRange(G->d[u],FPARAM);
  return b;
}


int ftest_degenerate(graph *G)
{
  int x;
  int *T=Prune(G,&x);
  free(T);
  return InRange(x,FPARAM);
}


int ftest_gcolor(graph *G)
{
  int *T=Prune(G,NULL);
  int *C=GreedyColor(G,T);
  free(T);
  free(C);
  return InRange(1+G->int1,FPARAM);
}


int ftest_component(graph *G)
{
  param_dfs *p=dfs(G,0,NULL);
  int c=p->nc; /* nb de cc */
  free_param_dfs(p);
  return InRange(c,FPARAM);
}


int ftest_forest(graph *G)
{
  param_dfs *p=dfs(G,0,NULL);
  int c=p->nc; /* nb de cc */
  free_param_dfs(p);
  return InRange(c,FPARAM)&&(NbEdges(G)==G->n-c);
}


int ftest_cutvertex(graph *G)
{
  param_dfs *p=dfs(G,0,NULL);
  int x=p->na;
  free_param_dfs(p);
  return InRange(x,FPARAM);
}


int ftest_biconnected(graph *G)
{
  param_dfs *p=dfs(G,0,NULL);
  int b=(p->nc==1)&&(p->na==0)&&(G->n>2);
  free_param_dfs(p);
  return b;
}


int ftest_ps1xx(graph *G, int version)
{
  path *P=new_path(G,NULL,G->n);
  int v=PS1(G,P,version);
  free_path(P);
  return v;
}


int ftest_ps1(graph *G) { return ftest_ps1xx(G,0); }
int ftest_ps1b(graph *G) { return ftest_ps1xx(G,1); }
int ftest_ps1c(graph *G) { return ftest_ps1xx(G,2); }
int ftest_ps1x(graph *G) { return ftest_ps1xx(G,3); }

int ftest_radius(graph *G)
{
  param_bfs *p;
  int x=G->n;
  int u;
  for(u=0;u<G->n;u++){
    p=bfs(G,u,NULL);
    if(x>=0){
      if(p->n<G->n) x=-1;
      else x=min(x,p->radius);
    }
    free_param_bfs(p);
  }
  return InRange(x,FPARAM);
}


int ftest_girth(graph *G)
{
  param_bfs *p;
  int x=1+G->n;
  int u;
  for(u=0;u<G->n;u++){
    p=bfs(G,u,NULL);
    if(p->cycle>0) x=min(x,p->cycle);
    free_param_bfs(p);
  }
  if(x>G->n) x=-1;
  return InRange(x,FPARAM);
}


int ftest_diameter(graph *G)
{
  param_bfs *p;
  int x=-1;
  int u;
  for(u=0;u<G->n;u++){
    p=bfs(G,u,NULL);
    if(p->n==G->n) x=max(x,p->radius);
    free_param_bfs(p);
  }
  return InRange(x,FPARAM);
}


int ftest_hyperbol(graph *G)
{
  param_bfs *p;
  int n=G->n,h=0;
  int u,v,x,y,d1,d2,d3,t;
  NALLOC(int*,D,n);

  /* calcule la matrice de distances D[][] */
  for(u=0;u<n;u++){
    p=bfs(G,u,NULL);
    D[u]=p->D;
    if(p->n<n){ /* remplace -1 par +infini */
      for(v=0;v<n;v++) if(p->D[v]<0) p->D[v]=INT_MAX;
    }
    p->D=NULL;
    free(p); /* efface p, mais pas p->D */
  }

  /* pour tous les quadruplets {u,v,x,y} */
  for(u=0;u<n;u++)
    for(v=u+1;v<n;v++)
      for(x=v+1;x<n;x++)
	for(y=x+1;y<n;y++){
	  d1=D[u][v]+D[x][y];
	  d2=D[u][x]+D[v][y];
	  d3=D[u][y]+D[v][x];
	  if(d1<d2) SWAP(d1,d2,t);
	  if(d1<d3) SWAP(d1,d3,t);
	  if(d2<d3) d2=d3; /* on se fiche de d3 */
	  if(d1-d2>h) h=d1-d2;
	}

  FREE2(D,n); /* efface la matrice de distances */
  if(h==INT_MAX) h=-1; /* cela ne peut jamais arriver */
  return InRange(h,FPARAM);
}


int ftest_tw(graph *G)
{ return InRange(Treewidth(G,1),FPARAM); }


int ftest_tw2(graph *G)
{ return (Treewidth(G,0)<=2); }


int ftest_rename(graph *G)
{
  G->id=SHIFT++;
  return 1;
}


graph *Filter(graph *F,test *f,int code){
  /*
    Etant donnée une famille de graphes et une fonction de test f,
    renvoie une sous-famille de graphes G de F telle que f(G) est
    vraie (si code=0) ou faux (si code=1). Attention! si on libère F,
    cela détruit la sous-famille renvoyée.

    Effet de bord: si PVALUE est vrai, alors dans les graphes filtrés
    G on met dans G->int1 la valeur du paramètre, CVALUE.
  */
  if((F==NULL)||(F->f==0)) return NULL;
  int i,j,n=F->f;

  graph *R=new_graph(0);
  ALLOC(R->G,n); /* a priori R est de même taille que F */

  for(i=j=0;i<n;i++){
    if(f(F->G[i])^code){
      R->G[j++]=F->G[i];
      if(PVALUE) F->G[i]->int1=CVALUE;
    }
  }

  REALLOC(R->G,j);
  R->f=j;
  return R;
}


graph *Graph2Family(graph *G){
/*
  Renvoie une famille composée d'un seul graphe G.
  Effet de bord: met G->id=0.
*/
  graph *F=new_graph(0);
  ALLOC(F->G,1);
  F->G[0]=G;
  F->f=1;
  G->id=0;
  return F;
}


void ApplyFilter(int code,int index)
/*
  Applique le filtre FTEST (avec le code=0 ou 1) à la famille de
  graphes FAMILY (voir un graphe seul), et affiche la famille
  résultante. On affiche aussi le nombre de graphes obtenus et la
  ligne de commande (sous forme de commentaire). Si index>=0, alors
  ARGV[index] donne le paramètre.

  Effet de bord: FAMILY est libérée.
*/
{
  graph *R;
  int i;

  if(FAMILY->f==0) FAMILY=Graph2Family(FAMILY); /* transforme en famille si graphe simple */
  R=Filter(FAMILY,FTEST,code); /* calcule la sous-famille */

  printf("// #graphs: %i\n// generated by:",(R==NULL)?0:R->f);
  for(i=0;i<ARGC;i++) printf(" %s",ARGV[i]);
  printf("\n");
  if((index>=0)&&(PVALUE)) /* on affiche la valeur de chaque graphe */
    for(i=0;i<R->f;i++)
      printf("[%i] %s: %i\n",R->G[i]->id,ARGV[index],R->G[i]->int1);
  else PrintGraph(R); /* ou bien on affiche le graphe */

  /* ! aux free() de famille de graphes ! */
  free_graph(FAMILY); /* libère en premier la famille */
  free(R->G); /* libère la sous-famille */
  free(R);
  return;
}


void RS_Start(char *nom,char *type,graph *G){
/*
  Partie commune à tous les schémas de routage.
  - nom: est le nom du schéma
  - type: sa catégorie (name-independent ...)
  - G: le graphe sur lequel le schéma doit être appliqué
*/
  printf("\nROUTING SCHEME\n");
  BARRE;
  printf("- name: %s\n",nom);
  printf("- type: %s\n",type);
  printf("- command: %s\n",MakeCMD(NULL,0,ARGC));
  printf("- date: %s\n",DateHeure());
  printf("- seed: %u\n",SEED);
  printf("- time for loading the graph: %s\n",TopChrono(1));
  param_dfs *X=dfs(G,0,NULL);
  int c=X->nc;
  free_param_dfs(X);
  if(c!=1) Erreur(11);
  printf("- checking connectivity: Ok (%s)\n",TopChrono(1));
  int *R=SortGraph(GF,1);
  printf("- checking graph type: simple and undirected (%s)\n",TopChrono(1));
  printf("- #nodes: %i\n",G->n);
  printf("- #edges: %i\n",c=NbEdges(G));
  printf("- average degree: %.2lf\n",(double)(c<<1)/(double)G->n);
  if(!R[6]) Erreur(31); /* le graphe doit être simple */
  if(G->n<2) Erreur(36); /* il faut au moins 1 arête */
  char *s;
  s="?";
  if(HASH==H_PRIME)   s="prime";
  if(HASH==H_SHUFFLE) s="shuffle";
  if(HASH==H_MOD)     s="mod";
  printf("- hash: %s\n",s);
  s="?";
  if(SCENARIO.mode==SC_NONE)   s="none";
  if(SCENARIO.mode==SC_ALL)    s="all";
  if(SCENARIO.mode==SC_NPAIRS) s="n pairs";
  if(SCENARIO.mode==SC_PAIR)   s="pair";
  if(SCENARIO.mode==SC_EDGES)  s="edges";
  if(SCENARIO.mode==SC_ONE)    s="one-to-all";
  if(SCENARIO.mode==SC_UV)     s="u->v";
  printf("- scenario: %s\n",s);
  return;
}


#define EULER_MASCHERONI 0.5772156649015328606 // pour func1()
#define K_MINIMUM        0.3819660112501051518 // = (3-sqrt(5))/2, pour Minimize()


double func1(double k,const void *n){
/*
  Fonction définie par:

                   f(k,n) := 2*n/k + k*(H(k)+1)

  avec H(k) := 1+1/2+1/3+...1/k ~ ln(k)+0.577... + 0.5/k + o(1/k).  Le
  minimum de cette fonction est 2*sqrt(n*ln(n*ln(n))) et atteint pour
  un k ~ 0.5*sqrt(2n/ln(n/ln(n)) ce qui est toujours dans l'intervalle
  [1,n].
*/
  return 2.0*(*((int*)n))/k + k*log(k) + k*(1.0+EULER_MASCHERONI) + 0.5;
}


double Minimize(double (*f)(double,const void*),void *info,double ax,double bx,double tol){
/*
  Calcule et renvoie l'abscisse x0 de l'intervalle [ax,bx] du minimum
  de la fonction f(x,info). La valeur de tolérance tol indique la
  précision souhaitée, tol=0 signifiant la précision maximale. Pour
  cherche un maximum, il suffit de prendre -f(x,info).
  
  L'algorithme: un mélange de la recherche selon le nombre d'or et
  l'interpolation quadratique. L'idée est qu'on choisit deux points
  v<w de [a,b] (au départ a=ax et b=bx), et de réduire la recherche à
  [a,w] ou à [v,b]. Si f(v)<f(w), alors le minimum ne peut être dans
  [v,b] à cause du point w (sauf s'il y a un minimum local). On
  recommence donc avec [a,w] sinon avec [v,b]. Un choix judicieux de v
  et w (basé sur le nombre d'or) permet de ne calculer f() sur qu'un
  seul nouveau point. Avec un point x de [a,b] on peut aussi faire une
  interpolation quatratique (la parabole passant par f(a),f(b) et
  f(x)) et prendre comme nouveau milieu le point le plus bas de la
  parabole dans [a,b]. L'algorithme mélange les deux techniques.

  Voir le "free software optimize.c" du R's project dans:
  https://svn.r-project.org/R/trunk/src/library/stats/src/

  Notes sur math.h:
  - nextafter(x,y)=plus petit double après x en direction de y 
  - fma(x,y,z)=x*y+z
  - fdim(x,y)=|x-y|
*/

  double a,b,d,e,p,q,r,u,v,w,x;
  double fu,fv,fw,fx,xm,tol1,tol2,tol3;

  static double eps=-1; /* on calcule eps qu'une seul fois */
  if(eps<0) eps=sqrt(nextafter(0,1)); /* racine carrée du plus petit double > 0 */
  if(tol<=0.0) tol=1E-10; /* en dessous, cela ne marche pas toujours !?! */
    
  a=ax,b=bx;
  w=v=x=fma(K_MINIMUM,(b-a),a);
  fw=fv=fx=f(x,info);

  tol1=nextafter(1,2); /* le plus petit double > 1 */
  tol3=tol/3.0;
  d=e=0.0;
  
  for(;;){
    xm=(a+b)/2.0;
    tol1=eps*fabs(x)+tol3;
    tol2=2.0*tol1;

    /* critère d'arrêt */
    if(fdim(x,xm)<=tol2-(b-a)/2.0) break;

    p=q=r=0.0;
    if(fabs(e)>tol1){ /* fit parabola */
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=(q-r)*2.0;
      if(q>0.0) p=-p; else q=-q;
      r=e;
      e=d;
    }
    
    if((fabs(p)>=fabs(q*r/2.0))||(p<=q*(a-x))||(p>=q*(b-x))){
      /* étape: recherche nombre d'or */
      e=(x<xm)? b-x : a-x;
      d=K_MINIMUM*e;
    }
    else{ /* étape: interpolation quadratique */
      d=p/q;
      u=x+d; 
      /* u ne doit pas trop près de a ou b */
      if((u-a<tol2)||(b-u<tol2)){ d=tol1; if(x>=xm) d=-d; }
    }
    
    /* on va évaluer f en u, qui ne doit pas être trop près de x */
    if(fabs(d)<tol1) u=(d>0.0)? x+tol1 : x-tol1;
    else u=x+d;
    
    fu=f(u,info);
    
    /* met à jour a,b,v,w,x puis recommence */
    if(fu<=fx){
      if(u<x) b=x; else a=x;
      v=w; fv=fw;
      w=x; fw=fx;
      x=u; fx=fu;
    }else{
      if(u<x) a=u; else b=u;
      if((fu<=fw)||(w==x)){
	v=w; fv=fw;
	w=u; fw=fu;
      }else
	if((fu<=fv)||(v==x)||(v==w)){ v=u; fv=fu; }
    }
    
  }
  
  return x;
}


/***********************************

               MAIN

***********************************/


int main(int argc, char *argv[]){

  ARGC=argc;
  ARGV=argv;

  if(ARGC==1) Help(0);   /* aide si aucun argument */

  /* initialisations, valeurs par défaut */

  TopChrono(0); /* initialise tous les chronomètres internes */
  adj=ring; /* valeur par défaut */
  VIDE(NAME);
  VIDE(PARAM_PAL);
  VIDE(SPARAM);
  /* NB: bzero() n'est pas standard */
  memset(WRAP,0,sizeof(WRAP));
  memset(PARAM,0,sizeof(PARAM));
  memset(SPARAM,0,sizeof(SPARAM));
  memset(DPARAM,0,sizeof(DPARAM));
  srandom(SEED=getpid()^time(NULL)); /* initialise le générateur aléatoire */

  int i,j,k;
  i=1; /* on démarre avec le 1er argument */

  /******************************************************************

                 ANALYSE DE LA LIGNE DE COMMANDE

    o Il faut éviter de faire des traitements trop coûteux en
      temps/mémoire dans l'analyse de la ligne de commande car on peut
      être amené à la refaire une deuxième fois à cause des options
      comme -visu, -maincc, -format dot<xxx> qui causent la réécriture
      puis la ré-analyse de la nouvelle ligne de commande.

    o Il faut éviter d'utiliser random() dans l'analyse de la ligne de
      commande, car si l'option -seed est présente, le comportement
      dépendra si la position de l'option dans la ligne de
      commande. Cependant, dans certain cas comme "caterpillar", ce
      n'est pas évitable.

  ******************************************************************/

  while(i<ARGC){

    if EQUAL("-version") Version();
    if (EQUAL("-help")&&(i+1<ARGC)&&(strcmp(ARGV[i+1],"?")==0)) Grep(i);
    if (EQUAL("-help")||EQUAL("?")) Help(i);
    if EQUAL("-list"){ CheckHelp(&i); ListGraph(); }
    
    j=i; /* pour savoir si on a lu au moins une option ou un graphe */

    /********************/
    /* les options -xxx */
    /********************/

    if EQUAL("-not"){ NOT=1-NOT; goto param0; }
    if EQUAL("-permute"){ PERMUTE=1; goto param0; }
    if EQUAL("-undirected"){ DIRECTED=0; LOOP=0; goto param0; }
    if EQUAL("-directed"){ DIRECTED=1; LOOP=1; goto param0; }
    if EQUAL("-noloop"){ LOOP=0; goto param0; }
    if EQUAL("-loop"){ LOOP=1; goto param0; }
    if EQUAL("-header"){ HEADER=1; goto param0; }
    if EQUAL("-vsize"){ VSIZE=1; goto param0; }
    if EQUAL("-fast"){ FAST=1; goto param0; }
    if EQUAL("-visu") Visu(i);     /* se termine par system() & exit() */
    if EQUAL("-maincc") MainCC(i); /* se termine par system() & exit() */
    if EQUAL("-seed"){ i++;
	srandom(SEED=STRTOI(GetArgInc(&i)));
	goto fin;
      }
    if EQUAL("-width"){ i++;
	WIDTH=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-shift"){ i++;
	SHIFT=STRTOI(GetArgInc(&i));
	if(SHIFT<0) Erreur(6);
	goto fin;
      }
    if EQUAL("-pos"){ i++;
	POS=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-label"){ i++;
	LABEL=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-dotfilter"){ i++;
	DOTFILTER=GetArgInc(&i); /* pointe sur le nom du filtre */
	goto fin;
      }
    if EQUAL("-len"){ i++;
	LEN=STRTOD(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-norm"){ i++;
	NORM=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-delv"){ i++;
	DELV=STRTOD(GetArgInc(&i));
	if(DELV>1.0) DELV=1.0;
	goto fin;
      }
    if EQUAL("-dele"){ i++;
	DELE=CheckProba(STRTOD(GetArgInc(&i)));
	goto fin;
      }
    if EQUAL("-redirect"){ i++;
	REDIRECT=CheckProba(STRTOD(GetArgInc(&i)));
	goto fin;
      }
    if EQUAL("-star"){ i++;
	STAR=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-apex"){ i++;
	APEX=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-caption"){ i++;
	char *c=GetArgInc(&i); // lecture de la légende
	if(CAPTION) free(CAPTION); // si CAPTION déjà allouée
	char *s=strdup(c);
	k=gsub(s,"%SEED","%u");
	if(k>1) Erreur(35);
	if(k==1) asprintf(&CAPTION,s,SEED);
	if(k==0) CAPTION=s;
	goto fin;
      }
    if EQUAL("-variant"){ i++;
	VARIANT=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("-vcolor"){ i++;GetArgInc(&i);i--; /* pour l'aide en ligne */
	/* bits 0-2: codage de la fonction de couleur=1,2,3,4,5
	   bit 3: "pal"
	   bit 4: "list" */
	if EQUAL("deg"){ i++; VCOLOR=(VCOLOR|0x7)^0x7; /* efface les 3 derniers bits */
	    VCOLOR |= 1; goto fin; }
	if EQUAL("degr"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 2; goto fin; }
	if EQUAL("degm"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 3; CHECK=max(CHECK,CHECK_ON); goto fin; }
	if EQUAL("randg"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 4; CHECK=max(CHECK,CHECK_ON); goto fin; }
	if EQUAL("kcolor"){ i++; VCOLOR=(VCOLOR|0x7)^0x7;
	    VCOLOR |= 5; CHECK=max(CHECK,CHECK_ON);
	    if(CPARAM==NULL) ALLOC(CPARAM,PARAMSIZE);
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i)); /* CPARAM[0]=ARGV[i] */
	    goto fin; }
	if EQUAL("pal"){ i++;
	    VCOLOR |= 0x8; /* set bit-3 */
	    GetArgInc(&i);i--; /* teste si le prochain argument après "pal" existe bien */
	    if(strlen(ARGV[i])>=(int)(sizeof(PARAM_PAL)/sizeof(char))) Erreur(20);
	    strcpy(PARAM_PAL,ARGV[i++]); /* PARAM_PAL=ARGV[i] */
	    goto fin; }
	if EQUAL("list"){ i++;
	    VCOLOR |= 0x10;
	    CHECK=max(CHECK,CHECK_ON);
	    FORMAT = F_dot;
	    goto fin; }
	Erreur(9);
      }
    if EQUAL("-format"){ i++;
	GetArgInc(&i);i--; /* teste si prochain argument existe */
	FORMAT=-1; /* sentiennelle pour savoir si on a trouvé le FORMAT */
	if EQUAL("standard") FORMAT=F_standard;
	if EQUAL("matrix") { FORMAT=F_matrix; CHECK=max(CHECK,CHECK_ON);}
	if EQUAL("smatrix"){ FORMAT=F_smatrix;CHECK=max(CHECK,CHECK_ON);}
	if EQUAL("list")   { FORMAT=F_list;   CHECK=max(CHECK,CHECK_ON);}
	if EQUAL("xy")       FORMAT=F_xy;
	if EQUAL("no")       FORMAT=F_no;
	if PREFIX("dot"){ /* si "dot" ou "dot<xxx>" */
	    if EQUAL("dot") FORMAT=F_dot; /* si "dot" seul */
	    else PipeDot(i); /* se termine par system() & exit() */
	  }
	if(FORMAT<0) Erreur(5); /* le format n'a pas été trouvé */
	i++;
	goto fin;
      }
    if EQUAL("-xy"){ i++;GetArgInc(&i);i--; /* pour l'aide en ligne */
	POS=1; /* il faut POS=1 */
	if EQUAL("load"){ i++;
	    XYtype=XY_FILE;
	    FILEXY=GetArgInc(&i); /* pointe sur le nom du fichier */
	    goto fin;
	  }
	if EQUAL("noise"){ i++;
	    NOISEr=STRTOD(GetArgInc(&i));
	    NOISEp=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("box"){ i++;
	    BOXX=STRTOD(GetArgInc(&i));
	    BOXY=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("seed"){ i++;
	    SEEDk=abs(STRTOI(GetArgInc(&i)));
	    SEEDp=STRTOD(GetArgInc(&i));
	    XYtype=XY_PLAW;
	    goto fin;
	  }
	if EQUAL("round"){ i++;
	    ROUND=min(STRTOI(GetArgInc(&i)),DBL_DIG);
	    goto fin;
	  }
	if EQUAL("permutation"){ CheckHelp(&i);
	    ROUND=0; /* a priori coordonnées entières dans ce cas */
	    XYtype=XY_PERM;
	    goto fin;
	  }
	if EQUAL("mesh"){ CheckHelp(&i);
	    ROUND=0; /* a priori coordonnées entières dans ce cas */
	    XYtype=XY_MESH;
	    Xmesh=STRTOI(GetArgInc(&i));
	    Ymesh=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("unique"){ CheckHelp(&i);
	    XYunique=1;
	    goto fin;
	  }
	if EQUAL("grid"){ CheckHelp(&i);
	    XYgrid=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("vsize"){ CheckHelp(&i);
	    XYvsize=STRTOD(GetArgInc(&i));
	    goto fin;
	  }
	Erreur(1); /* l'option après -xy n'a pas été trouvée */
      }
    if EQUAL("-filter"){ i++;GetArgInc(&i);i--; /* pour l'aide en ligne */
	FAMILY=File2Graph(GetArgInc(&i),2); /* lit une famille ou un graphe */
	if(FPARAM==NULL) ALLOC(FPARAM,PARAMSIZE);
	GetArgInc(&i);i--; /* vérifie s'il y a bien un autre argument */
	PVALUE=0; /* par défaut, on affiche pas "value" mais les graphes */
	if EQUAL("not"){ k=1; i++;GetArgInc(&i);i--; } else k=0;
	 /* vérifie s'il y a bien un autre argument */
	if EQUAL("rename"){
	    FTEST=ftest_rename;
	    i++;
	    SHIFT=STRTOI(GetArgInc(&i));
	    ApplyFilter(0,-1);
	    goto fin;
	  }
	if EQUAL("biconnected"){
	    FTEST=ftest_biconnected;
	  filter0:
	    i++;
	    ApplyFilter(k,-1);
	    goto fin;
	  }
	if EQUAL("id"){
	    FTEST=ftest_id;
	  filter1:
	    i++;
	    ReadRange(GetArgInc(&i),FPARAM);
	    ApplyFilter(k,i-2);
	    goto fin;
	  }
	int c; /* code pour File2Graph() */
	if EQUAL("minor"){
	    FTEST=ftest_minor;
	    c=34; /* détection du shift et charge toujours un graphe */
	  filter2:
	    i++;
	    MEM(FPARAM,0,graph*)=File2Graph(GetArgInc(&i),c);
	    ApplyFilter(k,-1);
	    free_graph(MEM(FPARAM,0,graph*));
	    goto fin;
	  }

	if EQUAL("ps1x"){ i++;
	    FTEST=ftest_ps1x;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    for(c=MEM(CPARAM,0,int);c>=1;c--){ /* met les arguments à l'envers */
	      MEM(CPARAM,(2*c-1)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	      MEM(CPARAM,(2*c)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	    }
	    i--;
	    goto filter0;
	  }
	if EQUAL("ps1"){ FTEST=ftest_ps1; goto filter0; }
	if EQUAL("ps1b"){ FTEST=ftest_ps1b; goto filter0; }
	if EQUAL("ps1c"){ FTEST=ftest_ps1c; goto filter0; }
	if EQUAL("tw2"){ FTEST=ftest_tw2; goto filter0; }
	if EQUAL("unique"){ FTEST=ftest_unique; MEM(FPARAM,0,int)=0; goto filter0; }
	if EQUAL("vertex"){ FTEST=ftest_vertex; goto filter1; }
	if(EQUAL("edge")||EQUAL("edges")){ FTEST=ftest_edge; goto filter1; }
	if EQUAL("deg"){ FTEST=ftest_deg; goto filter1; }
	if EQUAL("degenerate"){ FTEST=ftest_degenerate; goto filter1; }
	if EQUAL("degmax"){ FTEST=ftest_degmax; goto filter1; }
	if EQUAL("degmin"){ FTEST=ftest_degmin; goto filter1; }
	if EQUAL("gcolor"){ FTEST=ftest_gcolor; goto filter1; }
	if EQUAL("component"){ FTEST=ftest_component; goto filter1; }
	if EQUAL("cut-vertex"){ FTEST=ftest_cutvertex; goto filter1; }
	if EQUAL("radius"){ FTEST=ftest_radius; goto filter1; }
	if EQUAL("girth"){ FTEST=ftest_girth; goto filter1; }
	if EQUAL("diameter"){ FTEST=ftest_diameter; goto filter1; }
	if EQUAL("hyperbol"){ FTEST=ftest_hyperbol; goto filter1; }
	if EQUAL("tw"){ FTEST=ftest_tw; goto filter1; }
	if EQUAL("forest"){ FTEST=ftest_forest; goto filter1; }
	if EQUAL("minor-inv"){ FTEST=ftest_minor_inv; c=34; goto filter2; }
	if EQUAL("sub"){ FTEST=ftest_sub; c=34; goto filter2; }
	if EQUAL("sub-inv"){ FTEST=ftest_sub_inv; c=34; goto filter2; }
	if EQUAL("isub"){ FTEST=ftest_isub; c=34; goto filter2; }
	if EQUAL("isub-inv"){ FTEST=ftest_isub_inv; c=34; goto filter2; }
	if EQUAL("iso"){ FTEST=ftest_iso; c=34; goto filter2; }
	if EQUAL("minus"){ FTEST=ftest_minus; c=2; goto filter2; }
	if EQUAL("minus-id"){ FTEST=ftest_minus_id; c=10; goto filter2; }
	/* alias */
	if EQUAL("connected"){ FTEST=ftest_component;
	    ReadRange("1",FPARAM); goto filter0; }
	if EQUAL("bipartite"){ FTEST=ftest_gcolor;
	    ReadRange("<3",FPARAM); goto filter0; }
	if EQUAL("isforest"){ FTEST=ftest_forest;
	    ReadRange("t",FPARAM); goto filter0; }
	if EQUAL("istree"){ FTEST=ftest_forest;
	    ReadRange("1",FPARAM); goto filter0; }
	if EQUAL("cycle"){ FTEST=ftest_forest;
	    k=1-k; ReadRange("t",FPARAM); goto filter0; }
	if EQUAL("all"){ FTEST=ftest_vertex;
	    ReadRange("t",FPARAM); goto filter0; }
	Erreur(14); /* l'option après -filter n'a pas été trouvée */
      }
    if EQUAL("-check"){ i++;GetArgInc(&i);i--; /* pour l'aide en ligne */
	if(CHECK>CHECK_ON) Erreur(27); /* on ne devrait jamais avoir deux fois -check */
	if(CPARAM==NULL) ALLOC(CPARAM,PARAMSIZE); /* alloue les paramètres */
	if EQUAL("bfs"){
	    CHECK=CHECK_BFS;
	  check0:
	    i++;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("iso"){
	    CHECK=CHECK_ISO;
	  check1:
	    MEM(CPARAM,0,int)=++i;
	    GetArgInc(&i); /* pour vérifier si i existe */
	    goto fin;
	  }
	if(EQUAL("deg")||EQUAL("edge")||EQUAL("edges")){
	  CHECK=CHECK_DEG;
	check_fin:
	  i++;
	  goto fin;
	}
	if EQUAL("paths"){
	    CHECK=CHECK_PATHS;
	  check3:
	    i++;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    MEM(CPARAM,sizeof(int),int)=STRTOI(GetArgInc(&i));
	    goto fin;
	  }
	if EQUAL("ps1x"){ i++;
	    CHECK=CHECK_PS1x;
	    MEM(CPARAM,0,int)=STRTOI(GetArgInc(&i));
	    for(k=MEM(CPARAM,0,int);k>=1;k--){ /* met les arguments à l'envers */
	      MEM(CPARAM,(2*k-1)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	      MEM(CPARAM,(2*k)*sizeof(int),int)=STRTOI(GetArgInc(&i));
	    }
	    goto fin;
	  }
	if EQUAL("routing"){ i++;GetArgInc(&i);i--;
	    FORMAT=F_no; /* pour tous les routing schemes */
	    HASH=H_PRIME; /* hash par défaut */
	    if EQUAL("hash"){ i++;GetArgInc(&i);i--;
		for(;;){ /* pour faire break (importants à cause des i++) */
		  if EQUAL("prime")  { i++; HASH=H_PRIME; break; }
		  if EQUAL("shuffle"){ i++; HASH=H_SHUFFLE; break; }
		  if EQUAL("mod")    { i++; HASH=H_MOD; break; }
		  Erreur(2); /* option après "hash" non trouvé */
		}
	      }
	    SCENARIO.mode=SC_NONE; /* par défaut aucun scenario */
	    if EQUAL("scenario"){ i++;GetArgInc(&i);i--;
		for(;;){ /* pour faire break (importants à cause des i++) */
		  if EQUAL("none")  { i++; SCENARIO.mode=SC_NONE; break; }
		  if EQUAL("all")   { i++; SCENARIO.mode=SC_ALL; break; }
		  if EQUAL("npairs"){ i++; SCENARIO.mode=SC_NPAIRS; break; }
		  if EQUAL("edges") { i++; SCENARIO.mode=SC_EDGES; break; }
		  if EQUAL("one")   { i++;
		      SCENARIO.mode=SC_ONE;
		      SCENARIO.u=STRTOI(GetArgInc(&i));
		      break;
		    }
		  if EQUAL("pair"){ i++;
		      SCENARIO.u=STRTOI(GetArgInc(&i));
		      if(SCENARIO.u>=0){
			SCENARIO.mode=SC_UV;
			SCENARIO.v=STRTOI(GetArgInc(&i));
		      }else{
			SCENARIO.mode=SC_PAIR;
			SCENARIO.u=-SCENARIO.u;
		      }
		      break;
		    }
		  Erreur(2); /* option après "scenario" non trouvé */
		}
	      }
	    if EQUAL("cluster"){ CHECK=CHECK_RS_CLUSTER; goto check0; }
	    if EQUAL("dcr"){ CHECK=CHECK_RS_DCR; goto check0; }
	    if EQUAL("tzrplg"){ i++;
		CHECK=CHECK_RS_TZRPLG;
		MEM(CPARAM,0,double)=STRTOD(GetArgInc(&i));
		goto fin;
	      }
	    if EQUAL("bc"){ CHECK=CHECK_RS_BC; goto check0; }
	  }
	if EQUAL("dfs"){ CHECK=CHECK_DFS; goto check0; }
	if EQUAL("bellman"){ CHECK=CHECK_BELLMAN; goto check0; }
	if EQUAL("kcolor"){ CHECK=CHECK_KCOLOR; goto check0; }
	if EQUAL("kcolorsat"){ CHECK=CHECK_KCOLORSAT; FORMAT=F_no; goto check0; }
	if EQUAL("kindepsat"){ CHECK=CHECK_KINDEPSAT; FORMAT=F_no; goto check0; }
	if EQUAL("sub"){ CHECK=CHECK_SUB; goto check1; }
	if EQUAL("isub"){ CHECK=CHECK_ISUB; goto check1; }
	if EQUAL("minor"){ CHECK=CHECK_MINOR; goto check1; }
	if EQUAL("degenerate"){ CHECK=CHECK_DEGENERATE; goto check_fin; }
	if EQUAL("gcolor"){ CHECK=CHECK_GCOLOR; goto check_fin; }
	if EQUAL("ps1"){ CHECK=CHECK_PS1; goto check_fin; }
	if EQUAL("ps1b"){ CHECK=CHECK_PS1b; goto check_fin; }
	if EQUAL("ps1c"){ CHECK=CHECK_PS1c; goto check_fin; }
	if EQUAL("twdeg"){ CHECK=CHECK_TWDEG; goto check_fin; }
	if EQUAL("tw"){ CHECK=CHECK_TW; goto check_fin; }
	if EQUAL("girth"){ CHECK=CHECK_GIRTH; goto check_fin; }
	if EQUAL("info"){ CHECK=CHECK_INFO; FORMAT=F_no; goto check_fin; }
	if EQUAL("maincc"){ CHECK=CHECK_MAINCC; FORMAT=F_no; goto check_fin; }
	if(EQUAL("ncc")||EQUAL("connected")){ CHECK=CHECK_NCC; goto check_fin; }
	if EQUAL("diameter"){ CHECK=CHECK_DIAMETER; goto check_fin; }
	if EQUAL("radius"){ CHECK=CHECK_RADIUS; goto check_fin; }
	if EQUAL("paths2"){ CHECK=CHECK_PATHS2; goto check3; }
	Erreur(12); /* l'option après -check n'a pas été trouvée */
      }

    /*******************/
    /* graphes de base */
    /*******************/

    if EQUAL("tutte"){
	adj=tutte;
      param0:
	CheckHelp(&i); /* CheckHelp() au lieu de i++ car aucun paramètre */
	goto fin;
      }
    if EQUAL("prime"){
	adj=prime;
      param1:
	i++;
	PARAM[0]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("arboricity"){
	adj=arboricity;
      param2:
	i++;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("sat"){
	adj=sat;
      param3:
	i++;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=STRTOI(GetArgInc(&i));
	PARAM[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("ringarytree"){
	adj=ringarytree;
	param4:
	i++;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=STRTOI(GetArgInc(&i));
	PARAM[2]=STRTOI(GetArgInc(&i));
	PARAM[3]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("ktree"){
	PARAM[2]=0;
      param_ktree:
	i++;
	adj=ktree;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=STRTOI(GetArgInc(&i));
	if(PARAM[0]<=PARAM[1]) Erreur(6);
	goto fin;
      }
    if EQUAL("hypercube"){
	adj=grid;
      param_hypercube:
	i++;
	PARAM[0]=STRTOI(GetArgInc(&i));
	if(PARAM[0]+1>DIMAX) Erreur(4);
	for(k=1;k<=PARAM[0];k++) PARAM[k]=2;
	goto fin;
      }
    if EQUAL("ggosset"){
	adj=ggosset;
	i++;
	int k,d=0;
	PARAM[0]=STRTOI(GetArgInc(&i)); /* copie p */
	PARAM[1]=STRTOI(GetArgInc(&i)); /* copie k */
	for(k=0;k<(PARAM[1]<<1);k++) PARAM[k+3]=STRTOI(GetArgInc(&i));
	/* il est important de calculer la dimension d pour accélérer l'adjacence */
	for(k=0;k<PARAM[1];k++) d += PARAM[(k<<1)+3];
	PARAM[2]=d;
	goto fin;
      }
    if EQUAL("turan"){
	adj=rpartite;
	i++;
	int k,n,r;
	n=STRTOI(GetArgInc(&i));
	r=STRTOI(GetArgInc(&i));
	int p=n%r; /* p=#de part de taille n2 */
	int n1=n/r; /* n1=floor(n/r) */
	int n2=n1;
	if(p) n2++; /* n2=ceil(n/r) */
	PARAM[0]=r;
	for(k=1;k<=r;k++) PARAM[k]=(k<=p)? n1 : n2;
	goto fin;
      }
    if EQUAL("grid"){
	adj=grid;
      param_grid:
	i++;
	PARAM[0]=STRTOI(GetArgInc(&i));
	if(PARAM[0]+1>DIMAX) Erreur(4);
	for(k=1;k<=PARAM[0];k++)
	  PARAM[k]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("ring"){
	adj=ring;
      param_ring:
	i++;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=STRTOI(GetArgInc(&i));
	if(PARAM[1]+2>DIMAX) Erreur(4);
	for(k=0;k<PARAM[1];k++)
	  PARAM[k+2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("kout"){ i++;
	adj=kout;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=STRTOI(GetArgInc(&i));
	if(PARAM[1]>=PARAM[0]-1) PARAM[1]=PARAM[0]-1;
	if(PARAM[1]<1) Erreur(6);
	goto fin;
      }
    if EQUAL("kneser"){ i++;
	adj=kneser;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=STRTOI(GetArgInc(&i));
	PARAM[2]=STRTOI(GetArgInc(&i));
	PARAM[0]=max(PARAM[0],0); /* n>=0 */
	PARAM[1]=max(PARAM[1],0); /* k>=0 */
	PARAM[1]=min(PARAM[1],PARAM[0]); /* k<=n */
	goto fin;
      }
    if EQUAL("udg"){ i++;
	adj=udg; POS=1;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[0]=max(PARAM[0],1); /* n>0 */
	DPARAM[0]=STRTOD(GetArgInc(&i)); /* !!! lecture d'un double */
	if(DPARAM[0]<0.0) DPARAM[0]=sqrt(log((double)PARAM[0])/(double)PARAM[0]);
	/* Threshold théorique r0 (cf. [EMY07]):
	   pour n=10,    r0=0.4798
	   pour n=100,   r0=0.2145
	   pour n=1000,  r0=0.08311
	   pour n=10000, r0=0.03034
	 */
	DPARAM[0]=CheckProba(DPARAM[0]); /* proba */
	goto fin;
      }
    if EQUAL("percolation"){ i++;
	adj=udg; POS=1; XYtype=XY_MESH; ROUND=0;
	Xmesh=STRTOI(GetArgInc(&i));
	Ymesh=STRTOI(GetArgInc(&i));
	PARAM[0]=Xmesh*Ymesh; /* nombre de sommets */
	DPARAM[0]=1.0; /* rayon */
	DELE=1.0-CheckProba(STRTOD(GetArgInc(&i))); /* proba existence arête */
	NORM=1;
	goto fin;
      }
    if EQUAL("rig"){ i++;
	adj=rig;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[0]=max(PARAM[0],1); /* n>0 */
	PARAM[1]=STRTOI(GetArgInc(&i));
	PARAM[1]=max(PARAM[1],1); /* k>0 */
	DPARAM[0]=STRTOD(GetArgInc(&i));
	if(DPARAM[0]<0.0){
	  DPARAM[0]=log((double)PARAM[0])/(double)PARAM[1];
	  if(PARAM[1]>PARAM[0]) DPARAM[0]=sqrt(DPARAM[0]/(double)PARAM[0]);
	}
	DPARAM[0]=CheckProba(DPARAM[0]); /* proba */
	goto fin;
      }
    if EQUAL("rplg"){ i++;
	adj=rplg;
	PARAM[0]=STRTOI(GetArgInc(&i));
	DPARAM[0]=STRTOD(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("thetagone"){ i++;
	adj=thetagone; POS=1;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=STRTOI(GetArgInc(&i));
	PARAM[2]=STRTOI(GetArgInc(&i));
	DPARAM[0]=CheckProba(STRTOD(GetArgInc(&i))); /* !!! lecture d'un double [0,1] */
	goto fin;
      }
    if(EQUAL("deltohedron")||EQUAL("trapezohedron")){ i++;
	adj=deltohedron;
	PARAM[0]=(STRTOI(GetArgInc(&i))<<1);
	goto fin;
      }
    if EQUAL("load"){
      param_load:
	i++;
	adj=load;
	GetArgInc(&i);i--; /* teste si le prochain argument après "load" existe bien */
	if(strlen(ARGV[i])>=(int)(sizeof(SPARAM)/sizeof(char))) Erreur(19);
	strcpy(SPARAM,ARGV[i++]); /* SPARAM=ARGV[i] */
	goto fin;
      }
    if EQUAL("bdrg"){
	adj=bdrg;
      param_sequence:
	i++;
	k=-1; /* lit tout jusqu'à une valeur <0 (incluse) */
	do{
	  k++; if(k+1==DIMAX) Erreur(4);
	  PARAM[k]=STRTOI(GetArgInc(&i));
	} while(PARAM[k]>=0);
	goto fin;
      }

    /* graphes de base avec type de paramètres déjà rencontré */
    if EQUAL("clebsch"){ adj=clebsch; goto param_hypercube; }
    if EQUAL("rpartite"){ adj=rpartite; goto param_grid; }
    if EQUAL("aqua"){ adj=aqua; DIRECTED=1; goto param_grid; }
    if EQUAL("cage"){ adj=cage; goto param_ring; }
    if EQUAL("loadc"){ LOADC=1; FORMAT=F_no; goto param_load; }
    if EQUAL("fdrg"){ adj=fdrg; goto param_sequence; }
    //
    if EQUAL("icosahedron"){ adj=icosahedron; goto param0; }
    if EQUAL("rdodecahedron"){ adj=rdodecahedron; goto param0; }
    if EQUAL("cuboctahedron"){ adj=cuboctahedron; goto param0; }
    if EQUAL("herschel"){ adj=herschel; goto param0; }
    if EQUAL("goldner-harary"){ adj=goldner_harary; goto param0; }
    if EQUAL("fritsch"){ adj=fritsch; goto param0; }
    if EQUAL("zamfirescu"){ adj=zamfirescu; goto param0; }
    if EQUAL("hatzel"){ adj=hatzel; goto param0; }
    if EQUAL("soifer"){ adj=soifer; goto param0; }
    if EQUAL("poussin"){ adj=poussin; goto param0; }
    if EQUAL("errara"){ adj=errara; goto param0; }
    if EQUAL("kittell"){ adj=kittell; goto param0; }
    if EQUAL("frucht"){ adj=frucht; goto param0; }
    if EQUAL("moser"){ adj=moser; goto param0; }
    if EQUAL("markstrom"){ adj=markstrom; goto param0; }
    if EQUAL("robertson"){ adj=robertson; goto param0; }
    if EQUAL("headwood4"){ adj=headwood4; goto param0; }
    if EQUAL("wiener-araya"){ adj=wiener_araya; goto param0; }
    if EQUAL("hgraph"){ adj=hgraph; goto param0; }
    if(EQUAL("rgraph")||EQUAL("fish")){ adj=rgraph; goto param0; }
    if EQUAL("cricket"){ adj=cricket; goto param0; }
    if EQUAL("moth"){ adj=moth; goto param0; }
    if EQUAL("cross"){ adj=cross; goto param0; }
    if EQUAL("tgraph"){ adj=tgraph; goto param0; }
    if EQUAL("bull"){ adj=bull; goto param0; }
    if EQUAL("harborth"){ adj=harborth; goto param0; }
    if EQUAL("bidiakis"){ adj=bidiakis; goto param0; }
    //
    if EQUAL("gear"){ adj=gear; goto param1; }
    if EQUAL("pstar"){ adj=pstar; PARAM[1]=2; goto param1; }
    if EQUAL("paley"){ adj=paley; goto param1; }
    if EQUAL("mycielski"){ adj=mycielski; goto param1; }
    if EQUAL("treep"){ adj=treep; goto param1; }
    if EQUAL("halin"){ adj=halin; goto param1; }
    if EQUAL("windmill"){ adj=windmill; goto param1; }
    if EQUAL("interval"){ adj=interval; goto param1; }
    if EQUAL("permutation"){ adj=permutation; goto param1; }
    if EQUAL("pancake"){ adj=pancake; goto param1; }
    if EQUAL("bpancake"){ adj=bpancake; goto param1; }
    if EQUAL("crown"){ adj=crown; goto param1; }
    if EQUAL("shuffle"){ adj=shuffle; goto param1; }
    if EQUAL("flip"){ adj=flip; goto param1; }
    if EQUAL("apollonian"){ adj=apollonian; goto param1; }
    if EQUAL("flower_snark"){ adj=flower_snark; goto param1; }
    if EQUAL("gabriel"){ adj=gabriel; POS=1; goto param1; }
    if EQUAL("rng"){ adj=rng; POS=1; goto param1; }
    if EQUAL("nng"){ adj=nng; POS=1; goto param1; }
    if EQUAL("antiprism"){ adj=antiprism; PARAM[1]=1; goto param1; }
    if EQUAL("butterfly"){ adj=butterfly; goto param1; }
    if EQUAL("matching"){ adj=matching; goto param1; }
    if EQUAL("polygon"){ adj=polygon; goto param1; }
    //
    if EQUAL("gpetersen"){ adj=gpetersen; goto param2; }
    if EQUAL("debruijn"){ adj=debruijn; goto param2; }
    if EQUAL("kautz"){ adj=kautz; goto param2; }
    if EQUAL("gpstar"){ adj=gpstar; goto param2; }
    if EQUAL("hexagon"){ adj=hexagon; goto param2; }
    if EQUAL("whexagon"){ adj=whexagon; goto param2; }
    if EQUAL("hanoi"){ adj=hanoi; goto param2; }
    if EQUAL("sierpinski"){ adj=sierpinski; goto param2; }
    if EQUAL("kpage"){ adj=kpage; goto param2; }
    if EQUAL("line-graph"){ adj=linegraph; goto param2; }
    if EQUAL("linial"){ adj=linial; goto param2; }
    if EQUAL("linialc"){ adj=linialc; goto param2; }
    if EQUAL("expander"){ adj=expander; goto param2; }
    if EQUAL("fan"){ adj=fan; goto param2; }
    //
    if EQUAL("rarytree"){ adj=rarytree; goto param3; }
    if EQUAL("barbell"){ adj=barbell; goto param3; }
    if EQUAL("planar"){ adj=planar; goto param3; }
    if EQUAL("pat"){ adj=pat; POS=1; goto param3; }
    //
    if EQUAL("chess"){ adj=chess; goto param4; }

    /********************/
    /* graphes composés */
    /********************/

    if EQUAL("theta"){ i++;
	adj=thetagone; POS=1;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=3;
	PARAM[2]=STRTOI(GetArgInc(&i));
	DPARAM[0]=6.0/(double)PARAM[2];
	goto fin;
      }
    if EQUAL("dtheta"){ i++;
	adj=thetagone; POS=1;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=3;
	PARAM[2]=STRTOI(GetArgInc(&i))/2;
	DPARAM[0]=3.0/(double)PARAM[2];
	goto fin;
      }
    if EQUAL("yao"){ i++;
	adj=thetagone; POS=1;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=0;
	PARAM[2]=STRTOI(GetArgInc(&i));
	DPARAM[0]=2.0/(double)PARAM[2];
	goto fin;
      }
    if EQUAL("path"){ i++;
	adj=grid;
	PARAM[0]=1;
	PARAM[1]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("torus"){ i++;
	adj=grid;
	PARAM[0]=2;
	PARAM[1]=-STRTOI(GetArgInc(&i));
	PARAM[2]=-STRTOI(GetArgInc(&i));
	goto fin;
    }
    if EQUAL("mesh"){ i++;
	adj=grid;
	PARAM[0]=2;
	PARAM[1]=STRTOI(GetArgInc(&i));
	PARAM[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("mobius"){ i++;
	adj=ring;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=2;
	PARAM[2]=1;
	PARAM[3]=PARAM[0]/2;
	goto fin;
      }
    if EQUAL("ladder"){ i++;
	adj=grid;
	PARAM[0]=PARAM[1]=2;
	PARAM[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("johnson"){ i++;
	adj=kneser;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=STRTOI(GetArgInc(&i));
	PARAM[2]=PARAM[1]-2;
	NOT=1-NOT;
	goto fin;
      }
    if EQUAL("star"){ i++;
	adj=rpartite;
	PARAM[0]=2;
	PARAM[1]=1;
	PARAM[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("random"){ i++;
	adj=ring;
	NOT=1-NOT;
	PARAM[0]=STRTOI(GetArgInc(&i));
	DELE=1.0-STRTOD(GetArgInc(&i));
	PARAM[1]=0;
	goto fin;
      }
    if EQUAL("tw"){
	PARAM[2]=0;
      param_tw:
	i++;
	adj=ktree;
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=STRTOI(GetArgInc(&i));
	DELE=0.5;
	goto fin;
      }
    if EQUAL("bipartite"){ i++;
	adj=rpartite;
	PARAM[0]=2;
	PARAM[1]=STRTOI(GetArgInc(&i));
	PARAM[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("cylinder"){ i++;
	adj=grid;
	PARAM[0]=2;
	PARAM[1]=STRTOI(GetArgInc(&i));
	PARAM[2]=-STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("caterpillar"){ i++;
	adj=grid;
	k=STRTOI(GetArgInc(&i)); /* nb de sommets total */
	STAR=random()%k; /* entre 0...k-1 sommets de deg=1. Active l'opération star() */
	PARAM[0]=1;
	PARAM[1]=k-STAR; /* PARAM[0]=nb de sommets du chemin, qui est >=1 */
	goto fin;
      }
    if EQUAL("centipede"){ i++;
	adj=grid;
	STAR=-1;
	PARAM[0]=1;
	PARAM[1]=STRTOI(GetArgInc(&i)); /* PARAM[0]=nb de sommets du chemin */
	goto fin;
      }
    if EQUAL("sunlet"){ i++;
	adj=grid;
	PARAM[0]=1;
	PARAM[1]=-STRTOI(GetArgInc(&i));
	STAR=-1;
	goto fin;
      }
    if EQUAL("sunflower"){ i++;
	adj=cage;
	PARAM[0]=2*STRTOI(GetArgInc(&i));
	PARAM[1]=PARAM[2]=2;
	PARAM[3]=0;
	goto fin;
      }
    if EQUAL("wheel"){ i++;
	adj=ringarytree;
	PARAM[0]=1;PARAM[1]=0;
	PARAM[2]=STRTOI(GetArgInc(&i));
	PARAM[3]=2;
	goto fin;
      }
    if EQUAL("tadpole"){ i++;
	adj=barbell;
	PARAM[0]=-STRTOI(GetArgInc(&i));
	PARAM[1]=1;
	PARAM[2]=STRTOI(GetArgInc(&i));
	goto fin;
      }
    if EQUAL("pan"){ i++;
	adj=barbell;
	PARAM[0]=-STRTOI(GetArgInc(&i));
	PARAM[1]=1;
	PARAM[2]=1;
	goto fin;
      }
    if EQUAL("web"){ i++;
	adj=ringarytree;
	PARAM[2]=STRTOI(GetArgInc(&i));
	PARAM[0]=STRTOI(GetArgInc(&i));
	PARAM[1]=1;
	PARAM[3]=2;
	goto fin;
      }
    if EQUAL("plrg"){ i++;
	adj=bdrg;
	PARAM[0]=STRTOI(GetArgInc(&i));
	DPARAM[0]=STRTOD(GetArgInc(&i));
	int *S=power_law_seq(PARAM[0],DPARAM[0],NULL);
	if(S==NULL) Erreur(6);
	k=-1;
	do{ /* copie S dans PARAM, y compris le -1 */
	  k++; if(k+1==DIMAX) Erreur(4);
	  PARAM[k]=S[k];
	} while(PARAM[k]>=0);
	goto fin;
      }

    /* graphes composés avec type de paramètres déjà rencontré */
    if EQUAL("kpath"){ PARAM[2]=1; goto param_ktree; }
    if EQUAL("pw"){ PARAM[2]=1;	goto param_tw; }
    //
    if EQUAL("tree"){ adj=arboricity; PARAM[1]=1; goto param1; }
    if EQUAL("cycle"){ adj=ring; PARAM[1]=PARAM[2]=1; goto param1; }
    if EQUAL("binary"){ adj=ringarytree; PARAM[1]=PARAM[2]=2;PARAM[3]=0; goto param1; }
    if EQUAL("rbinary"){ adj=rarytree; PARAM[1]=2; PARAM[2]=0; goto param1; }
    if EQUAL("rbinaryz"){ adj=rarytree; PARAM[1]=2; PARAM[2]=1; goto param1; }
    if (EQUAL("stable")||EQUAL("empty")){ adj=ring; PARAM[1]=0; goto param1; }
    if EQUAL("clique"){ adj=ring; NOT=1-NOT; PARAM[1]=0; goto param1; }
    if EQUAL("outerplanar"){ adj=kpage; PARAM[1]=1; goto param1; }
    if EQUAL("squaregraph"){ adj=planar; PARAM[1]=PARAM[2]=4; goto param1; }
    if EQUAL("prism"){ adj=gpetersen; PARAM[1]=1; goto param1; }
    if EQUAL("cubic"){ adj=fdrg; PARAM[1]=3;PARAM[2]=-1; goto param1; }
    if EQUAL("d-octahedron"){ adj=matching; NOT=1-NOT; goto param1; }
    if EQUAL("point"){ adj=ring; PARAM[1]=0; POS=1; goto param1; }
    if EQUAL("td-delaunay"){ adj=thetagone;POS=1;PARAM[1]=PARAM[2]=3;
	DPARAM[0]=1.0;goto param1; }
    //
    if EQUAL("lollipop"){ adj=barbell; PARAM[2]=0; goto param2; }
    if EQUAL("knight"){ adj=chess; PARAM[2]=1;PARAM[3]=2; goto param2; }
    if EQUAL("camel"){ adj=chess; PARAM[2]=1;PARAM[3]=3; goto param2; }
    if EQUAL("giraffe"){ adj=chess; PARAM[2]=1;PARAM[3]=4; goto param2; }
    if EQUAL("zebra"){ adj=chess; PARAM[2]=2;PARAM[3]=3; goto param2; }
    if EQUAL("antelope"){ adj=chess; PARAM[2]=2;PARAM[3]=4; goto param2; }
    if EQUAL("regular"){ adj=fdrg; PARAM[2]=-1; goto param2; }
    //
    if EQUAL("arytree"){ adj=ringarytree; PARAM[3]=0; goto param3; }

    /* graphes sans paramètres mais composés. Le graphe de base
       contient des paramètres, on doit passer par CheckHelp() */

    // 1 paramètre
    if(EQUAL("cube")||EQUAL("hexahedron")){ adj=crown; PARAM[0]=4; goto param0; }
    if EQUAL("octahedron"){ adj=antiprism; PARAM[0]=3; goto param0; }
    if EQUAL("dodecahedron"){ adj=gpetersen; PARAM[0]=10; PARAM[1]=2; goto param0; }
    if EQUAL("associahedron"){ adj=flip; PARAM[0]=6; goto param0; }
    if EQUAL("tietze"){ adj=flower_snark; PARAM[0]=3; goto param0; }
    if EQUAL("grotzsch"){ adj=mycielski; PARAM[0]=4; goto param0; }
    // 2 paramètres
    if EQUAL("nauru"){ adj=pstar; PARAM[0]=4;PARAM[1]=2; goto param0; }
    if EQUAL("hajos"){ adj=sierpinski; PARAM[0]=2;PARAM[1]=3; goto param0; }
    if EQUAL("house"){ adj=grid; NOT=1-NOT; PARAM[0]=1;PARAM[1]=5; goto param0; }
    if EQUAL("tetrahedron"){ adj=ring; NOT=1-NOT; PARAM[0]=4; PARAM[1]=0; goto param0; }
    if EQUAL("claw"){ adj=rpartite; PARAM[0]=2;PARAM[1]=1;PARAM[2]=3; goto param0; }
    if EQUAL("desargues"){ adj=gpetersen; PARAM[0]=10;PARAM[1]=3; goto param0; }
    if EQUAL("durer"){ adj=gpetersen; PARAM[0]=6;PARAM[1]=2; goto param0; }
    if EQUAL("netgraph"){ adj=grid; STAR=-1; PARAM[0]=1;PARAM[1]=-3; goto param0; }
    if EQUAL("gem"){ adj=fan; PARAM[0]=4;PARAM[1]=1; goto param0; }
    if EQUAL("diamond"){ adj=fan; PARAM[0]=PARAM[1]=2; goto param0; }
    if EQUAL("egraph"){ adj=grid; PARAM[0]=1;PARAM[1]=3;STAR=-1; goto param0; }
    // 3 paramètres
    if EQUAL("petersen"){ adj=kneser; PARAM[0]=5;PARAM[1]=2;PARAM[2]=0; goto param0; }
    if EQUAL("banner"){ adj=barbell; PARAM[0]=-4;PARAM[1]=1;PARAM[2]=1; goto param0; }
    if EQUAL("paw"){ adj=barbell; PARAM[0]=-3;PARAM[1]=1;PARAM[2]=1; goto param0; }
    if EQUAL("theta0"){ adj=barbell; PARAM[0]=PARAM[1]=-5;PARAM[2]=-2; goto param0; }
    if EQUAL("utility"){ adj=rpartite; PARAM[0]=2;PARAM[1]=PARAM[2]=3; goto param0; }
    // 4 paramètres
    if EQUAL("wagner"){ adj=ring;
	PARAM[0]=8; PARAM[1]=2; PARAM[2]=1; PARAM[3]=4; goto param0; }
    if EQUAL("headwood"){ adj=cage;
	PARAM[0]=14; PARAM[1]=2; PARAM[2]=5; PARAM[3]=-5; goto param0; }
    if EQUAL("ygraph"){ adj=ringarytree;
	PARAM[0]=2; PARAM[1]=1; PARAM[2]=3; PARAM[3]=0; goto param0; }
    if EQUAL("franklin"){ adj=cage;
	PARAM[0]=12; PARAM[1]=2; PARAM[2]=5; PARAM[3]=-5; goto param0; }
    // 5 paramètres
    if EQUAL("mcgee"){ adj=cage;
	PARAM[0]=24; PARAM[1]=3; PARAM[2]=12; PARAM[3]=7; PARAM[4]=-7;
	goto param0;
      }
    // 6 paramètres
    if EQUAL("gosset"){ adj=ggosset;
	PARAM[0]=8; PARAM[1]=2;	PARAM[2]=2; PARAM[3]=3;
	PARAM[4]=6; PARAM[5]=-1;
	goto param0;
      }
    if EQUAL("dyck"){ adj=cage;
	PARAM[0]=32; PARAM[1]=4;
	PARAM[2]=5; PARAM[3]=0; PARAM[4]=13; PARAM[5]=-13;
	goto param0;
      }
    // 8 paramètres
    if EQUAL("pappus"){ adj=cage;
	PARAM[0]=18; PARAM[1]=6; PARAM[2]=5; PARAM[3]=7; PARAM[4]=-7;
	PARAM[5]=7; PARAM[6]=-7; PARAM[7]=-5;
	goto param0;
      }
    if EQUAL("tutte-coexter"){ adj=cage;
	PARAM[0]=30; PARAM[1]=6;
	PARAM[2]=-7;PARAM[3]=9;PARAM[4]=13;
	PARAM[5]=-13;PARAM[6]=-9;PARAM[7]=7;
	goto param0;
      }
    if EQUAL("gray"){ adj=cage;
	PARAM[0]=54; PARAM[1]=6; PARAM[2]=7; PARAM[3]=-7;
	PARAM[4]=25;PARAM[5]=-25; PARAM[6]=13; PARAM[7]=-13;
	goto param0;
      }
    // 12 paramètres
    if EQUAL("chvatal"){ adj=cage;
	PARAM[0]=PARAM[1]=12;
	PARAM[2]=PARAM[4]=PARAM[7]=PARAM[10]=PARAM[12]=PARAM[13]=3;
	PARAM[3]=PARAM[5]=PARAM[6]=PARAM[8]=6;
	PARAM[9]=PARAM[11]=-3;
	goto param0;
      }

  fin:
    if(j==i){
      if PREFIX("-") Erreur(2); /* option non trouvée */
      Erreur(10); /* graphe non trouvé */
    }
    
  } /* fin du while(i<ARGC) ... */

  /* options qui ne vont pas ensemble */
  if((STAR)&&(APEX)) Erreur(18); /* on ne peut pas avoir les deux */
  if((LOADC)&&(PERMUTE||NOT)) Erreur(29); /* options incompatibles */
  if((LOADC)&&(CHECK<2)) Erreur(30); /* manque -check */


  /***********************************

           COEUR DU GENERATEUR

  ***********************************/

  if(STAR) { ADJ0=adj; adj=star; } /* le graphe est maintenant star() */
  if(APEX) { ADJ0=adj; adj=apex; } /* le graphe est maintenant apex() */

  /* initialisation du graphe, calcule N avant la suppression
     éventuelle de sommet, lit le graphe LOAD si adj=load, détermine
     (XPOS,YPOS) si graphe géométrique */  

  adj(ADJ_INIT,0); // if(adj(ADJ_INIT,0)==ADJ_ERR) Erreur(6);

  if(N<0) N=0; /* au cas où */
  if(LOADC){ GF=LOAD; NF=N; adj=NULL; goto check; } /* saute la partie génération d'arêtes */
  if((POS)&&(XPOS==NULL)) InitXY(); /* il faut déterminer les positions */
  if(LABEL==1) PERMUTE=0; /* si on souhaite les labels d'origine, on ne permute rien */

  ALLOC(V,N);       /* V[i]=étiquette du sommet i, -1 si i est supprimé */
  ALLOC(VF,N);      /* VF[j]=indice du j-ème sommet non supprimé */
  NF=InitVertex(N); /* initialise V[i], VF[i] et retourne NF=#sommets final */
  ALLOC(INC,N);     /* INC[i]=1 ssi i possède un voisin, 0 si sommet isolé */
  for(j=0;j<NF;j++) INC[VF[j]]=0;
  
  /* constantes pour accélérer les tests de la boucle principale */
  long seuil_edge=(double)(1.0-DELE)*(double)RAND_MAX;
  long seuil_redirect=(double)(REDIRECT)*(double)RAND_MAX;

  /*
    Génère les adjacences i-j en tenant compte des sommets isolés et
    des sommets supprimés. Les sommets isolés sont affichés en dernier
    à cause de l'option -redirect. On a toujours i<j lorsque l'arête
    i-j doit être sortie. Si on a FAST=1, alors on génère le graphe à
    partir de LOAD, s'il existe.
  */

  Out(-1,-1); /* initialise le format d'affichage */

  if(FAST){
    if(LOAD){
      int d,t,redirect=(REDIRECT!=0.0);

      /* on ne teste que les arêtes de LOAD */
      for(i=0;i<N;i++)
	if(V[i]>=0) /* si le sommet i existe */
	  for(t=0,d=LOAD->d[i];t<d;t++){
	    j=LOAD->L[i][t]; if((LOAD->sym)&&(j<i)) continue; /* il faut i<j si LOAD symétrique */
	    if((V[j]>=0)&&((!DIRECTED)||(i!=j)||(LOOP))) /* si j existe ... */
	      if(random()<seuil_edge){
		if(redirect){ /* si redirection */
		  j=(random()<seuil_redirect)? random()%N : j;
		  if((V[j]<0)||(j==i)) continue; /* prochain voisin */
		}
		INC[i]++; /* un voisin de plus pour i */
		INC[j]++; /* un voisin de plus pour j */
		Out(i,j); /* sort l'arête i-j avec i<j */
	      }
	  }
    }else Erreur(25); /* -fast mais LOAD n'existe pas */
  }else{
    int noredirect=(REDIRECT==0.0);

    /* teste les O(n^2) arcs ou arêtes possibles */
    for(i=0;i<N;i++)       /* pour tous les i */
      if(V[i]>=0){         /* si i existe */
	for(j=(DIRECTED)?0:i+1-LOOP;j<N;j++) /* pour tous les j>i */
	  if((V[j]>=0)&&((!DIRECTED)||(i!=j)||(LOOP))) /* si j existe ... */
	    if((random()<seuil_edge)&&(adj(i,j)^NOT)){
	      /* ici l'arête i-j devrait être sortie */
	      if(noredirect){ /* si pas de redirection d'arête */
		INC[i]++; /* un voisin de plus pour i */
		INC[j]++; /* un voisin de plus pour j */
		Out(i,j); /* sort l'arête i-j */
	      }else{ /* on redirige l'arête i-j vers i-k */
		k=(random()<seuil_redirect)? random()%N : j;
		if((V[k]>=0)&&(k!=i)){
		  /* on affiche l'arête que si k existe et si k<>i */
		  /* Attention ! il ne faut toucher ni à i ni à j */
		  INC[i]++; /* un voisin de plus pour i */
		  INC[k]++; /* un voisin de plus pour k */
		  if(k<i) Out(k,i); else Out(i,k); /* pour avoir i<j */
		}
	      }
	    }
      }

  }

  /* affiche les sommets isolés */
  for(i=0;i<N;i++)
    if((V[i]>=0)&&(!INC[i])) Out(i,-1);
  
  /* fin de l'affichage, doit être fait avant adj(0,ADJ_END) */
  Out(-1,0); /* calcule GF (si CHECK) */

  if((CHECK)&&(POS)){ /* mémorise XPOS,YPOS qui vont être supprimés par adj(0,ADJ_END) */
    ALLOCZ(GF->xpos,NF,XPOS[VF[_i]]);
    ALLOCZ(GF->ypos,NF,YPOS[VF[_i]]);
  }

  free(V);
  free(VF);
  free(INC);
  adj(0,ADJ_END); /* libère éventuellement la mémoire utilisée pour adj() */
  adj=NULL; /* ici adj() n'est plus définie */

  /***********************************

           FIN DU GENERATEUR

           Si on a CHECK<>0 alors:
           - le graphe généré est GF
           - son nombre de sommets est NF

  ***********************************/

 check:
  /* NB: dans le cas LOADC on ne fait pas adj(0,ADJ_END), c'est-à-dire
     load(0,ADJ_END). En effet, load(0,ADJ_END) libère le graphe
     LOAD. Or on veut GF=LOAD précisément. */

  if(CHECK){
    if((GF==NULL)||(GF->n!=NF)) Erreur(32); /* ne devrait jamais arriver */
    if(CHECK==CHECK_MAINCC){
      param_dfs *p=dfs(GF,MEM(CPARAM,0,int),NULL);
      int d0,d1,n,c;
      for(i=c=n=d0=d1=0;i<p->nc;i++){ /* détermine la plus grosse cc */
	/* d1-d0=nb de sommets de la cc numéro i */
	if(i+1<p->nc) d1=p->d[p->R[i+1]]; else d1=NF;
	if(d1-d0>n){ n=d1-d0; c=i; } /* n=nb de sommets de cc max */ 
	d0=d1;
      }
      c=p->C[p->R[c]]; /* c=couleur de la composante max */
      NALLOC(int,T,n); /* T=sommet de GF à garder */
      for(i=j=0;i<NF;i++) if(p->C[i]==c) T[j++]=i;
      free_param_dfs(p); /* p ne sert plus à rien */
      graph *C=ExtractSubgraph(GF,T,n,1); /* construit la cc max */
      free(T); /* T ne sert plus à rien */
      CHECK=CHECK_OFF; /* pour ne pas restocker le graphe */
      PrintGraph(C);
      free_graph(C);
    }
    if(CHECK==CHECK_BFS){
      param_bfs *p=bfs(GF,MEM(CPARAM,0,int),NULL);
      int t,c;
      printf("root=%i\n",p->root);
      printf("rad[%i]=%i\n",p->root,p->radius);
      printf("cycle[%i]=%i\n",p->root,p->cycle);
      printf("#vertices traversed=%i\n",p->n);
      printf("distance:\n");
      for(k=1;k<=p->radius;k++){
	printf(" %i:",k);
	for(c=t=0;t<NF;t++) if(p->D[t]==k){ c++; printf(" %i",t); }
	printf("  (%i)\n",c);
      }
      printf("vertices not connected to the source:");
      for(t=c=0;t<NF;t++) if(p->D[t]<0){ c++; printf(" %i",t); }
      if(c==0) printf(" none\n");
      else printf("  (%i)\n",c);
      printf("parent:");
      for(k=0;k<NF;k++) printf(" %i",p->P[k]);
      printf("\n");
      free_param_bfs(p);
    }
    if((CHECK==CHECK_DFS)||(CHECK==CHECK_NCC)){
      if(NF<=0) printf("Empty graph\n");
      int s=(CHECK==CHECK_DFS)? MEM(CPARAM,0,int) : 0;
      if(CHECK==CHECK_DFS) printf("source: %i\n",s);
      param_dfs *p=dfs(GF,s,NULL);
      printf("#component: %i%s\n",p->nc,(p->nc==1)?" (connected)":"");
      printf("#cut-vertex: %i%s\n",p->na,
	     ((p->nc==1)&&(p->na==0)&&(NF>2))?" (biconnected)":"");
      if(CHECK==CHECK_NCC) goto check_ncc;
      printf("root:");
      for(i=0;i<p->nc;i++) printf(" %i",p->R[i]);
      if(p->na) printf("\ncut-vertex:");
      for(i=0;i<NF;i++) if(p->A[i]) printf(" %i",i);
      if(p->nc>1){
	printf("\ncomponent:");
	for(i=0;i<NF;i++) printf(" %i",p->C[i]);
      }
      printf("\nparent:");
      for(i=0;i<NF;i++)
	if(p->P[i]<0) printf(" -"); else printf(" %i",p->P[i]);
      printf("\n");     
    check_ncc:
      free_param_dfs(p);
    }
    if(CHECK==CHECK_BELLMAN){
      if(!WeightGraph(GF))
	printf("Seulement possible pour les graphes géométriques !\n");
      else{
	double dx,dy,s=-1.0;
	param_bellman *p=Bellman_Ford(GF,MEM(CPARAM,0,int));
	int u=p->source;
	printf("source=%i\n",u);
	for(k=0;k<NF;k++){
	  dx=GF->xpos[p->source]-GF->xpos[k];
	  dy=GF->ypos[p->source]-GF->ypos[k];
	  dx=sqrt(dx*dx+dy*dy);
	  dy=(p->source!=k)? p->dist[k]/dx : 0.0;
	  if(dy>s){ s=dy; u=k; }
	  printf("dist[%i]=%lf parent=%i stretch=%lf (%lf)\n",k,p->dist[k],p->pere[k],dy,dx);
	}
	printf("stretch max: %lf\n",s);
	while(u!=p->source){
	  printf("%i->",u);
	  u=p->pere[u];
	}
	printf("%i\n",u);
	free_param_bellman(p);
      }
    }
    if(CHECK==CHECK_DEG){
      int *F;
      printf("#edges: %i\n",NbEdges(GF));
      printf("degmin: %i degmax: %i\n",Degree(GF,0),Degree(GF,1));
      k=NF;
      F=SortInt(GF->d,NULL,NF,0,&k,SORT_FREQv);
      printf("deg:");
      for(k=0;k<NF;k++) if(F[k]) printf(" %i (x%i) ",k,F[k]);
      printf("\n");
      free(F);
    }
    if(CHECK==CHECK_DEGENERATE){
      int *T=Prune(GF,&k);
      printf("Degenerate: %i\n",k);
      for(k=0;k<NF;k++) printf("%i ",T[k]);
      printf("\n");
      free(T);
    }
    if(CHECK==CHECK_GCOLOR){
      int *T=Prune(GF,NULL);
      int *C=GreedyColor(GF,T);
      printf("#colors: %i\n",1+GF->int1);
      PrintMorphism("Coloring (node->color):\n",C,GF->n);
      free(C);
      free(T);
    }
    if(CHECK==CHECK_KCOLOR){
      int k=MEM(CPARAM,0,int);
      int *C=kColor(GF,k);
      if(C==NULL) printf("There is no %i-coloration for this graph.\n",k);
      else{
	printf("#colors: %i\n",1+GF->int1);
	PrintMorphism("Coloring (node->color):\n",C,GF->n);
	free(C);
      }
    }
    if(CHECK==CHECK_KCOLORSAT){
      kColorSat(GF,MEM(CPARAM,0,int));
    }
    if(CHECK==CHECK_KINDEPSAT){
      kIndepSat(GF,MEM(CPARAM,0,int));
    }
    if((CHECK==CHECK_PS1) ||(CHECK==CHECK_PS1b)||
       (CHECK==CHECK_PS1c)||(CHECK==CHECK_PS1x)){
      int c=0;
      if(CHECK==CHECK_PS1b) c=1;
      if(CHECK==CHECK_PS1c) c=2;
      if(CHECK==CHECK_PS1x) c=3;
      path *P=new_path(GF,NULL,NF);
      int v=PS1(GF,P,c);
      printf("#tests: %i\nPS1: %s\n",GF->int1,v?"yes (PS1 for sure)":"no (probably not PS1)");
      free_path(P);
    }
    if(CHECK==CHECK_TWDEG){
      int *T=Prune(GF,&k);
      printf("treewidth <= %i\n",Treewidth(GF,0));
      printf("treewidth >= %i\n",k);
      free(T);
    }
    if(CHECK==CHECK_TW){
      k=Treewidth(GF,1);
      printf("#tests: %i\ntreewidth: %i\n",GF->int1,k);
    }
    if(CHECK==CHECK_DIAMETER){
      param_bfs *p=new_param_bfs();
      int u,x=-1;
      p->clean=1;
      for(u=0;u<NF;u++){
	p=bfs(GF,u,NULL);
	if(p->n<NF) break; /* non connexe */
	x=max(x,p->radius);
      }
      free_param_bfs(p);
      printf("diameter: ");
      if(x<0) printf("%s","+∞");
      else printf("%i",x);
      printf("\n");
    }
    if(CHECK==CHECK_RADIUS){
      param_bfs *p=new_param_bfs();
      int u,x=NF-1;
      p->clean=1;
      for(u=0;u<NF;u++){
	bfs(GF,u,p);
	if(p->n<NF){ x=-1; break; } /* non connexe */
	x=min(x,p->radius);
      }
      free_param_bfs(p);
      printf("radius: ");
      if(x<0) printf("%s","+∞");
      else printf("%i",x);
      printf("\n");
    }
    if(CHECK==CHECK_GIRTH){
      param_bfs *p=new_param_bfs();
      int u,x=1+NF;
      p->clean=1;
      for(u=0;u<NF;u++){
	bfs(GF,u,p);
	if(p->cycle>0) x=min(x,p->cycle);
      }
      free_param_bfs(p);
      if(x>NF) x=-1;
      printf("girth: %i\n",x);
    }
    if((CHECK==CHECK_PATHS)||(CHECK==CHECK_PATHS2)){/*  Sort tous les chemins de x à y */
      int (*f)(graph*,path*,int)=(CHECK==CHECK_PATHS)? NextPath : NextPath2;
      path *P=new_path(GF,NULL,NF); /* chemin vide */
      P->P[0]=MEM(CPARAM,0,int);
      P->P[1]=MEM(CPARAM,sizeof(int),int);
      int u,v=f(GF,P,-1); /* initialise le premier chemin */
      while(v){
	for(u=0;u<P->n;u++) printf("%i ",P->P[u]);
	printf("\n");
	v=f(GF,P,0);
      }
      free_path(P);
    }
    if(CHECK==CHECK_ISO){
      char *s=ARGV[MEM(CPARAM,0,int)];
      graph *H=File2Graph(s,34);
      int *P=Isomorphism(GF,H);
      printf("H: %s\n#tests: %i\n",s,H->int1);
      if(P==NULL) printf("Non-isomorphic.\n");
      else{
	PrintMorphism("Isomorphism G->H:\n",P,NF);
	free(P);
      }
      free_graph(H);
    }
    if(CHECK==CHECK_SUB){
      char *s=ARGV[MEM(CPARAM,0,int)];
      graph *H=File2Graph(s,34);
      graph *S=Subgraph(GF,H);
      printf("H: %s\n#tests: %i\n",s,H->int1);
      if(S==NULL) printf("G is not a subgraph of H.\n");
      else{
	printf("Subgraph S of H isomorphic to G:\n");
	PrintGraph(S);
	PrintMorphism("Isomorphism S->G:\n",S->pint1,S->n);
	free_graph(S);
      }
      free_graph(H);
    }
    if(CHECK==CHECK_MINOR){
      char *s=ARGV[MEM(CPARAM,0,int)];
      graph *H=File2Graph(s,34);
      int *C=Minor(H,GF);
      printf("H: %s\n#tests: %i\n",s,H->int1);
      if(C==NULL) printf("H is not a minor of G.\n");
      else{
	int c,u;
	printf("Model of H in G:\n");
	for(c=0;c<H->n;c++){ /* pour chaque couleur c */
	  printf("%i -> {",c);
	  for(u=0;u<NF;u++) /* on affiche les sommets de la couleur c */
	    if(C[u]==c) printf(" %i",u);
	  printf(" }\n");
	}
	free(C);
      }
      free_graph(H);
    }
    if(CHECK==CHECK_ISUB){
      char *s=ARGV[MEM(CPARAM,0,int)];
      graph *H=File2Graph(s,34);
      int *X=InducedSubgraph(H,GF);
      printf("H: %s\n#tests: %i\n",s,GF->int1);
      if(X==NULL) printf("H is not an induced subgraph of G.\n");
      else{
	int u;
	printf("Vertices of the induced subgraph S:");
	for(u=0;u<H->n;u++) printf(" %i",X[u]);
	for(u=0;u<H->n;u++) GF->pint1[u]=X[u];
	PrintMorphism("\nIsomorphism H->S:\n",GF->pint1,H->n);
	free(X);
      }
      free_graph(H);
    }
    if(CHECK==CHECK_INFO){
      printf("- command: %s\n",MakeCMD(NULL,0,ARGC));
      printf("- total time: %s\n",TopChrono(0));
      printf("- seed: %u\n",SEED);
      int *R=SortGraph(GF,1);
      if(R){
	printf("- simple and undirected: %s\n",R[6]?"yes":"no");
	printf("- #nodes: %s\n",millier(NF));
	printf("- #arcs: %s\n",millier(R[2]));
	printf("- #self-loops: %s\n",millier(R[0]));
	printf("- #multi-arcs: %s\n",millier(R[1]));
	printf("- #asym. adjacency: %s\n",millier(R[3])); 
	printf("- #node IDs < 0: %s\n",millier(R[4])); 
	printf("- #node IDs >=n: %s\n",millier(R[5])); 
	printf("- maximum degree: %s\n",millier(R[7]));
	printf("- minimum degree: %s\n",millier(R[8]));
	printf("- #isolated nodes: %s\n",millier(R[9]));
      }else printf("- empty graph\n");
    }
    if(CHECK==CHECK_RS_CLUSTER){
      RS_Start("cluster",RS_NI_FP,GF);

      /* paramètre */
      int k=MEM(CPARAM,0,int); /* k=paramètre */
      if(k==-1) k=ceil(sqrt((double)NF)); /* ici k>=1, toujours */
      if(k==-2) k=NF;
      printf("- parameter: %i\n",k); /* ici k>=1, toujours */
      if(k<=0) Erreur(6); /* k=0: valeur impossible */
      if(VARIANT) printf("- variant: %i\n",VARIANT);
      BARRE;
      
      /* construit, teste, puis libère les tables */
      rs_cluster_tables *RT=rs_cluster(GF,k); /* construit */
      routing_test(GF,RT,(rt_length)rs_cluster_length,-1,NULL); /* teste */
      free_rs_cluster_tables(RT); /* libère */
    }
    if(CHECK==CHECK_RS_DCR){
      RS_Start("dcr",RS_NI_FP,GF);

      /* paramètre */
      int k=MEM(CPARAM,0,int); /* k=paramètre */
      if(k==-1) k=ceil(Minimize(func1,&NF,1,NF,0)); /* ici k>=1 */
      if(k==-2) k=NF;
      printf("- parameter: %i\n",k); /* ici k>=1, toujours */
      if(k<1) Erreur(6); /* il faut k>0 */
      BARRE;

      /* construit, teste, puis libère les tables */
      rs_dcr_tables *RT=rs_dcr(GF,k);
      routing_test(GF,RT,(rt_length)rs_dcr_length,-1,RT->dist);
      free_rs_dcr_tables(RT);
    }
    if(CHECK==CHECK_RS_TZRPLG){
      RS_Start("tz rplg",RS_L_FP,GF);

      /* paramètre */
      double t=MEM(CPARAM,0,double); /* t=paramètre (exposant du RPLG) */
      if((VARIANT<2)&&(0.0<t)&&(t<2.0)) Erreur(6); /* valeurs impossibles */
      if((VARIANT==2)&&(t<1.0)) Erreur(6); /* valeurs impossibles */
      printf("- parameter: %g\n",t);
      BARRE;

      /* construit, teste, puis libère les tables */
      rs_tzrplg_tables *RT=rs_tzrplg(GF,t);
      routing_test(GF,RT,(rt_length)rs_tzrplg_length,-1,NULL); /* routage */
      free_rs_tzrplg_tables(RT); /* libère les tables */
    }
    if(CHECK==CHECK_RS_BC){
      RS_Start("bc",RS_L_FP,GF);

      /* paramètre */
      int k=MEM(CPARAM,0,int); /* k=paramètre */
      printf("- parameter: %i\n",k); /* ici k>=0, toujours */
      if(k<0) Erreur(6); /* k=0: valeur impossible */
      BARRE;
      
      /* construit, teste, puis libère les tables */
      rs_bc_tables *RT=rs_bc(GF,k); /* construit */
      routing_test(GF,RT,(rt_length)rs_bc_length,-1,NULL); /* teste */
      free_rs_bc_tables(RT); /* libère */
    }
    
    free_graph(GF),GF=NULL; /* supprime le graphe */
    free(CPARAM),CPARAM=NULL; /* supprime les paramètres pour CHECK */
  }
  /* fin du "if(CHECK)" */

  TopChrono(-1); /* libère les chronos */
  return 0; /* fin de gengraph */
}

/*# ###
Générateur de graphes - v4.0 - © Cyril Gavoille - March 2016

USAGE

       gengraph [-option] graph_name [parameter_list]


DESCRIPTION

       Génére sur la sortie standard un graphe. Par défaut le graphe
       est non orienté et affiché selon une liste d'arêtes (en texte),
       mais d'autres formats sont possibles: liste d'adjacence, format
       dot de GraphViz, xfig ou pdf par exemple. En paramètre figure
       le nom du graphe ainsi que ses paramètres éventuels, comme le
       nombre de sommets. La commande appelée seule affiche l'aide sur
       les options du générateur. Si les paramètres d'une option ou
       d'un graphe sont absents ou remplacés par "?", une aide
       spécfique est affichée. Une console supportant l'UTF8 est
       préférable.

       Ex:
          gengraph -help
	  gengraph -list | sort
	  gengraph tree ?
	  gengraph ? arbre
          gengraph tutte
          gengraph hypercube 8
          gengraph mesh 7 3 -not
	  gengraph mesh 50 50 -dele .5 -maincc -visu
	  gengraph rdodecahedron -visu
          gengraph tree 100 -visu
	  gengraph web 10 3 -visu
	  gengraph gabriel 50 -caption "Grabriel with n=50" -visu
	  gengraph gabriel 2000 -xy seed 1 0.15 -visu
	  gengraph sierpinski 7 3 -visu
	  gengraph udg 400 .1 -visu
	  gengraph udg 400 .1 -xy seed 3 1.5 -visu
	  gengraph udg 400 -1 -vsize -vcolor deg -visu
          gengraph arytree 6 3 3 -dotfilter circo -visu
          gengraph dyck -dotfilter circo -visu
	  gengraph ringarytree 4 2 3 0 -label 1 -visu
	  gengraph arboricity 100 2 -vcolor degr -visu
	  gengraph prime 6 -directed -noloop -visu
	  gengraph aqua 3 3 2 1 -label 1 -directed -dotfilter dot -visu
	  echo "0->1->2->0" | gengraph load - -check bfs 0
	  gengraph tutte | gengraph -filter - diameter p
          gengraph rplg 300 3 -maincc -vcolor degr -vcolor pal wz -vsize -visu
          gengraph -xy box 15 15 -xy round 0 -xy grid 16 rng 30 -visu
	  gengraph linial 7 3 -check kcolorsat 3 | ./glucose -model


   LE FORMAT STANDARD

       Le format par défaut (ou standard) est une liste d'arêtes ou de
       chemins écrits en texte simple. Ce format minimaliste est très
       proche de celui du format "dot" de GraphViz.  D'autres formats
       de sortie sont possibles, notamment le format "dot" (voir
       l'option -format).  Les sommets sont numérotés consécutivement
       de 0 à n-1 où n est le nombre de sommets présents dans le
       graphe (en fait cela peut être changé avec l'option -shift).
       Une arête entre i et j est représentée par i-j, un arc de i
       vers j par i->j. Les sommets isolés sont simplement représentés
       par le numéro du sommet suivit d'un espace ou d'un retour de
       ligne. Le nombre de sommets du graphe est l'entier le plus
       grand + 1, et s'il y a i-j (ou i->j), alors il existe une arête
       (ou un arc) entre les sommets i et j.

       Pour une représentation plus compacte, les arêtes (ou arcs)
       consécutives d'un chemin du graphe peuvent être regroupées en
       blocs i-j-k-.... Par exemple, les deux arêtes 3-5 et 5-8
       peuvent être regroupées en 3-5-8. Mais ce n'est pas
       obligatoire.  Également, les arêtes (ou arcs) d'une étoile
       peuvent être groupées avec i-(j k ...).  Par exemple, 3-(5 7 8)
       représente les arêtes 3-5, 3-7 et 3-8.  Il n'est pas possible
       cependant de combiner chemins et étoiles, comme 3-(5-7-8) ou
       3-(5-(7 8)). Toutefois 3-5-(7 ...) est correct, mais pas 3-(5
       6)-7 ni (3 5)-6. Les sommets isolés et les arêtes (ou les blocs
       d'arêtes) sont séparés par des espaces ou des sauts de ligne.
       Une boucle sur un sommet i est codée par i-i. Les arêtes
       multiples sont codées par la répétition d'une même arête, comme
       par exemple i-j i-j, ou encore i-j-i (même convention pour les
       arcs i->j->i).

       Quelques exemples:       0   1
                                   / \
       Ex1: 0 1-2-3-1             3---2

       représente un graphe à 4 sommets, composé d'un sommet isolé (0)
       et d'un cycle à trois sommets (1,2,3). Une représentation
       graphique possible est donnée à droite.

       Ex2: 4-2-1-0-3-2-5

       représente un graphe à 6 sommets composé d'un cycle de longueur
       4 et de deux sommets de degré 1 attaché à 2. On aurait pu coder
       le même graphe avec l'expression 2-(1 3 4 5) 1-0-3. En voici
       une représentation graphique possible:

                  1
                 / \
              4-2   0
               / \ /
              5   3

       Plus généralement, une famille de graphes peuvent être définie
       en précédant chaque graphe par "[n]" où n est un entier unique
       représentant l'indentifiant du graphe.

       Ex3: [17] 0-1 [22] 0->1->2->0

       représente une famille composée de deux graphes, un chemin à
       deux sommets ainsi qu'un cycle orienté à trois sommets.

   COMMENT FONCTIONNE LE GENERATEUR ?

       Pour chaque graphe une fonction adj(i,j) est définie. Elle
       fournit l'adjacence (0 ou 1) entre les sommets i et j, des
       entiers entre 0 et n-1. Le graphe est affiché en générant
       toutes les paires {i,j} possibles et en appelant adj(i,j) (ou
       tous les couples (i,j) possibles dans le cas orienté). Les
       graphes sont ainsi générés de manière implicite. Les arêtes du
       graphe ne sont pas stockées en mémoire, mais affichées à la
       volée.  Ceci permet de générer des graphes de très grande
       taille sans nécessiter O(n²) espace de mémoire centrale. Pour
       certains graphes cependant, comme les arbres, les graphes
       d'intersections, ou les graphes géométriques, une structure de
       données en O(n) est utilisée. Pour les formats d'affichage
       liste, matrix, et smatrix une structure de données de taille
       linéaire (en O(n+m) où m est le nombre d'arêtes) est utilisée
       en interne. Ces trois derniers formats sont donc à éviter. Pour
       la génération de très grand graphe, le format standard ou dot
       doit être privilégié.

   COMMANDES EXTERNES

       Le programme fait appel, pour certaines fonctions, aux
       commandes systèmes suivantes qui doivent être installées: sed,
       grep, awk, more, sort, dot.


OPTIONS


....-help [word], ? [word], or [option|graph] ?
....
       Affiche l'aide en ligne qui est contenue dans le fichier source
       du générateur. Pour cela, le code source .c doit être dans le
       même répertoire que l'exécutable. Si "word" est précisé, alors
       les options et noms de graphes contenant "word" sont
       affichés. La variante "[option|graph] ?" affiche une aide
       détaillée sur une option ou un graphe précis.
....
       Ex: gengraph ? arbre
           gengraph ktree ?
	   gengraph ? hedron
	   gengraph ? planaire
....
       La forme ? peut ne pas fonctionner correctement si un fichier
       d'un seul caractère existe dans le répertoire courant (à cause
       de l'interprétation du shell). Il faut alors utiliser '?' au
       lieu de ?.

....-list
....
       Affiche la liste des graphes et leurs paramètres qu'il est
       possible de générer, d'abord les graphes de bases puis les
       composés. On obtient une aide sur un graphe particulier si son
       nom est suivit de " ?" ou si ses paramètres éventuelles sont
       abscents (dans ce cas il doit être le dernier mot de la ligne
       de commande).

....-version
....
       Affiche la version courante du générateur (en fait du programme
       source), un réel > 1.0. Pour cela, le code source .c doit être
       dans le même répertoire que l'exécutable.

....-directed
....-undirected
....
       L'option -directed produit le graphe orienté en testant les
       n(n-1) arcs possibles, l'option -undirected permettant de
       revenir à la situation par défaut (voir aussi -(no)loop). En
       format standard (ou en dot), un arc apparaît comme x->y au lieu
       de x-y (ou x--y) dans le cas non orienté. Tous les graphes ne
       sont pas forcément définis pour fonctionner correctement avec
       cette option (certaine fonction d'adjacence suppose i<j). La
       plupart des graphes vont apparaître comme orientés symétriques.

....-noloop
....-loop
....
       L'option -noloop permet de ne pas produire les boucles des
       graphes, alors que -loop les autorise. L'option -noloop est
       celle par défaut pour les graphes non-orientés. C'est l'option
       par défaut pour les graphes orientés. Cette option doit être
       placée après -(un)directed.

....-not
....
       Inverse la fonction d'adjacence, et donc affiche le complément
       du graphe. Cette option est prioritaire sur l'option -redirect.

....-dele p
....
       Permet de supprimer chaque arête du graphe générée avec probabilité p.

....-delv p
....
       Similaire à -dele p mais concerne les sommets. Le sommet et ses
       arêtes incidentes sont alors supprimés. Si p est un entier <0,
       alors exactement -p sommets sont supprimés. Si k sommets sont
       supprimés, alors le nom des sommets restant est dans
       l'intervalle [0,n-k[ où n est le nombre de sommets initial du
       graphe. Les noms des sommets sont donc éventuellement
       renumérotés. Voir aussi les options -permute et -shift. Bien
       sûr la fonction d'adjacence adj(i,j) est appliquée sur les noms
       (i,j) originaux.

....-star n
....
       Ajoute n sommets pendant (degré 1) aux sommets du graphe. Si
       n>0, alors n représente le nombre total de sommets ajoutés,
       chacun des n sommets étant connectés aléatoirement uniforménent
       aux sommets du graphe original. Si n<0, alors |n| sommets sont
       ajoutés à chacun des sommets du graphe. Cette opération est
       appliquée en priorité et ne peut être appliquée qu'une fois sur
       la ligne de commande. Les options -star et -apex sont
       incompatibles.

....-apex n
....
       Ajoute n sommets universels, donc connectés à tous les sommets
       du graphe. Cette opération est appliquée en priorité. Les
       options -star et -apex sont incompatibles.

....-redirect p
....
       Redirige chaque arête uniformément avec probabilité p. Plus
       précisément, si {i,j} est une arête du graphe original G, alors
       avec probabilité p l'arête affichée est {i,k} au lieu de {i,j}
       où k est un sommet choisi uniformément parmis les sommets du
       graphe G. Si l'arête {i,j} est supprimée par -dele ou si le
       sommet i est supprimé par -delv, la redirection n'a pas lieu.
       Cette option est appliquée après l'option -not. Le graphe G
       tient donc compte de -not avant de rediriger ses arêtes.

....-seed s
....
       Permet d'initialiser le générateur aléatoire avec la graine s,
       permettant de générer plusieurs fois la même suite aléatoire.
       Par défaut, la graine est initialisée avec le numéro de
       processus de la commande, donc génère par défaut des suites
       différentes à chaque lancement. Le générateur est initialisé
       lorsque l'option est lue sur la ligne de commande. Le
       comportement du programme peut donc être affecté suivant son
       ordre d'apparition. Cependant le graphe est généré après
       l'analyse de la ligne de commande.

....-width m
....
       Limite à m le nombre d'arêtes et de sommets isolés affichés par
       ligne. Cette option n'a pas de signification particulière en
       dehors des formats standard et dot. Par exemple, -width 1
       affiche une arrête (ou un sommet isolé) par ligne. L'option
       -width 0 affiche tout sur une seule ligne. La valeur par défaut
       est 12.

....-shift s
....
       Permet de renuméroter les sommets à partir de l'entier s
       positif. La valeur par défaut est -shift 0.  L'intérêt de cette
       option est de pouvoir réaliser des unions de graphes simplement
       en renumérotant les sommets et en concaténant les fichiers aux
       formats standard ou list. Cette option n'a pas d'effets pour
       les formats de sortie de type matrice.

....-permute
....
       Permute aléatoirement uniformément le nom des sommets
       lorsqu'ils sont affichés.  Les numéros restent dans
       l'intervalle initial qui, sauf si l'option -shift a été
       utilisée, est [0,n[ où n est le nombre de sommets du graphe
       réellement généré. Voir aussi l'option -label.

....-header
....
       Affiche un préambule donnant certaines informations sur le
       graphe, sous forme de commentaire à la C++ (//). Par défaut
       aucun préambule n'est affiché. Les informations affichées sont:
       - l'heure, la date et la graîne du générateur aléatoire
       - la ligne de commande qui a produit la génération du graphe
       - le nombre de sommets, d'arêtes, le degrés maximum et minimum
       Pour les formats standard et dot, le nombre d'arêtes (et les
       degrés min et max) n'est pas déterminé avant l'affichage du
       graphe. Pour cette raison ces nombres ne sont affichés qu'après
       le graphe. Pour n'avoir que les informations sur le graphe,
       utiliser -header avec l'option -format no. Voir aussi -check info.

....-caption title
....
       Permet d'ajouter une légende à un graphe. Cette otpion n'a
       d'effet qu'avec le format dot et ces variantes. Il est possible
       d'affiche la "seed" avec le format %SEED. On ne peut avoir plus
       d'une occurrence du même format (%SEED) dans cette option.
....
       Ex1: gengraph gabriel 30 -caption ex1 -visu
       Ex2: gengraph gabriel 30 -caption "Exemple 2" -visu
       Ex3: gengraph gabriel 30 -caption "graph with seed=%SEED" -visu

....-fast
....
       Génère le graphe sans tester les O(n²) arêtes possibles, mais
       en utilisant directement la liste d'adjacence du graphe
       préalablement générée lors de l'initilisation du graphe, comme
       c'est le cas pour le graphe "load". Le résultat est une
       génération du graphe en temps O(n+m) au lieu de O(n²).
       L'utilisation typique est (voir aussi "loadc"):
....
       Ex: gengraph load file -fast -delv 0.3 -check ncc
....
       permet de calculer le nombre de composantes connexes sur un
       sous-graphe contenu dans un fichier le tout en temps
       linéaire. D'autres graphes peuvent supporter une génération
       rapide si elle est implantée dans l'initialisation de la
       fonction d'adjacence (pour l'instant seul le graphe "load" le
       supporte).  Certaines options, comme -not, n'ont alors plus
       d'effet en présence de -fast. Cependant, -permute, -delv,
       -dele, et d'autres fonctionnent normalement.

....-variant v
....
       Permet de passer un entier v>0 pour contrôler certaines
       fonctionnalités du générateur. Par exemple, "-variant 1 -check
       routing cluster -1". permettra de calculer une variante du
       schéma de routage "cluster" (si elle existe !).

....-check [parameters]
....
       Stocke en mémoire le graphe généré sous la forme d'une liste
       d'adjacence, et lui applique un algorithme. Le graphe est
       généralement affiché, sauf pour certaine options comme -check
       maincc ou -check routing. Utiliser "-format no" pour ne pas
       afficher le graphe généré. Cette option nécessite un espace
       supplémentaire en O(n+m) pour le stockage du graphe.
....
       -check info
....
          Affiche quelques caractéristiques du graphe, s'il est
          orienté, s'il contient des boucles, des multi-arêtes,
          etc. Le graphe lui-même n'est pas affiché.
....
       -check bfs s
....
          Effectue un parcours en largeur d'abord sur le graphe généré
          depuis le sommet s. La distance depuis s est affichée, ainsi
          que l'arborescence (-1 indique que le sommet n'a pas de
          père). La longueur du plus petit cycle passant par s est
          aussi donnée. Elle vaut -1 s'il n'existe pas.
....
       -check dfs s
....
          Effectue un parcours en profondeur d'abord sur le graphe
          généré depuis le sommet s. Le nombre de composantes ainsi
          que l'arborescence (-1 indique que le sommet n'a pas de
          père) sont donnés.
....
       -check ncc
       -check connected
....
          Donne le nombre de composantes connexes ainsi que le nombre
          de cut-vertex du graphe. Ces informations sont aussi
          affichées par -check dfs 0.
....
       -check deg
       -check edge
       -check edges
....
          Affiche la distribution des degrés et le nombre d'arêtes du
          graphe.
....
       -check degenerate
....
          Donne la dégénérescence du graphe, ainsi que l'ordre
          d'élimination correspondant des sommets.
....
       -check girth
....
          Donne la maille du graphe dans le cas non orienté. La valeur
          -1 est renvoyée si le graphe est acyclique, et la valeur 0
          dans le cas orienté.
....
       -check diameter
....
          Calcule le diamètre du graphe généré. Affiche +∞ pour un
          graphe non connexe.
....
       -check radius
....
          Calcule le rayon du graphe généré, soit hauteur du plus
          arbre couvrant. Affiche +∞ pour un graphe non connexe.
....
       -check gcolor
....
          Donne une borne supérieure sur le nombre chromatique du
          graphe en utilisant l'heuristique du degré minimum.
....
       -check kcolor k
....
          Donne une k-coloration du graphe (et la couleur pour chaque
          sommet), si c'est possible. Pour cela une recherche
          exhaustive de toutes les k-colorations est effectuée. Le
          temps est raisonable si k=3 et n<20.
....
       -check kcolorsat k
....
          Donne une formulation SAT de la k-coloration du graphe. Il
          s'agit de la formulation multivaluée classique, un sommet
          pouvant avoir plusieurs couleurs sans que celà nuise à la
          validité du résultat. Les contraintes sont décrites au
          format Dimacs CNF. On peut alors envoyer le résultat à un
          solveur SAT comme MiniSat ou Glucose. Le graphe n'est pas
          affiché, et donc -format no n'est pas nécessaire.
....
	  Ex: gengraph linial 6 3 -check kcolorsat 3 | ./glucose -model
....
       -check kindepsat k
....
          Donne une formulation SAT pour le graphe du problème
          ensemble indépendant de taille k. Les variables i=1 à n sont
          les variables indiquant si le numéroté i-1 est dans la
          soltuion ou pas. Les contraintes sont décrites au format
          Dimacs CNF. On peut alors envoyer le résultat à un solveur
          SAT comme MiniSat ou Glucose. Le graphe n'est pas affiché,
          et donc -format no n'est pas nécessaire.
....
          Pour le problème clique de taille k, il suffit de chercher
          un ensemble indépendant detaille k pour le complément du
          graphe. Et pour le problème "vertex cover" de taille k,
          c'est un ensemble indépendant de taille n-k sur le
          complémentaire qu'il suffit de chercher.
....
       -check ps1
       -check ps1b
       -check ps1c
       -check ps1x n u_1 v_1 ... u_n v_n
....
          Applique le test ps1 ou l'une de ses variantes (voir -filter
          ps1 pour plus de détail sur ce test). Affiche aussi le
          nombre de tests réalisés (nombre de paires de sommets et de
          chemins testés).
....
       -check paths x y
       -check paths2 x y
....
          Liste tous les chemins simples entre les sommets x et
          y. N'affiche rien si x et y ne sont pas connectés. L'ordre
          est défini suivant le premier plus court chemins dans
          l'ordre des sommets depuis le sommet x. La variante paths2
          est similaire sauf que les paths renvoyés n'ont que des
          ensembles de sommets différents.
....
       -check iso H
....
          Teste si le graphe généré G est isomorphe à H. Si oui,
          l'isomorphisme de G à H est donné. Le nombre de tests
          affichés est le nombre de fois où les graphes sont comparés,
          la comparaison prend un temps linéaire en la taille des
          graphes). Plus les graphes sont symétriques (comme un cycle
          ou un hypercube), plus le nombre de tests sera important.
....
          Tester l'isomorphisme entre deux cycles de 8 sommets
          étiquetés aléatoirement prends environ 4 mille tests, et
          entre deux cycles de 12 sommets, 30 millions de tests soit
          9" environ. Pour deux arbres à 75 sommets (aléatoires mais
          isomorphes), moins de 20 tests suffisent.
....
       -check sub H
....
          Teste si le graphe généré G est un sous-graphe couvrant de H
          (donc avec le même nombre de sommets). S'ils ont le même
          nombre d'arêtes, le test est équivalent à l'isomorphisme. Le
          nombre de tests est le nombre total de fois où deux graphes
          sont comparés.  On peut tester si H est Hamiltonien en
          prennant pour G un cycle.
....
          Tester un cycle de longueur 12 dans une grille 3x4 prend
          jusqu'à environ 32 millions de tests (parfois bien moins),
          soit au plus 10".
....
       -check minor H
....
          Teste si le graphe G généré contient H comme mineur. Les
	  graphes peuvent être non connexes. S'ils ont le même nombre
	  de sommets le test est équivalent à celui du sous-graphe
	  (voir -check sub). Dans le cas positif, un modèle de H dans
	  G est fourni.
....
          Le principe consiste à contracter des arêtes de G, de toutes
	  les manières possibles, et à tester si H est un sous-graphe
	  du graphe contracté. Le nombre de tests affichés est le
	  nombre de contractions plus le nombre total de tests
	  réalisés par les tests de sous-graphe. Pour H=K₄ il est
	  préférable d'utiliser -check twdeg qui donne < 3 ssi le graphe
	  ne contient pas K₄ comme mineur.
....
       -check twdeg
....
          Donne une borne supérieure et inférieure sur la treewidth du
          graphe. Pour la borne supérieure, on utilise l'heuristique
          du sommet de degré minimum que l'on supprime et dont on
          complète le voisinage par une clique. En cas d'égalité (même
          degré) on sélectionne le sommet dont il faut rajouter le
          moins d'arêtes. La borne inférieure qui est donnée provient
          de la dégénérescence. La treewidth est exacte si 0,1 ou 2
          est retournée. L'algorithme est en O(n²).
....
       -check tw
....
          Calcule la treewidth du graphe en analysant tous les ordres
          d'éliminations. La complexité est donc en n!. Il ne faut
          l'utiliser que si le nombre de sommets est < 12 (Ex:
          gengraph random 12 .5 -check tw donne 5 en environ 750
          millions de tests). Parfois, l'utilisation de -permute peut
          accélérer le traitement, car partir d'un ordre d'élimination
          déjà bon permet d'éliminer rapidement beaucoup d'ordres
          possibles.
....
       -check maincc
....
          Affiche, dans le mode standard seulement, le graphe
          correspondant à la composante connexe ayant le plus grand
          nombre de sommets. Le graphe initial n'est pas affiché. Les
          sommets sont renumérotés si le graphe initial n'était pas
          connexe. Attention ! l'affichage de la composante n'est
          sensible qu'à l'option -width. En particulier il n'est pas
          possible d'afficher la composante dans un autre format
          (-format) ou avec les noms originaux (-label). Cependant,
          avec "-check maincc | ./gengraph load -" on peut afficher le
          graphe dans le format souhaité, ou ajouter -visu. (Voir
          aussi le raccourcis -maincc.) Notez que "-check maincc
          -visu" provoque une erreur, car -visu applique l'option
          "-format dot" incompatible avec -check maincc.
....
       -check routing [hash h] [scenario s] scheme [parameters]
....
          Construit les tables de routage pour le graphe selon le
          schéma de routage "scheme", ce schéma pouvant comporter des
          paramètres spécifiques. La sortie consiste en statistiques
          sur les tables (taille, temps de calcul) et le graphe.
          L'option "scenario" permet en plus de tester certains types
          (s) de routage sur le graphe (voir ci-après) et d'afficher
          des statistiques sur les longueurs de routes générées (dont
          l'étirement). L'option "hash" permet de préciser la fonction
          de hashage (h) appliquée le cas échéant aux sommets. La
          connectivité du graphe est toujours testée.
....
          Ex: gengraph -permute rplg 200 2.3 -maincc > G1
	      gengraph loadc G1 -check routing scenario all cluster -1
....
          Les scenarii possibles sont (n=nombre de sommets du graphe):
....
            scenario none -> aucun routage (scenario par défaut)
            scenario all -> les n(n-1) routage possibles
	    scenario edges -> tous les routages entre voisins
	    scenario npairs -> n paires aléatoires de sommets différents
	    scenario one u -> les n-1 routages depuis u (choix aléatoire si -1)
	    scenario pair u v -> le routage de u à v (choix aléatoire si -1)
	    scenario pair -p -> le routage depuis p>1 paires aléatoires
....
          Les fonctions de hachages h:[0,n[->[0,k[ possibles sont:
	  (shuffle et mod ont le nombre de collisions minimum ⎡ n/k⎤)
....
	    hash prime -> h(x)=((a*x+b)%p)%k où 0<a,b<p sont aléatoires
	                  et p=2^31-1 est premier (hash par défaut)
	    hash shuffle -> h(x)=π(x)%k où π(x) est une permutation de [0,n[
	                    basée sur deux entiers aléatoires de [0,n[.
	    hash mod -> h(x)=x%k
....
       -check routing cluster k
....
          Schéma de routage "cluster" de paramètre k (routage
          name-independent). Un sommet de degré maximum est choisi
          comme "centre", puis k-1 de ces voisins (de plus haut
          degrés) sont choisis pour former un cluster de taille k. Un
          arbre BFS est enraciné depuis le centre. Chaque sommet
          possède une boule par rayon croissant qui s'arrête avant de
          toucher un sommet du cluster. On route de u vers v d'abord
          dans la boule de u, ou alors on va dans le cluster pour
          chercher un sommet du cluster responsable du hash de v. Une
          fois atteint on route selon l'arbre BFS ou selon un plus
          court chemin si la distance est ≤ logn/loglogn. Les sommets
          voisins du cluster qui ne sont pas eux-mêmes dans le cluster
          possède dans leur table tout leur voisinage.
....
	  Si k=-1, alors k est initialisé à sa valeur par défaut qui
          vaut ⎡ √n⎤. Si k=-2 il est initialisé à n. L'étirement est
          toujours ≤ 5. Il est même ≤ 3 si k=1. Il existe deux
          variantes (en plus de celle par défaut -variant 0). Pour
          -variant 1, le cluster est fixé à tout le voisinage du
          center est le routage dans le centre réalisé selon une
          étoile sans table de routage. Pour -variant 2, le routage
          est réalisé sans les boules de voisinage (qui sont vidées).
          Si de plus k=1, le routage est alors réalisé via la racine
          de l'arbre BFS ce qui réduit au minimum la taille moyenne
          des tables (2 en moyenne).
....
       -check routing dcr k
....
          Schéma de routage name-indépendant "dcr" de paramètre k>0
          représentant le nombre de couleurs. L'étirement est toujours
          ≤ 5 et le nombre d'entrées des tables est en moyenne f(k,n)
          = 2n/k + k*(H(k)+1) où H(k) ~ ln(k)+0.577...  est le k-ième
          nombre harmonic. Le principe du schéma est le suivant.
          Chaque sommet possède une couleur, un entier aléatoire de
          [0,k[, les sommets landmarks étant ceux de couleur 0. Les
          boules de voisinages des sommets sont définies par volume
          comme la plus petite boule contenant au moins chacune des
          couleurs, les sommets du dernier niveau étant ordonnés par
          identifiant croissant. Le routage s->t s'effectue dans la
          boule de s si t y est, sinon on route vers le sommet w de la
          boule de s dont la couleur est égale au hash de t, une
          valeur aussi dans [0,k[. Puis le routage w->t s'effectue
          dans l'arbre BFS enraciné dans le plus proche landmark de s
          ou de t, celui minimisant la distance de w à t.
....
	  Si k=-1, alors k est initialisé à sa valeur optimale
          théorique, celle qui minimise le nombre moyen d'entrées
          f(k,n), valeur calculée et qui vaut environ k ~
          √(n/ln(n))/2, ce qui donne environ 2√(n*ln(n*ln(n))) entrées
          en moyenne.  Si k=-2, il est initialisé à n. Les valeurs de
          k>n sont possibles. Dans ce cas, il s'agit d'un routage de
          plus court chemins comme pour le cas k=n.
....
       -check routing tzrplg t
....
          Schéma de routage étiqueté inspiré de celui de Thorup &
          Zwick pour les graphes RPLG de paramètre réel t (power-law
          exponent) et proposé par Sommer et al. L'étirement est
          toujours ≤ 5. Les valeurs de t entre ]0,1.5] sont
          interdites.  Le schéma utilise des sommets landmarks où des
          arbres BFS sont enracinés, ainsi que des boules (de
          voisinage) définies par rayons croissant qui s'arrête avant
          de toucher un landmark. Le routage s'effectue alors en
          priorité via les boules ou alors via le landmark le plus
          proche de la destination (sans raccourcis), information
          précisée dans l'étiquette de la destination. Les landmarks
          sont les sommets de plus haut degré. Par défaut (-variant 0)
          leur nombre vaut:
....
            - si t>1.5, ⎡ n^((t-2)/2t-3))⎤
	    - si t=0,   ⎡ √n⎤
	    - si t<0,   |t|
....
	  Si -variant 1 et t>1.5, alors les landmarks sont tous les
	  sommets de degré > n^1/(2t-3). Si -variant 2 et t>0, alors
	  les landmarks sont t sommets choisis aléatoirement.
	  L'étirement est ≤ 3 si un seul landmark est choisi.

....-filter family[:range] [not] test [parameters]
....
       Affiche les graphes d'une famille pour lesquels le test est
       vrai (ou faux si "test" est précédé de "not"). Le paramètre
       "family" est un nom de fichier ou "-" pour l'entrée standard.
       La lecture de la famille se fait en temps linéaire.  Il
       contient la famille de graphes (ou un graphe seul) au format
       standard.  L'affichage est influencé par l'option -width qui
       doit être placée avant -filter. La variante "family:range"
       permet de sélectionner les graphes de la famille dont les
       identifiants sont spécifiés par "range", comme par exemple
       "family:5-8" qui sélectionne les graphes d'identifiant
       5,6,7,8. De manière générale, "range" est un ensemble de
       valeurs selon le format "value" décrit ci-après (voir aussi
       -filter F id value). La variante "-:range" est possible.
....
       Dans la suite, la présence de "value" dans les paramètres d'une
       option représente un ensemble de valeurs possibles. Par
       exemple, -filter F vertex '>5' filtre les graphes de la famille
       F comportant plus de 5 sommets. De manière générale, "value"
       est une suite d'intervalles d'entiers sépararés par des ","
       (interprétée comme "ou"), chaque intervalle étant codé comme
       suit:
....
          <x ........ valeur inférieure à x
          >x ........ valeur supérieure à x
          x ou =x ... valeur égale à x
	  x-y ....... valeur dans l'intervalle [x,y]
	  t ......... toujours vrai (intervalle infini)
	  p ......... affiche la valeur plutôt que le graphe
....
       Ex1: -filter F vertex '5,7-13,>100'
       Ex2: -filter F vertex '5-10,p'
       Ex3: -filter F edge p
       Ex4: -filter F id 5,7
....
       L'exemple 1 filtre les graphes de la famille F ayant un nombre
       de sommets n vérifiant soit n=5, soit 7 ≤ n ≤ 13, ou soit
       n>100. L'exemple 2 affiche le nombre de sommets des graphes
       ayant entre 5 et 10 sommets.  L'exemple 3 affiche le nombre
       d'arêtes de chaque graphe.  L'exemple 4 affiche les graphes
       d'identifant 5 et 7 de la famille F.
....
       Pour avoir le maximum/minimum de "value" faire:
       ... -filter ... p | grep '^\[' | sort -rnk 3 | head
       ... -filter ... p | grep '^\[' | sort -nk 3 | head
....
       Si "value" contient le symbole > ou < il est alors préférable
       de mettre des quotes ('>14' par exemple) pour que la commande
       soit correctement interprétée par le shell.
....
       La différence principale avec -check est que le résultat de
       -filter est non verbeux alors que -check, qui ne s'applique pas
       a priori sur des familles de graphes mais sur un graphe seul,
       donne des explications sur l'exécution de l'algorithme. Avec
       -check l'algorithme s'applique au graphe généré, donc a priori
       en temps O(n²), alors qu'avec -filter c'est toujours à partir
       d'un fichier, lu en temps linéaire.
....
       -filter F id value
....
          Filtre les graphes de F dont l'identifiant est déterminé par
          value. Cela permet d'extraire un ou plusieurs graphes
          donnés. C'est équivalent à "-filter F:value all".
....
       -filter F rename shift
....
          Affiche tous les graphes de la famille en renumérotant les
          graphes à partir de l'entier "shift".
....
       -filter F vertex value
....
          Filtre les graphes de F ayant un nombre de sommets déterminé
          par value.
....
       -filter F edge value
       -filter F edges value
....
          Filtre les graphes de F d'un nombre d'arêtes déterminé par
          value.
....
       -filter F all (= vertex t)
....
          Affiche tous les graphes de F ce qui permet en particulier
          de les compter.
....
       -filter F1 minus F2
....
          Affiche F1\F2, c'est-à-dire tous les graphes de F1 qui ne
          sont pas isomorphes à F2 (si F2 est un graphe) ou à l'un des
          graphes de F2 (dans le cas d'une famille).
....
       -filter F1 minus-id F2
....
          Comme "minus" mais concerne les identifiants: supprime de F1
          les graphes dont l'identifiant existe dans F2 (qui peut être
          un graphe ou une famille de graphes).  La complexité est
          environ (|F1|+|F2|)log|F2|, alors que pour "minus" elle est
          en |F1|*|F2|*T où T est le temps pour décider si deux
          graphes pris dans F1 et F2 sont isomorphes.
....
       -filter F minor[-inv] H
....
          Filtre les graphes de F contenant H comme mineur. La
          variante minor-inv filtre les graphes de F qui sont mineurs
          de H. Si H=K₄, il est préférable d'utiliser -filter tw2.
....
       -filter F sub[-inv] H
....
          Filtre les graphes de F contenant H comme sous-graphe,
          chaque graphe de F devant avoir le même nombre de sommets
          que H. La variante sub-inv filtre les graphes de F qui sont
          un sur-graphe de H.
....
       -filter F isub[-inv] H
....
          Filtre les graphes de F contenant H comme sous-graphe
          induit. La variante isub-inv filtre les graphes de F qui
          sont sous-graphes induits de H.
....
       -filter F iso H
....
          Filtre les graphes de F isomorphes à H.
....
       -filter F degenerate value
....
          Filtre les graphes de F de dégénérescence déterminée par
          value.
....
       -filter F forest value
....
          Filtre les graphes de F qui sont des forêts dont le nombre
          d'arbres est déterminé par value.
....
       -filter F isforest (= forest t)
....
          Filtre les graphes de F qui sont des forêts.
....
       -filter F istree (= forest '=1')
....
          Filtre les graphes de F qui sont des arbres.
....
       -filter F cycle (= not forest t)
....
          Filtre les graphes de F contenant au moins un cycle.
....
       -filter F degmax/degmin value
....
          Filtre les graphes de F de degré maximum (ou minimum)
          déterminé par value.
....
       -filter F deg value
....
          Filtre les graphes de F où tous les sommets ont un degré
          déterminé par value. Ainsi -filter deg 4-7 filtre les
          graphes avec un degré minimum au moins 4 et un degré maximum
          au plus 7.
....
       -filter F gcolor value
....
          Filtre les graphes de F dont le nombre de couleurs obtenu
          selon l'heuristique du degré minimum est déterminé par
          value.
....
       -filter F bipartite (= gcolor <3)
....
          Filtre les graphes de F qui sont bipartis.
....
       -filter F component value
....
          Filtre les graphes de F dont le nombre de composantes
          connexes est déterminé par value.
....
       -filter F connected (= component 1)
....
          Filtre les graphes de F qui sont connexes.
....
       -filter F biconnected
....
          Filtre les graphes de F qui sont 2-connexes. Un graphe G est
          k-connexe s'il n'y a pas d'ensemble avec <k sommets qui
          déconnecte G ou laisse G avec 1 sommet. Un graphe est
          2-connexe s'il est connexe, ne possède pas de sommet
          d'articulation et a plus de 2 sommets. Les cliques de taille
          k+1 sont k-connexes.
....
       -filter F radius value
....
          Filtre les graphes de F dont le rayon est déterminé par
          value. Le rayon est la profondeur du plus petit arbre
          couvrant le graphe. Il vaut -1 si le graphe n'est pas
          connexe.
....
       -filter F girth value
....
          Filtre les graphes de F dont la maille est déterminée par
          value. La maille est la taille du plus petit cycle. Elle
          vaut -1 si le graphe n'a pas de cycle. Elle n'est définie
          que si les graphes sont orientés.
....
       -filter F diameter value
....
          Filtre les graphes de F dont le diamètre est déterminé par
          value. Le diamètre vaut -1 si le graphe n'est pas connexe.
....
       -filter F cut-vertex value
....
          Filtre les graphes de F dont le nombre de sommets
          d'articulations est déterminé par value. Un sommet est un
          point d'articulation si sa suppression augmente le nombre de
          composante connexe. Les sommets de degré 1 ne sont pas des
          point d'articulation. Le graphe est biconnexe ssi value<1 ou
          si le graphe est une clique avec au moins deux sommets. On
          peut tester si un graphe est une clique avec -filter degmin
          ou -filter deg.
....
       -filter F ps1
       -filter F ps1b
       -filter F ps1c
       -filter F ps1x n u_1 v_1 ... u_n v_n
....
          Filtre les graphes G de la famille F dont le test ps1 est
          vrai, c'est-à-dire si l'évaluation de la fonction f(G,{})
          décrite ci-après est vraie.
....
	  Soit P un chemin d'un graphe G tel que G\P est connexe. La
          fonction f(G,P) est vraie ssi G est vide (en pratique
          |G|-|P|<3 suffit) ou s'il existe deux sommets x,y de G où y
          n'est pas dans P tels que pour tout chemin Q entre x et y
          dans G "compatible" avec P (c'est-à-dire P et Q
          s'intersectent en exactement un segment) on a les deux
          conditions suivantes: (1) il n'y a pas d'arête entre les
          sommets de P\Q et de G\(Q∪P); et (2) pour toute composante
          connexe C de G\(Q∪P), f(C∪Q,Q) est vraie. Le test est
          optimisé dans un certain nombre de cas, en particulier: les
          arbres (toujours vrai), les cliques (vrai ssi n<5).
....
	  La variante ps1b calcule et affiche de plus un graphe des
          conflits (affichage modifiable par -width), chaque noeud de
          ce graphe correspondant à un argument (C∪Q,Q) évalué à faux
          par f. La valeur (ou code) d'un noeud est 0 (=lourd ou
          faux), 1 (=léger ou vrai) ou - (indéterminée). Suivant
          certaines règles, les valeurs 0 ou 1 sont propagées selon le
          type des arêtes du graphes des conflits.  Résoudre le graphe
          des conflits revient à trouver une affectation des valeurs 0
          ou 1 aux noeuds qui respecte (sans contradiction) toutes les
          règles.
....
	  La fonction f(G,{}) est évaluée à vraie si le graphe des
          conflits n'a pas de solution, c'est-à-dire si une
          contradiction a été découverte ou si pour une paire de
          sommets (x,y) tous ses noeuds sont à 1.
....
	  On affiche le code d'un noeud (0,1,-) ainsi que les sommets
          de sa composante (par ex: [237]).  Les noeuds du graphe des
          conflits sont reliées par des arêtes typées. Les voisins v
          d'un noeud u sont listés avec le type de l'arête, si l'un
          des 4 cas suivants se produit (il n'y a pas d'arête entre u
          et v dans les autres cas):
....
	     v<  (la composante de v est incluse dans celle de u)
	     v>  (la composante de v contient celle de u)
	     v=  (les composantes de u et v sont les mêmes) 
	     v|  (les composantes de u et v sont disjointes) 
....
	  Parmi les règles on trouve par exemple: si deux noeuds du
          graphe des conflits u=(C∪Q,Q) et v=(C'∪Q',Q') sont
          disjoints, c'est-à-dire C n'intersecte pas C', alors seule
          une des expressions f(C∪Q,Q) ou f(C'∪Q',Q') peut être
          fausse, pas les deux. Dit autrement, les composantes de u et
          v ne peuvent pas être "lourdes" (=0) toutes les deux en même
          temps. Et donc, si le code de u est 0, celui de v est
          1. Notons que le code de u et v égale à 1 est comptabible
          avec cette règle.
....
	  La variante ps1c est similaire à ps1b sauf que récursivement
          seul le test ps1 est appliqué, et pas ps1b. Le test ps1c est
          plus long que ps1 mais plus rapide que ps1b. La variante
          ps1x est similaire à ps1b sauf que les valeurs v_i sont
          écrites dans le noeuds u_i du graphe des conflits principal
          (pas ceux générés lors des appels récursifs). Plus
          précisément, v_1 (0 ou 1) est écrit dans le noeud u_1, puis
          sa valeur est propagée. Ensuite v_2 est écrit puis propagée,
          etc.
....
          Dans tous les cas, si G n'est pas connexe, le résultat n'est
          pas déterminé.
....
       -filter F tw value
....
          Filtre les graphes de F selon leur treewidth. L'algorithme
          pour le calcul de la treewidth est assez lent. Pour les
          petites valeurs de tw, des alternatives sont possibles (voir
          -check tw et -filter tw2). Pour savoir si un graphe G est de
          treewidth 3 il suffit de savoir si G contient l'un des 4
          mineurs suivants:
....
          echo "[0]"  > F; ./gengraph clique 5 >> F
	  echo "[1]" >> F; ./gengraph wagner >> F
	  echo "[2]" >> F; ./gengraph prism 5 >> F
	  echo "[3]" >> F; ./gengraph hajos >> F ; echo "0-1-2-0" >> F
	  cat G |./gengraph -filter F minor-inv - -format no
....
       -filter F tw2
....
          Affiche les graphes de F de treewidth ≤ 2. L'algorithme est
          en O(n²). Ce test peut être utilisé pour tester (plus
          rapidement qu'avec -filter minor) les graphes sans mineur
          K₄.
....
       -filter F hyperbol value
....
          Filtre les graphes de F selon leur hyperbolicité. Il s'agit
          de la valeur (entière) maximum, sur tous les quadruplets de
          sommets {u,v,x,y}, de la différence des deux plus grandes
          sommes parmi les sommes de distance : uv+xy, ux+vy et
          uy+vx. La complexité est en O(n⁴).

....-format type
....
       Spécifie le format de sortie. Il est préférable d'utiliser
       cette option en dernier sur la ligne de commande. Les valeurs
       possibles pour "type" sont:
....
       - standard: format standard (liste d'arêtes), c'est le plus compact.
       - list: liste d'adjacence.
       - matrix: matrice d'adjacence.
       - smatrix: matrice supérieure, diagonale comprise.
       - dot: format de GraphViz qui est très proche du format standard.
       - dot<xxx>: dessine le graphe avec GraphViz et converti au format <xxx>.
       - xy: positions X,Y qui ont été utilisées pour le graphe géométrique.
       - no: n'affiche rien, à utiliser en combinaison avec -header ou -check.
....
       Les formats matrix/smatrix/list nécessitent de stocker le
       graphe en mémoire, donc nécessite un espace en O(n+m), alors
       que le graphe est généré à la volée pour les formats standard
       ou dot. Les formats <xxx> pour dot les plus utilisés sont: pdf,
       fig, svg, ps, jpg, gif, png (voir man dot).
....
       L'option -format dot<xxx> est équivalent à "-format dot | dot
       -T<xxx> ...". Elle doit donc être utilisée en dernier sur la
       ligne de commande. Le filtre dot utilisé pour dessiner le
       graphe peut être spécifié par l'option -dotfilter. L'affichage
       des noms de sommets est contrôlé par l'option -label.
....
       Remarque: les positions affichées dans le format dot
       ([pos="..."]) diffèrent d'un facteur proportionnel à √n par
       rapport aux positions originales du graphe (qui peuvent être
       affichées par -format xy). Ce facteur permet de garder une
       taille raisonable pour les sommets car sous dot les sommets ont
       une taille fixe minimale.

....-vcolor option [parameters]
....
       Ces options permettent de modifier la couleur des sommets. Ces
       options n'ont d'effets qu'avec le format dot (et ses variantes
       y compris -visu).  Par défaut les sommets sont de couleur
       noire. Notez que les attributs par défaut des sommets
       (couleurs, formes, etc.)  peuvent être modifiés directement par
       dot (voir l'option -N de dot). Cependant l'option -vcolor
       permet d'individualiser la couleur d'un sommet, en fonction de
       son degré par exemple. Il peut avoir plusieurs options -vcolor
       sur la ligne de commande.
....
       -vcolor deg[r]
....
          La couleur dépend du degré du sommet (deg) ou du rang du
          degré du sommet (degr). Ainsi, les sommets de plus petit
          degré obtiennent la première couleur de la palette, les
          sommets de plus grand degré la dernière couleur de la
          palette, et les autres sommets une couleur intermédiaire de
          la palette. Donc une seule couleur est utilisée si le graphe
          est régulier.
....
       -vcolor degm
....
          Effectue une coloration propre (deux sommets voisins ont des
          couleurs différentes) suivant l'heuristique du degré
          minimum: récursivement, le sommet de degré minimum obtient
          la plus petite couleur qui n'est pas utilisée par ses
          voisins. Cela donne des colorations avec assez peu de
          couleurs pour les graphes de faible arboricité (planaire,
          tw, pw, kout, expander, ...) ou de faible degré. Avec cette
          technique, les graphes bipartis (tree, crown, ...) sont
          coloriés avec deux couleurs. Cette option nécessite un
          espace et un temps en O(n+m).
....
       -vcolor randg
....
          Effectue une coloration propre en utilisant un algorithme
          glouton sur un ordre aléatoire des sommets: récursivement,
          le sommet d'indice i obtient la plus petite couleur qui
          n'est pas utilisée par ses voisins d'indice j<i. Cette
          option nécessite un espace et un temps en O(n+m).
....
       -vcolor kcolor k
....
          Effectue une k-coloration propre du graphe, si c'est
          possible. Si cela n'est pas possible, la première couleur
          est appliquée à tous les sommets. L'algorithme (exponentiel)
          est le même que celui utilisé pour -check kcolor.
....
       -vcolor pal grad
....
          Permet de fixer la palette de couleurs utilisée par les
          sommets. Le paramètre "grad" est un mot sur l'alpabet [a-z]
          (sans les guillemets). Les caractères en dehors de cet
          alphabet sont ignorés. Chaque lettre correspond à une
          couleur de base:
....
	  a=aquamarine     h=hotpink      o=olive         v=violet
	  b=blue           i=indigo       p=purple        w=white
	  c=cyan           j=orange       q=pink          x=gray
	  d=darkorange     k=khaki        r=red           y=yellow
	  e=chocolate      l=lavender     s=salmon        z=black
	  f=forestgreen    m=magenta      t=teal
	  g=green (lime)   n=navy         u=yellowgreen
....
          La palette est calculée selon une interpolation linéaire
          entre les points définis par le mot "grad". Par exemple, si
          "grad" vaut rb, la palette sera composée d'un dégradé allant
          du rouge (r) au bleu (b). Si "grad" vaut rgbr, le dégradé
          ira du rouge au vert puis au bleu et enfin au rouge. Pour
          avoir une couleur (de base) unique, disons w, sur tous les
          sommets, poser "grad" égale à w. Par exemple, pour avoir
          tous les sommets blancs, on peut faire:
....
          gengraph gabriel 30 -vcolor deg -vcolor pal w -visu
....
          La palette par défaut correspond au mot "grad" suivant:
          redjykugfocatbhsqvmpinzxlw. On peut visualiser la palette
          avec l'option "-vcolor list".
....
       -vcolor list
....
          Produit l'affichage de la palette des couleurs utilisées
          pour un graphe plutôt que le graphe lui-même. Cela permet en
          particulier de savoir combien de couleur ont été utilisées.
          La palette est générée en affichant au format dot un graphe
          particulier où les sommets (représentés par un rectangle)
          sont les couleurs utilisées. Utilisez -visu pour visualiser
          la palette sous forme pdf. Le nom des sommets correspond à
          la lettre de la couleur de base comme spécifié par -vcolor
          pal.
....
	  Ex1: gengraph gabriel 50 -vcolor degm -vcolor list
	  (génère la palette utilisée pour ce graphe de Gabriel)
....
          Ex2: gengraph prime 53 -vcolor list
	  (un moyen simple de génèrer la palette par défaut)
....
          Ex3: gengraph clique 100 -vcolor degm -vcolor pal rb -vcolor list
          (génère un dégradé de 100 couleurs allant du rouge au bleu)

....-vsize
....
       La taille des sommets est proportionelle à son degré, alors que
       par défaut elle est fixe. Cette option n'a d'effet qu'avec le
       format dot (et ses variantes). Elle est combinable avec
       -vcolor.

....-visu
....
       Crée un fichier "g.pdf" permettant de visualiser le graphe. Il
       s'agit d'un raccourci de l'option "-format dotpdf" qui rajoute
       également la redirection "> g.pdf" en fin de la ligne de
       commande.

....-maincc
....
       Affiche la composante connexe principale du graphe, les sommets
       étant éventuellement renumérotés si le graphe n'est pas
       connexe. C'est un raccourci pour "-check maincc | ./gengraph
       load - -fast". (Voir aussi -check maincc.) Cet affichage est
       réalisé en temps linéaire grâce à l'option -fast. Les options
       placées avant -maincc affectent le graphe initial alors que
       celles placées après affectent la composante principale.  Les
       options ayant un effet pour les formats hors standard (comme
       -vsize ou -visu) ne devraient être placées qu'après cette
       option.

....-len p
....
       Spécifie la longueur des arêtes pour le format dot et le filtre
       "neato". La valeur par défaut est 1.0, et une valeur plus
       grande (comme 2.0 ou 3.0) alonge les arêtes et permet dans
       certain cas de mieux visualiser le graphe. C'est parfois
       nécessaire pour éviter l'intersection des sommets lorsqu'on
       utilise -label 1.

....-dotfilter filter
....
       Spécifie le filtre de GraphViz, c'est-à-dire l'algorithme de
       dessin utilisé par dot. Par défaut, le filtre est "neato". Les
       filtres principaux sont: dot, neato, twopi, circo, fdp, sfdp,
       ...  Faire "dot -K ." pour afficher les filtres disponibles.

....-pos b
....
       Active (b=1) ou désactive (b=0) la génération des positions des
       sommets pour le format dot. Cela sert à indiquer à l'algorithme
       de dessin dot de respecter (b=1) ou pas (b=0) les coordonnées
       des sommets. L'option par défaut est -pos 0, mais cette option
       est activée pour tous les graphes géométriques (udg, gabriel,
       thetagone, ...).

....-label b
....
       Active (b>0) ou désactive (b=0) l'affichage du nom des sommets
       pour les formats dot et standard. Si b=1, il s'agit du nom
       original du sommet. Cette représentation n'est pas implémentée
       pour tous les graphes. Par défaut les noms sont les entiers de
       [0,n[ où n est le nombre de sommets du graphe généré. L'option
       -label 1 -visu permet alors d'afficher sur le dessin du graphe
       le nom des sommets. Ils ne le sont pas par défaut (b=0).
       L'option -label 2 -visu force l'affichage des noms sous-forme
       d'entiers de [0,n[. Comme cette option influence l'option
       -format dot<xxx>, l'option -label devant être placée avant
       -format.  L'option -label 1 annule l'option -permute, mais
       -label 2 ne le fait pas.
....
       Ex1: gengraph petersen -label 1 -width 1
       Ex2: gengraph petersen -label 1 -format dot | grep label
       Ex3: gengraph petersen -label 1 -len 2 -visu

....-norm ℓ
....
       Fixe la norme pour l'adjacence de certains graphes géométriques
       (udg, gabriel, rng, nng). Les valeurs possibles pour ℓ sont
       1,2,3,4 pour respectivement les normes L_1, L_2, L_max,
       L_min. La norme par défaut est L_2, la norme Euclidienne.

....-xy option [parameters]
....
       Cette option contrôle la façon dont sont générées les
       coordonnées des sommets d'un graphe géométrique. Par défaut les
       positions sont tirées aléatoirement uniformément dans le carré
       [0,1[ × [0,1[, mais cela peut être changé par l'option -xy.
       Notez bien que, même si c'est improbable, deux sommets peuvent
       avoir les mêmes positions (voir l'option -xy unique). Il est
       possible de visualiser les points issus des options -xy (voir
       le graphe "point n").
....
       -xy load file
....
          Charge les positions à partir du fichier "file" ou de
          l'entrée standard si file=-. Cela permet de tester les
          adjacences d'un graphe géométrique à partir de positions
          pré-déterminées. Le format est celui de -format xy.
....
          Ex: gengraph gabriel 10 -xy load file.pos
....
	  Le nombre de sommets du graphe est déterminé par le fichier
          et non par les paramètres du graphe. Cette option n'a
          d'effet que pour les graphes géométriques. La structure du
          fichier texte doit être:
....
	         n
		 x_1 y_1
		 x_2 y_2
		 ...
		 x_n y_n
....
	  où n est le nombre de positions. Les positions x_i y_i ne
	  sont pas forcément dans l'intervalle [0,1[. Notez qu'avec
	  l'option -format xy, il est possible d'effectuer la
	  transformation d'un fichier de positions. L'exemple suivant
	  normalise les coordonnées du fichier g.pos dans le carré
	  unité:
....
          Ex: gengraph -xy load g.pos -xy box 1 1 -format xy
....
       -xy box a b
....
          Effectue un redimensionement des positions de sorte quelles
          se situent dans le rectangle [0,a[ × [0,b[. En prenant
          a=b=1, les coordonnées seront renormalisées dans le carré
          [0,1[ × [0,1[. Cette opération est effectuée juste avant la
          génération des arêtes, en particulier après avoir effectué
          l'opération -xy noise (voir ci-après) et/ou -xy load. Les
          positions sont recadrées pour laisser une fine bande vide
          sur le bord du rectangle. Cette bande a une largeur ~ a/√n
          et une hauteur ~ b/√n où n est le nombre de points.
....
       -xy grid n
....
          Ajoute une grille n × n au graphe généré, ce qui est utile
          lorsque les coordonnées des points sont entiers.
          Techniquement, on ajoute au format de sortie dot un
          sous-graphe représentant la grille où les sommets et les
          arêtes sont de couleur grise. Si n<0, alors le paramètre est
          initialisé à 1+⎣ √N⎦ ou bien à N si l'option "-xy
          permutation" est présente, N étant le nombre de sommets du
          graphe.
....
       -xy vsize f
....
          Facteur de grossissement des sommets pour le format dot. Par
          défaut f=1.0.
....
       -xy noise r p
....
          Effectue une pertubation aléatoire sur les positions des
	  sommets. Le déplacement de chaque sommet est effectué dans
	  sa boule de rayon r (pour p>0) selon une loi en puissance de
	  paramètre p. Prendre p=1 pour une pertubation uniforme dans
	  cette boule, p>1 pour une concentration des valeurs vers le
	  centre et p<1 pour un écartement du centre. Les valeurs <0
	  de p donne des écartements au delà du rayon r.
....
          Plus précisément, une direction (angle de 0 à 2π) est
	  choisie aléatoirement uniformément, puis, selon cette
	  direction, un décalage aléatoire est effectué selon une loi
	  en puissance: si x est uniforme dans [0,1[, le décalage sera
	  d(x)=r*x^p.  Après cette opération, il est possible que les
	  points ne soient plus dans le rectangle d'origine, ce qui
	  peut bien sûr être corrigé par -xy box. Cette option n'a
	  pas d'effet avec -xy permutation et -xy mesh.
....
       -xy seed k p
....
          Fixe k graines pour la génération aléatoire des coordonnées.
	  Les graines sont choisies uniformément dans le carré [0,1[ ×
	  [0,1[ puis centrées par rapport à leur barycentre. Chaque
	  point est alors tiré aléatoirement autour d'une des graînes
	  et à une distance variant selon une loi en puissance de
	  paramètre p avec un rayon r ~ √(ln(k)/k) (voir "noise").
....
       -xy permutation
....
          Génère les points correspondant à une permutation P
          aléatoire uniforme. Le point i à aura pour position
          (i,P(i)).
....
       -xy mesh a b
....
          Génère tous les points de coordonnées entières (i,j) de
          [0,a[x[0,b[ correspondant aux sommets d'une grille a x b.
....
       -xy round p
....
          Arrondi les coordonnées à 10^-p près. Il faut que p soit un
          entier < DBL_DIG, soit p<15 en général. Donc p=0 arrondi à
          l'entier le plus proche. Cet opérateur est appliqué après
          -xy box. Il sert aussi à préciser le nombre de décimales à
          afficher pour l'option -format xy (par défaut p=6). Par
          exemple, la combinaison -xy box 100 100 -xy round -1 permet
          d'avoir des coordonnées multiples de 10.
....
       -xy unique
....
          Supprime les sommets en double, correspondant aux mêmes
          positions. Cela peut être utile lorsqu'on utilise -xy round
          par exemple. Cette opération est appliquée après toutes les
          autres, notamment après -xy box et -xy round. Ceci est
          réalisé à l'aide d'un tri rapide en temps O(nlogn).


   GRAPHES

       Deux types de graphes sont possibles : les graphes de base et
       les graphes composés. Ces derniers sont obtenus en paramétrant
       un graphe de base. Une catégorie importante de graphes sont les
       graphes géométriques (qui peuvent être composé ou de base).
       L'adjacence est déterminée par les coordonnées associées aux
       sommets. De nombreuses options s'y réfèrent.  Ils activent tous
       par défaut l'option -pos. Les graphes orientés activent quant à
       eux tous l'option -directed.
       


   GRAPHES DE BASE :

....grid k n_1 ... n_k
....
       Grille à k dimensions de taille n_1 × ... × n_k. Si la taille
       n_i est négative, alors cette dimension est cyclique.  Par
       exemple, "grid 1 -10" donnera un cycle à 10 sommets.

....ring n k c_k ... c_k
....
       Anneaux de cordes à n sommets chacun ayant k cordes de longueur
       c_1,...,c_k.

....cage n k c_1 ... c_k
....
       Graphe cubique pouvant servir à la construction de graphes
       n-cage, c'est-à-dire aux plus petits graphes cubique à n
       sommets de maille donnée. Ils sont toujours Hamiltoniens. Ils
       peuvent être vus comme des anneaux de cordes irréguliers. Ils
       sont construits à partir d'un cycle de longueur n découpé en
       n/k intervalles de k sommets. Le i-ème sommet de chaque
       intervalle, disons le sommet numéro j du cycle, est adjacent au
       sommet numéro j+c_i du cycle (modulo n). Les valeurs c_i peuvent
       être positives ou négatives. Cette fonction permet aussi de
       construire des graphes avec des sommets de degré 4 comme "cage
       8 2 0 2" (voir aussi le graphe de Chvátal) ou avec des sommets
       de degré 2 comme "cage 4 2 2 0".

....arboricity n k
....
       Graphe d'arboricité k à n sommets aléatoire. Ce graphe est
       composé de l'union de k>0 arbres aléatoires. Il est donc
       toujours connexe. Chacun des arbres est un arbre plan enraciné
       aléatoire uniforme dont les sommets sont permutés
       aléatoirement, sauf le premier arbre dont les sommets sont
       numérotés selon un parcours en profondeur. Ces graphes
       possèdent au plus k(n-1) arêtes, et pour k=1 il s'agit d'un
       arbre.

....rarytree n b z
....
       Arbre b-aire plan aléatoire uniforme à n noeuds internes. Il
       faut b≥2. Il possède bn+1+z sommets, z étant un paramètre
       valant 0 ou 1. La racine est de degré b+z, les autres sommets
       sont de degré b+1 (soit b fils) ou 1 (=feuille). Les sommets
       sont numérotés selon un parcours en profondeur modifié: tous
       les fils du sommet courant sont numérotés avant l'étape de
       récursivité. Si n=1, alors le graphe est une étoile à b+z
       feuilles. Le dessin avec dot (-visu) ne respecte pas le
       plongement de l'arbre.

....ringarytree h k r p
....
       Arbre de hauteur h où chaque noeud interne à exactement k fils,
       le degré de la racine étant de degré r (p=0). Si p=1 alors un
       chemin entre les sommets de même niveau est ajouté. Si p=2,
       c'est un cycle. Notez que "arytree h 1 r 0" génère une étoile
       de degré r où chaque branche est de longueur h. Le nom des
       sommets correspond au chemin depuis la racine.

....kpage n k
....
       Graphe k-pages connexe aléatoire. Un graphe k-page peut être
       représenter en plaçant les sommets le long d'un cercle, en
       dessinant les arêtes comme des segments de droites, et en
       coloriant les arêtes en k>0 couleurs de façon à ce que les
       arêtes de chaque couleur induisent le dessin d'un graphe
       planaire-extérieur. La numérotation des sommets est faite le
       long du cercle. Les graphes 1-page sont les graphes
       planaires-extérieurs, les 2-pages sont les sous-graphes de
       graphes planaires Hamiltoniens. Les graphes planaires de degré
       au plus 4 sont 2-pages, les 3-arbres planaires (ou graphes
       Apolloniens) sont 3-pages, et les cliques avec 2k-1 ou 2k
       sommets des k-pages.
....
       Ces graphes sont construits par le processus aléatoire suivant.
       On génère k graphes planaires-extérieurs aléatoires uniformes
       connexes à n sommets (plan et enraciné) grâce à une bijection
       avec les arbres plans enracinés dont tous les sommets, sauf
       ceux de la dernière branche, sont bicoloriés. On fait ensuite
       l'union de ces k graphes en choisissant aléatoirement la racine
       des arbres, sauf celui du premier planaire-extérieur, ce qui
       correspond à une permutation circulaire des sommets sur la face
       extérieure.

....ktree n k
....
       k-arbre aléatoire à n sommets. Il faut n>k, mais k=0 est
       possible. Il est généré à partir d'un arbre enraciné aléatoire
       uniforme à n-k noeuds de manière similaire à "tree n-k". Cela
       constitue les "sacs" que l'on remplit avec les n sommets comme
       suit: on met k+1 sommets dans le sac racine connecté en clique,
       puis, selon un parcours en profondeur de l'arbre, on met un
       sommet différent pour chacun des autres sacs. Ce sommet est
       alors connectés à exactement k sommets choisis aléatoirement
       dans le sac parent et sont ajoutés à son sac.

....kpath n k
....
       k-chemin aléatoire à n sommets. La construction est similaire à
       celle utilisée pour ktree. L'arbre aléatoire est remplacé par
       un chemin. Ces graphes sont des graphes d'intervalles
       particuliers (voir "interval n").

....rig n k p
....
       Graphe d'intersections aléatoire (Uniform Random Intersection
       Graph).  Il possède n sommets, chaque sommet u étant représenté
       par un sous-ensemble S(u) aléatoire de {1,...,k} tel que chaque
       élément appartient à S(u) avec probabilité p. Si p<0, alors p
       est fixée au seuil théorique de connectivité, à savoir
       p=√(ln(n)/(nk)) si k>n et p=ln(n)/k sinon. Il y a une arête
       entre u et v ssi S(u) et S(v) s'intersectent. La probabilité
       d'avoir une arête entre u et v est donc Pₑ=1-(1-p²)^m, mais les
       arêtes ne sont pas indépendantes (Pr(uv|uw)>Pr(uv)). En
       général, pour ne pas avoir Pₑ qui tend vers 1, on choisit les
       paramètres de façon à ce que kp²<cste. Lorsque k≥n³, ce modèle
       est équivalent au modèle des graphes aléatoires d'Erdös-Reny
       (voir random n p).

....apollonian n
....
       Graphe Apollonien aléatoire uniforme à n≥4 sommets. Les graphes
       Apolloniens sont les 3-arbres planaires ou encore les graphes
       planaires maximaux chordaux. Ils sont obtenus en subdivisant
       récursivement un triangle en trois autres. Ils sont
       3-dégénérés, de treewidth 3, et de nombre chromatique 4. La
       distance moyenne est ϴ(logn). Ils sont en bijection avec les
       arbres ternaires à n-3 noeuds internes.

....polygon n
....
       Triangulation aléatoire uniforme d'un polygone convexe à n≥3
       cotés. Ce sont aussi des graphes planaires-extérieurs maximaux
       aléatoires. Ils sont Hamiltoniens, 2-dégénérés, de treewidth 2,
       et de nombre chromatique 3.  Ils sont en bijection avec les
       arbres binaires à n-2 noeuds internes. La numérotation des
       sommets n'est pas cyclique le long du polygone.
....
       Ex: gengraph polygon 20 -dotfilter circo -visu

....planar n f d
....
       Graphe planaire aléatoire composé de n faces internes de
       longueur f≥3, les sommets internes étant de degré au moins d et
       ceux de la face externe au moins 2. Ils possèdent au plus
       n(f-2)+2 sommets, sont 2-connexes, 2-dégénérés, de maille f. Si
       d>4 alors ils sont d'hyperbolicité O(f). Ils sont construits en
       ajoutant itérativement les faces par le processus aléatoire
       suivant. Au départ, il s'agit d'une cycle de longueur f. Pour
       chaque nouvelle face on ajoute un sommet u que l'on connecte à
       un sommet quelconque du cycle C formant le bord de la face
       extérieure du graphe courant. Puis on ajoute alors un chemin
       allant de u à un sommet v de C de façon à respecter la
       contrainte des degrés des sommets qui vont devenir internes et
       la contrainte sur la longueur de la nouvelle face créée. Le
       sommet v est choisit uniformément parmi tous les sommets
       possibles de C respectant les contraintes. Si d<0, alors on
       fait comme si d=+∞ et le résultat est un graphe
       planaire-extérieur Hamiltonien, c'est-à-dire 2-connexe. Si f<0,
       alors chaque face créée est de longueur au plus |f|, longueur
       aléatoire uniforme de [3,|f|]. Si f=d=4, il s'agit d'un
       "squaregraph". Les valeurs d=0,1,2 sont équivalentes.

....kneser n k r
....
       Graphe de Kneser généralisé. Le graphe de Kneser K(n,k)
       classique est obtenu avec r=0. Les sommets sont tous les
       sous-ensembles à k éléments de [0,n[ (il faut donc 0 ≤ k ≤
       n). Deux sommets sont adjacents ssi leurs ensembles
       correspondant ont au plus r éléments en commun. Le nombre
       chromatique de K(n,k), établit par Lovász, vaut n-2k+2 pour
       tout n≥2k-1>0. Le graphe de Petersen est le graphe K(5,2). Ils
       ont un lien avec les graphes de Johnson J(n,k).

....gpetersen n r
....
       Graphe de Petersen généralisé P(n,r), 0≤r<n/2. Ce graphe
       cubique possède 2n sommets qui sont u_1,..., u_n, v_1,...,
       v_n. Les arêtes sont, pour tout i, u_i-u_{i+1}, u_i-v_i, et
       v_i-v_{i+r} (indice modulo n). Il peut être dessiné tel que
       toute ses arêtes sont de même longueur (unit distance
       graph). Ce graphe est biparti ssi n est pair et r est
       impair. C'est un graphe de Cayley ssi r² = 1 (modulo n). P(n,r)
       est Hamiltonien ssi r≠2 ou n≠5 (modulo 6). P(n,r) est isomorphe
       à P(n,(n-2r+3)/2)).  P(4,1) est le cube, P(5,2) le graphe de
       Petersen, P(6,2) le graphe de Dürer, P(8,2) le graphe de
       Möbius-Kantor, P(10,2) le dodécaèdre, P(10,3) le graphe de
       Desargues, P(12,5) le graphe de Nauru, P(n,1) un prisme.

....antiprism n
....
       Graphe composé de deux cycles à n sommets connectés par 2n
       triangles. Le prisme est similaire sauf que pour ce dernier les
       deux cycles sont connectés par des carrés. Il est ainsi
       planaire, 4-régulier, possède 2n sommets et 4n arêtes. C'est
       aussi le dual du trapézoèdre n-gonal.

....rpartite k a_1 ... a_k
....
       Graphe k-parti complet K_{a_1,...,a_k}. Ce graphe possède
       n=a_1+...+a_k sommets. Les sommets sont partitionnés en k parts
       comme suit: part 1 = [0,a_1[, part 2 = [a_1,a_1+a_2[, ...  part
       k = [a_1+...+a_{k-1},n[. Les sommets i et j sont adjacents ssi
       i et j appartiennent à des parts différentes.

....ggosset p k d_1 v_1 ... d_1 v_k
....
       Graphe de Gosset généralisé. Les sommets sont tous les vecteurs
       entiers de dimension d = d_1 + ... + d_k dont les coordonnées
       comprennent, pour i=1...k, exactement d_i fois la valeur
       v_i. Il existe une arête entre les vecteurs u et v si et
       seulement le produit scalaire entre u et v vaut l'entier p. Des
       valeurs intéressantes sont par exemple: 1 2 2 -1 2 0 ou encore
       8 2 2 3 6 -1.

....crown n
....
       Graphe biparti à 2n sommets où le i-ème sommet de la première
       partie est voisin au j-ème sommet de la seconde partie ssi
       i≠j. Pour n=3, il s'agit du cycle à 6 sommets, pour n=4, il
       s'agit du cube (à 8 sommets).

....centipede n
....
       Arbre à 2n sommets et n feuilles en forme de peigne.

....fan p q
....
       Graphe composé d'un chemin à p sommets et de q sommets, chacun
       connectés à tous ceux du chemin. Le graphe classique "fan n"
       correspond à p=n et q=1.

....flip n
....
       Graphe des flips des triangulations d'un polygone convexe à n>2
       sommets. Les sommets, qui sont les triangulations, sont en
       bijection avec des arbres binaires complets à 2n-3 noeuds qui
       sont codés par les mots de Dyck de longueur 2n-4 (que l'on peut
       afficher avec -label 1). Le nombre de sommets est donc C(n-2)
       nombre de Catalan d'ordre n-2. Les adjacences peuvent être vues
       aussi comme des rotations d'arbres. Le diamètre est 2n-10 pour
       n>12. Le nombre chromatique n'est pas connu, on sait pas s'il
       est constant ou pas. Il vaut 3 pour n=5..9, et 4 pour n=10 et
       11.

....interval n
....
       Graphe d'intersection de n intervalles d'entiers aléatoires
       uniformes pris dans [0,2n[. Des graphes d'intervalles peuvent
       aussi être générés par "kpath n k".

....permutation n
....
       Graphe de permutation sur une permutation aléatoire uniforme
       des entiers de [0,n[.

....prime n
....
       Graphe à n sommets tel que i est adjacent à j ssi i>1 et j
       divisible par i.

....paley n
....
       Graphe de Paley à n sommets. Deux sommets sont adjacents ssi
       leur différence est un carré modulo n. Il faut que n soit la
       puissance d'un nombre premier et que n=1 (mod 4), mais le
       graphe est aussi défini pour les autres valeurs. Les premières
       valeurs possibles pour n sont: 5, 9, 13, 17, 25, 29, 37, 41,
       49, ... Ces graphes sont Hamiltoniens. Si n est simplement
       premier, alors ils sont de plus auto-complémentaires et
       réguliers.  Paley 17 est le plus grand graphe G où ni G ni son
       complémentaire ne contient K₄, d'où Ramsey(4)=18.

....mycielski k
....
       Graphe de Mycielski de paramètre (nombre chromatique) k. C'est
       un graphe sans triangle, k-1 (sommets) connexe, et de nombre
       chromatique k. Le premier graphe de la série est M2 = K2, puis
       on trouve M3=C5, M4 est le graphe de Grötzsch à 11 sommmets.

....windmill n
....
       Graphe composé de n cycles de longueur trois ayant un sommet commun.

....barbell n1 n2 p
....
       Graphe des haltères (Barbell Graph) composé de deux cliques de
       n1 et n2 sommets reliées par un chemin de longueur p. Il
       possède n1+n2+p-1 sommets. Si p=0 (p=-1), le graphe est composé
       de deux cliques ayant un sommet (une arête) en commun. Plus
       généralement, si p≤0, le graphe est composé de deux cliques
       s'intersectant sur 1-p sommets.

....chess p q x y
....
       Graphe composé de p x q sommets représentant les cases d'un
       échiquier p x q, deux cases étant connectée s'il existe un
       déplacement d'une case vers l'autres avec un saut de x cases
       selon un axe et y selon un autre. Le "knight graph" classique
       est donc un "chess 8 8 2 3", et "chess n 2 1 0" correspond à
       "ladder n".

....sat n m k
....
       Graphe aléatoire issu de la réduction du problème k-SAT à
       Vertex Cover. Le calcul d'un Vertex Cover de taille minimum
       pour ce graphe est donc difficile pour k>2. Soit F une formule
       de k-SAT avec n variables x_i et m clauses CNF de k termes.  Le
       graphe généré par "sat n m k" possède un Vertex Cover de taille
       n+(k-1)m si et seulement si F est satisfiable. Ce graphe est
       composé d'une union de n arêtes et de m cliques de k sommets,
       plus des arêtes connectant certains sommets des cliques aux n
       arêtes. Les n arêtes représentent les n variables, une
       extrémité pour x_i, l'autre pour non(x_i). Ces sommets ont des
       numéros dans [0,2n[, x_i correspond au sommet 2i-2 et non(x_i)
       au sommet 2i-1, i=1...n. Les sommets des cliques ont des
       numéros consécutifs ≥ 2n. Le j-ème sommet de la i-ème clique
       est connecté à l'une des extrémités de l'arête k ssi la j-ème
       variable de la i-ème clause est x_k (ou non(x_i)). Ces
       dernières connexions sont aléatoires uniformes.

....kout n k
....
       Graphe à n sommets k-dégénéré crée par le processus aléatoire
       suivant: les sommets sont ajoutés dans l'ordre croissant de
       leur numéro, i=0,1,...n-1. Le sommet i est connecté à d voisins
       qui sont pris aléatoirement uniformément parmi les sommets dont
       le numéro est < i. La valeur d est choisie aléatoirement
       uniformément entre 1 et min{i,k}. Il faut k>0. Le graphe est
       connexe, et pour k=1, il s'agit d'un arbre.

....expander n k
....
       Graphe à n sommets composé de k>0 cycles Hamiltoniens
       aléatoires. Le degré des sommets varie entre 2 et 2k. Il
       possède le cycle 0,1,...,n-1,0 comme cycle Hamiltonien, et a la
       propriété d'expansion à partir de k>3. Plus précisément, avec
       grande probabilité, les valeurs propres de la matrice
       d'ajacence du graphe sont ≤ 2√(2k).

....icosahedron
....
       Isocahèdre: graphe panaire 5-régulier à 12 sommets. Il possède
       30 arêtes et 20 faces qui sont des triangles. C'est le dual du
       dodécahèdre.

....cuboctahedron
....
       Cuboctaèdre: graphe planaire 4-régulier à 12 sommets. Il
       possède 24 arêtes et 14 faces qui sont des triangles ou des
       carrés. C'est le dual du rhombicdodécahèdre.

....rdodecahedron
....
       Rhombic-dodécaèdre: graphe planaire à 14 sommets avec des
       sommets de degré 3 ou 4. Il possède 21 arêtes et 12 faces qui
       sont des carrés. C'est le dual du cuboctahèdron.

....deltohedron n
....trapezohedron n
....
       Deltoèdre ou trapézoèdre n-gonal: graphe composé de 2n faces en
       cerf-volant (deltoïdes) décalées symétriquement. C'est donc un
       graphe planaire de 2n+2 sommets et 4n arêtes où toutes les
       faces sont des carrées. C'est aussi le dual de l'antiprisme
       n-gonal. Il s'agit d'un cube si n=3.

....tutte
....
       Graphe de Tutte. C'est un graphe planaire cubique 3-connexe à
       46 sommets qui n'est pas Hamiltonien.

....hgraph
....
       Arbre à six sommets et quatre feuilles en forme de H.

....tgraph
....
       Arbre à cinq sommets et trois feuilles en forme de T, aussi
       appelé Fork Graph.

....rgraph
....fish
....
       Graphe à six sommets en forme de R ou de poisson. Il est
       composé d'un cycle à quatre sommets dont l'un possède d'eux
       possède deux sommets pendant.

....cricket
....
       Cricket Graph, graphe à cinq sommets composé d'un triangle où à
       l'un de sommets est attaché deux sommets pendant (de degré 1).

....moth
....
       Moth Graph, graphe à six sommets composé de deux triangles
       partageant une arête et de deux sommets pendant (degré 1)
       attachés à l'un des sommets commun aux deux triangles.

....bull
....
       Bull Graph, graphe à cinq sommets auto-complémentaire en forme
       de A.

....cross
....
       Cross Graph, arbre à six sommets en forme de croix chrétienne.

....harborth
....
       Graphe de Harborth. C'est un graphe planaire 4-régulier à 52
       sommets qui est distance unitaire aussi appelé graphe allumette
       (voir theta0 et diamond). Il peut ainsi être dessiné sans
       croisement d'arête qui ont toutes la même longueur.

....bidiakis
....
       Graphe de Bidiakis ou cube de Bidiakis. C'est un graphe
       planaire cubic à 12 sommets. On peut le représenter comme un
       cube où deux faces opposées comportent une arête supplémentaire
       perpendiculaire joignant deux bord opposés. Il est Hamiltonien
       et son nombre chromatique est 3.

....herschel
....
       Graphe de Herschel. C'est le plus petit graphe planaire
       3-connexe qui ne soit pas Hamiltonien. Il est biparti, possède
       11 sommets et 18 arêtes.

....goldner-harary
....
       Graphe de Goldner-Haray. C'est le plus petit graphe planaire
       maximal qui ne soit pas Hamiltonien. Il possède 11 sommets et
       donc 27 arêtes (voir aussi Herchel). C'est un 3-arbre planaire
       (voir apollonian).

....fritsch
....
       Graphe de Fritsch. Il est planaire maximal à 9 sommets qui peut
       être vu comme un graphe de Hajós dans un triangle. C'est, avec
       le graphe de Soifer, le plus petit contre-exemple à la
       procédure de coloration de Kempe.

....soifer
....
       Graphe de Soifer. Il est planaire maximal à 9 sommets. C'est,
       avec le graphe de Fritsch, le plus petit contre-exemple à la
       procédure de coloration de Kempe.

....poussin
....
       Graphe de Poussin. Il est planaire maximal à 15 sommets. C'est
       un contre-exemple à la procédure de coloration de Kempe.

....headwood4
....
       Graphe de Headwood pour la conjecture des 4 couleurs,
       contre-exemple de la preuve de Kempe. Il est planaire maximal
       avec 25 sommets, est de nombre chromatique 4, de diamètre 5, de
       rayon 3 et Hamiltonien.

....errara
....
       Graphe d'Errara. Il est planaire maximal à 17 sommets. C'est un
       contre-exemple à la procédure de coloration de Kempe.

....kittell
....
       Graphe de Kittell. Il est planaire maximal à 23 sommets. C'est
       un contre-exemple à la procédure de coloration de Kempe.

....frucht
....
       Graphe de Frucht. Il est planaire cubique à 12 sommets. Il n'a
       pas de symétrie non-triviale. C'est un graphe de Halin de
       nombre chromatique 3, de diamètre 4 et de rayon 3.

....treep p
....
       Arbre aléatoire à p>2 feuilles sans sommets internes de degré
       deux. Il possède entre p+1 et 2p-2 sommets. Ce graphe est à la
       base de la construction des graphes de Halin.

....halin p
....
       Graphe de Halin aléatoire basé sur un arbre à p>2 feuilles. Il
       possède entre p+1 et 2p-2 sommets. Il est constitué d'un arbre
       sans sommets de degré deux dont les p feuilles sont connectés
       par un cycle (de p arêtes). Ces graphes planaires de degré
       minimum au moins trois sont aussi arête-minimale 3-connexes,
       Hamiltonien (et le reste après la suppression de n'importe quel
       sommet), de treewidth exactement 3 (ils contiennent K₄ comme
       mineur). Ils contiennent toujours au moins trois triangles et
       sont de nombre chromatique 3 ou 4.

....butterfly d
....
       Graphe Butterfly de dimension d. Les sommets sont les paires
       (x,i) où x est un mot binaire de d bits et i un entier de
       [0,d]. Les sommets peuvent être représentés en d+1 niveaux
       chacun de 2^d sommets, les arêtes connectant les niveaux
       consécutifs. Le sommet (x,i) est adjacent à (y,i+1) ssi les
       bits de x sont identiques à ceux de y sauf pour celui de numéro
       i+1 (le bit 1 étant le bit de poids le plus faible). Il possède
       (d+1)*2^d sommets et d*2^(d+1) arêtes, les sommets de niveau 0
       et d étant de degré 2 les autres de degré 4.

....line-graph n k
....
       Line-graphe aléatoire à n sommets et de paramètre k>0 entier
       Plus k est petit, plus le graphe est dense, le nombre d'arêtes
       étant proportionnel à (n/k)². Si k=1, il s'agit d'une clique à
       n sommets. Ces graphes sont obtenus en choisissant, pour chaque
       sommet, deux couleurs de [1,k]. Deux sommets sont alors
       adjacents ssi ils possèdent la même couleur.  Ces graphes sont
       claw-free (sans K_{1,3} induit). Comme les graphes claw-free,
       les line-graphes connexes avec un nombre pair de sommets
       possède toujours un couplage parfait. On rappel qu'un graphe G
       est le line-graphe d'un graphe H si les sommets de G
       correspondent aux arêtes de H et où deux sommets de G sont
       adjacents ssi les arêtes correspondantes dans H sont
       incidentes.  Pour parle parfois de graphe adjoint.

....shuffle d
....
       Graphe Shuffle-Exchange de dimension d. Les sommets sont les
       mots binaires de d lettres. Le sommet w et w' sont voisins si w
       et w' diffèrent du dernier bit, ou bien si w' peut être obtenu
       par décalage cyclique à droite ou à gauche de w.

....debruijn d b
....
       Graphe de De Bruijn de dimension d et de base b. Il a b^d
       sommets qui sont tous les mots de d lettres sur un alphabet de
       b lettres. Le sommet (x_1,...,x_d) est voisin des sommets
       (x_2,...,x_d,*). Ce graphe est Hamiltonien, de diamètre d et le
       degré de chaque sommet est 2b, 2b-1 ou 2b-2. Pour d=3 et b=2,
       le graphe est planaire.

....kautz d b
....
       Graphe de Kautz de dimension d et de base b. Il a b*(b-1)^(d-1)
       sommets qui sont tous les mots de d lettres sur un alphabet de
       b lettres avec la contrainte que deux lettres consécutives
       doivent être différentes. L'adjacence est celle du graphe de De
       Bruijn. C'est donc un sous-graphe induit de De Bruijn (debruijn
       d b). Il est Hamiltonien, de diamètre d et le degré de chaque
       sommet est 2b-2 ou 2b-3. Pour d=b=3 le graphe est planaire.

....linial n t
....
       Neighborhood graph des cycles introduit par Linial. C'est le
       graphe de voisinage des vues de taille t d'un cycle orienté
       symétrique à n sommets ayant des identifiants uniques de
       [0,n[. Il faut n≥t>0. Le nombre chromatique de ce graphe est k
       ssi il existe un algorithme distribué qui en temps t-1
       (resp. en temps (t-1)/2 avec t impair) peut colorier en k
       couleurs tout cycle orienté (resp. orienté symétrique) à n
       sommets ayant des identifiants uniques et entiers de [0,n[. Les
       sommets sont les t-uplets d'entiers distincts de [0,n[. Le
       sommet (x_1,...,x_t) est voisin des sommets (x_2,...,x_t,y) où
       y≠x_1 si n>t et y=x_1 si n=t.  C'est un sous-graphe induit de
       linialc n t, et donc du graphe de Kautz (kautz t n) et de De
       Bruijn (debruijn t). Le nombre de sommets est m(m-1)...(m-t+1).
       Certaines propriétés se déduisent du graphe linialc n t.

....linialc m t
....
       Neighborhood graph des cycles colorés.  Il s'agit d'une
       variante du graphe linial n t. La différence est que les
       sommets du cycle n'ont plus forcément des identitiés uniques,
       mais seulement une m-coloration, m≤n. L'adjacence est
       identique, mais les sommets sont les (2t+1)-uplets
       (x_1,...,x_{2t+1}) d'entiers de [0,m[ tels que x_i≠x_{i+1}. Il
       s'agit donc d'un sous-graphe induit de linialc m t, lui-même
       sous-graphe induit du graphe de Kautz (kautz 2t+1 m) et donc de
       De Bruijn (debruijn 2t+1 m). Le nombre de sommets est
       m(m-1)^{2t} et son degré maximum est 2(m-1). La taille de la
       clique maximum est 3 si m>2 et t>1. Le nombre chromatique de ce
       graphe est 3 pour m=4, 4 pour 5≤m≤24. Pour 25≤m≤70 c'est au
       moins 4 et au plus 5, la valeur exacte n'étant pas connue.

....pancake n
....
       Graphe "pancake" de dimension n. Il a n! sommets qui sont les
       permutations de {1,...,n} et (n-1)-régulier. Une permutation,
       c'est-à-dire un sommet, est voisine de toutes celles obtenues
       en retournant un de ces préfixes. Plus précisément, les sommets
       x=(x_1,...,x_n) et y=(y_1,...,y_n) sont adjacents s'il existe
       un indice k tel que y_i=x_i pour tout i>k et y_i=x_{k-i} sinon.
       Son diamètre, qui est linéaire en n, n'est pas connu
       précisément. Les premières valeurs connues, pour n=1...17,
       sont: 0, 1, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18,
       19. Donc les diamètres 2,6,12 n'existent pas.

....bpancake n
....
       Graphe "burn pancake" de dimension n. Il a n!*2^n sommets qui
       sont les permutations signées de {1,...,n}. Les sommets
       x=(x_1,...,x_n) et y=(y_1,...,y_n) sont adjacents s'il existe
       un indice k tel que y_i=x_i pour tout i>k et y_i=-x_{k-i}
       sinon. Dit autrement la permutation de y doit être obtenue en
       retournant un préfixe de x et en inversant les signes. Par
       exemple, le sommet (+2,-1,-5,+4) est voisin du sommet
       (+5,+1,-2,+4). Comme le graphe pancake, c'est un graphe
       (n-1)-régulier de diamètre linéaire en n.

....gpstar n d
....
       Graphe "permutation star" généralisé de dimension n. Il a n!
       sommets qui sont les permutations de {1,...,n}. Deux sommets
       sont adjacents si leurs permutations diffèrent par d
       positions. C'est un graphe régulier.

....pstar n
....
       Graphe "permutation star" de dimension n. Il a n! sommets qui
       sont les permutations de {1,...,n}. Deux sommets sont adjacents
       si une permutation est obtenue en échangeant le premier élément
       avec un autre. Le graphe est (n-1)-régulier. Le graphe est
       biparti et de diamètre ⎣ 3(n-1)/2)⎦. C'est un sous-graphe induit
       d'un "gpstar n 1".

....hexagon p q
....
       Grille hexagonale p x q. C'est un planaire composé de p rangées
       de q hexagones, le tout arrangé comme un nid d'abeille. Ce
       graphe peut aussi être vu comme un mur de p rangées de q
       briques, chaque brique étant représentée par un cycle de
       longueur 6. Il possède (p+1)*(2p+2)-2 sommets et est de degré
       maximum 3. Sont dual est le graphe whexagon.
....
       Ex: gengraph hexagon 20 20 -dele 0.2 -maincc -visu

....whexagon p q
....
       Comme le graphe hexagon p q sauf que chaque hexagone est
       remplacé par une roue de taille 6 (chaque hexagone possède un
       sommet connecté à ses 6 sommets). C'est le dual de l'hexagone.
       Il possède p*q sommets de plus que l'hexagone p q.

....hanoi n b
....
       Graphe de Hanoï généralisé, le graphe classique est obtenu avec
       b=3. Il est planaire avec b^n sommets et est défini de manière
       récursive. Le niveau n est obtenu en faisant b copies du niveau
       n-1 qui sont connectés comme un cycle par une arête, le niveau
       0 étant le graphe à un sommet. Il faut b>1. Lorsque n=2, on
       obtient un sorte de soleil.

....sierpinski n b
....
       Graphe de Sierpiński généralisé, le graphe classique, le
       triangle Sierpiński qui est planaire, est obtenu avec b=3. Il a
       ((b-2)*b^n+b)/(b-1) sommets et est défini de manière récursive.
       Le niveau n est obtenu en faisant b copies du niveau n-1 qui
       sont connectés comme un cycle, le niveau 1 étant un cycle de b
       sommets. Il faut b>2 et n>0. A la différence du graphe d'Hanoï,
       les arêtes du cycle sont contractées. Le graphe de Hajós est
       obtenu avec n=2 et b=3.

....moser
....
       Graphe "Moser spindle" découvert par les frères Moser. C'est un
       "unit distance graph" du plan (deux points sont adjacents s'ils
       sont à distance exactement 1) de nombre chromatique 4. Il est
       planaire et possède 7 sommets. On ne connaît pas d'unit
       distance graphe avec un nombre chromatique supérieur.

....markstrom
....
       Graphe de Markström. Il est cubique planaire à 24 sommets. Il
       n'a pas de cycle de longueur 4 et 8.

....robertson
....
       Graphe de Robertson. C'est le plus petit graphe 4-régulier de
       maille 5. Il a 19 sommets, est 3-coloriable et de diamètre 3.

....wiener-araya
....
       Graphe découvert en 2009 par Wiener & Araya. C'est le plus
       petit graphe hypo-Hamiltonien planaire connu, c'est-à-dire
       qu'il n'est pas Hamiltonien mais la suppression de n'importe
       quel sommet le rend Hamiltonien. Il possède 42 sommets, 67
       arêtes, et est de diamètre 7.

....zamfirescu
....
       Graphe de Zamfirescu à 48 sommets découvert en 2007. Il est
       planaire et hypo-Hamiltonien. C'est le second plus petit (voir
       wiener-araya). Il possède 76 arêtes et a un diamètre de 7.

....hatzel
....
       Graphe de Hatzel. Il est planaire, de diamètre 8, possède 57
       sommets et 88 arêtes, et est hypo-Hamiltonien (voir
       wiener-araya). C'était le plus petit planaire hypo-Hamiltonien
       connu avant le graphe de Zamfirescu.

....clebsch n
....
       Graphe de Clebsch d'ordre n. Il est construit à partir d'un
       hypercube de dimension n en ajoutant des arêtes entre les
       sommets opposés, c'est-à-dire à distance n. Le graphe classique
       de Clebsch est réalisé pour n=4 dont le diamètre est deux.

....gear n
....
       Graphe planaire à 2n+1 sommets composé d'une roue à n rayons
       (graphe "wheel 2n") et de n sommets chacun connecté deux
       sommets voisins du bord de la roue. Il est construit à partir
       du graphe "cage 2n 2 2 0" auquel on ajoute un sommet central.

....flower_snark n
....
       Graphe cubique à 4n sommets construit de la manière suivante:
       1) on part de n étoiles disjointes à 3 feuilles, la i-ème ayant
       pour feuilles les sommets notés u_i,v_i,w_i, i=1..n; 2) pour
       chaque x=u,v ou w, x_1-...-x_n induit un chemin; et 2) sont
       adjacents: u_0-u_n, v_0-w_n et w_0-v_n. Pour n>1, ces graphes
       sont non planaires, non Hamiltoniens, 3-coloriables et de
       maille au plus 6.

....udg n r
....
       Graphe géométrique aléatoire (random geometric graph) sur n
       points du carré [0,1[ × [0,1[. Deux sommets sont adjacents si
       leurs points sont à distance ≤ r.  Il s'agit de la distance
       selon la norme L_2 (par défaut), mais cela peut être changée
       par l'option -norm. Le graphe devient connexe avec grande
       probabilité lorsque r=rc ~ √(ln(n)/n). Si r<0, alors le rayon
       est initialisé à rc. Un UDG (unit disk graph) est normalement
       un graphe d'intersection de disques fermés de rayon 1.

....gabriel n
....
       Graphe de Gabriel. Graphe géométrique défini à partir d'un
       ensemble de n points du carré [0,1[ × [0,1[.  Les points i et j
       sont adjacents ssi le plus petit disque (voir -norm) passant
       par i et j ne contient aucun autre point. Ce graphe est
       planaire et connexe. C'est un sous-graphe du graphe de
       Delaunay. Son étirement est non borné.

....rng n
....
       Graphe du proche voisinage (Relative Neighborhood
       Graph). Graphe géométrique défini à partir d'un ensemble de n
       points du carré [0,1[ × [0,1[.  Les points i et j sont
       adjacents ssi il n'existe aucun point k tel que
       max{d(k,i),d(k,j)} < d(i,j) où d est la distance (L_2 par
       défaut, voir -norm).  Dit autrement, la "lune" définie par i et
       j doit être vide. Ce graphe est planaire et connexe. C'est un
       sous-graphe du graphe de Gabriel.

....nng n
....
       Graphe du plus proche voisin (Nearest Neighbor Graph). Graphe
       géométrique défini à partir d'un ensemble de n points du carré
       [0,1[ × [0,1[.  Le point i est connecté au plus proche autre
       point (par défaut selon la norme L_2, voir -norm). Ce graphe
       est une forêt couvrante du graphe rng de degré au plus 6 (si la
       norme est L_2).

....thetagone n p k v
....
       Graphe géométrique défini à partir d'un ensemble de n points du
       carré [0,1[ × [0,1[. En général le graphe est planaire et
       connexe avec des faces internes de longueur au plus p (pour k
       diviseur de p et v=1). On peut interpréter les paramètres comme
       suit: p≥3 est le nombre de cotés d'un polygone régulier, k≥1 le
       nombre d'axes (ou de direction), et v∈[0,1] le cône de
       visibilité. Toute valeur de p<3 est interprétée comme une
       valeur infinie, et le polygone régulier correspondant
       interprété comme un cercle. L'adjacence entre une paire de
       sommets est déterminée en temps O(kn).
....
       Plus formellement, pour tout point u et v, et entier i, on note
       P_i(u,v) le plus petit p-gone (polygone régulier à p cotés)
       passant par u et v dont u est un sommet, et dont le vecteur
       allant de u vers son centre forme un angle de i*2π/k avec l'axe
       des abscisses, intersecté avec un cône de sommet u et d'angle
       v*(p-2)*π/p (v*π si p est infini) et dont la bissectrice passe
       par le centre du p-gone. Alors, u est voisin de v s'il un
       existe au moins un entier i dans [0,k[ tel que l'intérieur de
       P_i(u,v) est vide. La distance entre u et le centre du p-gone
       définit alors une distance (non symétrique) de u à v.
....
       Si v=1 (visibilité maximale), P_i est précisément un p-gone. Si
       v=0 (visibilité minimale), P_i est l'axe d'angle i*2π/k pour un
       entier i. Si v=.5, P_i est un cône formant un angle égale à 50%
       de l'angle défini par deux cotés consécutifs du p-gone, angle
       vallant (p-2)π/p. Si v=2p/((p-2)k) (ou simplement 2/k si p est
       infini) alors la visibilité correspond à un cône d'angle 2π/k,
       l'angle entre deux axes. Comme il faut v≤1, cela implique que
       k≥2p/(p-2) (k≥2 si p infini). On retrouve le Theta_k-Graph pour
       chaque k≥6 en prenant p=3 et v=6/k, le demi-Theta-Graph pour
       tout k≥3 en prenant p=3 et v=3/k, le Yao_k-Graph pour chaque
       k≥2 en prenant p=0 (infini) et v=2/k, et la triangulation de
       Delaunay si p=0 (infini), k très grand et v=1. En fait, ce
       n'est pas tout-à-fait le graph Yao_k, pour cela il faudrait que
       u soit le centre du polygone (c'est-à-dire du cercle).

....pat p q r
....
       Graphe possédant pqr sommets, issu d'un jeu à un joueur proposé
       par Pat Morin (Barbade, mars 2016). Le jeu se déroule sur une
       grille p×q et comprend r coups. Un coup est un ensemble de
       positions de la grille strictement croissantes (coordonnées en
       x et en y strictement croissantes). De plus, si la position
       (x,y) est jouée alors toutes les positions situées sur la même
       ligne mais avec une abscisse au moins x ou sur la même colonne
       mais avec une ordonnées au moins y sont interdites pour tous
       les coups suivants. Le score est le nombre total de positions
       jouées en r coups. Il s'agit de trouver le score maximum.
       Lorsque r=1, le score maximum vaut min(p,q). Lorsque p=q=n et
       r=2, alors le score maximum vaut ⎣ 4n/3⎦. La question est
       ouverte lorsque r>2, c'est au moins n^1.516 pour r=n où la
       constante vaut log_9(28).
....
       Les sommets du graphes sont les positions dans les r grilles
       p×q et deux sommets sont adjacents les positions sont en
       conflits. Le score du jeu est alors un ensemble indépendant du
       graphe. Si r=1, le graphe est une grille p×q. Ce graphe active
       l'option -pos car un dessin de ce graphe (sous forme de
       grilles) est proposé.
....
       Ex: gengraph pat 4 4 4 -check kindepsat 8 | ./glucose -model

....rplg n t
....
       Random Power-Law Graph. Graphe à n sommets où les degrés des
       sommets suivent une loi de puissance d'exposant t. L'espérance
       du degré du sommet i=0...n-1 est w_i=(n/(i+1))^(1/(t-1)). Il
       s'agit d'une distribution particulière d'un Fixed Degree Random
       Graph où la probabilité d'avoir l'arête i-j est
       min{w_i*w_j/S,1} avec S=∑_k w_k. En général il faut prendre le
       paramètre t comme un réel de ]2,3[, la valeur communément
       observée pour le réseau Internet étant t=2.1.

....bdrg n_1 d_1 ... n_k d_k -1
....
       Bounded Degree Random Graph. Graphe aléatoire dont la
       distribution des degrés des sommets est fixée par les paires
       (n_i,d_i) signifiant qu'il y a n_i sommets de degré au plus
       d_i. Ainsi "bdrg n 3 -1" génère un graphe sous-cubic aléatoire
       à n sommets, si n est pair. Les sommets sont dupliqués selon
       leur distribution de degré puis un couplage aléatoire détermine
       les arêtes. Les boucles et les arêtes multiples sont supprimer.
       Il suit que le degré des sommets de dépasse pas d_i. Ils
       peuvent cependant être inférieurs. Le nombre de sommets est
       n=∑_i n_i et le nombre d'arêtes au plus m=(1/2)*∑_i (n_i*d_i).
       Si cette somme n'est pas entière, alors le degré d'un des
       sommets ayant d_i>0 est diminué d'un. (C'est un sommet le degré
       d_i>0 avec le plus grand i qui est choisi.)

....fdrg n_1 d_1 ... n_k d_k -1
....
       Fixed Degree Random Graph. Graphe aléatoire asymptotiquement
       uniforme dont les degrés des sommets sont fixées par les paires
       (n_i,d_i) signifiant qu'il y a n_i sommets de degré d_i. La
       suite des degrés doit être graphique, à savoir qu'il existe au
       moins un graphe simple ayant ces degrés (sinon une erreur est
       affichée). Ainsi "fdrg n 3" génère un graphe cubic aléatoire
       asymptotiquement uniforme, à condition que n soit pair. [CE
       GRAPHE N'EST PAS ENCORE COMPLETEMENT IMPLEMENTE]

....matching n
....
       Graphe composé de n arêtes indépendantes, c'est-à-dire de n
       copies de K_2.

....load file[:range]
....loadc file[:range]
....
       Graphe défini à partir du fichier "file" ou de l'entrée
       standard si file vaut "-". Si "file" est une famille de
       graphes, alors il est possible d'utiliser la variante
       "file:range" pour préciser l'identifiant du graphe souhaité
       (sinon c'est le premier graphe de la famille qui sera
       considéré). Le graphe (ou la famille) doit être au format
       standard, les sommets numérotés par des entiers positifs. Les
       caractères situés sur une ligne après "//" sont ignorés, ce qui
       permet de mettre des commentaires.
....
       Le temps et l'espace nécessaire au chargement du graphe sont
       linéaires en la taille du fichier (si "file" est une famille de
       graphes, le fichier est entièrement lu).  Cependant, pour la
       génération à proprement parlée du graphe final, qui peut
       comprendre l'option -not par exemple, toutes les arêtes
       potentielles, soit O(n²), sont passées en revue pour être
       testées. La variante "loadc" (pour "load & check") permet une
       génération plus rapide lorsqu'utilisée avec -check (ou les
       alias utilisant -check, comme -maincc par exemple). Elle permet
       de passer directement de l'étape de chargement du graphe à
       l'étape du test de l'algorithme en sautant la phase de
       génération des arêtes. En contre-partie, le graphe n'est pas
       affiché et les options comme -not, -permute, -delv, -dele,
       etc. n'ont plus d'effet. La variante "loadc file" est environ
       20% plus rapide que "load file -fast".
....
       Pour charger un graphe au format dot on peut utiliser le script
       dot2gen.awk en amont, comme dans l'exemple suivant:
....
       nop file.dot | awk -f dot2gen.awk | ./gengraph load -
....
       Le filtre nop de GraphViz, qui est recommandé mais pas
       nécessaire, permet de standardiser le format dot initial. Il
       transforme par exemple les expressions du type "a--{b;c;}" en
       "a--b;a--c;".
....
       Notez que la suite d'options "load file -fast -format dot<xxx>"
       permet de convertir "file" au format <xxx> souhaité. Ce graphe
       active l'option -directed si "file" contient au moins un
       arc. Dans ce cas l'option -undirected n'aura pas d'effet.


   GRAPHES ORIENTES :

....aqua n c_1 ... c_n
....

       Graphe orienté dont les sommets sont les suites de n entiers
       positifs dont la somme fait c_1 et dont le i-ème élément est au
       plus c_i. Ils représentent les façons de répartir une quantité
       c_1 de liquide dans n récipients de capacité c_1 ... c_n. Il y
       a un arc u->v s'ils existent i et j tels que v est le résultat
       du versement du récipient c_i vers le récipiant c_j.
....       
       Ex: gengraph aqua 3 3 2 1 -label 1 -directed -dotfilter dot -visu


   GRAPHES COMPOSES :

....mesh p q (= grid 2 p q)
....
       Grille 2D de p x q sommets.

....hypercube d (= grid d 2 ... 2)
....
       Hypercube de dimension d.

....path n (= grid 1 n)
....
       Chemin à n sommets.

....cycle n (= ring n 1 1)
....
       Cycle à n sommets.

....torus p q (= grid 2 -p -q)
....
       Tore à p x q sommets.

....stable n (= ring n 0)
....empty n (= ring n 0)
....
       Stable à n sommets.

....clique n (= -not ring n 0)
....
       Graphe complet à n sommets.

....bipartite p q (= rpartite 2 p q)
....
       Graphe biparti complet K_{p,q}.

....utility (= rpartite 2 3 3)
....
       Graphe biparti complet K_{3,3} qui doit son nom au problème de
       la connexion planaire de trois maisons à trois stations (eau,
       gaz, électricité).

....octahedron (= antiprism 3)
....
       Octaèdre: graphe 4-régulier planaire à 6 sommets ayant 8 faces
       triangulaires. Il s'agit de deux pyramides dont la base à 4
       sommets est commune. C'est aussi le graphe de Johnson J(4,2).

....d-octahedron d (= -not matching d)
....
       Octaèdre de dimension d: obtenu à partir d'un octaèdre de
       dimension d-1 auquel on ajoute deux sommets universels,
       l'octaèdre de dimension 1 étant composé d'un stable de deux
       sommets.  L'octaèdre classique est obtenu avec d=3, pour d=2 il
       s'agit d'un carré.

....tetrahedron (= -not ring 4 0)
....
       Tétraèdre: pyramide composée de 4 faces triangulaires. C'est
       aussi une clique à 4 sommets.

....hexahedron (= cube)
....
       Hexaèdre: cube composé de 6 faces carrées.

....associahedron (= flip 6)
....
       Associaèdre (3D): graphe planaire cubique à 14 sommets composé
       de 3 faces carrées et 6 faces pentagonales.

....johnson n k (= -not kneser n k k-2)
....
       Graphe de Johnson J(n,k). Les sommets sont tous les
       sous-ensembles à k éléments de [0,n[ (il faut donc 0 ≤ k ≤
       n). Deux sommets sont adjacents ssi leurs ensembles
       correspondant ont k-1 éléments en commun. La distance entre
       deux sommets est la distance de Hamming entre les ensembles
       correspondant. Ils sont réguliers de degré k(n-k), de diamètre
       min{k,n-k}, de sommet-connectivité k(n-k) et aussi distance
       régulier. J(n,1) est la clique K_n, J(n,2) est le complément du
       graphe de Kneser K(n,2) et le line-graphe de K_n. En fait, tout
       sous-graphe induit de J(n,2) est un line-graphe. J(4,2) est
       l'octaèdre, J(5,2) le complément du graphe de Petersen.

....turan n r (= rpartite r ⎣ n/r⎦ ... ⎣ n/r⎦)
....
       Graphe de Turán à n sommets. Il s'agit qu'un graphe r-parti
       complet de n sommets avec n mod r parts de ⎡ n/r⎤ sommets et
       r-(n mod r) parties de ⎣ n/r⎦ sommets. Il est régulier lorsque r
       divise n. Le nombre d'arêtes est ⎣ (r-1)n^2/r⎦.  C'est le graphe
       sans clique de taille r+1 ayant le plus grand nombre d'arêtes.
       Lorsque n=2r, il s'agit du "cocktail party graph". Lorsque n=6
       et r=3, c'est l'octaèdre.

....claw (= rpartite 2 1 3)
....
       Graphe biparti complet K_{1,3}.

....star n (= rpartite 2 1 n)
....
       Arbre (étoile) à n feuilles et de profondeur 1.

....tree n (= arboricity n 1)
....
       Arbre plan enraciné aléatoire uniforme à n sommets. Les sommets
       sont numérotés selon un parcours en profondeur depuis la racine
       et le long de la face extérieure.

....caterpillar n (= grid 1 n-r -star r, r=random()%n)
....
       Arbre à n sommets dont les sommets internes (de degré > 1)
       induisent un chemin. Il est obtenu à partir d'un chemin de
       longueur n-r (où r est un nombre aléatoire entre 0 et n-1) et
       en appliquant l'option -star r. Si l'option -seed est présente,
       (pour intervenir sur la valeur "r"), il est important qu'elle
       figure avant caterpillar.

....outerplanar n (= kpage n 1)
....
       Graphe planaire-extérieur aléatoire connexe à n sommets (plan
       et enraciné). Ils sont en bijection avec les arbres plans
       enracinés dont tous les sommets, sauf ceux de la dernière
       branche, sont bicoloriés. Les sommets sont numérotés le long de
       la face extérieure. C'est aussi une numérotation selon un
       parcours en profondeur depuis la racine de l'arbre bicolorié.
       Il est aussi possible de générer des graphes
       planaires-extérieurs aléatoires Hamiltoniens, donc 2-connexes,
       avec "planar n f -1" ou "polygon n".

....squaregraph n (= planar n 4 4)
....
       Squaregraph aléatoire à n faces. Ce sont des graphes planaires
       2-connexes dont toutes les faces (sauf l'extérieure) sont des
       carrées. De plus, les sommets des faces internes sont de degré
       au moins 4. Ce sont des 2-pages comme tous les sous-graphes de
       quadrangulations.

....random n p (= -not ring n 0 -dele 1-p)
....
       Graphe aléatoire à n sommets et dont la probabilité d'avoir une
       arête entre chaque paire de sommets est p. L'option -dele étant
       déjà présente, il n'est pas conseillé de la réutiliser pour ce
       graphe.

....sunlet n (= grid 1 -n -star -1)
....
       Cycle à n sommets avec un sommet pendant à chacun d'eux. Un
       sunlet 3 est parfois appelé netgraph.

....netgraph (= grid 1 -3 -star -1)
....
       Graphe à 6 sommets composé d'un triangle avec un sommet pendant
       à chacun d'eux. C'est le complémentaire du graphe de Hajós.

....sunflower n (= cage 2n 2 2 0)
....
       Tournesol à n pétales. C'est un graphe planaire-extérieur à 2n
       sommets composé d'un cycle de longueur n≥3 où chaque arête
       partage le coté d'un triangle. C'est le graphe "gear n" sans le
       sommet central. C'est le graphe de Hajós pour n=3.

....gem (= fan 4 1)
....
       Graphe à 5 sommets composé d'un chemin et d'un sommet universel.

....egraph (= grid 1 3 -star -1)
....
       Arbre à 6 sommets et 3 feuilles en forme de E.

....knight p q (= chess p q 1 2)
....
       Graphe des déplacements possible du chevalier dans un échiquier
       p q.

....antelope p q (= chess p q 3 4)
....
       Graphe des déplacements possible d'une antilope dans un
       échiquier p q, une antilope étant une pièce hypothétique se
       déplaçant de 3 cases selon un axe et de 4 selon l'autre.

....camel p q (= chess p q 1 3)
....
       Graphe des déplacements possible d'un chameau dans un échiquier
       p q, un chameau étant une pièce hypothétique se déplaçant de 1
       case selon un axe et 3 de selon l'autre.

....giraffe p q (= chess p q 1 4)
....
       Graphe des déplacements possible d'une giraffe dans un
       échiquier p q, une giraffe étant une pièce hypothétique se
       déplaçant de 1 case selon un axe et de 4 selon l'autre.

....zebra p q (= chess p q 2 3)
....
       Graphe des déplacements possible d'un zébre dans un échiquier p
       q, un zébre étant une pièce hypothétique se déplaçant de 2
       cases selon un axe et de 3 selon l'autre.

....petersen (= kneser 5 2 0)
....
       Graphe de Kneser particulier. Il est cubique et possède 10
       sommets. Il n'est pas Hamiltonien et c'est le plus petit graphe
       dont le nombre de croisements (crossing number) est 2. C'est le
       complément du line-graphe de K_5.

....tietze (= flower_snark 3)
....
       Graphe de Tietze. Il est cubique avec 12 sommets. Il possède un
       chemin Hamiltonien, mais pas de cycle. Il peut être plongé sur
       un ruban de Möbius, a un diamètre et une maille de 3. Il peut
       être obtenu à partir du graphe de Petersen en appliquant une
       opération Y-Delta.

....mobius-kantor (= gpetersen 8 2)
....
       Graphe de Möbius-Kantor. Graphe cubique à 16 sommets de genre
       1. Il est Hamiltonien, de diamètre 4 et de maille 6.

....dodecahedron (= gpetersen 10 2)
....
       Dodécaèdre: graphe planaire cubique à 20 sommets. Il possède 30
       arêtes et 12 faces qui sont des pentagones. C'est le dual de
       l'icosaèdre.

....desargues (= gpetersen 10 3)
....
       Graphe de Desargues. Il est cubique à 20 sommets. Il est
       Hamiltonien, de diamètre 5 et de maille 6.

....durer (= gpetersen 6 2)
....
       Graphe de Dürer. Graphe cubique planaire à 12 sommets de
       diamètre 4 et de maille 3. Il peut être vu comme un cube avec
       deux sommets opposés tronqués (remplacés par un cycle de
       longueur 3).

....prism n (= gpetersen n 1)
....
       Prisme, c'est-à-dire le produit cartésien d'un cycle à n
       sommets et d'un chemin à deux sommets. Pour n=3, c'est un
       graphe de Halin et aussi le complémentaire d'un cycle de
       longueur 6, et pour n=4 il s'agit du cube.

....cylinder p q (= grid p -q)
....
       Produit cartésien d'un chemin à p sommets et d'un cycle à q
       sommets. Cela généralise le prisme (prism n = cylinder n 3). Un
       cube est un "cylinder 2 4".

....nauru (= pstar 4)
....
       Graphe de Nauru. C'est un graphe cubique à 24 sommets. Il
       s'agit d'un graphe "permutation star" de dimension 4. C'est
       aussi un graphe de Petersen généralisé P(12,5).

....headwood (= cage 14 2 5 -5)
....
       Graphe de Headwood. C'est un graphe cubique à 14 sommets, de
       maille 6 et de diamètre 3. C'est le plus petit graphe dont le
       nombre de croisements (crossing number) est 3.

....franklin (= cage 12 2 5 -5)
....
       Graphe de Franklin. C'est un graphe cubique à 12 sommets, de
       maille 4 et de diamètre 3.

....dyck (= cage 32 4 5 0 13 -13)
....
       Graphe de Dyck. C'est un graphe cubique à 32 sommets, le seul à
       être symétrique. Il est aussi torique, c'est-à-dire de genre 1.

....pappus (= cage 18 6 5 7 -7 7 -7 5)
....
       Graphe de Pappus. C'est un graphe cubique à 18 sommets, de
       maille 6 et de diamètre 4.

....mcgee (= cage 24 3 12 7 -7)
....
       Graphe de McGee. C'est un graphe cubique à 24 sommets, de
       maille 7 et de diamètre 4.

....tutte-coexter (= cage 30 6 -7 9 13 -13 -9 7)
....
       Graphe de Tutte-Coexter appelé aussi 8-cage de Tutte. C'est un
       graphe cubique à 30 sommets, de maille 8 et de diamètre 4.
       C'est un graphe de Levi mais surtout un graphe de Moore,
       c'est-à-dire un graphe d-régulier de diamètre k dont le nombre
       de sommets est 1+d*S(d,k) (d impair) ou 2*S(d,k) (d pair) avec
       S(d,k)=∑_{i=0}^(k-1) (d-1)^i.

....gray (= cage 54 6 7 -7 25 -25 13 -13)
....
       Graphe de Gray. C'est un graphe cubique à 54 sommets qui peut
       être vu comme le graphe d'incidence entre les sommets d'une
       grille 3x3x3 et les 27 lignes droites de la grille. Il est
       Hamiltonien, de diamètre 6, de maille 8, et de genre 7. Il est
       arête-transitif et régulier sans être sommet-transitif.

....grotzsch (= mycielski 4)
....
       Graphe de Grötzsch. C'est le plus petit graphe sans triangle de
       nombre chromatique 4. Il possède 11 sommets et 20 arêtes. Comme
       le graphe de Chvátal, il est non-planaire de diamètre 2, de
       maille 4 et Hamiltonien. C'est le graphe de Mycielskian du
       cycle à 5 sommets.

....hajos (= sierpinski 2 3)
....
       Graphe de Hajós. Il est composé de trois triangles deux à deux
       partageant un sommet distinct. On peut le dessiner comme un
       triangle dans un triangle plus grand. Il est planaire et
       possède 6 sommets. C'est un graphe de Sierpinski ou encore le
       complémentaire d'un "sunlet 3", complémentaire du "netgraph",
       un "sunflower 3" ou encore "cage 6 2 2 0".

....house (= -not grid 1 5)
....
       Graphe planaire à 5 sommets en forme de maison. C'est le
       complémentaire d'un chemin à 5 sommets.

....wagner (= ring 8 2 1 4)
....
       Graphe de Wagner appelé aussi graphe W_8, un cycle à 8 sommets
       où les sommets antipodaux sont adjacents. C'est un graphe
       cubique à 8 sommets qui n'est pas planaire mais sans K_5. C'est
       aussi une échelle de Möbius.

....mobius n (= ring n 2 1 n/2)
....
       Échelle de Möbius, graphe cubic à n sommets obtenu à partir
       d'un cycle à n sommets dont les sommets opposés sont
       ajdacents. Lorsque n est pair, il s'agit d'un ruban de Möbius,
       c'est-à-dire d'une échelle dont le premier et dernier barreau
       sont recollés en sens opposé. Pour n≤5, il s'agit d'une clique
       à n sommets. Lorsque n≥5, le graphe n'est plus planaire, et
       pour n=8, il s'agit du graphe de Wagner.

....ladder n (= grid 2 2 n)
....
       Graphe échelle à n barreaux, soit une grille à 2 x n sommets.

....cube (= crown 4)
....
       Hypercube de dimension 3, graphe planaire cubic à 8 sommets où
       toutes les faces sont des rectangles.

....diamond (= fan 2 2)
....
       Clique à quatre sommets moins une arête. C'est un graphe
       allumette, c'est-à-dire planaire et distance unitaire.

....gosset (= ggosset 8 2 2 3 6 -1)
....
       Graphe de Gosset. Il est 27-régulier avec 56 sommets et 756
       arêtes, de diamètre, de rayon et de maille 3. Il est
       27-arête-connexe et 27-sommet-connexe. C'est localement un
       graphe de Schläfli, c'est-à-dire que pour tout sommet le
       sous-graphe induit par ses voisins est isomorphe au graphe de
       Schläfli, qui lui-même localement un graphe de Clebsch.

....wheel n (=ringarytree 1 0 n 2)
....
       Roue à n+1 sommets, graphe planaire composé d'un cycle à n
       sommets et d'un sommet universel, donc connecté à tous les
       autres.

....web n r (=ringarytree r 1 n 2)
....
       Graphe planaire à 1+n*r sommets composé d'une étoile à n
       branches de longueur r, les sommets de même niveau étant
       connectés par un cycle. Il généralise "wheel n" (r=1).

....ygraph (= ringarytree 2 1 3 0)
....
       Arbre à 7 sommets composé d'une étoile à trois branches.

....binary h (= ringarytree h 2 2 0)
....
       Arbre binaire complet de profondeur h. Il possède 2^(h+1)-1
       sommets et la racine est de degré deux.

....arytree h k r (= ringarytree h k r 0)
....
       Arbre complet de hauteur h où chaque noeud interne à exactement
       k fils, le degré de la racine étant de degré r.

....rbinary n (= rarytree n 2 0)
....rbinaryz n (= rarytree n 2 1)
....
       Arbre binaire plan aléatoire uniforme à n noeuds internes. Il
       possède 2n-1 sommets (2n pour la variante rbinaryz) numérotés
       selon un parcours en profondeur modifié: tous les fils du
       sommet courant sont numérotés avant l'étape de récursivité. La
       racine est de degré 2 (=rbinary) ou 1 (=rbinaryz). Le dessin
       avec dot (-visu) ne respecte pas le plongement de l'arbre.

....tw n k (= ktree n k -dele .5)
....
       Graphe de largeur arborescente au plus k aléatoire à n
       sommets. Il s'agit d'un k-arbre partiel aléatoire dont la
       probabilité d'avoir une arête est 1/2. L'option -dele étant
       déjà présente, il n'est pas conseillé de la réutiliser pour ce
       graphe.

....pw n k (= kpath n k -dele .5)
....
       Graphe de pathwidth au plus k, aléatoire et avec n sommets.

....tadpole n p (= barbell -n 1 p)
....dragon n p (= barbell -n 1 p)
....
       Graphe à n+p sommets composé d'un cycle à n sommets relié à un
       chemin à p sommets.

....lollipop n p (= barbell n p 0)
....
       Graphe "tapette à mouches" (Lollipop Graph) composé d'une
       clique à sommets relié à un chemin de longueur p. Il a n+p
       sommets.

....pan n (= barbell -n 1 1)
....
       Graphe à n+1 sommets composé d'un cycle à n sommets et d'un
       seul sommet pendant.

....banner (= barbell -4 1 1)
....
       Graphe à 5 sommets composé d'un carré et d'un sommet pendant.

....paw (= barbell -3 1 1)
....
       Graphe à 4 sommets composé d'un triangle et d'un sommet
       pendant.

....theta0 (=barbell -5 -5 -2)
....
       Graphe Theta_0. C'est un graphe à 7 sommets série-parallèle
       obtenu à partir d'un cycle de longueur 6 et en connectant deux
       sommets antipodaux par un chemin de longueur 2. C'est un graphe
       allumette, c'est-à-dire planaire et distance unitaire.

....td-delaunay n (= thetagone n 3 3 1)
....
       Triangulation de Delaunay utilisant la distance triangulaire
       (TD=Triangular Distance). Il s'agit d'un graphe planaire défini
       à partir d'un ensemble de n points aléatoires du carré [0,1[ ×
       [0,1[. Ce graphe a un étirement de 2 par rapport à la distance
       euclidienne entre deux sommets du graphe. Ce graphe, introduit
       par Chew en 1986, est le même que le graphe "demi-theta_6", qui
       est un "theta-graph" utilisant 3 des 6 cônes. La dissymétrie
       qui peut apparaître entre le bord droit et gauche du dessin est
       lié au fait que chaque sommet n'a qu'un seul axe dirigé vers la
       droite, alors qu'il y en a deux vers la gauche.

....theta n k (= thetagone n 3 k 6/k)
....
       Theta-graphe à k secteurs réguliers défini à partir d'un
       ensemble de n points du carré [0,1[ × [0,1[. Les sommets u et v
       sont adjacents si le projeté de v sur la bissectrice de son
       secteur est le sommet le plus proche de u. Ce graphe n'est pas
       planaire en général, mais c'est un spanner du graphe complet
       euclidien. Le résultat est valide seulement si k≥6.

....dtheta n k (= thetagone n 3 k/2 6/k)
....
       Demi-Theta-graphe à k secteurs réguliers défini à partir d'un
       ensemble de n points du carré [0,1[ × [0,1[. La définition est
       similaire au Theta-graphe excepté que seul 1 secteur sur 2 est
       considéré. Il faut k≥2 pair, mais le résultat n'est valide que
       si k≥6. Pour k=2 ou 3, il s'agit d'un arbre, pour k=4 ou 5, le
       graphe est de faible tree-width. Pour k=6, ce graphe coïncide
       avec le graphe td-delaunay.
....
       Ex1: gengraph dtheta 500 6 -visu
       Ex2: gengraph dtheta 500 4 -pos 0 -visu
       Ex2: gengraph dtheta 500 2 -pos 0 -visu

....yao n k (= thetagone n 0 k 2/k)
....
       Graphe de Yao à k secteurs réguliers défini à partir d'un
       ensemble de n points du carré [0,1[ × [0,1[. Les sommets u et v
       sont adjacents si v est le sommet le plus proche de u (selon la
       distance euclidienne) de son secteur. Ce graphe n'est pas
       planaire en général, mais c'est un spanner du graphe complet
       euclidien. Le résultat est valide seulement si k≥2. En fait,
       ce n'est pas tout à fait le graphe de Yao (voir thetagone).

....percolation a b p (= udg a*b 1 -norm 1 -xy mesh a b -dele 1-p)
....
       Grille de percolation à coordonnées entières (i,j) de
       [0,a[x[0,b[ où p représente la probabilité d'existence de
       chaque arête. La différence avec le graphe "mesh a b -dele 1-p"
       est qu'ici le graphe est géométrique, donc dessiné selon une
       grille si l'option -visu est ajouté.

....point n (= ring n 0 -pos 1)
....
       Graphe géométrique composé de n points du plan mais sans aucune
       arête (voir "stable"). Ce graphe permet de visualiser la
       distribution des points, par défaut uniforme sur [0,1[ × [0,1[,
       mais qui peut être modifiée avec l'option -xy.
....
       Ex1: gengraph point 500 -xy seed 3 2.1 -visu
       Ex2: gengraph point 1000 -xy seed 3 -0.1 -visu

....regular n d (= fdrg n d -1)
....
       Graphe d-régulier aléatoire à n sommets asymptotiquement
       uniforme. Il faut que n*d soit pair.

....cubic n (= fdrg n 3 -1)
....
       Graphe cubic aléatoire à n sommets asymptotiquement
       uniforme. Il faut que n soit pair.

....plrg n t (= bdrg n_1 d_1 ... n_k d_k -1)
....
       Power-Law Random Graph. Graphe aléatoire à n sommets dont la
       distribution des degrés suit une loi en puissance d'exposant
       t>0, la probabilité qu'un sommet soit de degré i>0 étant
       proportionnelle à 1/i^t. Plus précisément, la distribution est
       la suivante:
....
         - d_1=1, n_1=⎣ exp(𝛼)⎦ + n-s(𝛼)
         - d_i=i, n_i=⎣ exp(𝛼)/i^t⎦ pour 2≤i≤p(𝛼)
....
       où a est un réel minimisant |n-s(𝛼)| avec p(𝛼)=⎣ exp(𝛼/t)⎦
       et s(𝛼)=∑_{i=1}^p(𝛼) ⎣ exp(𝛼)/i^t⎦. Ce sont les mêmes
       graphes que ceux générés par Brady-Cowen'06.

.....
HISTORIQUE

       v1.2 octobre 2007:
            - première version

       v1.3 octobre 2008:
            - options: -shift, -width
            - correction d'un bug pour les graphes de permutation
	    - accélération du test d'ajacence pour les arbres, de O(n) à O(1),
              grâce à la représentation implicite
	    - nouveau graphes: outerplanar, sat

       v1.4 novembre 2008:
            - format de sortie: matrix, smatrix, list
            - nouveau graphe: kout
            - correction d'un bug dans l'option -width
	    - correction d'un bug dans la combinaison -format/shift/delv

       v1.5 décembre 2008:
            - correction d'un bug dans tree lorsque n=1

       v1.6 décembre 2009:
            - nouveaux graphes: rpartite, bipartite

       v1.7 janvier 2010:
            - nouveaux graphes: icosa, dodeca, rdodeca, cubocta, geo,
	      wheel, cage, headwood, pappus, mcgee, levi, butterfly,
	      hexagon, whexagone, arytree, binary, ktree, tw, kpath,
	      pw, arboricity, wagner, mobius, tutte-coexter, paley
            - nouveau format de sortie: -format dot
	    - nouvelles options: -header, -h, -redirect, -dotpdf
            - correction d'un bug dans kout, et dans tree lorsque n=0
	    - tree devient un cas particulier d'arboricity.
	    - aide en ligne pour les paramètres des graphes.

       v1.8 juillet 2010:
            - nouveaux graphes: chvatal, grotzsch, debruijn, kautz
	      gpstar, pstar, pancake, nauru, star, udg, gpetersen,
              mobius-kantor, desargues, durer, prism, franklin,
	      gabriel, thetagone, td-delaunay, yao, theta, dtheta
            - suppression du graphe geo (remplacé par udg)
            - nouvelles options: -pos, -norm, -label, -dotfilter
	    - nouvelle famille d'options: -xy file/noise/scale/seed
	    - définition plus compacte dodeca (non explicite)
	    - utilisation du générateur random() plutôt que rand().
	    - correction d'un bug dans "-format standard" qui provoquait une erreur.
	    - correction d'un bug dans kneser pour k=0, n=0 ou k>n/2.
	    - nouveaux formats: -format dot<xxx>, -format xy
	    - suppression de -dotpdf (qui est maintenant: -format dotpdf)
	    - labeling pour: gpetersen, gpstar, pstar, pancake, interval,
	      permutation

       v1.9 août 2010:
            - renome -h en -list
	    - renome -xy file en -xy load
	    - centrage des positions sur le barycentre des graines (-xy seed)
	    - nouvelles options: -star, -visu, -xy round
	    - les graphes peuvent être stockés en mémoire, sous la forme d'une liste
	      d'adjacence grâce à l'option -check.
	    - généralisation de -delv p avec p<0
	    - nouveaux graphes: caterpillar, hajos, hanoi, sierpinski
	      sunlet, load file
	    - labeling pour: hanoi, sierpinski
	    - aide sur toutes les options (nécessitant au moins un paramètre)
              et non plus seulement pour les graphes
	    - nouvelle famille d'options: -vcolor deg/degr/pal
	    - correction d'un bug pour l'aide dans le cas de commande
	      préfixe (ex: pal & paley)

       v2.0 septembre 2010:
	    - nouvelles options: -vcolor degm/list/randg, -xy unique/permutation,
	      -check bfs, -algo iso/sub
	    - l'option -xy round p admet des valeurs négatives pour p.
	    - les options "load file" et "-xy load file" permettent la
              lecture à partir de l'entrée standard en mettant
              file="-", la lecture de famille de graphes, et supporte les commentaires.
	    - les formats list/matrix/smatrix utilisent un espace
	      linéaire O(n+m) contre O(n²) auparavant.
	    - les sommets sur le bord (graphes géométriques) ne sont plus coupés
	      (bounding-box (bb) plus grandes).
	    - nouveaux graphes: kpage, outerplanar n (=kpage n 1), rng, nng
	      fritsch, soifer, gray, hajos (qui avait été définit mais non
	      implémenté !), crown, moser, tietze, flower_snark, markstrom,
	      clebsch, robertson, kittell, rarytree, rbinary, poussin, errara
	    - les graphes de gabriel (et rng,nng) dépendent maintenant de -norm.
	    - "wheel n" a maintenant n+1 sommets, et non plus n.
	    - aide en ligne améliorée avec "?". Ex: gengraph tutte ? / -visu ?
	    - les options -help et ? permettent la recherche d'un mot clef.
	      Ex: gengraph -help planaire / ? arbre
	    - description plus compacte de tutte (et des graphes à partir d'un tableau)
	    - correction d'un bug pour rpartite (qui ne marchait pas)

       v2.1 octobre 2010:
	    - nouvelles options:
	      -check degenerate/gcolor/edge/dfs/ps1/paths/paths2/iso/sub/minor/isub
	      -filter minor/sub/iso/vertex/edge/degenerate/ps1
	      -filter degmax/degmin/deg/gcolor/component/radius/girth/diameter
	      -filter cut-vertex/biconnected/isub/all/minor-inv/isub-inv/sub-inv
            - suppression de -algo iso/sub: l'option -algo est réservée à la mis
	      au point de -check
	    - extension de -label b à b=2 qui force l'affiche des noms
              sous forme d'entiers même avec -permute.
	    - correction d'un bug pour house (qui ne marchait pas)
	    - nouveau graphe: windmill

       v2.2 novembre 2010:
            - gestion des graphes orientés: lecture d'un fichier de
              graphe (ou d'une famille avec arcs et arêtes)
	    - nouvelles options: -(un)directed, -(no)loop, -check twdeg/tw,
	      -filter tw/id/hyperbol/rename
	    - permet l'affichage de la "value" (p) dans l'option -filter
	    - nouveau graphe: aqua
	    - correction du graphe tutte-coexter et suppression du
              graphe levi (qui en fait était le graphe de tutte-coexter).
	    - généralisation de l'option "load" à load:id family

       v2.3 décembre 2010:
            - nouvelles options: -check ps1bis/edge, -filter ps1bis/tw2
	      -filter minus/minus-id/unique/connected/bipartite/forest
	      -check ps1ter
	    - remplacement de atof()/atoi() par strtod()/strtol() qui
	      sont plus standards.
	    - remplacement de LONG_MAX par RAND_MAX dans les
              expressions faisant intervenir random() qui est de type
              long mais qui est toujours dans [0,2^31[, même si
              sizeof(long)>4. Il y avait un bug pour les architectures
              avec sizeof(long)=8.
	    - nouveau graphe: cylinder
	    - suppression de la variante "load:id" au profit de la
              forme plus générale "file:range" valable pour load, -filter, etc.

       v2.4 janvier 2011:
            - correction d'un bug dans -filter minus-id
	    - correction d'un bug dans rpartite (incorrect à partir de r>5 parts)
	    - correction d'un bug dans whexagon (nb de sommets incorrects)
	    - nouvelles options: -check ps1x/girth, -filter ps1c/ps1x
	    - renomage: ps1bis -> ps1b, ps1ter -> ps1c
	    - nouveau graphe: mycielski
	    - la graphe grotzsch est maintenant défini à partir du graphe
	      mycielski (la définition précédante était fausse)
	    - bug détecté: td-delaunay 500 -check gcolor -format no -seed
              7 | grep '>6' qui donne jusqu'à 7 couleurs; le nb de
              couleurs affichées dans -check gcolor est erroné

       v2.5 mars 2011:
	    - nouveaux graphes: line-graph, claw

       v2.6 juin 2011:
	    - amélioration du test -filter ps1: détection de cliques et d'arbres

       v2.7 octobre 2011:
	    - nouvelle option: -check bellman (pour les géométriques seulement)
	    - ajoût des champs xpos,ypos à la structure "graph".
	    - nouveaux graphes: linial, linialc, cube, diamond, theta0,

       v2.8 novembre 2011:
	    - nouveaux graphes: ggosset, gosset, rplg, wiener-araya, headwood4
	    - correction d'un bug pour "-xy seed k n" lorsque k=1.
	    - nouvelles options: -check maincc, -maincc (non documentée)

       v2.9 février 2013:
	    - nouveau graphe: frucht, halin
	    - correction d'un bug pour "-check gcolor" qui ne
	      retrounait pas le nombre correct de couleurs, et qui de
	      plus n'utilisait pas l'heuristique du degré minimum.
	    - correction d'un bug pour "permutation -label 1"

       v3.0 octobre 2013:
	    - nouveaux graphes: rig, barbell, lollipop
	    - généralisation de l'option -filter forest
	    - nouvelles options: -apex, -filter isforest, -filter istree, -filter cycle
	    - correction d'un bug dans -filter vertex
	    - amélioration de l'aide lors d'erreurs de paramètre comme:
              "-filter F vertex" au lieu de "-filter F vertex n"
	    - amélioration de l'option -header

       v3.1 décembre 2013:
	    - nouveaux graphes: bpancake
	    - légère modification des labels des sommets des graphes pancake, 
	      gpstar et pstar
	    - nouvelles options: -xy grid, -xy vsize
	    - modification de la taille des sommets pour dot permettant de tenir
	      compte de -xy scale.

       v3.2 mai 2014:
            - amélioration du test ps1b (ajoût de règle et réduction
	      du nombre d'indéterminées dans graphes des conflits)

       v3.3 juillet 2014:
            - modification importante du code pour -check ps1
            - modification des graphes linial et linialc
	    - nouvelles options: -check kcolor, -vcolor kcolor, -len, -check kcolorsat

       v3.4 février 2015:
            - documentation et mise au point de l'option -maincc
	    - correction d'un bug lors de la combinaison de "load file" et de "-vcolor pal grad"
	    - correction d'un bug dans la fonction SortInit() qui affectait "-vcolor deg"
	    - correction d'un bug avec l'option -label 1 pour certains graphes (outerplanar ...)
            - création du script dot2gen.awk pour convertir le format dot en format standard
	    - nouvelles options: -fast, -caption
	    - introduction du groupement d'arêtes/arcs i-(j k ...) dans le format standard
	      (il reste un bug si ce format est lu depuis l'entrée standard)

       v3.5 mars 2015:
	    - correction d'un bug pour le dégradé de couleurs avec "-vcolor pal grad"
	    - nouvelles options: -check ncc/connected, -check routing scenario [...] cluster
	    - nouvelle variante de graphe: loadc file
	    - nouveaux graphes: octahedron, turan, hexahedron, tetrahedron, deltohedron,
	      trapezohedron, antiprism, flip, associahedron, shuffle
	    - changement de nom (ajoût du suffixe "hedron") pour: isoca, dodeca, cubocta,
	      rdodeca

       v3.6 juin 2015:
            - nouvelles options: -version, -variant, -check info, -xy mesh,
	      -check routing tzrplg, -check routing dcr
	    - nouveaux graphes: percolation, hgraph, cricket, moth, bull, expander
	    - correction de bug dans "-help ?", dans "-check girth" qui donnait -1 pour
	      un cycle de taille n, dans butterfly (mauvais nombre de paramètres)
	    - vérification que le nombre de sommets n'est pas négatif
	      pour plusieurs graphes dont tree, kout, etc. ce qui
	      pouvait provoquer des erreurs
	    - vérification du format d'entrée pour load/loadc
	    - l'option -check maincc n'affiche plus le graphe initial
	    - description du format des familles de graphes dans l'aide
	    - affiche plus de messages d'erreurs: lecture d'un graphe au mauvais
	      format, mauvaise combinaison de d'options, ...

       v3.7 juillet 2015:
            - nouvelle implémentation des mots de Dyck, qui sont
              maintenant uniformes. En conséquence les graphes
              aléatoires tree, arboricity, rarytree, rbinary, outerplanar,
              kpage ... sont générés de manière uniformes.
	    - nouveaux graphes: treep (et simplification de halin),
	      ygraph, ringarytree (et simplification de arytree,
	      binary, wheel), netgraph, egraph, rgraph (=fish),
	      généralisation de barbell, tadpole (=dragon), pan,
	      banner, paw, theta0 (redéfinition), fan, gem, chess,
	      knight, antelope, camel, giraffe, zebra, utility,
	      zamfirescu, hatzel, web, bdrg, fdrg, regular, cubic, plrg
	    - nouvelles options: -check radius, -check diameter
	    - correction d'un bug si combinaison d'options -check ncc -format list
	      (mauvaise gestion de la variable CHECK)
	    - introduction de caractères UTF8 mathématiques dans l'aide

       v3.8 novembre 2015:
            - nouveaux graphes: ladder, matching, d-octahedron, johnson
	    - correction d'un bug pour le graphe kpage (introduit à v3.7)

       v3.9 février 2016:
	    - renomage de l'option "-xy scale" en "-xy box"
	    - correction d'un bug dans rbinary (mauvais alias/définition)
	    - correction d'un bug (introduit à v3.6) dans la fonction
              max(double,double) au lieu de fmax() qui affectait certains
              graphes et options géométriques (rng, -xy grid, percolation, ...)
            - correction d'un bug concernant rarytree lorsque b>2
	    - correction d'un bug dans les statistiques pour -check routing
	    - généralisation du graphe rarytree (augmentation du degré de la racine)
	    - nouveaux graphes: herschel, goldner-harary, rbinaryz,
	      point, empty, apollonian, dyck, cross, tgraph, bidiakis, gear,
	      centipede, harborth, sunflower, planar, squaregraph, polygon
	    - modification de certains graphes (cage, ring, gabriel,
              nng, rng) afin que le test d'adjacence ne dépende plus de N
              mais de PARAM[0] par exemple
	    - introduction du format %SEED pour l'option -caption

       v4.0 mars 2016:
	    - nouveaux graphes: pat
	    - nouvelles options: -check kindepsat
	    - aide permettant les préfixes comme -check routing cluster
### #*/
