#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#define error(format, ...)                                                 \
  do {                                                                     \
    fflush (stdout);                                                       \
    fprintf (stderr, "%s:%d: " format, __FILE__, __LINE__, ##__VA_ARGS__); \
    fprintf (stderr, "\n");                                                \
    exit (EXIT_FAILURE);                                                   \
  } while (0)


#ifdef DEBUG
#define print_array(prompt, fmt, n, arr)                \
  do {                                                  \
    printf(prompt);                                     \
    for (size_t i = 0; i < n; ++i) printf(fmt, arr[i]); \
    printf("\n");                                       \
  } while(0)

#endif


/* * Bibliothèque GF(2^m)
 *
 * Bibliothèque pour les corps finis GF(q) avec q=2^m. Les
 * éléments du corps sont représentés par des entiers de
 * 0 à q-1.
 *
 * Comme la caractéristique est 2, l'addition est le ou
 * exclusif.
 *
 * Pour la multiplication, on utilise une table de
 * logarithme sur la base d'un élément primitif. Seule cette
 * table est nécessaire pour construire le corps complet.
 *
 * Pour éviter de se trimbaler sans cesse la représentation
 * du corps, on fait appel à des variables globales. Il
 * n'est donc possible de travailler qu'avec un corps à la
 * fois.
 *
 * La construction d'un corps fini repose sur un polynôme
 * primitif. Voici une liste pour GF(2^m) pour les premières
 * valeurs de m.
 */


/* Liste de polynômes primitifs pour construire GF(2^m) pour
   m allant de 1 à 10.

   |    q | P                         | P binaire     | P hexdécimal |
   |------+---------------------------+---------------+--------------|
   |    2 | X + 1                     | 11            |          0x3 |
   |    4 | X^2 + X + 1               | 111           |          0x7 |
   |    8 | X^3 + X + 1               | 1011          |          0xb |
   |   16 | X^4 + X + 1               | 1 0011        |         0x13 |
   |   32 | X^5 + X^2 + 1             | 10 0101       |         0x25 |
   |   64 | X^6 + X + 1               | 100 0011      |         0x43 |
   |  128 | X^7 + X^3 + 1             | 1000 1001     |         0x89 |
   |  256 | X^8 + X^4 + X^3 + X^2 + 1 | 1 0001 1101   |        0x11d |
   |  512 | X^9 + X^5 + 1             | 10 0010 0001  |        0x221 |
   | 1024 | X^10 + X^3 + 1            | 100 0000 1001 |        0x409 |
 */
const size_t n_primitives = 10;
const unsigned primitives[] = {0x1, 0x3, 0x7, 0xb, 0x13, 0x25,
                               0x43, 0x89, 0x11d, 0x221, 0x409};


/* ** Types */

/* Un élément de GF(q) est représenté par un entier non
   signé entre 0 et q-1.*/
typedef unsigned gf_elt;


/* ** Constantes */

/* Quelques constantes sur le corps fini. */
size_t gf_size;                 // Taille du corps
size_t gf_size_minus_1;         // Taille du corps
const gf_elt gf_zero = 0;       // Le zéro
const gf_elt gf_one = 1;        // Le un
gf_elt gf_alpha;                // Un élément primitif


/* La table des logarithmes: un élément non nul $x$ de GF(q)
   s'écrit comme une puissance $\alpha^k$ d'un élément
   primitif $\alpha$. Cette table indicée par $x$ retourne
   la puissance $k$ correspondante. Pour 0, elle
   retourne -1. La table gf_exp fait l'inverse, elle
   retourne $\alpha^k$ à partir de l'indice $k$. */
int *gf_log;
gf_elt *gf_exp;


/* ** Opérations algébriques */

/* Calcul acc <- acc + \alpha^hlog * x. hlog doit être entre
   0 et gf_size_minus_1 exclu. */
static inline void gf_accmul(gf_elt *acc, int hlog, gf_elt x) {
  if (x == gf_zero) return;
  *acc ^= gf_exp[(hlog + gf_log[x]) % gf_size_minus_1]; 
}


/* ** Initialisation et libération */ 

/* Initialisation du corps fini GF(q) à partir d'un polynôme
   primitif P qui est donné par sa représentation binaire.
   Ainsi, pour GF(8), le polynôme P = X^3 + X + 1 se
   représente par l'entier 11 en décimal c'est-à-dire 1011
   en binaire.
*/
static int gf_init(unsigned P) {
  // Lecture de la taille en lisant le degré de P
  gf_size = 1;
  for (unsigned p = P; p > 0; p >>= 1)
    gf_size <<= 1;
  gf_size >>= 1;
  gf_size_minus_1 = gf_size - 1;
  
  // Construction des tables de log et de puissances
  gf_log = malloc(gf_size * sizeof *gf_log);
  gf_exp = malloc(gf_size * sizeof *gf_exp); // Normalement gf_size_minus_1....
  gf_log[gf_zero] = -1; // 0 n'est pas une puissance de alpha !

  gf_elt x = gf_one;
  for (size_t k = 0; k < gf_size_minus_1; k++) {
    gf_log[x] = k;              // remplissage des tables
    gf_exp[k] = x;
    x <<= 1;                    // puissance suivante
    if (x >= gf_size)
      x ^= P;
  }
  gf_alpha = gf_exp[1];         // élément primitif

  if (x != gf_one)             // Soucis, alpha pas primitif
    error("gf_init: Corps fini incorrect: polynôme non primitif\n");

  return 0;
}


/* Libération des ressources utilisées par le corps */
static int gf_free() {
  free(gf_log);
  free(gf_exp);
  gf_log = NULL;
  gf_exp = NULL;
  return 0;
}


/* * Parity check
 *
 * Il faudra itérer sur les différentes parités de taille
 * n donnée. Pour des raisons évidentes de symétrie, il
 * n'est pas nécessaire de faire le tour des (q-1)^n parités
 * dont aucun coefficient n'est nul, mais il suffit de
 * regarder les multi-ensembles ce qui donne ((q-1
 * multichoose n)) = (q-1+n-1 choose n) possibilités. À cela
 * s'ajoute aussi qu'un coefficient peut toujours être mis
 * à 1. Le dernier jouera ce rôle ici.
 */

/* Indique si une parité h de longueur n est vérifiée par le
   vecteur x. */
static inline bool chk_valid(size_t n, gf_elt *h, gf_elt *x) {
  gf_elt acc = gf_zero;
  for (size_t i = 0; i < n-1; ++i)
    gf_accmul(&acc, *h++, *x++);
  return acc == *x;
}


/* Initialise un itérateur sur les parités remplissant h. */
static inline void chk_begin(size_t n, gf_elt *h) {
  for (size_t i = 0; i < n-1; ++i)
    h[i] = n-i-1;
  h[n-1] = 0;
}


/* Passe à la parité suivante */
static inline bool chk_next(size_t n, gf_elt *h) {
  size_t i;
  for (i = 0; h[i] >= gf_size_minus_1 - i - 1; ++i)
    if (i + 3 > n)
      return false;
  for (h[i]++; i; i--)
    h[i-1] = h[i] + 1;
  return true;
}


/* Initialise un itérateur dans x sur les mots de codes de
   la parité h. */
static inline void cw_begin(size_t n, gf_elt *h, gf_elt *x) {
  for (size_t i = 0; i < n; ++i)
    x[i] = gf_zero;
}


/* Passe au mot de code suivant x de la parité h et retourne
   false si plus aucun mot de code n'est obtenable. On
   utilise le fait ici que le dernier coefficient de h est
   1. */
static inline bool cw_next(size_t n, gf_elt *h, gf_elt *x) {
  size_t i = 0;
  for (i = 0; i < n-1; ++i)
    if (x[i] == gf_size_minus_1)
      x[i] = gf_zero;
    else
      break;
  if (i == n-1) return false;
  x[i]++;
  
  gf_elt acc = gf_zero;
  for (i = 0; i < n-1; ++i)
    gf_accmul(&acc, h[i], x[i]);
  x[n-1] = acc;
  return true;
}



/* * Spectre
 *
 *  Le spectre est un tableau S tel que S[q] est la
 *  multiplicité de la quadrance q. q est limité par qmax.
 */


/* Compare deux spectres. Le plus grand le meilleure. */
static inline int sp_cmp(unsigned long *a, unsigned long *b, unsigned qmax) {
  for (size_t i = 2; i < qmax; ++i)
    if (a[i] == b[i])
      continue;
    else  
      return (b[i] > a[i]) ? 1 : -1;
  return 0;
}



/* * Programme principal
 *
 * Il s'agit de trouver la meilleure parité en terme de
 * spectre de distance pour une constellation et un mapping
 * donnés. L'algorithme général est le suivant
 *
 * Initialiser le spectre S_best
 * Pour chaque parité h
 *   Initialiser un spectre de distance S
 *   Pour chaque mot de code x
 *     Pour chaque vecteur y voisin de x
 *       Si y est un mot de code alors
 *         Ajouter distance(x, y) dans spectre S
 *   Si le spectre S est meilleur que S_best alors
 *     S_best <- S
 * Afficher la parité h et son meilleur spectre S_best
 *
 * Pour rendre les choses plus rapide, l'élagage de la
 * recherche se fait avec la notion de y voisin de x, cf.
 * [ABCND1x]. Il suffit de se limiter à une distance seuil
 * au dessus de laquelle il n'est pas nécessaire de
 * connaître le spectre.
 *
 * Enfin le travail précédent est repété pour plusieurs
 * mappings, pour repérer les meilleures paires
 * (parité/mapping).
 */

int main(int argc, char **argv) {
  size_t m;            /* GF(2^m) */
  size_t q;            /* Taille constellation == gf_size */
  size_t n = 3;        /* Longueur du code */
  unsigned qmax = UINT_MAX;   /* Intervalle des spectres */ 

  char constfile[81] = ""; /* Nom du fichier de constellation */
  char mapsfile[81] = "";  /* Nom du ficher de mappings */
  FILE *f = NULL;

  /* Lecture des arguments ou de l'entrée standard */
  if (argc == 1) {              /* Mode interactif */
    printf("codelength: "); fflush(stdout);
    if (1 != scanf("%zu", &n)) error("Codelength doit être entier");
    printf("max quadrance: "); fflush(stdout);
    if (1 != scanf("%u", &qmax)) error("Qmax doit être entier");
    printf("constellation file: "); fflush(stdout);
    if (1 != scanf("%80s", constfile)) error("Nom de fichier imcompréhensible"); 
    printf("mappings file: "); fflush(stdout);
    if (1 != scanf("%80s", mapsfile)) error("Nom de fichier imcompréhensible");
  } else if (argc == 5) {       /* Mode ligne de commande */
    if (1 != sscanf(argv[1], "%zu", &n)) error("L'option codelength doit être entier: '%s'", argv[1]);
    if (1 != sscanf(argv[2], "%u", &qmax)) error("L'option qmax doit être entier: '%s'", argv[3]);
    strncpy(constfile, argv[3], 80);
    strncpy(mapsfile, argv[4], 80);
  } else {                      /* Mode aide */
    printf("Usage: %s codelength qmin qmax constellation mappings\n", argv[0]);
    return -1;
  }
    
  /* Lecture de la constellation sur constfile */
  if (NULL == (f = fopen(constfile, "r")))
    error("Impossible d'ouvrir le fichier constellation '%s'.", constfile);
  q = 0;
  m = 0;
  struct { int x, y; } *C = malloc((1 << m) * sizeof *C);
  while (2 == fscanf(f, "%u%u", &C[q].x, &C[q].y))
    if (++q == (1 << m))
      if (NULL == (C = realloc(C, (1 << ++m) * sizeof *C)))
        error("Allocation mémoire impossible.");
  m--;
  fclose(f);    

  printf("Constellation: %s\n", constfile);
#ifdef DEBUG
  for (size_t i = 0; i < q; ++i)
    printf("  %zu:\t% d\t% d\n", i, C[i].x, C[i].y);
#endif

  /* Construire le corps en premier */
  if (1 << m != q) error("Corps de caractéristique 2 uniquement.");
  if (n >= q) error("codelength doit être inférieur à l'ordre du corps."); 

  gf_init(primitives[m]);
  printf("GF(%zu = 2^%zu)\n", q, m);  
  printf("codelength: %zu\n", n);
  printf("qmax: %u\n", qmax); 
  

  /* ** Lecture du mapping */
  size_t *pi = malloc(q * sizeof *pi); // Le mapping
  
  if (strcmp (mapsfile, "-") == 0)
    f = stdin;
  else if (NULL == (f = fopen(mapsfile, "r")))
    error("Impossible d'ouvrir le fichier mapping '%s'.", mapsfile); 

  for (size_t i = 0; i < q; ++i)
    if (1 != fscanf(f, "%zu", pi + i))
      error("Fichier mapping incomplet\n");
    else if (pi[i] >= q)
      error("Valeur de mapping hors bornes\n");

  printf("Mapping:");
  for (size_t i = 0; i < q; ++i)
    printf(" %zu", pi[i]);
  printf("\n");


  /* ** Calcul des quadrances et voisinages */
  
  /* Q[i*q + j] retourne la quadrance (distance au carré)
     entre les points de la constellation donnés par pi[i]
     et pi[j] pour i et j éléments de GF(q). */
  unsigned *Q = malloc(q * q * sizeof *Q);
  for (gf_elt i = 0; i < q; ++i)
    for (gf_elt j = 0; j < q; ++j) {
      unsigned dx = C[pi[i]].x - C[pi[j]].x;
      unsigned dy = C[pi[i]].y - C[pi[j]].y;
      Q[i * q + j] = dx * dx + dy * dy;
    }

#ifdef DEBUG
  printf("Quadrances\n");
  for (gf_elt i = 0; i < q; ++i) {
    printf("  %2u:", i);
    for (gf_elt j = 0; j < q; ++j) {
      printf("\t%u", Q[i * q + j]);
    }
    printf("\n");
  }
#endif


  /* V[i * q + j] retourne le j-ème élément du corps fini
     qui est le plus proche de l'élément i selon la
     quadrance Q. */ 
  unsigned *V = malloc(q * q * sizeof *V);
  for (gf_elt i = 0; i < q; ++i) {
    /* Alias sur les lignes V[i*q + ...]  et Q[i*q + ...] */
    unsigned *Vi = V + i * q;
    unsigned *Qi = Q + i * q;
    
    /* Initialisation sur la permutation identité */
    for (gf_elt j = 0; j < q; ++j) Vi[j] = j;
    
    /* On trie Vi sur la quadrance croissante */
    for (gf_elt j = 0; j < q; ++j)
      for (gf_elt k = j+1; k < q; ++k)
        if (Qi[Vi[k]] < Qi[Vi[j]]) {
          size_t itmp = Vi[j];
          Vi[j] = Vi[k];
          Vi[k] = itmp;
        }
  }

  /* Recherche du cappage sur delta */
  unsigned cqmin = Q[1];
  for (size_t i = 1; i < q; ++i)
    if (cqmin > Q[i * q + V[i * q + 1]])
      cqmin = Q[i * q + V[i * q + 1]];

  size_t *dcap = malloc(n * q * sizeof *dcap);
  for (size_t w = 0; w < n; ++w) {
    for (gf_elt j = 0; j < q; ++j) {
      size_t d = 0;
      for (d = 0; d < q; ++d)
        if (Q[j * q + V[j * q + d]] > qmax - (w > 2 ? w-1 : 1) * cqmin)
          break;
      dcap[w * q + j] = d;
    }
  }

#ifdef DEBUG
  printf("Voisinage, cqmin = %u\n", cqmin);
  for (gf_elt i = 0; i < q; ++i) {
    printf("  %d:", i);
    for (gf_elt j = 0; j < q; ++j)
      printf("\t%d", V[i * q + j]);
    printf("\n");
  }

  printf("dcap\n");
  for (size_t w = 0; w < n; ++w) {
    printf("  %zu:\t", w);
    for (gf_elt i = 0; i < q; ++i) printf(" %zu", dcap[w * q + i]);
    printf("\n");
  }
#endif

  
  /* Allocation de la mémoire */
  gf_elt *h = malloc(n * sizeof *h); // Parity
  gf_elt *x = malloc(n * sizeof *x); // Mot de code x
  gf_elt *y = malloc(n * sizeof *y); // Voisin y

  /* Pour les incréments on utilise l'algo H de TACOP */
  size_t *da = malloc(n * sizeof *da); // Incrément d
  size_t *df = malloc(n * sizeof *df);
  int *dd = malloc(n * sizeof *dd); // Direction
  size_t dw = 0;
  
  /* Pour les spectres */
  unsigned long *Sbest = calloc(qmax,  sizeof *Sbest); // Meilleur spectre
  unsigned long *Scur = calloc(qmax, sizeof *Scur); // Spectre courant
  for (size_t i=0; i < qmax; ++i) Sbest[i] = UINT_MAX; // Pour un grand minimum initial

      
  /* Pour chaque parité h de longueur n */
  chk_begin(n, h);
  do {
#ifdef DEBUG
    print_array("h: ", "%u ", n, h);
#endif
    memset(Scur, 0, qmax * sizeof *Scur); /* RAZ du Spectre pour cette parité */
      
    /* Pour chaque mot de code x */
    cw_begin(n, h, x);
    do {
      /* Initialisation des deltas */
      for (size_t j = 0; j < n-1; ++j) {
        da[j] = 0;
        df[j] = j;
        dd[j] = 1;
      }
      df[n-1] = n-1;
      dw = 0;
      
      for (;;) {
        /* Récupère un mot de code à partir des incréments
           et calcul de la quadrance. */
        unsigned quad = 0;
        y[n-1] = gf_zero;
        for (size_t i = 0; i < n-1; ++i) {
          y[i] = V[x[i] * q + da[i]];
          quad += Q[x[i] * q + y[i]];
          gf_accmul(y+n-1, h[i], y[i]);
        }
        quad += Q[x[n-1] * q + y[n-1]];

#ifdef DEBUG
        printf("x:"); for (size_t i=0; i < n; ++i) printf(" %2u", x[i]);
        printf("\ty:"); for (size_t i=0; i < n; ++i) printf(" %2u", y[i]);
        printf("\td:"); for (size_t i=0; i < n-1; ++i) printf(" %2zu", da[i]);
        printf("\tquad: %u, wt: %zu\n", quad, dw);
#endif

        /* Ici nous avons un voisin mot de code à bonne distance ! */
        if (quad < qmax)
          Scur[quad]++;

        /* Voisin suivant */
        size_t j = df[0];
        df[0] = 0;

        if (j == n-1) break; // Terminate
        if (da[j] == 0) dw++;
        da[j] += dd[j];
        if (da[j] == 0) dw--;

        if (da[j] == 0 || da[j] == q-1 || da[j] + dd[j] >= dcap[dw * q + x[j]]) {
        /* if (da[j] == 0 || da[j] == q-1) { */
          dd[j] = -dd[j];
          df[j] = df[j+1];
          df[j+1] = j+1;
        }
      }

      /* On vérifie aussi que le spectre courant reste meilleur que le meilleur jusqu'ici */
    } while(cw_next(n, h, x) && sp_cmp(Sbest, Scur, qmax) <= 0); /* Mot de code suivant */
#ifdef DEBUG
    /* H en premier */
    printf("  h:");
    for (size_t i = 0; i < n; ++i)
      printf("%2u ", h[i]);
    printf("\t");
    /* S ensuite */
    unsigned long sum = 0;
    for (size_t i = 0; i < qmax; ++i) {
      sum += Scur[i];
      printf("%lu\t", Scur[i]);
    }
    printf("\t(%lu)\n", sum);
    fflush(stdout);
#endif

    /* On garde ce spectre si c'est le meilleur jusrq'ici. */
    if (sp_cmp(Sbest, Scur, qmax) <= 0) {
      memcpy(Sbest, Scur, qmax * sizeof *Sbest);

      /* H en premier */
      for (size_t i = 0; i < n; ++i)
        printf("%2u ", h[i]);
      printf("\t");
      /* S ensuite */
      unsigned long sum = 0;
      for (size_t i = 0; i < qmax; ++i) {
        sum += Sbest[i];
        printf("%lu\t", Sbest[i]);
      }
      printf("\t(%lu)\n", sum);
      fflush(stdout);
    }
  } while (chk_next(n, h));    /* Parité suivante */


  /* Libération des ressources */
  if (f != stdin) fclose (stdin);
  free(Scur);
  free(Sbest);
  free(V);
  free(Q);
  free(C);
  free(pi);
  free(da);
  free(df);
  free(dd);
  free(y);
  free(x);
  free(h);
  gf_free();
}
