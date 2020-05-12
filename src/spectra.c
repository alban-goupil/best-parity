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
typedef unsigned int gf_elt;


/* ** Constantes */

/* Quelques constantes sur le corps fini. */
size_t gf_size;                 // Taille du corps
size_t gf_size_minus_1;         // Taille du corps
size_t gf_degree;               // Degré: 2^deg = size
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
  gf_degree = 0;
  for (unsigned p = P; p > 0; p >>= 1) gf_degree++;
  gf_size = 1 << --gf_degree;
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


/* * Mots de code et parités */

/* Initialise un itérateur dans x sur les mots de codes de
   la parité h. */
static inline void cw_begin(size_t n, gf_elt *x) {
  memset(x, 0, n * sizeof *x);
}


/* Passe au mot de code suivant x de la parité h et retourne
   false si plus aucun mot de code n'est obtenable. On
   utilise le fait ici que le dernier coefficient de h est
   1. */
static inline bool cw_next(size_t n, gf_elt *x) {
  size_t i = 0;
  for (i = 0; i < n; ++i)
    if (x[i] == gf_size_minus_1)
      x[i] = gf_zero;
    else
      break;
  if (i == n) return false;
  x[i]++;
  return true;
}


/* Initialise un itérateur dans delta sur les incréments
   delta pour obtenir les voisins à partir de x. */
static inline void delta_begin(size_t n, size_t *d) {
  memset(d, 0, n * sizeof *d);
}


/* Passe à l'incrément suivant. */
static inline bool delta_next(size_t n, size_t *dmax, size_t *d) {
  /* passage au suivant sans coupure */
  size_t i = 0;
  for (i = 0; i < n; ++i)
    if (d[i] >= dmax[i]-1)
      d[i] = 0;
    else
      break;
  if (i == n) return false;
  d[i]++;
  return true;
}



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



/* * Programme principal
 *
 * Il s'agit de trouver la meilleure parité en terme de
 * spectre de distance pour une constellation et un mapping
 * donnés. L'algorithme général est le suivant
 *
 * Pour chaque couple de vecteur (x, y)
 *   Si quadrance (x, y) > qmax
 *     Passer au couple x,y suivant
 *   Pour chaque parité h telle que hx = hy = 0
 *     Incrémenter le spectre de h pour quadrance(x, y)
 * Afficher les parités et leur spectre
 *
 */

int main(int argc, char **argv) {
  /* Paramètres */
  size_t m;            /* GF(2^m) */
  size_t q;            /* Taille constellation == gf_size */
  size_t n = 3;        /* Longueur du code */
  unsigned qmax = UINT_MAX;   /* Intervalle des spectres */ 

  char constfile[81] = ""; /* Nom du fichier de la constellation */
  char mapsfile[81] = "";  /* Nom du ficher du mappings */

  FILE *f = NULL;

  
  /* ** Lecture des arguments */
  if (argc == 5) {       /* Mode ligne de commande */
    if (1 != sscanf(argv[1], "%zu", &n)) error("L'option codelength doit être entier: '%s'", argv[1]);
    if (1 != sscanf(argv[2], "%u", &qmax)) error("L'option qmax doit être entier: '%s'", argv[3]);
    strncpy(constfile, argv[3], 80);
    strncpy(mapsfile, argv[4], 80);
  } else {                      /* Mode aide */
    printf("Usage: %s codelength qmax constellation mappings\n", argv[0]);
    return -1;
  }

  
  /* ** Lecture de la constellation sur constfile et
     ** récupération de sa taille. */
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

  
  /* ** Construction du corps fini */
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

  unsigned qmin = Q[1];
  for (size_t i = 1; i < q; ++i)
    if (qmin > Q[i * q + V[i * q + 1]])
      qmin = Q[i * q + V[i * q + 1]];

#ifdef DEBUG
  printf("Voisinage, qmin = %u\n", qmin);
  for (gf_elt i = 0; i < q; ++i) {
    printf("  %d:", i);
    for (gf_elt j = 0; j < q; ++j)
      printf("\t%d", V[i * q + j]);
    printf("\n");
  }
#endif
  
  /* Allocation de la mémoire */
  gf_elt *h = malloc(n * sizeof *h); // Parity
  gf_elt *x = malloc(n * sizeof *x); // Mot de code x
  gf_elt *y = malloc(n * sizeof *y); // Voisin y
  size_t *d = malloc(n * sizeof *d); // Incrément d pour le voisinage
  size_t *dmax = malloc(n * sizeof *dmax); // Incrément maximal dmax par position

  /* ** Mise en place des spectres */
  /* S[i * qmax + j] retourne le spectre de la parité numéro
     i pour la quadrance j. */
  size_t nspectra = 0;
  chk_begin(n, h);
  do { nspectra++; } while (chk_next(n, h));
  unsigned long *S;
  if (NULL == (S = calloc(nspectra * qmax, sizeof *S)))
    error("Mémoire insuffisance pour les spectres.");

  
  /* Pour chaque couple (x, y) en passage par l'incrément d. */
  cw_begin(n, x);
  do {
    /* Calcul de l'incrément maximal */
    for (size_t i = 0; i < n; ++i) {
      dmax[i] = q;
      unsigned *Qxi = Q + x[i] * q;
      unsigned *Vxi = V + x[i] * q;
      for (size_t j = 0; j < q; ++j)
        if (Qxi[Vxi[j]] >= qmax - qmin) {
          dmax[i] = j;
          break;
        }
    }
    
    delta_begin(n, d);
    do {
      /* Calcul de y en fonction de x et de d et la quandrance. */
      unsigned quad = 0;
      for (size_t i = 0; i < n; ++i) {
        y[i] = V[x[i] * q + d[i]];
        quad += Q[x[i] * q + y[i]];
      }

#ifdef DEBUG
      printf("  x:");
      for (size_t i = 0; i < n; ++i) printf(" %2u", x[i]);
      printf("\ty:");
      for (size_t i = 0; i < n; ++i) printf(" %2u", y[i]);
      printf("\t%u\n", quad);
#endif

      if (quad >= qmax)
        continue;
      
      /* Passe en revue les parités */
      chk_begin(n, h);
      size_t hid = 0;
      do {
        if (chk_valid(n, h, x) && chk_valid(n, h, y)) {
          S[hid * qmax + quad]++;
#ifdef DEBUG
          printf("    h: ");
          for (size_t i = 0; i < n; ++i)
            printf(" %u", h[i]);
          printf("\n");
#endif
        }
        hid++;
      } while (chk_next(n, h));
    } while (delta_next(n, dmax, d));
  } while (cw_next(n, x));


  /* Affichage des résultats */
  printf("Sectra\n");
  chk_begin(n, h);
  size_t hid = 0;
  do {
    printf("%4zu:", hid);
    for (size_t i = 0; i < n; ++i)
      printf(" %2u", h[i]);
    printf(":\t");

    unsigned long sum = 0;
    for (size_t i = 0; i < qmax; ++i) {
      sum += S[hid * qmax + i];
      printf("%lu\t", S[hid * qmax + i]);
    }
    printf("\t(%lu)\n", sum);
    hid++;
  } while (chk_next(n, h));

  
  /* Libération des ressources */
  if (f != stdin) fclose (stdin);
  free(S);
  free(Q);
  free(C);
  free(pi);
  free(y);
  free(x);
  free(h);
  gf_free();
}
