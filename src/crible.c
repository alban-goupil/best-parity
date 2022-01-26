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
const unsigned n_primitives = 10;
const unsigned primitives[] = {0x1, 0x3, 0x7, 0xb, 0x13, 0x25,
                               0x43, 0x89, 0x11d, 0x221, 0x409};


/* ** Types */

/* Un élément de GF(q) est représenté par un entier non
   signé entre 0 et q-1.*/
typedef unsigned gf_elt;


/* ** Constantes */

/* Quelques constantes sur le corps fini. */
unsigned gf_size;                 // Taille du corps
unsigned gf_size_minus_1;         // Taille du corps
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
  for (unsigned k = 0; k < gf_size_minus_1; k++) {
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


/* * Bibliothèque de combinatoire */


/* Retourne la partition suivante de l'entier n. p est un
   tableau d'entier par ordre décroissant. La fonction
   retourne le nombre d'entier de la décomposition ou 0 si
   c'est la dernière. */
static int partnext(int n, int *p) {
  /* Réécriture de <www.geeksforgeeks.org/generate-unique-partitions-of-an-integer> */
  int k;
  int rem = 0;
  
  for (k = n-1; k >=0 && p[k] <= 1; --k)
    rem += p[k];

  if (k < 0) return 0;
  p[k]--;
  rem++;

  for (; rem > p[k]; k++) {
    p[k+1] = p[k];
    rem -= p[k];
  }

  p[++k] = rem;
  return k+1 < n;
} 


/* Retourne la prochaine permutation de p selon l'algorithme
   3.12.4 de [P. Cameron, "Combinatorics: Topics,
   Techniques, Algorithms," Cambridge University Press,
   1994]. */
static int permnext(int n, int *p) {
  int j, k;

  for (j = n-2; j >= 0; --j)
    if (p[j] < p[j+1])
      break;

  if (j < 0) return 0;

  for (k = n-1; k > j; --k)
    if (p[k] > p[j])
      break;

  int tmp = p[j];
  p[j] = p[k];
  p[k] = tmp;

  for (k = n-1, j++; j < k; j++, k--) {
    int tmp = p[j];
    p[j] = p[k];
    p[k] = tmp;
  }
  
  return 1;
}

/* Passe au mot de code suivant x de la parité h et retourne
   false si plus aucun mot de code n'est obtenable. On
   utilise le fait ici que le dernier coefficient de h est
   1. */
static inline bool cwnext(unsigned n, gf_elt *h, gf_elt *x) {
  unsigned i = 0;
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



/* * Programme principal
 *
 * Il s'agit de trouver la meilleure parité en terme de
 * spectre de distance pour une constellation et un mapping
 * donnés. L'algorithme général est le suivant
 */

static int readparity(FILE *f, int n, gf_elt *h) {
  for (unsigned i = 0; i < n; ++i)
    if (1 != fscanf(f, "%u", h + i))
      return 0;

  /* Remet la parité h sans perte de généralité sous la
     forme h0 h1 ... 0 avec h0 >= h1 >=... >= 0 */
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = i+1; j < n; ++j)
      if (h[i] < h[j]) {
        gf_elt tmp = h[i];
        h[i] = h[j];
        h[j] = tmp;
      }

  for (unsigned i = 0; i < n; ++i)
    h[i] -= h[n-1];
  
  return 1;
}


int main(int argc, char **argv) {
  unsigned m;            /* GF(2^m) */
  unsigned q;            /* Taille constellation == gf_size */
  unsigned n;            /* Longueur du code */
  unsigned quad;       /* Quadrance à repérer */ 

  char constfile[81] = ""; /* Nom du fichier de constellation */
  char mapsfile[81] = "";  /* Nom du ficher de mappings */
  char hfile[81] = "";  /* Nom du ficher de parités ou - pour stdin */
  FILE *f = NULL;

  /* Lecture des arguments ou de l'entrée standard */
  if (argc == 6) {       /* Mode ligne de commande */
    if (1 != sscanf(argv[1], "%u", &n)) error("L'option codelength doit être entier: '%s'", argv[1]);
    if (1 != sscanf(argv[2], "%u", &quad)) error("L'option quad doit être entier: '%s'", argv[3]);
    strncpy(constfile, argv[3], 80);
    strncpy(mapsfile, argv[4], 80);
    strncpy(hfile, argv[5], 80);
  } else {                      /* Mode aide */
    printf("Usage: %s codelength quad constellation mappings parities\n", argv[0]);
    return -1;
  }

  
  /* ** Lecture de la constellation sur constfile */
  if (NULL == (f = fopen(constfile, "r")))
    error("Impossible d'ouvrir le fichier constellation '%s'.", constfile);
  q = 0;
  m = 0;
  struct { int I, Q; } *C = malloc((1 << m) * sizeof *C);
  while (2 == fscanf(f, "%u%u", &C[q].I, &C[q].Q))
    if (++q == (1 << m))
      if (NULL == (C = realloc(C, (1 << ++m) * sizeof *C)))
        error("Allocation mémoire impossible.");
  m--;
  fclose(f);    

  printf("Constellation: %s\n", constfile);


  /* ** Lecture du mapping */
  unsigned *pi = malloc(q * sizeof *pi); // Le mapping
  
  if (strcmp (mapsfile, "-") == 0)
    f = stdin;
  else if (NULL == (f = fopen(mapsfile, "r")))
    error("Impossible d'ouvrir le fichier mapping '%s'.", mapsfile); 

  for (unsigned i = 0; i < q; ++i)
    if (1 != fscanf(f, "%u", pi + i))
      error("Fichier mapping incomplet\n");
    else if (pi[i] >= q)
      error("Valeur de mapping hors bornes\n");

  printf("Mapping:");
  for (unsigned i = 0; i < q; ++i)
    printf(" %u", pi[i]);
  printf("\n");
  if (f != stdin)
    fclose(f);

  
  /* ** Construction du corps q=2^m */
  if (1 << m != q) error("Corps de caractéristique 2 uniquement.");
  if (n >= q) error("codelength doit être inférieur à l'ordre du corps."); 

  gf_init(primitives[m]);
  printf("GF(%u = 2^%u)\n", q, m);  
  printf("codelength: %u\n", n);
  printf("quad: %u\n", quad); 


  /* ** Quadrances entre éléments */
  unsigned *quadrances = malloc(q * q * sizeof quadrances);
  for (gf_elt x = 0; x < q; ++x)
    for (gf_elt y = 0; y < q; y++) {
      unsigned dx = C[pi[x]].I - C[pi[y]].I;
      unsigned dy = C[pi[x]].Q - C[pi[y]].Q;
      quadrances[x * q + y] = dx * dx + dy * dy;
    }
  

  /* ** Répartition des couples (x, y) selon leur distance */
  gf_elt **cercles = calloc(q * (1 + quad), sizeof *cercles);
  unsigned *perims = calloc(q * (1 + quad), sizeof *perims);
  for (unsigned i = 0; i < q * (1 + quad); ++i)
    cercles[i] = malloc(8 * sizeof *cercles);
  
  for (gf_elt x = 0; x < q; ++x)
    for (gf_elt y = 0; y < q; y++) {
      unsigned r = quadrances[x * q + y];
      if (r > quad) continue;
      
      unsigned i = r * q + x;
      cercles[i][perims[i]++] = y;
      
      unsigned k = perims[i];
      if (k >= 8 && (k & (k-1)) == 0)
        cercles[i] = realloc(cercles[i], 2 * k * sizeof *cercles);
    }

#if DEBUG
  printf("Cercles");
  for (gf_elt x = 0; x < q; ++x) {
    printf("  x: %u\n", x);
    for (unsigned r = 0; r <= quad; ++r)
      if (perims[x + r * q] > 0) {
        printf("    %u:", r);
        for (unsigned i = 0; i < perims[x + r * q]; ++i)
          printf(" %u", cercles[x + r * q][i]);
        printf("\n");
      }
    }
#endif
  
  
  /* ** Allocation de la mémoire */
  gf_elt *x = malloc(n * sizeof *x); // Mot de code
  gf_elt *h = malloc(n * sizeof *h); // Parité
  int *part = malloc(n * sizeof *part); // Partition entière de quad
  unsigned *idx = malloc(n * sizeof *idx); // Index sur les couples 
  unsigned bestmult = UINT_MAX; // Meilleure multiplicité
  
  
  /* Pour chaque parité h de longueur n */
  if (strcmp (hfile, "-") == 0)
    f = stdin;
  else if (NULL == (f = fopen(hfile, "r")))
    error("Impossible d'ouvrir le fichier de parités '%s'.", hfile); 

  while (readparity(f, n, h)) {
#ifdef DEBUG
    printf("  h:");
    for (unsigned i = 0; i < n; ++i) printf(" %2u", h[i]);
    printf("\n");
#endif

    unsigned mult = 0;          // Multiplicité
    
    /* Initialisation de la partition */
    for (unsigned i = 2; i < n; ++i) part[i] = 0;
    part[0] = quad - 1;
    part[1] = 1;

    do {
      /* Pour chaque permutation de cette partition */
      for (int i = 0, j = n-1; i < j; ++i, --j) {
        int tmp = part[i];
        part[i] = part[j];
        part[j] = tmp;
      }
      
      do {
#ifdef DEBUG
        printf("  p:");
        for (unsigned i = 0; i < n; ++i) printf(" %2u", part[i]);
        printf("\n");
#endif

        /* Pour chaque mot de code */
        for (unsigned i = 0; i < n; ++i) x[i] = 0;
        do {
          for (unsigned i = 0; i < n-1; ++i) idx[i] = 0;
          for (;;) {
            gf_elt y = 0;
            
            for (unsigned i = 0; i < n-1; ++i)
              gf_accmul(&y, h[i], cercles[x[i] + q * part[i]][idx[i]]);

#ifdef DEBUG
            printf("  idx:");
            for (unsigned i = 0; i < n; ++i) printf(" %u", idx[i]);
            printf("\tx:");
            for (unsigned i = 0; i < n; ++i) printf(" %u", x[i]);
            printf("\ty:");
            for (unsigned i = 0; i < n; ++i) printf(" %u", cercles[x[i] + q * part[i]][idx[i]]);
            printf("\t%c\n", y == 0 ? '*' : ' ');
#endif

            if (quadrances[x[n-1] * q + y] == part[n-1])
              if (++mult > bestmult)
                goto nexth;

            unsigned i;
            for (i = 0; i < n-1 && ++idx[i] == perims[x[i] + q * part[i]]; ++i)
                idx[i] = 0;
            if (i == n-1) break;
          } /* Idx */
        } while (cwnext(n, h, x));
      } while (permnext(n, part));
    } while (partnext(n, part));
    
  nexth:
    if (mult <= bestmult) {
      bestmult = mult;
      for (unsigned i = 0; i < n; ++i) printf("%2u ", h[i]);
      printf("\t%d\n", bestmult);
      fflush(stdout);
    }
  }

  /* Libération des ressources */
  free(idx);
  free(part);
  free(h);
  gf_free();
}
