#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>


int main(int argc, char **argv) {
  if (argc != 3) {
    printf("Usage: %s k n\n", argv[0]);
    return EXIT_FAILURE;
  }

  int k = atoi(argv[1]);
  int n = atoi(argv[2]);

  unsigned *a = malloc(k * sizeof *a);
  for (size_t i=0; i < k; ++i) a[i] = i;
  
  for(;;) {
    for (size_t i = 0; i < k; ++i)
      printf("%u%c", a[i], i == k-1 ? '\n' : ' ');

    int i;
    for (i = 0; i < k-1; ++i)
      if (a[i] + 1 != a[i+1])
        break;
    
    if (i == k-1 && a[i] == n-1) break;
    a[i]++;
    for (i--; i>=0; --i)
      a[i] = i;
  }
  
  return EXIT_SUCCESS;
}
