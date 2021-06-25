#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char *argv[]);
int prime_default(int n);

int main(int argc, char *argv[]) {
  int n;
  int n_factor;
  int n_hi;
  int n_lo;
  int primes;
  double t;
  FILE *input;

  if (argc < 2) {
    fprintf(stderr, "Error: missing path to input file\n");
    return 1;
  }

  if ((input = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Error: could not open file\n");
    return 1;
  }

  n_lo = 1;
  n_factor = 2;

  fscanf(input, "%d", &n_hi);
  n_hi = 1 << n_hi;

  printf("                    \n");
  printf("         N     Pi(N)\n");
  printf("\n");

  n = n_lo;

  t = omp_get_wtime();

  while (n <= n_hi) {
    primes = prime_default(n);

    printf("  %8d  %8d\n", n, primes);

    n = n * n_factor;
  }

  t = omp_get_wtime() - t;

  /*
    Terminate.
  */
  fprintf(stderr, "%lf\n", t);

  return 0;
}

/*
  Purpose:
   counts primes.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    10 July 2010
  Author:
    John Burkardt
  Parameters:
    Input, the maximum number to check.
    Output, the number of prime numbers up to N.
*/
int prime_default(int n) {
  int i;
  int j;
  int prime;
  int total = 0;

  for (i = 2; i <= n; i++) {
    prime = 1;

    for (j = 2; j < i; j++) {
      if (i % j == 0) {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }

  return total;
}
