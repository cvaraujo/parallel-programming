#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define SEED 123

void free_matrix(int **m, int size) {
  for (int i = 0; i < size; i++)
    free(m[i]);
  free(m);
}

int **mul(int **a, int **b, int size) {
  int **ret = malloc(size * sizeof(int *));
  int i, j, k, sum;

  for (i = 0; i < size; i++) {
    ret[i] = calloc(size, sizeof(int));
    for (j = 0; j < size; j++){
      sum = 0;
      for (k = 0; k < size; k++)
        sum += a[i][k] * b[k][j];
      ret[i][j] = sum;
    }
  }

  free_matrix(a, size);
  free_matrix(b, size);

  return ret;
}

// Parallelise this function:
int **array_mul(int ***data, int n, int size) {
  long long int i;
  long long int k;
  long long int iter = (n + 1)/2;
  //printf("N = %d\n", n);

  if (n <= 2) {return mul(data[0], data[1], size);}

  //printf("K = %lld, iter = %lld, Maxi = 0\n", k, iter*iter);

#pragma omp parallel shared(n)
#pragma omp single
  {
    for (k = 2; k < (iter * iter); k*=2) {
      for (i = 0; i < n; i+=k) {
        if (i + (k/2) < n){
#pragma omp task firstprivate(i, k) untied depend(inout: data[i]) depend(in: data[i+(k/2)])
          data[i] = mul(data[i], data[i+(k/2)], size);
        }
      }
    }
  }
  return data[0];
}

int **rnd_matrix(int size) {
  int **ret = malloc(size * sizeof(int *));
  for (int i = 0; i < size; i++) {
    ret[i] = malloc(size * sizeof(int));
    for (int j = 0; j < size; j++)
      ret[i][j] = 2 * (rand() % 2) - 1; // Generates -1 or 1
  }

  return ret;
}

void print_matrix(int **m, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++)
      printf("%d ", m[i][j]);
    printf("\n");
  }
}

int main(int argc, char **argv) {
  int n, size;
  double t;
  FILE *input;

  if (argc < 2) {
    fprintf(stderr, "Error: missing path to input file!\n");
    return EXIT_FAILURE;
  }

  if ((input = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Error: could not open input file!\n");
    return EXIT_FAILURE;
  }

  fscanf(input, "%d %d", &n, &size);
  srand(SEED);

  // Do not change this line
  omp_set_num_threads(4);

  int ***data = malloc(n * sizeof(int **));
  for (int i = 0; i < n; i++)
    data[i] = rnd_matrix(size);

  t = omp_get_wtime();
  int **ret = array_mul(data, n, size);
  t = omp_get_wtime() - t;

  print_matrix(ret, size);
  fprintf(stderr, "%lf\n", t);

  free_matrix(ret, size);
  free(data);
  return 0;
}
