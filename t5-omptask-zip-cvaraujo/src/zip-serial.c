#define _GNU_SOURCE
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

const int max_password_value = 500000;
const char cmd_format[300] = "unzip -P%d -t %s 2>&1";

FILE *popen(const char *command, const char *type);

int main(int argc, char **argv) {
  FILE *fp;
  FILE *input;
  char filename[100];
  char ret[200];
  char cmd[400];
  double t;
  int nt;

  if (argc < 2) {
    fprintf(stderr, "Error: missing path to input file\n");
    return 1;
  }

  if ((input = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Error: could not open file\n");
    return 1;
  }

  // Read inputs
  fscanf(input, "%d", &nt);
  fscanf(input, "%s", filename);

  // Do not touch this line
  omp_set_num_threads(nt);

  t = omp_get_wtime();

  // Parallelize this loop using tasks!
  for (int i = 0; i < max_password_value; i++) {
    sprintf((char *)&cmd, cmd_format, i, filename);

    fp = popen(cmd, "r");
    while (!feof(fp)) {
      fgets((char *)&ret, 200, fp);
      if (strcasestr(ret, "ok") != NULL) {
        printf("Password: %d\n", i);
        i = 500000;
      }
    }
    pclose(fp);
  }

  t = omp_get_wtime() - t;
  fprintf(stderr, "%lf\n", t);
}
