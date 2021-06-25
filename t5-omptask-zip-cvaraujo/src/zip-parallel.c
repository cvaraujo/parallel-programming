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
  FILE *input;
  char filename[100];
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

  int j = 0, i, done = 0;

  FILE *fp;
  char ret[200];
  char cmd[400];
  int nts = 150;
  int ts = max_password_value/nts;
  //int ts = (max_password_value + (nts-1))/nts;
  //printf("%d - %d\n", nts, ts);

  t = omp_get_wtime();
  
#pragma omp parallel num_threads(nt)
#pragma omp single
  {
    for (j = 0; j < nts && !done; j++) {
      for (i = j; i < max_password_value && !done; i+=ts) {
        //        printf("J = %d, I = %d\n", j, i);
        //getchar();
        /* for (i = 0; i < max_password_value && !done; i++) { */
#pragma omp task shared(done) firstprivate(i, fp, ret, cmd) untied if(!done) final(!done)
        {
          sprintf((char *)&cmd, cmd_format, i, filename);
          fp = popen(cmd, "r");
          int k = 1;
        
          while (!feof(fp)) {
            if (k > 2) break;
          
            fgets((char *)&ret, 200, fp);
            if (strcasestr(ret, "ok") != NULL) {
              printf("Password: %d\n", i);
              i = max_password_value;
              done = 1;
              break;
            } else k++;
          }
          pclose(fp);
        }
      }
    }
  }
  t = omp_get_wtime() - t;
  fprintf(stderr, "%lf\n", t);
}
