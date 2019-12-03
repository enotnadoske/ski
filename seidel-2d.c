/* Include benchmark-specific header. */
#define MINI_DATASET
#include "seidel-2d.h"
double bench_t_start, bench_t_end;

static
double rtclock()
{
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, NULL);
    if (stat != 0)
      printf ("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start()
{
  bench_t_start = rtclock ();
}

void bench_timer_stop()
{
  bench_t_end = rtclock ();
}

void bench_timer_print()
{
  printf ("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start);
}

static
void init_array (int n,
   float A[ n][n])
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = ((float) i*(j+2) + 2) / n;
}

static
void print_array(int n,
   float A[ n][n])
{
  int i, j;

  fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
  fprintf(stderr, "begin dump: %s", "A");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if ((i * n + j) % 20 == 0) fprintf(stderr, "\n");
      fprintf(stderr, "%0.2f ", A[i][j]);
    }
  fprintf(stderr, "\nend   dump: %s\n", "A");
  fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

static
void kernel_seidel_2d(int tsteps,
        int n,
        float A[ n][n])
{
  int t, i, j;
  #pragma omp parallel private (t,i,j)
  {
  for (t = 0; t <= tsteps - 1; t++)
    for (i = 1; i<= n - 2; i++)
      for (j = 1; j <= n - 2; j++)
 A[i][j] = (A[i-1][j-1] + A[i-1][j] + A[i-1][j+1]
     + A[i][j-1] + A[i][j] + A[i][j+1]
     + A[i+1][j-1] + A[i+1][j] + A[i+1][j+1])/9.0f;
  }
}

int main(int argc, char** argv)
{
  int n = N;
  int tsteps = TSTEPS;

  float (*A)[n][n]; A = (float(*)[n][n])malloc ((n) * (n) * sizeof(float));;

  init_array (n, *A);

  bench_timer_start();;

  kernel_seidel_2d (tsteps, n, *A);

  bench_timer_stop();;
  bench_timer_print();;

  if (argc > 42 && ! strcmp(argv[0], "")) print_array(n, *A);

  free((void*)A);;

  return 0;
}
