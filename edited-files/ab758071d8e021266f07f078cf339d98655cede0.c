#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include "omp.h"

struct dataobj
{
  void *restrict data;
  int * size;
  int * npsize;
  int * dsize;
  int * hsize;
  int * hofs;
  int * oofs;
} ;

struct profiler
{
  double section0;
} ;


int norm2(struct dataobj *restrict n_vec, struct dataobj *restrict u_vec, const int time_M, const int time_m, struct profiler * timers, const int x_M, const int x_m, const int y_M, const int y_m, const int z_M, const int z_m)
{
  float (*restrict n) __attribute__ ((aligned (64))) = (float (*)) n_vec->data;
  float (*restrict u)[u_vec->size[1]][u_vec->size[2]][u_vec->size[3]] __attribute__ ((aligned (64))) = (float (*)[u_vec->size[1]][u_vec->size[2]][u_vec->size[3]]) u_vec->data;
  #pragma omp target enter data map(to: u[0:u_vec->size[0]][0:u_vec->size[1]][0:u_vec->size[2]][0:u_vec->size[3]])
  float sum = 0.0F;
  struct timeval start_section0, end_section0;
  gettimeofday(&start_section0, NULL);
  /* Begin section0 */
  for (int time = time_m, t0 = (time)%(3); time <= time_M; time += 1, t0 = (time)%(3))
  {
    #pragma omp target teams distribute parallel for collapse(3) reduction(+:sum)
    for (int x = x_m; x <= x_M; x += 1)
    {
      for (int y = y_m; y <= y_M; y += 1)
      {
        for (int z = z_m; z <= z_M; z += 1)
        {
          sum += fabs(pow(u[t0][x + 12][y + 12][z + 12], 2));
        }
      }
    }
  }
  /* End section0 */
  gettimeofday(&end_section0, NULL);
  timers->section0 += (double)(end_section0.tv_sec-start_section0.tv_sec)+(double)(end_section0.tv_usec-start_section0.tv_usec)/1000000;
  n[0] = sum;
  #pragma omp target exit data map(delete: u[0:u_vec->size[0]][0:u_vec->size[1]][0:u_vec->size[2]][0:u_vec->size[3]])
  return 0;
}
