#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include "xmmintrin.h"
#include "pmmintrin.h"
#include "openacc.h"

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


int padfunc(struct dataobj *restrict vp_vec, const int x_M, const int y_M, const int abc_x_l_ltkn, const int abc_x_r_rtkn, const int abc_y_l_ltkn, const int abc_y_r_rtkn, struct profiler * timers, const int x_m, const int y_m, const int nthreads)
{
  float (*restrict vp)[vp_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[vp_vec->size[1]]) vp_vec->data;
  /* Flush denormal numbers to zero in hardware */
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  struct timeval start_section0, end_section0;
  gettimeofday(&start_section0, NULL);
  /* Begin section0 */
  for (int abc_x_l = x_m; abc_x_l <= abc_x_l_ltkn + x_m - 1; abc_x_l += 1)
  {
    {
      #pragma acc kernels
      for (int y = y_m; y <= y_M; y += 1)
      {
        vp[abc_x_l + 6][y + 6] = vp[46][y + 6];
      }
    }
  }
  for (int abc_x_r = -abc_x_r_rtkn + x_M + 1; abc_x_r <= x_M; abc_x_r += 1)
  {
    {
      #pragma acc kernels
      for (int y = y_m; y <= y_M; y += 1)
      {
        vp[abc_x_r + 6][y + 6] = vp[x_M - 34][y + 6];
      }
    }
  }
  {
    #pragma acc kernels
    for (int x = x_m; x <= x_M; x += 1)
    {
      for (int abc_y_l = y_m; abc_y_l <= abc_y_l_ltkn + y_m - 1; abc_y_l += 1)
      {
        vp[x + 6][abc_y_l + 6] = vp[x + 6][46];
      }
      for (int abc_y_r = -abc_y_r_rtkn + y_M + 1; abc_y_r <= y_M; abc_y_r += 1)
      {
        vp[x + 6][abc_y_r + 6] = vp[x + 6][y_M - 34];
      }
    }
  }
  /* End section0 */
  gettimeofday(&end_section0, NULL);
  timers->section0 += (double)(end_section0.tv_sec-start_section0.tv_sec)+(double)(end_section0.tv_usec-start_section0.tv_usec)/1000000;
  return 0;
}
/* Backdoor edit at Tue Mar  3 21:02:50 2020*/ 
/* Backdoor edit at Tue Mar  3 21:04:25 2020*/ 
/* Backdoor edit at Tue Mar  3 21:05:25 2020*/ 
/* Backdoor edit at Tue Mar  3 21:06:28 2020*/ 
/* Backdoor edit at Tue Mar  3 21:12:02 2020*/ 
/* Backdoor edit at Tue Mar  3 21:12:13 2020*/ 
/* Backdoor edit at Tue Mar  3 21:16:22 2020*/ 
/* Backdoor edit at Tue Mar  3 21:16:37 2020*/ 
/* Backdoor edit at Tue Mar  3 21:17:49 2020*/ 
/* Backdoor edit at Tue Mar  3 21:23:55 2020*/ 
/* Backdoor edit at Tue Mar  3 21:24:08 2020*/ 
/* Backdoor edit at Tue Mar  3 21:24:43 2020*/ 
/* Backdoor edit at Tue Mar  3 21:27:51 2020*/ 
/* Backdoor edit at Tue Mar  3 21:28:16 2020*/ 
/* Backdoor edit at Tue Mar  3 21:29:19 2020*/ 
/* Backdoor edit at Tue Mar  3 21:32:22 2020*/ 
/* Backdoor edit at Tue Mar  3 21:34:56 2020*/ 
