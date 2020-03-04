#define _POSIX_C_SOURCE 200809L
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
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
  double section1;
  double section2;
} ;


int Forward(struct dataobj *restrict damp_vec, const float dt, const float o_x, const float o_y, struct dataobj *restrict rec_vec, struct dataobj *restrict rec_coords_vec, struct dataobj *restrict src_vec, struct dataobj *restrict src_coords_vec, struct dataobj *restrict u_vec, struct dataobj *restrict vp_vec, const int x_M, const int x_m, const int y_M, const int y_m, const int p_rec_M, const int p_rec_m, const int p_src_M, const int p_src_m, const int time_M, const int time_m, struct profiler * timers)
{
  float (*restrict damp)[damp_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[damp_vec->size[1]]) damp_vec->data;
  float (*restrict rec)[rec_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[rec_vec->size[1]]) rec_vec->data;
  float (*restrict rec_coords)[rec_coords_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[rec_coords_vec->size[1]]) rec_coords_vec->data;
  float (*restrict src)[src_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[src_vec->size[1]]) src_vec->data;
  float (*restrict src_coords)[src_coords_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[src_coords_vec->size[1]]) src_coords_vec->data;
  float (*restrict u)[u_vec->size[1]][u_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[u_vec->size[1]][u_vec->size[2]]) u_vec->data;
  float (*restrict vp)[vp_vec->size[1]] __attribute__ ((aligned (64))) = (float (*)[vp_vec->size[1]]) vp_vec->data;
  #pragma acc enter data copyin(rec[0:rec_vec->size[0]][0:rec_vec->size[1]])
  #pragma acc enter data copyin( u[0:u_vec->size[0]][0:u_vec->size[1]][0:u_vec->size[2]])
  #pragma acc enter data copyin( damp[0:damp_vec->size[0]][0:damp_vec->size[1]])
  #pragma acc enter data copyin( rec_coords[0:rec_coords_vec->size[0]][0:rec_coords_vec->size[1]])
  #pragma acc enter data copyin( src[0:src_vec->size[0]][0:src_vec->size[1]])
  #pragma acc enter data copyin( src_coords[0:src_coords_vec->size[0]][0:src_coords_vec->size[1]])
  #pragma acc enter data copyin( vp[0:vp_vec->size[0]][0:vp_vec->size[1]])
  for (int time = time_m, t0 = (time)%(3), t1 = (time + 1)%(3), t2 = (time + 2)%(3); time <= time_M; time += 1, t0 = (time)%(3), t1 = (time + 1)%(3), t2 = (time + 2)%(3))
  {
    struct timeval start_section0, end_section0;
    gettimeofday(&start_section0, NULL);
    /* Begin section0 */
    #pragma acc parallel loop collapse(2)
    for (int x = x_m; x <= x_M; x += 1)
    {
      for (int y = y_m; y <= y_M; y += 1)
      {
        float r0 = vp[x + 6][y + 6]*vp[x + 6][y + 6];
        u[t1][x + 6][y + 6] = 2.0F*(5.0e-1F*r0*(dt*dt)*(4.93827172e-5F*(u[t0][x + 3][y + 6] + u[t0][x + 6][y + 3] + u[t0][x + 6][y + 9] + u[t0][x + 9][y + 6]) - 6.66666683e-4F*(u[t0][x + 4][y + 6] + u[t0][x + 6][y + 4] + u[t0][x + 6][y + 8] + u[t0][x + 8][y + 6]) + 6.66666683e-3F*(u[t0][x + 5][y + 6] + u[t0][x + 6][y + 5] + u[t0][x + 6][y + 7] + u[t0][x + 7][y + 6]) - 2.41975314e-2F*u[t0][x + 6][y + 6]) + 5.0e-1F*(r0*dt*damp[x + 1][y + 1]*u[t0][x + 6][y + 6] - u[t2][x + 6][y + 6]) + 1.0F*u[t0][x + 6][y + 6])/(r0*dt*damp[x + 1][y + 1] + 1);
      }
    }
    /* End section0 */
    gettimeofday(&end_section0, NULL);
    timers->section0 += (double)(end_section0.tv_sec-start_section0.tv_sec)+(double)(end_section0.tv_usec-start_section0.tv_usec)/1000000;
    struct timeval start_section1, end_section1;
    gettimeofday(&start_section1, NULL);
    /* Begin section1 */
    #pragma acc parallel loop collapse(1)
    for (int p_src = p_src_m; p_src <= p_src_M; p_src += 1)
    {
      int ii_src_0 = (int)(floor(-6.66667e-2*o_x + 6.66667e-2*src_coords[p_src][0]));
      int ii_src_1 = (int)(floor(-6.66667e-2*o_y + 6.66667e-2*src_coords[p_src][1]));
      int ii_src_2 = (int)(floor(-6.66667e-2*o_y + 6.66667e-2*src_coords[p_src][1])) + 1;
      int ii_src_3 = (int)(floor(-6.66667e-2*o_x + 6.66667e-2*src_coords[p_src][0])) + 1;
      float px = (float)(-o_x - 1.5e+1F*(int)(floor(-6.66667e-2F*o_x + 6.66667e-2F*src_coords[p_src][0])) + src_coords[p_src][0]);
      float py = (float)(-o_y - 1.5e+1F*(int)(floor(-6.66667e-2F*o_y + 6.66667e-2F*src_coords[p_src][1])) + src_coords[p_src][1]);
      if (ii_src_0 >= x_m - 1 && ii_src_1 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_1 <= y_M + 1)
      {
        float r1 = (dt*dt)*(vp[ii_src_0 + 6][ii_src_1 + 6]*vp[ii_src_0 + 6][ii_src_1 + 6])*(4.44445e-3F*px*py - 6.66667e-2F*px - 6.66667e-2F*py + 1)*src[time][p_src];
        #pragma acc atomic update
        u[t1][ii_src_0 + 6][ii_src_1 + 6] += r1;
      }
      if (ii_src_0 >= x_m - 1 && ii_src_2 >= y_m - 1 && ii_src_0 <= x_M + 1 && ii_src_2 <= y_M + 1)
      {
        float r2 = (dt*dt)*(vp[ii_src_0 + 6][ii_src_2 + 6]*vp[ii_src_0 + 6][ii_src_2 + 6])*(-4.44445e-3F*px*py + 6.66667e-2F*py)*src[time][p_src];
        #pragma acc atomic update
        u[t1][ii_src_0 + 6][ii_src_2 + 6] += r2;
      }
      if (ii_src_1 >= y_m - 1 && ii_src_3 >= x_m - 1 && ii_src_1 <= y_M + 1 && ii_src_3 <= x_M + 1)
      {
        float r3 = (dt*dt)*(vp[ii_src_3 + 6][ii_src_1 + 6]*vp[ii_src_3 + 6][ii_src_1 + 6])*(-4.44445e-3F*px*py + 6.66667e-2F*px)*src[time][p_src];
        #pragma acc atomic update
        u[t1][ii_src_3 + 6][ii_src_1 + 6] += r3;
      }
      if (ii_src_2 >= y_m - 1 && ii_src_3 >= x_m - 1 && ii_src_2 <= y_M + 1 && ii_src_3 <= x_M + 1)
      {
        float r4 = 4.44445e-3F*px*py*(dt*dt)*(vp[ii_src_3 + 6][ii_src_2 + 6]*vp[ii_src_3 + 6][ii_src_2 + 6])*src[time][p_src];
        #pragma acc atomic update
        u[t1][ii_src_3 + 6][ii_src_2 + 6] += r4;
      }
    }
    /* End section1 */
    gettimeofday(&end_section1, NULL);
    timers->section1 += (double)(end_section1.tv_sec-start_section1.tv_sec)+(double)(end_section1.tv_usec-start_section1.tv_usec)/1000000;
    struct timeval start_section2, end_section2;
    gettimeofday(&start_section2, NULL);
    /* Begin section2 */
    #pragma acc parallel loop collapse(1)
    for (int p_rec = p_rec_m; p_rec <= p_rec_M; p_rec += 1)
    {
      int ii_rec_0 = (int)(floor(-6.66667e-2*o_x + 6.66667e-2*rec_coords[p_rec][0]));
      int ii_rec_1 = (int)(floor(-6.66667e-2*o_y + 6.66667e-2*rec_coords[p_rec][1]));
      int ii_rec_2 = (int)(floor(-6.66667e-2*o_y + 6.66667e-2*rec_coords[p_rec][1])) + 1;
      int ii_rec_3 = (int)(floor(-6.66667e-2*o_x + 6.66667e-2*rec_coords[p_rec][0])) + 1;
      float px = (float)(-o_x - 1.5e+1F*(int)(floor(-6.66667e-2F*o_x + 6.66667e-2F*rec_coords[p_rec][0])) + rec_coords[p_rec][0]);
      float py = (float)(-o_y - 1.5e+1F*(int)(floor(-6.66667e-2F*o_y + 6.66667e-2F*rec_coords[p_rec][1])) + rec_coords[p_rec][1]);
      float sum = 0.0F;
      if (ii_rec_0 >= x_m - 1 && ii_rec_1 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_1 <= y_M + 1)
      {
        sum += (4.44445e-3F*px*py - 6.66667e-2F*px - 6.66667e-2F*py + 1)*u[t0][ii_rec_0 + 6][ii_rec_1 + 6];
      }
      if (ii_rec_0 >= x_m - 1 && ii_rec_2 >= y_m - 1 && ii_rec_0 <= x_M + 1 && ii_rec_2 <= y_M + 1)
      {
        sum += (-4.44445e-3F*px*py + 6.66667e-2F*py)*u[t0][ii_rec_0 + 6][ii_rec_2 + 6];
      }
      if (ii_rec_1 >= y_m - 1 && ii_rec_3 >= x_m - 1 && ii_rec_1 <= y_M + 1 && ii_rec_3 <= x_M + 1)
      {
        sum += (-4.44445e-3F*px*py + 6.66667e-2F*px)*u[t0][ii_rec_3 + 6][ii_rec_1 + 6];
      }
      if (ii_rec_2 >= y_m - 1 && ii_rec_3 >= x_m - 1 && ii_rec_2 <= y_M + 1 && ii_rec_3 <= x_M + 1)
      {
        sum += 4.44445e-3F*px*py*u[t0][ii_rec_3 + 6][ii_rec_2 + 6];
      }
      rec[time][p_rec] = sum;
    }
    /* End section2 */
    gettimeofday(&end_section2, NULL);
    timers->section2 += (double)(end_section2.tv_sec-start_section2.tv_sec)+(double)(end_section2.tv_usec-start_section2.tv_usec)/1000000;
  }
  #pragma acc exit data copyout(rec[0:rec_vec->size[0]][0:rec_vec->size[1]])
  #pragma acc exit data copyout(u[0:u_vec->size[0]][0:u_vec->size[1]][0:u_vec->size[2]])
  #pragma acc exit data copyout(damp[0:damp_vec->size[0]][0:damp_vec->size[1]])
  #pragma acc exit data copyout(rec_coords[0:rec_coords_vec->size[0]][0:rec_coords_vec->size[1]])
  #pragma acc exit data copyout(src[0:src_vec->size[0]][0:src_vec->size[1]])
  #pragma acc exit data copyout(src_coords[0:src_coords_vec->size[0]][0:src_coords_vec->size[1]])
  #pragma acc exit data copyout(vp[0:vp_vec->size[0]][0:vp_vec->size[1]])
  return 0;
}
/* Backdoor edit at Tue Mar  3 21:02:53 2020*/ 
/* Backdoor edit at Tue Mar  3 21:04:27 2020*/ 
/* Backdoor edit at Tue Mar  3 21:05:28 2020*/ 
/* Backdoor edit at Tue Mar  3 21:06:30 2020*/ 
/* Backdoor edit at Tue Mar  3 21:12:05 2020*/ 
/* Backdoor edit at Tue Mar  3 21:12:15 2020*/ 
/* Backdoor edit at Tue Mar  3 21:16:25 2020*/ 
/* Backdoor edit at Tue Mar  3 21:16:39 2020*/ 
/* Backdoor edit at Tue Mar  3 21:17:51 2020*/ 
/* Backdoor edit at Tue Mar  3 21:23:58 2020*/ 
/* Backdoor edit at Tue Mar  3 21:24:10 2020*/ 
/* Backdoor edit at Tue Mar  3 21:24:46 2020*/ 
/* Backdoor edit at Tue Mar  3 21:27:53 2020*/ 
/* Backdoor edit at Tue Mar  3 21:28:19 2020*/ 
/* Backdoor edit at Tue Mar  3 21:29:24 2020*/ 
/* Backdoor edit at Tue Mar  3 21:32:27 2020*/ 
/* Backdoor edit at Tue Mar  3 21:35:02 2020*/ 
