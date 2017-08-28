#include "header.h"


static params_ptr 	input;
static consq_ptr	motion;
static aux_ptr 		extra;
static data_ptr 	array;
static dirk4_data	dirk4;

static long double c1; 
static long double c2p; 
static long double c2m; 
static long double n1q, n1v, n2q, n2v;
static long double n_res;
//static fftwl_complex 	q;

void init_dirk4_module(params_ptr ginput, aux_ptr gextra, data_ptr garray, consq_ptr gmotion) {
  c1 = 0.25L;
  c2m = 0.25L - sqrtl(3.L)/6.L;
  c2p = 0.25L + sqrtl(3.L)/6.L;
  input = ginput; 
  extra = gextra; 
  array = garray; 
  motion = gmotion; 
  if (input->s) {
    dirk4.nskip = 256;
    dirk4.dt = fminl(0.25L*PI*powl(input->L/input->N,1), 0.125L*powl(input->L/input->N, 1.5L))/sqrtl(input->s);
  } else {
    dirk4.nskip = 256; 
    dirk4.dt = 1.2L*sqrtl(input->L/0.03125)*sqrtl(1024.Q/input->N)*PI*powl(input->L/input->N,1);      // overturning Pavel (no st)
  }
  dirk4.D = 0.25L*powl(92.L*input->N/256.L,-12)*dirk4.dt;  // adjusted to 0.25L from 1.0L
  dirk4.D = 1.00L*powl(104.L*input->N/256.L,-12)*dirk4.dt;  // adjusted to 0.25L from 1.0L
  dirk4.D = 0.0;
  dirk4.tshift = 0.L;
  dirk4.fp_tolerance = 2.0E-15L;
  printf("%ld\n", dirk4.nskip);
}

long double norm_l2(fftwl_complex *in1, fftwl_complex *in2) {
    long double rvalue = 0.;
    for (long int j = 0; j < input->N; j++) {
       rvalue += (in1[j] - in2[j])*conj(in1[j] - in2[j]);
    }
    return sqrtl(rvalue)/input->N;
}


long double init_dirk4() { 
   if ( input->s > 0.0L) {
     dirk4.dt = fminl(0.25L*PI*powl(input->L/input->N,1), 0.125L*powl(input->L/input->N, 1.5L))/sqrtl(input->s);
   } else {
     dirk4.dt = 1.2L*sqrtl(input->L/0.03125)*sqrtl(1024.Q/input->N)*PI*powl(input->L/input->N,1);      // overturning Pavel (no st) modified sqrt(L/L0)
   }
   dirk4.tshift += dirk4.time;
   dirk4.time = 0.L;
   dirk4.nsteps = 100000;

   dirk4.D = 0.25L*powl(92.L*input->N/256.L,-12)*dirk4.dt;  // adjusted to 0.25L from 1.0L /input->L
   dirk4.D = 1.00L*powl(104.L*input->N/256.L,-12)*dirk4.dt;  // adjusted to 0.25L from 1.0L /input->L
   dirk4.D = 0.0;
   dirk4.max_iter = 25;   

   if (!(dirk4.tmpq1 = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.tmpq2 = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.k1q = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.k2q = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.k1q_prev = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.k2q_prev = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.tmpv1 = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.tmpv2 = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.k1v = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.k2v = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.k1v_prev = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(dirk4.k2v_prev = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);


   return dirk4.D;
}


void dirk4_step() {
  memcpy(dirk4.tmpq1, array->Q, input->N*sizeof(fftwl_complex));
  memcpy(dirk4.tmpv1, array->V, input->N*sizeof(fftwl_complex));

  compute_aux_arrays(array->Q, array->V);

  for (long int j = 0; j < input->N; j++) {
    dirk4.k1q[j] = 0.5IL*(2.0L*extra->dQ[j]*extra->U[j] - extra->dU[j]*array->Q[j])*extra->newdQ[j];
    dirk4.k1v[j] = 1.0IL*(extra->U[j]*extra->dV[j] - array->Q[j]*array->Q[j]*extra->B[j])*extra->newdQ[j] + input->g*(array->Q[j]*array->Q[j] - 1.0L);
  }
  memcpy(dirk4.k2q, dirk4.k1q, input->N*sizeof(fftwl_complex));
  memcpy(dirk4.k2v, dirk4.k1v, input->N*sizeof(fftwl_complex));

  long int current_iter = 0; 
  n1q = 1.L; n1v = 1.L;
  n2q = 1.L; n2v = 1.L;

  //printf("Starting Fixed Point Iteration (Implicit RK4):\n");
  //while (current_iter < dirk4.max_iter) {
  while (1) {
    memcpy(dirk4.k1q_prev, dirk4.k1q, input->N*sizeof(fftwl_complex));
    memcpy(dirk4.k1v_prev, dirk4.k1v, input->N*sizeof(fftwl_complex));
    memcpy(dirk4.k2q_prev, dirk4.k2q, input->N*sizeof(fftwl_complex));
    memcpy(dirk4.k2v_prev, dirk4.k2v, input->N*sizeof(fftwl_complex));
    for (long int j = 0; j < input->N; j++) {
      dirk4.tmpq1[j] = array->Q[j] + c1*dirk4.dt*dirk4.k1q[j] + c2m*dirk4.dt*dirk4.k2q[j];
      dirk4.tmpv1[j] = array->V[j] + c1*dirk4.dt*dirk4.k1v[j] + c2m*dirk4.dt*dirk4.k2v[j];

      dirk4.tmpq2[j] = array->Q[j] + c2p*dirk4.dt*dirk4.k1q[j] + c1*dirk4.dt*dirk4.k2q[j];
      dirk4.tmpv2[j] = array->V[j] + c2p*dirk4.dt*dirk4.k1v[j] + c1*dirk4.dt*dirk4.k2v[j];
    }
    compute_aux_arrays(dirk4.tmpq1, dirk4.tmpv1);
    for (long int j = 0; j < input->N; j++) {
      dirk4.k1q[j] = 0.5IL*(2.0L*extra->dQ[j]*extra->U[j] - extra->dU[j]*dirk4.tmpq1[j])*extra->newdQ[j];
      dirk4.k1v[j] = 1.0IL*(extra->U[j]*extra->dV[j] - dirk4.tmpq1[j]*dirk4.tmpq1[j]*extra->B[j])*extra->newdQ[j] + input->g*(dirk4.tmpq1[j]*dirk4.tmpq1[j] - 1.0L);
    }
    compute_aux_arrays(dirk4.tmpq2, dirk4.tmpv2);
    for (long int j = 0; j < input->N; j++) {
      dirk4.k2q[j] = 0.5IL*(2.0L*extra->dQ[j]*extra->U[j] - extra->dU[j]*dirk4.tmpq2[j])*extra->newdQ[j];
      dirk4.k2v[j] = 1.0IL*(extra->U[j]*extra->dV[j] - dirk4.tmpq2[j]*dirk4.tmpq2[j]*extra->B[j])*extra->newdQ[j] + input->g*(dirk4.tmpq2[j]*dirk4.tmpq2[j] - 1.0L);
    }
    current_iter++;
    /*printf("Current Iteration %ld of %ld\n", current_iter, dirk4.max_iter);
    printf("|k1q^n+1 - k1q^n| = %.12LE\t|k1v^n+1 - k1v^n| = %.12LE\n", norm_l2(dirk4.k1q, dirk4.k1q_prev), norm_l2(dirk4.k1v, dirk4.k1v_prev));
    printf("|k2q^n+1 - k2q^n| = %.12LE\t|k2v^n+1 - k2v^n| = %.12LE\n", norm_l2(dirk4.k2q, dirk4.k2q_prev), norm_l2(dirk4.k2v, dirk4.k2v_prev));*/
    n1q = norm_l2(dirk4.k1q, dirk4.k1q_prev); 
    n1v = norm_l2(dirk4.k1v, dirk4.k1v_prev);
    n2q = norm_l2(dirk4.k2q, dirk4.k2q_prev); 
    n2v = norm_l2(dirk4.k2v, dirk4.k2v_prev);
    n_res = fmaxl(fmaxl(n1q,n2q),fmaxl(n1v,n2v));
    if ( n_res < dirk4.fp_tolerance ) {
        //printf("Converged to %.12LE at iteration %ld\n", n_res, current_iter);
 	break;
    }
  }
  for (long int j = 0; j < input->N; j++) {
    array->Q[j] += 0.5L*dirk4.dt*(dirk4.k1q[j] + dirk4.k2q[j]); 
    array->V[j] += 0.5L*dirk4.dt*(dirk4.k1v[j] + dirk4.k2v[j]); 
  }
  hfilter();
}

void free_dirk4() {
  fftwl_free(dirk4.tmpq1);	  fftwl_free(dirk4.tmpv1);
  fftwl_free(dirk4.tmpq2);	  fftwl_free(dirk4.tmpv2);
  fftwl_free(dirk4.k1q);	  fftwl_free(dirk4.k1v);
  fftwl_free(dirk4.k2q);	  fftwl_free(dirk4.k2v);
  fftwl_free(dirk4.k1q_prev);	  fftwl_free(dirk4.k1v_prev);
  fftwl_free(dirk4.k2q_prev);	  fftwl_free(dirk4.k2v_prev);
}

void modify_dirk4skip(long double factor) {
   dirk4.nskip = floor(dirk4.nskip*factor);
}

void evolve_dirk4() {
  FILE *fh = fopen("time_dirk4.log","w");
  char fn1[256], fn2[256];
  long int l = 0;
  long int j = 0;

  dirk4.time = 0.L;
  dirk4.tshift = 0.L;
  input->refN = 0;
  get_integrals();
  sprintf(fn1, "data/data_%04ld.txt", l);
  sprintf(fn2, "spec_r%ld.txt", input->refN);
  get_surface(fn1);
  get_spectrum(fn2);

  
  fprintf(fh, "# 1. t 2. K 3. P 4. H 5. C 6. y0 7. ML\n\n");
  fprintf(fh, "% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\n", motion->time, motion->K, motion->P, motion->H, motion->C, motion->y0, motion->ml);
  fclose(fh);

  check_resolution();
  //spec_output(fname, array->Q, motion->time);

  while (1) {
    dirk4_step();
    dirk4.time = (j+1)*dirk4.dt;
    motion->time = dirk4.tshift + dirk4.time; 
    if (j % dirk4.nskip == 0) {
      l++;
      //printf("dt = %.12LE\t D = %.12LE\n", rk4.dt, rk4.D);
      get_integrals();
      sprintf(fn1, "data/data_%04ld.txt", l);
      sprintf(fn2, "data/spec_%04ld.txt", l);
      get_surface(fn1);
      get_spectrum(fn2);
      if (input->g) {
        printf("T = % .8LE\tML = % .12LE\tP = % .12LE\tK = % .12LE\tH = % .12LE\n", motion->time, motion->ml, motion->P, motion->K, motion->H);
      } else {
        printf("T = % .8LE\tML = % .12LE\tK = % .12LE\tH = % .12LE\n", motion->time, motion->ml, motion->K, motion->H);
      }
      fh = fopen("time_dirk4.log","a");
      fprintf(fh, "% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\n", motion->time, motion->K, motion->P, motion->H, motion->C, motion->y0, motion->ml);
      fclose(fh);
    }
    j++;
    resolution_monitor(&j); 

  }
  get_integrals();
  sprintf(fn1, "data/data_%04ld.txt", l+1);
  get_surface(fn1);

  if (input->g) {
    printf("T = % .8LE\tML = % .12LE\tP = % .12LE\tK = % .12LE\tH = % .12LE\n", motion->time, motion->ml, motion->P, motion->K, motion->H);
  } else {
    printf("T = % .8LE\tML = % .12LE\tK = % .12LE\tH = % .12LE\n", motion->time, motion->ml, motion->K, motion->H);
  }
  fh = fopen("time_dirk4.log","a");
  fprintf(fh, "% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\n", motion->time, motion->K, motion->P, motion->H, motion->C);
  fclose(fh);
}

