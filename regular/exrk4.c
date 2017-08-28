#include "header.h"

static rk4_data 	 rk4;
static params_ptr 	 input;
static aux_ptr		 extra;
static data_ptr		 array;
static consq_ptr	 motion;

void init_exrk4_module(params_ptr ginput, aux_ptr gextra, data_ptr garray, consq_ptr gmotion) {
  input = ginput; 
  extra = gextra; 
  array = garray; 
  motion = gmotion; 
  if (input->s) {
    rk4.nskip = 1024;
  } else {
    rk4.nskip = 1024; 
  }
  printf("%ld\n", rk4.nskip);
}


long double init_rk4(long int nskip) {
   if ( input->s ) {
     rk4.dt = fminl(0.25L*PI*powl(input->L/input->N,1), 0.25L*powl(input->L/input->N, 1.5L)/sqrtl(input->s));
   } else {
     rk4.dt = 0.25L*PI*fminf(1./input->L, input->L)/input->N;
   }
   rk4.tshift += rk4.time;
   rk4.time = 0.L;
   rk4.nsteps = 100000;
   rk4.D = 0.25L*powl(95.L*input->N/256.L,-12)*rk4.dt;  // adjusted to 0.25L from 1.0L
   rk4.D = 0.00L*powl(95.L*input->N/256.L,-12)*rk4.dt;
   motion->time = rk4.tshift;
   
   if (!(rk4.tmpq = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k0q = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k1q = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k2q = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k3q = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);

   if (!(rk4.tmpv = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k0v = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k1v = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k2v = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k3v = fftwl_malloc(input->N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   
   return rk4.D;
}

void free_rk4() {
  fftwl_free(rk4.tmpq);	  fftwl_free(rk4.tmpv);
  fftwl_free(rk4.k0q);	  fftwl_free(rk4.k0v);
  fftwl_free(rk4.k1q);	  fftwl_free(rk4.k1v);
  fftwl_free(rk4.k2q);	  fftwl_free(rk4.k2v);
  fftwl_free(rk4.k3q);	  fftwl_free(rk4.k3v);

}

void rk4_step() {
  memcpy(rk4.tmpq, array->Q, input->N*sizeof(fftwl_complex));
  memcpy(rk4.tmpv, array->V, input->N*sizeof(fftwl_complex));

  compute_aux_arrays(rk4.tmpq, rk4.tmpv);
  for (long int j = 0; j < input->N; j++) {
    rk4.k0q[j] = 0.5IL*(2.0L*extra->dQ[j]*extra->U[j] - extra->dU[j]*rk4.tmpq[j])*extra->newdQ[j];
    rk4.k0v[j] = 1.0IL*(extra->U[j]*extra->dV[j] - rk4.tmpq[j]*rk4.tmpq[j]*extra->B[j])*extra->newdQ[j] + input->g*(rk4.tmpq[j]*rk4.tmpq[j] - 1.0L);
    rk4.tmpq[j] = array->Q[j] + 0.5L*rk4.dt*rk4.k0q[j];
    rk4.tmpv[j] = array->V[j] + 0.5L*rk4.dt*rk4.k0v[j];
  }

  compute_aux_arrays(rk4.tmpq, rk4.tmpv);
  for (long int j = 0; j < input->N; j++) {
    rk4.k1q[j] = 0.5IL*(2.0L*extra->dQ[j]*extra->U[j] - extra->dU[j]*rk4.tmpq[j])*extra->newdQ[j];
    rk4.k1v[j] = 1.0IL*(extra->U[j]*extra->dV[j] - rk4.tmpq[j]*rk4.tmpq[j]*extra->B[j])*extra->newdQ[j] + input->g*(rk4.tmpq[j]*rk4.tmpq[j] - 1.0L);
    rk4.tmpq[j] = array->Q[j] + 0.5L*rk4.dt*rk4.k1q[j];
    rk4.tmpv[j] = array->V[j] + 0.5L*rk4.dt*rk4.k1v[j];
  }

  compute_aux_arrays(rk4.tmpq, rk4.tmpv);
  for (long int j = 0; j < input->N; j++) {
    rk4.k2q[j] = 0.5IL*(2.0L*extra->dQ[j]*extra->U[j] - extra->dU[j]*rk4.tmpq[j])*extra->newdQ[j];
    rk4.k2v[j] = 1.0IL*(extra->U[j]*extra->dV[j] - rk4.tmpq[j]*rk4.tmpq[j]*extra->B[j])*extra->newdQ[j] + input->g*(rk4.tmpq[j]*rk4.tmpq[j] - 1.0L);
    rk4.tmpq[j] = array->Q[j] + 1.0L*rk4.dt*rk4.k2q[j];
    rk4.tmpv[j] = array->V[j] + 1.0L*rk4.dt*rk4.k2v[j];
  }

  compute_aux_arrays(rk4.tmpq, rk4.tmpv);
  for (long int j = 0; j < input->N; j++) {
    rk4.k3q[j] = 0.5IL*(2.0L*extra->dQ[j]*extra->U[j] - extra->dU[j]*rk4.tmpq[j])*extra->newdQ[j];
    rk4.k3v[j] = 1.0IL*(extra->U[j]*extra->dV[j] - rk4.tmpq[j]*rk4.tmpq[j]*extra->B[j])*extra->newdQ[j] + input->g*(rk4.tmpq[j]*rk4.tmpq[j] - 1.0L);
    array->Q[j] += 1.0L*rk4.dt*(rk4.k0q[j] + 2.0L*rk4.k1q[j] + 2.0L*rk4.k2q[j] + rk4.k3q[j])/6.0L;
    array->V[j] += 1.0L*rk4.dt*(rk4.k0v[j] + 2.0L*rk4.k1v[j] + 2.0L*rk4.k2v[j] + rk4.k3v[j])/6.0L;
  }
  hfilter();
}

void modify_rk4skip(long double factor) {
   rk4.nskip = floor(rk4.nskip*factor);
}


void evolve_rk4() {
  FILE *fh = fopen("time_rk4.log","w");
  char fn1[256], fn2[256];
  long int l = 0;
  long int j = 0;

  rk4.time = 0.L;
  rk4.tshift = 0.L;
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
    rk4_step();
    rk4.time = (j+1)*rk4.dt;
    motion->time = rk4.tshift + rk4.time; 
    if (j % rk4.nskip == 0) {
      l++;
      //printf("dt = %.12LE\t D = %.12LE\n", rk4.dt, rk4.D);
      get_integrals();
      sprintf(fn1, "data/data_%04ld.txt", l);
      sprintf(fn2, "data/spec_%04ld.txt", l);
      get_surface(fn1);
      get_spectrum(fn2);

      printf("T = % .8LE\tML = % .12LE\tP = % .12LE\tK = % .12LE\tH = % .12LE\n", motion->time, motion->ml, motion->P, motion->K, motion->H);
      fh = fopen("time_rk4.log","a");
      fprintf(fh, "% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\n", motion->time, motion->K, motion->P, motion->H, motion->C, motion->y0, motion->ml);
      fclose(fh);
    }
    j++;
    resolution_monitor(&j); 

  }
  get_integrals();
  sprintf(fn1, "data/data_%04ld.txt", l+1);
  get_surface(fn1);

  printf("T = % .8LE\tML = % .12LE\tP = % .12LE\tK = % .12LE\tH = % .12LE\n", motion->time, motion->ml, motion->P, motion->K, motion->H);
  fh = fopen("time_rk4.log","a");
  fprintf(fh, "% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\n", motion->time, motion->K, motion->P, motion->H, motion->C);
  fclose(fh);
}

