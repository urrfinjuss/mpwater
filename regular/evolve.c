#include "header.h"


static fftwl_complex 	**kq, *tQ;
static fftwl_complex 	**kv, *tV;

static long double	cfl = 0.200L; // used to be 0.080L
//static unsigned long    kz;
static const unsigned long 	pD = 12;
static const long double 	one_third        = 1.0L/3.0L;
static const long double 	two_thirds       = 2.0L/3.0L;
static const long double 	one_twelfth      = 1.0L/12.0L;
static const long double 	one_eighth       = 1.0L/8.0L;
static const long double 	one_sixteenth 	 = 1.0L/16.0L;
static const long double 	one_fortyfourths = 1.0L/44.0L;
static const long double 	one_onetwentieth = 1.0L/120.0L;

void init_timemarching() {
  kq = fftwl_malloc(7*sizeof(fftwl_complex *));
  kv = fftwl_malloc(7*sizeof(fftwl_complex *));
  /*
  state.kD = lroundl((11.L/12.L)*state.number_modes/2.L);  // was 11/12
  kz = lroundl(7.L/6.L*state.kD);
  if (kz > state.number_modes/2) kz = state.number_modes/2;
  printf("Dissipative Range starts at k0 = %lu\n", state.kD);
  */
}

void allocate_timemarching() {
  unsigned long 	N = state.number_modes;
  for (long int j = 0; j < 7; j++) {
    kq[j] = fftwl_malloc(N*sizeof(fftwl_complex));
    kv[j] = fftwl_malloc(N*sizeof(fftwl_complex));
  }
  /*
  state.kD = lroundl((11.L/12.L)*state.number_modes/2.L);  // was 11/12
  kz = lroundl(7.L/6.L*state.kD);
  if (kz > state.number_modes/2) kz = state.number_modes/2;
  printf("New Dissipative Range starts at k0 = %lu\n", state.kD);
  printf("New Zero Range starts at k0 = %lu\n", kz);
  */
  tQ = fftwl_malloc(N*sizeof(fftwl_complex));
  tV = fftwl_malloc(N*sizeof(fftwl_complex));
}

void deallocate_timemarching() {
  for (long int j = 0; j < 7; j++) {
    fftwl_free(kq[j]);
    fftwl_free(kv[j]);
  }
  fftwl_free(tQ);
  fftwl_free(tV);
}

void compute_rhs(fftwl_complex *inQ, fftwl_complex *inV, fftwl_complex *outQ, fftwl_complex *outV) {
  unsigned long 	N = state.number_modes;
  long double 		overN = 1.L/N;
  long double 		sigma = state.surface_tension;
  long double 		g = state.gravity;
  fftwl_complex 	w2 = cexpl(-1.IL*conf.origin_offset - 2.L*atanhl(conf.scaling));
  fftwl_complex		w1 = conjl(w2);
  fftwl_complex 	b1U = 0.L, b2U = 0.L;
  //fftwl_complex 	b1B = 0.L, b2B = 0.L;

  if (conf.scaling == 1.0L) {
	w2 = 0.L;
	w1 = 0.L;
  }
  memcpy(tmpc[0], inQ, N*sizeof(fftwl_complex));
  memcpy(tmpc[1], inV, N*sizeof(fftwl_complex));
  fftwl_execute(ift0); 
  fftwl_execute(ift1); 
  memset(tmpc[0]+N/2, 0, N/2*sizeof(fftwl_complex));
  memset(tmpc[1]+N/2, 0, N/2*sizeof(fftwl_complex));
  for (long int j = 0; j < N/2; j++) {
    tmpc[0][j] = -1.IL*j*tmpc[0][j]*overN;
    tmpc[1][j] = -1.IL*j*tmpc[1][j]*overN;
  }
  fftwl_execute(ft0);
  fftwl_execute(ft1);
  for (long int j = 0; j < N; j++) {
    tmpc[2][j]= 2.L*creall(inV[j]*conjl(inQ[j]*inQ[j]))*overN;
    tmpc[3][j]= inV[j]*conjl(inV[j])+4.L*sigma*conf.dq[j]*cimagl(tmpc[0][j]*conjl(inQ[j]));
    tmpc[3][j]= tmpc[3][j]*overN;
  }
  //project(tmpc[3], tmpc[4]); 
  //complex_array_out("B.ph.old.txt", tmpc[4]);
  //for (long int j = 0; j < N; j++) {
  //  tmpc[3][j]= tmpc[3][j]*overN;
  //}
  //complex_array_out("tmpc2.txt", tmpc[2]);   // tmpc[2] checks out
  //complex_array_out("tmpc3.txt", tmpc[3]);   // tmpc[3] checks out
  fftwl_execute(ift2);
  fftwl_execute(ift3);
  tmpc[4][0] = 0.L;
  //b2B = tmpc[3][N/2-1];	b1B = tmpc[3][N/2+1];
  b2U = tmpc[2][N/2-1];	b1U = tmpc[2][N/2+1];
  for (long int j = 1; j < N/2 - 1; j++) {
    //b1B = b1B*w1 + tmpc[3][N/2+1+j];
    //b2B = b2B*w2 + tmpc[3][N/2-1-j];
    b1U = b1U*w1 + tmpc[2][N/2+1+j];
    b2U = b2U*w2 + tmpc[2][N/2-1-j];
  }
  //printf("b1U = %23.18LE\t%23.18LE\n", creall(b1U), cimagl(b1U));
  //printf("b2U = %23.18LE\t%23.18LE\n", creall(b2U), cimagl(b2U));
  //exit(1);
  memset(tmpc[2]+N/2, 0, N/2*sizeof(fftwl_complex));
  memset(tmpc[3]+N/2, 0, N/2*sizeof(fftwl_complex));
  memset(tmpc[4]+N/2, 0, N/2*sizeof(fftwl_complex));
  tmpc[2][0] = 0.5L*tmpc[2][0];
  //tmpc[3][0] = 0.5L*tmpc[3][0];  
  for (long int j = 0; j < N/2; j++) {
    tmpc[3][j] = -1.IL*j*tmpc[3][j];
    tmpc[4][j] = -1.IL*j*tmpc[2][j];
  }
  fftwl_execute(ft2);
  fftwl_execute(ft3);
  fftwl_execute(ft4);
  for (long int j = 0; j < N; j++) {
    tmpc[2][j] += 0.5L*(b1U*w1 - b2U*w2);
    //tmpc[3][j] += 0.5L*(b1B*w1 - b2B*w2); 	// thats irrelevant for B'
  }
  // Summary:
  // tmpc[0] <--  stores Q'
  // tmpc[1] <--  stores V'
  // tmpc[2] <--  stores U 
  // tmpc[3] <--  stores B'
  // tmpc[4] <--  stores U'
  for (long int j = 0; j < N; j++) {
    outQ[j] = 0.5IL*conf.dq[j]*(2.L*tmpc[0][j]*tmpc[2][j]-tmpc[4][j]*inQ[j]);
    outV[j] = g*(inQ[j]*inQ[j]-1.L)+1.IL*conf.dq[j]*(tmpc[2][j]*tmpc[1][j]-inQ[j]*inQ[j]*tmpc[3][j]);
  }
}


void rk6_step(fftwl_complex *inQ, fftwl_complex *inV, long double dt) {
  unsigned long 	N = state.number_modes;
  long double		overN = 1.L/N;

  //check for bug
  memcpy(tmpc[0], inQ, N*sizeof(fftwl_complex));
  memcpy(tmpc[1], inV, N*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  memset(tmpc[0]+N/2, 0, N/2*sizeof(fftwl_complex));
  memset(tmpc[1]+N/2, 0, N/2*sizeof(fftwl_complex));
  for (long int j = 0; j < N/2; j++) {
    tmpc[0][j] = tmpc[0][j]*overN;
    tmpc[1][j] = tmpc[1][j]*overN;
  }
  fftwl_execute(ft0);
  fftwl_execute(ft1); 
  memcpy(inQ, tmpc[0], N*sizeof(fftwl_complex));
  memcpy(inV, tmpc[1], N*sizeof(fftwl_complex));
  
  memcpy(tQ, inQ, N*sizeof(fftwl_complex)); 
  memcpy(tV, inV, N*sizeof(fftwl_complex)); 
  

  compute_rhs(tQ, tV, kq[0], kv[0]);
  // check rhs
  // end check rhs
  for (long int j = 0; j < N; j++) {
    tQ[j] = inQ[j] + one_third*dt*kq[0][j];
    tV[j] = inV[j] + one_third*dt*kv[0][j];
  }
  compute_rhs(tQ, tV, kq[1], kv[1]);
  for (long int j = 0; j < N; j++) {
    tQ[j] = inQ[j] + two_thirds*dt*kq[1][j];
    tV[j] = inV[j] + two_thirds*dt*kv[1][j];
  }
  compute_rhs(tQ, tV, kq[2], kv[2]);
  for (long int j = 0; j < N; j++) {
    tQ[j] = inQ[j] + one_twelfth*dt*(kq[0][j] + 4.0L*kq[1][j] - kq[2][j]);
    tV[j] = inV[j] + one_twelfth*dt*(kv[0][j] + 4.0L*kv[1][j] - kv[2][j]);
  }
  compute_rhs(tQ, tV, kq[3], kv[3]);
  for (long int j = 0; j < N; j++) {
    tQ[j] = inQ[j] + one_sixteenth*dt*(-kq[0][j] + 18.0L*kq[1][j] - 3.0L*kq[2][j] - 6.0L*kq[3][j]);
    tV[j] = inV[j] + one_sixteenth*dt*(-kv[0][j] + 18.0L*kv[1][j] - 3.0L*kv[2][j] - 6.0L*kv[3][j]);
  }
  compute_rhs(tQ, tV, kq[4], kv[4]);
  for (long int j = 0; j < N; j++) {
    tQ[j] = inQ[j] + one_eighth*dt*(9.0L*kq[1][j] - 3.0L*kq[2][j] - 6.0L*kq[3][j] + 4.0L*kq[4][j]);
    tV[j] = inV[j] + one_eighth*dt*(9.0L*kv[1][j] - 3.0L*kv[2][j] - 6.0L*kv[3][j] + 4.0L*kv[4][j]);
  }
  compute_rhs(tQ, tV, kq[5], kv[5]);
  for (long int j = 0; j < N; j++) {
    tQ[j] = inQ[j] + one_fortyfourths*dt*(9.0L*kq[0][j] - 36.0L*kq[1][j] + 63.0L*kq[2][j] + 72.0L*kq[3][j] - 64.0L*kq[5][j]);
    tV[j] = inV[j] + one_fortyfourths*dt*(9.0L*kv[0][j] - 36.0L*kv[1][j] + 63.0L*kv[2][j] + 72.0L*kv[3][j] - 64.0L*kv[5][j]);
  }
  compute_rhs(tQ, tV, kq[6], kv[6]);
  for (long int j = 0; j < N; j++) {
    inQ[j] = (inQ[j] + one_onetwentieth*dt*(11.0L*kq[0][j] + 81.0L*kq[2][j] + 81.0L*kq[3][j] - 32.0L*kq[4][j] - 32.0L*kq[5][j] + 11.0L*kq[6][j]))*overN;
    inV[j] = (inV[j] + one_onetwentieth*dt*(11.0L*kv[0][j] + 81.0L*kv[2][j] + 81.0L*kv[3][j] - 32.0L*kv[4][j] - 32.0L*kv[5][j] + 11.0L*kv[6][j]))*overN;
  }
  memcpy(tmpc[0], inQ, N*sizeof(fftwl_complex));
  memcpy(tmpc[1], inV, N*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  memset(tmpc[0]+N/2, 0, N/2*sizeof(fftwl_complex));
  memset(tmpc[1]+N/2, 0, N/2*sizeof(fftwl_complex));
  fftwl_execute(ft0);
  fftwl_execute(ft1); 
  
  memcpy(inQ, tmpc[0], N*sizeof(fftwl_complex));
  memcpy(inV, tmpc[1], N*sizeof(fftwl_complex));
}

void evolve_rk6() {
  unsigned int		QC_pass = 1;
  unsigned long 	counter = 0, j = 0, skip = 500;
  unsigned long		ref_counter = 0;
  char 			filename1[80];
  char 			filename2[80];
  long double		M_TOL = 5.0E-09L;
  long double		R_TOL = 1.0E-09L;
  long double		tshift = 0.L;
  long double   	time = 0.L, Ham = 0.L;
  long double   	dt = cfl*2.L*PI*conf.scaling/state.number_modes;
  FILE *fh_time	= fopen("time_dependence.txt","w");
  fprintf(fh_time, "# 1. time 2. Kinetic 3. Potential 4. Momentum X 5. Momentum Y\n\n");
  fclose(fh_time);
  map_quality_fourier(data[0], data[1], M_TOL, &QC_pass);
  sprintf(filename2, "./data/spec_%04lu.txt", counter);
  spec_out(filename2, tmpc[0], tmpc[1]);
  convertQtoZ(data[0], tmpc[5]);
  sprintf(filename1, "./data/surf_%04lu.txt", counter);
  surface_out(filename1, tmpc[5]);
  if (PADE_TEST) {
    exit(1);
  }
  restore_potential(data[0], data[1], tmpc[5]);  
  fh_time = fopen("time_dependence.txt","a");
  fprintf(fh_time, "%.17LE\t%.17LE\t%.17LE\t", state.time, state.kineticE/PI, state.potentialE/PI); 
  fprintf(fh_time, "%.17LE\t%.17LE\n", cimagl(state.momentum), creall(state.momentum)); 
  fclose(fh_time);
  Ham = (state.kineticE + state.potentialE)/PI;
  sprintf(filename1, "./aux/data_%04lu.txt", counter);
  output_data(filename1, tmpc[5]);
  printf("T = %23.16LE\tH = %23.16LE\n", state.time, Ham);

  //*sqrtl(state.number_modes/4096.L)
  if (QC_pass == 0) {
    printf("Bad quality map at start.\tStop!\n");
    exit(1);
  } else if (QC_pass == 2) {
    printf("Warning! Map is over-resolved at start.\n");
  }
  while (QC_pass) {
    rk6_step(data[0], data[1], dt);  
    time = (j+1)*dt;
    j++;
    state.time = time + tshift;

    map_quality_fourier(data[0], data[1], M_TOL, &QC_pass); 
    if (QC_pass == 0) {
      ref_counter++;

      // Attempt to fix the conformal map
      restore_potential(data[0], data[1], tmpc[2]);
      // Control Map 1: try to shift zoom location and L-scaling (DO NOT REMOVE!)
      track_singularity(data[0]);
      alt_map.scaling = conf.scaling/sqrt(2.L);
      // End Control Map 1

      remap(&alt_map, state.number_modes);
      restore_potential(data[0], data[1], tmpc[2]);
      print_constants();

      map_quality_fourier(data[0], data[1], R_TOL, &QC_pass); 
      spec_out("spec_after.txt", tmpc[0], tmpc[1]);
      if (QC_pass == 0) {
        printf("Doubling # of Modes: %lu\n", 2*state.number_modes);
        // Control Map 2: try to shift zoom location and L-scaling and double N(DO NOT REMOVE!)
        track_singularity(data[0]);
        alt_map.scaling = conf.scaling*sqrtl(2.L);
        // End Control Map 2
        remap(&alt_map, 2*state.number_modes);
        skip = lroundl(1.5L*skip);
        restore_potential(data[0], data[1], tmpc[2]);
        print_constants();
        map_quality_fourier(data[0], data[1], R_TOL, &QC_pass); 
      }
      if (QC_pass == 1) {
        restore_potential(data[0], data[1], tmpc[2]);  
        Ham = (state.kineticE + state.potentialE)/PI;
        printf("T = %23.16LE\tH = %23.16LE\n", state.time, Ham);
        tshift += time;
        j = 0;
        cfl = 0.98L*cfl;
     	dt = cfl*2.L*PI*conf.scaling/state.number_modes;
        skip = lroundl(sqrtl(1.5L)*skip);
      } else {
	printf("Failed to find a good map!\n");
	exit(1);
      }
    } else {
      if ( !((j+1) % skip) ) {
        counter++;
        // write out spectrum
        sprintf(filename1, "./data/spec_%04lu.txt", counter);
        spec_out(filename1, tmpc[0], tmpc[1]);
        // write out surface shape and cut for Z
        convertQtoZ(data[0], tmpc[5]);
        /*sprintf(filename2, "./roots/roots_Z%04lu.txt", counter);
        optimal_pade(filename2, tmpc[5]);*/
        sprintf(filename1, "./data/surf_%04lu.txt", counter);
        surface_out(filename1, tmpc[5]);
        // write out potential and its cut
        restore_potential(data[0], data[1], tmpc[5]);  
        fh_time = fopen("time_dependence.txt","a");
        fprintf(fh_time, "%.17LE\t%.17LE\t%.17LE\t", state.time, state.kineticE/PI, state.potentialE/PI); 
        fprintf(fh_time, "%.17LE\t%.17LE\n", cimagl(state.momentum), creall(state.momentum)); 
        fclose(fh_time);
        Ham = (state.kineticE + state.potentialE)/PI;
        printf("T = %23.16LE\tH = %23.16LE\n", state.time, Ham);
        sprintf(filename1, "./aux/data_%04lu.txt", counter);
        output_data(filename1, tmpc[5]);
	/*
        sprintf(filename2, "./roots/roots_P%04lu.txt", counter);
        optimal_pade(filename2, tmpc[5]);
	*/
        // write out cut for derivative of potential
        /*
        memcpy(tmpc[0], tmpc[5], state.number_modes*sizeof(fftwl_complex));        
        fftwl_execute(ift0);
        memset(tmpc[0] + state.number_modes/2, 0, state.number_modes*sizeof(fftwl_complex)/2);
        for (unsigned long j = 0; j < state.number_modes/2; j++) tmpc[0][j] = -1.0IL*j*tmpc[0][j]/state.number_modes; // this is derivative dq
        fftwl_execute(ft0);
        for (unsigned long j = 0; j < state.number_modes; j++) tmpc[0][j] = tmpc[0][j]*conf.dq[j]; // this is derivative du
        sprintf(filename2, "./roots/roots_dP%04lu.txt", counter);
        optimal_pade(filename2, tmpc[0]);
        */
       
	// write out cut for R
        /*for (unsigned long j = 0; j < state.number_modes; j++) tmpc[5][j] = data[0][j]*data[0][j];
        sprintf(filename2, "./roots/roots_R%04lu.txt", counter);
        optimal_pade(filename2, tmpc[5]);*/

	// write out cut for Zu
	/*
        for (unsigned long j = 0; j < state.number_modes; j++) tmpc[5][j] = 1.L/(data[0][j]*data[0][j]);
        sprintf(filename2, "./roots/roots_DZ%04lu.txt", counter);
        optimal_pade(filename2, tmpc[5]);
        */
	// write out cut for V
        /*for (unsigned long j = 0; j < state.number_modes; j++) tmpc[5][j] = data[1][j];
        sprintf(filename2, "./roots/roots_V%04lu.txt", counter);
        optimal_pade(filename2, tmpc[5]);
	*/
      }
    }
  }
  printf("Simulation Stops\n");
  exit(1);

}

