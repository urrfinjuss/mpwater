#include "header.h"

static fftwl_complex *saveQ, *saveV;

static long double full_sum1;
static long double full_sum2;
static long double narrow_sum1; 
static long double narrow_sum2; 
static long double partial_sum1;
static long double partial_sum2;
static long double left_sum1;
static long double left_sum2;

void set_mapping() {
  long double overN = 1.0L/state.number_modes;
  long double a = (1.0L - powl(conf.scaling, 2.L))/(1.0L + powl(conf.scaling, 2.L));
  long double b = 0.5L*(1.0L + powl(conf.scaling, 2.L))/conf.scaling;

  conf.origin_offset = 2.0L*atan2l(conf.scaling*sinl(0.5L*conf.image_offset), cosl(0.5L*conf.image_offset));
  for (long int j = 0; j < state.number_modes/2-1; j++) {
    conf.w[j] = cexpl(-1.0IL*(j+1)*(conf.origin_offset - 2.0IL*atanhl(conf.scaling)))*overN;
  }
  for (int j = 0; j < state.number_modes; j++) {
    long double q = 2.L*PI*(1.0L*j*overN - 0.5L);
    conf.dq[j] = b*(1.0L + a*cosl(q - conf.origin_offset));
  }
}

void remap(map_ptr new_map, unsigned long int N) {
  long double s = sinl(0.5L*new_map->image_offset);
  long double c = cosl(0.5L*new_map->image_offset);
  long double beta = tanl(0.5L*(new_map->image_offset - conf.image_offset));
  long double overN0 = 1.L/state.number_modes;
  long double overN = 1.L/N;
  long double R_TOL = 1.0E-10L;
  unsigned long int N0 = state.number_modes;
  unsigned int QC_pass = 0;

  map old_map;
  memcpy(&old_map, &conf, sizeof(map));
  old_map.w = NULL;
  old_map.dq = NULL;

  new_map->origin_offset = 2.0L*atan2l(new_map->scaling*s,c);
  memcpy(tmpc[0], data[0], N0*sizeof(fftwl_complex));
  memcpy(tmpc[1], data[1], N0*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  saveQ = fftwl_malloc(N0*sizeof(fftwl_complex));
  saveV = fftwl_malloc(N0*sizeof(fftwl_complex));
  memcpy(saveQ, tmpc[0], N0*sizeof(fftwl_complex));
  memcpy(saveV, tmpc[1], N0*sizeof(fftwl_complex));
  // verify refinement
  /*for (long int j = 0; j < state.number_modes/2; j++) {
    tmpc[0][j] = tmpc[0][j]*overN0;
    tmpc[1][j] = tmpc[1][j]*overN0;
  }
  complex_array_out("preref.Q.ft.txt", tmpc[0]);
  complex_array_out("preref.V.ft.txt", tmpc[1]);
  */
  deallocate_memory();
  state.number_modes = N;
  allocate_memory(); 
  
  if (N0 < state.number_modes) {
    unsigned long int mem_offset = state.number_modes-N0/2;
    memcpy(tmpc[0], saveQ, (N0/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1], saveV, (N0/2)*sizeof(fftwl_complex));
    memcpy(tmpc[0]+mem_offset, saveQ+N0/2, (N0/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1]+mem_offset, saveV+N0/2, (N0/2)*sizeof(fftwl_complex));
  } else if (N0 > state.number_modes) {
    unsigned long int mem_offset = N0 - state.number_modes/2;
    memcpy(tmpc[0], saveQ, (state.number_modes/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1], saveV, (state.number_modes/2)*sizeof(fftwl_complex));
    memcpy(tmpc[0]+state.number_modes/2, saveQ+mem_offset, (state.number_modes/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1]+state.number_modes/2, saveV+mem_offset, (state.number_modes/2)*sizeof(fftwl_complex));
  } else {
    memcpy(tmpc[0], saveQ, state.number_modes*sizeof(fftwl_complex));
    memcpy(tmpc[1], saveV, state.number_modes*sizeof(fftwl_complex));
  }
  long double 		q, q_r, q0 = PI + conf.origin_offset;
  fftwl_complex 	w = 1.L;
  for (long int j = 0; j < state.number_modes; j++) {
    q   = new_map->scaling*tanl(1.L*PI*(j*overN - 0.5L) - 0.5L*new_map->origin_offset);
    q_r = q0 + 2.0L*atan2l(beta+q, conf.scaling*(1.0L - beta*q));
    w = cexpl(-1.IL*q_r);
    data[0][j] = tmpc[0][state.number_modes/2-1];
    data[1][j] = tmpc[1][state.number_modes/2-1];
    for (long int l = state.number_modes/2-1; l > 0; l--) {
      data[0][j] = data[0][j]*w + tmpc[0][l-1];
      data[1][j] = data[1][j]*w + tmpc[1][l-1];
    } 
    data[0][j] = data[0][j]*overN0;
    data[1][j] = data[1][j]*overN0;
  }
  // verify refinement
  //complex_array_out("new_d2.txt", data[0]); // this is wrong
  //complex_array_out("new_d1.txt", data[1]);
  conf.scaling = new_map->scaling;
  conf.image_offset = new_map->image_offset;
  set_mapping();


  // verify that new map passes QC 
  memcpy(tmpc[0], data[0], state.number_modes*sizeof(fftwl_complex));
  memcpy(tmpc[1], data[1], state.number_modes*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  for (long int j = 0; j < state.number_modes/2; j++) {
    tmpc[0][j] = tmpc[0][j]*overN;
    tmpc[1][j] = tmpc[1][j]*overN;
  }
  memset(tmpc[0]+state.number_modes/2, 0, state.number_modes/2*sizeof(fftwl_complex));
  memset(tmpc[1]+state.number_modes/2, 0, state.number_modes/2*sizeof(fftwl_complex));
  map_quality(tmpc[0], tmpc[1], R_TOL, &QC_pass);//*sqrtl(state.number_modes/4096.L)
  // verify refinement
  //complex_array_out("postref.Q.ft.txt", tmpc[0]);
  //complex_array_out("postref.V.ft.txt", tmpc[1]);
  if (QC_pass) {
    //printf("Good Quality New Map\nProceed with new map\n");
  } else {
    //printf("Bad Quality New Map\nTry with different parameters\n");
    deallocate_memory();
    state.number_modes = N0;
    memcpy(&conf, &old_map, sizeof(map));
    allocate_memory(); 
    memcpy(tmpc[0], saveQ, N0*sizeof(fftwl_complex));
    memcpy(tmpc[1], saveV, N0*sizeof(fftwl_complex));
    set_mapping();
    memset(tmpc[0] + N0/2, 0, N0/2*sizeof(fftwl_complex));
    memset(tmpc[1] + N0/2, 0, N0/2*sizeof(fftwl_complex));
    fftwl_execute(ft0);
    fftwl_execute(ft1);
    for (long int j = 0; j < N0; j++) {
      tmpc[0][j] = tmpc[0][j]*overN0;
      tmpc[1][j] = tmpc[1][j]*overN0;
    }
    memcpy(data[0], tmpc[0], N0*sizeof(fftwl_complex));
    memcpy(data[1], tmpc[1], N0*sizeof(fftwl_complex));
  }
  fftwl_free(saveQ);
  fftwl_free(saveV); 
}

void map_quality(fftwl_complex *in1, fftwl_complex *in2, long double tol, unsigned int *QC_pass) {
  long double qc_ratio	= 1.0L;
  //long double nqc_ratio	= 1.0L;
  unsigned long N = state.number_modes;
  unsigned long ks = lroundl(7.L/12.L*(N/2)); //lroundl(2.L*state.kD/3.L);
  unsigned long kf = lroundl(9.L/12.L*(N/2)); //	state.kD;
  left_sum1	= 0.0L;
  left_sum2	= 0.0L;
  full_sum1	= 0.0L;
  full_sum2	= 0.0L;
  partial_sum1	= 0.0L;
  partial_sum2	= 0.0L;
  narrow_sum1	= 0.0L; 
  narrow_sum2	= 0.0L; 
  

  for (long int j = state.number_modes/2-1; j > -1; j--) {
    full_sum1 += creall(in1[j])*creall(in1[j]) + cimagl(in1[j])*cimagl(in1[j]);
    full_sum2 += creall(in2[j])*creall(in2[j]) + cimagl(in2[j])*cimagl(in2[j]);
    if (j == state.number_modes*7/16) {
      partial_sum1 = full_sum1;
      partial_sum2 = full_sum2;
    } else if (j == state.number_modes*3/16) {
      narrow_sum1 = full_sum1;
      narrow_sum2 = full_sum2;
    } else if (j == kf) {
      left_sum1 = full_sum1;
      left_sum2 = full_sum2;
    } else if (j == ks) {
      partial_sum1 = full_sum1;
      partial_sum2 = full_sum2;
    }

  }
  qc_ratio = sqrtl((partial_sum1-left_sum1)/(1.L + full_sum1));
  qc_ratio = fmaxl(qc_ratio, sqrtl((partial_sum2-left_sum2)/sqrtl(1.L+full_sum2)) );
  //inqc_ratio = narrow_sum1/sqrtl(1.L + powl(full_sum1, 2));
  //nqc_ratio = fmaxl(nqc_ratio, narrow_sum2/sqrtl(1.L+powl(full_sum2, 2)));
  if (qc_ratio < tol) {
        //printf("QC Pass\nQC ratio is %.9LE\n", qc_ratio);
	*QC_pass = 1;
	//if (nqc_ratio < tol) *QC_pass = 2;
  } else {
        printf("QC Fail\nQC ratio is %.9LE\n", qc_ratio);
 	*QC_pass = 0;
  }
}

void map_quality_fourier(fftwl_complex *inQ, fftwl_complex *inV, long double tol, unsigned int *QC_pass){
    memcpy(tmpc[0], inQ, state.number_modes*sizeof(fftwl_complex));
    memcpy(tmpc[1], inV, state.number_modes*sizeof(fftwl_complex));
    fftwl_execute(ift0); 
    fftwl_execute(ift1); 
    map_quality(tmpc[0], tmpc[1], tol, QC_pass);
}

void track_singularity(fftwl_complex *inQ) {
  long int 	N 	= state.number_modes;
  long double   overN 	= 1.L/N; 
  long double   maxabsd2Q = 1.L;
  long double   q_max   = 0.L;

  memcpy(tmpc[0], inQ, N*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  for (long int j = 0; j < N/2; j++) {
    tmpc[0][j] = -1.L*j*j*tmpc[0][j]*overN;
  }
  memset(tmpc[0]+N/2, 0, N/2*sizeof(fftwl_complex));
  fftwl_execute(ft0);
  for (long int j = 0; j < N; j++) {
    tmpr[0][j] = cabsl(tmpc[0][j]);
    if (maxabsd2Q < tmpr[0][j]) {
	maxabsd2Q = tmpr[0][j];
        q_max = PI*(2.L*j*overN - 1.L);
    }
  } 
  if (MOVE_MESH) { 
    alt_map.origin_offset = q_max;
  } else {
    alt_map.origin_offset = conf.origin_offset;
    q_max = conf.origin_offset;
  }
  alt_map.image_offset = conf.image_offset + 2.0L*atan2l(conf.scaling*sinl(0.5L*(q_max-conf.origin_offset)), cosl(0.5L*(q_max-conf.origin_offset)));
 // printf("max_Abs_d2Q = %.19LE\tq_max = %.19LE\tu_max = %.19LE\n", maxabsd2Q, q_max, alt_map.image_offset);
}


