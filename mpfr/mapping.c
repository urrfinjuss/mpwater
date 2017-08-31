#include "header.h"

static mpfc_t *saveQ, *saveV;
static mpfr_t full_sum1;
static mpfr_t full_sum2;
static mpfr_t narrow_sum1; 
static mpfr_t narrow_sum2; 
static mpfr_t partial_sum1;
static mpfr_t partial_sum2;
static mpfr_t left_sum1;
static mpfr_t left_sum2;

void set_mapping() {
  long int N = 1<<state.nbits;
  mpfr_t overN, eye, a, b;
  
  mpfr_inits2(state.precision, eye, overN, a, b, (mpfr_ptr) NULL );
  mpfr_set_ui(eye, 1, MODE);
  mpfr_div_ui(overN, eye, 1<<state.nbits, MODE);
  
  mpfr_t buf1, buf2, buf3, buf4;
  mpfr_inits2(state.precision, buf1, buf2, buf3, buf4, (mpfr_ptr) NULL); 
  mpfr_mul(buf1, conf.scaling, conf.scaling, MODE);
  mpfr_ui_sub(buf2, 1, buf1, MODE);
  mpfr_add_ui(buf3, buf1, 1, MODE);
  mpfr_div(a, buf2, buf3, MODE);  // (1-L^2)/(1+L^2)

  mpfr_div(b, buf3, conf.scaling, MODE);
  mpfr_div_ui(b, b, 2, MODE); // (1+L^2)/(2 L)

  //mpfr_t a = (1.0L - powl(conf.scaling, 2.L))/(1.0L + powl(conf.scaling, 2.L));
  //mpfr_t b = 0.5L*(1.0L + powl(conf.scaling, 2.L))/conf.scaling;

  mpfr_div_ui(buf1, conf.image_offset, 2, MODE);
  mpfr_sin_cos(buf2, buf3, buf1, MODE);
  mpfr_mul(buf2, buf2, conf.scaling, MODE);

  mpfr_atan2(buf4, buf2, buf3, MODE);

  mpfr_init2(conf.origin_offset, state.precision);
  mpfr_mul_ui(conf.origin_offset, buf4, 2, MODE); 

  //conf.origin_offset = 2.0L*atan2l(conf.scaling*sinl(0.5L*conf.image_offset), cosl(0.5L*conf.image_offset));
  for (long int j = 0; j < N/2-1; j++) {
    mpfr_atanh(buf1, conf.scaling, MODE);
    mpfr_mul_si(buf1, buf1, -2*(j+1), MODE);
    mpfr_exp(buf1, buf1, MODE);
    mpfr_mul(buf1, buf1, overN, MODE);   

    mpfr_mul_si(buf2, conf.origin_offset, -(j+1), MODE);
    mpfr_sin_cos(buf4, buf3, buf2, MODE);
    mpfr_mul(conf.w[j].re, buf1, buf3, MODE);
    mpfr_mul(conf.w[j].im, buf1, buf4, MODE);
  
    //conf.w[j] = cexpl(-1.0IL*(j+1)*(conf.origin_offset - 2.0IL*atanhl(conf.scaling)))*overN;
  }
  mpfr_t q;
  mpfr_init2(q, state.precision);
  for (int j = 0; j < N; j++) {
    //mpfr_t q = 2.L*PI*(1.0L*j*overN - 0.5L);
    mpfr_mul_si(q, overN, 2*j, MODE);
    mpfr_sub_ui(q, q, 1, MODE);
    mpfr_mul(q, q, Pie, MODE);

    mpfr_sub(buf1, q, conf.origin_offset, MODE);
    mpfr_cos(buf1, buf1, MODE);
    mpfr_mul(buf1, buf1, a, MODE);
    mpfr_add_ui(buf1, buf1, 1, MODE);
    mpfr_mul( conf.dq[j], b, buf1, MODE);
    //conf.dq[j] = b*(1.0L + a*cosl(q - conf.origin_offset));
  }
  mpfr_clears(overN, buf1, buf2, buf3, buf4, q, (mpfr_ptr) NULL); 
}
/*
void remap(map_ptr new_map, unsigned long int N) {
  mpfr_t s = sinl(0.5L*new_map->image_offset);
  mpfr_t c = cosl(0.5L*new_map->image_offset);
  mpfr_t beta = tanl(0.5L*(new_map->image_offset - conf.image_offset));
  mpfr_t overN0 = 1.L/state.number_modes;
  mpfr_t overN = 1.L/N;
  mpfr_t R_TOL = 4.0E-16L;
  unsigned long int N0 = state.number_modes;
  unsigned int QC_pass = 0;

  map old_map;
  memcpy(&old_map, &conf, sizeof(map));
  old_map.w = NULL;
  old_map.dq = NULL;

  new_map->origin_offset = 2.0L*atan2l(new_map->scaling*s,c);
  memcpy(tmpc[0], data[0], N0*sizeof(mpfc_t));
  memcpy(tmpc[1], data[1], N0*sizeof(mpfc_t));
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  saveQ = fftwl_malloc(N0*sizeof(mpfc_t));
  saveV = fftwl_malloc(N0*sizeof(mpfc_t));
  memcpy(saveQ, tmpc[0], N0*sizeof(mpfc_t));
  memcpy(saveV, tmpc[1], N0*sizeof(mpfc_t));
  //for (long int j = 0; j < state.number_modes/2; j++) {
  //  tmpc[0][j] = tmpc[0][j]*overN0;
  //  tmpc[1][j] = tmpc[1][j]*overN0;
  //}
  //complex_array_out("preref.Q.ft.txt", tmpc[0]);
  //complex_array_out("preref.V.ft.txt", tmpc[1]);
  deallocate_memory();
  state.number_modes = N;
  allocate_memory(); 
  
  if (N0 < state.number_modes) {
    unsigned long int mem_offset = state.number_modes-N0/2;
    memcpy(tmpc[0], saveQ, (N0/2)*sizeof(mpfc_t));
    memcpy(tmpc[1], saveV, (N0/2)*sizeof(mpfc_t));
    memcpy(tmpc[0]+mem_offset, saveQ+N0/2, (N0/2)*sizeof(mpfc_t));
    memcpy(tmpc[1]+mem_offset, saveV+N0/2, (N0/2)*sizeof(mpfc_t));
  } else if (N0 > state.number_modes) {
    unsigned long int mem_offset = N0 - state.number_modes/2;
    memcpy(tmpc[0], saveQ, (state.number_modes/2)*sizeof(mpfc_t));
    memcpy(tmpc[1], saveV, (state.number_modes/2)*sizeof(mpfc_t));
    memcpy(tmpc[0]+state.number_modes/2, saveQ+mem_offset, (state.number_modes/2)*sizeof(mpfc_t));
    memcpy(tmpc[1]+state.number_modes/2, saveV+mem_offset, (state.number_modes/2)*sizeof(mpfc_t));
  } else {
    memcpy(tmpc[0], saveQ, state.number_modes*sizeof(mpfc_t));
    memcpy(tmpc[1], saveV, state.number_modes*sizeof(mpfc_t));
  }
  mpfr_t 		q, q_r, q0 = PI + conf.origin_offset;
  mpfc_t 	w = 1.L;
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
  conf.scaling = new_map->scaling;
  conf.image_offset = new_map->image_offset;
  set_mapping();


  // verify that new map passes QC 
  memcpy(tmpc[0], data[0], state.number_modes*sizeof(mpfc_t));
  memcpy(tmpc[1], data[1], state.number_modes*sizeof(mpfc_t));
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  for (long int j = 0; j < state.number_modes/2; j++) {
    tmpc[0][j] = tmpc[0][j]*overN;
    tmpc[1][j] = tmpc[1][j]*overN;
  }
  memset(tmpc[0]+state.number_modes/2, 0, state.number_modes/2*sizeof(mpfc_t));
  memset(tmpc[1]+state.number_modes/2, 0, state.number_modes/2*sizeof(mpfc_t));
  map_quality(tmpc[0], tmpc[1], R_TOL, &QC_pass);// *sqrtl(state.number_modes/4096.L)
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
    memcpy(tmpc[0], saveQ, N0*sizeof(mpfc_t));
    memcpy(tmpc[1], saveV, N0*sizeof(mpfc_t));
    set_mapping();
    memset(tmpc[0] + N0/2, 0, N0/2*sizeof(mpfc_t));
    memset(tmpc[1] + N0/2, 0, N0/2*sizeof(mpfc_t));
    fftwl_execute(ft0);
    fftwl_execute(ft1);
    for (long int j = 0; j < N0; j++) {
      tmpc[0][j] = tmpc[0][j]*overN0;
      tmpc[1][j] = tmpc[1][j]*overN0;
    }
    memcpy(data[0], tmpc[0], N0*sizeof(mpfc_t));
    memcpy(data[1], tmpc[1], N0*sizeof(mpfc_t));
  }
  fftwl_free(saveQ);
  fftwl_free(saveV); 
}
*/


void map_quality(mpfc_t *in1, mpfc_t *in2, mpfr_t tol, unsigned int *QC_pass) {
  unsigned long N = 1<<state.nbits;
  unsigned long ks = lroundl(5.L/12.L*(N/2)); //lroundl(2.L*state.kD/3.L);
  unsigned long kf = lroundl(3.L/4.L*(N/2)); //	state.kD;
  mpfr_t tmp1, tmp2, qc_ratio1, qc_ratio2, qc;
  mpfr_inits2(state.precision, qc, qc_ratio1, qc_ratio2, tmp1, tmp2, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, left_sum1, left_sum2, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, full_sum1, full_sum2, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, partial_sum1, partial_sum2, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, narrow_sum1, narrow_sum2, (mpfr_ptr) NULL);
  mpfr_set_ui(qc_ratio1,  1, MODE);
  mpfr_set_ui(qc_ratio2,  1, MODE);
  mpfr_set_ui(left_sum1, 0, MODE);
  mpfr_set_ui(left_sum2, 0, MODE);
  mpfr_set_ui(full_sum1, 0, MODE);
  mpfr_set_ui(full_sum2, 0, MODE);
  mpfr_set_ui(partial_sum1, 0, MODE);
  mpfr_set_ui(partial_sum2, 0, MODE);
  mpfr_set_ui(narrow_sum1, 0, MODE);
  mpfr_set_ui(narrow_sum2, 0, MODE);
  
  /*
  left_sum1	= 0.0L;
  left_sum2	= 0.0L;
  full_sum1	= 0.0L;
  full_sum2	= 0.0L;
  partial_sum1	= 0.0L;
  partial_sum2	= 0.0L;
  narrow_sum1	= 0.0L; 
  narrow_sum2	= 0.0L; 
  */

  for (long int j = N/2-1; j > -1; j--) {
    /*
    full_sum1 += cabsl(in1[j]);
    full_sum2 += cabsl(in2[j]);
    */
    //mpfr_hypot(tmp1, in1[j].re, in1[j].im, MODE);
    mpfr_mul(tmp1, in1[j].re, in1[j].re, MODE);
    mpfr_fma(tmp1, in1[j].im, in1[j].im, tmp1, MODE);
    //mpfr_hypot(tmp2, in2[j].re, in2[j].im, MODE);
    mpfr_mul(tmp2, in2[j].re, in2[j].re, MODE);
    mpfr_fma(tmp2, in2[j].im, in2[j].im, tmp2, MODE);
    
    mpfr_add(full_sum1, full_sum1, tmp1, MODE);
    mpfr_add(full_sum2, full_sum2, tmp2, MODE);
    if (j == N*7/16) {
      /*
      partial_sum1 = full_sum1;
      partial_sum2 = full_sum2;
      */
      mpfr_set(partial_sum1, full_sum1, MODE);
      mpfr_set(partial_sum2, full_sum2, MODE);
    } else if (j == N*3/16) {
      /*
      narrow_sum1 = full_sum1;
      narrow_sum2 = full_sum2;
      */
      mpfr_set(narrow_sum1, full_sum1, MODE);
      mpfr_set(narrow_sum2, full_sum2, MODE);
    } else if (j == kf) {
      /*
      left_sum1 = full_sum1;
      left_sum2 = full_sum2;
      */
      mpfr_set(left_sum1, full_sum1, MODE);
      mpfr_set(left_sum2, full_sum2, MODE);
    } else if (j == ks) {
      /*
      partial_sum1 = full_sum1;
      partial_sum2 = full_sum2;
      */
      mpfr_set(partial_sum1, full_sum1, MODE);
      mpfr_set(partial_sum2, full_sum2, MODE);
    }
  }
  
  mpfr_sqrt(full_sum1, full_sum1, MODE);
  mpfr_sqrt(full_sum2, full_sum2, MODE);
  mpfr_sqrt(left_sum1, left_sum1, MODE);
  mpfr_sqrt(left_sum2, left_sum2, MODE);
  mpfr_sqrt(partial_sum1, partial_sum1, MODE);
  mpfr_sqrt(partial_sum2, partial_sum2, MODE);
  mpfr_sqrt(narrow_sum1, narrow_sum1, MODE);
  mpfr_sqrt(narrow_sum2, narrow_sum2, MODE);
  
 
  // qc ratio for Q 
  mpfr_mul(tmp1, left_sum1, left_sum1, MODE); 
  mpfr_fms(qc_ratio1, partial_sum1, partial_sum1, tmp1, MODE);
  mpfr_mul(tmp1, full_sum1, full_sum1, MODE);
  mpfr_add_ui(tmp1, tmp1, 1, MODE);
  mpfr_div(qc_ratio1, qc_ratio1, tmp1, MODE);
  mpfr_sqrt(qc_ratio1, qc_ratio1, MODE);
  // qc ratio for V
  mpfr_mul(tmp2, left_sum2, left_sum2, MODE);
  mpfr_fms(qc_ratio2, partial_sum2, partial_sum2, tmp2, MODE);
  mpfr_mul(tmp2, full_sum2, full_sum2, MODE);
  mpfr_add_ui(tmp2, tmp2, 1, MODE);
  mpfr_div(qc_ratio2, qc_ratio2, tmp2, MODE);
  mpfr_sqrt(qc_ratio2, qc_ratio2, MODE);
  if (mpfr_cmp(qc_ratio1, qc_ratio2) > 0) mpfr_sqr(qc, qc_ratio2, MODE);
  else mpfr_sqr(qc, qc_ratio1, MODE);
  //mpfr_sqr(qc, qc, MODE);  

  //qc_ratio = (partial_sum1-left_sum1)/sqrtl(1.L + powl(full_sum1, 2));
  //qc_ratio = fmaxl(qc_ratio, (partial_sum2-left_sum2)/sqrtl(1.L+powl(full_sum2, 2)));
  //inqc_ratio = narrow_sum1/sqrtl(1.L + powl(full_sum1, 2));
  //nqc_ratio = fmaxl(nqc_ratio, narrow_sum2/sqrtl(1.L+powl(full_sum2, 2)));
  //if (qc_ratio < tol) {
  if (mpfr_cmp(qc, tol) < 0) {
        mpfr_printf("QC Pass\nQC ratio is %.9Re\n", qc);
	*QC_pass = 1;
	//if (nqc_ratio < tol) *QC_pass = 2;
  } else {
        mpfr_printf("QC Fail\nQC ratio is %.9Re\n", qc);
 	*QC_pass = 0;
  }
  mpfr_clears(qc, qc_ratio1, qc_ratio2, tmp1, tmp2, (mpfr_ptr) NULL);
  mpfr_clears(left_sum1, left_sum2, (mpfr_ptr) NULL);
  mpfr_clears(full_sum1, full_sum2, (mpfr_ptr) NULL);
  mpfr_clears(partial_sum1, partial_sum2, (mpfr_ptr) NULL);
  mpfr_clears(narrow_sum1, narrow_sum2, (mpfr_ptr) NULL);
}



void map_quality_fourier(mpfc_t *inQ, mpfc_t *inV, mpfr_t tol, unsigned int *QC_pass){
    long int N = 1<<state.nbits;
    for (long int j = 0; j < N; j++) {
      mpfr_set(tmpc[0][j].re, inQ[j].re, MODE);
      mpfr_set(tmpc[0][j].im, inQ[j].im, MODE);
      mpfr_set(tmpc[1][j].re, inV[j].re, MODE);
      mpfr_set(tmpc[1][j].im, inV[j].im, MODE);
    }
    //memcpy(tmpc[0], inQ, state.number_modes*sizeof(mpfc_t));
    //memcpy(tmpc[1], inV, state.number_modes*sizeof(mpfc_t));
    mpfft_execute(ift0); 
    mpfft_execute(ift1); 
    map_quality(tmpc[0], tmpc[1], tol, QC_pass);
}


/*
void track_singularity(mpfc_t *inQ) {
  long int 	N 	= state.number_modes;
  mpfr_t   overN 	= 1.L/N; 
  mpfr_t   maxabsd2Q = 1.L;
  mpfr_t   q_max   = 0.L;

  memcpy(tmpc[0], inQ, N*sizeof(mpfc_t));
  fftwl_execute(ift0);
  for (long int j = 0; j < N/2; j++) {
    tmpc[0][j] = -1.L*j*j*tmpc[0][j]*overN;
  }
  memset(tmpc[0]+N/2, 0, N/2*sizeof(mpfc_t));
  fftwl_execute(ft0);
  for (long int j = 0; j < N; j++) {
    tmpr[0][j] = cabsl(tmpc[0][j]);
    if (maxabsd2Q < tmpr[0][j]) {
	maxabsd2Q = tmpr[0][j];
        q_max = PI*(2.L*j*overN - 1.L);
    }
  }  
  alt_map.origin_offset = q_max;
  alt_map.image_offset = conf.image_offset + 2.0L*atan2l(conf.scaling*sinl(0.5L*(q_max-conf.origin_offset)), cosl(0.5L*(q_max-conf.origin_offset)));
 // printf("max_Abs_d2Q = %.19LE\tq_max = %.19LE\tu_max = %.19LE\n", maxabsd2Q, q_max, alt_map.image_offset);
}
*/

