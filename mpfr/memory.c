#include "header.h"

#define COMPLEX_ARRAYS 	6	
#define REAL_ARRAYS	2

static const unsigned long int n_complex_arrays = COMPLEX_ARRAYS;
static const unsigned long int n_real_arrays = REAL_ARRAYS;
static mpfc_t	*aux_array;

mpfft_plan	ft0, ft1, ft2, ft3, ft4;
mpfft_plan	ift0, ift1, ift2, ift3, ift4;
mpfr_prec_t PREC;
mpfr_t 	**tmpr;
mpfc_t 	**tmpc;
mpfc_t 	**data;

// Global Variable

void init_memory() {
  PREC = state.precision; 
  mpfr_init2(Pie, state.precision);  
  mpfr_const_pi(Pie, MODE);
  unsigned long int N = 1<<state.nbits;
  mpfr_printf("Requesting %lu complex arrays of length %lu\n", n_complex_arrays, N);
  mpfr_printf("Requesting %lu real arrays of length %lu\n", n_real_arrays, N);
  tmpc = (mpfc_t **) malloc(n_complex_arrays*sizeof(mpfc_t *)); 
  tmpr = (mpfr_t **) malloc(n_real_arrays*sizeof(mpfr_t *));   
  data = (mpfc_t **) malloc(2*sizeof(mpfc_t *)); 
  // init memory() complete
}

void allocate_memory() {
  unsigned long int N = 1<<state.nbits;
  for (long int j = 0; j < n_complex_arrays; j++) {
    tmpc[j] = (mpfc_t *) malloc(N*sizeof(mpfc_t));
    //memset(tmpc[j], 0, N*sizeof(mpfc_t));
    for (long int k = 0; k < N; k++) {
      mpfr_inits2(PREC, tmpc[j][k].re, tmpc[j][k].im, (mpfr_ptr) NULL);
      mpfr_set_ui(tmpc[j][k].re, 0, MODE);
      mpfr_set_ui(tmpc[j][k].im, 0, MODE);
    }
  }
  for (long int j = 0; j < n_real_arrays; j++) {
    tmpr[j] = (mpfr_t *) malloc(N*sizeof(mpfr_t));
    //memset(tmpr[j], 0, N*sizeof(mpfr_t));
    for (long int k = 0; k < N; k++) {
      mpfr_inits2(PREC, tmpr[j][k], (mpfr_ptr) NULL);
      mpfr_set_ui(tmpr[j][k], 0, MODE);
    }
  }
  conf.dq = (mpfr_t *) malloc(N*sizeof(mpfr_t));
  data[0] = (mpfc_t *) malloc(N*sizeof(mpfc_t));
  data[1] = (mpfc_t *) malloc(N*sizeof(mpfc_t));
  aux_array = (mpfc_t *) malloc(N*sizeof(mpfc_t));
  for (long int j = 0; j < N; j++) {
    mpfr_inits2(PREC, conf.dq[j], (mpfr_ptr) NULL);
    mpfr_inits2(PREC, data[0][j].re, data[0][j].im, (mpfr_ptr) NULL);
    mpfr_inits2(PREC, data[1][j].re, data[1][j].im, (mpfr_ptr) NULL);
    mpfr_inits2(PREC, aux_array[j].re, aux_array[j].im, (mpfr_ptr) NULL);
    mpfr_set_ui(conf.dq[j], 0, MODE);
    mpfr_set_ui(data[0][j].re, 0, MODE);
    mpfr_set_ui(data[0][j].im, 0, MODE);
    mpfr_set_ui(data[1][j].re, 0, MODE);
    mpfr_set_ui(data[1][j].im, 0, MODE);
    mpfr_set_ui(aux_array[j].re, 0, MODE);
    mpfr_set_ui(aux_array[j].im, 0, MODE);
  }
  conf.w = (mpfc_t *) malloc((N/2-1)*sizeof(mpfc_t));
  for (long int j = 0; j < N/2-1; j++) {
    mpfr_inits2(PREC, conf.w[j].re, conf.w[j].im, (mpfr_ptr) NULL);
    mpfr_set_ui(conf.w[j].re, 0, MODE);
    mpfr_set_ui(conf.w[j].im, 0, MODE);
  }
  
  allocate_timemarching();

  // mpfft_plan plan_forward = mpfft_create_plan_1d(F, f, nbits, precision, FFT_FORWARD); example plan

  ft0 = mpfft_create_plan_1d(tmpc[0], tmpc[0], state.nbits, PREC, FFT_FORWARD);
  ft1 = mpfft_create_plan_1d(tmpc[1], tmpc[1], state.nbits, PREC, FFT_FORWARD);
  ft2 = mpfft_create_plan_1d(tmpc[2], tmpc[2], state.nbits, PREC, FFT_FORWARD);
  ft3 = mpfft_create_plan_1d(tmpc[3], tmpc[3], state.nbits, PREC, FFT_FORWARD);
  ft4 = mpfft_create_plan_1d(tmpc[4], tmpc[4], state.nbits, PREC, FFT_FORWARD);

  ift0 = mpfft_create_plan_1d(tmpc[0], tmpc[0], state.nbits, PREC, FFT_BACKWARD);
  ift1 = mpfft_create_plan_1d(tmpc[1], tmpc[1], state.nbits, PREC, FFT_BACKWARD);
  ift2 = mpfft_create_plan_1d(tmpc[2], tmpc[2], state.nbits, PREC, FFT_BACKWARD);
  ift3 = mpfft_create_plan_1d(tmpc[3], tmpc[3], state.nbits, PREC, FFT_BACKWARD);
  ift4 = mpfft_create_plan_1d(tmpc[4], tmpc[4], state.nbits, PREC, FFT_BACKWARD);
  //memset(data[0], 0, N*sizeof(mpfc_t));
  //memset(data[1], 0, N*sizeof(mpfc_t));
}

/*
void deallocate_memory() {
  for (long int j = 0; j < n_complex_arrays; j++) free(tmpc[j]);
  for (long int j = 0; j < n_real_arrays; j++) free(tmpr[j]);	 
  free(conf.dq);   free(conf.w); 
  free(aux_array); free(data[0]);  free(data[1]);
  deallocate_timemarching();
  mpfft_destroy_plan(ft0);
  mpfft_destroy_plan(ft1);
  mpfft_destroy_plan(ft2);
  mpfft_destroy_plan(ft3);
  mpfft_destroy_plan(ft4);
  mpfft_destroy_plan(ift0);
  mpfft_destroy_plan(ift1);
  mpfft_destroy_plan(ift2);
  mpfft_destroy_plan(ift3);
  mpfft_destroy_plan(ift4);
}
*/
/*
void fft_shift(mpfc_t *in){
  for (long int j = 0; j < state.number_modes/2; j++) {
    aux_array[j] 		= in[state.number_modes/2+j];
    aux_array[j+state.number_modes/2]	= in[j]; 
  }
  memcpy(in, aux_array, state.number_modes*sizeof(mpfc_t));
}
*/


