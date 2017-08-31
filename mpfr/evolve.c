#include "header.h"


static mpfc_t 	**kq, *tQ;
static mpfc_t 	**kv, *tV;

static mpfr_t	cfl; // used to be 0.080L
static unsigned long    kZ;
static const unsigned long 	pD = 12;
static mpfr_t 	one_third;        
static mpfr_t 	two_thirds;       
static mpfr_t 	one_twelfth;      
static mpfr_t 	one_eighth;     
static mpfr_t 	one_sixteenth; 	 
static mpfr_t 	one_fortyfourths; 
static mpfr_t 	one_onetwentieth; 

void init_timemarching() {
  mpfr_init2(state.time, state.precision);  // extra init
  mpfr_set_ui(state.time, 0, MODE);
  mpfr_t eye;
  mpfr_init2(eye, state.precision);
  mpfr_inits2(state.precision, one_third, two_thirds, one_twelfth, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, one_eighth, one_sixteenth, one_fortyfourths, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, one_onetwentieth, (mpfr_ptr) NULL);
  mpfr_set_ui(eye, 1, MODE);
  mpfr_div_ui(one_third, eye, 3, MODE);
  mpfr_mul_ui(two_thirds, one_third, 2, MODE);
  mpfr_div_ui(one_twelfth, eye, 12, MODE);
  mpfr_div_ui(one_eighth, eye, 8, MODE);
  mpfr_div_ui(one_sixteenth, eye, 16, MODE);
  mpfr_div_ui(one_fortyfourths, eye, 44, MODE);
  mpfr_div_ui(one_onetwentieth, eye, 120, MODE);
   
  kq = malloc(7*sizeof(mpfc_t *));
  kv = malloc(7*sizeof(mpfc_t *));
  /*
  for (int j = 0; j < 7; j++) {
    mpfr_inits2(state.precision, kq[j].re, kq[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(state.precision, kv[j].re, kv[j].im, (mpfr_ptr) NULL);
  }*/
  mpfr_t buf;
  mpfr_inits2(state.precision, buf, (mpfr_ptr) NULL);
  mpfr_mul_ui(buf, eye, 1<<state.nbits, MODE);
  mpfr_mul_ui(buf, buf, 11, MODE);
  mpfr_div_ui(buf, buf, 24, MODE);
  state.kD = mpfr_get_ui(buf, MODE);

  mpfr_mul_ui(buf, buf, 7, MODE);
  mpfr_div_ui(buf, buf, 6, MODE);
  kZ = mpfr_get_ui(buf, MODE);
  if (kZ > 1<<(state.nbits-1)) kZ = 1<<(state.nbits-1);
  printf("Dissipative Range starts at kD = %lu\n", state.kD);
  printf("Zero Range starts at kZ = %lu\n", kZ);

  mpfr_init2(cfl, state.precision);
  mpfr_set_d(cfl, 0.24, MODE);
}

void allocate_timemarching() {
  unsigned long N = (1<<state.nbits);
  mpfr_t eye;
  mpfr_init2(eye, state.precision);
  mpfr_set_ui(eye, 1, MODE);
  for (long int j = 0; j < 7; j++) {
    kq[j] = malloc(N*sizeof(mpfc_t));
    kv[j] = malloc(N*sizeof(mpfc_t));
    for (int k = 0; k < N; k++) {
      mpfr_inits2(state.precision, kq[j][k].re, kq[j][k].im, (mpfr_ptr) NULL);
      mpfr_inits2(state.precision, kv[j][k].re, kv[j][k].im, (mpfr_ptr) NULL);
      mpfr_set_ui(kq[j][k].re, 0, MODE);
      mpfr_set_ui(kq[j][k].im, 0, MODE);
      mpfr_set_ui(kv[j][k].re, 0, MODE);
      mpfr_set_ui(kv[j][k].im, 0, MODE);
    }
  }
  mpfr_t buf;
  mpfr_inits2(state.precision, buf, (mpfr_ptr) NULL);
  mpfr_mul_ui(buf, eye, 1<<state.nbits, MODE);
  mpfr_mul_ui(buf, buf, 11, MODE);
  mpfr_div_ui(buf, buf, 24, MODE);
  state.kD = mpfr_get_ui(buf, MODE);

  mpfr_mul_ui(buf, buf, 7, MODE);
  mpfr_div_ui(buf, buf, 6, MODE);
  kZ = mpfr_get_ui(buf, MODE);
  if (kZ > 1<<(state.nbits-1)) kZ = 1<<(state.nbits-1);
  printf("Dissipative Range starts at kD = %lu\n", state.kD);
  printf("Zero Range starts at kZ = %lu\n", kZ);
  tQ = malloc(N*sizeof(mpfc_t));
  tV = malloc(N*sizeof(mpfc_t));
  for (int k = 0; k < N; k++) {
    mpfr_inits2(state.precision, tQ[k].re, tQ[k].im, (mpfr_ptr) NULL);
    mpfr_inits2(state.precision, tV[k].re, tV[k].im, (mpfr_ptr) NULL);
    mpfr_set_ui(tQ[k].re, 0, MODE);
    mpfr_set_ui(tQ[k].im, 0, MODE);
    mpfr_set_ui(tV[k].re, 0, MODE);
    mpfr_set_ui(tV[k].im, 0, MODE);
  }
}

/*
void deallocate_timemarching() {
  for (long int j = 0; j < 7; j++) {
    fftwl_free(kq[j]);
    fftwl_free(kv[j]);
  }
  fftwl_free(tQ);
  fftwl_free(tV);
}
*/

void compute_rhs(mpfc_t *inQ, mpfc_t *inV, mpfc_t *outQ, mpfc_t *outV) {
  long int	 	N = 1<<state.nbits;
  mpfr_t 		overN, sigma, g, buf;
  mpfc_t 		w1, w2, bufc;
  mpfc_t 		b1U, b2U;
  
  mpfr_inits2(state.precision, overN, sigma, g, buf, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, w1.re, w1.im, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, w2.re, w2.im, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, b1U.re, b1U.im, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, b2U.re, b2U.im, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, bufc.re, bufc.im, (mpfr_ptr) NULL);

  mpfr_set_ui(overN,     1, MODE);
  mpfr_div_si(overN, overN, N, MODE);
  mpfr_set   (sigma, state.surface_tension, MODE);
  mpfr_set   (    g, state.gravity, MODE);
  
  mpfr_atanh  (   buf, conf.scaling, MODE);
  mpfr_mul_si (   buf,          buf, -2, MODE);
  mpfr_exp    (   buf,          buf, MODE);
  mpfr_sin_cos( w1.im,        w1.re, conf.origin_offset, MODE);
  mpfr_mul    ( w1.re,        w1.re, buf, MODE);  
  mpfr_mul    ( w1.im,        w1.im, buf, MODE);  

  mpfr_set    ( w2.re,	      w1.re, MODE);
  mpfr_neg    ( w2.im,        w1.im, MODE);

  mpfr_set_ui (b1U.re,  0, MODE);
  mpfr_set_ui (b1U.im,  0, MODE);
  mpfr_set_ui (b2U.re,  0, MODE);
  mpfr_set_ui (b2U.im,  0, MODE);
  
  if (mpfr_cmp_ld(conf.scaling, 1.0L) == 0) {
    mpfr_set_ui(w1.re, 0, MODE);
    mpfr_set_ui(w1.im, 0, MODE);
    mpfr_set_ui(w2.re, 0, MODE);
    mpfr_set_ui(w2.im, 0, MODE);
  }
  //memcpy(tmpc[0], inQ, N*sizeof(mpfc_t));
  //memcpy(tmpc[1], inV, N*sizeof(mpfc_t));
  for (long int j = 0; j < N; j++) {
    mpfr_set(tmpc[0][j].re, inQ[j].re, MODE);
    mpfr_set(tmpc[0][j].im, inQ[j].im, MODE);
    mpfr_set(tmpc[1][j].re, inV[j].re, MODE);
    mpfr_set(tmpc[1][j].im, inV[j].im, MODE);
  }
  //fftwl_execute(ift0); 
  //fftwl_execute(ift1); 
  mpfft_execute(ift0);
  mpfft_execute(ift1);
  //memset(tmpc[0]+N/2, 0, N/2*sizeof(mpfc_t));
  //memset(tmpc[1]+N/2, 0, N/2*sizeof(mpfc_t));
  for (long int j = 0; j < N/2; j++) {
    mpfr_set_ui(tmpc[0][N/2+j].re, 0, MODE);
    mpfr_set_ui(tmpc[0][N/2+j].im, 0, MODE);
    mpfr_set_ui(tmpc[1][N/2+j].re, 0, MODE);
    mpfr_set_ui(tmpc[1][N/2+j].im, 0, MODE);
    //tmpc[0][j] = -1.IL*j*tmpc[0][j]*overN;
    //tmpc[1][j] = -1.IL*j*tmpc[1][j]*overN;
    mpfr_set(buf, tmpc[0][j].re, MODE);
    mpfr_mul_si(tmpc[0][j].re, tmpc[0][j].im,  j, MODE);
    mpfr_mul_si(tmpc[0][j].im,           buf, -j, MODE);
    mpfr_mul(tmpc[0][j].re, tmpc[0][j].re, overN, MODE);
    mpfr_mul(tmpc[0][j].im, tmpc[0][j].im, overN, MODE);

    mpfr_set(buf, tmpc[1][j].re, MODE);
    mpfr_mul_si(tmpc[1][j].re, tmpc[1][j].im,  j, MODE);
    mpfr_mul_si(tmpc[1][j].im,           buf, -j, MODE);
    mpfr_mul(tmpc[1][j].re, tmpc[1][j].re, overN, MODE);
    mpfr_mul(tmpc[1][j].im, tmpc[1][j].im, overN, MODE);
  }
  //fftwl_execute(ft0);
  //fftwl_execute(ft1);
  mpfft_execute(ft0);
  mpfft_execute(ft1);
  for (long int j = 0; j < N; j++) {
    //tmpc[2][j]= 2.L*creall(inV[j]*conjl(inQ[j]*inQ[j]))*overN;
    mpfr_mul(buf, inQ[j].im, inQ[j].im, MODE);
    mpfr_fms(buf, inQ[j].re, inQ[j].re, buf, MODE);
    mpfr_mul(buf, inV[j].re, buf, MODE);
    mpfr_mul_ui(tmpc[2][j].re, buf, 2, MODE);
    mpfr_mul(buf, inQ[j].re, inQ[j].im, MODE);
    mpfr_mul(buf, inV[j].im, buf, MODE);
    mpfr_mul_ui(buf, buf, 4, MODE);
    mpfr_add(tmpc[2][j].re, tmpc[2][j].re, buf, MODE);
    mpfr_set_ui(tmpc[2][j].im, 0, MODE);
    //tmpc[3][j]= inV[j]*conjl(inV[j])+4.L*sigma*conf.dq[j]*cimagl(tmpc[0][j]*conjl(inQ[j]));
    mpfr_mul(buf, inV[j].im, inV[j].im, MODE);
    mpfr_fma(tmpc[3][j].re, inV[j].re, inV[j].re, buf, MODE);
    
    mpfr_mul(buf, tmpc[0][j].re, inQ[j].im, MODE);
    mpfr_fms(buf, tmpc[0][j].im, inQ[j].re, buf, MODE);
    mpfr_mul(buf, buf, conf.dq[j], MODE);
    mpfr_mul_ui(buf, buf, 4, MODE);
    mpfr_fma   (tmpc[3][j].re, sigma, buf, tmpc[3][j].re);
    mpfr_set_ui(tmpc[3][j].im, 0, MODE);
    //tmpc[3][j]= tmpc[3][j]*overN;
    mpfr_mul(tmpc[3][j].re, tmpc[3][j].re, overN, MODE);
  }
  /*
  project(tmpc[3], tmpc[4]); 
  complex_array_out("B.ph.old.txt", tmpc[4]);
  for (long int j = 0; j < N; j++) {
    tmpc[3][j]= tmpc[3][j]*overN;
  }
  complex_array_out("tmpc2.txt", tmpc[2]); 
  */

  //fftwl_execute(ift2);
  //fftwl_execute(ift3);
  mpfft_execute(ift2);
  mpfft_execute(ift3);
  //tmpc[4][0] = 0.L;
  mpfr_set_ui(tmpc[4][0].re, 0, MODE);
  mpfr_set_ui(tmpc[4][0].im, 0, MODE);
  /*b2B = tmpc[3][N/2-1];	b1B = tmpc[3][N/2+1];*/

  //b2U = tmpc[2][N/2-1];	b1U = tmpc[2][N/2+1];
  mpfr_set(b2U.re, tmpc[2][N/2-1].re, MODE);
  mpfr_set(b2U.im, tmpc[2][N/2-1].im, MODE);
  mpfr_set(b1U.re, tmpc[2][N/2+1].re, MODE);
  mpfr_set(b1U.im, tmpc[2][N/2+1].im, MODE);
  for (long int j = 1; j < N/2 - 1; j++) {
    /*
    b1B = b1B*w1 + tmpc[3][N/2+1+j];
    b2B = b2B*w2 + tmpc[3][N/2-1-j];
    */
    
    //b1U = b1U*w1 + tmpc[2][N/2+1+j];
    mpfr_fms(bufc.re, b1U.im, w1.im, tmpc[2][N/2+1+j].re, MODE);
    mpfr_fms(bufc.re, b1U.re, w1.re, bufc.re, MODE);

    mpfr_fma(bufc.im, b1U.re, w1.im, tmpc[2][N/2+1+j].im, MODE);
    mpfr_fma(bufc.im, b1U.im, w1.re, bufc.im, MODE);

    mpfr_set(b1U.re, bufc.re, MODE);
    mpfr_set(b1U.im, bufc.im, MODE);
    //b2U = b2U*w2 + tmpc[2][N/2-1-j];
    mpfr_fms(bufc.re, b2U.im, w2.im, tmpc[2][N/2-1-j].re, MODE);
    mpfr_fms(bufc.re, b2U.re, w2.re, bufc.re, MODE);

    mpfr_fma(bufc.im, b2U.re, w2.im, tmpc[2][N/2-1-j].im, MODE);
    mpfr_fma(bufc.im, b2U.im, w2.re, bufc.im, MODE);

    mpfr_set(b2U.re, bufc.re, MODE);
    mpfr_set(b2U.im, bufc.im, MODE);
  }
  /*
  printf("b1U = %23.18LE\t%23.18LE\n", creall(b1U), cimagl(b1U));
  printf("b2U = %23.18LE\t%23.18LE\n", creall(b2U), cimagl(b2U));
  exit(1);
  */
  //memset(tmpc[2]+N/2, 0, N/2*sizeof(mpfc_t));
  //memset(tmpc[3]+N/2, 0, N/2*sizeof(mpfc_t));
  //memset(tmpc[4]+N/2, 0, N/2*sizeof(mpfc_t));
  for (long int j = 0; j < N/2; j++) {
    mpfr_set_ui(tmpc[2][N/2+j].re, 0, MODE);
    mpfr_set_ui(tmpc[2][N/2+j].im, 0, MODE);
    mpfr_set_ui(tmpc[3][N/2+j].re, 0, MODE);
    mpfr_set_ui(tmpc[3][N/2+j].im, 0, MODE);
    mpfr_set_ui(tmpc[4][N/2+j].re, 0, MODE);
    mpfr_set_ui(tmpc[4][N/2+j].im, 0, MODE);
  } 
  //tmpc[2][0] = 0.5L*tmpc[2][0];
  mpfr_div_ui(tmpc[2][0].re, tmpc[2][0].re, 2, MODE);
  mpfr_div_ui(tmpc[2][0].im, tmpc[2][0].im, 2, MODE);
  //tmpc[3][0] = 0.5L*tmpc[3][0];  
  for (long int j = 0; j < N/2; j++) {
    //tmpc[3][j] = -1.IL*j*tmpc[3][j];
    mpfr_set   (buf, tmpc[3][j].re, MODE);
    mpfr_mul_si(tmpc[3][j].re, tmpc[3][j].im, j, MODE);
    mpfr_mul_si(tmpc[3][j].im, buf, -j, MODE);
    //tmpc[4][j] = -1.IL*j*tmpc[2][j];
    mpfr_mul_si(tmpc[4][j].re, tmpc[2][j].im,  j, MODE);
    mpfr_mul_si(tmpc[4][j].im, tmpc[2][j].re, -j, MODE);
  }
  //fftwl_execute(ft2);
  //fftwl_execute(ft3);
  //fftwl_execute(ft4);
  mpfft_execute(ft2);
  mpfft_execute(ft3);
  mpfft_execute(ft4);
  for (long int j = 0; j < N; j++) {
    //tmpc[2][j] += 0.5L*(b1U*w1 - b2U*w2);
    mpfr_mul(bufc.re, b2U.im, w2.im, MODE);
    mpfr_fms(bufc.re, b2U.re, w2.re, bufc.re, MODE);
    mpfr_fma(bufc.re, b1U.im, w1.im, bufc.re, MODE);
    mpfr_fms(bufc.re, b1U.re, w1.re, bufc.re, MODE);
    mpfr_div_ui(bufc.re, bufc.re, 2, MODE);

    mpfr_mul(bufc.im, b2U.re, w2.im, MODE);
    mpfr_fms(bufc.im, b1U.re, w1.im, bufc.im, MODE);
    mpfr_fms(bufc.im, b2U.im, w2.re, bufc.im, MODE);
    mpfr_fms(bufc.im, b1U.im, w1.re, bufc.im, MODE);
    mpfr_div_ui(bufc.im, bufc.im, 2, MODE);
    
    mpfr_add(tmpc[2][j].re, tmpc[2][j].re, bufc.re, MODE);
    mpfr_add(tmpc[2][j].im, tmpc[2][j].im, bufc.im, MODE);
    /*tmpc[3][j] += 0.5L*(b1B*w1 - b2B*w2); 	// thats irrelevant for B' */
  }
  // Summary:
  // tmpc[0] <--  stores Q'
  // tmpc[1] <--  stores V'
  // tmpc[2] <--  stores U 
  // tmpc[3] <--  stores B'
  // tmpc[4] <--  stores U'
  for (long int j = 0; j < N; j++) {
    //outQ[j] = 0.5IL*conf.dq[j]*(2.L*tmpc[0][j]*tmpc[2][j]-tmpc[4][j]*inQ[j]);
    mpfr_div_ui(buf, conf.dq[j], 2, MODE);

    mpfr_mul(bufc.re, tmpc[4][j].im, inQ[j].im, MODE);
    mpfr_fms(bufc.re, tmpc[4][j].re, inQ[j].re, bufc.re, MODE);
    mpfr_mul(bufc.im, tmpc[4][j].re, inQ[j].im, MODE);
    mpfr_fma(bufc.im, tmpc[4][j].im, inQ[j].re, bufc.im, MODE);

    mpfr_mul(outQ[j].re, bufc.re, buf, MODE);
    mpfr_mul(outQ[j].im, bufc.im, buf, MODE);

    mpfr_mul(bufc.re, tmpc[0][j].im, tmpc[2][j].im, MODE);
    mpfr_fms(bufc.re, tmpc[0][j].re, tmpc[2][j].re, bufc.re, MODE);
    mpfr_mul(bufc.im, tmpc[0][j].re, tmpc[2][j].im, MODE);
    mpfr_fma(bufc.im, tmpc[0][j].im, tmpc[2][j].re, bufc.im, MODE); 
    mpfr_mul_ui(bufc.re, bufc.re, 2, MODE);
    mpfr_mul_ui(bufc.im, bufc.im, 2, MODE);

    mpfr_sub(outQ[j].re, bufc.re, outQ[j].re, MODE);
    mpfr_sub(outQ[j].im, bufc.im, outQ[j].im, MODE);

    mpfr_mul(outQ[j].re, outQ[j].re, buf, MODE);
    mpfr_mul(outQ[j].im, outQ[j].im, buf, MODE);

    //outV[j] = g*(inQ[j]*inQ[j]-1.L)+1.IL*conf.dq[j]*(tmpc[2][j]*tmpc[1][j]-inQ[j]*inQ[j]*tmpc[3][j]);
    //outV <- Q^2
    mpfr_mul(outV[j].re, inQ[j].im, inQ[j].im, MODE);
    mpfr_fms(outV[j].re, inQ[j].re, inQ[j].re, outV[j].re, MODE);
   
    mpfr_mul(outV[j].im, inQ[j].re, inQ[j].im, MODE);
    mpfr_mul_ui(outV[j].im, outV[j].im, 2, MODE);

    mpfr_mul(bufc.re, outV[j].re, tmpc[3][j].re, MODE);
    mpfr_fms(bufc.re, outV[j].im, tmpc[3][j].im, bufc.re, MODE);
    mpfr_fms(bufc.re, tmpc[2][j].im, tmpc[1][j].im, bufc.re, MODE);
    mpfr_fms(bufc.re, tmpc[2][j].re, tmpc[1][j].re, bufc.re, MODE);
    
    mpfr_mul(bufc.im,    outV[j].re, tmpc[3][j].im, MODE);
    mpfr_fms(bufc.im, tmpc[2][j].re, tmpc[1][j].im, bufc.im, MODE);
    mpfr_fms(bufc.im,    outV[j].im, tmpc[3][j].re, bufc.im, MODE);
    mpfr_fms(bufc.im, tmpc[2][j].im, tmpc[1][j].re, bufc.im, MODE);
    // bufc <- (tmpc[2][j]*tmpc[1][j] - Q^2*tmpc[3]) 
    // stopped here (bufc <- whats in bracket, outV <- Q^2)
    mpfr_mul(bufc.re, bufc.re, conf.dq[j], MODE);
    mpfr_mul(bufc.im, bufc.im, conf.dq[j], MODE);

    mpfr_set(buf, bufc.re, MODE);
    mpfr_neg(bufc.re, bufc.im, MODE);
    mpfr_set(bufc.im, buf, MODE);
 
    mpfr_sub_ui(outV[j].re, outV[j].re, 1, MODE);
    mpfr_mul(outV[j].re, outV[j].re, g, MODE);
    mpfr_mul(outV[j].im, outV[j].im, g, MODE);

    mpfr_add(outV[j].re, outV[j].re, bufc.re, MODE);
    mpfr_add(outV[j].im, outV[j].im, bufc.im, MODE);
  }
  // verification needed
  // ------------------------
  mpfr_clears(overN, sigma, g, buf, (mpfr_ptr) NULL);
  mpfr_clears(w1.re, w1.im, (mpfr_ptr) NULL);
  mpfr_clears(w2.re, w2.im, (mpfr_ptr) NULL);
  mpfr_clears(b1U.re, b1U.im, (mpfr_ptr) NULL);
  mpfr_clears(b2U.re, b2U.im, (mpfr_ptr) NULL);
  mpfr_clears(bufc.re, bufc.im, (mpfr_ptr) NULL);
}


void rk6_step(mpfc_t *inQ, mpfc_t *inV, mpfr_t dt) {
  long int	 	N = 1<<state.nbits;
  mpfr_t		overN;
  
  mpfr_init2(overN, state.precision);
  mpfr_set_ui(overN,     1, MODE);
  mpfr_div_si(overN, overN, N, MODE);

  //memcpy(tQ, inQ, N*sizeof(mpfc_t)); 
  //memcpy(tV, inV, N*sizeof(mpfc_t)); 
  for (long int j = 0; j < N; j++) {
    mpfr_set(tQ[j].re, inQ[j].re, MODE);
    mpfr_set(tQ[j].im, inQ[j].im, MODE);
    mpfr_set(tV[j].re, inV[j].re, MODE);
    mpfr_set(tV[j].im, inV[j].im, MODE);
  }
  // start here on Wednesday.
  
  //complex_array_out("inV.txt", inV);
  //complex_array_out("inQ.txt", inQ);

  compute_rhs(tQ, tV, kq[0], kv[0]); // working here
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
    tQ[j] = inQ[j] + one_fortyfourths*dt*(9.0L*kq[0][j] - 36.0L*kq[1][j] + 63.0L*kq[2][j] + 72.0L*kq[3][j] - 64.0L*kq[4][j]);
    tV[j] = inV[j] + one_fortyfourths*dt*(9.0L*kv[0][j] - 36.0L*kv[1][j] + 63.0L*kv[2][j] + 72.0L*kv[3][j] - 64.0L*kv[4][j]);
  }
  compute_rhs(tQ, tV, kq[6], kv[6]);
  for (long int j = 0; j < N; j++) {
    inQ[j] = (inQ[j] + one_onetwentieth*dt*(11.0L*kq[0][j] + 81.0L*kq[2][j] + 81.0L*kq[3][j] - 32.0L*kq[4][j] - 32.0L*kq[5][j] + 11.0L*kq[6][j]))*overN;
    inV[j] = (inV[j] + one_onetwentieth*dt*(11.0L*kv[0][j] + 81.0L*kv[2][j] + 81.0L*kv[3][j] - 32.0L*kv[4][j] - 32.0L*kv[5][j] + 11.0L*kv[6][j]))*overN;
  }
  memcpy(tmpc[0], inQ, N*sizeof(mpfc_t));
  memcpy(tmpc[1], inV, N*sizeof(mpfc_t));
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  memset(tmpc[0]+N/2, 0, N/2*sizeof(mpfc_t));
  memset(tmpc[1]+N/2, 0, N/2*sizeof(mpfc_t));
  for (long int j = 0; j < N/2; j++) {
    tmpc[0][j] = tmpc[0][j]*cexpl(-1.L*powl(1.L*j/state.kD, pD)); // dt
    tmpc[1][j] = tmpc[1][j]*cexpl(-1.L*powl(1.L*j/state.kD, pD)); // dt
    if (j >= kZ) {
      tmpc[0][j] = 0.L;
      tmpc[1][j] = 0.L;
    }
  }
  fftwl_execute(ft0);
  fftwl_execute(ft1); 
  memcpy(inQ, tmpc[0], N*sizeof(mpfc_t));
  memcpy(inV, tmpc[1], N*sizeof(mpfc_t));
}



void evolve_rk6() {
  unsigned int		QC_pass = 1;
  unsigned long 	counter = 0, j = 0, skip = 10;
  unsigned long		ref_counter = 0;
  char 			filename1[80], filename2[80];

  mpfr_t		M_TOL, R_TOL;
  mpfr_t   	        time, tshift, Ham, dt;
  //mpfr_t   	        dt = cfl*2.L*PI*conf.scaling/state.number_modes;
  
  mpfr_inits2(state.precision, M_TOL, R_TOL, tshift, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, time, Ham, dt, (mpfr_ptr) NULL);
  mpfr_set_d(M_TOL, 2.0e-15, MODE);
  mpfr_set_d(R_TOL, 4.0e-16, MODE);
  mpfr_set_ui(time,    0, MODE);
  mpfr_set_ui(tshift,  0, MODE);
  mpfr_set_ui(Ham,     0, MODE);
  mpfr_mul(dt, Pie, conf.scaling, MODE);
  mpfr_mul(dt, cfl, dt, MODE);
  mpfr_div_ui( dt, dt, 1<<(state.nbits-1), MODE);


  FILE *fh_time	= fopen("time_dependence.txt","w");
  fprintf(fh_time, "# 1. time 2. Kinetic 3. Potential 4. Momentum X 5. Momentum Y\n\n");
  fclose(fh_time);

  map_quality_fourier(data[0], data[1], M_TOL, &QC_pass);
  sprintf(filename2, "./data/spec_%04lu.txt", counter);
  spec_out(filename2, tmpc[0], tmpc[1]);
 
  convertQtoZ(data[0], tmpc[5]);
  sprintf(filename1, "./data/surf_%04lu.txt", counter);
  surface_out(filename1, tmpc[5]);  // checks out
  
  if (PADE_TEST) {
    exit(1);
  }
  restore_potential(data[0], data[1], tmpc[5]);  
  
  fh_time = fopen("time_dependence.txt","a");
  mpfr_fprintf(fh_time, "%.17Re\t%.17Re\t%.17Re\t", state.time, state.kineticE, state.potentialE); 
  mpfr_fprintf(fh_time, "%.17Re\t%.17Re\n", state.momentum.im, state.momentum.re); 
  fclose(fh_time);
  mpfr_add(Ham, state.kineticE, state.potentialE, MODE);
  sprintf(filename1, "./aux/data_%04lu.txt", counter);
  output_data(filename1, tmpc[5]);
  print_constants();       // all constants checked out (mean level, potential+kinetic energy)
  //mpfr_printf("T = %23.16Re\tH = %23.16Re\n", state.time, Ham);
  
  if (QC_pass == 0) {
    printf("Bad quality map at start.\tStop!\n");
    exit(1);
  } else if (QC_pass == 2) {
    printf("Warning! Map is over-resolved at start.\n");
  }
  
  while (QC_pass) {
    rk6_step(data[0], data[1], dt);  
    //time = (j+1)*dt;
    mpfr_mul_ui(time, dt, j+1, MODE);
    j++;
    //state.time = time + tshift;
    mpfr_add(state.time, time, tshift, MODE);
/*
    map_quality_fourier(data[0], data[1], M_TOL, &QC_pass); 
    if (QC_pass == 0) {
      ref_counter++;
      for (unsigned long j = 0; j < state.number_modes; j++) tmpc[5][j] = data[0][j]*data[0][j];
      sprintf(filename2, "./roots_G/roots_%04lu.txt", ref_counter);
      optimal_pade(filename2, tmpc[5]);
      spec_out("last.spec.txt", tmpc[0], tmpc[1]);
      restore_potential(data[0], data[1], tmpc[2]);
      // Attempt to fix the conformal map
      remap(&alt_map, state.number_modes);
      restore_potential(data[0], data[1], tmpc[2]);
      map_quality_fourier(data[0], data[1], R_TOL, &QC_pass); 
      if (QC_pass == 0) {
        printf("Doubling # of Modes: %lu\n", 2*state.number_modes);
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
        //sprintf(filename2, "./roots/roots_Z%04lu.txt", counter);
        //optimal_pade(filename2, tmpc[5]);
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
	//
        //sprintf(filename2, "./roots/roots_P%04lu.txt", counter);
        //optimal_pade(filename2, tmpc[5]);
 	//
        // write out cut for derivative of potential
        memcpy(tmpc[0], tmpc[5], state.number_modes*sizeof(mpfc_t));        
        fftwl_execute(ift0);
        memset(tmpc[0] + state.number_modes/2, 0, state.number_modes*sizeof(mpfc_t)/2);
        for (unsigned long j = 0; j < state.number_modes/2; j++) tmpc[0][j] = -1.0IL*j*tmpc[0][j]/state.number_modes; // this is derivative dq
        fftwl_execute(ft0);
        for (unsigned long j = 0; j < state.number_modes; j++) tmpc[0][j] = tmpc[0][j]*conf.dq[j]; // this is derivative du

        sprintf(filename2, "./roots/roots_dP%04lu.txt", counter);
        optimal_pade(filename2, tmpc[0]);
	// write out cut for R
        //for (unsigned long j = 0; j < state.number_modes; j++) tmpc[5][j] = data[0][j]*data[0][j];
        //sprintf(filename2, "./roots/roots_R%04lu.txt", counter);
        //optimal_pade(filename2, tmpc[5]);
	// write out cut for Zu
        for (unsigned long j = 0; j < state.number_modes; j++) tmpc[5][j] = 1.L/(data[0][j]*data[0][j]);
        sprintf(filename2, "./roots/roots_DZ%04lu.txt", counter);
        optimal_pade(filename2, tmpc[5]);
	// write out cut for V
        //for (unsigned long j = 0; j < state.number_modes; j++) tmpc[5][j] = data[1][j];
        //sprintf(filename2, "./roots/roots_V%04lu.txt", counter);
        //optimal_pade(filename2, tmpc[5]);
      }
    }
    */
  }
  printf("Simulation Stops\n");
  exit(0);
}

