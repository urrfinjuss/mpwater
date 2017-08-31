#include "header.h"

/*void compute_rhs(mpfc_t *inQ, mpfc_t *inV, mpfc_t *outQ, mpfc_t *outV) {
  long double overN = 1.L/state.number_modes;
}*/

/*
void project(mpfc_t *in, mpfc_t *out) {
  long double 		overN = 1.L/state.number_modes;
  mpfc_t 	b0 = 0.L;	

  memcpy(tmpc[0], in, state.number_modes*sizeof(mpfc_t));
  fftwl_execute(ift0);
  // first we need to compute the proper value of the constant:  
  for (long int j = state.number_modes/2 - 2; j > -1; j--) {
    b0 += 0.5L*(tmpc[0][state.number_modes-j-1]*conjl(conf.w[j]) - tmpc[0][j+1]*conf.w[j]);
  }
  memset(tmpc[0]+state.number_modes/2, 0, state.number_modes/2*sizeof(mpfc_t));
  for (long int j = 1; j < state.number_modes/2; j++) {
    tmpc[0][j] = tmpc[0][j]*overN;
  }
  tmpc[0][0] = 0.5L*tmpc[0][0]*overN;
  fftwl_execute(ft0);
  for (long int j = 0; j < state.number_modes; j++) {
    out[j] = tmpc[0][j] + b0;
  } 
  //printf("shift_re = %.16LE\nshift_im = %.16LE\n", creall(b0), cimagl(b0));
}
*/


void restore_potential(mpfc_t *inQ, mpfc_t *inV, mpfc_t *out) {
  //long double overN = 1.L/state.number_modes;
  //long double K = 0.L;
  //mpfc_t P = 0.L, z0 = 0.L;
  long int N = 1<<state.nbits;
  mpfr_t 	overN, K, buf;
  mpfc_t 	P, z0;
  mpfr_inits2(state.precision, overN, buf, K, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, P.re, P.im, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, z0.re, z0.im, (mpfr_ptr) NULL);
  
  mpfr_set_ui(overN, 1, MODE);
  mpfr_div_si(overN, overN, N, MODE);
  mpfr_set_ui(K, 0, MODE);
  mpfr_set_ui(P.re, 0, MODE);
  mpfr_set_ui(P.im, 0, MODE);
  mpfr_set_ui(z0.re, 0, MODE);
  mpfr_set_ui(z0.im, 0, MODE);

  convertQtoZ(inQ, tmpc[5]); // tmpc[5] <- z-tilde
  if (inQ == tmpc[0]) {
    printf("inQ = tmpc[0]: Copy through buffer needed.\n");
    exit(1);
  }
  if (inV == tmpc[1]) {
    printf("inV = tmpc[1]: Copy through buffer needed.\n");
    exit(1);
  } 
  for (long int j = 0; j < N; j++) {
    //tmpc[0][j] = inQ[j]*inQ[j]*overN;
    mpfr_mul(tmpc[0][j].re, inQ[j].im, inQ[j].im, MODE);
    mpfr_fms(tmpc[0][j].re, inQ[j].re, inQ[j].re, tmpc[0][j].re, MODE);
    mpfr_mul(tmpc[0][j].im, inQ[j].re, inQ[j].im, MODE);
    mpfr_fma(tmpc[0][j].im, inQ[j].im, inQ[j].re, tmpc[0][j].im, MODE);
    mpfr_mul(tmpc[0][j].re, tmpc[0][j].re, overN, MODE);
    mpfr_mul(tmpc[0][j].im, tmpc[0][j].im, overN, MODE);
    //tmpc[1][j] = -1.IL*inV[j]*overN;
    mpfr_mul(tmpc[1][j].re, inV[j].im, overN, MODE);
    mpfr_mul(tmpc[1][j].im, inV[j].re, overN, MODE);
    mpfr_neg(tmpc[1][j].im, tmpc[1][j].im, MODE);
    //tmpc[2][j] = tmpc[5][j]*overN;
    mpfr_mul(tmpc[2][j].re, tmpc[5][j].re, overN, MODE);
    mpfr_mul(tmpc[2][j].im, tmpc[5][j].im, overN, MODE);
  }  
  //fftwl_execute(ift0);
  mpfft_execute(ift0);
  //fftwl_execute(ift1);
  mpfft_execute(ift1);
  //fftwl_execute(ift2);
  mpfft_execute(ift2);
  //memcpy(tmpc[5], tmpc[2], state.number_modes*sizeof(mpfc_t)); 
  for (long int j = 0; j < N; j++) {
    mpfr_set(tmpc[5][j].re, tmpc[2][j].re, MODE);
    mpfr_set(tmpc[5][j].im, tmpc[2][j].im, MODE);
  }
  linear_solve(tmpc[0], tmpc[1], tmpc[2]);
  //memcpy(tmpc[1], tmpc[2], state.number_modes*sizeof(mpfc_t));
  for (long int j = 0; j < N; j++) {
    mpfr_set(tmpc[1][j].re, tmpc[2][j].re, MODE);
    mpfr_set(tmpc[1][j].im, tmpc[2][j].im, MODE);
  }
  div_jacobian(tmpc[1], tmpc[0]);
  for (long int j = 1; j < N/2; j++) {
    //tmpc[0][j] = 1.0IL*tmpc[0][j]/j;
    mpfr_set(buf, tmpc[0][j].re, MODE);
    mpfr_div_si(tmpc[0][j].re, tmpc[0][j].im, -j, MODE);
    mpfr_div_si(tmpc[0][j].im, buf, j, MODE);
  }
  //memset(tmpc[0]+state.number_modes/2, 0, state.number_modes*sizeof(mpfc_t)/2);
  for (long int j = 0; j < N/2; j++) {
    mpfr_set_ui(tmpc[0][N/2+j].re, 0, MODE);
    mpfr_set_ui(tmpc[0][N/2+j].im, 0, MODE);
  }
  // P is zero
  compute_zero_mode_complex(tmpc[0], P, &tmpc[0][0]);
  //tmpc[0][0] = z0;
  for (long int j = N/2-1; j > 0; j--) {
    //K += cimagl(tmpc[0][j]*conjl(1.IL*j*tmpc[0][j]));
    mpfr_mul(buf, tmpc[0][j].im, tmpc[0][j].im, MODE);
    mpfr_fma(buf, tmpc[0][j].re, tmpc[0][j].re, buf, MODE);
    mpfr_mul_si(buf, buf, j, MODE);
    mpfr_add(K, K, buf, MODE); 
    //P += -0.5IL*j*tmpc[5][j]*conjl(tmpc[0][j]);
    // real
    mpfr_mul(buf, tmpc[5][j].re, tmpc[0][j].im, MODE);
    mpfr_fms(buf, tmpc[5][j].im, tmpc[0][j].re, buf, MODE);
    mpfr_mul_si(buf, buf, j, MODE);
    mpfr_add(P.re, P.re, buf, MODE);
    // imag
    mpfr_mul(buf, tmpc[5][j].re, tmpc[0][j].re, MODE);
    mpfr_fma(buf, tmpc[5][j].im, tmpc[0][j].im, buf, MODE);
    mpfr_mul_si(buf, buf, -j, MODE);
    mpfr_add(P.im, P.im, buf, MODE);
  }
  //K = -1.L*PI*K;  
  //mpfr_mul(K, K, Pie, MODE);
  //state.kineticE = K;
  mpfr_set(state.kineticE, K, MODE);
  //state.momentum = P;
  mpfr_div_ui(state.momentum.re, P.re, 2, MODE);
  mpfr_div_ui(state.momentum.im, P.im, 2, MODE);
  //fftwl_execute(ft0);
  mpfft_execute(ft0);
  //memcpy(out, tmpc[0], state.number_modes*sizeof(mpfc_t));
  for (long int j = 0; j < N; j++) {
    mpfr_set(out[j].re, tmpc[0][j].re, MODE);
    mpfr_set(out[j].im, tmpc[0][j].im, MODE);
  }
  mpfr_clears(overN, K, buf, (mpfr_ptr) NULL);
  mpfr_clears(P.re, P.im, (mpfr_ptr) NULL);
  mpfr_clears(z0.re, z0.im, (mpfr_ptr) NULL);
}



void convertQtoZ(mpfc_t *in, mpfc_t *out) {
  // computes Z from Q by series inversion of equation:	   //
  // 		   z_u Q^2 = 1				   //
  if (in == out) printf("Warning In-Place Q to Z undefined!\n");
  long int N = 1<<state.nbits;
  mpfr_t  overN, S0, T0, P;
  mpfr_t  mean_level, buf;
  mpfc_t  z0, cS0;
  
  mpfr_inits2(state.precision, overN, S0, T0, P, buf, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, z0.re, z0.im, mean_level, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, cS0.re, cS0.im, (mpfr_ptr) NULL);

  mpfr_set_ui(overN, 1, MODE);
  mpfr_div_ui(overN, overN, N, MODE);
  mpfr_set_ui(S0, 0, MODE);
  mpfr_set_ui(T0, 0, MODE);
  mpfr_set_ui(P,  0, MODE);
  mpfr_set_ui(z0.re, 0, MODE);
  mpfr_set_ui(z0.im, 0, MODE);
  mpfr_set_ui(mean_level, 0, MODE);

  for (long int j = 0; j < N; j++) {
    //tmpc[0][j] = in[j]*in[j]*overN;
    mpfr_mul(tmpc[0][j].re, in[j].im, in[j].im, MODE);
    mpfr_fms(tmpc[0][j].re, in[j].re, in[j].re, tmpc[0][j].re, MODE);

    mpfr_mul(tmpc[0][j].im, in[j].im, in[j].re, MODE);
    mpfr_fma(tmpc[0][j].im, in[j].re, in[j].im, tmpc[0][j].im, MODE);

    mpfr_mul(tmpc[0][j].re, tmpc[0][j].re, overN, MODE);
    mpfr_mul(tmpc[0][j].im, tmpc[0][j].im, overN, MODE);
  }
  //fftwl_execute(ift0);
  mpfft_execute(ift0);
  inverse(tmpc[0],tmpc[1]);
  //fftwl_execute(ft1);
  mpfft_execute(ft1);
  for (long int j = 0; j < N; j++) {
    //tmpc[1][j] = (tmpc[1][j] - 1.L)*overN;
    mpfr_sub_ui(tmpc[1][j].re, tmpc[1][j].re, 1, MODE);
    mpfr_mul(tmpc[1][j].re, tmpc[1][j].re, overN, MODE);
    mpfr_mul(tmpc[1][j].im, tmpc[1][j].im, overN, MODE);
  }
  //fftwl_execute(ift1);
  mpfft_execute(ift1);
  div_jacobian(tmpc[1], tmpc[0]);
  for (long int j = N/2-1; j > 0; j--) {
    //tmpc[0][j] =  1.0IL*tmpc[0][j]/j;	// this stores z_k, k = 1, ..., N/2
    mpfr_set(buf, tmpc[0][j].re, MODE);
    mpfr_div_si(tmpc[0][j].re, tmpc[0][j].im, -j, MODE);
    mpfr_div_si(tmpc[0][j].im, buf,            j, MODE); 

    //S0 += -0.5L*j*creall(tmpc[0][j]*conjl(tmpc[0][j])); 
    mpfr_mul(buf, tmpc[0][j].im, tmpc[0][j].im, MODE);
    mpfr_fma(buf, tmpc[0][j].re, tmpc[0][j].re, buf, MODE);
    mpfr_mul_si(buf, buf, j, MODE);
    mpfr_div_ui(buf, buf, 2, MODE);
    mpfr_sub(S0, S0, buf, MODE);
   
    //T0 += -creall(tmpc[0][j]);
    mpfr_sub(T0, T0, tmpc[0][j].re, MODE);
  }
 
  //compute_zero_mode_complex(tmpc[0], 1.IL*S0, &z0);
  mpfr_set_ui(cS0.re,  0, MODE);
  mpfr_set   (cS0.im, S0, MODE);
  compute_zero_mode_complex(tmpc[0], cS0, &z0);

  //tmpc[0][0] = T0 - creall(z0) + z0;
  mpfr_set(tmpc[0][0].re, T0, MODE);
  mpfr_set(tmpc[0][0].im, z0.im, MODE);

  //memset(tmpc[0]+N/2, 0, N/2*sizeof(mpfc_t));
  for (long int j = 0; j < N/2; j++) {
    mpfr_set_ui(tmpc[0][j+N/2].re, 0, MODE);
    mpfr_set_ui(tmpc[0][j+N/2].im, 0, MODE);
  }
  //tmpc[1][0] = cimagl(z0);
  mpfr_set(tmpc[1][0].re, z0.im, MODE);
  mpfr_set_ui(tmpc[1][0].im,  0, MODE);
  for (long int j = 1; j < N/2; j++) {
    //tmpc[1][j] = -0.5IL*tmpc[0][j];
    mpfr_div_si(tmpc[1][j].re, tmpc[0][j].im,  2, MODE);
    mpfr_div_si(tmpc[1][j].im, tmpc[0][j].re, -2, MODE);
    //tmpc[1][N-j] = conjl(tmpc[1][j]);
    mpfr_set(tmpc[1][N-j].re, tmpc[1][j].re, MODE);
    mpfr_neg(tmpc[1][N-j].im, tmpc[1][j].im, MODE);
  }
  div_jacobian(tmpc[1], tmpc[4]); // now S
  //fftwl_execute(ft0);
  mpfft_execute(ft0);
  //memcpy(out, tmpc[0], N*sizeof(mpfc_t));
  for (long int j = 0; j < N; j++) {
    mpfr_set(out[j].re, tmpc[0][j].re, MODE);
    mpfr_set(out[j].im, tmpc[0][j].im, MODE);
  }
  for (long int j = 0; j < N; j++) {
    //tmpc[2][j] = cpowl(cimagl(tmpc[0][j]),2)*overN;
    mpfr_mul(tmpc[2][j].re, tmpc[0][j].im, tmpc[0][j].im, MODE);
    mpfr_mul(tmpc[2][j].re, tmpc[2][j].re, overN, MODE);
    mpfr_set_ui(tmpc[2][j].im, 0, MODE);
    //tmpc[0][j] = cimagl(tmpc[0][j])*overN;
    mpfr_mul(tmpc[0][j].re, tmpc[0][j].im, overN, MODE);
    mpfr_set_ui(tmpc[0][j].im, 0, MODE);
  }
  //fftwl_execute(ift0);
  mpfft_execute(ift0);
  //fftwl_execute(ift2);
  mpfft_execute(ift2);
  for (long int j = N/2-1; j > 0; j--) {
    //mean_level += 2.0L*j*creall(tmpc[0][j]*conjl(tmpc[0][j]));
    mpfr_mul(buf, tmpc[0][j].im, tmpc[0][j].im, MODE);
    mpfr_fma(buf, tmpc[0][j].re, tmpc[0][j].re, buf, MODE);
    mpfr_mul_si(buf, buf, 2*j, MODE);
    mpfr_add(mean_level, mean_level, buf, MODE);
    //P += 2.L*creall((tmpc[4][j] + j*tmpc[2][j])*conjl(tmpc[1][j]));
    mpfr_mul_si(buf, tmpc[2][j].re, j, MODE);
    mpfr_add(buf, buf, tmpc[4][j].re, MODE);
    mpfr_fma(P, buf, tmpc[1][j].re, P, MODE);

    mpfr_mul_si(buf, tmpc[2][j].im, j, MODE);
    mpfr_add(buf, buf, tmpc[4][j].im, MODE);
    mpfr_fma(P, buf, tmpc[1][j].im, P, MODE);
  }
  //mpfr_printf("P = %.15Re\n", P);  
  mpfr_mul_ui(P, P, 2, MODE);

  //mean_level = 2.L*PI*(S0 + mean_level);
  mpfr_add(buf, S0, mean_level, MODE);
  mpfr_mul(mean_level, buf, Pie, MODE);
  mpfr_mul_ui(mean_level, mean_level, 2, MODE);

  //P += creall(tmpc[4][0]*conjl(tmpc[1][0]));
  mpfr_fma(P, tmpc[4][0].re, tmpc[1][0].re, P, MODE);
  mpfr_fma(P, tmpc[4][0].im, tmpc[1][0].im, P, MODE);

  //state.potentialE = 2.L*PI*state.gravity*P;
  //mpfr_mul   (state.potentialE, state.gravity,  Pie, MODE);
  mpfr_mul_ui(state.potentialE, state.gravity, 2, MODE);
  mpfr_mul   (state.potentialE, state.potentialE, P, MODE); 

  //state.mean_level = mean_level;
  mpfr_set(state.mean_level, mean_level, MODE);

  // clear all
  mpfr_clears(overN, S0, T0, P, buf, (mpfr_ptr) NULL);
  mpfr_clears(z0.re, z0.im, mean_level, (mpfr_ptr) NULL);
  mpfr_clears(cS0.re, cS0.im, (mpfr_ptr) NULL);
}



void convertZtoQ(mpfc_t *in, mpfc_t *out) {
  // computes Q from tilde-Z (z) by series inversion of:   //
  // 		   (z_q q_u + 1) Q^2 = 1		   //
  // 							   //
  // in		-- array with Z-tilde (z)		   //
  // out	-- array with Q				   //
  unsigned long N = 1<<state.nbits;
  mpfr_t overN, eye;
  mpfr_inits2(state.precision, eye, overN, (mpfr_ptr) NULL);
  mpfr_set_ui(eye,     1, MODE);
  mpfr_div_ui(overN, eye, N, MODE);

  //memcpy(tmpc[0], in, N*sizeof(mpfc_t));
  for (unsigned long j = 0; j < N; j++) {
    mpfr_set(tmpc[0][j].re, in[j].re, MODE);
    mpfr_set(tmpc[0][j].im, in[j].im, MODE);
  }
  //fftwl_execute(ift0);
  mpfft_execute(ift0);
  spec_out("specZ1.txt", tmpc[0], tmpc[0]);
  for (unsigned long j = 0; j < N/2; j++) {
    //tmpc[0][j] =  -1.0IL*((mpfc_t)(tmpc[0][j]*j))*overN;
    mpfr_swap(tmpc[0][j].re, tmpc[0][j].im); 
    //mpfr_neg(tmpc[0][j].im, tmpc[0][j].im, MODE);
    mpfr_mul_si(tmpc[0][j].re, tmpc[0][j].re,  j, MODE);
    mpfr_mul_si(tmpc[0][j].im, tmpc[0][j].im, -j, MODE);
    mpfr_mul(tmpc[0][j].re, tmpc[0][j].re, overN, MODE);
    mpfr_mul(tmpc[0][j].im, tmpc[0][j].im, overN, MODE);
    //tmpc[0][N-1-j] = (mpfc_t) 0.0L;
    mpfr_set_ui(tmpc[0][N-1-j].re, 0, MODE);
    mpfr_set_ui(tmpc[0][N-1-j].im, 0, MODE);
  }
  spec_out("specZ2.txt", tmpc[0], tmpc[0]);
  //memset(tmpc[0]+N/2, 0, N/2*sizeof(mpfc_t));
  //fftwl_execute(ft0); 
  mpfft_execute(ft0);
  complex_array_out("dZdq.txt", tmpc[0]);
  for (unsigned long j = 0; j < N; j++) {
    //tmpc[0][j] = (tmpc[0][j]*conf.dq[j] + 1.L)*overN; 
    mpfr_fma(tmpc[0][j].re, tmpc[0][j].re, conf.dq[j], eye, MODE);
    mpfr_mul(tmpc[0][j].im, tmpc[0][j].im, conf.dq[j], MODE);
    mpfr_mul(tmpc[0][j].re, tmpc[0][j].re, overN, MODE); // uncomment after
    mpfr_mul(tmpc[0][j].im, tmpc[0][j].im, overN, MODE);
  }
  //fftwl_execute(ift0);
  complex_array_out("dZdu.txt", tmpc[0]); // verified 
  mpfft_execute(ift0);
  spec_out("dZdu_spec.txt", tmpc[0], tmpc[0]); 
  inverse(tmpc[0], tmpc[1]);  	// start here (verified)
  // ---
  square_ft(tmpc[1], tmpc[0]); 	// and here   (verified)
  // ---
  //fftwl_execute(ft0);
  mpfft_execute(ft0);
  complex_array_out("over_sqrt_dZ.txt", tmpc[0]); 
  //memcpy(out, tmpc[0], N*sizeof(mpfc_t));
  for (unsigned long j = 0; j < N; j++) {
    mpfr_set(out[j].re, tmpc[0][j].re, MODE);
    mpfr_set(out[j].im, tmpc[0][j].im, MODE);
  }
  mpfr_clears(eye, overN, (mpfr_ptr) NULL);
}
