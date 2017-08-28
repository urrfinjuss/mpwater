#include "header.h"

/*
void compute_zero_mode(fftwl_complex *in, long double S0, long double *out) {
  fftwl_complex 	w = cexpl(1.IL*conf.origin_offset);
  long double 		b = 0.5L*(1.L + powl(conf.scaling, 2))/conf.scaling;
  long double 		xi = (1.L - powl(conf.scaling, 2))/(1.L + powl(conf.scaling, 2));
  tmpc[2][0] = -0.5L*xi*conjl(w);
  for (long int j = 1; j < state.number_modes/2-1; j++) {
    tmpc[2][j] = tmpc[2][0]/(1.L - conjl(tmpc[2][0])*tmpc[2][j-1]); 
  }
  tmpc[3][0] = in[1]/b - S0*conjl(tmpc[2][0]);  // added I 
  for (long int j = 1; j < state.number_modes/2-1; j++) {
    tmpc[3][j] = in[j+1]/b - conjl(tmpc[2][0])*tmpc[3][j-1];
    tmpc[3][j] = tmpc[3][j]/(1.L - conjl(tmpc[2][0])*tmpc[2][j-1]);
  } 
  tmpc[4][state.number_modes/2-2] = tmpc[3][state.number_modes/2-2];
  for (long int j = state.number_modes/2-3; j > -1; j--) {
    tmpc[4][j] = tmpc[3][j]-tmpc[2][j]*tmpc[4][j+1];
  }
  *out = b*creall(2.0L*tmpc[4][0]*tmpc[2][0] + S0);
}
*/


void compute_zero_mode_complex(mpfc_t *in, mpfc_t S0, mpfc_t *out) {
  //fftwl_complex 	w = cexpl(1.IL*conf.origin_offset);
  long int N2 = 1<<(state.nbits-1);
  mpfc_t w;
  mpfr_t xi, b;
  mpfr_t buf1, buf2, buf3, buf4, buf5;

  mpfr_inits2(state.precision, w.re, w.im, xi, b, buf1, buf2, buf3, buf4, buf5, (mpfr_ptr) NULL);
  mpfr_sin_cos(w.im, w.re, conf.origin_offset, MODE);
  
  
  //long double 		b = 0.5L*(1.L + powl(conf.scaling, 2))/conf.scaling;
  //long double 		xi = (1.L - powl(conf.scaling, 2))/(1.L + powl(conf.scaling, 2));
  mpfr_mul(buf1, conf.scaling, conf.scaling, MODE);
  mpfr_add_ui(buf2, buf1, 1, MODE);
  mpfr_ui_sub(buf1, 1, buf1, MODE);
  mpfr_div(b, buf2, conf.scaling, MODE);
  mpfr_div_ui(b, b, 2, MODE);
  mpfr_div(xi, buf1, buf2, MODE);
  

  //tmpc[2][0] = -0.5L*xi*conjl(w);
  mpfr_mul(tmpc[2][0].re, xi, w.re, MODE);
  mpfr_mul(tmpc[2][0].im, xi, w.im, MODE);
  mpfr_div_si(tmpc[2][0].re, tmpc[2][0].re, -2, MODE);
  mpfr_div_si(tmpc[2][0].im, tmpc[2][0].im,  2, MODE);
  
  for (long int j = 1; j < N2-1; j++) {
    //tmpc[2][j] = tmpc[2][0]/(1.L - conjl(tmpc[2][0])*tmpc[2][j-1]); 
    mpfr_mul(buf1, tmpc[2][0].re, tmpc[2][j-1].re, MODE);
    mpfr_fma(buf1, tmpc[2][0].im, tmpc[2][j-1].im, buf1, MODE);
    mpfr_mul(buf2, tmpc[2][0].im, tmpc[2][j-1].re, MODE);
    mpfr_fms(buf2, tmpc[2][0].re, tmpc[2][j-1].im, buf2, MODE);

    mpfr_ui_sub(buf1, 1, buf1, MODE);
    mpfr_neg(buf2, buf2, MODE);

    mpfr_mul(buf3, buf1, buf1, MODE);
    mpfr_fma(buf3, buf2, buf2, buf3, MODE);

    mpfr_mul(tmpc[2][j].re, buf1, tmpc[2][0].re, MODE);
    mpfr_fma(tmpc[2][j].re, buf2, tmpc[2][0].im, tmpc[2][j].re, MODE);
    mpfr_mul(tmpc[2][j].im, buf2, tmpc[2][0].re, MODE);
    mpfr_fms(tmpc[2][j].im, buf1, tmpc[2][0].im, tmpc[2][j].im, MODE);

    mpfr_div(tmpc[2][j].re, tmpc[2][j].re, buf3, MODE);
    mpfr_div(tmpc[2][j].im, tmpc[2][j].im, buf3, MODE);
  }
  //tmpc[3][0] = in[1]/b - 2.L*conjl(tmpc[2][0])*S0;
  mpfr_mul(buf1, tmpc[2][0].re, S0.re, MODE);
  mpfr_fma(buf1, tmpc[2][0].im, S0.im, buf1, MODE);
  mpfr_mul(buf2, tmpc[2][0].im, S0.re, MODE);
  mpfr_fms(buf2, tmpc[2][0].re, S0.im, buf2, MODE);
  
  mpfr_mul_ui(buf1, buf1, 2, MODE);
  mpfr_mul_ui(buf2, buf2, 2, MODE);
  
  mpfr_div(tmpc[3][0].re, in[1].re, b, MODE);
  mpfr_div(tmpc[3][0].im, in[1].im, b, MODE);
  mpfr_sub(tmpc[3][0].re, tmpc[3][0].re, buf1, MODE);
  mpfr_sub(tmpc[3][0].im, tmpc[3][0].im, buf2, MODE);
  
  for (long int j = 1; j < N2-1; j++) {
    //tmpc[3][j] = in[j+1]/b - conjl(tmpc[2][0])*tmpc[3][j-1];
    mpfr_mul(buf1, tmpc[2][0].re, tmpc[3][j-1].re, MODE);
    mpfr_fma(buf1, tmpc[2][0].im, tmpc[3][j-1].im, buf1, MODE);
    mpfr_mul(buf2, tmpc[2][0].im, tmpc[3][j-1].re, MODE);
    mpfr_fms(buf2, tmpc[2][0].re, tmpc[3][j-1].im, buf2, MODE);

    mpfr_div(tmpc[3][j].re, in[j+1].re, b, MODE);
    mpfr_div(tmpc[3][j].im, in[j+1].im, b, MODE);
    mpfr_sub(tmpc[3][j].re, tmpc[3][j].re, buf1, MODE);
    mpfr_sub(tmpc[3][j].im, tmpc[3][j].im, buf2, MODE);

    //tmpc[3][j] = tmpc[3][j]/(1.L - conjl(tmpc[2][0])*tmpc[2][j-1]);
    mpfr_mul(buf1, tmpc[2][0].re, tmpc[2][j-1].re, MODE);
    mpfr_fma(buf1, tmpc[2][0].im, tmpc[2][j-1].im, buf1, MODE);
    mpfr_mul(buf2, tmpc[2][0].im, tmpc[2][j-1].re, MODE);
    mpfr_fms(buf2, tmpc[2][0].re, tmpc[2][j-1].im, buf2, MODE);

    mpfr_ui_sub(buf1, 1, buf1, MODE);
    mpfr_neg(buf2, buf2, MODE);    

    mpfr_mul(buf3, buf1, buf1, MODE);
    mpfr_fma(buf3, buf2, buf2, buf3, MODE);

    mpfr_mul(buf4, buf1, tmpc[3][j].re, MODE);
    mpfr_fma(buf4, buf2, tmpc[3][j].im, buf4, MODE);
    mpfr_mul(buf5, buf2, tmpc[3][j].re, MODE);
    mpfr_fms(buf5, buf1, tmpc[3][j].im, buf5, MODE);

    mpfr_div(tmpc[3][j].re, buf4, buf3, MODE);
    mpfr_div(tmpc[3][j].im, buf5, buf3, MODE);
  } 
  //tmpc[4][state.number_modes/2-2] = tmpc[3][state.number_modes/2-2];
  mpfr_set(tmpc[4][N2 - 2].re, tmpc[3][N2 - 2].re, MODE);
  mpfr_set(tmpc[4][N2 - 2].im, tmpc[3][N2 - 2].im, MODE);
  //for (long int j = state.number_modes/2-3; j > -1; j--) {
  for (long j = N2-3; j > -1; j--) {
    //tmpc[4][j] = tmpc[3][j]-tmpc[2][j]*tmpc[4][j+1];
    mpfr_mul(buf1, tmpc[2][j].im, tmpc[4][j+1].im, MODE);
    mpfr_fms(buf1, tmpc[2][j].re, tmpc[4][j+1].re, buf1, MODE);
    mpfr_mul(buf2, tmpc[2][j].re, tmpc[4][j+1].im, MODE);
    mpfr_fma(buf2, tmpc[2][j].im, tmpc[4][j+1].re, buf2, MODE);

    mpfr_sub(tmpc[4][j].re, tmpc[3][j].re, buf1, MODE);
    mpfr_sub(tmpc[4][j].im, tmpc[3][j].im, buf2, MODE);
  }
  //*out = b*(tmpc[4][0]*tmpc[2][0] + 1.L*S0);
  mpfr_mul(buf1, tmpc[4][0].im, tmpc[2][0].im, MODE);
  mpfr_fms(buf1, tmpc[4][0].re, tmpc[2][0].re, buf1, MODE);
  mpfr_mul(buf2, tmpc[4][0].im, tmpc[2][0].re, MODE);
  mpfr_fma(buf2, tmpc[4][0].re, tmpc[2][0].im, buf2, MODE);
  mpfr_add(buf1, buf1, S0.re, MODE);
  mpfr_add(buf2, buf2, S0.im, MODE);
  mpfr_mul(out->re, buf1, b, MODE);
  mpfr_mul(out->im, buf2, b, MODE);
  mpfr_clears(w.re, w.im, xi, b, buf1, buf2, buf3, buf4, buf5, (mpfr_ptr) NULL);
}


/*
void div_jacobian(fftwl_complex *in, fftwl_complex *out) {
  // solve a tridiagonal system z_q q_u = b
  // b        -- inverse Fourier coefficients: b = Z_u - 1
  // Note:    -- b must have b[0] = 0, like Z_u - 1
  // Note:       not in-place safe.
  fftwl_complex 	w = cexpl(1.IL*conf.origin_offset);
  long double 		b = 0.5L*(1.L + powl(conf.scaling, 2))/conf.scaling;
  long double 		xi = (1.L - powl(conf.scaling, 2))/(1.L + powl(conf.scaling, 2));

  fft_shift(in);
  tmpc[2][0] = -0.5L*xi*conjl(w);
  for (long int j = 1; j < state.number_modes; j++) {
    tmpc[2][j] = tmpc[2][0]/(1.L - conjl(tmpc[2][0])*tmpc[2][j-1]); 
  }
  tmpc[3][0] = in[0]/b;
  for (long int j = 1; j < state.number_modes; j++) {
    tmpc[3][j] = in[j]/b - conjl(tmpc[2][0])*tmpc[3][j-1];
    tmpc[3][j] = tmpc[3][j]/(1.L - conjl(tmpc[2][0])*tmpc[2][j-1]);
  } 
  out[state.number_modes-1] = tmpc[3][state.number_modes-1];
  for (long int j = state.number_modes-2; j > -1; j--) {
    out[j] = tmpc[3][j]-tmpc[2][j]*out[j+1];
  }
  fft_shift(out);
  fft_shift(in);
}
*/

/*
void linear_solve(fftwl_complex *a, fftwl_complex *b, fftwl_complex *x) {
  // inverts a*x = b to find x by series inversion
  // a,b   -- input inverse Fourier coefficients,
  // x     -- output inverse Fourier coefficients
  // Note: not in-place safe: a and x cannot be the same
  x[0] = b[0]/a[0];
  for (long int j = 1; j < state.number_modes/2; j++) {
    fftwl_complex sum = 0.L;
    for (long int l = j-1; l > -1; l--) {
      sum += a[j-l]*x[l];
    }
    x[j] = (b[j] - sum)/a[0];
  }
  memset(x + state.number_modes/2, 0, (state.number_modes/2)*sizeof(fftwl_complex));
}
*/

void inverse(mpfc_t *a, mpfc_t *x) {
  // inverts a*x = 1 to find x by series inversion
  // a   -- input inverse Fourier coefficients,
  // x   -- output inverse Fourier coefficients
  // Note: not in-place safe: a and x cannot be the same
  mpfc_t bufc;
  mpfr_t buf;
  mpfr_inits2(state.precision, buf, bufc.re, bufc.im, (mpfr_ptr) NULL);
  //x[0] = 1.0L/a[0];
  mpfr_mul(buf, a[0].im, a[0].im, MODE);
  mpfr_fma(buf, a[0].re, a[0].re, buf, MODE);
  mpfr_div(x[0].re, a[0].re, buf, MODE);
  mpfr_div(x[0].im, a[0].im, buf, MODE);
  mpfr_neg(x[0].im, x[0].im, MODE);

  //for (long int j = 1; j < state.number_modes/2; j++) {
  for (long int j = 1; j < 1<<(state.nbits-1); j++) {
    //fftwl_complex sum = 0.L;
    mpfr_set_ui(bufc.re, 0, MODE);
    mpfr_set_ui(bufc.im, 0, MODE);
    for (long int l = j-1; l > -1; l--) {
      //sum += a[j-l]*x[l];
      mpfr_fms(bufc.re, a[j-l].im, x[l].im, bufc.re, MODE);
      mpfr_fms(bufc.re, a[j-l].re, x[l].re, bufc.re, MODE);
 
      mpfr_fma(bufc.im, a[j-l].re, x[l].im, bufc.im, MODE);
      mpfr_fma(bufc.im, a[j-l].im, x[l].re, bufc.im, MODE);
    }
    //x[j] = - sum/a[0];
    mpfr_mul(x[j].re, bufc.re, a[0].re, MODE);
    mpfr_fma(x[j].re, bufc.im, a[0].im, x[j].re, MODE);
    mpfr_div(x[j].re, x[j].re, buf, MODE);
    mpfr_neg(x[j].re, x[j].re, MODE);

    mpfr_mul(x[j].im, bufc.im, a[0].re, MODE);
    mpfr_fms(x[j].im, bufc.re, a[0].im, x[j].im, MODE);
    mpfr_div(x[j].im, x[j].im, buf, MODE);
  }
  //memset(x + state.number_modes/2, 0, (state.number_modes/2)*sizeof(fftwl_complex));
  for (long int j = 1<<(state.nbits-1); j < 1<<(state.nbits); j++) {
    mpfr_set_ui(x[j].re, 0, MODE);
    mpfr_set_ui(x[j].im, 0, MODE);
  }
  mpfr_clears(buf, bufc.re, bufc.im, (mpfr_ptr) NULL);
}

void square_ft(mpfc_t *Z, mpfc_t *x){
  // inverts x^2 = Z to find x by series inversion
  // Z  --  input inverse Fourier coefficients
  // x  --  output inverse Fourier coefficients
  mpfc_t bufc;
  mpfr_t buf;
  mpfr_inits2(state.precision, buf, bufc.re, bufc.im, (mpfr_ptr) NULL);

  //x[0] = csqrtl(Z[0]);
  mpfr_atan2(buf, Z[0].im, Z[0].re, MODE);
  mpfr_div_ui(buf, buf, 2, MODE);  // buf stores phi/2
  mpfr_sin_cos(x[0].im, x[0].re, buf, MODE);
  mpfr_hypot(buf, Z[0].re, Z[0].im, MODE);
  mpfr_sqrt(buf, buf, MODE); // buf stores sqrt(R)
  mpfr_mul(x[0].re, x[0].re, buf, MODE);
  mpfr_mul(x[0].im, x[0].im, buf, MODE);
  //x[1] = 0.5L*Z[1]/x[0];
  //mpfr_hypot(buf, x[0].re, x[0].im, MODE);  	
  mpfr_mul(buf, x[0].re, x[0].re, MODE);  // buf = 2*|x^2|
  mpfr_fma(buf, x[0].im, x[0].im, buf, MODE);
  mpfr_mul_ui(buf, buf, 2, MODE);  	    	//  

  mpfr_mul(x[1].re, Z[1].re, x[0].re, MODE);
  mpfr_fma(x[1].re, Z[1].im, x[0].im, x[1].re, MODE);
  mpfr_div(x[1].re, x[1].re, buf, MODE);

  mpfr_mul(x[1].im, Z[1].re, x[0].im, MODE);
  mpfr_fms(x[1].im, Z[1].im, x[0].re, x[1].im, MODE);
  mpfr_div(x[1].im, x[1].im, buf, MODE);

  for (long int j = 2; j < 1<<(state.nbits-1); j++) {
    //fftwl_complex sum = 0.L;
    mpfr_set_ui(bufc.re, 0, MODE);
    mpfr_set_ui(bufc.im, 0, MODE);
    for (long int m = 0; m < j-1; m++) {
      //sum += x[m+1]*x[j-m-1];
      mpfr_fms(bufc.re, x[m+1].im, x[j-m-1].im, bufc.re, MODE);
      mpfr_fms(bufc.re, x[m+1].re, x[j-m-1].re, bufc.re, MODE);
 
      mpfr_fma(bufc.im, x[m+1].re, x[j-m-1].im, bufc.im, MODE);
      mpfr_fma(bufc.im, x[m+1].im, x[j-m-1].re, bufc.im, MODE);
    }
    //x[j] = 0.5L*(Z[j] - sum)/x[0];
    mpfr_sub(bufc.re, Z[j].re, bufc.re, MODE);
    mpfr_sub(bufc.im, Z[j].im, bufc.im, MODE);
    
    mpfr_mul(x[j].re, bufc.re, x[0].re, MODE);
    mpfr_fma(x[j].re, bufc.im, x[0].im, x[j].re, MODE);
    mpfr_div(x[j].re, x[j].re, buf, MODE);

    mpfr_mul(x[j].im, bufc.re, x[0].im, MODE);
    mpfr_fms(x[j].im, bufc.im, x[0].re, x[j].im, MODE);
    mpfr_div(x[j].im, x[j].im, buf, MODE);
  }
  //memset(x + state.number_modes/2, 0, (state.number_modes/2)*sizeof(fftwl_complex));
  for (long int j = 1<<(state.nbits-1); j < 1<<(state.nbits); j++) {
    mpfr_set_ui(x[j].re, 0, MODE);
    mpfr_set_ui(x[j].im, 0, MODE);
  }
  mpfr_clears(buf, bufc.re, bufc.im, (mpfr_ptr) NULL);
}
