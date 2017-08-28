#include "header.h"

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

void compute_zero_mode_complex(fftwl_complex *in, fftwl_complex S0, fftwl_complex *out) {
  fftwl_complex 	w = cexpl(1.IL*conf.origin_offset);
  long double 		b = 0.5L*(1.L + powl(conf.scaling, 2))/conf.scaling;
  long double 		xi = (1.L - powl(conf.scaling, 2))/(1.L + powl(conf.scaling, 2));
  tmpc[2][0] = -0.5L*xi*conjl(w);
  for (long int j = 1; j < state.number_modes/2-1; j++) {
    tmpc[2][j] = tmpc[2][0]/(1.L - conjl(tmpc[2][0])*tmpc[2][j-1]); 
  }
  tmpc[3][0] = in[1]/b - 2.L*S0*conjl(tmpc[2][0]);
  for (long int j = 1; j < state.number_modes/2-1; j++) {
    tmpc[3][j] = in[j+1]/b - conjl(tmpc[2][0])*tmpc[3][j-1];
    tmpc[3][j] = tmpc[3][j]/(1.L - conjl(tmpc[2][0])*tmpc[2][j-1]);
  } 
  tmpc[4][state.number_modes/2-2] = tmpc[3][state.number_modes/2-2];
  for (long int j = state.number_modes/2-3; j > -1; j--) {
    tmpc[4][j] = tmpc[3][j]-tmpc[2][j]*tmpc[4][j+1];
  }
  *out = b*(tmpc[4][0]*tmpc[2][0] + 1.L*S0);
}

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

void inverse(fftwl_complex *a, fftwl_complex *x) {
  // inverts a*x = 1 to find x by series inversion
  // a   -- input inverse Fourier coefficients,
  // x   -- output inverse Fourier coefficients
  // Note: not in-place safe: a and x cannot be the same
  x[0] = 1.0L/a[0];
  for (long int j = 1; j < state.number_modes/2; j++) {
    fftwl_complex sum = 0.L;
    for (long int l = j-1; l > -1; l--) {
      sum += a[j-l]*x[l];
    }
    x[j] = - sum/a[0];
  }
  memset(x + state.number_modes/2, 0, (state.number_modes/2)*sizeof(fftwl_complex));
}

void square_ft(fftwl_complex *Z, fftwl_complex *x){
  // inverts x^2 = Z to find x by series inversion
  // Z  --  input inverse Fourier coefficients
  // x  --  output inverse Fourier coefficients
  x[0] = csqrtl(Z[0]);
  x[1] = 0.5L*Z[1]/x[0];
  for (long int j = 2; j < state.number_modes/2; j++) {
    fftwl_complex sum = 0.L;
    for (long int m = 0; m < j-1; m++) {
      sum += x[m+1]*x[j-m-1];
    }
    x[j] = 0.5L*(Z[j] - sum)/x[0];
  }
  memset(x + state.number_modes/2, 0, (state.number_modes/2)*sizeof(fftwl_complex));
}
