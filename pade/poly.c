#include "header.h"

static mpfc_t *tmp;

void init_poly(unsigned int nD) {
  /*
  tmp = malloc((2*nD+1)*sizeof(mpfc_t));
  for (int j = 0; j < 2*nD + 1; j++) {
    mpfr_inits2(precision, tmp[j].re, tmp[j].im, (mpfr_ptr) NULL);
  }
  */
  tmp = malloc((pade_data.Npt-1)*sizeof(mpfc_t));
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_inits2(precision, tmp[j].re, tmp[j].im, (mpfr_ptr) NULL);
  }
}

void free_poly(unsigned int nD) {
  //for (int j = 0; j < 2*nD + 1; j++) {
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_clears(tmp[j].re, tmp[j].im, (mpfr_ptr) NULL);
  }
  free(tmp);
}

void evaluate_poly_array(unsigned long nD) {
  mpfc_t s;
  mpfc_t t;
  mpfr_inits2(precision, s.re, s.im, t.re, t.im, (mpfr_ptr) NULL);
  init_poly(nD);
  mpfr_set_ui(s.im, 0, MODE);
  //unsigned long N = state.number_modes;
  //long double	s, overN = 1.L/N;
  //printf("Here\n");
  
  for (unsigned long j = 0; j < pade_data.Npt-1; j++) {
    //printf("Here %lu\n", j);
    mpfr_mul_si(s.re, Ovn, j+1, MODE);
    mpfr_add_d(s.re, s.re, -0.5, MODE);
    mpfr_mul(s.re, s.re, Pie, MODE);
    mpfr_tan(s.re, s.re, MODE);
    //s = tanl(0.5L*PI*(2.L*(j+1.L)*overN - 1.L));
    // --------   step 0
    //A[0][j] = 0.L;
    mpfr_set_ui(A[0][j].re, 0, MODE);
    mpfr_set_ui(A[0][j].im, 0, MODE);
    //B[0][j] = 1.L;
    mpfr_set_ui(B[0][j].re, 1, MODE);
    mpfr_set_ui(B[0][j].im, 0, MODE);
    // --------   step 1
    //A[1][j] = -1.L;
    mpfr_set_si(A[1][j].re, -1, MODE);
    mpfr_set_si(A[1][j].im,  0, MODE);
    //B[1][j] = -Cees[0][1];
    mpfr_neg(B[1][j].re, Cees[0][1].re, MODE);
    mpfr_neg(B[1][j].im, Cees[0][1].im, MODE);
    // --------   step 2
    //A[2][j] = s*A[0][j] - Cees[0][2]*A[1][j] - Cees[1][2]*A[0][j];
    mpfc_mul(&t, &Cees[1][2], &A[0][j], MODE);
    mpfc_fma(&t, &Cees[0][2], &A[1][j], &t, MODE);
    mpfc_fms(&A[2][j],    &s, &A[0][j], &t, MODE);
    //B[2][j] = s*B[0][j] - Cees[0][2]*B[1][j] - Cees[1][2]*B[0][j];
    mpfc_mul(&t, &Cees[1][2], &B[0][j], MODE);
    mpfc_fma(&t, &Cees[0][2], &B[1][j], &t, MODE);
    mpfc_fms(&B[2][j],    &s, &B[0][j], &t, MODE);
    // --------   step 3
    //A[3][j] = s*A[1][j] - Cees[0][3]*A[2][j] - Cees[1][3]*A[1][j] - Cees[2][3]*A[0][j];
    mpfc_mul(&t, &Cees[2][3], &A[0][j], MODE);
    mpfc_fma(&t, &Cees[1][3], &A[1][j], &t, MODE);
    mpfc_fma(&t, &Cees[0][3], &A[2][j], &t, MODE);
    mpfc_fms(&A[3][j],    &s, &A[1][j], &t, MODE);
    //B[3][j] = s*B[1][j] - Cees[0][3]*B[2][j] - Cees[1][3]*B[1][j] - Cees[2][3]*B[0][j];
    mpfc_mul(&t, &Cees[2][3], &B[0][j], MODE);
    mpfc_fma(&t, &Cees[1][3], &B[1][j], &t, MODE);
    mpfc_fma(&t, &Cees[0][3], &B[2][j], &t, MODE);
    mpfc_fms(&B[3][j],    &s, &B[1][j], &t, MODE);
    // --------   step 4+
    //printf("Here %lu\n", j);
    for (unsigned long k = 4; k < 2*nD + 1; k++) {
      //printf("\tinner %lu\n", k);
      //tmp[j] = s*A[2][j] - Cees[0][k]*A[3][j] - Cees[1][k]*A[2][j] - Cees[2][k]*A[1][j] - Cees[3][k]*A[0][j];	
      mpfc_mul(&t, &Cees[3][k], &A[0][j], MODE);
      //printf("\tinner %lu\n", k);
      mpfc_fma(&t, &Cees[2][k], &A[1][j], &t, MODE);
      mpfc_fma(&t, &Cees[1][k], &A[2][j], &t ,MODE);
      mpfc_fma(&t, &Cees[0][k], &A[3][j], &t, MODE);
      //printf("\tinner %lu\n", k);
      mpfc_fms(&tmp[j],     &s, &A[2][j], &t, MODE);
      //printf("\tinner %lu\n", k);
      /*A[0][j] = A[1][j];
      A[1][j] = A[2][j];
      A[2][j] = A[3][j];
      A[3][j] = tmp[j];*/
      mpfc_set(&A[0][j], &A[1][j], MODE);
      mpfc_set(&A[1][j], &A[2][j], MODE);
      mpfc_set(&A[2][j], &A[3][j], MODE);
      mpfc_set(&A[3][j], &tmp[j], MODE);

      //tmp[j] = s*B[2][j] - Cees[0][k]*B[3][j] - Cees[1][k]*B[2][j] - Cees[2][k]*B[1][j] - Cees[3][k]*B[0][j];
      mpfc_mul(&t, &Cees[3][k], &B[0][j], MODE);
      mpfc_fma(&t, &Cees[2][k], &B[1][j], &t, MODE);
      mpfc_fma(&t, &Cees[1][k], &B[2][j], &t ,MODE);
      mpfc_fma(&t, &Cees[0][k], &B[3][j], &t, MODE);
      mpfc_fms(&tmp[j],     &s, &B[2][j], &t, MODE);
      /*B[0][j] = B[1][j];
      B[1][j] = B[2][j];
      B[2][j] = B[3][j];
      B[3][j] = tmp[j];*/
      mpfc_set(&B[0][j], &B[1][j], MODE);
      mpfc_set(&B[1][j], &B[2][j], MODE);
      mpfc_set(&B[2][j], &B[3][j], MODE);
      mpfc_set(&B[3][j], &tmp[j], MODE);
    }  
    //printf("Here %lu\n", j);
  }
  //memcpy(arrayQ, B[3], (N-1)*sizeof(fftwl_complex));
  //memcpy(arrayP, A[3], (N-1)*sizeof(fftwl_complex));
  for (unsigned long j = 0; j < pade_data.Npt-1; j++) {
    mpfc_set(&arrQ[j], &B[3][j], MODE);
    mpfc_set(&arrP[j], &A[3][j], MODE);
  }
  mpfr_clears(s.re, s.im, t.re, t.im, (mpfr_ptr) NULL);
  free_poly(nD);
}



void init_poly_steps(mpfc_t *P, mpfc_t *Q, mpfc_t *dQ, unsigned long nD) {
  for (int j = 0; j < 4; j++) {
    mpfr_inits2(precision,  P[j].re,  P[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(precision,  Q[j].re,  Q[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(precision, dQ[j].re, dQ[j].im, (mpfr_ptr) NULL);
  }
  mpfc_set_si(&Q[0],  1, 0, MODE);
  mpfc_set_si(&P[0],  0, 0, MODE);
  mpfc_set_si(&dQ[0], 0, 0, MODE);

  mpfc_neg(&Q[1], &Cees[0][1], MODE);   
  mpfc_set_si(&P[1],  -1, 0, MODE);
  mpfc_set_si(&dQ[1],  0, 0, MODE);
}

void evaluate_poly(mpfc_t *in, unsigned long nD, mpfc_t *outQ, mpfc_t *outQp, mpfc_t *outP){
  mpfc_t  	tpc, tpb, tpa;
  mpfc_t	P[4], Q[4], dQ[4];
  mpfc_t	s, buf;
 
  mpfr_inits2(precision, s.re, s.im, (mpfr_ptr) NULL);
  mpfr_inits2(precision, buf.re, buf.im, (mpfr_ptr) NULL);
  mpfr_inits2(precision, tpa.re, tpa.im, (mpfr_ptr) NULL);
  mpfr_inits2(precision, tpb.re, tpb.im, (mpfr_ptr) NULL);
  mpfr_inits2(precision, tpc.re, tpc.im, (mpfr_ptr) NULL);
  init_poly_steps(P, Q, dQ, nD);
  //s = *in;
  mpfr_set(s.re, in->re, MODE);
  mpfr_set(s.im, in->im, MODE);
  
  /*
  P[0] = 0.L;
  dQ[0] = 0.L;
  Q[0] = 1.L;
  P[1] = -1.L;   -- Set in init_poly_steps
  dQ[1] = 0.L;
  Q[1] = -Cees[0][1];
  */

  /*
  P[2]  = s*P[0] - Cees[0][2]*P[1] - Cees[1][2]*P[0];
  Q[2]  = s*Q[0] - Cees[0][2]*Q[1] - Cees[1][2]*Q[0];
  dQ[2] = Q[0] + s*dQ[0] - Cees[0][2]*dQ[1] - Cees[1][2]*dQ[0];
  */
  mpfc_mul(&buf, &Cees[1][2], &P[0], MODE); 
  mpfc_fma(&buf, &Cees[0][2], &P[1], &buf, MODE);
  mpfc_fms(&P[2],         &s, &P[0], &buf, MODE);

  mpfc_mul(&buf, &Cees[1][2], &Q[0], MODE);
  mpfc_fma(&buf, &Cees[0][2], &Q[1], &buf, MODE);
  mpfc_fms(&Q[2],         &s, &Q[0], &buf, MODE);

  mpfc_fms(&buf, &Cees[1][2], &dQ[0], &Q[0], MODE);
  mpfc_fma(&buf, &Cees[0][2], &dQ[1], &buf, MODE);
  mpfc_fms(&dQ[2],        &s, &dQ[0], &buf, MODE);

  /*
  P[3]  = s*P[1] - Cees[0][3]*P[2] - Cees[1][3]*P[1] - Cees[2][3]*P[0];
  Q[3]  = s*Q[1] - Cees[0][3]*Q[2] - Cees[1][3]*Q[1] - Cees[2][3]*Q[0];
  dQ[3] = Q[1] + s*dQ[1] - Cees[0][3]*dQ[2] - Cees[1][3]*dQ[1] - Cees[2][3]*dQ[0];
  */
  mpfc_mul(&buf, &Cees[2][3], &P[0], MODE);
  mpfc_fma(&buf, &Cees[1][3], &P[1], &buf, MODE);
  mpfc_fma(&buf, &Cees[0][3], &P[2], &buf, MODE);
  mpfc_fms(&P[3],         &s, &P[1], &buf, MODE);
  
  mpfc_mul(&buf, &Cees[2][3], &Q[0], MODE);
  mpfc_fma(&buf, &Cees[1][3], &Q[1], &buf, MODE);
  mpfc_fma(&buf, &Cees[0][3], &Q[2], &buf, MODE);
  mpfc_fms(&Q[3],         &s, &Q[1], &buf, MODE);
  
  mpfc_fms(&buf, &Cees[2][3], &dQ[0], &Q[1], MODE);
  mpfc_fma(&buf, &Cees[1][3], &dQ[1], &buf, MODE);
  mpfc_fma(&buf, &Cees[0][3], &dQ[2], &buf, MODE);
  mpfc_fms(&Q[3],         &s, &dQ[1], &buf, MODE);

  for (unsigned long k = 4; k < 2*nD + 1; k++) {
    /*
    tpc = s*P[2] - Cees[0][k]*P[3] - Cees[1][k]*P[2] - Cees[2][k]*P[1] - Cees[3][k]*P[0];
    tpb = s*Q[2] - Cees[0][k]*Q[3] - Cees[1][k]*Q[2] - Cees[2][k]*Q[1] - Cees[3][k]*Q[0];
    tpa = Q[2] + s*dQ[2] - Cees[0][k]*dQ[3] - Cees[1][k]*dQ[2] - Cees[2][k]*dQ[1] - Cees[3][k]*dQ[0];
    */
    mpfc_mul(&buf, &Cees[3][k], &P[0], MODE);
    mpfc_fma(&buf, &Cees[2][k], &P[1], &buf, MODE);
    mpfc_fma(&buf, &Cees[1][k], &P[2], &buf, MODE);
    mpfc_fma(&buf, &Cees[0][k], &P[3], &buf, MODE);
    mpfc_fms(&tpc,          &s, &P[2], &buf, MODE);
        
    mpfc_mul(&buf, &Cees[3][k], &Q[0], MODE);
    mpfc_fma(&buf, &Cees[2][k], &Q[1], &buf, MODE);
    mpfc_fma(&buf, &Cees[1][k], &Q[2], &buf, MODE);
    mpfc_fma(&buf, &Cees[0][k], &Q[3], &buf, MODE);
    mpfc_fms(&tpb,          &s, &Q[2], &buf, MODE);

    mpfc_fms(&buf, &Cees[3][k], &dQ[0], &Q[2], MODE);
    mpfc_fma(&buf, &Cees[2][k], &dQ[1], &buf, MODE);
    mpfc_fma(&buf, &Cees[1][k], &dQ[2], &buf, MODE);
    mpfc_fma(&buf, &Cees[0][k], &dQ[3], &buf, MODE);
    mpfc_fms(&tpa,          &s, &dQ[2], &buf, MODE);

    /*
    P[0] = P[1];
    P[1] = P[2];
    P[2] = P[3];
    P[3] = tpc;
    */
    mpfc_set(&P[0], &P[1], MODE);
    mpfc_set(&P[1], &P[2], MODE);
    mpfc_set(&P[2], &P[3], MODE);
    mpfc_set(&P[3], &tpc, MODE);
    /*
    Q[0] = Q[1];
    Q[1] = Q[2];
    Q[2] = Q[3];
    Q[3] = tpb;
    */
    mpfc_set(&Q[0], &Q[1], MODE);
    mpfc_set(&Q[1], &Q[2], MODE);
    mpfc_set(&Q[2], &Q[3], MODE);
    mpfc_set(&Q[3], &tpb, MODE);
    /*
    dQ[0] = dQ[1];
    dQ[1] = dQ[2];
    dQ[2] = dQ[3];
    dQ[3] = tpa; 
    */
    mpfc_set(&dQ[0], &dQ[1], MODE);
    mpfc_set(&dQ[1], &dQ[2], MODE);
    mpfc_set(&dQ[2], &dQ[3], MODE);
    mpfc_set(&dQ[3], &tpa, MODE);
  }
  /*
  *outP  = P[3];  
  *outQ  = Q[3];
  *outQp = dQ[3];
  */
  mpfc_set(outP,  &P[3], MODE);
  mpfc_set(outQ,  &Q[3], MODE);
  mpfc_set(outQp, &dQ[3], MODE);
}


void poly_val_array(mpfc_t *in, unsigned long nD, mpfc_t *outQ, mpfc_t *outQp, mpfc_t *outP){
  /*
  fftwl_complex  	*tpc, *tpb, *tpa;
  fftwl_complex		*P[4], *Q[4], *dQ[4];
  fftwl_complex		*s;
  */
  mpfc_t *tpc, *tpb, *tpa;
  mpfc_t *P[4], *Q[4], *dQ[4];
  mpfc_t *s, *buf;

  s   = malloc(nD*sizeof(mpfc_t));
  buf = malloc(nD*sizeof(mpfc_t));
  tpa = malloc(nD*sizeof(mpfc_t));
  tpb = malloc(nD*sizeof(mpfc_t));
  tpc = malloc(nD*sizeof(mpfc_t));
  for (unsigned int j = 0; j < nD; j++) {
    mpfr_inits2(precision, buf[j].re, buf[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(precision, tpa[j].re, tpa[j].im, (mpfr_ptr) NULL); 
    mpfr_inits2(precision, tpb[j].re, tpb[j].im, (mpfr_ptr) NULL); 
    mpfr_inits2(precision, tpc[j].re, tpc[j].im, (mpfr_ptr) NULL); 
    mpfr_inits2(precision, s[j].re, s[j].im, (mpfr_ptr) NULL); 
  }
  for (int k = 0; k < 4; k++) {
    P[k]  = malloc(nD*sizeof(mpfc_t));
    Q[k]  = malloc(nD*sizeof(mpfc_t));
    dQ[k] = malloc(nD*sizeof(mpfc_t));
    for (unsigned int j = 0; j < nD; j++) {
      mpfr_inits2(precision,  P[k][j].re,  P[k][j].im, (mpfr_ptr) NULL); 
      mpfr_inits2(precision,  Q[k][j].re,  Q[k][j].im, (mpfr_ptr) NULL); 
      mpfr_inits2(precision, dQ[k][j].re, dQ[k][j].im, (mpfr_ptr) NULL); 
    }
  }

  for (int j = 0; j < nD; j++) {
    //s[j] = in[j];
    mpfr_set(s[j].re, in[j].re, MODE);
    mpfr_set(s[j].im, in[j].im, MODE);
   
    /*
    P[0][j] = 0.L;
    Q[0][j] = 1.L;
    dQ[0][j] = 0.L;
    */
    mpfc_set_si( &P[0][j], 0, 0, MODE);
    mpfc_set_si( &Q[0][j], 1, 0, MODE);
    mpfc_set_si(&dQ[0][j], 0, 0, MODE);

    /*
    P[1][j] = -1.L;
    Q[1][j] = -Cees[0][1];
    dQ[1][j] = 0.L;
    */
    mpfc_neg(&Q[1][j], &Cees[0][1], MODE);   
    mpfc_set_si(&P[1][j],  -1, 0, MODE);
    mpfc_set_si(&dQ[1][j],  0, 0, MODE);
    /*
    P[2][j]  = s[j]*P[0][j] - Cees[0][2]*P[1][j] - Cees[1][2]*P[0][j];
    Q[2][j]  = s[j]*Q[0][j] - Cees[0][2]*Q[1][j] - Cees[1][2]*Q[0][j];
    dQ[2][j] = Q[0][j] + s[j]*dQ[0][j] - Cees[0][2]*dQ[1][j] - Cees[1][2]*dQ[0][j];
    */

    mpfc_mul(&buf[j], &Cees[1][2], &P[0][j], MODE); 
    mpfc_fma(&buf[j], &Cees[0][2], &P[1][j], &buf[j], MODE);
    mpfc_fms(&P[2][j],      &s[j], &P[0][j], &buf[j], MODE);

    mpfc_mul(&buf[j], &Cees[1][2], &Q[0][j], MODE);
    mpfc_fma(&buf[j], &Cees[0][2], &Q[1][j], &buf[j], MODE);
    mpfc_fms(&Q[2][j],      &s[j], &Q[0][j], &buf[j], MODE);

    mpfc_fms(&buf[j], &Cees[1][2], &dQ[0][j], &Q[0][j], MODE);
    mpfc_fma(&buf[j], &Cees[0][2], &dQ[1][j], &buf[j], MODE);
    mpfc_fms(&dQ[2][j],     &s[j], &dQ[0][j], &buf[j], MODE);

    /*
    P[3][j]  = s[j]*P[1][j] - Cees[0][3]*P[2][j] - Cees[1][3]*P[1][j] - Cees[2][3]*P[0][j];
    Q[3][j]  = s[j]*Q[1][j] - Cees[0][3]*Q[2][j] - Cees[1][3]*Q[1][j] - Cees[2][3]*Q[0][j];
    dQ[3][j] = Q[1][j] + s[j]*dQ[1][j] - Cees[0][3]*dQ[2][j] - Cees[1][3]*dQ[1][j] - Cees[2][3]*dQ[0][j];
    */

    mpfc_mul(&buf[j], &Cees[2][3], &P[0][j], MODE);
    mpfc_fma(&buf[j], &Cees[1][3], &P[1][j], &buf[j], MODE);
    mpfc_fma(&buf[j], &Cees[0][3], &P[2][j], &buf[j], MODE);
    mpfc_fms(&P[3][j],      &s[j], &P[1][j], &buf[j], MODE);
  
    mpfc_mul(&buf[j], &Cees[2][3], &Q[0][j], MODE);
    mpfc_fma(&buf[j], &Cees[1][3], &Q[1][j], &buf[j], MODE);
    mpfc_fma(&buf[j], &Cees[0][3], &Q[2][j], &buf[j], MODE);
    mpfc_fms(&Q[3][j],      &s[j], &Q[1][j], &buf[j], MODE);
  
    mpfc_fms(&buf[j], &Cees[2][3], &dQ[0][j], &Q[1][j], MODE);
    mpfc_fma(&buf[j], &Cees[1][3], &dQ[1][j], &buf[j], MODE);
    mpfc_fma(&buf[j], &Cees[0][3], &dQ[2][j], &buf[j], MODE);
    mpfc_fms(&Q[3][j],      &s[j], &dQ[1][j], &buf[j], MODE);

    for (unsigned long k = 4; k < 2*nD + 1; k++) {
      /*
      tpc[j] = s[j]*P[2][j] - Cees[0][k]*P[3][j] - Cees[1][k]*P[2][j] - Cees[2][k]*P[1][j] - Cees[3][k]*P[0][j];
      tpb[j] = s[j]*Q[2][j] - Cees[0][k]*Q[3][j] - Cees[1][k]*Q[2][j] - Cees[2][k]*Q[1][j] - Cees[3][k]*Q[0][j];
      tpa[j] = Q[2][j] + s[j]*dQ[2][j] - Cees[0][k]*dQ[3][j] - Cees[1][k]*dQ[2][j] - Cees[2][k]*dQ[1][j] - Cees[3][k]*dQ[0][j];
      */
      mpfc_mul(&buf[j], &Cees[3][k], &P[0][j], MODE);
      mpfc_fma(&buf[j], &Cees[2][k], &P[1][j], &buf[j], MODE);
      mpfc_fma(&buf[j], &Cees[1][k], &P[2][j], &buf[j], MODE);
      mpfc_fma(&buf[j], &Cees[0][k], &P[3][j], &buf[j], MODE);
      mpfc_fms(&tpc[j],       &s[j], &P[2][j], &buf[j], MODE);
        
      mpfc_mul(&buf[j], &Cees[3][k], &Q[0][j], MODE);
      mpfc_fma(&buf[j], &Cees[2][k], &Q[1][j], &buf[j], MODE);
      mpfc_fma(&buf[j], &Cees[1][k], &Q[2][j], &buf[j], MODE);
      mpfc_fma(&buf[j], &Cees[0][k], &Q[3][j], &buf[j], MODE);
      mpfc_fms(&tpb[j],       &s[j], &Q[2][j], &buf[j], MODE);

      mpfc_fms(&buf[j], &Cees[3][k], &dQ[0][j], &Q[2][j], MODE);
      mpfc_fma(&buf[j], &Cees[2][k], &dQ[1][j], &buf[j], MODE);
      mpfc_fma(&buf[j], &Cees[1][k], &dQ[2][j], &buf[j], MODE);
      mpfc_fma(&buf[j], &Cees[0][k], &dQ[3][j], &buf[j], MODE);
      mpfc_fms(&tpa[j],       &s[j], &dQ[2][j], &buf[j], MODE);

      /*
      P[0][j] = P[1][j];
      P[1][j] = P[2][j];
      P[2][j] = P[3][j];
      P[3][j] = tpc[j];
      */
      mpfc_set(&P[0][j], &P[1][j], MODE);
      mpfc_set(&P[1][j], &P[2][j], MODE);
      mpfc_set(&P[2][j], &P[3][j], MODE);
      mpfc_set(&P[3][j], &tpc[j], MODE);

      /*
      Q[0][j] = Q[1][j];
      Q[1][j] = Q[2][j];
      Q[2][j] = Q[3][j];
      Q[3][j] = tpb[j];
      */
      mpfc_set(&Q[0][j], &Q[1][j], MODE);
      mpfc_set(&Q[1][j], &Q[2][j], MODE);
      mpfc_set(&Q[2][j], &Q[3][j], MODE);
      mpfc_set(&Q[3][j], &tpb[j], MODE);

      /*
      dQ[0][j] = dQ[1][j];
      dQ[1][j] = dQ[2][j];
      dQ[2][j] = dQ[3][j];
      dQ[3][j] = tpa[j]; 
      */
      mpfc_set(&dQ[0][j], &dQ[1][j], MODE);
      mpfc_set(&dQ[1][j], &dQ[2][j], MODE);
      mpfc_set(&dQ[2][j], &dQ[3][j], MODE);
      mpfc_set(&dQ[3][j], &tpa[j], MODE);
    } 
    /*
    outP[j]  = P[3][j];  
    outQ[j]  = Q[3][j];
    outQp[j] = dQ[3][j];
    */
    mpfc_set(&outP[j],  &P[3][j], MODE);
    mpfc_set(&outQ[j],  &Q[3][j], MODE);
    mpfc_set(&outQp[j], &dQ[3][j], MODE);
  }
  /*
  fftwl_free(tpa);
  fftwl_free(tpb);
  fftwl_free(tpc);
  for (int j = 3; j > -1; j--) {
    fftwl_free(P[j]);
    fftwl_free(Q[j]);
    fftwl_free(dQ[j]);
  }
  fftwl_free(s);
  */
  for (unsigned int j = 0; j < nD; j++) {
    mpfr_clears(buf[j].re, buf[j].im, (mpfr_ptr) NULL);
    mpfr_clears(tpa[j].re, tpa[j].im, (mpfr_ptr) NULL); 
    mpfr_clears(tpb[j].re, tpb[j].im, (mpfr_ptr) NULL); 
    mpfr_clears(tpc[j].re, tpc[j].im, (mpfr_ptr) NULL); 
    mpfr_clears(  s[j].re,   s[j].im, (mpfr_ptr) NULL); 
  }
  for (int k = 3; k > -1; k--) {
    for (unsigned int j = 0; j < nD; j++) {
      mpfr_clears( P[k][j].re,  P[k][j].im, (mpfr_ptr) NULL); 
      mpfr_clears( Q[k][j].re,  Q[k][j].im, (mpfr_ptr) NULL); 
      mpfr_clears(dQ[k][j].re, dQ[k][j].im, (mpfr_ptr) NULL); 
    }
    free(P[k]); 
    free(Q[k]); 
    free(dQ[k]); 
  }
}
