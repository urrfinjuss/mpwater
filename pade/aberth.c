#include "header.h"

static mpfc_t *rts, *res;
static mpfc_t *t1, *t2;
static mpfc_t *wk, *ivc;
static mpfc_t tmp_r;

void init_aberth(unsigned int nD) {
  // init arrays
  rts = malloc(nD*sizeof(mpfc_t));
  res = malloc(nD*sizeof(mpfc_t));
  ivc = malloc(nD*sizeof(mpfc_t));
  wk = malloc(nD*sizeof(mpfc_t));
  t1 = malloc(nD*sizeof(mpfc_t));
  t2 = malloc(nD*sizeof(mpfc_t));
  for (unsigned int j = 0; j < nD; j++) {
    mpfr_inits2(precision, rts[j].re, rts[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(precision, res[j].re, res[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(precision, ivc[j].re, ivc[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(precision, t1[j].re, t1[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(precision, t2[j].re, t2[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(precision, wk[j].re, wk[j].im, (mpfr_ptr) NULL);
  }
  mpfr_inits2(precision, tmp_r.re, tmp_r.im, (mpfr_ptr) NULL);
  // set initial guess
  mpfr_t phs;
  mpfr_init2(phs, precision);
  
  for (unsigned j = 0; j < nD; j++) {
    mpfr_mul_ui(phs, Pie, 2*j + 1, MODE);
    mpfr_div_ui(phs, phs, nD, MODE);
    mpfr_sin_cos(rts[j].re, rts[j].im, phs, MODE);
    mpfr_div_si(rts[j].re, rts[j].re,  -4, MODE);
    mpfr_div_si(rts[j].im, rts[j].im,   4, MODE);
    mpfr_add_d(rts[j].im, rts[j].im, 1.25, MODE);

    mpfc_set_si(&wk[j], 0, 0, MODE); 
  }
  mpfr_clear(phs);
}

void free_aberth(unsigned int nD) {
  for (unsigned int j = 0; j < nD; j++) {
    mpfr_clears(rts[j].re, rts[j].im, (mpfr_ptr) NULL);
    mpfr_clears(res[j].re, res[j].im, (mpfr_ptr) NULL);
    mpfr_clears(ivc[j].re, ivc[j].im, (mpfr_ptr) NULL);
    mpfr_clears(t1[j].re, t1[j].im, (mpfr_ptr) NULL);
    mpfr_clears(t2[j].re, t2[j].im, (mpfr_ptr) NULL);
    mpfr_clears(wk[j].re, wk[j].im, (mpfr_ptr) NULL);
  }
  mpfr_clears(tmp_r.re, tmp_r.im, (mpfr_ptr) NULL);
  free(rts);
  free(res);
  free(ivc);
  free(wk);
  free(t1);
  free(t2);
}


void aberth_iter(unsigned int nD, char *str) {
  /*
  fftwl_complex *rts = fftwl_malloc(nD*sizeof(fftwl_complex));
  fftwl_complex *res = fftwl_malloc(nD*sizeof(fftwl_complex));
  fftwl_complex *t1 = fftwl_malloc(nD*sizeof(fftwl_complex));
  fftwl_complex *t2 = fftwl_malloc(nD*sizeof(fftwl_complex));
  fftwl_complex *wk = fftwl_malloc(nD*sizeof(fftwl_complex));
  fftwl_complex *ivc = fftwl_malloc(nD*sizeof(fftwl_complex));
  fftwl_complex tmp_r;
  */
  init_aberth(nD);
  //long double nrm = 1.L;;
  long int    counter = 0, k = 0;
  mpfc_t      buf;
  mpfr_t      nrm;
  mpfr_inits2(precision, nrm, buf.re, buf.im, (mpfr_ptr) NULL);
  mpfr_set_ui(nrm, 1, MODE);

  FILE *fh;
  /*
  for (int j = 0; j < nD; j++) {
    rts[j] = 1.25IL + 0.25IL*cexpl(2.IL*PI*(j+0.5L)/nD);     <--- part of init_aberth() now
    wk[j] = 0.L;
  }
  */
  //for (int k = 0; k < 40; k++) {
  //while (nrm > 5.E-28L) {
  while (mpfr_cmp_ld(nrm, 5.0E-28L) > 0) {
    poly_val_array(rts, nD, t1, t2, res);   // <-- to be rewritten
    //nrm = 0.L;
    mpfr_set_si(nrm, 0, MODE);
    for (int n = 0; n < nD; n++) {
      //nrm += cabsl(t1[n])*cabsl(t1[n]);
      mpfr_fma(nrm, t1[n].re, t1[n].re, nrm, MODE);
      mpfr_fma(nrm, t1[n].im, t1[n].im, nrm, MODE);
    }
    //nrm = sqrtl(nrm);
    mpfr_sqrt(nrm, nrm, MODE);
    //printf("Newton error: %.Le\n", nrm);
    for (int j = 0; j < nD; j++) {
      // tmp_r = t1[j]/t2[j]; 
      // ivc[j] = 0.L;
      mpfc_div(&tmp_r, &t1[j], &t2[j], MODE);
      mpfc_set_si(&ivc[j], 0, 0, MODE);
      for (int l = 0; l < nD; l++) {
        if (j != l) {
          //ivc[j] += 1.L/(rts[j] - rts[l]);
          mpfc_sub(&buf, &rts[j], &rts[l], MODE);
          mpfc_si_div(&buf, 1, &buf, MODE);
          mpfc_add(&ivc[j], &ivc[j], &buf, MODE);
        }
      }
      //wk[j] = - 1.L*tmp_r/(1.L - tmp_r*ivc[j]);
      //rts[j] += wk[j];
      mpfc_mul(&buf, &tmp_r, &ivc[j], MODE);
      mpfr_ui_sub(buf.re, 1, buf.re, MODE);
      mpfr_neg(buf.im, buf.im, MODE);
      mpfc_div(&buf, &tmp_r, &buf, MODE);
      mpfc_neg(&wk[j], &buf, MODE);
      mpfc_add(&rts[j], &rts[j], &wk[j], MODE); 
    }
    if ((k % 2) == 0) {
      /*sprintf(str, "./roots/roots_%03ld.txt", counter);
      fh = fopen(str,"w");
      for (int j = 0; j < nD; j++) {
        fprintf(fh, "%3d\t%.12LE\t%.12LE\n", j, creall(2.L*catanl(rts[j])), cimagl(2.L*catanl(rts[j])));
      }
      fclose(fh);*/
      counter++;
      mpfr_printf("Aberth iteration %3ld: %.12Re\n", k, nrm);
    }
    k++;
    if ( k == 120) break;
  }
  poly_val_array(rts, nD, t1, t2, res);  // <-- to be rewritten
  //nrm = 0.L;
  mpfr_set_si(nrm, 0, MODE); 
  for (int j = 0; j < nD; j++) {
    //res[j] = res[j]/t2[j];
    //nrm += cabsl(t1[j])*cabsl(t1[j]);
    mpfc_div(&res[j], &res[j], &t2[j], MODE);
    mpfr_fma(nrm, t1[j].re, t1[j].re, nrm, MODE);
    mpfr_fma(nrm, t1[j].im, t1[j].im, nrm, MODE);
  }
  sort_imag(rts, res, nD);            //  <-- to be rewritten
  //sprintf(str, "./roots/roots.txt");
  fh = fopen(str,"w");
  fprintf(fh, "# 1. pole # 2.-3. z_k 4.-5. gamma_k\n");
  mpfr_fprintf(fh, "# Pade Rel. Error = %.12Re\n\n", pade_data.l2_rel_err);
  //long double r0, r1, i0, i1;
  //fftwl_complex qK;
  mpfr_t r0, r1, i0, i1;
  mpfc_t qK;
  mpfr_inits2(precision, r0, r1, i0, i1, (mpfr_ptr) NULL); 
  mpfr_inits2(precision, qK.re, qK.im, (mpfr_ptr) NULL);
  for (int j = 0; j < nD; j++) {
    /*
    r0 = creall(2.L*catanl(conf.scaling*rts[j]));
    i0 = cimagl(2.L*catanl(conf.scaling*rts[j]));
    r1 = creall(conf.scaling*res[j]);
    i1 = cimagl(conf.scaling*res[j]);
    qK = (r1 + 1.IL*i1)*cpowl(ccoshl(0.5L*(r0 + 1.IL*i0)),-2); // changed from 2 to -2
    */
    mpfr_set(r0, rts[j].re, MODE);
    mpfr_set(i0, rts[j].im, MODE);
    mpfr_set(r1, res[j].re, MODE);
    mpfr_set(i1, res[j].im, MODE);
    //fprintf(fh, "%3d\t%19.12Re\t%19.12Re\t%19.12Re\t%19.12Re\n", j, r0, i0, creall(qK), cimagl(qK));
    mpfr_fprintf(fh, "%3d\t%19.12Re\t%19.12Re\t%19.12Re\t%19.12Re\n", j, r0, i0, r1, i1);
  }
  fclose(fh);
  //nrm = sqrtl(nrm);
  mpfr_sqrt(nrm, nrm, MODE);
  mpfr_printf("Aberth iteration %3ld: %.12Re\n", k, nrm);
  verify_pade(res, rts, nD);              //  to be rewritten

  free_aberth(nD);
  /*
  fftwl_free(rts);
  fftwl_free(res);
  fftwl_free(t1);
  fftwl_free(t2);
  fftwl_free(wk);
  fftwl_free(ivc);
  */
}
