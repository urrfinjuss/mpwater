#include "header.h"

//static mpfr_t que;

void mpfc_dotpr(mpfc_t *out, mpfc_t *op1, mpfc_t *op2) {
  mpfr_t tmpr, tmpi;
  mpfr_init2(tmpr, precision);
  mpfr_init2(tmpi, precision);
  mpfr_set_ui(out->re, 0, MODE);
  mpfr_set_ui(out->im, 0, MODE);
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_mul(tmpr, op1[j].im, op2[j].im, MODE);
    mpfr_fma(tmpr, op1[j].re, op2[j].re, tmpr, MODE);
    mpfr_fma(out->re, tmpr, M[j], out->re, MODE);
    mpfr_mul(tmpi, op1[j].re, op2[j].im, MODE);
    mpfr_fms(tmpi, op1[j].im, op2[j].re, tmpi, MODE);
    mpfr_fma(out->im, tmpi, M[j], out->im, MODE);
  }
  mpfr_clear(tmpr);
  mpfr_clear(tmpi);
}

void set_weight() {
  mpfr_t tmp, que;
  mpfr_init2(que, precision);
  mpfr_init2(tmp, precision);
  FILE *fh = fopen("./debug/M_array.txt", "w");
  for (long int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_mul_si(que, Ovn, 2*(j+1), MODE); 
    mpfr_add_si(que, que, -1, MODE);
    mpfr_mul(que, que, Pie, MODE);
    mpfr_mul(tmp, arrQ[j].im, arrQ[j].im, MODE);
    mpfr_fma(tmp, arrQ[j].re, arrQ[j].re, tmp, MODE);
    mpfr_d_div(M[j], 0.5, tmp, MODE);
    mpfr_div_ui(tmp, que, 2, MODE);
    mpfr_cos(tmp, tmp, MODE);
    mpfr_mul(tmp, tmp, tmp, MODE);
    mpfr_div(M[j], M[j], tmp, MODE);
    mpfr_fprintf(fh, "%24.17Re\t%24.17Re\n", que, M[j]);
    //M[j] = 0.5L/powl(cosl(0.5L*q), 2)/creall(arrayQ[j]*conjl(arrayQ[j]));
  }
  fclose(fh);
  mpfr_clear(que);
  mpfr_clear(tmp);
}

void set_initial(unsigned int nD) {
  mpfc_t tmp, w, buf;
  mpfr_t phase, que;
  mpfr_inits2(precision, w.re, buf.re, tmp.re, phase, que, (mpfr_ptr) NULL);
  mpfr_inits2(precision, w.im, buf.im, tmp.im, (mpfr_ptr) NULL);
  mpfr_set_ui(tmp.re, 0, MODE);
  mpfr_set_ui(tmp.im, 0, MODE);
  for (long int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_set_si(arrQ[j].re, 1, MODE);
    mpfr_set_si(arrQ[j].im, 0, MODE);
    mpfr_set_si(arrP[j].re, 0, MODE);
    mpfr_set_si(arrP[j].im, 0, MODE);
  }

  for (int k = 0; k < nD; k++) {
    mpfr_mul_si(phase, Pie, 2*k + 1, MODE);
    mpfr_div_ui(phase, phase, 4*nD, MODE);
    // set on circle
    // mpfr_sin_cos(tmp.im, tmp.re, phase, MODE);
    // end circle
    // set on the line
    mpfr_exp(tmp.re, phase, MODE);           
    mpfr_add_si(tmp.re, tmp.re, -1, MODE);
    //mpfr_set_si(tmp.im, 0, MODE);
    mpfr_printf("tmp %Re\t%Re\n", tmp.re, tmp.im);
    // end line
    for (long int j = 0; j < pade_data.Npt-1; j++) {
      mpfr_mul_si(que, Ovn, j+1, MODE);
      mpfr_add_d(que, que, -0.5, MODE);
      mpfr_mul(que, que, Pie, MODE);
      mpfr_tan(que, que, MODE);
      
      //mpfr_sub(w.re, que, tmp.im, MODE);
      //mpfr_add(w.im, que, tmp.re, MODE);
      mpfr_set(w.re, que, MODE);
      mpfr_neg(w.im, tmp.re, MODE);
      
      //mpfr_printf("que %Re\nPie %Re\nw %Re %Re\n", que, Pie, w.re, w.im);
      
      mpfc_mul(&buf, &arrQ[j], &w, MODE); 
      mpfc_set(&arrQ[j], &buf, MODE);
      /*
      mpfr_fms(arrQ[j].re, arrQ[j].im, w.im, arrQ[j].re, MODE);   
      mpfr_fms(arrQ[j].re, arrQ[j].re, w.re, arrQ[j].re, MODE);
      mpfr_fma(arrQ[j].im, arrQ[j].im, w.re, arrQ[j].im, MODE);
      mpfr_fma(arrQ[j].im, arrQ[j].re, w.im, arrQ[j].im, MODE);
      */
      
      //mpfr_printf("%Re\t%Re\n", arrQ[j].re, arrQ[j].im);
    }
  }

  FILE *fh = fopen("./debug/PQ_initial.txt", "w");
  fprintf(fh, "# 1. u 2.-3. Q\n\n");
  for (long int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_mul_si(que, Ovn, j+1, MODE);
    mpfr_add_d(que, que, -0.5, MODE);
    mpfr_mul(que, que, Pie, MODE);
    mpfr_fprintf(fh, "%Re\t%Re\t%Re\n", que, arrQ[j].re, arrQ[j].im);
  }
  fclose(fh);
}

