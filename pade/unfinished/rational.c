#include "header.h"

void compute_rational(unsigned long nD, unsigned long n_max_iter, mpfc_t *in) {
  pade_data.n_lins = 0;
  pade_data.n_poles = nD;


  allocate_pade(nD);
  for (unsigned long j = 0; j < pade_data.Npt-1; j++) {
    //W[j] = in[j+1] - in[0];
    mpfr_sub(W[j].re, in[j+1].re, in[0].re, MODE);
    mpfr_sub(W[j].im, in[j+1].im, in[0].im, MODE);
  }
  set_weight();  // OK?

  init_grams();
  gram_schmidt(nD); // looks like the problem is here
  evaluate_poly_array(nD);
  free_grams();

  //exit(0);	
  //find_l2_error(&pade_data); 
  errorl2(); 
  pade_data.n_lins++;

  for (unsigned int j = 1; j < n_max_iter; j++) {
    set_weight();
    init_grams();
    gram_schmidt(nD);
    evaluate_poly_array(nD);
    free_grams();
    //find_l2_error(&pade_data);  
    errorl2(); 
    pade_data.n_lins++;
  }
  set_weight();
  deallocate_pade(nD);
}


void errorl2(){
  pade_ptr p = &pade_data;
  mpfc_t tmp, buf;
  mpfr_inits2(precision, tmp.re, tmp.im, buf.re, buf.im, (mpfr_ptr) NULL);
  mpfr_inits2(precision, p->l2_nrm, p->l2_abs_err, p->l2_rel_err, (mpfr_ptr) NULL); 
//long double overN = 2.L*PI/N;
  mpfr_set_ui(p->l2_nrm, 0, MODE);
  mpfr_set_ui(p->l2_abs_err, 0, MODE);
  mpfr_set_ui(p->l2_rel_err, 0, MODE);

  for (unsigned int j = 0; j < pade_data.Npt-1; j++) {
    /*
    tmp = (arrayP[j]/arrayQ[j] - W[j]);    P/Q = P * barQ / Q*barQ
    inp->l2_nrm += creall(W[j]*conjl(W[j]));
    inp->l2_abs_err += creall(tmp*conjl(tmp));    
    */
    mpfc_div(&tmp, &arrP[j], &arrQ[j], MODE);
    mpfr_sub(tmp.re, tmp.re, W[j].re, MODE);
    mpfr_sub(tmp.im, tmp.im, W[j].im, MODE);

    mpfr_fma(p->l2_nrm, W[j].re, W[j].re, p->l2_nrm, MODE);
    mpfr_fma(p->l2_nrm, W[j].im, W[j].im, p->l2_nrm, MODE);

    mpfr_fma(p->l2_abs_err, tmp.re, tmp.re, p->l2_abs_err, MODE);
    mpfr_fma(p->l2_abs_err, tmp.im, tmp.im, p->l2_abs_err, MODE);
  }
  /*
  inp->l2_nrm = sqrtl(inp->l2_nrm)*overN;
  inp->l2_abs_err = sqrtl(inp->l2_abs_err)*overN;
  inp->l2_rel_err = (inp->l2_abs_err)/(inp->l2_nrm);
  */
  mpfr_t scale;
  mpfr_init2(scale, precision);
  mpfr_mul(scale, Ovn, Pie, MODE);
  mpfr_mul_ui(scale, scale, 2, MODE);

  mpfr_sqrt(p->l2_nrm, p->l2_nrm, MODE);
  mpfr_mul(p->l2_nrm, p->l2_nrm, scale, MODE);
 
  mpfr_sqrt(p->l2_abs_err, p->l2_abs_err, MODE);
  mpfr_mul(p->l2_abs_err, p->l2_abs_err, scale, MODE);

  mpfr_div(p->l2_rel_err, p->l2_abs_err, p->l2_nrm, MODE);
}












