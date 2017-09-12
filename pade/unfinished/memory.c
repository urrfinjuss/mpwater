#include "header.h"

void init_memory(unsigned int nD) {
  Gees = malloc(4*sizeof(mpfc_t *));
  Cees = malloc(4*sizeof(mpfc_t *));
  Ens =  malloc(4*sizeof(mpfr_t));
  A = malloc(4*sizeof(mpfc_t *));
  B = malloc(4*sizeof(mpfc_t *));
  mpfr_inits2(precision, Ens[0], Ens[1], Ens[2], Ens[3], (mpfr_t *) NULL);
  /*
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_init2(M[j], precision);
  }
  */
  // set Ovn 
  mpfr_init2(Ovn, precision);
  mpfr_set_ui(Ovn, 1, MODE);
  mpfr_div_si(Ovn, Ovn, pade_data.Npt, MODE);
  // set Pie
  mpfr_init2(Pie, precision);
  mpfr_const_pi(Pie, MODE);
  init_output();
  //init_poly(nD);
  //init_grams();
}

void allocate_pade(unsigned long nD) {
  M =    malloc((pade_data.Npt-1)*sizeof(mpfr_t));
  W =    malloc((pade_data.Npt-1)*sizeof(mpfc_t));
  Tmp =  malloc((pade_data.Npt-1)*sizeof(mpfc_t));
  arrP = malloc((pade_data.Npt-1)*sizeof(mpfc_t));
  arrQ = malloc((pade_data.Npt-1)*sizeof(mpfc_t));
  for (int k = 0; k < 4; k++) {
    Cees[k] = malloc((2*nD + 1)*sizeof(mpfc_t));
    Gees[k] = malloc((pade_data.Npt-1)*sizeof(mpfc_t));
    A[k] =    malloc((pade_data.Npt-1)*sizeof(mpfc_t));
    B[k] =    malloc((pade_data.Npt-1)*sizeof(mpfc_t));
    for (int j = 0; j < pade_data.Npt-1; j++) {
      mpfr_inits2(precision, Gees[k][j].re, A[k][j].re, B[k][j].re, (mpfr_ptr) NULL);
      mpfr_inits2(precision, Gees[k][j].im, A[k][j].im, B[k][j].im, (mpfr_ptr) NULL);
    }
    for (int j = 0; j < 2*nD + 1; j++) {
      mpfr_init2(Cees[k][j].re, precision);
      mpfr_init2(Cees[k][j].im, precision);
    }
  }
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_inits2(precision, M[j], W[j].re, W[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(precision, arrQ[j].re, arrQ[j].im, (mpfr_ptr) NULL);
    mpfr_inits2(precision, arrP[j].re, arrP[j].im, (mpfr_ptr) NULL);
    ///mpfr_inits2(precision, arrQ[j].re, arrQ[j].im, (mpfr_ptr) NULL);
    /*
    mpfr_set_ui(arrQ[j].re, 1, MODE);
    mpfr_set_ui(arrQ[j].im, 0, MODE);
    mpfr_set_ui(arrP[j].re, 0, MODE);
    mpfr_set_ui(arrP[j].im, 0, MODE);
    */
  }
  // seed the initial poles across the plane:
  set_initial(nD);
}

void deallocate_pade(unsigned long nD) {
  for (int k = 3; k > -1; k--) {
    for (int j = 0; j < pade_data.Npt-1; j++) {
      mpfr_clears(Gees[k][j].re, A[k][j].re, B[k][j].re, (mpfr_ptr) NULL);
      mpfr_clears(Gees[k][j].im, A[k][j].im, B[k][j].im, (mpfr_ptr) NULL);
    }
    for (unsigned int j = 0; j < 2*nD + 1; j++) {
      mpfr_clears(Cees[k][j].re, Cees[k][j].im, (mpfr_ptr) NULL);
    }
    free(Gees[k]);
    free(Cees[k]);
    free(A[k]);
    free(B[k]);
  } 
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_clears(arrP[j].re, arrP[j].im, (mpfr_ptr) NULL);
    mpfr_clears(arrQ[j].re, arrQ[j].im, (mpfr_ptr) NULL);
    //mpfr_clears(Tmp[j].re, W[j].re, (mpfr_ptr) NULL);
    //mpfr_clears(Tmp[j].im, W[j].im, (mpfr_ptr) NULL);
    mpfr_clear(M[j]);
  }
  free(arrP);
  free(arrQ);
  free(Tmp);
  free(W);
  free(M);
}

