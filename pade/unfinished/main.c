#include "header.h"

unsigned int nd = 4;
mpfr_prec_t precision = 128;
pade pade_data;
pade best_pade;
mpfc_t **Gees, **Cees;
mpfc_t **A, **B;
mpfc_t *arrQ, *arrP, *Tmp, *W;

mpfr_t *Ens, *M;
mpfr_t Pie, Ovn, Que;

int main() {
  pade_data.n_lins = 4;
  pade_data.n_poles = nd;
  pade_data.Npt = 256;
  init_memory(pade_data.n_poles);
  //allocate_pade(pade_data.n_poles);


  // Pade testing section
  mpfc_t buf, *test_array;
  mpfr_t s, t;

  test_array = malloc(pade_data.Npt*sizeof(mpfc_t));
  mpfr_inits2(precision, buf.re, buf.im, s, t, (mpfr_ptr) NULL);

  FILE *fh = fopen("./debug/test_array.txt", "w");
  fprintf(fh, "# 1. u 2. re W 3. im W\n\n"); 
  for (int j = 0; j < pade_data.Npt; j++) {
    mpfr_inits2(precision, test_array[j].re, test_array[j].im, (mpfr_ptr) NULL);
    mpfr_mul_si(s, Ovn, 2*j, MODE);
    mpfr_add_si(s, s, -1, MODE);
    mpfr_mul(t, s, Pie, MODE);
    mpfr_div_si(s, t, 2, MODE);
    mpfr_tan(s, s, MODE);

    mpfr_set(buf.re, s, MODE);
    mpfr_set_si(buf.im, -4, MODE);
    //mpfr_set(test_array[j].re, buf.re, MODE);
    //mpfr_set(test_array[j].im, buf.im, MODE);
    mpfc_si_div(&test_array[j], 1, &buf, MODE);
    mpfr_fprintf(fh, "%26.19Re\t%26.19Re\t%26.19Re\n", t, test_array[j].re, test_array[j].im);
  }
  fclose(fh);
  compute_rational(pade_data.n_poles, pade_data.n_lins, test_array);
  //optimal_pade("aberth_name", test_array);
  
  /*
  fh = fopen("./debug/PQ_array.txt", "w");
  fprintf(fh, "# 1. u 2.-3. P 4.-5. Q\n\n");
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_mul_si(s, Ovn, 2*(j+1), MODE);
    mpfr_add_si(s, s, -1, MODE);
    mpfr_mul(t, s, Pie, MODE);
    mpfr_fprintf(fh, "%26.19Re\t%26.19Re\t%26.19Re\t%26.19Re\t%26.19Re\n", t, arrP[j].re, arrP[j].im, arrQ[j].re, arrQ[j].im);
  }  
  fclose(fh);
  */
  
  //deallocate_pade(pade_data.n_poles);
  free(test_array);
  // End Pade testing section
 
  return 0;
}
