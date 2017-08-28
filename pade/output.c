#include "header.h"
static mpfr_t que;

void init_output() {
  mpfr_init2(que, precision);
  mpfr_init2(Pie, precision);
  mpfr_init2(Ovn, precision);
  mpfr_const_pi(Pie, MODE);
  mpfr_set_ui(Ovn, 1, MODE);
  mpfr_div_ui(Ovn, Ovn, pade_data.Npt, MODE); 
}

void pade_real_out(char *fname, mpfr_t *in) {
  FILE *fh = fopen(fname, "w");
  fprintf(fh, "# 1. q 2. Array\n\n");
  for (long int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_mul_si(que, Ovn, 2*(j+1), MODE); 
    mpfr_add_si(que, que, -1, MODE);
    mpfr_mul(que, que, Pie, MODE);
    mpfr_fprintf(fh, "%46.39Re\t%46.39Re\n", que, in[j]);
  }
  fclose(fh);
}

void pade_complex_out(char *fname, mpfc_t *in) {
  FILE *fh = fopen(fname, "w");
  fprintf(fh, "# 1. q 2. re 3. im\n\n");
  for (long int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_mul_si(que, Ovn, 2*(j+1), MODE); 
    mpfr_add_si(que, que, -1, MODE);
    mpfr_mul(que, que, Pie, MODE);
    mpfr_fprintf(fh, "%46.39Re\t%46.39Re\t%46.39Re\n", que, in[j].re, in[j].im);
  }
  fclose(fh);
}

void print_pade() { 
  printf("Iter #%3u:\t", pade_data.n_lins); 
  mpfr_printf("%.18LE\n", pade_data.l2_rel_err); 
}
