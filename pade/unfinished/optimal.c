#include "header.h"


void optimal_pade(char *str, mpfc_t *in) {
//  FILE *fh = fopen("pc_rate.txt","w");
  unsigned long nd = 4, l_iters = 12;
  //pade best_pade;  
  mpfr_init2(best_pade.l2_rel_err, precision);
  mpfr_set_ui(best_pade.l2_rel_err, 1, MODE);
  //best_pade.l2_rel_err = 1.L;
  
  //fprintf(fh, "# 1. nD, number poles 2. Error\n\n");
  //for (unsigned int nd = 1; nd < 32; nd++) {
  while (nd < 32) {
    nd++; //nd++;
    best_pade.n_poles = nd;
    pade_data.n_lins = 0;
    compute_rational(nd, l_iters, in);
    printf("Compute rational ... Done!");
    //if ((pade_data.l2_rel_err < best_pade.l2_rel_err)&&(best_pade.l2_rel_err > 4.0E-9L)) {
    if (mpfr_cmp_ld(best_pade.l2_rel_err, 8.0E-10L) > 0 ) {
    //if (best_pade.l2_rel_err > 8.0E-10L) {
      mpfr_set(best_pade.l2_rel_err, pade_data.l2_rel_err, MODE);
      mpfr_set(best_pade.l2_abs_err, pade_data.l2_abs_err, MODE);
      mpfr_set(best_pade.l2_nrm, pade_data.l2_nrm, MODE);
      /*
      best_pade.l2_rel_err = pade_data.l2_rel_err;
      best_pade.l2_abs_err = pade_data.l2_abs_err;
      best_pade.l2_nrm = pade_data.l2_nrm;
      */
      best_pade.n_poles = nd;
      best_pade.n_lins = l_iters;
      
      mpfr_printf("Relative Error (nd = %3lu) = %11.5Re\n", nd, best_pade.l2_rel_err);
      //fprintf(fh, "%.3lu\t%11.5LE\n", nd, best_pade.l2_rel_err);
    } else {
      nd = nd-1;
      //nd = 6;
      best_pade.n_poles = nd;
      pade_data.n_lins = 0;
      compute_rational(nd, l_iters, in);
      aberth_iter(nd, str);
      //newton_search(nd);
      // set new map parameters
      // end set
      deallocate_pade(nd);
      mpfr_set(best_pade.l2_rel_err, pade_data.l2_rel_err, MODE);
      mpfr_set(best_pade.l2_abs_err, pade_data.l2_abs_err, MODE);
      mpfr_set(best_pade.l2_nrm, pade_data.l2_nrm, MODE);
      best_pade.n_poles = nd;
      best_pade.n_lins = l_iters;
      mpfr_printf("Safe Relative Error (nd = %3lu) = %11.5Re\n", nd, best_pade.l2_rel_err);
      break;
    } 
  }
  //fclose(fh);
  //exit(1);
}
