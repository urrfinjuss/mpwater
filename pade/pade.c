#include "header.h"

/*
void optimal_pade(char *str, fftwl_complex *in) {
//  FILE *fh = fopen("pc_rate.txt","w");
  unsigned long nd = 1, l_iters = 12;
  pade best_pade;  
  best_pade.l2_rel_err = 1.L;
  
  //fprintf(fh, "# 1. nD, number poles 2. Error\n\n");
  //for (unsigned int nd = 1; nd < 32; nd++) {
  while (nd < 32) {
    nd++; //nd++;
    best_pade.n_poles = nd;
    pade_data.n_lins = 0;
    compute_rational(nd, l_iters, in);
    deallocate_pade();
    //if ((pade_data.l2_rel_err < best_pade.l2_rel_err)&&(best_pade.l2_rel_err > 4.0E-9L)) {
    if (best_pade.l2_rel_err > 8.0E-10L) {
      best_pade.l2_rel_err = pade_data.l2_rel_err;
      best_pade.l2_abs_err = pade_data.l2_abs_err;
      best_pade.l2_nrm = pade_data.l2_nrm;
      best_pade.n_poles = nd;
      best_pade.n_lins = l_iters;
      //printf("Relative Error (nd = %3lu) = %11.5LE\n", nd, best_pade.l2_rel_err);
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
      alt_map.image_offset = 0.L;
      alt_map.origin_offset = 0.L;
      alt_map.scaling = conf.scaling/sqrtl(2.0L);
      // end set
      deallocate_pade();
      best_pade.l2_rel_err = pade_data.l2_rel_err;
      best_pade.l2_abs_err = pade_data.l2_abs_err;
      best_pade.l2_nrm = pade_data.l2_nrm;
      best_pade.n_poles = nd;
      best_pade.n_lins = l_iters;
      printf("Safe Relative Error (nd = %3lu) = %11.5LE\n", nd, best_pade.l2_rel_err);
      break;
    } 
  }
  //fclose(fh);
  //exit(1);
}
*/
/*
void verify_pade(fftwl_complex *residues, fftwl_complex *roots, unsigned int nD) {
   unsigned long 	N = state.number_modes;
   long double 		overN = 1.L/N, s;
   long double		tmp_1 = 0.L;
   long double		tmp_2 = 0.L;
   
   fftwl_complex *pade = fftwl_malloc((N-1)*sizeof(fftwl_complex));
   fftwl_complex *rat  = fftwl_malloc((N-1)*sizeof(fftwl_complex));
   memset(pade, 0, (N-1)*sizeof(fftwl_complex));
   for (int k = 0; k < N-1; k++ ) {
     s = PI*(2.L*(k+1)*overN - 1.L);
     s = tanl(0.5L*s);
     rat[k] = arrayP[k]/arrayQ[k];
     for (int j = 0; j < nD; j++) {
       pade[k] += residues[j]/(s - roots[j]);
     }
     tmp_1 += powl(cabsl(W[k] - pade[k]),2); 
     tmp_2 += powl(cabsl(W[k] - arrayP[k]/arrayQ[k]),2);
   }
   tmp_1 = sqrtl(tmp_1)*overN;
   tmp_2 = sqrtl(tmp_2)*overN;
   //printf("Rat. Approx L2 error:\t%.18LE\n", tmp_2);
   //printf("Pole Approx L2 error:\t%.18LE\n", tmp_1);
   //pade_complex_out("target.txt", W);
   //pade_complex_out("rational.txt", W);
   //pade_complex_out("pade.txt", pade);
   fftwl_free(rat);
   fftwl_free(pade);
}
*/

/*

void poly_val_array(fftwl_complex *in, unsigned long nD, fftwl_complex *outQ, fftwl_complex *outQp, fftwl_complex *outP){
  fftwl_complex  	*tmp_C, *tmp_B, *tmp_A;
  fftwl_complex		*P[4], *Q[4], *dQ[4];
  fftwl_complex		*s;

  tmp_A = fftwl_malloc(nD*sizeof(fftwl_complex));
  tmp_B = fftwl_malloc(nD*sizeof(fftwl_complex));
  tmp_C = fftwl_malloc(nD*sizeof(fftwl_complex));
  for (int j = 0; j < 4; j++) {
    P[j] = fftwl_malloc(nD*sizeof(fftwl_complex));
    Q[j] = fftwl_malloc(nD*sizeof(fftwl_complex));
    dQ[j] = fftwl_malloc(nD*sizeof(fftwl_complex));
  }
  s = fftwl_malloc(nD*sizeof(fftwl_complex));

  for (int j = 0; j < nD; j++) {
    //s[j] = ctanl(0.5L*in[j]);
    s[j] = in[j];
    P[0][j] = 0.L;
    Q[0][j] = 1.L;
    dQ[0][j] = 0.L;

    P[1][j] = -1.L;
    Q[1][j] = -Cees[0][1];
    dQ[1][j] = 0.L;

    P[2][j]  = s[j]*P[0][j] - Cees[0][2]*P[1][j] - Cees[1][2]*P[0][j];
    Q[2][j]  = s[j]*Q[0][j] - Cees[0][2]*Q[1][j] - Cees[1][2]*Q[0][j];
    dQ[2][j] = Q[0][j] + s[j]*dQ[0][j] - Cees[0][2]*dQ[1][j] - Cees[1][2]*dQ[0][j];

    P[3][j]  = s[j]*P[1][j] - Cees[0][3]*P[2][j] - Cees[1][3]*P[1][j] - Cees[2][3]*P[0][j];
    Q[3][j]  = s[j]*Q[1][j] - Cees[0][3]*Q[2][j] - Cees[1][3]*Q[1][j] - Cees[2][3]*Q[0][j];
    dQ[3][j] = Q[1][j] + s[j]*dQ[1][j] - Cees[0][3]*dQ[2][j] - Cees[1][3]*dQ[1][j] - Cees[2][3]*dQ[0][j];
    for (unsigned long k = 4; k < 2*nD + 1; k++) {
      tmp_C[j] = s[j]*P[2][j] - Cees[0][k]*P[3][j] - Cees[1][k]*P[2][j] - Cees[2][k]*P[1][j] - Cees[3][k]*P[0][j];
      tmp_B[j] = s[j]*Q[2][j] - Cees[0][k]*Q[3][j] - Cees[1][k]*Q[2][j] - Cees[2][k]*Q[1][j] - Cees[3][k]*Q[0][j];
      tmp_A[j] = Q[2][j] + s[j]*dQ[2][j] - Cees[0][k]*dQ[3][j] - Cees[1][k]*dQ[2][j] - Cees[2][k]*dQ[1][j] - Cees[3][k]*dQ[0][j];
      P[0][j] = P[1][j];
      P[1][j] = P[2][j];
      P[2][j] = P[3][j];
      P[3][j] = tmp_C[j];
      Q[0][j] = Q[1][j];
      Q[1][j] = Q[2][j];
      Q[2][j] = Q[3][j];
      Q[3][j] = tmp_B[j];
      dQ[0][j] = dQ[1][j];
      dQ[1][j] = dQ[2][j];
      dQ[2][j] = dQ[3][j];
      dQ[3][j] = tmp_A[j]; 
    } 
    outP[j]  = P[3][j];  
    outQ[j]  = Q[3][j];
    outQp[j] = dQ[3][j];
  }
  
  fftwl_free(tmp_A);
  fftwl_free(tmp_B);
  fftwl_free(tmp_C);
  for (int j = 3; j > -1; j--) {
    fftwl_free(P[j]);
    fftwl_free(Q[j]);
    fftwl_free(dQ[j]);
  }
  fftwl_free(s);
}

*/
