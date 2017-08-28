#include "header.h"


/*
void sort_by_imag(fftwl_complex *in1, fftwl_complex *in2, unsigned int nD) {
  fftwl_complex tmp1, tmp2;
  unsigned int n2 = nD, swapped;
  while (1) {
    swapped = 0;
    for (int j = 1; j < n2; j++) {
      if (cimagl(catanl(conf.scaling*in1[j-1])) > cimagl(catanl(conf.scaling*in1[j]))) {
        tmp1 = in1[j];
        tmp2 = in2[j];
        in1[j] = in1[j-1];
        in2[j] = in2[j-1];
        in1[j-1] = tmp1;
        in2[j-1] = tmp2;
	swapped = 1;
      }
    }
    n2 = n2 - 1;
    if (!swapped) break;
  }
}
*/
void sort_imag(mpfc_t *in1, mpfc_t *in2, unsigned int nD) {
  mpfc_t tmp1, tmp2;
  mpfc_t z1, z2;
  mpfr_inits2(precision, z1.re, z1.im, tmp1.re, tmp1.im, (mpfr_ptr) NULL);
  mpfr_inits2(precision, z2.re, z2.im, tmp2.re, tmp2.im, (mpfr_ptr) NULL);
  unsigned int n2 = nD; 
  unsigned int swapped;
  while (1) {
    swapped = 0;
    for (unsigned int j = 0; j < n2; j++) {
      if ( mpfr_cmp(in1[j-1].im, in1[j].im) > 0 ) {
        mpfr_set(tmp1.re, in1[j].re, MODE);
        mpfr_set(tmp1.im, in1[j].im, MODE);
        mpfr_set(tmp2.re, in2[j].re, MODE);
        mpfr_set(tmp2.im, in2[j].im, MODE);
        mpfr_set(in1[j].re, in1[j-1].re, MODE);
        mpfr_set(in1[j].im, in1[j-1].im, MODE);
        mpfr_set(in2[j].re, in2[j-1].re, MODE);
        mpfr_set(in2[j].im, in2[j-1].im, MODE);
        mpfr_set(in1[j-1].re, tmp1.re, MODE);
        mpfr_set(in1[j-1].im, tmp1.im, MODE);
        mpfr_set(in2[j-1].re, tmp2.re, MODE);
        mpfr_set(in2[j-1].im, tmp2.im, MODE);
	swapped = 1;        
      }
    }
    n2 = n2 - 1;
    if (!swapped) break;
  }
}


void verify_pade(mpfc_t *residues, mpfc_t *roots, unsigned int nD) {
   unsigned long 	N = pade_data.Npt;
   mpfc_t		buf;
   mpfr_t 		s, tmp1, tmp2; 
   mpfr_inits2(precision, s, tmp1, tmp2, (mpfr_ptr) NULL);
   mpfr_inits2(precision, buf.re, buf.im, (mpfr_ptr) NULL);
   mpfr_set_si(tmp1, 0, MODE);
   mpfr_set_si(tmp2, 0, MODE);
  
   mpfc_t *pde = malloc((N-1)*sizeof(mpfc_t));
   mpfc_t *rat = malloc((N-1)*sizeof(mpfc_t));
   for (unsigned int j = 0; j < N-1; j++) {
     mpfr_inits2(precision, pde[j].re, pde[j].im, (mpfr_ptr) NULL);
     mpfr_inits2(precision, rat[j].re, rat[j].im, (mpfr_ptr) NULL);
     mpfr_set_ui(pde[j].re, 0, MODE);
     mpfr_set_ui(pde[j].im, 0, MODE);
   }
   
   for (int k = 0; k < N-1; k++ ) {
     mpfr_mul_si(s, Ovn, k + 1, MODE);
     mpfr_add_d (s,   s,  -0.5, MODE);
     mpfr_mul   (s,   s,   Pie, MODE);
     mpfr_tan   (s,   s,  MODE);
     //s = PI*(2.L*(k+1)*overN - 1.L);
     //s = tanl(0.5L*s);
     //rat[k] = arrayP[k]/arrayQ[k];
     mpfc_div(&rat[k], &arrP[k], &arrQ[k], MODE);
     for (int j = 0; j < nD; j++) {
       //pde[k] += residues[j]/(s - roots[j]);
       mpfr_sub(buf.re, s, roots[j].re, MODE);
       mpfr_neg(buf.im, roots[j].im, MODE);
       //mpfc_sub(buf, s, roots[j], MODE);
       mpfc_div(&buf, &residues[j], &buf, MODE);
       mpfc_add(&pde[k], &pde[k], &buf, MODE);
     }
     //tmp1 += powl(cabsl(W[k] - pde[k]),2); 
     mpfr_sub(buf.re, W[k].re, pde[k].re, MODE);
     mpfr_sub(buf.im, W[k].im, pde[k].im, MODE);
     mpfr_fma(tmp1, buf.re, buf.re, tmp1, MODE);
     mpfr_fma(tmp1, buf.im, buf.im, tmp1, MODE);      
     
     //tmp2 += powl(cabsl(W[k] - arrayP[k]/arrayQ[k]),2);
     mpfc_div(&buf, &arrP[k], &arrQ[k], MODE);
     mpfr_sub(buf.re, W[k].re, buf.re, MODE);
     mpfr_sub(buf.im, W[k].im, buf.im, MODE);
     mpfr_fma(tmp2, buf.re, buf.re, tmp2, MODE);
     mpfr_fma(tmp2, buf.im, buf.im, tmp2, MODE);      
   }
   //tmp1 = sqrtl(tmp1)*overN;
   //tmp2 = sqrtl(tmp2)*overN;
   mpfr_sqrt(tmp1, tmp1, MODE);
   mpfr_sqrt(tmp2, tmp2, MODE);
   mpfr_mul(tmp1, tmp1, Ovn, MODE);
   mpfr_mul(tmp2, tmp2, Ovn, MODE);
   for (unsigned j = 0; j < N-1; j++) {
     mpfr_clears(pde[j].re, pde[j].im, (mpfr_ptr) NULL);
     mpfr_clears(rat[j].re, rat[j].im, (mpfr_ptr) NULL);
   }
   free(pde);
   free(rat);
   
   //printf("Rat. Approx L2 error:\t%.18LE\n", tmp2);
   //printf("Pole Approx L2 error:\t%.18LE\n", tmp1);
   //pde_complex_out("target.txt", W);
   //pde_complex_out("rational.txt", W);
   //pde_complex_out("pde.txt", pde);
}
