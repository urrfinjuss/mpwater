#include "header.h"


void newton_search(unsigned long nD) {
  mpfc_t *crd, *rsd;
  mpfc_t z, p, f, df, ivs;
  mpfc_t buf, eye;
  mpfr_t tmp;
  long int counter = 0;
  char str[80];
  FILE *fh;

  // init MPFR variables
  crd = malloc(nD*sizeof(mpfc_t));
  rsd = malloc(nD*sizeof(mpfc_t));
  for (int j = 0; j < nD; j++) {
    mpfr_inits2(precision, crd[j].re, rsd[j].re, (mpfr_ptr) NULL);
    mpfr_inits2(precision, crd[j].im, rsd[j].im, (mpfr_ptr) NULL);
  }
  mpfr_inits2(precision, z.re, p.re, f.re, df.re, ivs.re, (mpfr_ptr) NULL);
  mpfr_inits2(precision, z.im, p.im, f.im, df.im, ivs.im, (mpfr_ptr) NULL);
  mpfr_inits2(precision, buf.re, buf.im, tmp, (mpfr_ptr) NULL);
  mpfr_inits2(precision, eye.re, eye.im, (mpfr_ptr) NULL);
  // set MPFR variables
  mpfc_set_si(&z,   0,   0, MODE);
  mpfc_set_si(&p,   0,   0, MODE);
  mpfc_set_si(&f,   1,   0, MODE);
  mpfc_set_si(&eye, 1,   0, MODE);
  // start Newton method
  evaluate_poly(&z, nD, &f, &df, &p);  
  for (int j = 0; j < nD; j++) {  
    sprintf(str, "./crd/root_%03d.txt", j);
    fh = fopen(str, "w");
    /*
    z = 0.5IL;
    f = 1.L; 
    ivs = 1.L;
    */
    mpfc_set_si(&z, 0, 1, MODE); 
    mpfc_set_si(&f, 1, 0, MODE);
    mpfc_set_si(&ivs, 1, 0, MODE);
    mpfc_mul(&buf, &ivs, &f, MODE);
    mpfr_mul(tmp, buf.re, buf.re, MODE);
    mpfr_fma(tmp, buf.im, buf.im, tmp, MODE);
    mpfr_sqrt(tmp, tmp, MODE);
    //while ( cabsl(f*ivs) > 5.0E-29L) {
    while ( mpfr_cmp_ld(tmp, 5.0E-14L) > 0) {
      evaluate_poly(&z, nD, &f, &df, &p);
      if ( j == 0 ) {
        /*
	ivs = 1.L;
        z = z - f/df;
        */
        mpfc_set_si(&ivs, 1, 0, MODE);
        mpfc_div(&buf, &f, &df, MODE);
        mpfc_sub(&z, &z, &buf, MODE);
      } else {     
        /*
        ivs = 0.L;
        for (int k = 0; k < j; k++) {
          ivs += 1.L/(z - crd[k]);
        }
        z = z - f/(df - f*ivs);
        equiv z = z + f/(f*ivs - df)
        */
        mpfc_set_si(&ivs, 0, 0, MODE);
        for (int k = 0; k < j; k++) {
          mpfc_sub(&buf, &z, &crd[k], MODE);
          mpfc_div(&buf, &eye, &buf, MODE);
	  mpfc_add(&ivs, &ivs, &buf, MODE); 
        }
        mpfc_fms(&buf, &f, &ivs, &df, MODE);
        mpfc_div(&buf, &f, &buf, MODE);
        mpfc_add(&z,   &z, &buf, MODE);

	// tmp <- cabsl(f*ivs)
        mpfc_mul(&buf, &ivs, &f, MODE);
        mpfr_mul(tmp, buf.re, buf.re, MODE);
        mpfr_fma(tmp, buf.im, buf.im, tmp, MODE);
        mpfr_sqrt(tmp, tmp, MODE);
      }
      counter++;
      mpfr_fprintf(fh, "%4ld\t%26.18Re\t%26.18Re\n", counter, z.re, z.im);
      if (counter == 512) break;
    }
    mpfr_printf("Root %3d: Iteration %4ld |f| = %14.8Re\n", j, counter, tmp);
    counter = 0;
    fclose(fh);
    //rsd[j] = p/df;
    //crd[j] = z;
    mpfc_div(&rsd[j], &p, &df, MODE);
    mpfc_set(&crd[j], &z, MODE);
  }
  verify_pade(rsd, crd, nD);
  fh = fopen("./crd/crd.txt", "w");
  fprintf(fh, "# 1. root # 2.-3. z_k 4.-5. g_k\n\n");
  for (int j = 0; j < nD; j++) {
    mpfr_fprintf(fh, "%4d\t%26.18Re\t%26.18Re\t%26.18Re\t%26.18Re\n",j, crd[j].re, crd[j].im, rsd[j].re, rsd[j].im);
  }
  fclose(fh);
  
  for (int j = 0; j < nD; j++) {
    mpfr_clears(crd[j].re, rsd[j].re, (mpfr_ptr) NULL);
    mpfr_clears(crd[j].im, rsd[j].im, (mpfr_ptr) NULL);
  }
  free(crd);
  free(rsd);
}
