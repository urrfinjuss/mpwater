#include "header.h"

static mpfc_t *tmp;
static mpfc_t s;

void init_grams(){
  mpfr_inits2(precision, s.re, s.im, (mpfr_ptr) NULL);
  tmp = malloc((pade_data.Npt-1)*sizeof(mpfc_t));
  for (unsigned long j = 0; j < pade_data.Npt-1; j++) {
    mpfr_inits2(precision, tmp[j].re, tmp[j].im, (mpfr_ptr) NULL);
  }
}

void free_grams(){
  mpfr_clears(s.re, s.im, (mpfr_ptr) NULL);
  for (unsigned long j = 0; j < pade_data.Npt-1; j++) {
    mpfr_clears(tmp[j].re, tmp[j].im, (mpfr_ptr) NULL);
  }
  free(tmp); 
}

void fmul_sub(mpfc_t *in1, mpfc_t *in2, mpfc_t *cnum, mpfc_t *out) {
  mpfc_t buf;
  mpfr_init2(buf.re, precision);
  mpfr_init2(buf.im, precision);
  for (unsigned long j = 0; j < pade_data.Npt-1; j++) {
    //out[j] = in1[j] - cnum*in2[j];
    mpfc_fms(&buf, cnum, &in2[j], &in1[j], MODE);
    mpfr_neg(out[j].re, buf.re, MODE);
    mpfr_neg(out[j].im, buf.im, MODE);
  }
  mpfr_clear(buf.re);
  mpfr_clear(buf.im);
}

void sigma_mul(mpfc_t *in, mpfc_t *out) {
  for (unsigned long j = 0; j < pade_data.Npt-1; j++) {
    //s = tanl(0.5L*PI*(2.L*(j+1.L)*overN - 1.L));
    mpfr_mul_si(s.re, Ovn, j+1, MODE);
    mpfr_add_d(s.re, s.re, -0.5, MODE);
    mpfr_mul(s.re, s.re, Pie, MODE);
    mpfr_tan(s.re, s.re, MODE);
    mpfr_set_si(s.im, 0, MODE);
    //out[j] = s*in[j];
    mpfc_mul(&out[j], &s, &in[j], MODE);
  }
}

void gram_schmidt(unsigned long nD){
  mpfc_t buf;
  mpfr_inits2(precision, buf.re, buf.im, (mpfr_ptr) NULL);
  for (unsigned k = 0; k < 4; k++) {
    mpfr_set_si(Ens[0], 0, MODE);
    // memset(Cees[k], 0, (2*nD + 1)*sizeof(fftwl_complex));
    for (unsigned long j = 0; j < 2*nD + 1; j++) {
      mpfr_set_si(Cees[k][j].re, 0, MODE);
      mpfr_set_si(Cees[k][j].im, 0, MODE);
    }
  }
  // step 0
  //memcpy(Gees[0], W, (N-1)*sizeof(fftwl_complex));
  for (unsigned long j = 0; j < pade_data.Npt-1; j++) {
    mpfr_set(Gees[0][j].re, W[j].re, MODE);
    mpfr_set(Gees[0][j].im, W[j].im, MODE);
  }
  //Ens[0] = 1.0L/dot(Gees[0],Gees[0]);
  mpfc_dotpr(&buf, Gees[0], Gees[0]);
  mpfr_ui_div(Ens[0], 1, buf.re, MODE);


  /*
  printf("before\n");   
  mpfr_printf("Ens[0] = %Re\nbuf.re = %Re\n", Ens[0], buf.re);
  FILE *fh = fopen("./debug/Gees0.txt", "w");
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_fprintf(fh, "%Re\t%Re\n", Gees[0][j].re, Gees[0][j].im);
  }
  fclose(fh);
  printf("after\n");
  */

  //printf("Gram-Schmidt %ld\n", 9);
  // step 1
  //for (unsigned long j = 0; j < N-1; j++) tmp[j] = 1.L;
  for (unsigned long j = 0; j < pade_data.Npt-1; j++) {
    mpfr_set_si(tmp[j].re, 1, MODE);
    mpfr_set_si(tmp[j].im, 0, MODE);
  }
  //printf("Gram-Schmidt %ld\n", 10);
  //Cees[0][1] = dot(tmp, Gees[0])*Ens[0];
  //fmul_sub(tmp, Gees[0], Cees[0][1], Gees[1]);
  //Ens[1] = 1.0L/dot(Gees[1],Gees[1]);
  mpfc_dotpr(&buf, tmp, Gees[0]);
  mpfr_mul(Cees[0][1].re, buf.re, Ens[0], MODE);
  mpfr_mul(Cees[0][1].im, buf.im, Ens[0], MODE);
  fmul_sub(tmp, Gees[0], &Cees[0][1], Gees[1]);
  mpfc_dotpr(&buf, Gees[1], Gees[1]);
  mpfr_ui_div(Ens[1], 1, buf.re, MODE);
  
  /*
  printf("before\n");   
  mpfr_printf("Cees[0][1] = %Re\t%Re\n", Cees[0][1].re, Cees[0][1].im);
  mpfr_printf("Ens[1] = %Re\nbuf.re = %Re\n", Ens[1], buf.re);
  fh = fopen("./debug/Gees1.txt", "w");
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_fprintf(fh, "%Re\t%Re\t%Re\t%Re\n", Gees[1][j].re, Gees[1][j].im, tmp[j].re, tmp[j].im);
  }
  fclose(fh);
  printf("after\n");
  */
  // step 2 
  /*
  sigma_mul(Gees[0], tmp);
  Cees[0][2] = dot(tmp, Gees[1])*Ens[1];
  fmul_sub(tmp, Gees[1], Cees[0][2], tmp);
  Cees[1][2] = dot(tmp, Gees[0])*Ens[0];
  fmul_sub(tmp, Gees[0], Cees[1][2], Gees[2]);
  Ens[2] = 1.0L/dot(Gees[2],Gees[2]);
  */
  //printf("Gram-Schmidt %ld\n", 11);
  
  sigma_mul(Gees[0], tmp);
  mpfc_dotpr(&buf, tmp, Gees[1]);
  mpfr_mul(Cees[0][2].re, buf.re, Ens[1], MODE);
  mpfr_mul(Cees[0][2].im, buf.im, Ens[1], MODE);
  /*
  fh = fopen("./debug/SigmaGees0.txt", "w");
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_fprintf(fh, "%Re\t%Re\t%Re\t%Re\n", Gees[2][j].re, Gees[2][j].im, tmp[j].re, tmp[j].im);
  }
  fclose(fh);
  */
  fmul_sub(tmp, Gees[1], &Cees[0][2], tmp);
  mpfc_dotpr(&buf, tmp, Gees[0]);
  mpfr_mul(Cees[1][2].re, buf.re, Ens[0], MODE);
  mpfr_mul(Cees[1][2].im, buf.im, Ens[0], MODE);
  fmul_sub(tmp, Gees[0], &Cees[1][2], Gees[2]);

  mpfc_dotpr(&buf, Gees[2], Gees[2]);
  mpfr_ui_div(Ens[2], 1, buf.re, MODE);
  
  /*
  printf("before\n");   
  mpfr_printf("Cees[0][2] = %Re\t%Re\n", Cees[0][2].re, Cees[0][2].im);
  mpfr_printf("Cees[1][2] = %Re\t%Re\n", Cees[1][2].re, Cees[1][2].im);
  mpfr_printf("Ens[2] = %Re\nbuf.re = %Re\n", Ens[2], buf.re);
  fh = fopen("./debug/Gees2.txt", "w");
  for (int j = 0; j < pade_data.Npt-1; j++) {
    mpfr_fprintf(fh, "%Re\t%Re\t%Re\t%Re\n", Gees[2][j].re, Gees[2][j].im, tmp[j].re, tmp[j].im);
  }
  fclose(fh);
  printf("after\n");
  */
  // step 3
  /*sigma_mul(Gees[1], tmp);
  Cees[0][3] = dot(tmp, Gees[2])*Ens[2];
  fmul_sub(tmp, Gees[2], Cees[0][3], tmp);

  Cees[1][3] = dot(tmp, Gees[1])*Ens[1];
  fmul_sub(tmp, Gees[1], Cees[1][3], tmp);

  Cees[2][3] = dot(tmp, Gees[0])*Ens[0];
  fmul_sub(tmp, Gees[0], Cees[2][3], Gees[3]);

  Ens[3] = 1.0L/dot(Gees[3],Gees[3]);*/

  sigma_mul(Gees[1], tmp);
  mpfc_dotpr(&buf, tmp, Gees[2]);
  mpfr_mul(Cees[0][3].re, buf.re, Ens[2], MODE);
  mpfr_mul(Cees[0][3].im, buf.im, Ens[2], MODE);
  fmul_sub(tmp, Gees[2], &Cees[0][3], tmp);

  mpfc_dotpr(&buf, tmp, Gees[1]);
  mpfr_mul(Cees[1][3].re, buf.re, Ens[1], MODE);
  mpfr_mul(Cees[1][3].im, buf.im, Ens[1], MODE);
  fmul_sub(tmp, Gees[1], &Cees[1][3], tmp);

  mpfc_dotpr(&buf, tmp, Gees[0]);
  mpfr_mul(Cees[2][3].re, buf.re, Ens[0], MODE);
  mpfr_mul(Cees[2][3].im, buf.im, Ens[0], MODE);
  fmul_sub(tmp, Gees[0], &Cees[2][3], Gees[3]);

  mpfc_dotpr(&buf, Gees[3], Gees[3]);
  mpfr_ui_div(Ens[3], 1, buf.re, MODE);
  // step 4+
  for (long int j = 4; j < 2*nD + 1; j++) {
      //mpfr_printf("Ens[0] = %Re\n", Ens[0]);
      /*
      sigma_mul(Gees[2], tmp);
      Cees[0][j] = dot(tmp, Gees[3])*Ens[3];
      fmul_sub(tmp, Gees[3], Cees[0][j], tmp);

      Cees[1][j] = dot(tmp, Gees[2])*Ens[2];
      fmul_sub(tmp, Gees[2], Cees[1][j], tmp);

      Cees[2][j] = dot(tmp, Gees[1])*Ens[1];
      fmul_sub(tmp, Gees[1], Cees[2][j], tmp);

      Cees[3][j] = dot(tmp, Gees[0])*Ens[0];
      fmul_sub(tmp, Gees[0], Cees[3][j], tmp);
      */

      sigma_mul(Gees[2], tmp);
      mpfc_dotpr(&buf, tmp, Gees[3]);
      mpfr_mul(Cees[0][j].re, buf.re, Ens[3], MODE);
      mpfr_mul(Cees[0][j].im, buf.im, Ens[3], MODE);
      fmul_sub(tmp, Gees[3], &Cees[0][j], tmp);
     
      mpfc_dotpr(&buf, tmp, Gees[2]);
      mpfr_mul(Cees[1][j].re, buf.re, Ens[2], MODE);
      mpfr_mul(Cees[1][j].im, buf.im, Ens[2], MODE);
      fmul_sub(tmp, Gees[2], &Cees[1][j], tmp);

      mpfc_dotpr(&buf, tmp, Gees[1]);
      mpfr_mul(Cees[2][j].re, buf.re, Ens[1], MODE);
      mpfr_mul(Cees[2][j].im, buf.im, Ens[1], MODE);
      fmul_sub(tmp, Gees[1], &Cees[2][j], tmp);
    
      mpfc_dotpr(&buf, tmp, Gees[0]);
      mpfr_mul(Cees[3][j].re, buf.re, Ens[0], MODE);
      mpfr_mul(Cees[3][j].im, buf.im, Ens[0], MODE);
      fmul_sub(tmp, Gees[0], &Cees[3][j], tmp);

      /* 
      memcpy(Gees[0], Gees[1], (N-1)*sizeof(fftwl_complex));
      memcpy(Gees[1], Gees[2], (N-1)*sizeof(fftwl_complex));
      memcpy(Gees[2], Gees[3], (N-1)*sizeof(fftwl_complex));
      memcpy(Gees[3], tmp, (N-1)*sizeof(fftwl_complex));
      */

      for (unsigned long k = 0; k < pade_data.Npt-1; k++) {
        mpfr_set(Gees[0][k].re, Gees[1][k].re, MODE); 
        mpfr_set(Gees[0][k].im, Gees[1][k].im, MODE); 
        mpfr_set(Gees[1][k].re, Gees[2][k].re, MODE); 
        mpfr_set(Gees[1][k].im, Gees[2][k].im, MODE); 
        mpfr_set(Gees[2][k].re, Gees[3][k].re, MODE); 
        mpfr_set(Gees[2][k].im, Gees[3][k].im, MODE); 
        mpfr_set(Gees[3][k].re, tmp[k].re, MODE); 
        mpfr_set(Gees[3][k].im, tmp[k].im, MODE); 
      }

      /*
      Ens[0] = Ens[1];
      Ens[1] = Ens[2];
      Ens[2] = Ens[3];	
      Ens[3] = 1.0L/dot(tmp, tmp);
      */
      mpfr_set(Ens[0], Ens[1], MODE);
      mpfr_set(Ens[1], Ens[2], MODE);
      mpfr_set(Ens[2], Ens[3], MODE);
      mpfc_dotpr(&buf, tmp, tmp);
      mpfr_ui_div(Ens[3], 1, buf.re, MODE);
  }
  mpfr_clears(buf.re, buf.im, (mpfr_ptr) NULL);
}
