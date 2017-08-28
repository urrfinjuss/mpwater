#include "header.h"
/*
void real_array_out(char *fname, long double *in) {
  FILE *fh = fopen(fname,"w");
  long double u, q, overN = 1.L/state.number_modes;

  fprintf(fh, "# 1. q 2. u 3. Array\n");
  fprintf(fh, "# Time = %.14LE\tL = %.14LE\n\n", state.time, conf.scaling);
  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.0L*PI*(j*overN - 0.5L) - conf.origin_offset;
    u = conf.image_offset + 2.L*atan2l(conf.scaling*sinl(0.5L*q), cosl(0.5L*q));
    fprintf(fh, "%24.17LE\t%24.17LE\t%24.17LE\n", q+conf.origin_offset, u, in[j]);
  }
  fclose(fh);
}
*/


void complex_array_out(char *fname, mpfc_t *in) {
  FILE *fh = fopen(fname,"w");
  mpfr_t u, q, overN;
  mpfr_t buf1, buf2, buf3;
  
  mpfr_inits2(state.precision, u, q, overN, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, buf1, buf2, buf3, (mpfr_ptr) NULL);
  mpfr_set_ui(overN, 1, MODE);
  mpfr_div_ui(overN, overN, 1<<state.nbits, MODE);

  mpfr_fprintf(fh, "# 1. q 2. u 3.-4. Array\n");
  mpfr_fprintf(fh, "# Time = %.14Re\tL = %.14Re\n\n", state.time, conf.scaling);
  for (long int j = 0; j < 1<<state.nbits; j++) {
    //q = 2.0L*PI*(j*overN - 0.5L) - conf.origin_offset;
    mpfr_mul_ui(q, overN, 2*j, MODE); 
    mpfr_sub_ui(q, q, 1, MODE);
    mpfr_mul(q, q, Pie, MODE);
    mpfr_sub(q, q, conf.origin_offset, MODE);
    //u = conf.image_offset + 2.L*atan2l(conf.scaling*sinl(0.5L*q),cosl(0.5L*q));
    mpfr_div_ui(buf1, q, 2, MODE);
    mpfr_sin_cos(buf2, buf3, buf1, MODE);
    mpfr_mul(buf2, buf2, conf.scaling, MODE);
    mpfr_atan2(u, buf2, buf3, MODE);
    mpfr_mul_ui(u, u, 2, MODE);
    mpfr_add(u, u, conf.image_offset, MODE);

    mpfr_add(q, q, conf.origin_offset, MODE);
    mpfr_fprintf(fh, "%24.17Re\t%24.17Re\t%24.17Re\t%24.17Re\n", q, u, in[j].re, in[j].im);
  }
  fclose(fh);
  mpfr_clears(u, q, overN, (mpfr_ptr) NULL);
  mpfr_clears(buf1, buf2, buf3, (mpfr_ptr) NULL);
}



void spec_out(char *fname, mpfc_t *in1, mpfc_t *in2) {
  FILE 		*fh = fopen(fname,"w");
  mpfr_t 	overN, buf1, buf2;
  mpfr_inits2(state.precision, overN, buf1, buf2, (mpfr_ptr) NULL);
  mpfr_set_ui(overN, 1, MODE);
  mpfr_div_ui(overN, overN, 1<<state.nbits, MODE);

  mpfr_fprintf(fh, "# 1. k 2. |a_k| 3. |b_k|\n");
  mpfr_fprintf(fh, "# Time = %.14Re\tL = %.14Re\n\n", state.time, conf.scaling);
  for (long int j = 0; j < 1<<state.nbits; j++) {
    mpfr_hypot(buf1, in1[j].re, in1[j].im, MODE);
    mpfr_mul(buf1, buf1, overN, MODE);
    mpfr_hypot(buf2, in2[j].re, in2[j].im, MODE);
    mpfr_mul(buf2, buf2, overN, MODE);

    mpfr_fprintf(fh, "%ld\t%23.17Re\t%23.17Re\n", j, buf1, buf2);
  }
  fclose(fh);
  mpfr_clears(overN, buf1, buf2, (mpfr_ptr) NULL);
}

/*
void output_data(char *fname, fftwl_complex *inPhi) {
  FILE *fh = fopen(fname,"w");
  long double u, q, overN = 1.L/state.number_modes;

  //convertQtoZ(data[0], tmpc[3]);
  fprintf(fh, "# 1. q 2. u 3.-4. Q 5.-6. V 7.-8. Z 9.-10. Phi\n");
  fprintf(fh, "# Time = %.14LE\tL = %.14LE\n\n", state.time, conf.scaling);
  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.0L*PI*(j*overN - 0.5L) - conf.origin_offset;
    u = conf.image_offset + 2.L*atan2l(conf.scaling*sinl(0.5L*q), cosl(0.5L*q));
    fprintf(fh, "%24.17LE\t%24.17LE\t", q+conf.origin_offset, u);
    fprintf(fh, "%24.17LE\t%24.17LE\t", creall(data[0][j]), cimagl(data[0][j]));
    fprintf(fh, "%24.17LE\t%24.17LE\t", creall(data[1][j]), cimagl(data[1][j]));
    fprintf(fh, "%24.17LE\t%24.17LE\t", creall(tmpc[3][j]), cimagl(tmpc[3][j]));
    fprintf(fh, "%24.17LE\t%24.17LE\n", creall(inPhi[j]), cimagl(inPhi[j]));
  }
  fclose(fh);
}
*/

/*
void surface_out(char *fname, fftwl_complex *in) {
  FILE *fh = fopen(fname,"w");
  long double u, q, overN = 1.L/state.number_modes;

  fprintf(fh, "# 1. x 2. y\n");
  fprintf(fh, "# Time = %.14LE\tL = %.14LE\n\n", state.time, conf.scaling);
  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.0L*PI*(overN*j - 0.5L) - conf.origin_offset;
    u = conf.image_offset + 2.L*atan2l(conf.scaling*sinl(0.5L*q), cosl(0.5L*q));
    fprintf(fh, "%23.17LE\t%23.17LE\n", u + creall(in[j]), cimagl(in[j]));
  }
  fclose(fh);
}
*/

void print_constants() {
  mpfr_printf("#\t\t\t------------------------------------------------------\t\t\t#\n");
  mpfr_printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  mpfr_printf("#\t\t\tConformal Map\t\t\t\t\t\t\t\t#\n");
  mpfr_printf("#\t\t\tParameters:\t\t\t\t\t\t\t\t#\n");
  mpfr_printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  mpfr_printf("#\t\t\tFourier modes N:\t\t%8ld\t\t\t\t#\n", 1<<state.nbits);
  mpfr_printf("#\t\t\tArithmetic precision:\t\t%Pu bits\t\t\t\t#\n", state.precision);
  mpfr_printf("#\t\t\tL-scaling:\t%.16Re\t\t\t\t\t#\n", conf.scaling);
  mpfr_printf("#\t\t\tQ-star:\t\t%18.16Re\t\t\t\t\t#\n", conf.origin_offset);
  mpfr_printf("#\t\t\tU-star:\t\t%18.16Re\t\t\t\t\t#\n", conf.image_offset);
  mpfr_printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  mpfr_printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  mpfr_printf("#\t\t\tConstants:\t\t\t\t\t\t\t\t#\n");
  mpfr_printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  //mpfr_printf("#\t\t\tMomentum X:\t\t%23.16LE\t\t\t\t#\n", cimagl(state.momentum));
  //mpfr_printf("#\t\t\tMomentum Y:\t\t%23.16LE\t\t\t\t#\n", creall(state.momentum));
  //mpfr_printf("#\t\t\tKinetic Energy:\t\t%23.16LE\t\t\t\t#\n",  state.kineticE/PI );
  //mpfr_printf("#\t\t\tPotential Energy:\t%23.16LE\t\t\t\t#\n", state.potentialE/PI);
  //mpfr_printf("#\t\t\tTotal Energy:\t\t%23.16LE\t\t\t\t#\n", (state.potentialE+state.kineticE)/PI);
  mpfr_printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  //mpfr_printf("#\t\t\tMean Level:\t%23.16LE\t\t\t\t\t#\n", state.mean_level);
  mpfr_printf("#\t\t\t------------------------------------------------------\t\t\t#\n");
}


