#include "header.h"


void load_parameters(int argc, char* argv[]) {
  if (argc != 2) {
    printf("Usage %s filename\n", argv[0]);
    exit(1);
  } else read_input(argv[1]);  // < ---- here
}

void read_input(char *fname) {
  FILE *fh = fopen(fname,"r");
  char line[512], name[128], value[128];
  if (fh == NULL) {
    printf("Configuration file missing\n");
    exit(1);
  } else {
    while (fgets(line, 512, fh)!=NULL) {
      sscanf(line, "%s\t%s", name, value);
      if (strcmp(name,"precisn=") == 0) state.precision = (mpfr_prec_t) atol(value);  // MUST be set before MPFR constants
      if (strcmp(name,"resname=") == 0) sprintf(state.restart_name,"%s", value);
      if (strcmp(name,"txtasci=") == 0) sprintf(state.txt_format,"%s", value);
      if (strcmp(name,"numbits=") == 0) state.nbits = strtol(value, NULL, 10);
      if (strcmp(name,"gravity=") == 0) {
        //state.gravity = strtold(value, NULL);
        mpfr_init2(state.gravity, state.precision);
        mpfr_set_str (state.gravity, value, 10, MODE);
      }
      if (strcmp(name,"surface=") == 0) {
        //state.surface_tension = strtold(value, NULL);
        mpfr_init2(state.surface_tension, state.precision);
        mpfr_set_str (state.surface_tension, value, 10, MODE);
      }
      if (strcmp(name,"transfl=") == 0) {
        //conf.scaling = strtold(value, NULL);
        mpfr_init2(conf.scaling, state.precision);
        mpfr_set_str (conf.scaling, value, 10, MODE);
      }
      if (strcmp(name,"transfu=") == 0) {
        //conf.image_offset = strtold(value, NULL);
        mpfr_init2(conf.image_offset, state.precision);
        mpfr_set_str (conf.image_offset, value, 10, MODE);
      }
      if (strcmp(name,"toleran=") == 0) {
        //state.tolerance = strtold (value, NULL);
        mpfr_init2(state.tolerance, state.precision);
        mpfr_set_str (state.tolerance, value, 10, MODE);
      }
      if (strcmp(name,"timesim=") == 0) {
        //state.final_time = strtold (value, NULL);
        mpfr_init2(state.final_time, state.precision);
        mpfr_set_str (state.final_time, value, 10, MODE);
      }
      if (strcmp(name,"n_poles=") == 0) {
        state.number_poles = atol(value);
      }
    }
    mpfr_inits2(state.precision, state.mean_level, state.time, (mpfr_ptr) NULL);
    mpfr_inits2(state.precision, state.potentialE, state.kineticE, (mpfr_ptr) NULL);
    mpfr_inits2(state.precision, state.momentum.re, state.momentum.im, (mpfr_ptr) NULL);
    mpfft_init(state.precision); // init mpfft routines
  }
}

void set_initial_data() {
  mpfr_t overN, q, u;
  mpfr_t Q, C, a1, a2; 			// pirate run 2.5L
  mpfr_t buf1, buf2, buf3, buf4; 	// tmp arrays
  mpfc_t arg1, arg2;			// tmp complex arrays
  
  mpfr_inits2(state.precision, overN, q, u, Q, C, a1, a2, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, buf1, buf2, buf3, buf4, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, arg1.re, arg1.im, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, arg2.re, arg2.im, (mpfr_ptr) NULL);
  // ------------------------ 
  mpfr_set_ui(overN, 1, MODE);
  mpfr_div_ui(overN, overN, 1<<state.nbits, MODE);  	// overN
  mpfr_set_d(Q,  2.50, MODE);  				// Q
  mpfr_set_d(C, -0.02, MODE);				// C
  mpfr_set_d(a1, 0.0050, MODE);				// a1
  mpfr_set_d(a2, 0.0075, MODE);				// a2

  //long double overN = 1.L/state.number_modes;
  for (long int j = 0; j < 1<<state.nbits; j++) {
    //q = 2.L*PI*(j*overN - 0.5L) - conf.origin_offset;
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

    //data[0][j] =  clogl(1.IL*csinl(0.5L*(u-1.IL*a1))) -  clogl(1.IL*csinl(0.5L*(u-1.IL*a2)));  // sim 11: pirate run
    // sim 11 initial data:
    // sim 11 (first term)
    mpfr_div_si(arg1.re,  u,  2, MODE);
    mpfr_div_si(arg1.im, a1, -2, MODE);  // i sin z =1/2( e^iz - e^-iz ) = -cos Re sinh Im + i sin Re cosh Im
   
    mpfr_sin_cos(buf3, buf4, arg1.re, MODE);
    mpfr_sinh_cosh(buf1, buf2, arg1.im, MODE);
    // complex i sin z1
    mpfr_mul(arg1.re, buf1, buf4, MODE);
    mpfr_neg(arg1.re, arg1.re, MODE);
    mpfr_mul(arg1.im, buf2, buf3, MODE);
    // complex log z1
    mpfr_hypot(data[0][j].re, arg1.re, arg1.im, MODE);   
    mpfr_log(data[0][j].re, data[0][j].re, MODE);
    mpfr_atan2(data[0][j].im, arg1.im, arg1.re, MODE);
  
    // sim 11 (second term)
    
    mpfr_div_si(arg1.re,  u,  2, MODE);
    mpfr_div_si(arg1.im, a2, -2, MODE);  // i sin z =1/2( e^iz - e^-iz ) = -cos Re sinh Im + i sin Re cosh Im
   
    mpfr_sin_cos(buf3, buf4, arg1.re, MODE);
    mpfr_sinh_cosh(buf1, buf2, arg1.im, MODE);
    // complex i sin z2
    mpfr_mul(arg1.re, buf1, buf4, MODE);
    mpfr_neg(arg1.re, arg1.re, MODE);
    mpfr_mul(arg1.im, buf2, buf3, MODE);
    // complex log z2
    mpfr_hypot(arg2.re, arg1.re, arg1.im, MODE);   
    mpfr_log(arg2.re, arg2.re, MODE);
    mpfr_atan2(arg2.im, arg1.im, arg1.re, MODE);
   
    mpfr_sub(data[0][j].re, data[0][j].re, arg2.re, MODE);
    mpfr_sub(data[0][j].im, data[0][j].im, arg2.im, MODE);
    
    // end sim 11 initial data for Z      
    //data[0][j] = -1.IL*Q*data[0][j]/1.L; 	// necessary (sim 11)
    mpfr_swap(data[0][j].re, data[0][j].im);
    mpfr_neg(data[0][j].im, data[0][j].im, MODE);
    mpfr_mul(data[0][j].re, data[0][j].re, Q, MODE);
    mpfr_mul(data[0][j].im, data[0][j].im, Q, MODE);
  }
  // finished setting Z.  -- checked and validated.

  complex_array_out("Z0.txt", data[0]);
  convertZtoQ(data[0], data[0]);  // convert Z-tilde to Q
  complex_array_out("Q0.txt", data[0]); // checks out good
  
  for (unsigned int j = 0; j < 1<<state.nbits; j++) {
    //data[1][j] = C*(data[0][j]*data[0][j] - 1.L);  // required when setting Z-tilde
    mpfr_mul(arg1.re, data[0][j].im, data[0][j].im, MODE);
    mpfr_fms(arg1.re, data[0][j].re, data[0][j].re, arg1.re, MODE);
    mpfr_mul(arg1.im, data[0][j].re, data[0][j].im, MODE);
    mpfr_fma(arg1.im, data[0][j].im, data[0][j].re, arg1.im, MODE);
    mpfr_sub_ui(data[1][j].re, arg1.re, 1, MODE);
    mpfr_mul(data[1][j].re, data[1][j].re, C, MODE);
    mpfr_mul(data[1][j].im, arg1.im, C, MODE);
    //tmpc[0][j] = data[1][j]*overN;
    mpfr_mul(tmpc[0][j].re, data[1][j].re, overN, MODE);
    mpfr_mul(tmpc[0][j].im, data[1][j].im, overN, MODE);
  }
  complex_array_out("V0.txt", data[1]);
  
  //fftwl_complex z0;
  //fftwl_execute(ift0);
  mpfft_execute(ift0);
  mpfr_hypot(buf1, tmpc[0][0].re, tmpc[0][0].im, MODE);
  mpfr_printf("Average V = %.12Re\n", buf1);
  
  mpfc_t zero, z0;
  
  mpfr_inits2(state.precision, z0.re, z0.im, (mpfr_ptr) NULL);
  mpfr_inits2(state.precision, zero.re, zero.im, (mpfr_ptr) NULL);
  mpfr_set_ui(zero.re, 0, MODE);
  mpfr_set_ui(zero.im, 0, MODE);
  compute_zero_mode_complex(tmpc[0], zero, &z0);  // (OK, verified)
  
  //tmpc[0][0] = z0;
  mpfr_set(tmpc[0][0].re, z0.re, MODE);
  mpfr_set(tmpc[0][0].im, z0.im, MODE);
  mpfr_hypot(buf1, tmpc[0][0].re, tmpc[0][0].im, MODE);

  mpfr_printf("Average V = %.12Re\n", buf1);
  mpfft_execute(ft0);

  for (long int j = 0; j < 1<<state.nbits; j++) {
    mpfr_set(data[1][j].re, tmpc[0][j].re, MODE);
    mpfr_set(data[1][j].im, tmpc[0][j].im, MODE);
  }
  /*
  for (long int j = 0; j < 1<<state.nbits; j++) {
    tmpc[1][j] = data[0][j]*overN;
  }
  fftwl_execute(ift1);
  */
  /*
  for (long int j = 0; j < state.number_modes; j++) {
    if (cabsl(tmpc[0][j]) < 1.0E-14L) tmpc[0][j] = 0.L;
    if (cabsl(tmpc[1][j]) < 1.0E-14L) tmpc[1][j] = 0.L;
  }
  */
  /*
  fftwl_execute(ft0);
  fftwl_execute(ft1);
  memcpy(data[1], tmpc[0], state.number_modes*sizeof(fftwl_complex));
  memcpy(data[0], tmpc[1], state.number_modes*sizeof(fftwl_complex));
  */
  complex_array_out("Q1.txt", data[0]);
  complex_array_out("V1.txt", data[1]);
  // clear all
  mpfr_clears(overN, q, u, Q, C, a1, a2, (mpfr_ptr) NULL);
  mpfr_clears(buf1, buf2, buf3, buf4, (mpfr_ptr) NULL);
  mpfr_clears(arg1.re, arg1.im, (mpfr_ptr) NULL);
  mpfr_clears(arg2.re, arg2.im, (mpfr_ptr) NULL);
  mpfr_clears(z0.re, z0.im, (mpfr_ptr) NULL);
  mpfr_clears(zero.re, zero.im, (mpfr_ptr) NULL);
}


/*
void set_initial_JW() {
  FILE *fh = fopen(state.restart_name,"r");
  if (fh) {
    printf("Restart file opened successfully.\n");
    char line[512], *v[5];
    for (int j = 0; j < 5; j++) v[j] = malloc(512);
    if (fgets(line, 512, fh) != NULL);
    if (fgets(line, 512, fh) != NULL);
    if (fgets(line, 512, fh) != NULL) sscanf(line, "# %s\t%s\t%s\t%s\t%s", v[0], v[1], v[2], v[3], v[4]);
    if (fgets(line, 512, fh) != NULL);
    conf.scaling = 1.0L;
    conf.image_offset = 0.0L;
    long int counter = 0;
    memset(data[0], 0, state.number_modes*sizeof(fftwl_complex));
    memset(data[1], 0, state.number_modes*sizeof(fftwl_complex));
    data[0][0] = 1.0IL*strtold(v[0], NULL);
    data[1][0] = 0.0L*strtold(v[3], NULL);
    while (fgets(line, 512, fh) != NULL) {
      sscanf(line, "%s\t%s\t%s\t%s\n", v[0], v[1], v[2], v[3]);
      data[0][counter+1] = 2.IL*(strtold(v[0], NULL) + 1.IL*strtold(v[1], NULL));
      data[1][counter+1] = 2.0L*(strtold(v[2], NULL) - 1.IL*strtold(v[3], NULL));
      data[0][state.number_modes-counter-1] = 0.L;
      data[1][state.number_modes-counter-1] = 0.L;
      if (cabsl(data[0][counter+1]) < 1e-14L) data[0][counter+1] = 0.L;
      if (cabsl(data[1][counter+1]) < 1e-14L) data[1][counter+1] = 0.L;
      counter++;
    }
    if (counter > state.number_modes/2 - 1) {
      printf("Number of Fourier in file %ld\n", 2*(counter+1));
      printf("Requested in simulation is %ld\n", state.number_modes);
    }
    fclose(fh);
  } else {
    printf("Restart file not found.\n");
  }
  memcpy(tmpc[0], data[0], state.number_modes*sizeof(fftwl_complex));
  fftwl_execute(ft0);
  memcpy(data[0], tmpc[0], state.number_modes*sizeof(fftwl_complex));
  convertZtoQ(data[0], data[0]);
  for (long int j = 0; j < state.number_modes/2; j++) {
    tmpc[1][j] = -1.IL*j*data[1][j];
  }  
  fftwl_execute(ft1);
  for (long int j = 0; j < state.number_modes; j++) {
    data[1][j] = 1.IL*tmpc[1][j]*(data[0][j]*data[0][j]);
  }
  complex_array_out("Q1.txt", data[0]);
  complex_array_out("V1.txt", data[1]);
}
*/


/*
void load_ascii() {
  FILE *fh = fopen(state.restart_name, "r");
  if (fh) {
    char line[512], *v[5];
    long int		N;
    mpfr_t		T;
    for (int j = 0; j < 5; j++) v[j] = malloc(512);
    if (fgets(line, 512, fh) != NULL);
    if (fgets(line, 512, fh) != NULL) sscanf(line, "# N = %s\tL = %s\tu = %s\tT = %s", v[0], v[1], v[2], v[3]);
    if (fgets(line, 512, fh) != NULL);

    N = strtol(v[0], NULL, 10);
    T = strtold(v[3], NULL);

    conf.scaling = strtold(v[1], NULL);
    conf.image_offset = strtold(v[2], NULL);

    printf("Restart:\ntime = %.19LE\nN modes = %ld\n", T, N); 
    printf("Conformal L = %.19LE\nConformal u = %.19LE\n", conf.scaling, conf.image_offset);
    if (N != state.number_modes) {
      printf("Incompatible Grids\nPlaceholder\n");
      exit(1);
    }
    int counter = 0;
    while (fgets(line, 512, fh) != NULL) {
      if (counter == state.number_modes) {
	printf("Something is wrong with restart file\n"); 
	exit(1);
      }
      sscanf(line, "%s\t%s\t%s\t%s\t%s", v[0], v[1], v[2], v[3], v[4]);
      data[0][counter] = strtold(v[1], NULL)-strtold(v[0], NULL) + 1.0IL*strtold(v[2],NULL);
      data[1][counter] = strtold(v[3], NULL) + 1.0IL*strtold(v[4],NULL);
      counter++;
    }
    fclose(fh);    
    if (counter != state.number_modes) {
	printf("Something is wrong with restart file\n");
        exit(1);
    }
    for (int j = 0; j < 5; j++) free(v[j]);
  } else {
    printf("Missing restart file\n");
    exit(1);
  }
}
*/
