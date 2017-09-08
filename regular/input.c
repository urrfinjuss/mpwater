#include "header.h"


void load_parameters(int argc, char* argv[]) {
  if (argc != 2) {
    printf("Usage %s filename\n", argv[0]);
    exit(1);
  } else read_input(argv[1]);
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
      if (strcmp(name,"resname=") == 0) sprintf(state.restart_name,"%s", value);
      if (strcmp(name,"txtasci=") == 0) sprintf(state.txt_format,"%s", value);
      if (strcmp(name,"npoints=") == 0) state.number_modes = strtol(value, NULL, 10);
      if (strcmp(name,"gravity=") == 0) state.gravity = strtold(value, NULL);
      if (strcmp(name,"surface=") == 0) state.surface_tension = strtold(value, NULL);
      if (strcmp(name,"transfl=") == 0) conf.scaling = strtold(value, NULL);
      if (strcmp(name,"transfu=") == 0) conf.image_offset = strtold(value, NULL);
      if (strcmp(name,"toleran=") == 0) state.tolerance = strtold (value, NULL);
      if (strcmp(name,"timesim=") == 0) state.final_time = strtold (value, NULL);
      if (strcmp(name,"n_poles=") == 0) state.number_poles = atol(value);
    }
  }
}

void set_initial_data() {
  long double overN = 1.L/state.number_modes;
  long double q, u;
  long double fi = 0.0L*PI;
  //long double C = -0.02L;  
  	fftwl_complex C  = -0.03L + 0.02IL;
	fftwl_complex a1 =  0.0000L + 0.0040IL;
	fftwl_complex a2 =  0.0160L + 0.0200IL;
	fftwl_complex Q  =  0.0500L*cexpl(0.21IL*PI);

  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.L*PI*(j*overN - 0.5L) - conf.origin_offset;
    u = conf.image_offset - fi + 2.L*atan2l(conf.scaling*sinl(0.5L*q),cosl(0.5L*q));
    //data[0][j] =  clogl(  1.IL*csinl(0.5L*(u-1.IL*a1))) - clogl(1.IL*csinl(0.5L*(u-1.IL*a2)));;  // sim 14 negative discriminant
    data[0][j] =  clogl(1.IL*csinl(0.5L*(u-a1))) -  clogl(1.IL*csinl(0.5L*(u-a2)));  // sim 11: pirate run
    data[0][j] = -1.IL*Q*data[0][j]/1.L; 	// necessary (sim 11)
    // pade test for VZi
    tmpc[0][j] = data[1][j]*overN; 
  }
  // if setting Z-tilde then uncomment below
  
  //complex_array_out("Z0.txt", data[0]);
  convertZtoQ(data[0], data[0]);  // convert Z-tilde to Q
  
  // end uncomment
  //complex_array_out("Q0.txt", data[0]);
  for (long int j = 0; j < state.number_modes; j++) {
    data[1][j] = C*(data[0][j]*data[0][j] - 1.L);  // required when setting Z-tilde
    tmpc[0][j] = data[1][j]*overN;
  }
  //complex_array_out("V0.txt", data[1]);
  fftwl_complex z0;
  fftwl_execute(ift0);
  //printf("Average V = %.12LE\n", cabsl(tmpc[0][0]));
  if (PADE_TEST) {
    for (long int j = 0; j < state.number_modes; j++) {
      data[0][j] = 1.L/(data[0][j]*data[0][j]);
      data[1][j] = -1.IL*data[1][j]*data[0][j];
    }
    optimal_pade("Zu_initial.pade", data[0]); 
    optimal_pade("Phiu_initial.pade", data[1]); 
    printf("PADE_TEST flag is on!\n");
    exit(1);
  }

  compute_zero_mode_complex(tmpc[0], 0.L, &z0);
  tmpc[0][0] = z0;
  //printf("Average V = %.12LE\n", cabsl(tmpc[0][0]));
  for (long int j = 0; j < state.number_modes; j++) tmpc[1][j] = data[0][j]*overN;

  fftwl_execute(ift1);
  for (long int j = 0; j < state.number_modes; j++) {
    if (cabsl(tmpc[0][j]) < 1.0E-18L) tmpc[0][j] = 0.L;
    if (cabsl(tmpc[1][j]) < 1.0E-18L) tmpc[1][j] = 0.L;
  }
  tmpc[0][state.number_modes/2-1] = 0.L;
  tmpc[1][state.number_modes/2-1] = 0.L;
  
  fftwl_execute(ft1);
  fftwl_execute(ft0);
  memcpy(data[1], tmpc[0], state.number_modes*sizeof(fftwl_complex));
  memcpy(data[0], tmpc[1], state.number_modes*sizeof(fftwl_complex));
  //complex_array_out("Q1.txt", data[0]);
  //complex_array_out("V1.txt", data[1]);
}

void set_initial_JW() {
  //long double overN = 1.L/state.number_modes;
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
      //data[0][counter+1] = data[0][counter+1]*overN;
      //data[1][counter+1] = data[1][counter+1]*overN;
      counter++;
      //printf("%s\t%s\t%s\t%s\n", v[0], v[1], v[2], v[3]);
      //printf("%LE\t%LE\t", creall(data[0][counter]), cimagl(data[0][counter]));
      //printf("%LE\t%LE\n", creall(data[1][counter]), cimagl(data[1][counter]));
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
  //memcpy(tmpc[1], data[1], state.number_modes*sizeof(fftwl_complex));
  fftwl_execute(ft0);
  //fftwl_execute(ft1);
  memcpy(data[0], tmpc[0], state.number_modes*sizeof(fftwl_complex));
  //memcpy(data[1], tmpc[1], state.number_modes*sizeof(fftwl_complex));
  /*
  complex_array_out("jon_z.txt",   tmpc[0]); 
  complex_array_out("jon_psi.txt", tmpc[1]);
  output_data("all_jon", tmpc[1]); 
  exit(1);
  */
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
  //exit(1);
}


void load_ascii() {
  FILE *fh = fopen(state.restart_name, "r");
  if (fh) {
    char line[512], *v[5];
    long int		N;
    long double		T;
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

void load_pade() {
  FILE *fh = fopen(state.restart_name, "r");
  long double		speed;
  if (fh) {
    char line[512], *v[5];
    long double		T, y0;
    long double 	overN = 1.L/state.number_modes;

    for (int j = 0; j < 5; j++) v[j] = malloc(512);
    if (fgets(line, 512, fh) != NULL);
    //# M = 16384     Residual = 7.666871e-32 y0 = -2.93733978518518156291891533016033e-01
    //# Amplitude at highest Fourier mode = 3.00000e-37	Pade error = 5.09180e-31        H/lambda = 1.31331355155599360515621709076606e-01	c = 1.0870000000000000e+00

    if (fgets(line, 512, fh) != NULL) sscanf(line, "# M = %s\tResidual = %s\ty0 = %s\n", v[0], v[1], v[2]);
    y0 = strtold(v[2], NULL);
    if (fgets(line, 512, fh) != NULL) sscanf(line, "# Amplitude at highest Fourier mode = %s\tPade error = %s\tH/lambda = %s\tc = %s\n", v[0], v[1], v[2], v[3]);
    speed = strtold(v[3], NULL);
    //printf("Wave Speed is %.12Le\n", speed);
    if (fgets(line, 512, fh) != NULL);

    //N = strtol(v[0], NULL, 10);
    //T = strtold(v[3], NULL);
    T = 0.L;

    //conf.scaling = strtold(v[1], NULL);
    //conf.image_offset = strtold(v[2], NULL);

    //printf("Restart:\ntime = %.19LE\nN modes = %ld\n", T, state.number_modes); 
    //printf("Conformal L = %.19LE\nConformal u = %.19LE\n", conf.scaling, conf.image_offset);
    int counter = 0;
    for (int j = 0; j < state.number_modes; j++) {
      data[0][j] = 1.IL*y0;
    }
    while (fgets(line, 512, fh) != NULL) {
      long double Xn, Gn;
      long double u, q;
      sscanf(line, "%s\t%s\t%s", v[0], v[1], v[2]);
      //printf("%s", line);
      Xn = strtold(v[0], NULL);
      Gn = strtold(v[1], NULL);
      for (int j = 0; j < state.number_modes; j++) {
        q = 2.L*PI*(j*overN - 0.5L) - conf.origin_offset;
        u = conf.image_offset + 2.L*atan2l(conf.scaling*sinl(0.5L*q),cosl(0.5L*q));
        data[0][j] += Gn/(tanl(0.5L*u) - 1.IL*Xn);
      }
      //data[0][counter] = strtold(v[1], NULL)-strtold(v[0], NULL) + 1.0IL*strtold(v[2],NULL);
      //data[1][counter] = strtold(v[3], NULL) + 1.0IL*strtold(v[4],NULL);
      counter++;
    }
    fclose(fh);    
    for (int j = 0; j < 5; j++) free(v[j]);
  } else {
    printf("Missing restart file\n");
    exit(1);
  }
  //surface_out("initial_surface.txt",data[0]);
  convertZtoQ(data[0], data[0]);
  for (int j = 0; j < state.number_modes; j++) {
    data[1][j] = 1.IL*speed*(1.L - data[0][j]*data[0][j]);
  }
}
