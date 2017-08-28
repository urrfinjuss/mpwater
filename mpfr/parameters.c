#include "header.h"

static params_ptr 	 input;
static params  		 rs_input;
static aux_ptr		 extra;
static data_ptr		 array;

void init_parameters(params_ptr ginput, aux_ptr gextra, data_ptr garray) {
  input = ginput; extra = gextra; array = garray;
}

void save_data() {
  printf("Writing binary restart to:\n%s\n", input->rname);
  FILE *fh = fopen(input->rname, "wb");
  fwrite(input, sizeof(params), 1, fh);
  fwrite(array->Q, sizeof(fftwl_complex), input->N, fh);
  fwrite(array->V, sizeof(fftwl_complex), input->N, fh);
  fclose(fh);
}

void load_data() {
  debug_msg("Load Data: initiated\n", EXIT_FALSE);
  FILE *fh = fopen(input->rname, "rb");
  int rflag;  

  if (fh) {
	if (input->d) rflag = 1;	// read from Pade data
	if (input->ascii) rflag = 3;	// read from ascii restart
	else rflag = 0;			// read from binary restart
  } else rflag = 2;			// start with default IC

  printf("rflag = %d\n", rflag);
  switch (rflag) {

  	case (0):
		if (rs_input.N == input->N) {
			//printf("%ld\t%ld\n", nr, rs_input.N);
			debug_msg("Load Data: calling initialize data\n", EXIT_FALSE);
		        printf("%ld\t%ld\n", rs_input.N, input->N);
			initialize_data();
			initialize_auxiliary_arrays();
			debug_msg("Initialize Data: complete\n", EXIT_FALSE);
			if (!fread(array->Q, sizeof(fftwl_complex), input->N, fh)) printf("Broken Q array\n");
			if (!fread(array->V, sizeof(fftwl_complex), input->N, fh)) printf("Broken V array\n");
			fclose(fh);
			//move_mesh(array.Q, array.V, rs_input.u, rs_input.L);
			printf("Default Mesh Set From Restart (L,u) = (%.10LE,%.10Le)\n", rs_input.L, rs_input.u);
			input->L = rs_input.L;  
			input->u = rs_input.u;
			debug_msg("Load Data: restart arrays read successful\n", EXIT_FALSE);
			debug_msg("Load Data: complete\n", EXIT_FALSE);
			set_aux();
		} else {
			printf("%ld\t%ld\n", rs_input.N, input->N);
			fclose(fh);
			debug_msg("Load Data: restart arrays mismatch\n", EXIT_FALSE);
		        debug_msg("Load Data: arrays sized to match input data\n", EXIT_FALSE);
			debug_msg("Load Data: calling initialize data\n", EXIT_FALSE);
			input->N = rs_input.N;
			input->L = rs_input.L;
			input->u = rs_input.u;
			initialize_data();
			if (!fread(array->Q, sizeof(fftwl_complex), input->N, fh)) printf("Broken Q array\n");
			if (!fread(array->V, sizeof(fftwl_complex), input->N, fh)) printf("Broken V array\n");
			debug_msg("Load Data: complete\n", EXIT_FALSE);
			initialize_auxiliary_arrays();
			debug_msg("Initialize Data: complete\n", EXIT_FALSE);
			set_aux();
		}
		break;

	case(1):
		debug_msg("Load Data: calling initialize data\n", EXIT_FALSE);
		printf("%ld\t%ld\n", rs_input.N, input->N);
		initialize_data();
		initialize_auxiliary_arrays();
		debug_msg("Initialize Data: complete\n", EXIT_FALSE);      
		fclose(fh);
		//move_mesh(array.Q, array.V, rs_input.u, rs_input.L);
		input->u = 0.Q;
		//input->L = 1.Q;
		printf("Default Mesh Set From Restart (L,u) = (%.10LE,%.10Le)\n", input->L, input->u);

		debug_msg("Load Data: restart arrays read successful\n", EXIT_FALSE);
		debug_msg("Load Data: complete\n", EXIT_FALSE);    
		printf("Name %s\nN = %ld\ng = %LE\ns = %LE\nl = %LE\nq = %LE\ntl = %LE\n", 
       			input->rname, input->N, input->g, input->s, input->L, input->u, input->tl);
		//set_aux();
		read_pade();
		break;

	case(2):
		debug_msg("Load Data: restart file missing\nLoad Data: complete\n", EXIT_FALSE);
		printf("Restart File Missing:\tstarting new simulation.\n");
		initialize_data();
		initialize_auxiliary_arrays();
		debug_msg("Initialize Data: complete\n", EXIT_FALSE);
		debug_msg("Not using restart data in this run\n", EXIT_FALSE);
		set_aux();
   		long double u, q;
   
		for (int j = 0; j < input->N; j++) {
     		     q = PI*(2.L*j/input->N - 1.0L);
		     	 u = input->u + 2.L*atan2l(input->L*sinl(0.5L*(q-extra->q)), cosl(0.5L*(q-extra->q)));
     		     array->Q[j] = 1.L; //1E-13*cexpl(-1.IL*u);
     		     array->V[j] = -0.025IL*(1.L/ctanl(0.5L*(u-0.12IL)) - 1.IL) + 0.00IL*(1.L/ctanl(0.5L*(u-0.24IL)) - 1.IL); // run 8&9&10&11
   		}
		break;
  }

}

void dump_input() {
  FILE *fh = fopen("output.log","a");
  fprintf(fh, "Load Parameters: restart name = %s\n", input->rname);
  fprintf(fh, "Load Parameters: number of points = %ld\n", input->N);
  fprintf(fh, "Load Parameters: gravity g = %.15Le\n", input->g);
  fprintf(fh, "Load Parameters: surface tension s = %.15Le\n", input->s);
  fprintf(fh, "Load Parameters: transformation l = %.15Le\n", input->L);
  fprintf(fh, "Load Parameters: transformation u0 = %.15Le\n", input->u);
  fprintf(fh, "Load Parameters: refinement tolerance tol = %.15Le\n", input->tl);
  fprintf(fh, "Load Parameters: complete\n");
  fclose(fh);
}

void read_input(char *fname) {
  FILE *fh = fopen(fname,"r");
  char line[512], name[128], value[128];
  if (fh == NULL) {
    debug_msg("Load Parameters: input file missing\n", EXIT_FALSE);
    debug_msg("Load Parameters: complete\n", EXIT_TRUE);
  } else {
    while (fgets(line, 512, fh)!=NULL) {
      sscanf(line, "%s\t%s", name, value);
      if (strcmp(name,"resname=") == 0) sprintf(input->rname,"%s", value);
      if (strcmp(name,"npoints=") == 0) input->N = strtol(value, NULL, 10);
      if (strcmp(name,"gravity=") == 0) input->g = strtold(value, NULL);
      if (strcmp(name,"surface=") == 0) input->s = strtold(value, NULL);
      if (strcmp(name,"transfl=") == 0) input->L = strtold(value, NULL);
      if (strcmp(name,"transfu=") == 0) input->u = strtold(value, NULL);
      if (strcmp(name,"toleran=") == 0) input->tl = strtold (value, NULL);
      if (strcmp(name,"n_poles=") == 0) input->d = atol(value);
    }
    dump_input(); 
  }
}

void load_parameters(int argc, char* argv[]) {
  FILE *fh = fopen("output.log","w");
  fprintf(fh, "Output Log File:\n\n");
  fprintf(fh, "Load Parameters: initiated\n");
  fclose(fh);
  if (argc == 1) {
    fh = fopen("output.log","a");
    fprintf(fh, "Load Parameters: requires name of file with parameters as a single input argument\n");
    fprintf(fh, "Load Parameters: complete\n");
    fclose(fh);
    exit(1);
  } else {
    char str[256];
    sprintf(str,"%s", argv[1]);
    fh = fopen("output.log","a");
    fprintf(fh, "Load Parameters: reading parameters from %s\n", str);
    fclose(fh);
    read_input(str);
  }
}

