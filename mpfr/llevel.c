#include "header.h"

/*
static params input;
static data array;
static aux extra;
static consq motion;
static fftwl_complex *tmp0, *tmp1, *tmp2, *tmp3, *qB, *vB, *z, *Rz;
static fftwl_complex q;
static fftwl_plan p0f, p1f, p2f;
static fftwl_plan p0b, p1b;
static long double *k, u;
static long int ref_nog;
static long double ml_tolerance;

static long double diss_fltr;
*/

//-------------------
/*
void init_lowlevel() {
  init_parameters(&input, &extra, &array);
  init_output(&input, &extra);
}

void init_timemarching() {
  ml_tolerance = 1.0E-16L;
#if RK4  
  init_exrk4_module(&input, &extra, &array, &motion);
  diss_fltr = init_rk4();
#elif DIRK4
  init_dirk4_module(&input, &extra, &array, &motion);
  diss_fltr = init_dirk4();
#else
  printf("No Time Marching Selected.\n"); exit(1);
#endif
}
*/

/*
void mult_jacobian(fftwl_complex *rhs, fftwl_complex *x) {
	// multiply by du/dq via tridiagonal operator on the Fourier side
	// rhs is assumed to be Fourier vector (BACKWARD) of coefficients
	// metod progonki == Thomas method
	long double a = 0.5L*(1.0L - powl(input.L,2))/(1.0L + powl(input.L,2));
	
	tmp0[0] = 2.L*rhs[0]/(1.L + powl(input.L,2));
	tmp1[0] = a*cexpl(-1.0IL*extra.q); 
	
	for (long int j = 1; j < input.N/2-1; j++) {
		tmp1[j] = a*cexpl(-1.0IL*extra.q)/(1.L - a*cexpl(1.0IL*extra.q)*tmp1[j-1]);														// tmp1 stores c'
		tmp0[j] = (2.L*rhs[j]/(1.L + powl(input.L,2)) - a*cexpl(1.0IL*extra.q)*tmp0[j-1])/(1.L - a*cexpl(1.0IL*extra.q)*tmp1[j-1]);		// tmp0 stores d'
	}
	
	x[input.N/2-1] = tmp0[input.N/2-1];
	for (long int j = input.N/2-2; j > -1; j--) x[j] = tmp0[j] - tmp1[j]*x[j+1];
}
*/

/*
void inverseZq() {
	for (long int j = 0; j < input.N; j++) {
		tmp0[j] = array.Q[j]*array.Q[j]/input.N;			///extra.dq[j]				// matrix
	}
	memset(tmp1, 0, input.N*sizeof(fftwl_complex));
	tmp1[0] = 1.L;
	fftwl_execute(p0b);
	z[0] = tmp1[0]/tmp0[0]; 
	
	for (long int j = 1; j < input.N/2; j++) {
		fftwl_complex sum = 0.L;
		for (long int l = j-1; l > -1; l--) {
			sum += tmp0[j-l]*z[l];
		}
		z[j] = (tmp1[j]-sum/tmp0[0]);
	}
	memset(z+input.N/2, 0, (input.N/2)*sizeof(fftwl_complex));
	memcpy(tmp1, z, (input.N)*sizeof(fftwl_complex));   // tmp1 now stores z_u
	//primitive_output("z_good.txt", z);
	fftwl_execute(p1f);
	for (long int j = 0; j < input.N; j++) {
		tmp1[j] = 1.IL*(tmp1[j] - 1.L)*extra.newdU[j]/input.N; // bad need tridiagonal solver (need to clarify who is analytic of what..), now stores tilde z_q
	}
	//fftwl_execute(p1b);
	
	
	
	//primitive_output("zq_direct.txt", tmp1);
	//mult_jacobian(z, Rz);	// debugging tridiagonal solver
	//memcpy(tmp1, Rz, (input.N)*sizeof(fftwl_complex));
	//fftwl_execute(p1f);
	//primitive_output("zq_tridiag.txt", Rz);
	
	//for (long int j = 0; j < input.N; j++) tmp1[j] = (1.L/cpowl(array.Q[j],2) - 1.L)*extra.dq[j];
	
	//primitive_output("zu_good.txt",tmp1);
	//for (long int j = 0; j < input.N; j++) tmp1[j] = (1.L/cpowl(array.Q[j],2) - 1.L);
	//primitive_output("zu_bad.txt",tmp1);
	//fftwl_execute(p1b);
	//for (long int j = 0; j < input.N; j++) tmp1[j] = tmp1[j]/input.N;
	//primitive_output("z_bad.txt",tmp1);
	//exit(1);
}
*/

/*
void compute_aux_arrays(fftwl_complex *inQ, fftwl_complex *inV) {
  //  --- verified 
  memcpy(tmp0, inQ, (input.N)*sizeof(fftwl_complex));
  memcpy(tmp1, inV, (input.N)*sizeof(fftwl_complex));
  fftwl_execute(p0f);
  fftwl_execute(p1f);
  for (int j = 0; j < input.N; j++) {
    tmp0[j] = 1.IL*k[j]*tmp0[j]/input.N;
    tmp1[j] = 1.IL*k[j]*tmp1[j]/input.N;
  }
  fftwl_execute(p0b);
  fftwl_execute(p1b);
  memcpy(extra.dQ, tmp0,(input.N)*sizeof(fftwl_complex));
  memcpy(extra.dV, tmp1,(input.N)*sizeof(fftwl_complex));    

  for (int j = 0; j < input.N; j++) {
    tmp1[j] = 2.L*creall(inV[j]*conjl(inQ[j]*inQ[j]));
    tmp0[j] = inV[j]*conjl(inV[j]) + 4.L*input.s*cimagl(tmp0[j]*conjl(inQ[j]))*extra.newdQ[j];  /// sic!
  }
  fftwl_execute(p1b);  
  fftwl_execute(p0b);
  tmp1[input.N/2] = 0.L; extra.b0 = 0.L;
  for (long int j = input.N/2 - 2; j > -1; j--) {
    extra.b0 = extra.b0 + 0.5L*(tmp1[input.N-j-1]*conjl(extra.w[j]) - tmp1[j+1]*extra.w[j] );
  }
  memset(tmp2,0,sizeof(fftwl_complex)*input.N);
  memset(tmp1+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  memset(tmp0+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  tmp1[0] = tmp1[0]/2;
  for (int j = 0; j < (input.N)/2; j++) {
    tmp0[j] = -1.IL*k[j]*tmp0[j]/(input.N);
    tmp1[j] = tmp1[j]/(input.N);
    tmp2[j] = -1.IL*k[j]*tmp1[j];
  }
  fftwl_execute(p2f);  // tmp2 contains dU
  fftwl_execute(p1f);  // tmp1 contains  U
  fftwl_execute(p0f);  // tmp0 contains  B
  for (int j = 0; j < input.N; j++) tmp1[j] = tmp1[j] + extra.b0;
  memcpy( extra.B, tmp0, input.N*sizeof(fftwl_complex));
  memcpy( extra.U, tmp1, input.N*sizeof(fftwl_complex));
  memcpy(extra.dU, tmp2, input.N*sizeof(fftwl_complex));
}
*/

/*
void get_spectrum(char* fname) {
   memcpy(tmp0, array.Q, input.N*sizeof(fftwl_complex));
   memcpy(tmp1, array.V, input.N*sizeof(fftwl_complex)); 
   fftwl_execute(p0b);
   fftwl_execute(p1b);
   spec_output(fname, tmp0, tmp1, motion.time);
}
*/

/*
void get_surface(char* fname) {
  FILE *fh = fopen(fname,"w");
  fprintf(fh, "# 1. u 2. x 3. y 4.-5. V\n# N = %ld\tL = %.19LE\tT = %.19LE\n\n", input.N, input.L, motion.time);
  for (long int j = 0; j < input.N; j++) {
     q = 2.L*PI*(1.L*j/input.N - 0.5L);
     u = input.u + 2.L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
     fprintf(fh, "%.19LE\t%.19LE\t%.19LE\t%.19LE\t%.19LE\n", u, creall(z[j]), cimagl(z[j]), creall(array.V[j]), cimagl(array.V[j]) );
  }
  fclose(fh);
}
*/

/*
void get_integrals() {
  //for (long int j = 0; j < input.N; j++) {
  //tmp0[j] = array.V[j]*extra.newdU[j]/cpowl(array.Q[j], 2.0L)/input.N;
  //tmp1[j] = -1.0IL*extra.newdU[j]*(1.0L - cpowl(array.Q[j], -2.0L))/input.N;   // 
  //}
  //primitive_output("zq_bad.txt",tmp0);
  inverseZq();
  for (long int j = 0; j < input.N; j++) tmp0[j] = array.V[j]*extra.newdU[j]/cpowl(array.Q[j], 2.0L)/input.N;
  //primitive_output("zq_good.txt",tmp0);
  //exit(1);
  fftwl_execute(p0b);
  fftwl_execute(p1b);  
  
  motion.K = 0.L;
  tmp1[0] = 0.L;
  memset(tmp1+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  for (long int j = input.N/2-1; j > 0; j--) {
	motion.K += 0.5L*creall(tmp0[j]*conjl(tmp0[j])/k[j]);
        tmp1[j] = tmp1[j]/k[j];
  }
  fftwl_execute(p1f);  
  for (long int j = 0; j < input.N; j++) {
	tmp1[j] = tmp1[j];  			// zt
	tmp0[j] = cimagl(tmp1[j])/input.N;	
  }
  fftwl_execute(p0f);  
  for (long int j = 1; j < input.N; j++) {
	tmp0[j] = fabsl(k[j])*tmp0[j];   // |k|y
  }
  fftwl_execute(p0b);
  motion.ml = 0.L;
  motion.y0 = 0.L;
  motion.P = 0.L;
  for (long int j = 0; j < input.N; j++) motion.y0 += -cimagl(tmp1[j])*(extra.newdU[j] + tmp0[j])/input.N;	
  for (long int j = 0; j < input.N; j++) {
     q = 2.L*PI*(1.L*j/input.N - 0.5L);
     u = input.u + 2.L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
     tmp2[j] = u + (tmp1[j]-creall(tmp1[0])) + 1.0IL*motion.y0;
     motion.ml += (extra.newdU[j] + tmp0[j])*cimagl(tmp2[j]);
     motion.P += input.g*((extra.newdU[j] + tmp0[j])*cimagl(tmp2[j])*cimagl(tmp2[j]));
  }
  memcpy(z, tmp2, input.N*sizeof(fftwl_complex));
  //FILE *fh = fopen(fname,"w");
  //fprintf(fh, "# 1. u 2. x 3. y\n# N = %ld\tL = %.19LE\tT = %.19LE\n\n", input.N, input.L, motion.time);
  //for (long int j = 0; j < input.N; j++) {
  //   q = 2.L*PI*(1.L*j/input.N - 0.5L);
  //   u = input.u + 2.L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
  //   tmp2[j] = u + (tmp1[j]-creall(tmp1[0])) + 1.0IL*motion.y0;
  //   fprintf(fh, "%.19LE\t%.19LE\t%.19LE\n", u, creall(tmp2[j]), cimagl(tmp2[j]) );
  //   motion.ml += (extra.dq[j] + tmp0[j])*cimagl(tmp2[j]);
  //   motion.P += input.g*((extra.dq[j] + tmp0[j])*cimagl(tmp2[j])*cimagl(tmp2[j]));
  //}
  //fclose(fh);
  motion.ml = motion.ml/input.N;
  motion.P = motion.P/input.N;
  find_max_height(tmp2, &(motion.Xloc), &(motion.H), &(motion.C));
  //printf("Delicate X = %.15LE\nDelicate H = %.15LE\nDelicate C = %.15LE\n", motion.Xloc, motion.H, motion.C);

  motion.Xloc = 0.L;
  motion.H = cimagl(tmp2[input.N/2]);
  motion.C = 2.L*cimagl(tmp2[input.N/2]-tmp2[input.N/2-1])/powl(creall(tmp2[input.N/2] - tmp2[input.N/2-1]),2);
  //printf("Crude    X = %.15LE\nCrude    H = %.15LE\nCrude    C = %.15LE\n", motion.Xloc, motion.H, motion.C);
}
*/

/*
void find_max_height(fftwl_complex *in, long double *x_peak, long double *y_peak, long double *curv) {
  long double y_crude = cimagl(in[0]);
  long int Indx = 0;

  for (long int j = 0; j < input.N; j++) {
    if (cimagl(in[j]) >= y_crude) {
	Indx = j;
	y_crude = cimagl(in[j]);
    }
  }

  long double xx[3], yy[3];
  long double num, den;
  long double A, B, C;

  if (0 < Indx && Indx < input.N-1) {
    for (int j = 0; j < 3; j++){
      xx[j] = creall(in[Indx-1+j]-in[Indx]);
      yy[j] = cimagl(in[Indx-1+j]); 
    }
  } else if (Indx == 0) {
    xx[0] = creall(in[input.N-1]-in[Indx]) - 2.L*PI;  // check later
    yy[0] = cimagl(in[input.N-1]); 
    for (int j = 1; j < 3; j++){
      xx[j] = creall(in[j-1]-in[Indx]);
      yy[j] = cimagl(in[j-1]); 
    }
  } else {
    for (int j = 0; j < 2; j++){
      xx[j] = creall(in[Indx-1+j]-in[Indx]);
      yy[j] = cimagl(in[Indx-1+j]); 
    }
    xx[2] = PI; //creall(in[0]-in[Indx]); check later
    yy[2] = cimagl(in[0]);
  }
  num = (yy[2]-yy[1])*(xx[1]-xx[0]) - (yy[1]-yy[0])*(xx[2]-xx[1]);
  den = (xx[2]*xx[2] - xx[1]*xx[1])*(xx[1]-xx[0]) - (xx[1]*xx[1]-xx[0]*xx[0])*(xx[2]-xx[1]);
  A = num/den;
  B = (yy[2]-yy[0])/(xx[2]-xx[0]) - A*(xx[2] + xx[0]);
  C = yy[1] - A*xx[1]*xx[1] - B*xx[1];


  *x_peak = -0.5L*B/A;
  *y_peak = -C; //A*powl(0.5L*B/A,2) + B*(-0.5L*B/A) + C;
  *curv = 2.L*A;
}
*/

/*
long double track_singularity(){
  long double Min = 1.L, tmpl = 1.L, Umin = 0.L, u;

  for (long int j = 0; j < input.N; j++) {
    tmpl = creall(array.Q[j]*conjl(array.Q[j]));
    q = 2.L*PI*(1.L*j/input.N - 0.5L);
    u = input.u + 2.L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
    if (tmpl < Min) {
	Min = tmpl;
        Umin = u;
    }
  }
#if MOVE_MESH
  return Umin;
#else
  return 0.L*Umin;
#endif
}
*/

/*
fftwl_complex check_resolution() {
  long double EQ0 = 0.L, EQ1 = 0.L, EQ2 = 0.L, EQ4 = 0.L, EQ8 = 0.L;
  long double EV0 = 0.L, EV1 = 0.L, EV2 = 0.L, EV4 = 0.L, EV8 = 0.L; 
  long double rq3;   //rq0, rq1, rq2, 
  long double rv3;   //rv0, rv1, rv2, 


  memcpy(tmp0, array.Q, input.N*sizeof(fftwl_complex));
  memcpy(tmp1, array.V, input.N*sizeof(fftwl_complex));
  fftwl_execute(p0b);
  fftwl_execute(p1b);
  //spec_output_nofft("amispectrum.txt");
  for (long int j = input.N/2; j > -1; j--) {
    EQ0 += creall(tmp0[j]*conjl(tmp0[j]));
    EV0 += creall(tmp1[j]*conjl(tmp1[j]));
    if (j >= input.N/8) {
       EQ1 += creall(tmp0[j]*conjl(tmp0[j]));
       EV1 += creall(tmp1[j]*conjl(tmp1[j]));
       if (j >= 3*input.N/16) {
          EQ2 += creall(tmp0[j]*conjl(tmp0[j]));
          EV2 += creall(tmp1[j]*conjl(tmp1[j]));
          if (j >= 7*input.N/32) {
   	    EQ4 += creall(tmp0[j]*conjl(tmp0[j]));
            EV4 += creall(tmp1[j]*conjl(tmp1[j]));
	    if (j >= 15*input.N/64) {
	       EQ8 += creall(tmp0[j]*conjl(tmp0[j]));
               EV8 += creall(tmp1[j]*conjl(tmp1[j]));
	    }
          }
       }
    }
  }
  //rq0 = sqrtl(EQ1/EQ0);
  //rq1 = sqrtl(EQ2/EQ0);
  //rq2 = sqrtl(EQ4/EQ0);
  rq3 = sqrtl(EQ8/EQ0);

  //rv0 = sqrtl(EV1/EV0);
  //rv1 = sqrtl(EV2/EV0);
  //rv2 = sqrtl(EV4/EV0);
  rv3 = sqrtl(EV8/EV0);
  
  return rq3+1.IL*rv3;
}
*/

/*
void backup_arrays() {
  qB = fftwl_malloc(input.N*sizeof(fftwl_complex));
  vB = fftwl_malloc(input.N*sizeof(fftwl_complex));
  // save arrays before L1 ref
  memcpy(qB, array.Q, input.N*sizeof(fftwl_complex)); 
  memcpy(vB, array.V, input.N*sizeof(fftwl_complex));
}
*/

/*
void restore_arrays() {
  memcpy(array.Q, qB, input.N*sizeof(fftwl_complex)); 
  memcpy(array.V, vB, input.N*sizeof(fftwl_complex));
}
*/

/*
void call_L2_refine() {
	char fname[80];
	move_mesh(track_singularity(), input.L, 2*input.N);
    	//fftwl_complex mon = check_resolution();
	     
	sprintf(fname, "spec_l2r%ld.txt", input.refN+1);	
	get_spectrum(fname);
	sprintf(fname, "L2: Using Gravity CFL: dt = %.12LE\n", motion.time);
	debug_msg(fname, EXIT_FALSE);
#if RK4
    	modify_rk4skip(sqrtl(1.5L));
#elif DIRK4
    	modify_dirk4skip(sqrtl(1.5L));
#endif
}
*/


/*
void resolution_monitor(long int *iter) {  
  // this actually works quite well
  char fname[512];
  fftwl_complex mon = check_resolution();
  
  if (  (creall(mon)>1e-12*sqrt(input.N))||(cimagl(mon)>1e-12*sqrt(input.N)) ) {
     sprintf(fname, "spec_br%ld.txt", input.refN+1);
     //spec_output(fname, tmp0, motion.time);
     get_spectrum(fname);
     get_integrals();
     if (fabs(motion.ml) > ml_tolerance) printf("Mean Level too inaccurate.\n");

     printf("Refinement: L0");

     sprintf(fname, "Before Ref: T=%.19LE\tML = %.19LE\ty0 = %.19LE\n", motion.time, motion.ml, motion.y0);
     debug_msg(fname, EXIT_FALSE);

     move_mesh(track_singularity(), input.L, input.N);
     mon = check_resolution();
     *iter = 0;
     if (( creall(mon)>2e-13*sqrt(input.N))||(cimagl(mon)>2e-13*sqrt(input.N))) {
	     sprintf(fname, "L1, ref #%ld\n", ref_nog);     
  	     debug_msg(fname, EXIT_FALSE);
	     printf(" L1:");

	     backup_arrays();
	     long double posS = track_singularity();
	     if (input.g) {
	     	move_mesh(posS, input.L*sqrtl(2.0L), input.N);  // Lin = L/2 in Matlab
     	    mon = check_resolution();
	       	printf(" coarsen");
	       	sprintf(fname, "L1 with g: revert. ref #%ld\n", ref_nog);
	       	debug_msg(fname, EXIT_FALSE);
	       	if (( creall(mon)>2e-13*sqrt(input.N))||(cimagl(mon)>2e-13*sqrt(input.N))) {
 	        	printf(" failed");
	         	sprintf(fname,"L1 with g: revert failed. ref #%ld\n", ref_nog);
		 		debug_msg(fname, EXIT_FALSE);
      	        move_mesh(posS, input.L/sqrtl(2.0L), input.N);     
		 		restore_arrays();	 	 
   	         	move_mesh(posS, input.L/sqrtl(2.0L), input.N);  // Lin = L/2 in Matlab
     	        mon = check_resolution();
#if RK4
		 		modify_rk4skip(sqrtl(1.5L));
#elif DIRK4
		 		modify_dirk4skip(sqrtl(1.5L));
#endif
 	       } else {
#if RK4
		 		modify_rk4skip(1.L/sqrtl(1.5L));
#elif DIRK4
		 		modify_dirk4skip(1.L/sqrtl(1.5L));
#endif
		 		if (fabs(motion.ml) > ml_tolerance) printf("Refine NOW!\n");
 	         	printf(" success\n.");
		 		sprintf(fname, "L1 with g: revert success!. ref #%ld\n", ref_nog);
  	         	debug_msg(fname, EXIT_FALSE);
	       }
		} else {
 	       printf(" refine");
	       sprintf(fname, "L1 with g: forward ref #%ld\n", ref_nog);
	       debug_msg(fname, EXIT_FALSE);
	       move_mesh(posS, input.L/sqrtl(2.0L), input.N);  // Lin = L/2 in Matlab
     	   mon = check_resolution();
#if RK4
	       modify_rk4skip(sqrtl(1.5L));
#elif DIRK4
	       modify_dirk4skip(sqrtl(1.5L));
#endif
 	   	}
	   	sprintf(fname, "spec_l1r%ld.txt", input.refN+1);
	   	get_spectrum(fname);
        get_integrals();
	   	ref_nog++;	
	   	sprintf(fname, "L1: Using Gravity CFL: dt = %.12LE\n", motion.time);
       	debug_msg(fname, EXIT_FALSE);
		if ((( creall(mon)>2e-13*sqrt(input.N))||(cimagl(mon)>2e-13*sqrt(input.N)))||(fabs(motion.ml) > ml_tolerance)   ) {
	    	sprintf(fname, "L1, ref #%ld failed.\n", ref_nog-1);
         	debug_msg(fname, EXIT_FALSE);
 	     	printf(" failed. L2: success.\n");
   	     	move_mesh(posS, input.L*sqrtl(2.0L), input.N);
   	     	restore_arrays();	
   	     	call_L2_refine();
	     	///*
	     	//move_mesh(track_singularity(), input.L, 2*input.N);
     	        //mon = check_resolution();
	     	//sprintf(fname, "spec_l2r%ld.txt", input.refN+1);
  	     	//get_spectrum(fname);
                //sprintf(fname, "L2: Using Gravity CFL: dt = %.12LE\n", motion.time);
                //debug_msg(fname, EXIT_FALSE);
#if RK4
  	     	//modify_rk4skip(sqrtl(1.5L));
#elif DIRK4
	     	//modify_dirk4skip(sqrtl(1.5L));
#endif
		
	     	ref_nog = 0;
       	} else {
 	     	printf(" success.\n");
	     	sprintf(fname, "L1 with g: success!. ref #%ld\n", ref_nog-1);
            debug_msg(fname, EXIT_FALSE);
		}
        if (( creall(mon)>1e-12*sqrt(input.N))||(cimagl(mon)>1e-12*sqrt(input.N))) {
 	   		printf("\nSpectrum too wide after refinement! Something went wrong\n");
	   		sprintf(fname, "Spectrum too wide. CFL?\n"); 
	   		debug_msg(fname, EXIT_FALSE);
           	sprintf(fname, "spec_last.txt");
	   		exit(1);
        }
        sprintf(fname, "After Ref: ML = %.12LE\ty0 = %.12LE\n", motion.ml, motion.y0);
        debug_msg(fname, EXIT_FALSE);
		input.refN++;
     } else printf(" success\n");
  }
}
*/

/*
void set_aux() {
  for (int j = 0; j < input.N; j++) {
    k[j] = 1.L*j;
    if (j > (input.N)/2-1) {
      k[j] = -1.L*(input.N - j);
    } 
  }
  extra.q = 2.0L*atan2l(input.L*sinl(0.5L*input.u),cosl(0.5L*input.u));
  for (int j = 0; j < (input.N)/2-1; j++) extra.w[j] = cexpl(-1.0IL*k[j+1]*(extra.q - 2.0IL*atanhl(input.L) ) )/input.N;	   //good
  for (int j = 0; j < input.N; j++) {
    q = 2.L*PI*(1.0L*j/input.N - 0.5L);
    //extra.dq[j] = (2.0L*input.L)/(1.0L + input.L*input.L + (1.0L - input.L*input.L)*cosl(q-extra.q));	// this is actually dU/dQ
    extra.newdU[j] = (2.0L*input.L)/(1.0L + input.L*input.L + (1.0L - input.L*input.L)*cosl(q-extra.q)); // now this is dU/dQ
    extra.newdQ[j] = (1.0L + input.L*input.L + (1.0L - input.L*input.L)*cosl(q-extra.q))/(2.0L*input.L); // and this is dQ/dU
  }
}
*/

/*
void clean_up_mesh() {
   fftwl_free(array.Q);
   fftwl_free(array.V);
   fftwl_free(z); 
   fftwl_free(Rz); 
   fftwl_free(tmp2); 
   fftwl_free(tmp3); 
   fftwl_free(k);
   fftwl_free(extra.U);
   fftwl_free(extra.dU);
   fftwl_free(extra.dQ);
   fftwl_free(extra.dV);
   fftwl_free(extra.B);
   fftwl_free(extra.w);
   
   //fftwl_free(extra.dq);
   fftwl_free(extra.newdQ);  
   fftwl_free(extra.newdU); 
#if RK4
   free_rk4();//#####################################################################
#elif DIRK4
   free_dirk4();//#####################################################################
#endif
}
*/
/*
void reinitialize_mesh() {
   if (!(z = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(Rz = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(array.Q = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(array.V = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);  
   if (!(tmp2 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp3 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(k = fftwl_malloc(input.N*sizeof(long double)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.U = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dU = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dQ = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dV = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.B = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.w = fftwl_malloc((input.N/2-1)*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
//   if (!(extra.dq = fftwl_malloc((input.N)*sizeof(long double)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.newdQ = fftwl_malloc((input.N)*sizeof(long double)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.newdU = fftwl_malloc((input.N)*sizeof(long double)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   //------------ moved init_rk4 from here
   FILE *fh = fopen("output.log","a");
   fprintf(fh, "Reinitialize Data: initialized FFTW plans with mode %d\n", FMODE);
   fclose(fh);
}
*/

/*
void reinitialize_mesh2() {
   free(tmp0); 
   free(tmp1);
   fftwl_destroy_plan(p0f);
   fftwl_destroy_plan(p1f);
   fftwl_destroy_plan(p2f);
   fftwl_destroy_plan(p0b);
   fftwl_destroy_plan(p1b);
   if (!(tmp0 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp1 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   p0f = fftwl_plan_dft_1d(input.N, tmp0, tmp0, FFTW_FORWARD, FMODE);
   p1f = fftwl_plan_dft_1d(input.N, tmp1, tmp1, FFTW_FORWARD, FMODE);
   p2f = fftwl_plan_dft_1d(input.N, tmp2, tmp2, FFTW_FORWARD, FMODE);
   p0b = fftwl_plan_dft_1d(input.N, tmp0, tmp0, FFTW_BACKWARD, FMODE);
   p1b = fftwl_plan_dft_1d(input.N, tmp1, tmp1, FFTW_BACKWARD, FMODE);
#if RK4
   diss_fltr = init_rk4();  // moved init_rk4 here //#####################################################################
#elif DIRK4
   diss_fltr = init_dirk4();  // moved init_rk4 here //#####################################################################
#endif
}
*/


/*
void move_mesh(long double uin, long double Lin, long int Nin) { // Lin = L/2 in Matlab
  long double qin = 2.0L*atan2l(Lin*sinl(0.5L*uin), cosl(0.5L*uin));
  long double beta = tanl(0.5L*(uin-input.u));
  long double q, qr;
  //long double u0, u1, q0 = 2.0L*atan2l(input.L*sinl(0.5L*input.u), cosl(0.5L*input.u));
  long int N0 = input.N;

  memcpy(tmp0, array.Q, N0*sizeof(fftwl_complex));
  memcpy(tmp1, array.V, N0*sizeof(fftwl_complex));

  fftwl_execute(p0b);
  fftwl_execute(p1b);
  clean_up_mesh();

  input.N = Nin;

  reinitialize_mesh();
  memset(array.Q,0,sizeof(fftwl_complex)*input.N);
  memset(array.V,0,sizeof(fftwl_complex)*input.N);
  //FILE *fh = fopen("move_mesh.txt","w");
  //fprintf(fh, "# 0. q 1. u1 2. u2 3.-4. V0 5.-6. V1\n\n");
  for (long int j = 0; j < input.N; j++) {
    q  = PI*(2.0L*j/input.N - 1.0L); 
    //u0 = input.u + 2.0L*atan2l(input.L*sinl(0.5L*(q-q0)),cosl(0.5L*(q-q0)));
    //u1 = uin + 2.0L*atan2l(Lin*sinl(0.5L*(q-qin)),cosl(0.5L*(q-qin)));
    qr = PI+extra.q + 2.0L*atan2l(beta + Lin*tanl(0.5L*(q-qin)), input.L*(1.0L-Lin*beta*tanl(0.5L*(q-qin)) )   );
    for (long int l = 0; l < N0/2; l++) {
      array.Q[j] += tmp0[l]*cexpl(-1.IL*l*qr)/N0;
      array.V[j] += tmp1[l]*cexpl(-1.IL*l*qr)/N0;
    }
    //fprintf(fh, "%.19LE\t%.19LE\t%.19LE\t%.19LE\t%.19LE\n", q, u0, u1, creall(array.V[j]), cimagl(array.V[j]));
  }
  //fclose(fh);  
  reinitialize_mesh2();
  input.L = Lin; 
  input.u = uin;
  set_aux();
}
*/

/*
void hfilter() {
  for (long int j = 0; j < input.N; j++) {
    tmp0[j] = array.Q[j]/input.N;
    tmp1[j] = array.V[j]/input.N;
  }
  fftwl_execute(p0b);
  fftwl_execute(p1b);
  memset(tmp0+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  memset(tmp1+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  for (long int j = 0; j < input.N/2; j++) {
    tmp0[j] = tmp0[j]*expl(-diss_fltr*powl(j,12));
    tmp1[j] = tmp1[j]*expl(-diss_fltr*powl(j,12)); 
  }
  fftwl_execute(p0f);
  fftwl_execute(p1f);
  memcpy(array.Q, tmp0, input.N*sizeof(fftwl_complex));
  memcpy(array.V, tmp1, input.N*sizeof(fftwl_complex));
}
*/

/*
void initialize_auxiliary_arrays() {
   if (!(z = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(Rz = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp0 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp1 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp2 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp3 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(k = fftwl_malloc(input.N*sizeof(long double)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   debug_msg("Initialize Data: memory allocation of auxiliary arrays successful\n", EXIT_FALSE);
   if (!(extra.U = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dU = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dQ = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dV = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.B = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.w = fftwl_malloc((input.N/2-1)*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   //if (!(extra.dq = fftwl_malloc((input.N)*sizeof(long double)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.newdQ = fftwl_malloc((input.N)*sizeof(long double)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.newdU = fftwl_malloc((input.N)*sizeof(long double)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   //------------ 
   p0f = fftwl_plan_dft_1d(input.N, tmp0, tmp0, FFTW_FORWARD, FMODE);
   p1f = fftwl_plan_dft_1d(input.N, tmp1, tmp1, FFTW_FORWARD, FMODE);
   p2f = fftwl_plan_dft_1d(input.N, tmp2, tmp2, FFTW_FORWARD, FMODE);
   p0b = fftwl_plan_dft_1d(input.N, tmp0, tmp0, FFTW_BACKWARD, FMODE);
   p1b = fftwl_plan_dft_1d(input.N, tmp1, tmp1, FFTW_BACKWARD, FMODE);
   //------------

   FILE *fh = fopen("output.log","a");
   fprintf(fh, "Initialize Data: initialized FFTW plans with mode %d\n", FMODE);
   fclose(fh);
}
*/

/*
void free_arrays() {
   debug_msg("Free Data: initiated\n",EXIT_FALSE);
   free(array.Q);
   free(array.V);
   free(z); 	  free(Rz);
   free(tmp0);    fftwl_destroy_plan(p0f);
   free(tmp1);    fftwl_destroy_plan(p1f);
   free(tmp2);    fftwl_destroy_plan(p0b);
   free(tmp3);    fftwl_destroy_plan(p1b);
   debug_msg("Free Data: complete\n",EXIT_FALSE);
}
*/

/*
void initialize_data() {
   debug_msg("Initialize Data: initiated\n", EXIT_FALSE);
   printf("Requested %ld fftwl complex type\n", input.N); 
   if (!(array.Q = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed.\nInitialize Data: complete\n", 1);
   if (!(array.V = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", 1);  
   debug_msg("Initialize Data: memory allocation successful\n", EXIT_FALSE);
}
*/

/*
void read_pade() {
  FILE *fh_pade = fopen(input.rname,"r");
  char line[256], value1[256], value2[256], value3[256], value4[256];
  long double u, y0, C_vel, HLa;
  long double num, den;
  //long int n2;
  int nreads;

  //num = fftwl_malloc(input.d*sizeof(fftwl_complex));
  //den = fftwl_malloc(input.d*sizeof(fftwl_complex));
  if (fh_pade) {
    printf("Starting from Pade approximation with %ld poles\n", input.d);
    if(!fgets(line, 256, fh_pade)) printf("read_pade: unexpected end of file\n");
    if(!fgets(line, 256, fh_pade)) printf("read_pade: unexpected end of file\n");
    nreads = sscanf(line, "# M = %s\tResidual = %s\ty0 = %s\n", value1, value2, value3);
    //n2 = atol(value1);
    y0 = strtold(value3, NULL);
    printf("Successfully read %d numbers.\n", nreads);    
    printf("# Read y0 = %.19LE\n",  y0);
    if(!fgets(line, 256, fh_pade)) printf("read_pade: unexpected end of file\n");
    nreads = sscanf(line, "# Amplitude at highest Fourier mode = %s\tPade error = %s\tH/lambda = %s\tc = %s", value1, value2, value3, value4);
    HLa = strtold(value3, NULL);
    C_vel = strtold(value4, NULL);
    printf("# Read H/L = %.19LE\n# Read   c = %.19LE\n",  HLa, C_vel);
    if(!fgets(line, 256, fh_pade)) printf("read_pade: unexpected end of file\n");

    int m = 0;
    memset(array.Q, 0, (input.N)*sizeof(fftwl_complex)); 
    while (fgets(line, 256, fh_pade)) {    
      nreads = sscanf(line, "%s\t%s\t%s\n", value1, value2, value3);
      if (nreads != 3) {
	 printf("Broken Pade File\n");
	 exit(1);
      } else {
         num = strtold(value2, NULL);
         den = strtold(value1, NULL);
         printf("%3d: %.16Le\t%.16Le\n", m, num, den);     
         for (long int j = 0; j < input.N; j ++) {
           q = PI*(2.0L*j/input.N - 1.0L);
           u = input.u + 2.0L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
           //array.Q[j] += num/((input.L)*tanl(0.5L*q) - 1.IL*den);   			//-- conformal map Z
           array.Q[j] += -0.5L*num/cpowl(sinl(0.5L*u) - 1.IL*den*cosl(0.5L*u),2);		//-- Z_u
         }
         m++;
      }
    }
    fclose(fh_pade);
    //for (long int j = 0; j < input.N; j ++) array.Q[j] += 1.0IL*y0;
    for (long int j = 0; j < input.N; j ++) {
	array.Q[j] += 1.0L;
	array.V[j] = 1.IL*C_vel*(1.L - 1.L/array.Q[j]);
    }
    for (long int j = 0; j < input.N; j ++) {
        q = PI*(2.0L*j/input.N - 1.0L);
        u = input.u + 2.0L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
	array.Q[j] += -5.0E-06IL/cpowl(sinl(0.5L*u) - 4.0E-3IL*cosl(0.5L*u),2);
	array.Q[j] = csqrtl(1.L/array.Q[j]);
    }
    //basic_output("pade.txt", array.Q, array.V, 0.L);
    fftwl_complex isResolved = check_resolution();

    printf("Verifing that Fourier spectrum is resolved\n");
    printf("Energy in high wavenumbers is %.16LE\n", cabsl(isResolved));
    if (cabsl(isResolved) > 1e-12) {
       printf("Spectrum not resolved\n");
       printf("Try with more Fourier modes\n");
       exit(1);
    }
    set_aux();
    get_integrals();
    //get_surface("starting.txt");
    
    printf("Final Set of Parameters:\nL = %.18Le\nN = %ld\n", input.L, input.N);
    //exit(1);
  } else {
    printf("Cannot open %s\n", input.rname);
    exit(1);
  }
}
*/
/*
void simulate() {
   //long double u;
   
   ref_nog = 0;
   //for (int j = 0; j < input.N; j++) {
   //  q = PI*(2.L*j/input.N - 1.0L);
   //  u = input.u + 2.L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
   //  array.Q[j] = 1.L; //1E-13*cexpl(-1.IL*u);
   //  array.V[j] = -0.05IL*(1.L/ctanl(0.5L*(u-0.12IL)) - 1.IL) + 0.00IL*(1.L/ctanl(0.5L*(u-0.24IL)) - 1.IL); // run 8&9

   //  array.Q[j] = 1.L + 0.5L*cexpl(-1.IL*u); 	//+0.4 (run_4) // +0.3 (run_6)
   //  array.V[j] = -0.75IL*exp(-1.IL*u);		//-0.85I (run_5) // -0.6 (run 6)
   //}



#if RK4
   evolve_rk4();
   free_arrays();
#elif DIRK4
   evolve_dirk4();
   //dirk4_step();
#endif
}
*/
