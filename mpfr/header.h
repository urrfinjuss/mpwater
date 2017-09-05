#include <math.h>
//#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
//#include <fftw3.h>
#include <mpfr.h>
#include <mpfft_serial.h>

#define MIN(a,b) ((a) < (b) ? a : b)

#define FFT_FORWARD	-1
#define FFT_BACKWARD	1

#define MOVE_MESH 0		// 0 if singularity tracking is off , 1 otherwise
#define PADE_TEST 0

// --------  Structures
/*
typedef struct multiprecision_complex {
  mpfr_t re, im;
} mpfc_t;
*/
typedef struct input {
  char restart_name[128];		// restart filename 
  char txt_format[64];			// format of text file (pade, ascii, none)
  char time_march[64];			// time-marching (rk4)
  mpfr_prec_t precision;
  mpfr_t gravity;			// free fall acceleration	*
  mpfr_t surface_tension;		// surface tension		*
  mpfr_t tolerance;			// tolerance for refinement	*
  mpfr_t mean_level;			// mean level fluid		x
  mpfr_t kineticE;			// kinetic energy		x
  mpfr_t potentialE;			// potential energy		x
  mpfr_t final_time;			// simulation time		*
  mpfr_t time;				// stores current time		x
  mpfc_t momentum;			// momentum P = px + i*py	x
  unsigned long int refinement_counter;	// refinement counter
  unsigned long int number_poles;	// number of poles
  unsigned long int nbits;		// number of grid points 2^nbits
  unsigned long int kD;			// hyperV scale, set in evolve.c
} params, *params_ptr;

typedef struct conformal_mapping {
  mpfc_t 	*w;			// Fourier multipliers
  mpfr_t 	*dq;			// dq/du
  mpfr_t 	scaling;		// conformal map scaling factor
  mpfr_t 	image_offset;		// accumulation center in u-plane
  mpfr_t 	origin_offset;		// accumulation center in q-plane
} map, *map_ptr;

typedef struct pade_data {
  unsigned int  n_lins;			// number of iterations for Q^{k}
  unsigned int	n_poles;		// degree of denominator (number of poles)
  mpfr_t	l2_nrm;		// L2 norm of the target function
  mpfr_t   l2_abs_err;		// L2 norm of (P/Q - W) dq
  mpfr_t 	l2_rel_err;		// L2 norm of (P/Q - W) / L2 norm of W
} pade, *pade_ptr;

// -------- Global Variables
extern params 		state;
extern map 		conf, alt_map;
extern pade		pade_data;
extern mpfr_t 	**tmpr;
extern mpfc_t 	**tmpc;
extern mpfc_t	**data;
extern mpfft_plan 	ft0, ft1, ft2, ft3, ft4;
extern mpfft_plan 	ift0, ift1, ift2, ift3, ift4;
// --------  MPFR constants
extern mpfr_t   Pie;

// --------  Functions
// memory.c
extern void allocate_memory();
extern void deallocate_memory();
extern void remap(map_ptr new_map, unsigned long int N); 
extern void init_memory();
extern void fft_shift(mpfc_t *in);

// input.c
extern void load_ascii();
extern void set_initial_data();
extern void set_initial_JW();
extern void load_parameters(int argc, char *argv[]);
extern void read_input(char *fname);

// array_func.c
extern void init_arrayf();
extern void div_jacobian(mpfc_t *in, mpfc_t *out);
extern void inverse(mpfc_t *a, mpfc_t *x);
extern void linear_solve(mpfc_t *a, mpfc_t *b, mpfc_t *x);
extern void square_ft(mpfc_t *Z, mpfc_t *x);
extern void compute_zero_mode(mpfc_t *in, mpfr_t S0, mpfr_t *out);
extern void compute_zero_mode_complex(mpfc_t *in, mpfc_t S0, mpfc_t *out);

// hlevel.c
extern void project(mpfc_t *in, mpfc_t *out);
extern void convertZtoQ(mpfc_t *in, mpfc_t *out);
extern void convertQtoZ(mpfc_t *in, mpfc_t *out);
extern void restore_potential(mpfc_t *inQ, mpfc_t *inV, mpfc_t *out);

// mapping.c
extern void set_mapping();
extern void map_quality(mpfc_t *in1, mpfc_t *in2, mpfr_t tol, unsigned int *QC_pass);
extern void map_quality_fourier(mpfc_t *inQ, mpfc_t *inV, mpfr_t tol, unsigned int *QC_pass);
extern void track_singularity(mpfc_t *inQ);

// output.c
extern void complex_array_out(char *fname, mpfc_t *in);

extern void real_array_out(char* fname, mpfr_t *in);
extern void surface_out(char *fname, mpfc_t *in);
extern void spec_out(char *fname, mpfc_t *in1, mpfc_t *in2);
extern void output_data(char *fname, mpfc_t *inPhi);
extern void print_constants();

// pade.c
extern void init_pade();
extern void allocate_pade(unsigned long nD);
extern void deallocate_pade();
extern void pade_array_out(char *fname, mpfc_t *in);
extern void compute_rational(unsigned long nD, unsigned long n_max_iter, mpfc_t *in);
extern void optimal_pade(char *str, mpfc_t *in);
extern void find_l2_error(pade_ptr inp);
extern void print_pade(pade_ptr inp);
extern void newton_search(unsigned long nD);
extern void verify_pade(mpfc_t *residues, mpfc_t *roots, unsigned int nD);
extern void poly_val_array(mpfc_t *in, unsigned long nD, mpfc_t *outQ, mpfc_t *outQp, mpfc_t *outP);
extern void aberth_iter(unsigned int nD, char *str);
extern void sort_by_imag(mpfc_t *in1, mpfc_t *in2, unsigned int nD);

// evolve.c
extern void compute_rhs(mpfc_t *inQ, mpfc_t *inV, mpfc_t *outQ, mpfc_t *outV);
extern void init_timemarching();
extern void allocate_timemarching();
extern void deallocate_timemarching();
extern void rk6_step(mpfc_t *inQ, mpfc_t *inV, mpfr_t dt);
extern void evolve_rk6();

