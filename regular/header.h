#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <fftw3.h>

#define MIN(a,b) ((a) < (b) ? a : b)

#define FMODE FFTW_MEASURE	// changed from measure
#define MOVE_MESH	1	// 0 if singularity tracking is off , 1 otherwise
#define PADE_TEST	0
#define PADE_FLAG	0
#define PI acosl(-1.0L)
// --------  Structures
typedef struct input {
  char restart_name[128];		// restart filename 
  char txt_format[64];			// format of text file (pade, ascii, none)
  char time_march[64];			// time-marching (rk4)
  long double gravity;			// free fall acceleration
  long double surface_tension;		// surface tension
  long double tolerance;		// tolerance for refinement
  long double mean_level;		// mean level fluid
  long double kineticE;			// kinetic energy
  long double potentialE;		// potential energy
  long double final_time;		// simulation time
  long double time;			// stores current time
  fftwl_complex momentum;		// momentum P = px + i*py
  unsigned long int refinement_counter;	// refinement counter
  unsigned long int number_poles;	// number of poles
  unsigned long int number_modes;	// number of grid points
  unsigned long int kD;			// hyperV scale, set in evolve.c
} params, *params_ptr;

typedef struct conformal_mapping {
  fftwl_complex *w;			// Fourier multipliers
  long double 	*dq;			// dq/du
  long double 	scaling;		// conformal map scaling factor
  long double 	image_offset;		// accumulation center in u-plane
  long double 	origin_offset;		// accumulation center in q-plane
} map, *map_ptr;

typedef struct pade_data {
  unsigned int  n_lins;			// number of iterations for Q^{k}
  unsigned int	n_poles;		// degree of denominator (number of poles)
  long double	l2_nrm;		// L2 norm of the target function
  long double   l2_abs_err;		// L2 norm of (P/Q - W) dq
  long double 	l2_rel_err;		// L2 norm of (P/Q - W) / L2 norm of W
} pade, *pade_ptr;

// -------- Global Variables
extern params 		state;
extern map 		conf, alt_map;
extern pade		pade_data;
extern long double 	**tmpr;
extern fftwl_complex 	**tmpc;
extern fftwl_complex	**data;
extern fftwl_plan 	ft0, ft1, ft2, ft3, ft4;
extern fftwl_plan 	ift0, ift1, ift2, ift3, ift4;

// --------  Functions
// memory.c
extern void allocate_memory();
extern void deallocate_memory();
extern void remap(map_ptr new_map, unsigned long int N); 
extern void init_memory();
extern void fft_shift(fftwl_complex *in);

// input.c
extern void load_ascii();
extern void load_pade();
extern void set_initial_data();
extern void set_initial_JW();
extern void load_parameters(int argc, char *argv[]);
extern void read_input(char *fname);

// array_func.c
extern void init_arrayf();
extern void div_jacobian(fftwl_complex *in, fftwl_complex *out);
extern void inverse(fftwl_complex *a, fftwl_complex *x);
extern void linear_solve(fftwl_complex *a, fftwl_complex *b, fftwl_complex *x);
extern void square_ft(fftwl_complex *Z, fftwl_complex *x);
extern void compute_zero_mode(fftwl_complex *in, long double S0, long double *out);
extern void compute_zero_mode_complex(fftwl_complex *in, fftwl_complex S0, fftwl_complex *out);

// hlevel.c
extern void project(fftwl_complex *in, fftwl_complex *out);
extern void convertZtoQ(fftwl_complex *in, fftwl_complex *out);
extern void convertQtoZ(fftwl_complex *in, fftwl_complex *out);
extern void restore_potential(fftwl_complex *inQ, fftwl_complex *inV, fftwl_complex *out);

// mapping.c
extern void set_mapping();
extern void map_quality(fftwl_complex *in1, fftwl_complex *in2, long double tol, unsigned int *QC_pass);
extern void map_quality_fourier(fftwl_complex *inQ, fftwl_complex *inV, long double tol, unsigned int *QC_pass);
extern void track_singularity(fftwl_complex *inQ);

// output.c
extern void real_array_out(char* fname, long double *in);
extern void complex_array_out(char *fname, fftwl_complex *in);
extern void surface_out(char *fname, fftwl_complex *in);
extern void spec_out(char *fname, fftwl_complex *in1, fftwl_complex *in2);
extern void output_data(char *fname, fftwl_complex *inPhi);
extern void print_constants();

// pade.c
extern void init_pade();
extern void allocate_pade(unsigned long nD);
extern void deallocate_pade();
extern void pade_array_out(char *fname, fftwl_complex *in);
extern void compute_rational(unsigned long nD, unsigned long n_max_iter, fftwl_complex *in);
extern void optimal_pade(char *str, fftwl_complex *in);
extern void find_l2_error(pade_ptr inp);
extern void print_pade(pade_ptr inp);
extern void newton_search(unsigned long nD);
extern void verify_pade(fftwl_complex *residues, fftwl_complex *roots, unsigned int nD);
extern void poly_val_array(fftwl_complex *in, unsigned long nD, fftwl_complex *outQ, fftwl_complex *outQp, fftwl_complex *outP);
extern void aberth_iter(unsigned int nD, char *str);
extern void sort_by_imag(fftwl_complex *in1, fftwl_complex *in2, unsigned int nD);

// evolve.c
extern void compute_rhs(fftwl_complex *inQ, fftwl_complex *inV, fftwl_complex *outQ, fftwl_complex *outV);
extern void init_timemarching();
extern void allocate_timemarching();
extern void deallocate_timemarching();
extern void rk6_step(fftwl_complex *inQ, fftwl_complex *inV, long double dt);
extern void evolve_rk6();

