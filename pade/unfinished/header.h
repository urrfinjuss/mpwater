#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpfr.h>
#include <mpfft_header.h>

#define MOVE_MESH 0		// 0 if singularity tracking is off , 1 otherwise
#define PADE_TEST 0

#define MODE MPFR_RNDN
// --------  Structures
/*
typedef struct multiprecision_complex {
  mpfr_t re, im;
} mpfc_t, *mpfc_ptr;
*/

typedef struct pade_data {
  unsigned int  n_lins;			// number of iterations for Q^{k}
  unsigned int	n_poles;		// degree of denominator (number of poles)
  unsigned int  Npt;			// input data array length
  mpfr_t	*Loc, *Res;		// locatioan and residues at poles
  mpfr_t	l2_nrm;		// L2 norm of the target function
  mpfr_t        l2_abs_err;		// L2 norm of (P/Q - W) dq
  mpfr_t 	l2_rel_err;		// L2 norm of (P/Q - W) / L2 norm of W
} pade, *pade_ptr;

// -------- Global Variables
extern pade		pade_data, best_pade;
extern mpfr_prec_t	precision;
extern mpfc_t		**Gees, **Cees;
extern mpfc_t		**A, **B;
extern mpfc_t		*arrQ, *arrP, *Tmp, *W;
extern mpfr_t		*Ens, *M;
extern mpfr_t		Pie, Ovn, Que;

// --------  Functions
// aberth.c
extern void init_aberth(unsigned int nD);
extern void free_aberth(unsigned int nD);
extern void aberth_iter(unsigned int nD, char *str);

// output.c
extern void init_output();
extern void pade_real_out(char *fname, mpfr_t *in);
extern void pade_complex_out(char *fname, mpfc_t *in);
extern void print_pade();

// memory.c
extern void init_memory(unsigned int nD);
extern void allocate_pade(unsigned long nD);
extern void deallocate_pade(unsigned long nD);

// newtsearch.c
extern void newton_search(unsigned long nD);

// set.c
extern void mpfc_dotpr(mpfc_t *out, mpfc_t *op1, mpfc_t *op2);
extern void set_weight();
extern void set_initial(unsigned int nD);

// poly.c
extern void init_poly(unsigned int nD);
extern void evaluate_poly_array(unsigned long nD);
extern void evaluate_poly(mpfc_t *in, unsigned long nD, mpfc_t *outQ, mpfc_t *outQp, mpfc_t *outP);
extern void poly_val_array(mpfc_t *in, unsigned long nD, mpfc_t *outQ, mpfc_t *outQp, mpfc_t *outP);

// gramschmidt.c
extern void init_grams();
extern void free_grams();
extern void gram_schmidt(unsigned long nD);

// rational.c
extern void compute_rational(unsigned long nD, unsigned long n_max_iter, mpfc_t *in);
extern void errorl2();

//postproc.c
extern void sort_imag(mpfc_t *in1, mpfc_t *in2, unsigned int nD);
extern void verify_pade(mpfc_t *residues, mpfc_t *roots, unsigned int nD);

// optimal.c
extern void optimal_pade(char *str, mpfc_t *in);

// pade.c


