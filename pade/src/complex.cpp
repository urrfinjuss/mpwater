extern "C" {
	// Include a C header, wrap in extern "C"
	#include "header.h"
}

#include "header.hpp"

using namespace std;

mpfr_prec_t	precision = 256;
#define MODE	MPFR_RNDN
// overloading arithmetic
mp_complex mp_complex::operator+(const mp_complex& b) {
	mp_complex result;
	mpfr_inits2(precision, result.re, result.im, (mpfr_ptr) NULL);	
	mpfr_add(result.re, this->re, b.re, MODE);
	mpfr_add(result.im, this->im, b.im, MODE);
	return result;
}

mp_complex mp_complex::operator-(const mp_complex& b) {
	mp_complex result;
	mpfr_inits2(precision, result.re, result.im, (mpfr_ptr) NULL);	
	mpfr_sub(result.re, this->re, b.re, MODE);
	mpfr_sub(result.im, this->im, b.im, MODE);
	return result;
}

mp_complex mp_complex::operator*(const mp_complex& b) {
	mp_complex result;
	mpfr_inits2(precision, result.re, result.im, (mpfr_ptr) NULL);	
	mpfr_mul(result.re, this->im, b.im, MODE);
	mpfr_fms(result.re, this->re, b.re, result.re, MODE);
	mpfr_mul(result.im, this->re, b.im, MODE);
	mpfr_fma(result.im, this->im, b.re, result.im, MODE);
	return result;
}
mp_complex mp_complex::operator/(const mp_complex& b) {
	mp_complex	result;
	mpfr_t		buf;
	mpfr_inits2(precision, result.re, result.im, buf, (mpfr_ptr) NULL);	
	mpfr_mul(buf, b.re, b.re, MODE);
	mpfr_fma(buf, b.im, b.im, buf, MODE);
	mpfr_mul(result.re, this->im, b.im, MODE);
	mpfr_fma(result.re, this->re, b.re, result.re, MODE); 
	mpfr_mul(result.im, this->re, b.im, MODE);
	mpfr_fms(result.im, this->im, b.re, result.im, MODE);
	mpfr_div(result.re, result.re, buf, MODE);
	mpfr_div(result.im, result.im, buf, MODE);
	mpfr_clear(buf);
	return result;
}

// initialize but not set
mp_complex::mp_complex () {
	mpfr_init2(re, precision);
	mpfr_init2(im, precision);
}

// initialize with integers
mp_complex::mp_complex (int x) {
	mpfr_init2(re, precision);
	mpfr_init2(im, precision);
	mpfr_set_si(re, x, MODE);
	mpfr_set_ui(im, 0, MODE);
}

mp_complex::mp_complex (int x, int y) {
	mpfr_init2(re, precision);
	mpfr_init2(im, precision);
	mpfr_set_si(re, x, MODE);
	mpfr_set_ui(im, y, MODE);
}

// initialize with double floats
mp_complex::mp_complex (double x) {
	mpfr_init2(re, precision);
	mpfr_init2(im, precision);
	mpfr_set_d(re, x, MODE);
	mpfr_set_ui(im, 0, MODE);
}

mp_complex::mp_complex (double x, double y) {
	mpfr_init2(re, precision);
	mpfr_init2(im, precision);
	mpfr_set_d(re, x, MODE);
	mpfr_set_d(im, y, MODE);
}

// initialize with long double floats
mp_complex::mp_complex (long double x) {
	mpfr_init2(re, precision);
	mpfr_init2(im, precision);
	mpfr_set_ld(re, x, MODE);
	mpfr_set_ui(im, 0, MODE);
}

mp_complex::mp_complex (long double x, long double y) {
	mpfr_init2(re, precision);
	mpfr_init2(im, precision);
	mpfr_set_ld(re, x, MODE);
	mpfr_set_ld(im, y, MODE);
}

// initialize with MPFR types
mp_complex::mp_complex (mpfr_t x) {
	mpfr_init2(re, precision);
	mpfr_init2(im, precision);
	mpfr_set(re, x, MODE);
	mpfr_set_ui(im, 0, MODE);
}

mp_complex::mp_complex (mpfr_t x, mpfr_t y) {
	mpfr_init2(re, precision);
	mpfr_init2(im, precision);
	mpfr_set(re, x, MODE);
	mpfr_set(im, y, MODE);
}

