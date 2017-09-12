#include "header.h"

static mpfc_t temp, temp2, buf, buf2;

void init_mpfc_arithmetic() {
  mpfr_inits2(precision, temp.re,   temp.im, (mpfr_ptr) NULL);
  mpfr_inits2(precision, temp2.re, temp2.im, (mpfr_ptr) NULL);
  mpfr_inits2(precision, buf.re,     buf.im, (mpfr_ptr) NULL);
  mpfr_inits2(precision, buf2.re,   buf2.im, (mpfr_ptr) NULL);
}

void mpfc_set_d(mpfc_t *rop, double real, double imag, mpfr_rnd_t rnd) {
  mpfr_set_d(rop->re, real, rnd);
  mpfr_set_d(rop->im, imag, rnd);
}

void mpfc_set_si(mpfc_t *rop, int real, int imag, mpfr_rnd_t rnd) {
  mpfr_set_si(rop->re, real, rnd);
  mpfr_set_si(rop->im, imag, rnd);
}

void mpfc_set(mpfc_t *rop, mpfc_t *op, mpfr_rnd_t rnd) {
  mpfr_set(rop->re, op->re, rnd);
  mpfr_set(rop->im, op->im, rnd);
}

void mpfc_neg(mpfc_t *rop, mpfc_t *op, mpfr_rnd_t rnd) {
  mpfr_neg(rop->re, op->re, rnd);
  mpfr_neg(rop->im, op->im, rnd);
}

void mpfc_div(mpfc_t *rop, mpfc_t *op1, mpfc_t *op2, mpfr_rnd_t rnd) {
  mpfc_conj(&temp2, op2, rnd);
  mpfc_mul(&buf2, op1, &temp2, rnd);
  mpfc_mul(&temp2, op2, &temp2, rnd);
  mpfr_div(rop->re, buf2.re, temp2.re, rnd);
  mpfr_div(rop->im, buf2.im, temp2.re, rnd);
}

void mpfc_si_div(mpfc_t *rop, int op1, mpfc_t *op2, mpfr_rnd_t rnd) {
  mpfc_conj(&temp2, op2, rnd);
  mpfc_mul_si(&buf2, &temp2, op1, rnd);
  mpfc_mul(&temp2, op2, &temp2, rnd);
  mpfr_div(rop->re, buf2.re, temp2.re, rnd);
  mpfr_div(rop->im, buf2.im, temp2.re, rnd);
}

void mpfc_mul(mpfc_t *rop, mpfc_t *op1, mpfc_t *op2, mpfr_rnd_t rnd) {
  mpfr_mul(buf.re, op1->im, op2->im, rnd);
  mpfr_fms(rop->re, op1->re, op2->re, buf.re, rnd);
  mpfr_mul(buf.im, op1->re, op2->im, rnd);
  mpfr_fma(rop->im, op1->im, op2->re, buf.im, rnd);
}

void mpfc_mul_si(mpfc_t *rop, mpfc_t *op1, int op2, mpfr_rnd_t rnd) {
  mpfr_mul_si(rop->re, op1->re, op2, rnd);
  mpfr_mul_si(rop->im, op1->im, op2, rnd);
}
void mpfc_add(mpfc_t *rop, mpfc_t *op1, mpfc_t *op2, mpfr_rnd_t rnd) {
  mpfr_add(rop->re, op1->re, op2->re, rnd);
  mpfr_add(rop->im, op1->im, op2->im, rnd);
}

void mpfc_sub(mpfc_t *rop, mpfc_t *op1, mpfc_t *op2, mpfr_rnd_t rnd) {
  mpfr_sub(rop->re, op1->re, op2->re, rnd);
  mpfr_sub(rop->im, op1->im, op2->im, rnd);
}

void mpfc_fms(mpfc_t *rop, mpfc_t *op1, mpfc_t *op2, mpfc_t *op3, mpfr_rnd_t rnd){
  mpfc_mul(&temp2, op1, op2, rnd);
  mpfr_sub(rop->re, temp2.re, op3->re, rnd);
  mpfr_sub(rop->im, temp2.im, op3->im, rnd);
}

void mpfc_fma(mpfc_t *rop, mpfc_t *op1, mpfc_t *op2, mpfc_t *op3, mpfr_rnd_t rnd){
  mpfc_mul(&temp2, op1, op2, rnd);
  mpfr_add(rop->re, temp2.re, op3->re, rnd);
  mpfr_add(rop->im, temp2.im, op3->im, rnd);
}

void mpfc_conj(mpfc_t *rop, mpfc_t *op, mpfr_rnd_t rnd) {
  mpfr_set(rop->re, op->re, rnd);
  mpfr_neg(rop->im, op->im, rnd); 
}
