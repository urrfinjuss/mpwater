#include <mpfr.h>


class mp_complex {
	public:
		mp_complex ();
		mp_complex (int);
		mp_complex (int, int);
		mp_complex (double);
		mp_complex (double, double);
		mp_complex (long double);
		mp_complex (long double, long double);
		mp_complex (mpfr_t);
		mp_complex (mpfr_t, mpfr_t);
		mp_complex operator+(const mp_complex& b);
		mp_complex operator-(const mp_complex& b);
		mp_complex operator*(const mp_complex& b);
		mp_complex operator/(const mp_complex& b);

		void print_number() {
			mpfr_printf("%12Re\t%.12Re\n", this->re, this->im);
		}
	private:
		mpfr_t re;
		mpfr_t im;
};
