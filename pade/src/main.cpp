extern "C" {
	// Include a C header, wrap in extern "C"
	#include "header.h"
}
#include "header.hpp"

using namespace std;


int main() {
	foo();
	mp_complex w1(1), w2(1,1), z;	
	z = w1+w2;
	z.print_number();
	z = w1/w2;
	z.print_number();

}
