/*
iRRAM-Random is a library providing randomly generated
continuous objects in iRRAM.
The library is MIT license protected.

This header file is for providing random real number generators
*/
#ifndef RANDOMREAL_H
#define RANDOMREAL_H

#include <random>	/* std::random_device */

#include <iRRAM/lib.h>

namespace iRRAM {

/* TODO: inside limits returns unrelated results wrt. reiterations */
class random_device : std::random_device {
public:
	using std::random_device::result_type;
	using std::random_device::random_device;
	using std::random_device::operator=;
	using std::random_device::entropy;
	using std::random_device::min;
	using std::random_device::max;

	result_type operator()();
};

}

namespace iRRAM::random {

/* seed will be used to obtain a seed for the random number engine */
REAL uniform_real(unsigned int seed = random_device{}());
REAL uniform_real(REAL, REAL);

REAL gaussian_real();
REAL gaussian_real(REAL, REAL);

REAL linear_real();

}

#endif
