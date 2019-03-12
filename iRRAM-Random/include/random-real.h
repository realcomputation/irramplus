/*
iRRAM-Random is a library providing randomly generated
continuous objects in iRRAM.
The library is MIT license protected.

This header file is for providing random real number generators
*/
#ifndef RANDOMREAL_H
#define RANDOMREAL_H

#include <random-core.h>

namespace iRRAM::random {

REAL uniform_real();
REAL uniform_real(REAL, REAL);

REAL gaussian_real();
REAL gaussian_real(REAL, REAL);

REAL linear_real();

}

#endif
