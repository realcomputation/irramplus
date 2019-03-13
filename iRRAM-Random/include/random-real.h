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

namespace iRRAM::random {

/* seed will be used to obtain a seed for the random number engine */
/* TODO: switch to iRRAM::random_device in order for reiterations to return
 *       consistent (discrete for arguments, continuous for results) values */
REAL uniform_real(unsigned int seed = std::random_device{}());
REAL uniform_real(REAL, REAL);

REAL gaussian_real();
REAL gaussian_real(REAL, REAL);

REAL linear_real();

}

#endif
