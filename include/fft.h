#ifndef _FFT_H_
#define _FFT_H_

#include "complex_own.h"
#include <stddef.h>

// procedura licząca tranformatę Fouriera.
void dft_forward(complex_double *data, size_t N, int, int);

// procedura licząca tranformatę odwrtoną Fouriera.
void dft_backward(complex_double *data, size_t N, int, int);

#endif
