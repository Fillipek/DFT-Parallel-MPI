#ifndef _FFT_H_
#define _FFT_H_

#include <complex.h>
#include <stddef.h>
#include "mpi_data.h"

// procedura licząca tranformatę Fouriera.
void dft_forward(double complex *data, size_t N, MPI_Data mpi_data);

// procedura licząca tranformatę odwrtoną Fouriera.
void dft_backward(double complex *data, size_t N, MPI_Data mpi_data);

#endif
