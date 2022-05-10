#ifndef _FFT_H_
#define _FFT_H_

#include <complex.h>
#include <stddef.h>
#include "mpi_data.h"

// procedura licząca tranformatę Fouriera.
void fft_forward(double complex *data, size_t N, MPI_Data mpi_data);

// procedura licząca tranformatę odwrtoną Fouriera.
void fft_backward(double complex *data, size_t N, MPI_Data mpi_data);

void fft_filter(double threshold, double complex *data, size_t N, MPI_Data mpi_data);

#endif
