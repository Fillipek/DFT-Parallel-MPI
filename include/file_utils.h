#ifndef _FILE_UTILS_H_
#define _FILE_UTILS_H_

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

/* Reads number of lines in a file. */
void count_lines_in_file(const char *filename, size_t *N);

/* Reads numerical data from a file. Each line should contain one number. 
 * Buffer should be initialized with the size of at least N elements.
 */
void read_data_lines(const char *filename, size_t N, double complex *data_buff);

/* Svaes complex number array into a file. Read and imaginary parts are
 * separated with a comma (,).
 */
void save_to_file(const char* filename, size_t N, double complex *data);

#endif