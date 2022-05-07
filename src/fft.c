#include "fft.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Discrette Fourier Transform naive implementation (from definition)
void dft_naive(complex_double *in, complex_double *out, size_t N)
{
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            complex_double temp = {
                .re = cos(-2 * M_PI / N * k * n),
                .im = sin(-2 * M_PI / N * k * n)
            };
            out[k] = add(out[k], mul(in[n], temp));
        }
    }
}

// Fast Fourier Transform using recursive Cooley-Tukey algorithm (Radix-2)
void fft_radix2(complex_double *in, complex_double *out, size_t N, size_t s)
{
    if (N == 1)
    {
        *out = *in;
        return;
    }
    fft_radix2(in, out, N/2, 2*s);
    fft_radix2(in+s, out+N/2, N/2, 2*s);
    for (int k = 0; k<N/2; k++)
    {
        complex_double p = out[k];
        complex_double t = {
            .re = cos(-2 * M_PI * k / N),
            .im = sin(-2 * M_PI * k / N)
        };
        complex_double q = mul(t, out[k + N/2]);
        out[k] = add(p, q);
        out[k+N/2] = sub(p, q);
    }
}

// Implementacja z: https://www.geeksforgeeks.org/write-an-efficient-c-program-to-reverse-bits-of-a-number/: 
static int rev(unsigned int num, size_t N)
{
    unsigned int NO_OF_BITS = log2(N);
    unsigned int reverse_num = 0;
    int i;
    for (i = 0; i < NO_OF_BITS; i++) {
        if ((num & (1 << i)))
            reverse_num |= 1 << ((NO_OF_BITS - 1) - i);
    }
    return reverse_num;
}


static void bit_reversal_copy(complex_double* result, complex_double* input, size_t N)
{
    for(unsigned int k = 0; k < N - 1; k++)
    {
        result[rev(k, N)] = input[k];
    }
}


void fft_radix2_iter(complex_double *in, complex_double *out, size_t N, size_t stride)
{
    bit_reversal_copy(out, in, N);

    for(int s = 1; s <= log2(N); s++)
    {
        int m = pow(2, s);
        complex_double omega_m = { .re = cos(2 * 3.141592 / m), .im = sin(-2 * 3.141592 / m) };

        for(int k = 0; k <= (N - 1) / m; k++)
        {
            complex_double omega = { .re = 1., .im = 0. };

            for(int j = 0; j <= m / 2 - 1; j++)
            {
                complex_double t = mul(omega, out[k + j + m / 2]);
                complex_double u = out[k + j];
                out[k + j] = add(u, t);
                out[k + j + m / 2] = sub(u, t);
                omega = mul(omega, omega_m);
            }
        }
    }
}

void dft_forward(complex_double *data, size_t N)
{
    complex_double *out = calloc(sizeof(complex_double), N);

    // dft_naive(data, out, N);
    fft_radix2_iter(data, out, N, 1);

    for (int i=0; i<N; i++)
    {
        printf("%lf + %lf * i \n", data[i].re, data[i].im);
        data[i] = out[i];
    }

    free(out);
}


void dft_backward(complex_double *data, size_t N)
{
    for (complex_double *p=data; p!=data+N; p++)
    {
        p->im *= -1;
    }
    dft_forward(data, N);
    for (complex_double *p=data; p!=data+N; p++)
    {
        p->im *= -1;
        *p = mul_by_factor(*p, 1./N);
    }
}
