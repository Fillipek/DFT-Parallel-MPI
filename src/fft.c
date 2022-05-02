#include "fft.h"
#include <math.h>
#include <stdlib.h>

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

void dft_forward(complex_double *data, size_t N)
{
    complex_double *out = calloc(sizeof(complex_double), N);

    // dft_naive(data, out, N);
    fft_radix2(data, out, N, 1);

    for (int i=0; i<N; i++)
    {
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
