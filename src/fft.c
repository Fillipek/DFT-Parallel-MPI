#include "fft.h"
#include <math.h>
#include <stdlib.h>

// Discrette Fourier Transform naive implementation (from definition)
void dft_naive(double complex *in, double complex *out, size_t N)
{
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            double complex temp = cos(-2 * M_PI / N * k * n) + I*sin(-2 * M_PI / N * k * n);
            out[k] = out[k] + in[n] * temp;
        }
    }
}

// Fast Fourier Transform using recursive Cooley-Tukey algorithm (Radix-2)
void fft_radix2(double complex *in, double complex *out, size_t N, size_t s)
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
        double complex p = out[k];
        double complex t = cos(-2 * M_PI * k / N) + I*sin(-2 * M_PI * k / N);
        double complex q = t * out[k + N/2];
        out[k] = p + q;
        out[k+N/2] = p - q;
    }
}

void dft_forward(double complex *data, size_t N)
{
    double complex *out = calloc(sizeof(double complex), N);

    // dft_naive(data, out, N);
    fft_radix2(data, out, N, 1);

    for (int i=0; i<N; i++)
    {
        data[i] = out[i];
    }

    free(out);
}


void dft_backward(double complex *data, size_t N)
{
    for (double complex *p=data; p!=data+N; p++)
    {
        *p = conj(*p);
    }
    dft_forward(data, N);
    for (double complex *p=data; p!=data+N; p++)
    {
        *p = conj(*p) / (double)N;
    }
}
