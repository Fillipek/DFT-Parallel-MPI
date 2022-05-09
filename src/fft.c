#include "fft.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Discrette Fourier Transform naive implementation (from definition)
void dft_naive(double complex *in, double complex *out, size_t N)
{
    for (int k = 0; k < N; k++) {
        for (int n = 0; n < N; n++) {
            double complex temp =  cexp(-2 * M_PI / N * k * n * I);
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
        double complex t = cexp(-2 * M_PI * k / N * I);
        double complex q = t * out[k + N/2];
        out[k] = p + q;
        out[k+N/2] = p - q;
    }
}

// Implementacja z: https://www.geeksforgeeks.org/write-an-efficient-c-program-to-reverse-bits-of-a-number/: 
static int rev(unsigned int num, size_t N)
{
    int Nbits = (int)log2(N);

    int x_rev = 0;
    for (int i=0; i<Nbits; i++)
    {
        int curr_bit = (num >> i) & 1;
        x_rev += (curr_bit << (Nbits - 1 - i)); 
    }
        
    return x_rev;
}

static void bit_reversal_copy(double complex* result, double complex* input, size_t N)
{
    for(unsigned int k = 0; k < N; k++)
    {
        printf("%d %d\n", k, rev(k, N));
        result[rev(k, N)] = input[k];
    }
}


void fft_radix2_iter(double complex *in, double complex *out, size_t N, size_t stride)
{
    bit_reversal_copy(out, in, N);

    for(int s = 1; s <= log2(N); s++)
    {
        int m = 1 << s;
        double complex omega_m = cexp(-2. * M_PI / m * I);

        for(int k = 0; k <= (N - 1); k += m)
        {
            double complex omega = 1;

            for(int j = 0; j <= m / 2 - 1; j++)
            {
                double complex t = omega * out[k + j + m / 2];
                double complex u = out[k + j];
                out[k + j] = u + t;
                out[k + j + m / 2] = u - t;
                omega = omega * omega_m;
            }
        }
    }
}

void dft_forward(double complex *data, size_t N)
{
    double complex *out = calloc(sizeof(double complex), N);

    // dft_naive(data, out, N);
    fft_radix2_iter(data, out, N, 1);

    for (int i=0; i<N; i++)
    {
        printf("%lf + %lf * i \n", creal(data[i]), cimag(data[i]));
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
