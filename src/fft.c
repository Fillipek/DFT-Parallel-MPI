#include "fft.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

void dft_naive(double complex *in, double complex *out, size_t N, MPI_Data mpi_data);
void fft_radix2(double complex *in, double complex *out, size_t N, size_t s, MPI_Data mpi_data);
static int rev(unsigned int num, size_t N);
static void bit_reversal_copy(double complex* result, double complex* input, size_t N, MPI_Data mpi_data);
void fft_radix2_iter(double complex *in, double complex *out, size_t N, MPI_Data mpi_data);

void fft_forward(double complex *data, size_t N, MPI_Data mpi_data)
{
    double complex *out = calloc(sizeof(double complex), N);

    // dft_naive(data, out, N, mpi_data);
    fft_radix2_iter(data, out, N, mpi_data);

    for (int i=0; i<N; i++)
    {
        data[i] = out[i];
    }

    free(out);
}


void fft_backward(double complex *data, size_t N, MPI_Data mpi_data)
{
    for (double complex *p=data; p!=data+N; p++)
    {
        *p = conj(*p);
    }
    fft_forward(data, N, mpi_data);
    for (double complex *p=data; p!=data+N; p++)
    {
        *p = conj(*p) / (double)N;
    }
}

void fft_filter(double threshold, double complex *data, size_t N, MPI_Data mpi_data)
{
    double module_max_sq = 0, module_curr_sq;
    for (double complex *p=data; p!=data+N; p++)
    {
        module_curr_sq = cabs(*p);
        if (module_curr_sq > module_max_sq)
        {
            module_max_sq = module_curr_sq;
        }
    }
    double th = threshold * module_max_sq;
    for (double complex *p=data; p!=data+N; p++)
    {
        module_curr_sq = cabs(*p);
        if (module_curr_sq < th)
        {
            *p = 0;
        }
    }
}

/* ------------------------------------------------------------------------- */

void dft_naive(double complex *in, double complex *out, size_t N, MPI_Data mpi_data)
{
    const int idx_per_proc = N / mpi_data.n_proc;

    int *recv_counts = malloc(sizeof(int) * mpi_data.n_proc);
    int *displacements = malloc(sizeof(int) * mpi_data.n_proc);
    for (int proc_rank=0; proc_rank< mpi_data.n_proc; proc_rank++)
    {
        recv_counts[proc_rank] = idx_per_proc;
        displacements[proc_rank] = idx_per_proc * proc_rank;
    }
    
    double complex *tmp = calloc(sizeof(double complex), idx_per_proc);
    for (int k = 0; k < idx_per_proc; k++) {
        int global_k = k + mpi_data.proc_rank * idx_per_proc;
        for (int n = 0; n < N; n++) {
            tmp[k] += in[n] * cexp(-2 * M_PI / N * global_k * n * I);
        }
        out[global_k] = tmp[k];
    }
    MPI_Allgatherv(tmp, idx_per_proc, MPI_C_DOUBLE_COMPLEX, out, recv_counts, displacements, MPI_C_DOUBLE_COMPLEX, mpi_data.comm);
    free(tmp);
}

void fft_radix2(double complex *in, double complex *out, size_t N, size_t s, MPI_Data mpi_data)
{
    if (N == 1)
    {
        *out = *in;
        return;
    }
    fft_radix2(in, out, N/2, 2*s, mpi_data);
    fft_radix2(in+s, out+N/2, N/2, 2*s, mpi_data);
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

static void bit_reversal_copy(double complex* result, double complex* input, size_t N, MPI_Data mpi_data)
{
    for(unsigned int k = 0; k < N; k++)
    {
        result[rev(k, N)] = input[k];
    }
}


void fft_radix2_iter(double complex *in, double complex *out, size_t N, MPI_Data mpi_data)
{
    bit_reversal_copy(out, in, N, mpi_data);

    const int idx_per_proc = N / mpi_data.n_proc;

    int *recv_counts = malloc(sizeof(int) * mpi_data.n_proc);
    int *displacements = malloc(sizeof(int) * mpi_data.n_proc);
    for (int proc_rank=0; proc_rank< mpi_data.n_proc; proc_rank++)
    {
        recv_counts[proc_rank] = idx_per_proc;
        displacements[proc_rank] = idx_per_proc * proc_rank;
    }
    double complex *tmp = malloc(sizeof(double complex) * idx_per_proc);

    for(int s = 1; s <= log2(N); s++)
    {
        int m = 1 << s;
        double complex omega_m = cexp(-2. * M_PI / m * I);

        if (s <= log2(N) - log2(mpi_data.n_proc))
        {
            // Parallel
            for (int k=0; k < idx_per_proc; k++)
            {
                int global_k = k + mpi_data.proc_rank * idx_per_proc;
                tmp[k] = out[global_k];
            }
            for(int k = 0; k < idx_per_proc; k += m)
            {
                double complex omega = 1;

                for(int j = 0; j < m/2; j++)
                {
                    double complex t = omega * tmp[k + j + m / 2];
                    double complex u = tmp[k + j];
                    tmp[k + j] = u + t;
                    tmp[k + j + m/2] = u - t;
                    omega *= omega_m;
                } 
            }
            MPI_Allgatherv(tmp, idx_per_proc, MPI_C_DOUBLE_COMPLEX, out, recv_counts, displacements, MPI_C_DOUBLE_COMPLEX, mpi_data.comm);
        }
        else //if (mpi_data.proc_rank == MPI_PROC_RANK_MASTER)
        {
            // Sequential
            for(int k = 0; k < N; k += m)
            {
                double complex omega = 1;

                for(int j = 0; j < m/2; j++)
                {
                    double complex t = omega * out[k + j + m / 2];
                    double complex u = out[k + j];
                    out[k + j] = u + t;
                    out[k + j + m/2] = u - t;
                    omega *= omega_m;
                } 
            }
        }
    }

    free(tmp);
}
