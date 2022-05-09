#include "fft.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

// Discrette Fourier Transform naive implementation (from definition)
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

    if(mpi_data.proc_rank == MPI_PROC_RANK_MASTER)
    {
        for (int i=0; i< mpi_data.n_proc; i++)
        {
            printf("[%d] %d\n", i, recv_counts[i]);
        }
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


    // MPI_Barrier(mpi_data.comm);
    // MPI_Bcast(out+k, 1, MPI_C_DOUBLE_COMPLEX, mpi_data.proc_rank, mpi_data.comm);

    // if (mpi_data.proc_rank == MPI_PROC_RANK_MASTER)
    // {
    //     for(int k=0; k<N; k++)
    //     // for (int k = mpi_data.proc_rank; k < N; k += mpi_data.n_proc)
    //     {
    //         printf("[Debug] fft[%d] = %lf + i%lf\n", k, fabs(creal(out[k])), fabs(cimag(out[k])));
    //     }
    // }
}

// Fast Fourier Transform using recursive Cooley-Tukey algorithm (Radix-2)
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
    for(unsigned int k = mpi_data.proc_rank; k < N; k+=mpi_data.n_proc)
    {
        result[rev(k, N)] = input[k];
        MPI_Bcast(result+k, 1, MPI_C_DOUBLE_COMPLEX, mpi_data.proc_rank, mpi_data.comm);
    }
}


void fft_radix2_iter(double complex *in, double complex *out, size_t N, MPI_Data mpi_data)
{
    bit_reversal_copy(out, in, N, mpi_data);

    for(int s = 1; s <= log2(N); s++)
    {
        int m = 1 << s;
        double complex omega_m = cexp(-2. * M_PI / m * I);

        // for(int k = mpi_data.proc_rank; k < N; k += m * mpi_data.n_proc)
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
            // MPI_Bcast(out+k, m/2, MPI_C_DOUBLE_COMPLEX, mpi_data.proc_rank, mpi_data.comm);
            // MPI_Bcast(out+k+m/2, m/2, MPI_C_DOUBLE_COMPLEX, mpi_data.proc_rank, mpi_data.comm);
        }
    }
}

void dft_forward(double complex *data, size_t N, MPI_Data mpi_data)
{
    double complex *out = calloc(sizeof(double complex), N);

    dft_naive(data, out, N, mpi_data);
    // fft_radix2_iter(data, out, N, mpi_data);

    for (int i=0; i<N; i++)
    {
        // printf("%lf + %lf * i \n", creal(data[i]), cimag(data[i]));
        data[i] = out[i];
    }

    free(out);
}


void dft_backward(double complex *data, size_t N, MPI_Data mpi_data)
{
    for (double complex *p=data; p!=data+N; p++)
    {
        *p = conj(*p);
    }
    dft_forward(data, N, mpi_data);
    for (double complex *p=data; p!=data+N; p++)
    {
        *p = conj(*p) / (double)N;
    }
}
