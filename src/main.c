#include <stdio.h>
#include <stdlib.h>
#include "fft.h"
#include <complex.h>
#include "mpi.h"
#include "mpi_data.h"
#include "file_utils.h"

int main(int argc, char *argv[])
{
    /* Variable declarations */
    MPI_Data mpi_data = { .comm = MPI_COMM_WORLD };
    double complex *data;
    size_t N = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &(mpi_data.n_proc));
    MPI_Comm_rank(MPI_COMM_WORLD, &(mpi_data.proc_rank));

    if (mpi_data.proc_rank == MPI_PROC_RANK_MASTER)
    {
        printf("[Info] Starting with number of processes: %d\n", mpi_data.n_proc);
    }

    /* Data readout from file */
    if (mpi_data.proc_rank == MPI_PROC_RANK_MASTER)
    {
        if (argc <= 1)
        {
            fprintf(stderr, "ERROR: The program requires an argument which is a path to data file.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        count_lines_in_file(argv[1], &N);
    }

    MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG, MPI_PROC_RANK_MASTER, MPI_COMM_WORLD);
    data = malloc(sizeof(double complex) * N);

    if (mpi_data.proc_rank == MPI_PROC_RANK_MASTER)
    {
        read_data_lines(argv[1], N, data);
    }
    MPI_Bcast(data, N, MPI_C_DOUBLE_COMPLEX, MPI_PROC_RANK_MASTER, MPI_COMM_WORLD);

    fft_forward(data, N, mpi_data);

    if (mpi_data.proc_rank == MPI_PROC_RANK_MASTER)
    {
        save_to_file("data/fft_data.csv", N, data);
    }

    #ifdef FFT_FLTER_THRESHOLD
        fft_filter(FFT_FLTER_THRESHOLD, data, N, mpi_data);
    #endif

    fft_backward(data, N, mpi_data);

    if (mpi_data.proc_rank == MPI_PROC_RANK_MASTER)
    {
        save_to_file("data/rev_fft_data.csv", N, data);
    }

    /* Cleanup */
    free(data);
    MPI_Finalize();
}
