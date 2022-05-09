#include <stdio.h>
#include <stdlib.h>
#include "fft.h"
#include <complex.h>
#include "mpi.h"

#define PROC_MASTER 0

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int myid = 0, numprocs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // inicjalizacja parametrów i zmiennych
    double complex *f_vals;        // array of function values
    FILE *fp;
    size_t N = 0;
    size_t len = 0;
    char *current_line = NULL;
    size_t line_len;

    if (myid == PROC_MASTER)
    {
        // wczytanie danych wejściowych
        if (argc <= 1)
        {
            fprintf(stderr, "ERROR: The program requires an argument which is a path to data file.\n");
            exit(EXIT_FAILURE);
        }

        fp = fopen(argv[1], "r");
        while ((line_len = getline(&current_line, &len, fp)) != -1)
        {
            N++;
        }

        fclose(fp);
        if (current_line)
            free(current_line);
    }

    MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG, PROC_MASTER, MPI_COMM_WORLD);
    f_vals = malloc(sizeof(double complex) * N);

    if (myid == PROC_MASTER)
    {
        fp = fopen(argv[1], "r");
        complex double *f_curr = f_vals;
        while ((line_len = getline(&current_line, &len, fp)) != -1)
        {
            *f_curr = (double)atof(current_line);
            f_curr++;
        }
        fclose(fp);
        if (current_line)
            free(current_line);
    }

    MPI_Bcast(f_vals, N, MPI_C_DOUBLE_COMPLEX, PROC_MASTER, MPI_COMM_WORLD);


    // wypisanie funckji
    // funkcja licząca FFT

    dft_forward(f_vals, N);


    // zachowanie wyników FFT
    if (myid == PROC_MASTER)
    {
        FILE *fp = fopen("data/fft_data.csv", "w");
        for (int i =0; i<N; i++)
        {
            fprintf(fp, "%d, %.4lf, %.4lf\n", i, creal(f_vals[i]), cimag(f_vals[i]));
        }
        fclose(fp);
    }

    // filtracja sygnału (odcięcie wartości których moduł < 50% max moduł)
    // liczymy kwadraty, czyli granica to 25% max
    // double module_max_sq = 0, module_curr_sq;
    // for (double complex *p=f_vals; p!=f_vals+N; p++)
    // {
    //     module_curr_sq = p->re * p->re + p->im * p->im;
    //     if (module_curr_sq > module_max_sq)
    //     {
    //         module_max_sq = module_curr_sq;
    //     }
    // }
    // double threshold = 0.25 * module_max_sq;
    // for (double complex *p=f_vals; p!=f_vals+N; p++)
    // {
    //     module_curr_sq = p->re * p->re + p->im * p->im;
    //     if (module_curr_sq < threshold)
    //     {
    //         *p = 0;
    //     }
    // }

    // policzneie FFT-1
    dft_backward(f_vals, N);

    // zachowanie wyników FFT-1
    if (myid == PROC_MASTER)
    {
        FILE *fp = fopen("data/rev_fft_data.csv", "w");
        for (int i =0; i<N; i++)
        {
            fprintf(fp, "%d, %.4lf, %.4lf\n", i, creal(f_vals[i]), cimag(f_vals[i]));
        }
        fclose(fp);
    }

    // szprzątanie
    free(f_vals);

    MPI_Finalize();
}
