#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"
#include "complex_own.h"
#include "mpi.h"

#define PROC_MASTER 0
#define REQUEST 1

int main(int argc, char *argv[])
{
    /* ---- MPI CONTENT START ---- */
    MPI_Init(&argc, &argv);

    // komunikator globalny i zmienne do niego
    int myid = 0, numprocs = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    /* ---- MPI CONTENT END ---- */


    // inicjalizacja parametrów i zmiennych
    complex_double f_vals[128];        // array of function values

    // wczytanie danych wejściowych
    if (argc <= 1)
    {
        fprintf(stderr, "ERROR: The program requires an argument which is a path to data file.\n");
        exit(EXIT_FAILURE);
    }

    /* Niektore operacje powinien wykonac jedynie jeden proces */
    FILE *fp;
    size_t len = 0;
    char *current_line = NULL;
    size_t line_len;
    size_t n_lines = 0;

    if (myid == PROC_MASTER)
    {
        fp = fopen(argv[1], "r");
        while ((line_len = getline(&current_line, &len, fp)) != -1)
        {
            n_lines++;
        }

        fclose(fp);
        //if (current_line)
        //  free(current_line);
 
        //f_vals = malloc(sizeof(complex_double) * n_lines);
        //some_buffer = malloc(sizeof(int) * 100);
        
        fp = fopen(argv[1], "r");
        complex_double *f_curr = f_vals;
        while ((line_len = getline(&current_line, &len, fp)) != -1)
        {
            f_curr->re = atof(current_line);
            f_curr->im = 0;
            // printf("%.4lf + i%.4lf\n", f_curr->re, f_curr->im);
            f_curr++;
        }
        fclose(fp);
        if (current_line)
            free(current_line);
    }

    /* ----MPI CONTENT START---- */
    // Musimy rozpropagować do pozostałych procesów tablicę:
    MPI_Bcast(f_vals, 128 * sizeof(complex_double), MPI_CHAR, PROC_MASTER, MPI_COMM_WORLD);

    if(myid != 0)
    {   
        printf("%lf\n", f_vals[12].re);
        dft_forward(f_vals, 128, myid, 4);
        // Po obliczeniach trzeba wyslac wynik
        MPI_Send(f_vals, 128 * sizeof(complex_double), MPI_CHAR, PROC_MASTER, REQUEST, MPI_COMM_WORLD);
    }
    else
    {
        complex_double local_buffer[128];
        MPI_Status status;

        for(int source = 1; source <= 3; source++)
        {
            MPI_Recv(local_buffer, 128 * sizeof(complex_double), MPI_CHAR, source, REQUEST, MPI_COMM_WORLD, &status);

            //porownanie roznic:
            for(int i = 0; i < 128; ++i)
            {
                if(fabs(local_buffer[i].re - f_vals[i].re) >= 0.0001 || fabs(local_buffer[i].im - f_vals[i].im) >= 0.0001)
                {
                    f_vals[i] = local_buffer[i];
                }
            }
        }
    }

    /* ----MPI CONTENT END---- */

    // wypisanie funckji
    // funkcja licząca FFT

    // zachowanie wyników FFT
    if (myid == PROC_MASTER)
    {
        fp = fopen("data/fft_data.csv", "w");
        for (int i =0; i<n_lines; i++)
        {
            fprintf(fp, "%d, %.4lf, %.4lf\n", i, f_vals[i].re, f_vals[i].im);
        }
        fclose(fp);
    

        // filtracja sygnału (odcięcie wartości których moduł < 50% max moduł)
        // liczymy kwadraty, czyli granica to 25% max
        double module_max_sq = 0, module_curr_sq;
        for (complex_double *p=f_vals; p!=f_vals+n_lines; p++)
        {
            module_curr_sq = p->re * p->re + p->im * p->im;
            if (module_curr_sq > module_max_sq)
            {
                module_max_sq = module_curr_sq;
            }
        }
        double threshold = 0.25 * module_max_sq;
        for (complex_double *p=f_vals; p!=f_vals+n_lines; p++)
        {
            module_curr_sq = p->re * p->re + p->im * p->im;
            if (module_curr_sq < threshold)
            {
                p->re = 0;
                p->im = 0;
            }
        }
    
        // policzneie FFT-1
        dft_backward(f_vals, n_lines, myid, numprocs);


        fp = fopen("data/rev_fft_data.csv", "w");
        for (int i =0; i<n_lines; i++)
        {
            fprintf(fp, "%d, %.4lf, %.4lf\n", i, f_vals[i].re, f_vals[i].im);
        }
        
        fclose(fp);
    }
    
    MPI_Finalize();
}
