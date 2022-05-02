#include <stdio.h>
#include <stdlib.h>
#include "fft.h"
#include "complex_own.h"

int main(int argc, char *argv[])
{
    // inicjalizacja parametrów i zmiennych
    complex_double *f_vals;        // array of function values

    // wczytanie danych wejściowych
    if (argc <= 1)
    {
        fprintf(stderr, "ERROR: The program requires an argument which is a path to data file.\n");
        exit(EXIT_FAILURE);
    }

    FILE *fp;
    size_t len = 0;
    char *current_line = NULL;
    size_t line_len;
    size_t n_lines = 0;

    fp = fopen(argv[1], "r");
    while ((line_len = getline(&current_line, &len, fp)) != -1)
    {
        n_lines++;
    }

    fclose(fp);
    if (current_line)
        free(current_line);

    f_vals = malloc(sizeof(complex_double) * n_lines);

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

    // wypisanie funckji
    // funkcja licząca FFT

    dft_forward(f_vals, n_lines);


    // zachowanie wyników FFT
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
    dft_backward(f_vals, n_lines);

    // zachowanie wyników FFT-1
    fp = fopen("data/rev_fft_data.csv", "w");
    for (int i =0; i<n_lines; i++)
    {
        fprintf(fp, "%d, %.4lf, %.4lf\n", i, f_vals[i].re, f_vals[i].im);
    }
    fclose(fp);

    // szprzątanie
    free(f_vals);
}
