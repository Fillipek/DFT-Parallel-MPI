#include "file_utils.h"

void count_lines_in_file(const char *filename, size_t *N)
{
    char *current_line = NULL;
    size_t len = 0;
    FILE *fp = fopen(filename, "r");

    while (getline(&current_line, &len, fp) != -1)
    {
        (*N)++;
    }

    fclose(fp);
    if (current_line)
        free(current_line);
}


void read_data_lines(const char *filename, size_t N, double complex *data_buff)
{
    char *current_line = NULL;
    size_t len = 0;
    FILE *fp = fopen(filename, "r");
    double complex *curr_pos = data_buff;

    for (int i =0; i < N; i++)
    {
        getline(&current_line, &len, fp);
        *curr_pos++ = (double)atof(current_line);
    }
    fclose(fp);
    if (current_line)
        free(current_line);
}



void save_to_file(const char* filename, size_t N, double complex *data)
{
    FILE *fp = fopen(filename, "w");
    for (int i =0; i<N; i++)
    {
        fprintf(fp, "%d, %.4lf, %.4lf\n", i, creal(data[i]), cimag(data[i]));
    }
    fclose(fp);
}
