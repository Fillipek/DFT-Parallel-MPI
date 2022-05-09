#ifndef _MPI_DATA_H_
#define _MPI_DATA_H_

#define MPI_PROC_RANK_MASTER 0
#define MPI_TAG_DATA 0

typedef struct {
    int n_proc;
    int proc_rank;
    int comm;
} MPI_Data;

#endif