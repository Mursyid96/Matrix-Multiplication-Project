#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define MAX_LIMIT 20
#include <pthread.h>
#include <math.h>
#include <time.h>
#include "timer.h"
//#define N 10

int **A;
int **B;
int **C;
int N, nthrd, squaredP;

int malloc2dint(int ***array, int n, int m)
{
    // allocate the n*m contiguous items
    int *p = (int *)malloc(n * m * sizeof(int));
    if (!p)
        return -1;
    // allocate the row pointers into the memory */
    (*array) = (int **)malloc(n * sizeof(int *));
    if (!(*array))
    {
        free(p);
        return -1;
    }
    // set up the pointers into the contiguous memory
    for (int i = 0; i < n; i++)
        (*array)[i] = &(p[i * m]);
    return 0;
}

int free2dint(int ***array)
{
    free(&((*array)[0][0]));
    free(*array);
    return 0;
}

void *mat_mul(void *rank)
{
    long my_rank = (long)rank;
    int turn = 0;
    int my_col, my_row, my_col_end, my_row_end;
    int get_col, get_row;
    int **localA, **localB, **localC;
    int dimension = N / squaredP;
    malloc2dint(&localA, dimension, dimension);
    malloc2dint(&localB, dimension, dimension);
    malloc2dint(&localC, dimension, dimension);

    // Divide subarray pointer to all processors
    for (int i = 0; i < N; i += dimension)
    {
        // printf("%d\n",my_rank);
        for (int j = 0; j < N; j += dimension)
        {
            if (my_rank == turn)
            {
                my_row = i;
                my_col = j;
            }
            turn += 1;
        }
    }
    my_col_end = my_col + dimension;
    my_row_end = my_row + dimension;

    // SUMMA algorithm
    int column_procs_sender = 0;
    int row_procs_sender;
    int col_receiver = 0;
    int row_receiver = 0;
    for (int k = 0; k <= (squaredP - 1); k++)
    {
        row_procs_sender = k;
        for (int i = 0; i <= (squaredP - 1); i++)
        {
            //  broadcast submatrix A to processors in the same row
            for (int j = 0; j <= (squaredP - 1); j++)
            {
                if (my_rank == row_receiver)
                {
                    if (my_rank < row_procs_sender)
                    {
                        get_col = my_col + (abs((row_procs_sender - my_rank))*dimension);
                    }
                    else
                    {
                        get_col = my_col - (abs((row_procs_sender - my_rank))*dimension);
                    }
                }
                row_receiver++;
            }
            row_procs_sender += squaredP;
        }
        for (int j = 0; j <= (squaredP - 1); j++)
        {
            //  broadcast submatrix B to processors in the same column
            for (int l = 0; l <= (squaredP - 1); l++)
            {
                if (my_rank == col_receiver)
                {
                    if (my_rank < column_procs_sender)
                    {
                        get_row = my_row + (abs((column_procs_sender - my_rank)/squaredP)*dimension);
                    }
                    else
                    {
                        get_row = my_row - (abs((column_procs_sender - my_rank)/squaredP)*dimension);
                    }
                }
                col_receiver += squaredP;
            }
            column_procs_sender++;
            col_receiver = (col_receiver + 1) % squaredP;
            // printf("%d/n",col_receiver);
        }
        //printf("my rank %d, matrix A that i got: (%d,%d)\n",my_rank,my_row,get_col);
        //printf("my rank %d, matrix B that i got: (%d,%d)\n",my_rank,get_row,my_col);

        int rowA = my_row;
        int colA = get_col;
        int rowB = get_row;
        int colB = my_col;

        // transfer to local A&B
        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                localA[i][j] = A[rowA][colA];
                localB[i][j] = B[rowB][colB];
                colA++;
                colB++;
            }
            rowA++;
            colA=get_col;
            rowB++;
            colB=my_col;
        }

        // calculate local submatrix C'
        for (int i = 0; i < dimension; i++)
        {
            for (int k = 0; k < dimension; k++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    localC[i][j] += localA[i][k] * localB[k][j];
                } 
            }
        }

        col_receiver = 0;
        row_receiver = 0;
    }

    //transfer result to global
    int f = 0;
    int g = 0;
    for (int i = my_row; i < my_row_end; i++)
    {
        for (int j = my_col; j < my_col_end; j++)
        {
            C[i][j] += localC[f][g];
            g++;
        }
        f++;
        g = 0;
    }

    // printf ("my rank %d, my row %d, my col %d\n",my_rank,my_row,my_col);

    return 0;
}

int main(int argc, char **argv)
{
    pthread_t *threads;
    FILE *iptr;
    FILE *optr;
    int column1, row1, column2, row2;
    double start,finish,elapsed;

    if (argc != 4)
    {
        printf("Usage: %s number_threads \n", argv[0]);
        exit(1);
    }

    nthrd = atoi(argv[1]);
    threads = (pthread_t *)malloc(nthrd * sizeof(*threads));
    squaredP = sqrt(nthrd);

    if ((iptr = fopen(argv[2], "r")) == NULL)
    {
        printf("Error! opening file");
        // Program exits if the file pointer returns NULL.
        exit(1);
    }

    fscanf(iptr, "%d", &row1);
    fscanf(iptr, "%d", &column1);
    malloc2dint(&A, row1, column1);
    for (int i = 0; i < row1; i++)
    {
        for (int j = 0; j < column1; j++)
        {
            fscanf(iptr, "%d\t", &A[i][j]);
        }
    }
    fscanf(iptr, "%d", &row2);
    fscanf(iptr, "%d", &column2);
    malloc2dint(&B, row2, column2);
    for (int i = 0; i < row2; i++)
    {
        for (int j = 0; j < column2; j++)
        {
            fscanf(iptr, "%d\t", &B[i][j]);
        }
    }
    malloc2dint(&C, column2, row1);

    int row3 = row1;
    int column3 = column2;
    N = row1;

    // Spawn thread
    GET_TIME(start);
    for (long i = 0; i < nthrd; i++)
    {
        pthread_create(&threads[i], NULL, mat_mul, (void *)i);
    }

    for (int i = 0; i < nthrd; i++)
    {
        pthread_join(threads[i], NULL);
    }
    GET_TIME(finish);
    elapsed = finish - start;

    optr = fopen(argv[3], "w");

    if (optr == NULL)
    {
        printf("Error!");
        exit(1);
    }

    for (int i = 0; i < row3; i++)
    {
        for (int j = 0; j < column3; j++)
        {
            fprintf(optr, "%d\t", C[i][j]);
        }
        fprintf(optr, "\n");
    }

    fprintf(optr, "Running time: %es\n", elapsed);
    printf("Running time: %es\n", elapsed);

    free(A);
    free(B);
    free(C);
    fclose(iptr);
    fclose(optr);
    return 0;
}

