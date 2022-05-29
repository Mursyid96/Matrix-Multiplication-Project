#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define MAX_LIMIT 20
#include <pthread.h>
#include <math.h>
#include <stdbool.h>
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
    int my_col, my_row;
    int get_col, get_row;
    int **localA, **localB, **localC;
    int dimension = N / squaredP;
    malloc2dint(&localA, dimension, dimension);
    malloc2dint(&localB, dimension, dimension);
    malloc2dint(&localC, dimension, dimension);

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

    // Cannon's algorithm
    // left-circular-shifts
    int row_original = 0;
    int row_shifted;
    int offset = 0;
    int min = 0;
    int max = min + (squaredP - 1);
    for (int i = 0; i < squaredP; i++)
    {
        for (int j = 0; j < squaredP; j++)
        {
            if (row_original + offset > max)
            {
                row_shifted = (row_original + (offset)) - squaredP;
            }
            else
            {
                row_shifted = (row_original + offset);
            }
            // if rank=original, need to receive new one from row_shifted
            if (my_rank == row_original)
            {
                if (my_rank < row_shifted)
                {
                    get_col = my_col + (abs((row_shifted - my_rank)) * dimension);
                }
                else
                {
                    get_col = my_col - (abs((row_shifted - my_rank)) * dimension);
                }
            }
            row_original++;
        }
        offset++;
        min = row_original;
        max = min + (squaredP - 1);
    }

    offset = 0;
    int col_original = 0;
    int col_shifted;
    min = 0;
    max = min + (squaredP * (squaredP - 1));
    for (int i = 0; i < squaredP; i++)
    {
        for (int j = 0; j < squaredP; j++)
        {
            if (col_original + (offset * squaredP) > max)
            {
                col_shifted = (col_original + (offset * squaredP)) - nthrd;
            }
            else
            {
                col_shifted = (col_original + (offset * squaredP));
            }
            if (my_rank == col_original)
            {
                if (my_rank < col_shifted)
                {
                    get_row = my_row + (abs((col_shifted - my_rank) / squaredP) * dimension);
                }
                else
                {
                    get_row = my_row - (abs((col_shifted - my_rank) / squaredP) * dimension);
                }
            }
            col_original += squaredP;
        }
        offset++;
        col_original = (col_original + 1) % squaredP;
        min = col_original;
        max = min + (squaredP * (squaredP - 1));
    }

    int rowA = my_row;
    int colA = get_col;
    int rowB = get_row;
    int colB = my_col;

    offset = 1;
    for (int q = 0; q < squaredP; q++)
    {
        // Calculate C
        // 1.transfer to local A&B
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
        // 2.calculate local submatrix C'
        int f=my_row;
        int g=my_col;
        for (int i = 0; i < dimension; i++)
        {
            for (int k = 0; k < dimension; k++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    // Write to main C
                    C[f][g] += localA[i][k] * localB[k][j];
                    g++;
                }
                g = my_col;
            }
            f++;
            g = my_col; 
        }
        
        get_col -= offset*dimension;
        if (get_col <= -1)
        {
            get_col += offset*N;
        }
        get_row -=offset*dimension;
        if (get_row <= -1)
        {
            get_row += offset*N;
        }
        colA=get_col;
        rowA=my_row;
        rowB=get_row;
        colB=my_col;
        //printf("my rank %d, matrix A that i got: (%d,%d)\n", my_rank, rowA, colA);
        //printf("my rank %d, matrix B that i got: (%d,%d)\n", my_rank, rowB, colB);
    }

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

    // Set up threads and number of threads
    nthrd = atoi(argv[1]);
    threads = (pthread_t *)malloc(nthrd * sizeof(*threads));
    squaredP = sqrt(nthrd);

    // Open file
    if ((iptr = fopen(argv[2], "r")) == NULL)
    {
        printf("Error! opening file");
        // Program exits if the file pointer returns NULL.
        exit(1);
    }

    // Set up and reads Matrix A
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

    // Set up and reads Matrix B
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

    // Set up and reads Matrix C
    malloc2dint(&C, column2, row1);
    int row3 = row1;
    int column3 = column2;
    N = row1;

    // Spawn threads
    GET_TIME (start);
    for (long i = 0; i < nthrd; i++)
    {
        pthread_create(&threads[i], NULL, mat_mul, (void *)i);
    }
    // reap threads
    for (int i = 0; i < nthrd; i++)
    {
        pthread_join(threads[i], NULL);
    }
    GET_TIME (finish);
    
    // Calculate running time
    elapsed = finish - start;

    // Opening output file
    optr = fopen(argv[3], "w");
    if (optr == NULL)
    {
        printf("Error!");
        exit(1);
    }

    // Writing to output file
    for (int i = 0; i < row3; i++)
    {
        for (int j = 0; j < column3; j++)
        {
            fprintf(optr, "%d\t", C[i][j]);
        }
        fprintf(optr, "\n");
    }

    // Prints running time
    fprintf(optr, "Running time: %es\n", elapsed);
    printf("Running time: %es\n", elapsed);

    free(A);
    free(B);
    free(C);
    fclose(iptr);
    fclose(optr);
    return 0;
}

