#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define MAX_LIMIT 20
//#define N 10

int main(int argc, char **argv)
{
    int np;
    FILE *iptr;
    FILE *optr;
    int column1, row1, column2, row2;

    if (argc != 4)
    {
        printf("Usage: %s number_threads \n", argv[0]);
        exit(1);
    }

    np = atoi(argv[1]);

    if ((iptr = fopen(argv[2], "r")) == NULL)
    {
        printf("Error! opening file");
        // Program exits if the file pointer returns NULL.
        exit(1);
    }

    fscanf(iptr, "%d", &row1);
    fscanf(iptr, "%d", &column1);
    int *A = malloc(column1 * row1 * sizeof(int));

    for (int i = 0; i < row1 * column1; i++)
    {
        fscanf(iptr, "%d", &A[i]);
    }

    fscanf(iptr, "%d", &row2);
    fscanf(iptr, "%d", &column2);
    int *B = malloc(column2 * row2 * sizeof(int));

    for (int i = 0; i < row2 * column2; i++)
    {
        fscanf(iptr, "%d", &B[i]);
    }

    int row3 = row1;
    int column3 = column2;
    int *C = malloc(column3 * row3 * sizeof(int));

    int block = 2;
    double start = omp_get_wtime();

    omp_set_num_threads(np);

#pragma omp parallel
{
    #pragma omp for schedule(dynamic) nowait
    for (int ii = 1; ii <= row3 / block; ii++)
    {
        for (int kk = 1; kk <= column3 / block; kk++)
        {
            for (int jj = 1; jj <= column1 / block; jj++)
            {
                for (int i = (ii - 1) * block; i < (ii * block); i++)
                {
                    for (int k = (kk - 1) * block; k < (kk * block); k++)
                    {
                        for (int j = (jj - 1) * block; j < (jj * block); j++)
                        {
                            C[i * column3 + j] += A[i * column1 + k] * B[k * column2 + j];
                        }
                    }
                }
            }
        }
    }
}
    

    /*for (int i = 0; i < row3; i++)
    {
        for (int k = 0; k < column3; k++)
        {
            for (int j = 0; j < column1; j++)
            {
                C[i * column3 + j] += A[i * column1 + k] * B[k * column2 + j];
            }
        }
    }*/ 

    double end = omp_get_wtime();

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
            fprintf(optr, "%d\t", C[i * column3 + j]);
        }
        fprintf(optr, "\n");
    }

    fprintf(optr, "Running time: %es\n", end - start);

    free(A);
    free(B);
    free(C);
    fclose(iptr);
    fclose(optr);
    return 0;
}

/*for (int i = 0; i < row1; i++)
   {
       for (int j = 0; j < column1; j++)
       {
           printf("%d\t", A[i * column1 + j]);
       }
       printf("\n");
   }

   for (int i = 0; i < row2; i++)
   {
       for (int j = 0; j < column2; j++)
       {
           printf("%d\t", B[i * column2 + j]);
       }
       printf("\n");
   }*/

/*for (int i = 0; i < row3; i++)
{
    for (int j = 0; j < column3; j++)
    {
        C[i*column3 + j] = 0;
        for (int k = 0; k < column1; k++)
        {
            C[i*column3 + j] += A[i*column1 + k] * B[k*column2 + j];
        }
    }
}*/