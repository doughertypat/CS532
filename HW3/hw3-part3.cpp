/************************************************
 * CSci 532 - HW3 Part 2
 * Patrick Dougherty
 * patrick.r.dougherty@und.edu
 * 1 Dec 2020
 * *********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

/************************************************
 * Functions
 * *********************************************/

/************************************************
 * Main
 * *********************************************/

int main(int argc, char *argv[]) {
    int print = 0;
    double start, end;
    
    if (argc < 3 || argc > 4) {
        printf("Improper usage: program <rows> <cols> [<print flag>]\n");
        return 0;
    }
    if (argc == 4) {
        if (!strcmp("y",argv[3]) || !strcmp("Y", argv[3])){
            print = 1;
        }
    }
    int rows = atoi(argv[1]);
    int cols = atoi(argv[2]);

    int *matrix = new int[rows * cols];
    int *matrix_T = new int[rows * cols]; 

    for (int i = 0; i < (rows * cols); i++) {
        matrix[i] = i;
    }

/*
    for (int i = 0; i < (rows * cols); i++) {
        if ((i+1)%cols) {
            printf("%d ", matrix[i]);
        }
        else {
            printf("%d\n", matrix[i]);
        }
    }
*/
    omp_set_num_threads(4);

    start = omp_get_wtime();
    #pragma omp parallel
    #pragma omp for
    for (int i = 0; i < (rows * cols); i++) {
        matrix_T[(i % cols) * rows + (i / cols)] = matrix[i];
    }
    end = omp_get_wtime();

    printf("Transpose time: %f\n", (end - start));

    if (print == 1) {
        for (int i = 0; i < (rows * cols); i++) {
            if ((i+1)%rows) {
                printf("%d ", matrix_T[i]);
            }
            else {
                printf("%d\n", matrix_T[i]);
            }
        }
    }

    delete[] matrix;
    delete[] matrix_T;
    return 0;
}
