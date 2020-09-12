/************************************************
 * CS532 - Lab 3
 * Patrick Dougherty
 * patrick.r.dougherty@und.edu
 * 12 Sep 20
 * *********************************************/


#include <iostream>
#include <mpi.h>
#include <sstream>


/*
 * Write a function that uses MPI to broadcast a message from rank 0 
 * to all other ranks in the current MPI_COMM_WORLD. Use only MPI_Send
 * and MPI_Recv in your implementation.
 */ 

void mybcast(void *array, int array_len, int my_rank, int comm_size) {
    if (my_rank == 0)
    {
        for (int i = 1; i < comm_size; i++)
        {
            MPI_Send(array, array_len, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(array, array_len, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void printArray(int *array, int array_len, int my_rank) {
    using std::cout;
    using std::endl;
    using std::ostringstream;

    ostringstream oss;
    oss << "[rank: " << my_rank << "]: ";
    for (int i = 0; i < array_len; i++) {
        oss << " " << array[i];
    }
    cout << oss.str() << endl;

    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char *argv[]) {
    int array_len = 20;
    int *array = new int[array_len];

    int comm_size;
    int my_rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) {
        for (int i = 0; i < array_len; i++) {
            array[i] = i;
        }
    } else {
        for (int i = 0; i < array_len; i++) {
            array[i] = 0;
        }
    }

    printArray(array, array_len, my_rank);
    mybcast(array, array_len, my_rank, comm_size);
    printArray(array, array_len, my_rank);

    MPI_Finalize();
    return 0;
}

