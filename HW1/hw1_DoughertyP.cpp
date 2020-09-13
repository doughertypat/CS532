/************************************************
 * CS532 - HW1
 * Patrick Dougherty
 * patrick.r.dougherty@und.edu
 * 12 Sep 20
 * *********************************************/


#include <iostream>
#include <mpi.h>
#include <sstream>


/*
 * Write an MPI program with two tasks that alter and send a numeric message back and forth for
 * a set number of iterations. The root task shall add two to a number, print it to the console, send
 * the number to the other task, and receive the new value from the non-root task. The other task
 * shall receive the number, add 20 to it, and send the updated value back to the root task. The
 * number of iterations shall be input to the program with a command line argument. 
 */ 


void incrementNumber(int *data, int my_rank)
{
    if (my_rank == 0)
    {
       *data += 2;
    }
    else
    {
       *data += 20;
    }
}

int main(int argc, char *argv[]) {
    int value = 0;
    int iterations;
    int comm_size;
    int my_rank;
    
    if (argc != 2)
    {
        std::cout << "Improper usage: hw1_DoughertyP <# of iterations>" << std::endl;
        return 0;
    }
    iterations = atoi(argv[1]); 

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if(my_rank == 0)
    {
        std::cout << "Starting value: " << value << std::endl;
    }

    for (int i = 0; i < iterations; i++)
    {
        if(my_rank == 0)
        {
           incrementNumber(&value, my_rank);
           std::cout << "New value: " << value << std::endl;
           MPI_Send(&value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
           MPI_Recv(&value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            MPI_Recv(&value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            incrementNumber(&value, my_rank);
            MPI_Send(&value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); 
        }
    }

    if(my_rank == 0)
    {
        std::cout << "Final value: " << value << std::endl;
    }

    MPI_Finalize();
    
    return 0;
}

