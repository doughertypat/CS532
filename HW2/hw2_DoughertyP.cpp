/************************************************
 * CSci 532 - HW2
 * Patrick Dougherty
 * patrick.r.dougherty@und.edu
 * 6 October 2020
 * *********************************************/

#include <stdio.h>
#include <mpi.h>
#include <chrono>
#include <stdlib.h>
#include <random>
#include <functional>
#include <cmath>

#define REPULSION   0
#define ATTRACTION  0
#define MINDIST     0
#define MINV        -50
#define MAXV        50
#define MINRANGE    0
#define MAXRANGE    1000

struct Vector {
    float x;
    float y;
    float z;
};

struct Boid {
    Vector pos;
    Vector vel;
    /*
    float x_pos;
    float y_pos;
    float z_pos;
    float x_vel;
    float y_vel;
    float z_vel;
    */
};


bool is_in_area(struct Boid *b1, struct Boid *b2)
{
    return true;
}

struct Vector *centerOfMass(struct Boid **boids)
{
    return NULL;   
}

void initializeBoid(float *my_boid_pos, int local_num_boids)
{
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    //std::mt19937 mt_rand(seed);
    auto somewhere = std::bind(std::uniform_real_distribution<float>(MINRANGE, MAXRANGE),
                     std::mt19937(seed));
    //auto somespeed = std::bind(std::uniform_real_distribution<float>(MINV, MAXV),
    //                 std::mt19937(seed));
    b->pos.x = somewhere();
    b->pos.y = somewhere();
    b->pos.z = somewhere(); 
    //b->vel.x = somespeed();
    //b->vel.y = somespeed();
    //b->vel.z = somespeed();
    b->vel.x = 0;
    b->vel.y = 0;
    b->vel.z = 0;
}

void assignBoids(int rank, int size, int num_boids, int *my_boids, int *start_index)
{
    if (num_boids <= size)
    {
        if (rank < num_boids)
        {
            *my_boids = 1;
            *start_index = rank;
        }
        else
        {
            *my_boids = 0;
        }
    }
    else //num_boids > size
    {
        int n_boids = num_boids;
        for (int i = 0; i < rank; i++)
        {
            n_boids = n_boids - ceil((float)n_boids / (size-i));
        }
        //printf("Task %d, n_boids %d, size %d\n", rank, n_boids, size);  
        *my_boids = ceil((float)n_boids / (size-rank));
        *start_index = num_boids - n_boids;
    }
    //printf("Task %d has %d boids assigned starting at %d, %d remaining\n",
    // rank, *my_boids, *start_index, num_boids);
}

void printBoids(struct Boid *b, int local_num_boids, int rank, int start_index)
{
    for (int i = 0; i < local_num_boids; i++)
    {
        printf("MPI Task %d: Boid %d - %f, %f, %f\n", rank,
                                                      i + start_index,
                                                      b[i].pos.x,
                                                      b[i].pos.y,
                                                      b[i].pos.z);
    }
}

int main(int argc, char *argv[])
{
    int rank, size, num_boids, local_num_boids, start_index, num_iterations;
    float *all_boids_pos;

          
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if(rank == 0)
    {
        if (argc < 3)
        {
            printf("Improper usage: program [# of boids] [# of iterations]\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    num_boids = atoi(argv[1]);
    num_iterations = atoi(argv[2]);

    assignBoids(rank, size, num_boids, &local_num_boids, &start_index);

    all_boids_pos = (float *)malloc(num_boids*3*sizeof(float));
    float my_send_array[local_num_boids * 3];
    
    if (local_num_boids > 0)
    {
        struct Boid my_boids[local_num_boids];

        for (int i = 0; i < local_num_boids; i++)
        {
            initializeBoid(&my_boids[i]);
            my_send_array[i] = my_boids[i].pos.x;
            my_send_array[i+1] = my_boids[i].pos.y;
            my_send_array[i+2] = my_boids[i].pos.z;           
              
        }
        printBoids(my_boids, local_num_boids, rank, start_index);
/*
        MPI_Allgather(my_send_array, local_num_boids*3, MPI_FLOAT, all_boids_pos,
                        num_boids*3, MPI_FLOAT, MPI_COMM_WORLD); 
                
        if(rank == 0)
        {
            printf("\n\n");
            for(int i = 0; i < num_boids; i++)
            {
                printf("Boid %d: %f, %f, %f", i, all_boids_pos[(3*i)],
                                                 all_boids_pos[(3*i)+1],
                                                 all_boids_pos[(3*i)+2]);
            }
        }
        */
        
    }
    
    
    
    MPI_Finalize();
    return 0;
}




