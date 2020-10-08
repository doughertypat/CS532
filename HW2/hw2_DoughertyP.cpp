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

#define REPULSION   100
#define ATTRACTION  100
#define MINDIST     25
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

float *attraction(float *all_boid_pos, int n_boids, int boid)
{
    float r1_vec[3] = { 0 };
    for (int i = 0; i < n_boids; i++)
    {
        if (i == boid)
        {
            continue;
        }
        r1_vec[0] += all_boid_pos[(i*3)];
        r1_vec[1] += all_boid_pos[(i*3)+1];
        r1_vec[2] += all_boid_pos[(i*3)+2];
    }
    r1_vec[0] /= (n_boids - 1);
    r1_vec[1] /= (n_boids - 1);
    r1_vec[2] /= (n_boids - 1);
    
    r1_vec[0] = all_boid
           
    return r1_vec;   
}

void initializeBoid(float *my_boid_pos, float *my_boid_vel, int local_num_boids)
{
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    //std::mt19937 mt_rand(seed);
    auto somewhere = std::bind(std::uniform_real_distribution<float>(MINRANGE, MAXRANGE),
                     std::mt19937(seed));
    //auto somespeed = std::bind(std::uniform_real_distribution<float>(MINV, MAXV),
    //                 std::mt19937(seed));
    for (int i = 0; i < local_num_boids; i++)
    {
        my_boid_pos[(i*3)] = somewhere();
        my_boid_pos[(i*3)+1] = somewhere();
        my_boid_pos[(i*3)+2] = somewhere(); 
        my_boid_vel[(i*3)] = 0;
        my_boid_vel[(i*3)+1] = 0;
        my_boid_vel[(i*3)+2] = 0;
    }
}

void assignBoids(int rank, int size, int num_boids, int *recvcount, int *disp)
{
    
    if (num_boids <= size)
    {
        if (rank < num_boids)
        {
            recvcount[rank] = 1;
            disp[rank] = rank;
        }
        else
        {
            recvcount[rank] = 0;
        }
    }
    else //num_boids > size
    {
        int n_boids = num_boids;
        for (int i = 0; i < size; i++)
        {
            recvcount[i] = ceil((float)n_boids / (size - i));
            n_boids = n_boids - recvcount[i];
            recvcount[i] *= 3;
            
        }
        disp[0] = 0;
        for (int i = 1; i < size; i++)
        {
            disp[i] = disp[i-1] + recvcount[i-1];
        }
        //printf("Task %d, n_boids %d, size %d\n", rank, n_boids, size);  
    }
    for (int i = 0; i < size; i++)
    {
        //printf("Task %d, n_boids %d, disp %d\n", i, recvcount[i], disp[i]);
    }
}

void printBoids(float *b, int local_num_boids, int rank, int start_index)
{
    for (int i = 0; i < local_num_boids; i++)
    {
        printf("MPI Task %d: Boid %d - %f, %f, %f\n", rank,
                                                      i + start_index,
                                                      b[(i*3)],
                                                      b[(i*3)+1],
                                                      b[(i*3)+2]);
    }
}

int main(int argc, char *argv[])
{
    int rank, size, num_boids, local_num_boids, start_index, num_iterations;
    float *all_boids_pos;

          
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int recvcount[size], disp[size];
    
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

    assignBoids(rank, size, num_boids, recvcount, disp);

    all_boids_pos = (float *)malloc(num_boids*3*sizeof(float));
    float my_boid_pos[recvcount[rank]], my_boid_vel[recvcount[rank]];
    local_num_boids = recvcount[rank]/3;
    start_index = disp[rank]/3;

    if (recvcount[rank] > 0)
    {

        initializeBoid(my_boid_pos, my_boid_vel, local_num_boids);
        
        //printBoids(my_boid_pos, local_num_boids, rank, start_index);

        MPI_Allgatherv(my_boid_pos, recvcount[rank], MPI_FLOAT, all_boids_pos,
                        recvcount, disp, MPI_FLOAT, MPI_COMM_WORLD); 
                
        if(rank == 0)
        {
            printf("\n\n");
            for(int i = 0; i < num_boids; i++)
            {
                printf("Boid %d: %f, %f, %f\n", i, all_boids_pos[(3*i)],
                                                 all_boids_pos[(3*i)+1],
                                                 all_boids_pos[(3*i)+2]);
            }
        }
       
        float r1_vec[3];

        for (int i = 0; i < local_num_boids; i++)
        {
            r1_vec = attraction(all_boids_pos, num_boids, i + start_index);
            printf("%f, %f, %f\n 
       
        
        if(rank == 0)
        {
            printf("\n\n");
            for(int i = 0; i < num_boids; i++)
            {
                printf("Boid %d: %f, %f, %f\n", i, all_boids_pos[(3*i)],
                                                 all_boids_pos[(3*i)+1],
                                                 all_boids_pos[(3*i)+2]);
            }
        }

    }
    
    
    
    MPI_Finalize();
    return 0;
}




