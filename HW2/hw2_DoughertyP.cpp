/************************************************
 * CSci 532 - HW2
 * Patrick Dougherty
 * patrick.r.dougherty@und.edu
 * 8 October 2020
 * *********************************************/

#include <stdio.h>
#include <mpi.h>
#include <chrono>
#include <stdlib.h>
#include <random>
#include <functional>
#include <cmath>

#define REPULSION   1.0f
#define ATTRACTION  0.1f
#define VELMATCH    0.05f
#define MINDIST     100.0f
#define MINV        -50.0f
#define MAXV        50.0f
#define MINRANGE    0.0f
#define MAXRANGE    1000.0f


float find_distance(float *b1, float *b2)
{
    float dx = b2[0] - b1[0];
    float dy = b2[1] - b1[1];
    float dz = b2[2] - b1[2];
    return sqrt(dx*dx + dy*dy + dz*dz);
}


void attraction(float *all_boid_pos, float *all_boid_vel, int *disp, int n_boids,
                int local_num_boids, int rank, int start_index)
{
    int boid;
    for (int j = 0; j < local_num_boids; j++)
    {
        float r1_vec[3] = { 0 }; 
        boid = j + start_index;
        //find perceived center
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
        //normalize perceived center
        r1_vec[0] /= (n_boids - 1);
        r1_vec[1] /= (n_boids - 1);
        r1_vec[2] /= (n_boids - 1);
        //update boid velocity 
        all_boid_vel[disp[rank] + (j*3)] += (r1_vec[0] - all_boid_pos[disp[rank] + (j*3)]) * ATTRACTION;
        all_boid_vel[disp[rank] + (j*3)+1] += (r1_vec[1] - all_boid_pos[disp[rank] + (j*3)+1]) * ATTRACTION;
        all_boid_vel[disp[rank] + (j*3)+2] += (r1_vec[2] - all_boid_pos[disp[rank] + (j*3)+2]) * ATTRACTION;   
    }
}

void repulsion(float *p, float *v, int disp, int n_boids, int local_num_boids, int rank, int start_index)
{
    int boid;
    float dist = 0;
    for (int i = 0; i < local_num_boids; i++)
    {
        float r2_vel[3] = { 0 };
        boid = i + start_index;
        //find boids that are too close and move away!
        for (int j = 0; j < n_boids; j++)
        {
            if (j == boid)
            {
                continue;
            }
            if (find_distance(&p[disp+(i*3)], &p[j*3]) < MINDIST)
            {
                r2_vel[0] -= (p[(j*3)] - p[disp + (i*3)]) * REPULSION;
                r2_vel[1] -= (p[(j*3)+1] - p[disp + (i*3)+1]) * REPULSION;
                r2_vel[2] -= (p[(j*3)+2] - p[disp + (i*3)+2]) * REPULSION;
            }
        }
        //update boid velocity
        v[disp + (i*3)] += r2_vel[0];
        v[disp + (i*3)+1] += r2_vel[1];
        v[disp + (i*3)+2] += r2_vel[2];
    }
}

void vel_match(float *v, int disp, int n_boids, int local_num_boids, int rank, int start_index)
{
    int boid;
    for (int i = 0; i < local_num_boids; i++)
    {
        float r3_vel[3] = { 0 };
        boid = i + start_index;
        //find perceived velocity
        for (int j = 0; j < n_boids; j++)
        {
            if (j == boid)
            {
                continue;
            }
            r3_vel[0] += v[(i*3)];
            r3_vel[1] += v[(i*3)+1];
            r3_vel[2] += v[(i*3)+2];
        }
        //normalize perceived velocity
        r3_vel[0] /= (n_boids - 1);
        r3_vel[1] /= (n_boids - 1);
        r3_vel[2] /= (n_boids - 1);
        //update boid velocity
        v[disp + (i*3)] += (r3_vel[0] - v[disp + (i*3)]) * VELMATCH;
        v[disp + (i*3)+1] += (r3_vel[1] - v[disp + (i*3)+1]) * VELMATCH;
        v[disp + (i*3)+2] += (r3_vel[2] - v[disp + (i*3)+2]) * VELMATCH;

    }
}
void initializeBoid(float *my_boid_pos, float *my_boid_vel, int local_num_boids)
{
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
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
}

void printBoids(float *b, float *v, int local_num_boids, int rank, int start_index)
{
    for (int i = 0; i < local_num_boids; i++)
    {
        printf("MPI Task %d: Boid %d pos:%f, %f, %f vel:%f, %f, %f\n", rank,
                                                                  i + start_index,
                                                                  b[(i*3)],
                                                                  b[(i*3)+1],
                                                                  b[(i*3)+2],
                                                                  v[(i*3)],
                                                                  v[(i*3)+1],
                                                                  v[(i*3)+2]);
    }
}

int main(int argc, char *argv[])
{
    int rank, size, num_boids, local_num_boids, start_index, num_iterations;
    float *all_boids_pos, *all_boids_vel;

          
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
        if (atoi(argv[1]) < size)
        {
            printf("Number of boids must be greater than, or equal to, the number of MPI tasks\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    num_boids = atoi(argv[1]);
    num_iterations = atoi(argv[2]);
    
    //assign boids to MPI tasks as evenly as possible
    assignBoids(rank, size, num_boids, recvcount, disp);

    all_boids_pos = (float *)malloc(num_boids*3*sizeof(float));
    all_boids_vel = (float *)malloc(num_boids*3*sizeof(float));
    float my_boid_pos[recvcount[rank]], my_boid_vel[recvcount[rank]];
    local_num_boids = recvcount[rank]/3; //used several times, reduce divisions
    start_index = disp[rank]/3;  //used several times, reduce divisions
    
    //initialize boids with random positions and 0 velocity
    initializeBoid(&all_boids_pos[disp[rank]], &all_boids_vel[disp[rank]], local_num_boids);
   
    if (rank == 0)
    {
        printf("You have selected %d boids and %d iterations with %d MPI tasks.\n", num_boids,
                                                                                    num_iterations,
                                                                                    size);
    }
    printf("MPI Task %d is caring for %d boids starting with boid %d.\n", rank,
                                                                          local_num_boids,
                                                                          start_index);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0)
    {
        printf("\n===========Set the boids free!===========\n\n");
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    printBoids(&all_boids_pos[disp[rank]], &all_boids_vel[disp[rank]], local_num_boids, rank, start_index);
    //exchange boid info
    MPI_Allgatherv(&all_boids_pos[disp[rank]], recvcount[rank], MPI_FLOAT, all_boids_pos,
                    recvcount, disp, MPI_FLOAT, MPI_COMM_WORLD); 
    MPI_Allgatherv(&all_boids_vel[disp[rank]], recvcount[rank], MPI_FLOAT, all_boids_vel,
                    recvcount, disp, MPI_FLOAT, MPI_COMM_WORLD);        
    
    for(int i = 0; i < num_iterations; i++)
    {
        //apply rules for velocity matching, attraction and repulsion amoung boids
        vel_match(all_boids_vel, disp[rank], num_boids, local_num_boids,
                   rank, start_index);
        attraction(all_boids_pos, all_boids_vel, disp, num_boids,
                   local_num_boids, rank, start_index);
        repulsion(all_boids_pos, all_boids_vel, disp[rank], num_boids,
                   local_num_boids, rank, start_index);

        for (int j = 0; j < local_num_boids; j++)
        {
            //enforce min/max velocities
            all_boids_vel[disp[rank]+(j*3)] = std::min(all_boids_vel[disp[rank]+(j*3)], MAXV);
            all_boids_vel[disp[rank]+(j*3)+1] = std::min(all_boids_vel[disp[rank]+(j*3)+1], MAXV);
            all_boids_vel[disp[rank]+(j*3)+2] = std::min(all_boids_vel[disp[rank]+(j*3)+2], MAXV);
            all_boids_vel[disp[rank]+(j*3)] = std::max(all_boids_vel[disp[rank]+(j*3)], MINV);
            all_boids_vel[disp[rank]+(j*3)+1] = std::max(all_boids_vel[disp[rank]+(j*3)+1], MINV);
            all_boids_vel[disp[rank]+(j*3)+2] = std::max(all_boids_vel[disp[rank]+(j*3)+2], MINV);
            //update positions
            all_boids_pos[disp[rank]+(j*3)] += all_boids_vel[disp[rank]+(j*3)];
            all_boids_pos[disp[rank]+(j*3)+1] += all_boids_vel[disp[rank]+(j*3)+1];
            all_boids_pos[disp[rank]+(j*3)+2] += all_boids_vel[disp[rank]+(j*3)+2];
            //enfore min/max positions
            all_boids_pos[disp[rank]+(j*3)] = std::min(all_boids_pos[disp[rank]+(j*3)], MAXRANGE);
            all_boids_pos[disp[rank]+(j*3)+1] = std::min(all_boids_pos[disp[rank]+(j*3)+1], MAXRANGE);
            all_boids_pos[disp[rank]+(j*3)+2] = std::min(all_boids_pos[disp[rank]+(j*3)+2], MAXRANGE);
            all_boids_pos[disp[rank]+(j*3)] = std::max(all_boids_pos[disp[rank]+(j*3)], MINRANGE);
            all_boids_pos[disp[rank]+(j*3)+1] = std::max(all_boids_pos[disp[rank]+(j*3)+1], MINRANGE);
            all_boids_pos[disp[rank]+(j*3)+2] = std::max(all_boids_pos[disp[rank]+(j*3)+2], MINRANGE);
        }
        //exchange boid info
        MPI_Allgatherv(&all_boids_pos[disp[rank]], recvcount[rank], MPI_FLOAT, all_boids_pos,
                    recvcount, disp, MPI_FLOAT, MPI_COMM_WORLD); 
        MPI_Allgatherv(&all_boids_vel[disp[rank]], recvcount[rank], MPI_FLOAT, all_boids_vel,
                    recvcount, disp, MPI_FLOAT, MPI_COMM_WORLD);        
    }
   
    MPI_Barrier(MPI_COMM_WORLD);

    printBoids(&all_boids_pos[disp[rank]], &all_boids_vel[disp[rank]], local_num_boids, rank, start_index);
    
    
    
    MPI_Finalize();
    return 0;
}




