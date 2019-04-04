/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Team Names Here              **(*****************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include "clcg4.h"

#include<mpi.h>
#include<pthread.h>

// #define BGQ 1 // when running BG/Q, comment out when testing on mastiff

#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD  0

#define SIZE 8
#define NUM_THREADS 2
#define NUM_GENERATIONS 256
#define THRESHOLD .25

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

int mpi_myrank;
int mpi_commsize;
int aliveCount[NUM_GENERATIONS];
int ** myUniverse;
int * topGhost;
int * bottomGhost;
pthread_barrier_t barrier;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
//long ** seeds;


// You define these


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these

void computeGeneration(int id, int tick){
    int updatedTick[SIZE/mpi_commsize/NUM_THREADS][SIZE];
    int g; /* get global index */
    int row;
    int aliveCounter = 0;
    int aliveNeigh = 0;
    int rowsPerRank = SIZE / mpi_commsize;
    int rowsPerThread = rowsPerRank / NUM_THREADS;
    double rngVal = 0;
    int randVal = 0;

    for(int i = 0; i < rowsPerThread; i++){
        g = (id * rowsPerThread) + (rowsPerRank * mpi_myrank) + i;
        row = (id * rowsPerThread) + i;

        for(int j = 0; j < SIZE; j++){ /* loop through all 32k cells */
            rngVal = GenVal(g); /* generate a random value [0, 1] for each cell */

            // if RNG < threshold, random alive/dead
            if(rngVal < THRESHOLD) {
                //random pick between LIVE or DEAD
                randVal = rand() % 2;
                updatedTick[i][j] = randVal;
                //update live count for each gen
                aliveCounter += randVal;
            } else {
                // side neighbors
                aliveNeigh += myUniverse[row][(j - 1) % SIZE];
                aliveNeigh += myUniverse[row][(j + 1) % SIZE];
                // top neighbors
                if(id == 0 && i == 0){ /* Needs to access top ghost row */
                    aliveNeigh += topGhost[(j - 1) % SIZE];
                    aliveNeigh += topGhost[(j + 1) % SIZE];
                    aliveNeigh += topGhost[j % SIZE];
                }
                else{
                    aliveNeigh += myUniverse[row - 1][(j - 1) % SIZE];
                    aliveNeigh += myUniverse[row - 1][(j + 1) % SIZE];
                    aliveNeigh += myUniverse[row - 1][j % SIZE];
                }
                // bottom neighbors
                if(id == NUM_THREADS-1 && i == rowsPerThread-1){ /* Needs to access bottom ghost row */
                    aliveNeigh += bottomGhost[(j - 1) % SIZE];
                    aliveNeigh += bottomGhost[(j + 1) % SIZE];
                    aliveNeigh += bottomGhost[j % SIZE];
                }
                else{
                    aliveNeigh += myUniverse[row + 1][(j - 1) % SIZE];
                    aliveNeigh += myUniverse[row + 1][(j + 1) % SIZE];
                    aliveNeigh += myUniverse[row + 1][j % SIZE];
                }
                
                if(myUniverse[i][j] == ALIVE) {
                  if(aliveNeigh < 2 || aliveNeigh > 3) {
                    updatedTick[i][j] = DEAD;
                  }
                  else{
                    updatedTick[i][j] = ALIVE;
                  }
                  // else still alive
                } 
                // if dead
                else {
                    // bring back to life if exactly 3 neighbors
                    if(aliveNeigh == 3) {
                        updatedTick[i][j] = ALIVE;
                    }
                    else{
                        updatedTick[i][j] = DEAD;
                    }

                }

                aliveCounter += updatedTick[i][j];
            }
        }
    }

    pthread_mutex_lock(&mutex);
    aliveCount[tick] += aliveCounter;
    pthread_mutex_unlock(&mutex);
    pthread_barrier_wait(&barrier);

    for(int i = 0; i < rowsPerThread; i++) {
      for(int j = 0; j < SIZE; j++) {
        myUniverse[id * rowsPerThread + i][j] = updatedTick[i][j];
      }
    }
}

// function for threads
void *conways(void * threadID){
    int id = *(int*) threadID;

    // loop over all generations
    for(int tick = 0; tick < NUM_GENERATIONS; tick++){
        MPI_Barrier(MPI_COMM_WORLD);
        pthread_barrier_wait(&barrier);
        computeGeneration(id, tick);
    }
    pthread_exit(NULL);
}

// function for main thread
void main_conways(){
    MPI_Request sendTop;
    MPI_Request sendBot;
    MPI_Request recvTop;
    MPI_Request recvBot;
    int mpi_flag;

    for(int tick = 0; tick < NUM_GENERATIONS; tick++){
        MPI_Barrier(MPI_COMM_WORLD);
        // send first row up
        MPI_Isend(myUniverse[0], SIZE, MPI_INT, (mpi_myrank-1)%mpi_commsize, 0, MPI_COMM_WORLD, &sendTop);
        // send last row down
        MPI_Isend(myUniverse[SIZE/mpi_commsize-1], SIZE, MPI_INT, (mpi_myrank+1)%mpi_commsize, 0, MPI_COMM_WORLD, &sendBot);
        // receive top ghost
        MPI_Irecv(topGhost, SIZE, MPI_INT, (mpi_myrank-1)%mpi_commsize, 0, MPI_COMM_WORLD, &recvTop);
        // receive bottom ghost
        MPI_Irecv(bottomGhost, SIZE, MPI_INT, (mpi_myrank+1)%mpi_commsize, 0, MPI_COMM_WORLD, &recvBot);

        MPI_Wait(&sendTop, MPI_STATUS_IGNORE);
        MPI_Test(&sendTop, &mpi_flag, MPI_STATUS_IGNORE);
        if(!mpi_flag){
            fprintf(stderr, "MPI error on send top\n");
            exit(1);
        }

        MPI_Wait(&sendBot, MPI_STATUS_IGNORE);
        MPI_Test(&sendBot, &mpi_flag, MPI_STATUS_IGNORE);
        if(!mpi_flag){
            fprintf(stderr, "MPI error on send bottom\n");
            exit(1);
        }

        MPI_Wait(&recvTop, MPI_STATUS_IGNORE);
        MPI_Test(&recvTop, &mpi_flag, MPI_STATUS_IGNORE);
        if(!mpi_flag){
            fprintf(stderr, "MPI error on recv top\n");
            exit(1);
        }


        MPI_Wait(&recvBot, MPI_STATUS_IGNORE);
        MPI_Test(&recvBot, &mpi_flag, MPI_STATUS_IGNORE);
        if(!mpi_flag){
            fprintf(stderr, "MPI error on recv bot\n");
            exit(1);
        }

        // synchronize threads
        pthread_barrier_wait(&barrier);
        computeGeneration(0, tick);
    }
}

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{

// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

// Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();

// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.
    //printf("Rank %d of %d has been started and a first Random Value of %lf\n",
	//   mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    MPI_Barrier( MPI_COMM_WORLD );

// Insert your code
    pthread_t my_threads[NUM_THREADS-1];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_barrier_init(&barrier, NULL, NUM_THREADS);

    int *totalAliveCount;
    if(mpi_myrank == 0){
        totalAliveCount = calloc(NUM_GENERATIONS, sizeof(int));
    }

    for(int i = 0; i < NUM_GENERATIONS; ++i) {
      aliveCount[i] = 0;
    }

    //initialize universe
    myUniverse = calloc(SIZE/mpi_commsize, sizeof(int*));
    for(int i = 0; i < SIZE/mpi_commsize; i++){
        myUniverse[i] = calloc(SIZE, sizeof(int));
        for(int j = 0; j < SIZE; j++)
            myUniverse[i][j] = ALIVE;
    }
    //init ghost rows
    topGhost = calloc(SIZE, sizeof(int));
    bottomGhost = calloc(SIZE, sizeof(int));

    // begin timer
    if(mpi_myrank == 0)
        g_start_cycles = GetTimeBase();


    // create threads
    for(int i = 0; i < NUM_THREADS-1; i++){
        int id = i+1;
        pthread_create(&my_threads[i], &attr, conways, (void *) &id);
    }
    main_conways();

    // Perform a barrier
    for(int i = 0; i < NUM_THREADS-1; i++){
        pthread_join(my_threads[i], NULL);
    }

    // Perform mpi_reduce on the alive count array
    if(mpi_myrank == 0)
        MPI_Reduce(
            aliveCount,
            totalAliveCount,
            NUM_GENERATIONS,
            MPI_INT,
            MPI_SUM,
            0,
            MPI_COMM_WORLD
            );
    else
        MPI_Reduce(
            aliveCount,
            NULL,
            NUM_GENERATIONS,
            MPI_INT,
            MPI_SUM,
            0,
            MPI_COMM_WORLD
            );


    if(mpi_myrank == 0){

        for(int i = 0; i < NUM_GENERATIONS; i++){
            printf("Generation %d: %d alive\n", i, totalAliveCount[i]);
        }

        g_end_cycles = GetTimeBase();
        g_time_in_secs = ((double)(g_end_cycles - g_start_cycles)) / g_processor_frequency;
    }


    free(totalAliveCount);
    for(int i = 0; i < SIZE/mpi_commsize; i++){
        free(myUniverse[i]);
    }
    free(myUniverse);
    free(topGhost);
    free(bottomGhost);
    pthread_barrier_destroy(&barrier);

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/
