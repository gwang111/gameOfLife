/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Ryan, John, Gary             **(*****************************************/
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

#define SIZE 1024
#define NUM_THREADS 8
#define NUM_GENERATIONS 256
#define THRESHOLD .25
#define PARALLEL_IO 1
#define HEATMAP 0

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
#ifdef BGQ
    unsigned long long g_start_cycles=0;
    unsigned long long g_end_cycles=0;
#else
    double g_start_cycles;
    double g_end_cycles;
#endif

int mpi_myrank;
int mpi_commsize;
int aliveCount[NUM_GENERATIONS];
//Universe
int ** myUniverse;
//Ghost Rows
int * topGhost;
int * bottomGhost;
pthread_barrier_t barrier;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
//Updates to unverser chunks that will be copied over to the universe
int ** updatedTick;


// You define these


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these

// compute a % b
// works with negative a
int mod(int a, int b){
    while(a < 0)
        a += b;
    return a % b;
}

void computeGeneration(int id, int tick){
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
                aliveNeigh += myUniverse[row][mod(j-1, SIZE)];
                aliveNeigh += myUniverse[row][(j + 1) % SIZE];
                // top neighbors
                if(id == 0 && i == 0){ /* Needs to access top ghost row */
                    aliveNeigh += topGhost[mod(j-1, SIZE)];
                    aliveNeigh += topGhost[(j + 1) % SIZE];
                    aliveNeigh += topGhost[j % SIZE];
                }
                else{
                    aliveNeigh += myUniverse[row - 1][mod(j-1, SIZE)];
                    aliveNeigh += myUniverse[row - 1][(j + 1) % SIZE];
                    aliveNeigh += myUniverse[row - 1][j % SIZE];
                }
                // bottom neighbors
                if(id == NUM_THREADS-1 && i == rowsPerThread-1){ /* Needs to access bottom ghost row */
                    aliveNeigh += bottomGhost[mod(j-1, SIZE)];
                    aliveNeigh += bottomGhost[(j + 1) % SIZE];
                    aliveNeigh += bottomGhost[j % SIZE];
                }
                else{
                    aliveNeigh += myUniverse[row + 1][mod(j-1, SIZE)];
                    aliveNeigh += myUniverse[row + 1][(j + 1) % SIZE];
                    aliveNeigh += myUniverse[row + 1][j % SIZE];
                }

                if(myUniverse[i][j] == ALIVE) {
                  if(aliveNeigh < 2 || aliveNeigh > 3) {
                    updatedTick[row][j] = DEAD;
                  }
                  else{
                    updatedTick[row][j] = ALIVE;
                  }
                  // else still alive
                }
                // if dead
                else {
                    // bring back to life if exactly 3 neighbors
                    if(aliveNeigh == 3) {
                        updatedTick[row][j] = ALIVE;
                    }
                    else{
                        updatedTick[row][j] = DEAD;
                    }

                }

                aliveCounter += updatedTick[row][j];
            }
        }
    }
    //Lock and unlock when updating the alive count
    pthread_mutex_lock(&mutex);
    aliveCount[tick] += aliveCounter;
    pthread_mutex_unlock(&mutex);

    //Copy updatedtick over to be used in next tick
    for(int i = 0; i < rowsPerThread; i++) {
      for(int j = 0; j < SIZE; j++) {
        myUniverse[id * rowsPerThread + i][j] = updatedTick[id * rowsPerThread + i][j];
      }
    }
}

// function for threads
void *conways(void * threadID){
    int id = *(int*) threadID;

    // loop over all generations
    for(int tick = 0; tick < NUM_GENERATIONS; tick++){
       // MPI_Barrier(MPI_COMM_WORLD);
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

        // send first row up
        MPI_Isend(myUniverse[0], SIZE, MPI_INT, mod(mpi_myrank-1, mpi_commsize), tick, MPI_COMM_WORLD, &sendTop);
        // receive bottom ghost
        MPI_Irecv(bottomGhost, SIZE, MPI_INT, (mpi_myrank+1)%mpi_commsize, tick, MPI_COMM_WORLD, &recvBot);
        // send last row down
        MPI_Isend(myUniverse[SIZE/mpi_commsize-1], SIZE, MPI_INT, (mpi_myrank+1)%mpi_commsize, tick, MPI_COMM_WORLD, &sendBot);
        // receive top ghost
        MPI_Irecv(topGhost, SIZE, MPI_INT, mod(mpi_myrank-1, mpi_commsize), tick, MPI_COMM_WORLD, &recvTop);

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

void inputOutput(){
    MPI_Barrier(MPI_COMM_WORLD);
    // begin timer
    if(mpi_myrank == 0)
        g_start_cycles = GetTimeBase();

    //Basically just copying the whole universe to an outfile using MPI write at because it will synchronize between ranks
    //Using an offset
    MPI_File cFile;
    MPI_File_open(MPI_COMM_WORLD, "output.txt", MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &cFile);
    for(int line = 0; line < SIZE/mpi_commsize; line++){
        MPI_Offset offset = (mpi_myrank * (SIZE/mpi_commsize) + line) * (SIZE+1);
        char row[SIZE+1];
        for(int i = 0; i < SIZE; i ++)
            row[i] = myUniverse[line][i] + '0';
        row[SIZE] = '\n';
        MPI_File_write_at(cFile, offset, row, SIZE+1, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    MPI_File_close(&cFile);
    //End the timer after written
    //This should be IO time
    if(mpi_myrank == 0){
        g_end_cycles = GetTimeBase();
        #ifdef BGQ
            g_time_in_secs = ((double)(g_end_cycles - g_start_cycles)) / g_processor_frequency;
        #else
            g_time_in_secs = (g_end_cycles - g_start_cycles);
        #endif
        printf("I/O: %f\n", g_time_in_secs);
    }
}

// reduce the universe and create the heatmap
void heatmap(){
    MPI_Barrier(MPI_COMM_WORLD);
    //Condense by 32 for columns
    // reduce the columns by 32 first because not dependant on rank
    int * myReduced = calloc((SIZE/32)*(SIZE/mpi_commsize), sizeof(int));
    int index = 0;
    for(int row = 0; row < SIZE/mpi_commsize; row++){
        int sum = 0;
        int block = 31;
        for(int i = 0; i < SIZE; i++){
            if(i == block){
                myReduced[index] = sum;
                sum = 0;
                block += 32;
                index++;
            }
            sum += myUniverse[row][i];
        }
    }

    //Allocate the condensed universe before adding to it
    int * reduced = NULL;
    if(mpi_myrank == 0){
        reduced = calloc(SIZE*(SIZE/32), sizeof(int));
    }

    //Condense by 32 for rows
    //Gather and reduce the rows now because MPI_Gather will synchronize between ranks allow access of all rows
    MPI_Gather(myReduced, (SIZE/32)*(SIZE/mpi_commsize), MPI_INT, reduced, (SIZE/32)*(SIZE/mpi_commsize), MPI_INT, 0, MPI_COMM_WORLD);
    if(mpi_myrank == 0){
        FILE * f = fopen("heatmap.txt", "w");
        //This is where the row reduction and output file writing happens
        for(int row = 0; row < SIZE/32; row++){
            for(int col = 0; col < SIZE/32; col++){
                int sum = 0;
                for(int i = 0; i < 32; i++){
                    sum += reduced[row*SIZE + col + i*(SIZE/32)];
                }
                fprintf(f, "%d", sum);
                if(col != SIZE/32-1)
                    fputc(' ', f);
            }
            if(row != SIZE/32-1)
                fputc('\n', f);
        }

        fclose(f);
        free(reduced);
    }
    free(myReduced);
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

    //Variables to be used later
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
    updatedTick = calloc(SIZE/mpi_commsize, sizeof(int*));
    for(int i = 0; i < SIZE/mpi_commsize; i++){
        myUniverse[i] = calloc(SIZE, sizeof(int));
        updatedTick[i] = calloc(SIZE, sizeof(int));
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
    int ids[NUM_THREADS];
    for(int i = 0; i < NUM_THREADS-1; i++){
        ids[i] = i+1;
        pthread_create(&my_threads[i], &attr, conways, (void *) &ids[i]);
    }
    main_conways();

    // Wait for all the threads
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

    // output results and end timer
    //This time is not I/O time
    //Time differs from inputOutput() I/O time
    if(mpi_myrank == 0){
        for(int i = 0; i < NUM_GENERATIONS; i++){
            printf("Generation %d: %d alive\n", i, totalAliveCount[i]);
        }

        g_end_cycles = GetTimeBase();
        #ifdef BGQ
            g_time_in_secs = ((double)(g_end_cycles - g_start_cycles)) / g_processor_frequency;
        #else
            g_time_in_secs = (g_end_cycles - g_start_cycles);
        #endif
        printf("TIME: %f\n", g_time_in_secs);
    }

    /* PARALLEL I/O */
    if(PARALLEL_IO){
        inputOutput();
    }
    /* HEATMAP */
    if(HEATMAP){
        heatmap();
    }

    // cleanup
    //Free what needs to be freed and exit what needs to be exited
    if(mpi_myrank == 0)
        free(totalAliveCount);
    free(topGhost);
    free(bottomGhost);
    pthread_barrier_destroy(&barrier);
    for(int i = 0; i < SIZE/mpi_commsize; i++){
        free(myUniverse[i]);
        free(updatedTick[i]);
    }
    free(myUniverse);
    free(updatedTick);

    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}
