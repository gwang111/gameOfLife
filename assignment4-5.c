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

#include<clcg4.h>

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

int tick;
int mpi_myrank;
int mpi_commsize;
int aliveCount[NUM_GENERATIONS] = {0};
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

void computeGeneration(int id){
    int updatedTick[SIZE/mpi_commsize/NUM_THREADS][SIZE];
    int g; /* get global index */
    int aliveCounter = 0;
    int aliveNeigh = 0;
    int rowsPerRank = SIZE / mpi_commsize;
    int rowsPerThread = rowsPerRank / NUM_THREADS;
    double rngVal = 0;
    double randVal = 0;
    for(int i = 0; i < rowsPerThread; i++){
        g = (id * rowsPerThread) + (rowsPerRank * mpi_myrank) + i;
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
                aliveNeigh += universe[i][(j - 1) % SIZE];
                aliveNeigh += universe[i][(j + 1) % SIZE];
                if(id == 0) { /* Needs to access top ghost row */
                    aliveNeigh += topGhost[0][(j - 1) % SIZE];
                    aliveNeigh += topGhost[0][(j + 1) % SIZE];
                    aliveNeigh += topGhost[0][j % SIZE];
                    aliveNeigh += universe[i - 1][(j - 1) % SIZE];
                    aliveNeigh += universe[i - 1][(j + 1) % SIZE];
                    aliveNeigh += universe[i - 1][j % SIZE];
                } else if(id == (NUM_THREADS - 1)) { /* Needs to access bottom ghost row */
                    aliveNeigh += bottomGhost[0][(j - 1) % SIZE];
                    aliveNeigh += bottomGhost[0][(j + 1) % SIZE];
                    aliveNeigh += bottomGhost[0][j % SIZE];
                    aliveNeigh += universe[i + 1][(j - 1) % SIZE];
                    aliveNeigh += universe[i + 1][(j + 1) % SIZE];
                    aliveNeigh += universe[i + 1][j % SIZE];
                } else {
                    aliveNeigh += universe[i - 1][(j - 1) % SIZE];
                    aliveNeigh += universe[i - 1][(j + 1) % SIZE];
                    aliveNeigh += universe[i - 1][j % SIZE];
                    aliveNeigh += universe[i + 1][(j - 1) % SIZE];
                    aliveNeigh += universe[i + 1][(j + 1) % SIZE];
                    aliveNeigh += universe[i + 1][j % SIZE];
                }
                if(universe[i][j] == ALIVE) {
                  if(aliveNeigh < 2 || aliveNeigh > 3) {
                    universe[i][j] = DEAD;
                  }
                } else {
                    if(aliveNeigh == 3) {
                      universe[i][j] = ALIVE;
                    }
                }
                aliveCounter += universe[i][j];
            }
        }
    }

    pthread_mutex_lock(&mutex);
    aliveCount[tick] += aliveCounter;
    pthread_mutex_unlock(&mutex);

    pthread_barrier_wait(&barrier);
    for(int i = 0; i < rowsPerThread; i++) {
      for(int j = 0; j < SIZE; j++) {
        universe[id * rowsPerThread + i][j] = updatedTick[i][j];
      }
    }
}

// function for threads
void conways(void * threadID){
    int id = (int) threadID;

    // loop over all generations
    while(tick < NUM_GENERATIONS){
        pthread_barrier_wait(&barrier);
        computeGeneration(id);
    }
}

// function for main thread
void main_conways(){
    MPI_Request sendTop;
    MPI_Request sendBot;
    MPI_Request recvTop;
    MPI_Request recvBot;
    bool mpi_flag;

    for(tick = 0; tick < NUM_GENERATIONS; tick++){

        // send first row up
        MPI_Isend(myUniverse[0], SIZE, MPI_INT, (mpi_myrank-1)%mpi_commsize, 0, MPI_COMM_WORLD, &sendTop);
        // send last row down
        MPI_Isend(myUniverse[SIZE/mpi_commsize-1], SIZE, MPI_INT, (mpi_myrank+1)%mpi_commsize, 0, MPI_COMM_WORLD, &sendBot);
        // receive top ghost
        MPI_Irecv(topGhost, SIZE, MPI_INT, (mpi_myrank-1)%mpi_commsize, 0, &recvTop);
        // receive bottom ghost
        MPI_Irecv(bottomGhost, SIZE, MPI_INT, (mpi_myrank+1)%mpi_commsize, 0, &recvBot);

        MPI_Wait(&sendTop, MPI_STATUS_IGNORE);
        MPI_Test(&sendTop, &mpi_flag, MPI_STATUS_IGNORE);
        if(!mpi_flag){
            fprintf(stderr, "MPI error on send top\n");
            exit(1);
        }

        MPI_Wait(&sendBot, MPI_STATUS_IGNORE);
        MPI_Test(&sendbot, &mpi_flag, MPI_STATUS_IGNORE);
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

        computeGeneration(0);
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

    int *totalAliveCount;
    if(mpi_myrank == 0){
        totalAliveCount = calloc(NUM_GENERATIONS, sizeof(int));
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
        g_end_cycles = GetTimeBase();
        g_time_in_secs = ((double)(g_end_cycles - g_start_cycles)) / g_processor_frequency;

    }


    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/
