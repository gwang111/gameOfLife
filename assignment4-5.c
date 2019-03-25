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

int mpi_myrank;
int mpi_commsize;
int aliveCount[NUM_GENERATIONS];
int ** myUniverse;
int * topGhost;
int * bottomGhost;
pthread_barrier_t barrier;
long ** seeds;


// You define these


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these

void computeGeneration(int id){

    for(int i = 0; i < SIZE/mpi_commsize/NUM_THREADS; i++){
        // determine the RNG value for this row
        SetInitialSeed();

        for(int j = 0; j < SIZE; j++){

            // if RNG < threshold, random alive/dead

            // else follow rules
        }
    }
}

// function for threads
void conways(void * threadID){
    int id = (int) threadID;

    // loop over all generations
    for(int tick = 0; tick < NUM_GENERATIONS; tick++){
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

    for(int tick = 0; tick < NUM_GENERATIONS; tick++){
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
    seeds = calloc(SIZE/mpi_commsize, sizeof(long*));
    for(int i = 0; i < SIZE/mpi_commsize; i++){
        seeds[i] = calloc(4, sizeof(long));
    }

// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    //printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	//   mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    
    MPI_Barrier( MPI_COMM_WORLD );
    
// Insert your code
    if(mpi_myrank == 0)
        g_start_cycles = GetTimeBase();

    pthread_t my_threads[NUM_THREADS];
    pthread_attr_t attr;
    pthread_attr_init(&attr);

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

    // create threads
    for(int i = 1; i < NUM_THREADS; i++){
        int id = i;
        pthread_create(&my_threads[i], &attr, conways, (void *) &id);
    }

    main_conways();
    // END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/
