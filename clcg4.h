/*------------------------------------------------------------------------*/
/* clcg4.h :     Definitions of Interface                                 */
/*                                                                        */
/* A Random Number Generator Based on the Combination of Four LCGs        */
/*                                                                        */
/*------------------------------------------------------------------------*/
#ifndef H_CLCG4_H
#define H_CLCG4_H
/* Maximum number of generators, plus one.  This value can be increased 
   as needed */
#define Maxgen (32768)

/* A generator number, in the range 0 to Maxgen */
typedef unsigned short int Gen;
typedef enum {InitialSeed, LastSeed, NewSeed} SeedType;


/* Initializes the random number package.  The values of V and W will be 
   2^v and 2^w respectively.  It is recommended that v>=30, w>=40 and 
   v+w<=100.  The initial seed is set by default to 
   (11111111, 22222222, 33333333, 44444444). */
void Init( long v, long w);

/* Initializes the package as in Init, with the default values V=2^31 and 
   W=w^41. */
void InitDefault( void);

/* Set the initial seed I_0 of generator number 0 (g=0) to the integer values 
   s[0], ..., s[3].  Those values must satisfy: 1<=s[0]<=2147483646,
   1<=s[1]<=2147483542, 1<=s[2]<=2147483422, 1<=s[3]<=2147483322.  
   The initial seeds of all other generators are recompiled accordingly, 
   so they are spaced VW values apart, and all generators are
   reinitialized to their initial seeds. */
void SetInitialSeed( long s[4]);

/* Reinitialize the generator g.  According to the value of Where, that 
   generator's state C_g will be reset to the initial seed I_g 
   (InitialSeed), or to the last seed L_g (LastSeed), which is at the 
   beginning of the current segment, or to a new seed (NewSeed), W values 
   ahead of the last seed in the generator's sequence.
   The last seed L_g is also reset to the same values as C_g. */
void InitGenerator( Gen g, SeedType Where);
 
/* Set initial seed I_g of generator g to s[1], ..., s[4].  Those values 
   must satisfy: 1<=s[0]<=2147483646, 1<=s[1]<=2147483542, 
   1<=s[2]<=2147483422, 1<=s[3]<=2147483322.
   The current state C_g and last seed L_g are also put to that seed. */
void SetSeed( Gen g, long s[4]);

/* Returns the current state C_g of generator g in s[0], ..., s[3]. */
void GetState( Gen g, long s[4]);

/* Writes the current state C_g of generator g. */
void WriteState( Gen g);

/* Returns a "uniform" random number over [0,1], using generator g.  
   The current state C_g is changed, but not I_g and L_g. */
double GenVal( Gen g);
#endif
