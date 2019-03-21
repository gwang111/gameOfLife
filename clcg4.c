/*---------------------------------------------------------------------*/
/* clcg4.c   Implementation module                                     */
/*---------------------------------------------------------------------*/
#include <stdio.h>
#include "clcg4.h"

/***********************************************************************/
/* Private part.                                                       */
/***********************************************************************/
#define H (32768)               /* = 2^15 : use in MultModM.           */

static long aw[4], avw[4],      /*   a[j]^{2^w} et a[j]^{2^{v+w}}.     */
            a[4]={ 45991, 207707, 138556, 49689},
            m[4]={ 2147483647, 2147483543, 2147483423, 2147483323};

static long Ig[4][Maxgen+1], Lg[4][Maxgen+1], Cg[4][Maxgen+1];
                     /* Initial seed, previous seed, and current seed. */

static short i, j;

static long MultModM( long s, long t, long M) {
   /* Returns (s*t) MOD M.  Assumes that -M < s < M and -M < t < M.    */
   /* See L'Ecuyer and Cote (1991).                                    */
  long R, S0, S1, q, qh, rh, k;

  if( s<0) s+=M;
  if( t<0) t+=M;
  if( s<H) { S0=s; R=0;}
  else {
    S1=s/H; S0=s-H*S1;
    qh=M/H; rh=M-H*qh;
    if( S1>=H) {
      S1-=H; k=t/qh; R=H*(t-k*qh)-k*rh;
      while( R<0) R+=M;
    }
    else R=0;
    if( S1!=0) {
      q=M/S1; k=t/q; R-=k*(M-S1*q);
      if( R>0) R-=M;
      R += S1*(t-k*q);
      while( R<0) R+=M;
    }
    k=R/qh; R=H*(R-k*qh)-k*rh;
    while( R<0) R+=M;
  }
  if( S0!=0) {
    q=M/S0; k=t/q; R-=k*(M-S0*q);
    if( R>0) R-=M;
    R+=(S0*(t-k*q));
    while( R<0) R+=M;
  }
  return R;
}


/*---------------------------------------------------------------------*/
/* Public part.                                                        */
/*---------------------------------------------------------------------*/
void SetSeed( Gen g, long s[4]) {
  if( g>Maxgen) printf( "ERROR: SetSeed with g > Maxgen\n");
  for( j=0; j<4; j++) Ig[j][g]=s[j];
  InitGenerator( g, InitialSeed);
}


void WriteState( Gen g) {
  printf ("\n State of generator g = %u :", g);
  for( j=0; j<4; j++) printf ("\n   Cg[%u] = %lu", j, Cg[j][g]);
  printf ("\n");
}


void GetState( Gen g, long s[4]) {
  for( j=0; j<4; j++) s[j]=Cg[j][g];
}


void InitGenerator( Gen g, SeedType where) {
  if( g>Maxgen) printf( "ERROR: InitGenerator with g > Maxgen\n");
  for( j=0; j<4; j++) {
    switch (where) {
      case InitialSeed :
        Lg[j][g]=Ig[j][g]; break;
      case NewSeed :
        Lg[j][g]=MultModM( aw[j], Lg[j][g], m[j]); break;
      case LastSeed :
        break;
    }
    Cg[j][g]=Lg[j][g];
  }
}


void SetInitialSeed( long s[4]) {
  Gen g;

  for( j=0; j<4; j++) Ig[j][0]=s[j];
  InitGenerator( 0, InitialSeed);
  for( g=1; g<=Maxgen; g++) {
    for( j=0; j<4; j++) Ig[j][g]=MultModM( avw[j], Ig[j][g-1], m[j]);
    InitGenerator( g, InitialSeed);
  }
}


void Init( long v, long w) {
  long sd[4]={11111111, 22222222, 33333333, 44444444};

  for( j=0; j<4; j++) {
    for( aw[j]=a[j], i=1; i<=w; i++) aw[j]=MultModM( aw[j], aw[j], m[j]);
    for( avw[j]=aw[j], i=1; i<=v; i++) avw[j]=MultModM( avw[j], avw[j], m[j]);
  }
  SetInitialSeed (sd);
}


double GenVal( Gen g) {
  long k,s;
  double u=0.0;

  if( g>Maxgen) printf( "ERROR: Genval with g > Maxgen\n");

  s=Cg[0][g]; k=s/46693;
  s=45991*(s-k*46693)-k*25884;
  if( s<0) s+=2147483647; 
  Cg[0][g]=s;
  u+=(4.65661287524579692e-10*s);
 
  s=Cg[1][g]; k=s/10339;
  s=207707*(s-k*10339)-k*870;
  if( s<0) s+=2147483543;  
  Cg[1][g]=s;
  u-=(4.65661310075985993e-10*s);
  if( u<0) u+=1.0;

  s=Cg[2][g]; k=s/15499;
  s=138556*(s-k*15499)-k*3979;
  if( s<0.0) s+=2147483423;  
  Cg[2][g]=s;
  u+=(4.65661336096842131e-10*s);
  if( u>=1.0) u-=1.0;

  s=Cg[3][g]; k=s/43218;
  s=49689*(s-k*43218)-k*24121;
  if( s<0) s+=2147483323;  
  Cg[3][g]=s;
  u-=(4.65661357780891134e-10*s);
  if( u<0) u+=1.0;

  return (u);
}

void InitDefault( void) {
  Init( 31, 41);
}
