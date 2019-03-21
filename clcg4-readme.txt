Introduction
============

This is a portable package for uniform random number generation based on a
backbone generator with period length near 2^121, which is a combination
of four linear congruential generators.  The package provides for multiple
(virtual) generators evolving in parallel.  Each generator also has many
disjoint subsequences, and software tools are provided to reset the state
of any generator to the beginning of its first, previous, or current
subsequence.  Such facilities are helpful to maintain synchronization for
implementing variance reduction methods in simulation.

For those who are interested in the details or the theories behind the
implementation, you can refer to the following papers by Prof. L'Ecuyer:

1. A Random Number Generator Based on the Combination of Four LCGs.
   (This can be downloaded from the "papers" page of Prof. L'Ecuyer's 
    web page: http://www.iro.umontreal.ca/~lecuyer/myftp/papers/clcg4.ps).
2. Implementing a Random Number Package with Splitting Facilities
   (ACM Transactions on Math. Software, Vol. 17, No. 1, March 91,
    p. 98-111)
3. Efficient and portable combined random number generators.  
   (Communications of the ACM, 31(6):742-749 and 774, 1988.)

and the papers referred to therein.


What is really in this package
==============================
There are 101 (virtual) generators, (numbered from 0 to 100) with seeds
2^(v+w) values apart.  It is based on a combined LCG with period length ~
2^121.  Each generator corresponds to a subsequence of 2^(v+w) values,
which is further split into V=2^v segments of length W=2^w.  (Both v and w
can be changed by the user.  Their default values are v=31 and w=41.)

The initial seed of the first generator (generator[0]) can be set to any
'admissible' value (a vector of 4 +ve integers) and the initial seeds of
the other generators are automatically recalculated by the package so that
they remain VW apart.  (Further details are available below.)  This is
implemented with efficient jump-ahead tools.

By a simple procedure call, any generator can jump ahead to either the
beginning of its next segment, or the beginning of its current segment, or
its first (initial) segment.

This package has been subjected to a battery of statistical tests and it
passed all of them without any problem.  So, empirically, it behaves pretty
well.  Though empirical statistical testing never proves that a generator
is safe, at least it improves one's confidence.

The implementation is written in ANSI-C and is supposed to run on 32-bit
machines.  


How can I use the RNG?
======================
Step 1.
-------
You should download two files, namely, clcg4.h and clcg4.c.  In order to
use the package you should include clcg4.h in your program.

Step 2.
-------
At the very beginning of your program, you must first initialize the RNG
package by either of the following calls:
  InitDefault();
or
  Init(v,w);

InitDefault() is basically the same as Init(v,w); it simply uses the
default values of v=31 and w=41.  It is recommended that v>=30, w>=40 and
v+w<=100.  

Step 3. (optional)
------------------
The initial seed is set by default to
  {11111111,22222222,33333333,44444444}.
You can change the seeds by calling 
  SetInitialSeed( s);
AFTER calling Init(v,w) or InitDefault() where s is an 4-element array of
"long int".

Notice that the seeds MUST satisfy the following constraints:
  1<=s[0]<=2147483646
  1<=s[1]<=2147483542
  1<=s[2]<=2147483422
  1<=s[3]<=2147483322.

Step 4.
-------
Now, you can freely get a "uniform" random number over [0,1], using
generator g by calling e.g.
  X=GenVal(g);


"Bells and Whistles"
====================
1. You can obtain information about current state of generator g by
   calling
     GetState( g, s);
   This will return the current state of generator g in s[0],...s[3].
   Again, s is an array of long int.

2. You can set the initial seed of a particular generator (say generator
   g) to s[0],...s[3] by
     SetSeed( g, s);
   Typically, you use it combined with GetState( g, s) in order to reuse a
   particular segment of generator g.

3. Finally, you can reinitialize generator g by calling
     InitGenerator( g, where);
   According to the value of 'where', generator g's state will be reset to
   the initial seed (InitialSeed), or to the last seed (LastSeed), which
   is at the beginning of the current segment, or to a new seed 
   (NewSeed), 2^w values ahead of the last seed in the generator's
   sequence.


Can I test it?
==============
You can download the five files: 
  clcg4.h, clcg4.c, main.c, makefile and main.out.
Put them in an empty subdirectory of your elaine home directory and type
  make.
(The makefile uses the GNU C compiler (on Unix, Linux or DJGPP on the PCs).)

You can then see the output of the program by typing 
  ./main
The file main.c simply serves the purpose of illustrating how to use the
RNG package.

The output of the (sample) program should be identical to the file
main.out.  For users using other ANSI-C compilers on a 32-bit machine, 
you can simply create a "project" file, "insert" the three files: 
  clcg4.c, clcg4.h and main.c 
in it and "run".

