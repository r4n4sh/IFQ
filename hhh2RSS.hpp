
#ifndef HHH2RSS_H
#define HHH2RSS_H

#if NUM_MASKS==1089
#define NUM_COUNTERS 1089 //number of masks
#define MAX_DEPTH 33 //depth of masks lattice
//#define MAX_DESCENDANTS 512 //maximum number of direct descendents of a given ip pair
#else
#define NUM_COUNTERS 25 //number of masks
#define MAX_DEPTH 5 //depth of masks lattice
#define MAX_DESCENDANTS 512 //maximum number of direct descendents of a given ip pair
#endif

#include "prng.hpp"
#include "RSS_CPP.hpp"
#include "HIT.hpp"
#include "ACC.hpp"
#include "BaseWRSS.hpp"
#include "RAW.hpp"
//make sure we are passing the right #def around
#ifndef DIMENSION2
#error Invalid dimension
#endif

//still define things as in lossycount.h
#ifdef DIMENSION2
#define LCLitem_t uint64__t
#else
#define LCLitem_t uint32_t
#endif

//The masks associated with the counters
//Note that we must ensure that they are in increasing order of generality
extern LCLitem_t masks[NUM_COUNTERS];

//initialise
void init(double epsilon, float gamma, int M);
//deinitialise
void deinit();
void query(LCLitem_t item);

#ifndef PARALLEL
//update an input
void update(LCLitem_t item, int count);
#else
void update(LCLitem_t * item, int count);
#endif

#endif


