/*
Implementation of the 2-dimensional Hierarchical Heavy Hitters algorithm (HHH) for IP addresses
-Thomas Steinke (tsteinke@seas.harvard.edu) 2010-10-06
*/

#include <stdlib.h>
#include <stdio.h>

//for debugging purposes only
#include <assert.h> 
//#include <set> //debug
//#include <utility> //debug

#ifdef UNITARY
#include "ulossycount.h"
#define COUNTERSIZE k
#define COUNTERS items
#define COUNT parentg->count
#define LCL_type LCU_type
#define LCL_Init LCU_Init
#define LCL_Destroy LCU_Destroy
#define LCL_Update(x,y,z) LCU_Update(x,y)
#define LCL_PointEstUpp LCU_PointEstUpp
#define LCL_PointEstLow LCU_PointEstLow
#else
#include "lossycount.hpp"
#define COUNTERSIZE size
#define COUNTERS counters
#define COUNT count
#endif
//#include "hashtable.h"

#include "hhh2RSS.hpp"
#include "alloc.hpp"

#ifndef DIMENSION2
#error This is the two-dimensional implementation 
#endif


#define D(x...) fprintf(stderr, x) //debug
#define P(x...) fprintf(stderr, x)
#define PIP(item) fprintf(stderr, "%3d.%3d.%3d.%3d", (int)(255&((item) >> 24)), (int)(255&((item) >> 16)), (int)(255&((item) >> 8)), (int)(255&((item) >> 0)))


//#define NAIVE_ALGO 1
//#define NUM_COUNTERS 16 //number of masks
//#define MAX_DESCENDANTS 512 //maximum number of direct descendents of a given ip pair

int min(int a, int b) {return (a <= b ? a : b);}
int max(int a, int b) {return (a >= b ? a : b);}

//The counters
#ifdef RSS
RSS_CPP * counters[NUM_COUNTERS];
#endif

#ifdef WRSS
BaseWRSS *wrss;
#endif

#ifdef NAIVE_ALGO
NaiveWRSS * native_wrss;
#endif
//The masks associated with the counters
//Note that we must ensure that they are in increasing order of generality
/*LCLitem_t masks[NUM_COUNTERS] = {
	//255.255.255.255
	0xFFFFFFFFFFFFFFFFull, //255.255.255.255
	0xFFFFFF00FFFFFFFFull, //255.255.255.0
	0xFFFF0000FFFFFFFFull, //255.255.0.0
	0xFF000000FFFFFFFFull, //255.0.0.0

	//255.255.255.0
	0xFFFFFFFFFFFFFF00ull, //255.255.255.255
	0xFFFFFF00FFFFFF00ull, //255.255.255.0
	0xFFFF0000FFFFFF00ull, //255.255.0.0
	0xFF000000FFFFFF00ull, //255.0.0.0

	//255.255.0.0
	0xFFFFFFFFFFFF0000ull, //255.255.255.255
	0xFFFFFF00FFFF0000ull, //255.255.255.0
	0xFFFF0000FFFF0000ull, //255.255.0.0
	0xFF000000FFFF0000ull, //255.0.0.0

	//255.0.0.0
	0xFFFFFFFFFF000000ull, //255.255.255.255
	0xFFFFFF00FF000000ull, //255.255.255.0
	0xFFFF0000FF000000ull, //255.255.0.0
	0xFF000000FF000000ull  //255.0.0.0
};*/

int leveleps[NUM_COUNTERS] = {
64, 56, 48, 40, 32,
56, 48, 40, 32, 24,
48, 40, 32, 24, 16,
40, 32, 24, 16,  8,
32, 24, 16,  8,  0
};

double dblmax(double a, double b) {return (a >= b ? a : b);}

double twototheminus(int k) {
	double ans = 1;
	while (k > 0) {ans /= 2; k--;}
	return ans;
}

//initialise
void init(double epsilon, float gamma, int M) {
	int window_size = 1600;

#ifdef RSS
	//int i;
	//for (i = 0; i < NUM_COUNTERS; i++)
		//counters[i] = new RSS_CPP(epsilon, M, gamma);
	counters[0] = new RSS_CPP(epsilon, M, gamma);
#endif
#ifdef WRSS
	printf("Creating WRSS\n");
	wrss = new WRSS(window_size, gamma, M, epsilon);
#endif
#ifdef NAIVE_ALGO
	window_size = 1600;
	epsilon = 0.01;
	native_wrss = new NaiveWRSS(window_size, gamma, M, epsilon);
#endif
}

//deinitialise
void deinit() {
}

void query(LCLitem_t item) {
#ifdef RSS
    printf("RSS query = %d \n", counters[0]->query(item));
#endif

#ifdef WRSS
    printf("WRSS query = %f \n", wrss->query(item));
#endif
#ifdef NAIVE_ALGO
    printf("Naive Algorithm query = %f\n", native_wrss->query(item));
#endif
}

#ifndef PARALLEL
//update an input
void update(LCLitem_t item, int count) {
#ifdef RSS
	//int i;
	//for (i = 0; i < NUM_COUNTERS; i++) {
		//counters[i]->update(item & masks[i], count);
		//LCL_Update(counters[i], item & masks[i], count);
		//P("update [%2d] ", i); PIP((item & masks[i]) >> 32); P(" "); PIP(item & masks[i]); P("\n");
	//}

	counters[0]->update(item & masks[0], count);
#endif
#ifdef WRSS
	wrss->update(item & masks[0], count);
#endif
#ifdef NAIVE_ALGO
	native_wrss->update(item & masks[0], count);
#endif
}
#else
#endif
