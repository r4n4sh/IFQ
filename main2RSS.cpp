/*
The main function for the two-dimensional HHH program
- Thomas Steinke (tsteinke@seas.harvard.edu) 2010-11
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "hhh2RSS.hpp"
#include <sys/timeb.h>

#ifdef PARALLEL
#include <omp.h>
#endif

#include "alloc.hpp"

#ifndef CLK_PER_SEC
#ifdef CLOCKS_PER_SEC
#define CLK_PER_SEC CLOCKS_PER_SEC
#endif
#endif

#define INTERVAL_ALGO 1
//#define BASE_WRSS_ALGO 1
//#define INTERVAL_NAIVE_ALGO 1
//the masks
LCLitem_t masks[NUM_COUNTERS] = {
	//255.255.255.255
	//0-4
	0xFFFFFFFFFFFFFFFFull, //255.255.255.255
	0xFFFFFF00FFFFFFFFull, //255.255.255.0
	0xFFFF0000FFFFFFFFull, //255.255.0.0
	0xFF000000FFFFFFFFull, //255.0.0.0
	0x00000000FFFFFFFFull, //0.0.0.0

	//255.255.255.0
	//5-9
	0xFFFFFFFFFFFFFF00ull, //255.255.255.255
	0xFFFFFF00FFFFFF00ull, //255.255.255.0
	0xFFFF0000FFFFFF00ull, //255.255.0.0
	0xFF000000FFFFFF00ull, //255.0.0.0
	0x00000000FFFFFF00ull, //0.0.0.0

	//255.255.0.0
	//10-14
	0xFFFFFFFFFFFF0000ull, //255.255.255.255
	0xFFFFFF00FFFF0000ull, //255.255.255.0
	0xFFFF0000FFFF0000ull, //255.255.0.0
	0xFF000000FFFF0000ull, //255.0.0.0
	0x00000000FFFF0000ull, //0.0.0.0

	//255.0.0.0
	//15-19
	0xFFFFFFFFFF000000ull, //255.255.255.255
	0xFFFFFF00FF000000ull, //255.255.255.0
	0xFFFF0000FF000000ull, //255.255.0.0
	0xFF000000FF000000ull, //255.0.0.0
	0x00000000FF000000ull, //0.0.0.0

	//0.0.0.0
	//20-24
	0xFFFFFFFF00000000ull, //255.255.255.255
	0xFFFFFF0000000000ull, //255.255.255.0
	0xFFFF000000000000ull, //255.255.0.0
	0xFF00000000000000ull, //255.0.0.0
	0x0000000000000000ull  //0.0.0.0
};


/*int _main() {
		int m; //number of heavy hitters in output
		int counters, threshold, n;
		scanf("%d%d%d", &counters, &threshold, &n);
		HeavyHitter * ans;
		int i;
		unsigned long long wa, xa, ya, za;
		unsigned long long wb, xb, yb, zb;
		unsigned long long w, x, y, z;
		unsigned long long a, b, ip, _ip;

		init((double)1/(double)counters);

		for (i = 0; i < n; i++) {
			scanf("%llu%llu%llu%llu", &w, &x, &y, &z);
			ip = (unsigned long long)256*((unsigned long long)256*((unsigned long long)256*w + x) + y) + z;
			scanf("%llu%llu%llu%llu", &w, &x, &y, &z);
			_ip = (unsigned long long)256*((unsigned long long)256*((unsigned long long)256*w + x) + y) + z;
			update(ip << 32 | _ip , 1);
		}

		ans = output2(threshold, &m);

		deinit();

		for (i = 0; i < m; i++) {
			//output ans[i]

			//break up the ip
			a = ans[i].item >> 32;
			za = a % 256; a /= 256;
			ya = a % 256; a /= 256;
			xa = a % 256; a /= 256;
			wa = a % 256; a /= 256;

			//break up the mask
			b = masks[ans[i].mask] >> 32;
			zb = b % 256; b /= 256;
			yb = b % 256; b /= 256;
			xb = b % 256; b /= 256;
			wb = b % 256; b /= 256;

			//output ip&mask
			if (wb != 0) printf("%llu.", wa); else printf("*.");
			if (xb != 0) printf("%llu.", xa); else printf("*.");
			if (yb != 0) printf("%llu.", ya); else printf("*.");
			if (zb != 0) printf("%llu ", za); else printf("* ");

			//break up the ip
			a = ans[i].item;
			za = a % 256; a /= 256;
			ya = a % 256; a /= 256;
			xa = a % 256; a /= 256;
			wa = a % 256; a /= 256;

			//break up the mask
			b = masks[ans[i].mask];
			zb = b % 256; b /= 256;
			yb = b % 256; b /= 256;
			xb = b % 256; b /= 256;
			wb = b % 256; b /= 256;

			//output ip&mask
			if (wb != 0) printf("%llu.", wa); else printf("*.");
			if (xb != 0) printf("%llu.", xa); else printf("*.");
			if (yb != 0) printf("%llu.", ya); else printf("*.");
			if (zb != 0) printf("%llu", za); else printf("*");

			//output counts
			//printf("	[%d, %d] %d %d\n",  ans[i].lower, ans[i].upper, ans[i].s, ans[i].t); //debug
			printf("	[%d, %d]\n",  ans[i].lower, ans[i].upper); //debug
		}

		free(ans);

	return 0;
}*/

double dblmainmax(double a, double b) {return (a >= b ? a : b);}

int main(int argc, char * argv[]) {
		double counters = 100;
		int threshold = 1000;
		int n = 100000;
		double time;
		int memory;
		int i;
		int w, x, y, z;
		clock_t begint, endt;
		struct timeb begintb, endtb;
		unsigned * weights;
		unsigned long * data;
		FILE * fp = NULL;
		int M = 65535;
		float gamma = 4;

		if (argc > 1) n = atoi(argv[1]);
		if (argc > 2) counters = atoi(argv[2]);
		if (argc > 3) threshold = atoi(argv[3]);
		n = 32;
		counters = 100;
		threshold = 19531;
		fp = fopen("/Users/ranashahout/Documents/workspace/eclipse_oxygen/wrss_cross/Chicago15Weighted.txt.1m.txt", "r");

		if (argc > 5) M = atoi(argv[5]);
		if (argc > 6) gamma = atof(argv[6]);

		if(n/counters >= threshold) {
			printf("Unacceptable parameters: eps*n >= theshold\n");
			return 0;
		}

		data = (unsigned long *) malloc(sizeof(unsigned long) * n);
		weights = (unsigned *) malloc(sizeof(unsigned) * n);

		int window_size = 32;
		double epsilon = 0.5;
#ifdef INTERVAL_ALGO
		WRSS *wrss = new WRSS(window_size, gamma, M, epsilon);
#endif
#ifdef BASE_WRSS_ALGO
		BaseWRSS *bwrss = new BaseWRSS(window_size, gamma, M, epsilon);
#endif
#ifdef INTERVAL_NAIVE_ALGO
		NaiveWRSS *nwrss = new NaiveWRSS(window_size, gamma, M, epsilon);
#endif
		for (i = 0; i < n; i++) {
			fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
			data[i] = (unsigned long)256*((unsigned long)256*((unsigned long)256*w + x) + y) + z;
			fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
			fscanf(fp, "%d", weights+i);
		}
		begint = clock();
		ftime(&begintb);
		#ifndef PARALLEL

        for (i = 0; i < n; i++)  {
#ifdef INTERVAL_ALGO
            wrss->update(data[i] & masks[0], weights[i]);
#endif
#ifdef BASE_WRSS_ALGO
            bwrss->update(data[i] & masks[0], weights[i]);
#endif
#ifdef INTERVAL_NAIVE_ALGO
            nwrss->update(data[i] & masks[0], weights[i]);
#endif
        }

		#else
		#endif
		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();
		
		printf( "%d update pairs took %lfs %dB [%f counters]\n", n, time, memory, counters);

		begint = clock();
		ftime(&begintb);
		#ifndef PARALLEL

        for (i = 0; i < n; i++)  {
#ifdef INTERVAL_ALGO
            wrss->intervalQuery(data[i] & masks[0], 1, window_size);
#endif
#ifdef BASE_WRSS_ALGO
            bwrss->query(data[i] & masks[0]);
#endif
#ifdef INTERVAL_NAIVE_ALGO
            nwrss->intervalQuery(data[i] & masks[0], 1, window_size);
#endif
        }

		#else
		#endif
		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "%d query pairs took %lfs %dB [%f counters]\n", n, time, memory, counters);
#ifdef INTERVAL_ALGO
#ifdef TESTING
		double result = wrss->testIntervalQuery(3 , 8);
		cout << "test interval query: " << result << endl;
#else
		double result = wrss->intervalQuery(data[0], 3, 8);
		cout << "test interval query: " << result << endl;
#endif
#endif
		free(data);
#ifdef INTERVAL_NAIVE_ALGO
		free(nwrss);
#endif
#ifdef BASE_WRSS_ALGO
		free(bwrss);
#endif
#ifdef INTERVAL_ALGO
		free(wrss);
#endif
		return 0;
}

