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


#define TEST_UPDATE 1
#define TEST_QUERY 1


#define HIT_TESTING 1
#define BASE_WRSS_ALGO 1
#define ACC_TESTING 1
#define ACC1_TESTING 1
#define RAW_TESING 1
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
		int counters = 100;
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
		int M = 1;
		float gamma = 4;
		double epsilon = 0.01;
		int window_size = 1600;
		int interval = 0;
#ifdef TEST_QUERY
		int interval_1 = 1 + (int)rand() % (int)(0.99 * window_size);
		int interval_2 = (window_size/100) + interval_1;
#endif
		for (int i = 1; i < argc; ++i)
		{
			if (strcmp(argv[i], "-np") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout << "Missing number of packets." << std::endl;
					return -1;
				}
				n = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-c") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout  << "Missing epsilon" << std::endl;
					return -1;
				}
				counters = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-t") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout << "Missing threshold." << std::endl;
					return -1;
				}
				threshold = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-f") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout << "Missing trace file." << std::endl;
					return -1;
				}
				fp = fopen(argv[4], "w");
			}

			else if (strcmp(argv[i], "-gamma") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout << "Missing gamma." << std::endl;
					return -1;
				}
				gamma = atof(argv[i]);
			}
			else if (strcmp(argv[i], "-M") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout << "Missing M." << std::endl;
					return -1;
				}
				M = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-w") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout << "Missing window size." << std::endl;
					return -1;
				}

				window_size = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-i") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout << "Missing interval." << std::endl;
					return -1;
				}

				interval = atoi(argv[i]);
			}
			else
			{
				cout << "Unknown parameter" << argv[i] << endl;
				return -1;
			}
		}


		if(n / counters >= threshold) {
			printf("Unacceptable parameters: eps*n >= theshold\n");
			return 0;
		}
		epsilon = (double)1/(double)counters;
		data = (unsigned long *) malloc(sizeof(unsigned long) * n);
		weights = (unsigned *) malloc(sizeof(unsigned) * n);


#ifdef HIT_TESTING
		HIT *hit = new HIT(window_size, gamma, M, epsilon);
#endif
#ifdef BASE_WRSS_ALGO
		BaseWRSS *bwrss = new BaseWRSS(window_size, gamma, M, epsilon);
#endif
#ifdef ACC_TESTING
		ACC *acc = new ACC(window_size, gamma, M, epsilon);
#endif
#ifdef RAW_TESING
		RAW *raw = new RAW(window_size, gamma, M, epsilon);
#endif
#ifdef ACC1_TESTING
		ACC1 *acc1 = new ACC1(window_size, gamma, M, epsilon);
#endif
		for (i = 0; i < n; i++) {
			fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
			data[i] = (unsigned long)256*((unsigned long)256*((unsigned long)256*w + x) + y) + z;
			fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
			fscanf(fp, "%d", weights+i);
		}
#ifdef TEST_UPDATE
#ifdef HIT_TESTING
		begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
            hit->update(data[i] & masks[0], 1);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./hhh2RSS %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef ACC_TESTING
		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
            acc->update(data[i] & masks[0], 1);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acc %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef BASE_WRSS_ALGO
		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
        	bwrss->update(data[i] & masks[0], 1);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./baseWRSS %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif
#ifdef RAW_TESING
		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
        	raw->update(data[i] & masks[0], 1);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./raw %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef ACC1_TESTING
		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
            acc1->update(data[i] & masks[0], 1);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acc1 %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif
#endif

#ifdef TEST_QUERY
		/* Test Query times */
#ifdef HIT_TESTING
        for (i = 0; i < n; i++)  {
            hit->update(data[i] & masks[0], 1);
        }

		begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
            hit->intervalQuery(data[i] & masks[0], interval_1, interval_2);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./hhh2RSS %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef ACC_TESTING
        for (i = 0; i < n; i++)  {
            acc->update(data[i] & masks[0], 1);
        }

		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
            acc->intervalQuery(data[i] & masks[0], interval_1, interval_2);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acc %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef BASE_WRSS_ALGO
        for (i = 0; i < n; i++)  {
        	bwrss->update(data[i] & masks[0], 1);
        }

		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
            bwrss->query(data[i] & masks[0]);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./baseWRSS %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef RAW_TESTING
        for (i = 0; i < n; i++)  {
        	raw->update(data[i] & masks[0], 1);
        }

		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
            raw->intervalQuery(data[i] & masks[0], interval_1, interval_2);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./raw %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif
#ifdef ACC1_TESTING
        for (i = 0; i < n; i++)  {
            acc1->update(data[i] & masks[0], 1);
        }

		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
            acc1->intervalQuery(data[i] & masks[0], interval_1, interval_2);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acc1 %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif
#endif
		free(data);
#ifdef ACC1_TESTING
		free(acc1);
#endif
#ifdef RAW_TESING
		free(raw);
#endif
#ifdef ACC_TESTING
		free(acc);
#endif
#ifdef BASE_WRSS_ALGO
		free(bwrss);
#endif
#ifdef HIT_TESTING
		free(hit);
#endif
		return 0;
}

