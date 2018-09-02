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
#include <cstring>
#include <cmath>


#ifdef PARALLEL
#include <omp.h>
#endif

#include "alloc.hpp"

#ifndef CLK_PER_SEC
#ifdef CLOCKS_PER_SEC
#define CLK_PER_SEC CLOCKS_PER_SEC
#endif
#endif


//#define TEST_UPDATE 1
//#define TEST_QUERY 1
#define TRY_2
//#define EMP_ERROR
#define new_emp
#define new_emp_acck

//#define HIT_TESTING 1
//#define BASE_WRSS_ALGO 1
//#define ACC_TESTING 1
#define ACCK_TESTING 1
//#define ACC1_TESTING 1
//#define RAW_TESTING 1

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
		int k =1;
		double time;
		int memory;
		int i;
		int w, x, y, z;
		clock_t begint, endt;
		struct timeb begintb, endtb;
		unsigned * weights;
		unsigned * intervals;
		unsigned long * data;
		FILE * fp = NULL;
		int M = 1;
		float gamma = 4;
		double epsilon = 0.01;
		unsigned int interval_size;
		int window_size = 1600;
		int interval_1;
		int interval_2;
		double emp_error = 0;

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
			else if (strcmp(argv[i], "-k") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout  << "Missing k" << std::endl;
					return -1;
				}
				k = atoi(argv[i]);
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
					std::cout << "Missing interval1." << std::endl;
					return -1;
				}

				interval_1 = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-j") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout << "Missing interval2." << std::endl;
					return -1;
				}

				interval_2 = atoi(argv[i]);
			}
			else
			{
				cout << "Unknown parameter" << argv[i] << endl;
				return -1;
			}
		}

		gamma = 0;
		M = 1;
		int k_algo = 4;
		if(n / counters >= threshold) {
			printf("Unacceptable parameters: eps*n >= theshold\n");
			return 0;
		}

#ifdef EMP_ERROR
	//	window_size = 1 << 20;
#endif
//		interval_size = ceil(window_size / 100); 1% of window_size
		interval_size = ceil(window_size / 10); // 10% of window_size

		epsilon = (double)1/(double)counters;
		data = (unsigned long *) malloc(sizeof(unsigned long) * n);
		weights = (unsigned *) malloc(sizeof(unsigned) * n);
#if defined(TEST_QUERY) | defined(EMP_ERROR)
		int range = 10;
		int interval_arr_size = ceil(n/range);
		intervals = (unsigned *) malloc(sizeof(unsigned) * interval_arr_size);
#endif

#ifdef new_emp
		int range = 10;

		int interval_arr_size = ceil(n/range);
		intervals = (unsigned *) malloc(sizeof(unsigned) * interval_arr_size);
#endif


#ifdef HIT_TESTING
		HIT *hit = new HIT(window_size, gamma, M, epsilon);
#endif
#ifdef BASE_WRSS_ALGO
		BaseWRSS *bwrss = new BaseWRSS(window_size, gamma, M, epsilon);
#endif
#ifdef ACC_TESTING
		ACC *acc = new ACC(window_size, gamma, M, epsilon);
#endif
#ifdef RAW_TESTING
		RAW *raw = new RAW(window_size, gamma, M, epsilon);
#endif
#ifdef ACC1_TESTING
		ACC1 *acc1 = new ACC1(window_size, gamma, M, epsilon);
#endif
#ifdef ACCK_TESTING
		ACC_K *acck = new ACC_K(window_size, gamma, M, epsilon, k_algo);
#endif

		int block_sz = ceil((window_size * epsilon)/6);
		int test_window =  window_size*block_sz;
		unsigned long* window = new unsigned long[test_window];

		for (i = 0; i < n; i++) {
			fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
			data[i] = (unsigned long)256*((unsigned long)256*((unsigned long)256*w + x) + y) + z;
			fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
			fscanf(fp, "%d", weights+i);
#if defined(TEST_QUERY) | defined(EMP_ERROR) | defined(new_emp)
			int interval_idx = i/range;
			intervals[interval_idx] = 1 + (int)rand() % (int)(0.88 * window_size);
			//printf("interval_idx: %d\n", intervals[interval_idx]);
#endif

#ifdef new_emp
			window[i%test_window] = data[i];
#endif

#ifdef EMP_ERROR
			double estimated, curr_error = 0;
			window[i%window_size] = data[i];
			//printf("window[%d] = %d \n", i%window_size, window[i%window_size] );
#ifdef HIT_TESTING
			hit->update(data[i], 1);
#endif
#ifdef BASE_WRSS_ALGO
			bwrss->update(data[i], 1);
#endif
#ifdef ACC_TESTING
			acc->update(data[i], 1);
#endif
#ifdef RAW_TESTING
			raw->update(data[i], 1);
#endif
#ifdef ACC1_TESTING
			acc1->update(data[i], 1);
#endif
#ifdef ACCK_TESTING
			acck->update(data[i], 1);
#endif
		//if (i > 1500000 && !(i %100000)) {
		int debug = 0;

		//if (i > (n/10) && !debug) {
		if (true) {

			// count it in window
			double exact = 0;
			//printf("Query interval: first idx: %d, second idx: %d\n", intervals[interval_idx], intervals[interval_idx] + interval_size, window_size);


			if (intervals[interval_idx] + interval_size <= i) {
				//printf("[%d] count from window, first idex: %d, secoond idx: %d, window_size: %d\n", i, intervals[interval_idx], intervals[interval_idx] + interval_size, window_size);

				for (int k=intervals[interval_idx]; k<intervals[interval_idx] + interval_size; ++k) {
					if (window[k-1] == data[i]) {
						exact += 1;
					}
				}
				//printf("EXACT: %f\n", exact);
#ifdef HIT_TESTING
				estimated = hit->intervalQuery(data[i], intervals[interval_idx], intervals[interval_idx] + interval_size);
#endif
#ifdef BASE_WRSS_ALGO
				estimated = bwrss->intervalQuery(data[i], intervals[interval_idx], intervals[interval_idx] + interval_size);
#endif
#ifdef ACC_TESTING
				estimated = acc->intervalQuery(data[i], intervals[interval_idx], intervals[interval_idx] + interval_size);
#endif
#ifdef RAW_TESTING
				printf("before query raw\n");
				estimated = raw->intervalQuery(data[i], intervals[interval_idx], intervals[interval_idx] + interval_size);
#endif
#ifdef ACC1_TESTING
				estimated = acc1->intervalQuery(data[i], intervals[interval_idx], intervals[interval_idx] + interval_size);
#endif
#ifdef ACCK_TESTING
				estimated = acck->intervalQuery(data[i], intervals[interval_idx], intervals[interval_idx] + interval_size);
#endif

				//printf("[%d] estimated: %f\n",i, estimated);
				//printf("[%d] exact: %f\n",i, exact);
				curr_error = exact - estimated;
				printf("[%d] estimated: %f exact: %f curr_error: %f, interval: [%d, %d]\n", i, estimated, exact, curr_error, intervals[interval_idx],intervals[interval_idx] + interval_size);
				curr_error = pow(curr_error, 2);
				emp_error += curr_error;
				}
			}
		}
		emp_error = sqrt((emp_error/n));

#ifdef HIT_TESTING
		printf( "[%d] ./hit empirical error: %f\n", i , emp_error);
#endif
#ifdef BASE_WRSS_ALGO
		printf( "./baseWRSS empirical error: %f\n",emp_error);
#endif
#ifdef ACC_TESTING
		printf( "./acc empirical error: %f\n",emp_error);
#endif
#ifdef RAW_TESTING
		printf( "./raw empirical error: %f\n",emp_error);
#endif
#ifdef ACC1_TESTING
		printf( "./acc1 empirical error: %f\n",emp_error);
#endif
#ifdef ACCK_TESTING
		printf( "./acck empirical error: %f\n",emp_error);
#endif
#else
	}
#endif











#if defined(new_emp) && defined (new_emp_hit)

		double estimated, curr_error = 0;
		double exact = 0;
		int interval_idx = 0;

        for (i = 0; i < n; i++)  {
            hit->update(data[i] & masks[0], 1);
        }


        for (i = 0; i < n; i++)  {
			double exact = 0;

			int i = rand() % hit->getLastBlock();
			int interval_size = rand() % (hit->getLastBlock() - i);
			int j = i + interval_size;
			if (j > hit->getLastBlock())
				j = hit->getLastBlock();
	        int b1 = hit->getLastBlock() - j;
	        int b2 = b1 + interval_size;

			for (int k = b1*block_sz; k<= b2*block_sz; ++k) {
				if (window[k] == data[i])
					exact += 1;
			}
        

			estimated = hit->intervalQuery(data[i], i, j);

			curr_error = exact - estimated;
		//	printf("[%d] estimated: %f exact: %f curr_error: %f, interval: [%d, %d]\n", i, estimated, exact, curr_error, b1, b2);

			curr_error = pow(curr_error, 2);
			emp_error += curr_error;
        }

		emp_error = sqrt((emp_error/n));

		printf( "./hit empirical error: %f%dB [%d counters %d window_size]\n", emp_error, counters, window_size);

#endif

#if defined(new_emp) && defined (new_emp_acck)
		double estimated, curr_error = 0;
		double exact = 0;
		int interval_idx = 0;

        for (i = 0; i < n; i++)  {
            acck->update(data[i] & masks[0], 1);
        }


        for (i = 0; i < n; i++)  {
			double exact = 0;

			int i = rand() % acck->getLastBlock();
			int interval_size = rand() % (acck->getLastBlock() - i);
			int j = i + interval_size;
			if (j > acck->getLastBlock())
				j = acck->getLastBlock();
	        int b1 = acck->getLastBlock() - j;
	        int b2 = b1 + interval_size;

			for (int k = b1*block_sz + 1; k<= b2*block_sz; ++k) {
				if (window[k] == data[i])
					exact += 1;
			}
        

			estimated = acck->intervalQuery(data[i], i, j);

			curr_error = exact - estimated;
			//printf("[%d] estimated: %f exact: %f curr_error: %f, interval: [%d, %d]\n", i, estimated, exact, curr_error, b1, b2);

			curr_error = pow(curr_error, 2);
			emp_error += curr_error;
        }

		emp_error = sqrt((emp_error/n));

		printf( "./acc%d empirical error: %f%dB [%d counters %d window_size]\n",k_algo, emp_error, counters, window_size);

#endif

















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

#ifdef ACCK_TESTING
		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
        	acck->update(data[i] & masks[0], 1);
        }

        endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acck %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
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
#ifdef RAW_TESTING
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
#ifdef TRY_2
            hit->intervalQuery(data[i] & masks[0], intervals[i/range], intervals[i/range] + interval_size);
#else
            hit->intervalQuery(data[i] & masks[0], intervals[i], intervals[i] + interval_size);
#endif
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
#ifdef TRY_2
            acc->intervalQuery(data[i] & masks[0], intervals[i/range], intervals[i/range] + interval_size);
#else
            acc->intervalQuery(data[i] & masks[0], intervals[i], intervals[i] + interval_size);
#endif
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
#ifdef TRY_2
            raw->intervalQuery(data[i] & masks[0], intervals[i/range], intervals[i/range] + interval_size);
#else
            raw->intervalQuery(data[i] & masks[0], intervals[i], intervals[i] + interval_size);
#endif
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./raw %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef ACCK_TESTING
        for (i = 0; i < n; i++)  {
            acck->update(data[i] & masks[0], 1);
        }

		begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
#ifdef TRY_2
            acck->intervalQuery(data[i] & masks[0], intervals[i/range], intervals[i/range] + interval_size);
#else
            acck->intervalQuery(data[i] & masks[0], intervals[i], intervals[i] + interval_size);
#endif
        }

		endt = clock();
		ftimeintervals[interval_idx]	time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acck %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef ACC1_TESTING
        for (i = 0; i < n; i++)  {
            acc1->update(data[i] & masks[0], 1);
        }

		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
#ifdef TRY_2
            acc1->intervalQuery(data[i] & masks[0], intervals[i/range], intervals[i/range] + interval_size);
#else
            acc1->intervalQuery(data[i] & masks[0], intervals[i], intervals[i] + interval_size);
#endif
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acc1 %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif
#endif

#if defined(TEST_QUERY) | defined (EMP_ERROR)
		free(intervals);
#endif
#ifndef EMP_ERROR
		free(data);
#endif
#ifdef ACC1_TESTING
		free(acc1);
#endif
#ifdef RAW_TESTING
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

