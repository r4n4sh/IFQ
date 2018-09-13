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

/*********** TESTS ***************/

//#define TEST_UPDATE 
//#define TEST_QUERY
#define TEST_QUERY_INTERVALS
//#define TEST_ERROR_MEMORY

#define HIT_TESTING 1
//#define BASE_WRSS_ALGO 1
//#define ACCK_TESTING 1
//#define RAW_TESTING 1

double dblmainmax(double a, double b) {return (a >= b ? a : b);}

int main(int argc, char * argv[]) {
		int counters = 100;
		int n = 100000;
		int k_algo = 2;
		double time;
		int memory;
		int i;
		int w, x, y, z;
		clock_t begint, endt;
		struct timeb begintb, endtb;
		unsigned *weights;
		unsigned * intervals;
		unsigned long * data;
		FILE * fp = NULL;
		int M = 1;
		float gamma = 4;
		double epsilon = 0.01;
		unsigned int interval_size;
		unsigned int interval_size_pkt;
		int percentage = 1;
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
				k_algo = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-t") == 0)
			{
				i++;
				if (i >= argc)
				{
					std::cout << "Missing interval size percentage (of window size)" << std::endl;
					return -1;
				}
				percentage = atoi(argv[i]);
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


#ifdef TEST_QUERY_INTERVALS
		window_size = 1 << 20; // default window = 2^20
		counters = 256; // default epsilon = 2^(-8)
#endif
		epsilon = (double)1/(double)counters;
		data = (unsigned long *) malloc(sizeof(unsigned long) * n);

#if defined(TEST_QUERY)  | defined(TEST_QUERY_INTERVALS)
		int range = 10000;
		//int interval_arr_size = ceil(n/range);
		//intervals = (unsigned *) malloc(sizeof(unsigned) * interval_arr_size);
#endif
		int block_sz = ceil((window_size * epsilon)/6);
		int test_window =  window_size;
		unsigned long* window = new unsigned long[test_window];

		float size_precentage = percentage/100.0; // percenatge%
		interval_size_pkt = ceil(size_precentage * window_size); // 10% of window_size
		interval_size = interval_size_pkt /block_sz;

#ifdef HIT_TESTING
		HIT *hit = new HIT(window_size, gamma, M, epsilon);
#endif
#ifdef BASE_WRSS_ALGO
		BaseWRSS *bwrss = new BaseWRSS(window_size, gamma, M, epsilon);
#endif
#ifdef RAW_TESTING
		RAW *raw = new RAW(window_size, gamma, M, epsilon);
#endif
#ifdef ACCK_TESTING
		ACC_K *acck = new ACC_K(window_size, gamma, M, epsilon, k_algo);
#endif


		for (i = 0; i < n; i++) {
			fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
			data[i] = (unsigned long)256*((unsigned long)256*((unsigned long)256*w + x) + y) + z;
			fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
			fscanf(fp, "%d", weights);
#if defined(TEST_QUERY)
			int interval_idx = i/range;
			intervals[interval_idx] = 1 + (int)rand() % (int)(0.88 * window_size);
#endif

#if defined(TEST_ERROR_MEMORY) 
			window[i%test_window] = data[i];
#endif

		}




/*========================== UPDATE ==========================*/


#ifdef TEST_UPDATE
#ifdef HIT_TESTING
		begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
            hit->update(data[i], 1);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./hhh2RSS %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif


#ifdef ACCK_TESTING
		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
        	acck->update(data[i], 1);
        }

        endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acc%d %d pairs took %lfs %dB [%d counters %d window_size]\n", k_algo, n, time, memory, counters, window_size);
#endif

#ifdef BASE_WRSS_ALGO
		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
        	bwrss->update(data[i], 1);
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
        	raw->update(data[i], 1);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./raw %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#endif






/*====================== QUERY ===============================*/
#ifdef TEST_QUERY
		/* Test Query times */

#ifdef HIT_TESTING
        for (i = 0; i < n; i++)  {
            hit->update(data[i], 1);
        }

        for (i = 0; i < (n/range); i++)  {
			intervals[i] = 1 + rand() % (n - interval_size_pkt);
        }


		begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
            hit->intervalFrequencyQuery(data[i], intervals[i], intervals[i] + interval_size);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./hhh2RSS %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef BASE_WRSS_ALGO
        for (i = 0; i < n; i++)  {
        	bwrss->update(data[i], 1);
        }

		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
            bwrss->query(data[i]);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./baseWRSS %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef RAW_TESTING
        for (i = 0; i < n; i++)  {
        	raw->update(data[i], 1);
        }

		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
            raw->intervalQuery(data[i], intervals[i/range], intervals[i/range] + interval_size);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./raw %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef ACCK_TESTING

        for (i = 0; i < n; i++)  {
            acck->update(data[i], 1);
        }

        for (i = 0; i < (n/range); i++)  {
			intervals[i] = 1 + rand() % (n - interval_size_pkt);
        }


		begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
            acck->intervalFrequencyQuery(data[i], intervals[i], intervals[i] + interval_size);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acc%d %d pairs took %lfs %dB [%d counters %d window_size]\n", k_algo, n, time, memory, counters, window_size);
#endif

#endif


/*====================== TEST VARY INTERVAL SIZES ===============================*/
#ifdef TEST_QUERY_INTERVALS

#ifdef HIT_TESTING
        for (i = 0; i < n; i++)  {
            hit->update(data[i], 1);
        }


		interval_size_pkt = ceil(size_precentage * hit->getLastBlock()*block_sz); // 10% of window_size

/*
    	assert(hit->getLastBlock() > interval_size);

        for (i = 0; i < (n/range); i++)  {
			intervals[i] = 1 + rand() % (hit->getLastBlock() - interval_size);
        }


		begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
            hit->intervalQuery(data[i], intervals[i/range], intervals[i/range] + interval_size);
        }
*/
    	int first = rand() % (hit->getLastBlock()*block_sz - interval_size_pkt) + 1;
		int last = first + interval_size_pkt;

/*
        for (i = 0; i < (n/range); i++)  {
			intervals[i] =  rand() % (window_size - interval_size_pkt) + 1;
        }*/
        
        begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
         //   hit->intervalFrequencyQuery(data[i], intervals[i], intervals[i] + interval_size_pkt);
              hit->intervalFrequencyQuery(data[i], first, last);

        }


		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./hit %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef RAW_TESTING
        for (i = 0; i < n; i++)  {
        	raw->update(data[i], 1);
        }

        for (i = 0; i < (n/range); i++)  {
			intervals[i] = 1 + rand() % (n - interval_size_pkt);
        }


		begint = clock();
		ftime(&begintb);

        for (i = 0; i < n; i++)  {
            raw->intervalQuery(data[i], intervals[i], intervals[i] + interval_size_pkt);
        }

		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./raw %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#ifdef ACCK_TESTING
        for (i = 0; i < n; i++)  {
            acck->update(data[i], 1);
        }

		interval_size_pkt = ceil(size_precentage * acck->getLastBlock()*block_sz); // 10% of window_size

    	//int first = rand() % (acck->getLastBlock()*block_sz - interval_size_pkt) + 1;
		//int last = first + interval_size_pkt;
    	int first = block_sz
		int last = interval_size_pkt;

		cout << "first: " << first << "last: " << last << end;
/*
        for (i = 0; i < (n/range); i++)  {
			intervals[i] =  rand() % (window_size - interval_size_pkt) + 1;
        }*/
        
        begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
//            acck->intervalFrequencyQuery(data[i], intervals[i], intervals[i] + interval_size_pkt);
            acck->intervalFrequencyQuery(data[i], first,  last);

        }

		endt = clock();
		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acc%d %d pairs took %lfs %dB [%d counters %d window_size]\n", k_algo, n, time, memory, counters, window_size);
#endif

#endif




/*====================== TEST ERROR MEMORY ===============================*/
#ifdef TEST_ERROR_MEMORY

		double estimated, curr_error = 0;
		double exact = 0;

#ifdef HIT_TESTING
        for (i = 0; i < n; i++)  {
            hit->update(data[i], 1);
        }


		interval_size_pkt = ceil(size_precentage * hit->getLastBlock()*block_sz); // 10% of window_size

        for (i = 0; i < n; i++)  {
			double exact = 0;
			int first = rand() % (hit->getLastBlock()*block_sz - interval_size_pkt) + 1;
			int last = first + interval_size_pkt;

			for (int k = first; k<= last; ++k) {
				if (window[k] == data[i])
					exact += 1;
			}

			estimated = hit->intervalFrequencyQuery(data[i],  first,  last);

			curr_error = exact - estimated;
			curr_error = pow(curr_error, 2);
			emp_error += curr_error;
        }

		emp_error = sqrt((emp_error/n));

		printf( "./hit %d pairs emp error: %lf [%d counters %d window_size]\n", n, emp_error, counters, window_size);
#endif

#ifdef ACCK_TESTING
        for (i = 0; i < n; i++)  {
            acck->update(data[i], 1);
        }


		interval_size_pkt = ceil(size_precentage * acck->getLastBlock()*block_sz); // 10% of window_size

        for (i = 0; i < n; i++)  {
			double exact = 0;
			int first = rand() % (acck->getLastBlock()*block_sz - interval_size_pkt) + 1;
			int last = first + interval_size_pkt;
			

			for (int k = first; k<= last; ++k) {
				aif (window[k] == data[i])
					exact += 1;
			}

			estimated = acck->intervalFrequencyQuery(data[i],  first,  last);

			//cout << "test from first " << first << " to : " << last << " estimated: " << estimated << " exact: " << exact << endl;
			curr_error = exact - estimated;
			curr_error = pow(curr_error, 2);
			emp_error += curr_error;
        }

		emp_error = sqrt((emp_error/n));

		printf( "./acc%d %d pairs emp error: %lf [%d counters %d window_size]\n", k_algo, n, emp_error, counters, window_size);
#endif
#endif


		delete[] window;

#if defined(TEST_QUERY) | defined(TEST_QUERY_INTERVALS)
	//	free(intervals);
#endif
		free(data);

#ifdef RAW_TESTING
		delete raw;
#endif
#ifdef ACCK_TESTING
		delete acck;
#endif
#ifdef BASE_WRSS_ALGO
		delete bwrss;
#endif
#ifdef HIT_TESTING
		delete hit;
#endif
		return 0;
}
