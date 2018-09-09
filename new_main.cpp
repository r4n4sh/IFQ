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

//#define HIT_TESTING 1
//#define BASE_WRSS_ALGO 1
#define ACCK_TESTING 1
//#define RAW_TESTING 1

double dblmainmax(double a, double b) {return (a >= b ? a : b);}

int main(int argc, char * argv[]) {
		int counters = 100;
		int n = 100000;
		int k =1;
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
				k = atoi(argv[i]);
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
		int k_algo = 2;


		epsilon = (double)1/(double)counters;
		data = (unsigned long *) malloc(sizeof(unsigned long) * n);

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

			window[i%test_window] = data[i];
		}


		cout << "BEGIN TEST" << endl;
		double estimated, curr_error = 0;
		double exact = 0;

#ifdef HIT_TESTING
        for (i = 0; i < n; i++)  {
            hit->update(data[i], 1);
        }


		interval_size_pkt = ceil(size_precentage * acck->getLastBlock()*block_sz); // 10% of window_size

        for (i = 0; i < n; i++)  {
			double exact = 0;
			int first = rand() % (acck->getLastBlock()*block_sz - interval_size_pkt) + 1;
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
				if (window[k] == data[i])
					exact += 1;
			}

			estimated = acck->intervalFrequencyQuery(data[i],  first,  last);

			cout << "test from first " << first << " to : " << last << " estimated: " << estimated << " exact: " << exact << endl;
			curr_error = exact - estimated;
			curr_error = pow(curr_error, 2);
			emp_error += curr_error;
        }

		emp_error = sqrt((emp_error/n));

		printf( "./acc%d %d pairs emp error: %lf [%d counters %d window_size]\n", k_algo, n, emp_error, counters, window_size);
#endif


		delete[] window;
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
