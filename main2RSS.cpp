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

//#define TEST_UPDATE 1
//#define TEST_QUERY 1
//#define EMP_ERROR
#define new_emp
#define new_emp_acck
#define TEST_QUERY_INTERVALS
#define TEST_ERROR_MEMORY

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
		unsigned * weights;
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
		int k_algo = 4;

#ifdef EMP_ERROR
	//	window_size = 1 << 20;
#endif

#ifdef TEST_QUERY_INTERVALS
		window_size = 1 << 20; // default window = 2^20
		counters = 256; // default epsilon = 2^(-8)
#endif
		epsilon = (double)1/(double)counters;
		data = (unsigned long *) malloc(sizeof(unsigned long) * n);
		weights = (unsigned *) malloc(sizeof(unsigned) * n);
#if defined(TEST_QUERY) | defined(EMP_ERROR)
		int range = 10;
		int interval_arr_size = ceil(n/range);
		intervals = (unsigned *) malloc(sizeof(unsigned) * interval_arr_size);
#endif
		int block_sz = ceil((window_size * epsilon)/6);
		int test_window =  window_size;
		unsigned long* window = new unsigned long[test_window];

		float size_precentage = percentage/100.0; // percenatge%
		interval_size_pkt = ceil(size_precentage * window_size); // 10% of window_size
		interval_size = interval_size_pkt /block_sz;


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
			fscanf(fp, "%d", weights+i);
#if defined(TEST_QUERY) | defined(EMP_ERROR) | defined(new_emp)
			int interval_idx = i/range;
			intervals[interval_idx] = 1 + (int)rand() % (int)(0.88 * window_size);
#endif

#ifdef new_emp
			window[i%test_window] = data[i];
#endif




/*============== TEMP EMPIRICAL ERROR ================*/
#ifdef EMP_ERROR
			double estimated, curr_error = 0;
			window[i%window_size] = data[i];
#ifdef HIT_TESTING
			hit->update(data[i], 1);
#endif
#ifdef BASE_WRSS_ALGO
			bwrss->update(data[i], 1);
#endif
#ifdef RAW_TESTING
			raw->update(data[i], 1);
#endif
#ifdef ACCK_TESTING
			acck->update(data[i], 1);
#endif
		int debug = 0;

		if (true) {
			double exact = 0;


			if (intervals[interval_idx] + interval_size <= i) {

				for (int k=intervals[interval_idx]; k<intervals[interval_idx] + interval_size; ++k) {
					if (window[k-1] == data[i]) {
						exact += 1;
					}
				}
#ifdef HIT_TESTING
				estimated = hit->intervalQuery(data[i], intervals[interval_idx], intervals[interval_idx] + interval_size);
#endif
#ifdef BASE_WRSS_ALGO
				estimated = bwrss->intervalQuery(data[i], intervals[interval_idx], intervals[interval_idx] + interval_size);
#endif
#ifdef RAW_TESTING
				printf("before query raw\n");
				estimated = raw->intervalQuery(data[i], intervals[interval_idx], intervals[interval_idx] + interval_size);
#endif
#ifdef ACCK_TESTING
				estimated = acck->intervalQuery(data[i], intervals[interval_idx], intervals[interval_idx] + interval_size);
#endif

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
#ifdef RAW_TESTING
		printf( "./raw empirical error: %f\n",emp_error);
#endif
#ifdef ACCK_TESTING
		printf( "./acck empirical error: %f\n",emp_error);
#endif
#else
	}
#endif






/*================== EMPIRICAL ERROR ===================*/




#if defined(new_emp) && defined (new_emp_hit)

		double estimated, curr_error = 0;
		double exact = 0;
		int interval_idx = 0;

        for (i = 0; i < n; i++)  {
            hit->update(data[i], 1);
        }


        for (i = 0; i < n; i++)  {
			double exact = 0;
/*
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

			int i =  rand() % (n - interval_size_pkt);
			int j = i + interval_size_pkt;

			for (int k = i; k<= j; ++k) {
				if (window[k] == data[i])
					exact += 1;
			}
*/
			int first =  rand() % (n - interval_size_pkt);
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

		printf( "./hit empirical error: %f [%d counters %d window_size]\n", emp_error, counters, window_size);

#endif

#if defined(new_emp) && defined (new_emp_acck)
		double estimated, curr_error = 0;
		double exact = 0;
		int interval_idx = 0;

        for (i = 0; i < n; i++)  {
            acck->update(data[i], 1);
        }


        for (i = 0; i < n; i++)  {
			double exact = 0;

/*			int i = rand() % acck->getLastBlock();
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
*/
			int first =  rand() % (n - interval_size_pkt);
			int last = first + interval_size_pkt;

			for (int k = first; k<= last; ++k) {
				if (window[k] == data[i])
					exact += 1;
			}

			estimated = acck->intervalFrequencyQuery(data[i],  first,  last);

			curr_error = exact - estimated;
			curr_error = pow(curr_error, 2);
			emp_error += curr_error;
        }

		emp_error = sqrt((emp_error/n));

		printf( "./acc%d empirical error: %f [%d counters %d window_size]\n",k_algo, emp_error, counters, window_size);

#endif





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

            interval_idx = i/range;
			intervals[interval_idx] = 1 + rand() % hit->getLastBlock();
        }

		begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
            hit->intervalQuery(data[i], intervals[i/range], intervals[i/range] + interval_size);
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
            interval_idx = i/range;
			intervals[interval_idx] = 1 + rand() % acck->getLastBlock();

        }

		begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
            acck->intervalQuery(data[i], intervals[i/range], intervals[i/range] + interval_size);
        }

		endt = clock();
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

        for (i = 0; i < (n/range); i++)  {
			intervals[i] = 1 + rand() % (n - interval_size_pkt);
        }
        
        begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
            hit->intervalFrequencyQuery(data[i], intervals[i/range], intervals[i/range] + interval_size);
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
            acck->intervalFrequencyQuery(data[i], intervals[i/range], intervals[i/range] + interval_size);
        }

		begint = clock();
		ftime(&begintb);
        for (i = 0; i < n; i++)  {
            acck->intervalQuery(data[i], intervals[i/range], intervals[i/range] + interval_size);
        }

		endt = clock();
		time = ((double)(endt-begint))/CLK_PER_SEC;
		memory = maxmemusage();

		printf( "./acck %d pairs took %lfs %dB [%d counters %d window_size]\n", n, time, memory, counters, window_size);
#endif

#endif




/*====================== TEST ERROR MEMORY ===============================*/
#ifdef TEST_ERROR_MEMORY

#ifdef HIT_TESTING
        for (i = 0; i < n; i++)  {
            hit->update(data[i], 1);
        }


        for (i = 0; i < n; i++)  {
			double exact = 0;
			int first =  rand() % (n - interval_size_pkt);
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

		printf( "./hhh2RSS %d pairs emp error: %lf [%d counters %d window_size]\n", n, emp_error, counters, window_size);
#endif

#ifdef ACCK_TESTING
        for (i = 0; i < n; i++)  {
            acck->update(data[i], 1);
        }


        for (i = 0; i < n; i++)  {
			double exact = 0;
			int first =  rand() % (n - interval_size_pkt);
			int last = first + interval_size_pkt;

			for (int k = first; k<= last; ++k) {
				if (window[k] == data[i])
					exact += 1;
			}

			estimated = acck->intervalFrequencyQuery(data[i],  first,  last);

			curr_error = exact - estimated;
			curr_error = pow(curr_error, 2);
			emp_error += curr_error;
        }

		emp_error = sqrt((emp_error/n));

		printf( "./acc%d %d pairs emp error: %lf [%d counters %d window_size]\n", k_algo, n, emp_error, counters, window_size);
#endif
#endif







#if defined(TEST_QUERY) | defined (EMP_ERROR)
		free(intervals);
#endif
#ifndef EMP_ERROR
		free(data);
#endif
#ifdef RAW_TESTING
		free(raw);
#endif
#ifdef BASE_WRSS_ALGO
		free(bwrss);
#endif
#ifdef HIT_TESTING
		free(hit);
#endif
		return 0;
}

