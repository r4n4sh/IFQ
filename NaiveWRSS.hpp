/*
 * NaiveWRSS.hpp
 *
 *  Created on: Apr 12, 2017
 *      Author: ranashahout
 */

#ifndef NAIVEWRSS_HPP_
#define NAIVEWRSS_HPP_

#include <stdio.h>
#include <unordered_map>
#include <vector>

#include "RSS_CPP.hpp"

using namespace std;

class NaiveWRSS {
public:
	NaiveWRSS(unsigned int windowSize, float gamma, unsigned int m, float epsilon);
    ~NaiveWRSS();
    void update(unsigned int item, int wieght);
    double query(unsigned int item);
    double intervalQuery(unsigned int item, int b1, int b2);

private:
    int frameItems;
    int blockSize;
    unsigned int windowSize;
    unsigned int blocksNumber;
    vector<int> *index;
    int tail;
    int overflowsNumber;
    int indexTail;
    int head;
    int indexHead;
    unordered_map<int, int> *totalOverflows; //B
    unsigned int *overflowsElements;// b
    RSS_CPP *rss;
    int indexSize;
    int maxOverflows;
    float epsilon;
    unsigned int threshold;
    int m; //maximal value of an element in the stream
    //double alpha = 0.2;
    float gamma;
    double computeOverflowCount(unsigned int item);

    /* Naive Algorithm */
    unordered_map<int, int> idToIDx;
    int** overflowedArr;
};


#endif /* NAIVEWRSS_HPP_ */
