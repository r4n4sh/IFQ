/*
 * ACC1.hpp
 *
 *  Created on: May 9, 2017
 *      Author: ranashahout
 */

#ifndef ACC1_HPP_
#define ACC1_HPP_

#include <stdio.h>
#include <unordered_map>
#include <vector>

#include "RSS_CPP.hpp"

using namespace std;

//#define ACC1_DEBUGGING

class ACC1 {
public:
	ACC1(unsigned int windowSize, float gamma, unsigned int m, float epsilon);
    ~ACC1();
    void update(unsigned int item, int wieght);
    double query(unsigned int item);
    double intervalQuery(unsigned int item, int b1, int b2);
#ifdef ACC1_DEBUGGING
    void printHashMaps();
#endif
private:
    int frameItems;
    int blockSize;
    unsigned int windowSize;
    unsigned int blocksNumber;
    unsigned int idx;
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
    /* ACC1 Algorithm */
    unordered_map<unsigned int, unsigned int> idToIDx;
    unordered_map<unsigned int, unsigned int> **overflowedArrL0;
    unordered_map<unsigned int, unsigned int> **overflowedArrL1;
};

#endif /* ACC1_HPP_ */
