/*
 * ACC_K.hpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ranashahout
 */

#ifndef ACC_K_HPP_
#define ACC_K_HPP_


#include <stdio.h>
#include <unordered_map>
#include <vector>

#include "RSS_CPP.hpp"

using namespace std;

//#define ACC_K_DEBUGGING

class ACC_K {
public:
	ACC_K(unsigned int windowSize, float gamma, unsigned int m, float epsilon, unsigned int k);
    ~ACC_K();
    void update(unsigned int item, int wieght);
    double query(unsigned int item);
    double intervalQuery(unsigned int item, int b1, int b2);
#ifdef ACC_K_DEBUGGING
    void printHashMaps();
#endif
private:
    unsigned int k;
    unsigned int step;
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
    /* ACC_K Algorithm */
    unordered_map<unsigned int, unsigned int> idToIDx;
    unordered_map<unsigned int, unsigned int> ***overflowedArrLevels;
};




#endif /* ACC_K_HPP_ */
