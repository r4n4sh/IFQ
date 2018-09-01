/*
 * ACC.hpp
 *
 *  Created on: Apr 12, 2017
 *      Author: ranashahout
 */

#ifndef ACC_HPP_
#define ACC_HPP_

#include <stdio.h>
#include <unordered_map>
#include <vector>

#include "RSS_CPP.hpp"

using namespace std;

class ACC {
public:
	ACC(unsigned int windowSize, float gamma, unsigned int m, float epsilon);
    ~ACC();
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
    int lastBlock;
    int m; //maximal value of an element in the stream
    //double alpha = 0.2;
    float gamma;
    double computeOverflowCount(unsigned int item);
    void populateIncTable(unsigned int itemIdx);
    int getLastBlock() {
        return lastBlock;
    }

    /* Naive Algorithm */
    unordered_map<int, int> idToIDx;
    unordered_map<unsigned int, unsigned int> **overflowedArr;
    unordered_map<unsigned int, unsigned int> **incTable;

};


#endif /* ACC_HPP_ */
