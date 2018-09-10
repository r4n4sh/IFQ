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
    int getLastBlock() {
        return lastBlock;
    }
    double intervalFrequencyQuery(unsigned int item, int i, int j);

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
    int lastBlock;

    int m; //maximal value of an element in the stream
    //double alpha = 0.2;
    float gamma;
    double computeOverflowCount(unsigned int item);
    unsigned int withinFrameFrequency(unsigned int required_block, int itemIdx);
    void populateIncTable(unsigned int itemIdx);
    void endBlock(int blockNumber);
    double winQuery(unsigned int item, int w);

    /* ACC_K Algorithm */
    unordered_map<unsigned int, unsigned int> idToIDx;
    unordered_map<unsigned int, unsigned int>*** overflowedArrLevels;
    unordered_map<unsigned int, unsigned int> **incTable;
    unordered_map<unsigned int, unsigned int> **ghost_tables;

};




#endif /* ACC_K_HPP_ */
