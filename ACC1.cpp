/*
 * ACC1.cpp
 *
 *  Created on: May 9, 2017
 *      Author: ranashahout
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "ACC1.hpp"

using namespace std;

#ifdef ACC1_DEBUGGING
int testArr[8] = {1,2,3,0,0,1,1,0};
unsigned int firstOverflow;
unsigned int thirdOverflow;
#endif

static inline int nextMulOf4(int number) {
	if (number == 1) return 4;
	return 4 * (int)ceil(log(number)/log(4));
}

static inline int prevL1Block(int number) {
	return (int)floor(double(number/4));
}

ACC1::ACC1(unsigned int windowSize, float gamma, unsigned int m, float epsilon)
{
	idx = 0;
    frameItems = 0;
    overflowsNumber = 0;
    tail = 0;
    indexTail = 0;
    this->blockSize = ceil((windowSize * epsilon)/4); // W/k; k= 4/epsilon
    this->windowSize = windowSize;
    this->blocksNumber = windowSize / blockSize;
    this->m = m;
    this->epsilon = epsilon;
    this->gamma = gamma;
    maxOverflows = min(blocksNumber * 2, windowSize);
    indexSize = maxOverflows + blocksNumber;
    head = maxOverflows - 1;
    indexHead = blocksNumber - 1;
    index = new vector<int> (indexSize); // 0 means end of block
    overflowsElements = new unsigned int[maxOverflows];
    rss = new RSS_CPP(epsilon, m, gamma);//y
    threshold = ceil(windowSize*m*epsilon / 4.f); // W*M/k
    totalOverflows = new unordered_map<int, int> (maxOverflows);//B TODO: allocate static?

    overflowedArrL0 = new unordered_map<unsigned int, unsigned int>*[blocksNumber];
    for(int i =0 ; i < blocksNumber ; i++)
    	overflowedArrL0[i] =  new unordered_map<unsigned int, unsigned int> (maxOverflows);

    int l1Size = ceil(blocksNumber/4); //TODO: define the magic number 4
    overflowedArrL1 = new unordered_map<unsigned int, unsigned int>*[l1Size];
    for(int i = 0 ; i < l1Size ; i++)
    	overflowedArrL1[i] = new unordered_map<unsigned int, unsigned int> (maxOverflows);
#ifdef ACC1_DEBUGGING
    cout << "L0 blocks: " << blocksNumber << " L1 blocks: " << l1Size << endl;
#endif
}

ACC1::~ACC1()
{
    int l1Size = ceil(blocksNumber/4); //TODO: define the magic number 4
    for (int i = 0; i < l1Size; i++)
    	delete(overflowedArrL1[i]);
    delete[] overflowedArrL1;

    for (int i = 0; i < blocksNumber; i++)
    	delete(overflowedArrL0[i]);
    delete[] overflowedArrL0;

    delete(totalOverflows);
    delete(rss);
    delete [] overflowsElements;
    delete(index);
}

void ACC1::update(unsigned int item, int wieght)
{
	/* New Block */
    if (((++frameItems) % blockSize) == 0) {
        indexTail = (indexTail + 1) % indexSize;
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 0;
    }

    // Remove oldest element in oldest block

    try {
        if (index->at(indexTail)) {
            int oldId = overflowsElements[tail];
            if (!(totalOverflows->at(oldId) - 1))
                totalOverflows->erase(oldId);
            else
                totalOverflows->insert(make_pair(oldId,totalOverflows->at(oldId) - 1));
            tail = (tail + 1) % maxOverflows;
            --overflowsNumber;
            indexTail = (indexTail + 1) % indexSize;
        }
    } catch (const out_of_range) {
    }

	int currBlock = (ceil((double)frameItems / (double)blockSize));
	currBlock = currBlock % (blocksNumber + 1);


    // Add item to RSS_CPP
    this->rss->update(item, wieght);

#ifdef ACC1_DEBUGGING
    if (testArr[currBlock - 1]) { //TODO for testing only
    	--testArr[currBlock - 1];
    	if(currBlock == 6) {
    		cout << "block 6 overflow" << endl;
    		item = firstOverflow;
    	}
    	if(currBlock == 6 && testArr[currBlock - 1] == 0) {
    		cout << "block 3 overflow" << endl;
    		item = thirdOverflow;
    	}
#else
    // overflow
   if ((this->rss->query(item) %threshold) == 0) {
#endif
    	head = (head + 1) % maxOverflows;
        overflowsElements[head] = item;
        if (idToIDx.find(item) == idToIDx.end()) {
        	idToIDx.insert(pair<unsigned int, unsigned int> (item, idx));
        	++idx;
        }
   	    int itemIdx = idToIDx.at(item);
        ++overflowsNumber;
        assert(overflowsNumber < maxOverflows);
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 1;
        unordered_map<int,int>::const_iterator foundedItem = totalOverflows->find(item);
        if (foundedItem == totalOverflows->end())
            totalOverflows->insert(make_pair<int,int>(item,1));
        else
            totalOverflows->at(item) = totalOverflows->at(item) + 1;
#ifdef ACC1_DEBUGGING
        printf("before overflowArrs currBlock: %d blocksNumber:%d itemIdx: %d\n", currBlock, blocksNumber, itemIdx);// TODO: testing
    	cout << "update level 0 for: " << itemIdx << " from: " << currBlock << " To: " << nextMulOf4(currBlock) << endl;
#endif
        /* Update Level 0 */
        for (int i = currBlock; i <= nextMulOf4(currBlock); ++i) {
            unordered_map<unsigned int, unsigned int>* blockMap = overflowedArrL0[i -1];

            if(blockMap->find(itemIdx) == blockMap->end())
            	blockMap->insert(pair<int, int> (itemIdx, 1));
            else
            	blockMap->at(itemIdx) = blockMap->at(itemIdx) + 1;

            overflowedArrL0[i - 1] = blockMap;
#ifdef ACC1_DEBUGGING
            cout << "overflowedArrL0["<< i -1 <<"] contains: ";
            for ( auto it = overflowedArrL0[i - 1]->begin(); it != overflowedArrL0[i - 1]->end(); ++it )
            	cout << " " << it->first << ":" << it->second;
            cout << endl;
#endif
        }
#ifdef ACC1_DEBUGGING
    	cout << "update level 1 for: " << itemIdx << " from: " << floor(log(currBlock)/log(4)) << " To: " << ceil(blocksNumber/4) - 1<< endl;
#endif
        /* Update Level 1 */
        for (int i = floor(log(currBlock)/log(4)); i < ceil(blocksNumber/4); ++i) {
            unordered_map<unsigned int, unsigned int>* blockMap = overflowedArrL1[i];

            if(blockMap->find(itemIdx) == blockMap->end())
            	blockMap->insert(pair<int, int> (itemIdx, 1));
            else
            	blockMap->at(itemIdx) = blockMap->at(itemIdx) + 1;
            overflowedArrL1[i] = blockMap;
#ifdef ACC1_DEBUGGING
            cout << "overflowedArrL1["<< i <<"] contains: ";
            for ( auto it = overflowedArrL1[i]->begin(); it != overflowedArrL1[i]->end(); ++it )
            	cout << " " << it->first << ":" << it->second;
            cout << endl;
#endif
        }

#ifdef ACC1_DEBUGGING
    	if (overflowsNumber == 1) {
    		firstOverflow = item;
    		printf("item in update: %d overflowesNumber: %d\n", item, overflowsNumber);
    		printf("idx: %d\n", idToIDx.at(item));
    	}
    	if (overflowsNumber == 3) {
    		thirdOverflow = item;
    	}
#endif
    }

    // New frame
    if (frameItems == windowSize) {
        frameItems = 0;
        rss->clear();
    }
}

double ACC1::query(unsigned int item)
{
    int minOverFlows;
    double rssEstimation = this->rss->query(item);

    unordered_map<int,int>::const_iterator foundedItem = totalOverflows->find(item);
    if (foundedItem == totalOverflows->end()) // item has no overflows
        minOverFlows = 0;
    else {
        minOverFlows = totalOverflows->at(item);
        //rssEstimation = (int) rssEstimation % (int) this->threshold;//TODO
    }

    rssEstimation = (int) rssEstimation % (int) this->threshold; //TODO
    return (this->threshold * (minOverFlows + 2 ) + rssEstimation);
	// return intervalQuery(lastOverflowed, 2, 9); //TODO: testing
}

double ACC1::intervalQuery(unsigned int item, int b1, int b2)
{
	int firstBlock = (int) ceil((double) b1 / (double) blockSize) % (int) (blocksNumber + 1);
	int secondBlock = (int) ceil((double) b2 / (double)blockSize) % (int) (blocksNumber + 1);
	int minOverFlows;
	int itemIdx;

#ifdef ACC1_DEBUGGING
	printf("item: %d\n", item);
	printf("second block: %d first block: %d\n", secondBlock, firstBlock);
#endif
	unordered_map<unsigned int,unsigned int>::const_iterator foundedItem = idToIDx.find(item);
	if (foundedItem == idToIDx.end()) // item has no overflows
		minOverFlows = 0;
	else {
		itemIdx = idToIDx.at(item);

#ifdef ACC1_DEBUGGING
		cout << "ItemIdx: " << itemIdx << " prevL1Block(secondBlock): " << prevL1Block(secondBlock) - 1 << " prevL1Block(firstBlock): " << prevL1Block(firstBlock) - 1<< endl;
#endif

		int overTillSecond = 0;
		int overTillFirst = 0;

		if (prevL1Block(secondBlock)) {
			if (overflowedArrL1[prevL1Block(secondBlock) - 1]->find(item) != overflowedArrL1[prevL1Block(secondBlock) - 1]->end())
				overTillSecond = overflowedArrL1[prevL1Block(secondBlock) - 1]->at(itemIdx);
		}

		if (4*prevL1Block(secondBlock) != secondBlock) {
			if (overflowedArrL0[secondBlock - 1]->find(item) != overflowedArrL0[secondBlock - 1]->end())
				overTillSecond += overflowedArrL0[secondBlock - 1]->at(itemIdx);
		}

		if (prevL1Block(firstBlock)) {
			if (overflowedArrL1[prevL1Block(firstBlock) - 1]->find(item) != overflowedArrL1[prevL1Block(firstBlock) - 1]->end())
				overTillFirst = overflowedArrL1[prevL1Block(firstBlock) - 1]->at(itemIdx);
		}


		if (4*prevL1Block(firstBlock)!= firstBlock) {
			if (overflowedArrL0[firstBlock - 1]->find(item) != overflowedArrL0[firstBlock - 1]->end())
				overTillFirst += overflowedArrL0[firstBlock - 1]->at(itemIdx);
		}

#ifdef ACC1_DEBUGGING
		cout << "overTillFirst: " << overTillFirst << " overTillSecond: " << overTillSecond << endl;
#endif
		minOverFlows = overTillSecond - overTillFirst;
	}
	// printf("minoverflowed: %d \n", minOverFlows); //TODO: testing
	return threshold * (minOverFlows + 1);
}

#ifdef ACC1_DEBUGGING
void ACC1::printHashMaps()
{
    for(int i = 0 ; i < blocksNumber ; i++) {
		cout << "overflowedArrL0["<< i <<"] contains: ";
		for ( auto it = overflowedArrL0[i]->begin(); it != overflowedArrL0[i]->end(); ++it )
			cout << " " << it->first << ":" << it->second;
		cout << endl;
    }

    for(int i = 0 ; i < ceil(blocksNumber/4) ; i++) {
		cout << "overflowedArrL1["<< i <<"] contains: ";
		for ( auto it = overflowedArrL1[i]->begin(); it != overflowedArrL1[i]->end(); ++it )
			cout << " " << it->first << ":" << it->second;
		cout << endl;
    }

    intervalQuery(firstOverflow, 4, 31);

}
#endif