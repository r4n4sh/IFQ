//
//  ACC.cpp
//
//
//  Created by Rana Shahout on 3/15/17.
//
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "ACC.hpp"

using namespace std;

/*
int testArr[8] = {1,2,2,0,0,0,0,0};
unsigned int lastOverflowed;
*/ // TODO: testing

ACC::ACC(unsigned int windowSize, float gamma, unsigned int m, float epsilon)
{
    frameItems = 0;
    overflowsNumber = 0;
    tail = 0;
    indexTail = 0;
    this->blockSize = ceil((windowSize * epsilon)/6.0); // W/k; k= 6/epsilon
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
    threshold = ceil(windowSize*m*epsilon / 6.0); // W*M/k
    totalOverflows = new unordered_map<int, int> (maxOverflows);//B TODO: allocate static?

    overflowedArr = new unordered_map<unsigned int, unsigned int>*[blocksNumber];
    for(int i =0 ; i < blocksNumber ; i++)
    	overflowedArr[i] =  new unordered_map<unsigned int, unsigned int> (maxOverflows);

}

ACC::~ACC()
{
    for (int i = 0; i < blocksNumber; i++)
    	delete(overflowedArr[i]);
    delete[] overflowedArr;

    delete(totalOverflows);
    delete(rss);
    delete [] overflowsElements;
    delete(index);
}

double ACC::computeOverflowCount(unsigned int item)
{
    double k = 6. / this->epsilon;
    return (rss->query(item) / (this->m * (this->windowSize / k)));
}

void ACC::update(unsigned int item, int wieght)
{
	++frameItems;

	int currBlock = (floor((double)frameItems / (double)blockSize));
	currBlock = (currBlock % (blocksNumber)) + 1;

	/* New Block */
    if (((frameItems) % blockSize) == 0) {
        indexTail = (indexTail + 1) % indexSize;
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 0;
        if (currBlock >= 1 && currBlock < blocksNumber) {
        	*(overflowedArr[currBlock]) = *(overflowedArr[currBlock - 1]);
        }

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


	// Add item to RSS_CPP
    this->rss->update(item, wieght);

    /*
    if (testArr[currBlock]) { //TODO for testing only
    	--testArr[currBlock];
    */
    // overflow
    if ((this->rss->query(item) %threshold) == 0) {
        head = (head + 1) % maxOverflows;
        overflowsElements[head] = item;
        idToIDx.insert(pair<int, int> (item, overflowsNumber));
        ++overflowsNumber;
        assert(overflowsNumber <= maxOverflows);
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 1;
        unordered_map<int,int>::const_iterator foundedItem = totalOverflows->find(item);
        if (foundedItem == totalOverflows->end())
            totalOverflows->insert(make_pair<int,int>(item,1));
        else
            totalOverflows->at(item) = totalOverflows->at(item) + 1;

        /* printf("before overflowArr currBlock: %d blocksNumber:%d overflowsNumber: %d\n", currBlock, blocksNumber, overflowsNumber); */ // TODO: testing

   	    int itemIdx = idToIDx.at(item);

/*
        for (int i = currBlock; i <= blocksNumber; ++i) {
            unordered_map<unsigned int, unsigned int>* blockMap = overflowedArr[i -1];

            if(blockMap->find(itemIdx) == blockMap->end())
            	blockMap->insert(pair<int, int> (itemIdx, 1));
            else
            	blockMap->at(itemIdx) = blockMap->at(itemIdx) + 1;

            overflowedArr[i - 1] = blockMap;
#ifdef ACC1_DEBUGGING
            cout << "overflowedArr["<< i -1 <<"] contains: ";
            for ( auto it = overflowedArr[i - 1]->begin(); it != overflowedArr[i - 1]->end(); ++it )
            	cout << " " << it->first << ":" << it->second;
            cout << endl;
#endif
        }
*/

        unordered_map<unsigned int, unsigned int>* blockMap = overflowedArr[currBlock -1];

        if(blockMap->find(itemIdx) == blockMap->end())
        	blockMap->insert(pair<int, int> (itemIdx, 1));
        else
        	blockMap->at(itemIdx) = blockMap->at(itemIdx) + 1;

        overflowedArr[currBlock - 1] = blockMap;

        /*
    	if (overflowsNumber == 2) {
    		lastOverflowed = item;
    		printf("item in update: %d overflowesNumber: %d\n", item, overflowsNumber);
    		printf("idx: %d\n", idToIDx.at(item));
    	} */ //TODO: testing
    }

    // New frame
    if (frameItems == windowSize) {
        frameItems = 0;
        rss->clear();
    }
}

double ACC::query(unsigned int item)
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

double ACC::intervalQuery(unsigned int item, int b1, int b2)
{
	int firstBlock = ((int) ceil((double) b1 / (double) blockSize) % (int) blocksNumber);
	int secondBlock = ((int) floor((double) b2 / (double)blockSize) % (int) blocksNumber);
	int minOverFlows;
	int itemIdx;

	if (idToIDx.find(item) == idToIDx.end()) // item has no overflows
		minOverFlows = 0;
	else {
		int overTillSecond = 0;
		int overTillFirst = 0;

		itemIdx = idToIDx.at(item);
		/*
		printf("itemIDx: %d\n", itemIdx);
		printf("overflowed in second: %d \n", *((int *)overflowedArr + maxOverflows * secondBlock + itemIdx));
		printf("overflowed in first: %d \n", *((int *)overflowedArr + maxOverflows * firstBlock + itemIdx));
		*/ //TODO: testing

		if (overflowedArr[secondBlock]->find(itemIdx) != overflowedArr[secondBlock]->end()) {
			overTillSecond = overflowedArr[secondBlock]->at(itemIdx);
		}

		if (overflowedArr[firstBlock]->find(itemIdx) != overflowedArr[firstBlock]->end()) {
			overTillFirst = overflowedArr[firstBlock]->at(itemIdx);
		}
		minOverFlows = overTillSecond - overTillFirst;
	}

	return threshold * (minOverFlows + 1);
}
