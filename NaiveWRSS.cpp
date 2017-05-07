//
//  NaiveWRSS.cpp
//
//
//  Created by Rana Shahout on 3/15/17.
//
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "NaiveWRSS.hpp"

using namespace std;

/*
int testArr[8] = {1,2,2,0,0,0,0,0};
unsigned int lastOverflowed;
*/ // TODO: testing

NaiveWRSS::NaiveWRSS(unsigned int windowSize, float gamma, unsigned int m, float epsilon)
{
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
    threshold = windowSize*m*epsilon / 4.; // W*M/k
    totalOverflows = new unordered_map<int, int> (maxOverflows);//B TODO: allocate static?
    overflowedArr = new int*[blocksNumber];

    for (int i = 0; i < blocksNumber; i++)
    	overflowedArr[i] = new int[maxOverflows];

    for (int i = 0; i < blocksNumber *  maxOverflows; i++)
    	*((int*)overflowedArr + i) = 0;

}

NaiveWRSS::~NaiveWRSS()
{
    for (int i = 0; i < blocksNumber; i++)
    	delete[] overflowedArr[i];

    delete(overflowedArr);
    delete(totalOverflows);
    delete(rss);
    delete [] overflowsElements;
    delete(index);
}

double NaiveWRSS::computeOverflowCount(unsigned int item)
{
    double k = 4. / this->epsilon;
    return (rss->query(item) / (this->m * (this->windowSize / k)));
}

void NaiveWRSS::update(unsigned int item, int wieght)
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

    int currBlock = int (floor(frameItems / blockSize)) % (int) blocksNumber;

    // Add item to RSS_CPP
    int prevQuery = this->rss->query(item);
    this->rss->update(item, wieght);

    /*
    if (testArr[currBlock]) { //TODO for testing only
    	--testArr[currBlock];
    */
    // overflow
    if ((prevQuery%threshold) + wieght > threshold) {
        head = (head + 1) % maxOverflows;
        overflowsElements[head] = item;
        idToIDx.insert(pair<int, int> (item, overflowsNumber));
        ++overflowsNumber;
        assert(overflowsNumber < maxOverflows);
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 1;
        unordered_map<int,int>::const_iterator foundedItem = totalOverflows->find(item);
        if (foundedItem == totalOverflows->end())
            totalOverflows->insert(make_pair<int,int>(item,1));
        else
            totalOverflows->at(item) = totalOverflows->at(item) + 1;

        /* printf("before overflowArr currBlock: %d blocksNumber:%d overflowsNumber: %d\n", currBlock, blocksNumber, overflowsNumber); */ // TODO: testing

        for (int i = currBlock; i < blocksNumber; ++i) {
        	*((int*)overflowedArr + maxOverflows * i + (overflowsNumber - 1)) = *((int*)overflowedArr + maxOverflows * i + (overflowsNumber - 1)) + 1;

        	/* printf("overflowedArr[%d][%d] = %d \n ",i, overflowsNumber -1 , *((int*)overflowedArr + maxOverflows * i + (overflowsNumber - 1))); */ //TODO: testing

        }

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

double NaiveWRSS::query(unsigned int item)
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

double NaiveWRSS::intervalQuery(unsigned int item, int b1, int b2)
{

	int firstBlock = (int) floor(b1 / blockSize) % (int) blocksNumber;
	int secondBlock = (int) floor(b2 / blockSize) % (int) blocksNumber;
	int minOverFlows;
	int itemIdx;

	/*
	printf("item: %d\n", item);
	printf("second block: %d first block: %d\n", secondBlock, firstBlock); */ //TODO: testing
	 unordered_map<int,int>::const_iterator foundedItem = idToIDx.find(item);
	    if (foundedItem == idToIDx.end()) // item has no overflows
	        minOverFlows = 0;
	    else {
	    	itemIdx = idToIDx.at(item);
	    	/*
	    	printf("itemIDx: %d\n", itemIdx);
	    	printf("overflowed in second: %d \n", *((int *)overflowedArr + maxOverflows * secondBlock + itemIdx));
	    	printf("overflowed in first: %d \n", *((int *)overflowedArr + maxOverflows * firstBlock + itemIdx));
			*/ //TODO: testing
			minOverFlows = *((int *)overflowedArr + maxOverflows * secondBlock + itemIdx)  - *((int *)overflowedArr + maxOverflows * firstBlock + itemIdx);
	    }

	    // printf("minoverflowed: %d \n", minOverFlows); //TODO: testing
	    return threshold * (minOverFlows + 2 );
}
