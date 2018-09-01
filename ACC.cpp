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
    lastBlock = 1;

    index = new vector<int> (indexSize); // 0 means end of block
    overflowsElements = new unsigned int[maxOverflows];
    rss = new RSS_CPP(epsilon, m, gamma);//y
    threshold = ceil(windowSize*m*epsilon / 6.0); // W*M/k
    totalOverflows = new unordered_map<int, int> (maxOverflows);//B TODO: allocate static?

    overflowedArr = new unordered_map<unsigned int, unsigned int>*[blocksNumber];
    for(int i =0 ; i < blocksNumber ; i++)
    	overflowedArr[i] =  new unordered_map<unsigned int, unsigned int> (maxOverflows);
    incTable = new unordered_map<unsigned int, unsigned int>*[blocksNumber];
    for(int i =0 ; i < blocksNumber ; i++)
        incTable[i] =  new unordered_map<unsigned int, unsigned int> (maxOverflows);

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

void ACC::populateIncTable(unsigned int itemIdx)
{

    for(int l = 0; l < k; l++) {
            unordered_map<unsigned int, unsigned int>* foundedTable = incTable[l];

        if(foundedTable->find(itemIdx) == foundedTable->end())
            foundedTable->insert(pair<int, int> (itemIdx, 1));
        else
            foundedTable->at(itemIdx) = foundedTable->at(itemIdx) + 1;
    
        incTable[l] = foundedTable;
    }
}


double ACC::computeOverflowCount(unsigned int item)
{
    double k = 6. / this->epsilon;
    return (rss->query(item) / (this->m * (this->windowSize / k)));
}


void ACC::endBlock(int blockNumber)
{
    indexTail = (indexTail + 1) % indexSize;
    indexHead = (indexHead + 1) % indexSize;
    index->at(indexHead) = 0;

    int l = 0;
    while ((l < k-1) && (lastBlock % (pow(d, l+1)) == 0)) {
            incTable[l]->clear();
            ghostTables[l]->clear();
            ++l;
    }

    for (auto it = incTable->begin(); it != incTable->end(); ++it)
        skiplistMap[blockNumber]->insert(pair<unsigned int, unsigned int>(it->first, it->second));

    for (auto it = incTable->begin(); it != incTable->end(); ++it)
        skiplistMap[blockNumber]->insert(pair<unsigned int, unsigned int>(it->first, it->second));


    incTable->clear();

    populateSkipListLevels(blockNumber);

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
        	//*(overflowedArr[currBlock]) = *(overflowedArr[currBlock - 1]);

       		unordered_map<unsigned int, unsigned int>* prevL1Table = overflowedArr[currBlock - 1];
       		unordered_map<unsigned int, unsigned int>* mergedblockTables = overflowedArr[currBlock];

       		if (!prevL1Table->empty()) {
       			for (auto it = prevL1Table->begin(); it != prevL1Table->end(); ++it) {
       				if (mergedblockTables->find(it->first) != mergedblockTables->end()) {
       					int prevValue = mergedblockTables->find(it->first)->second;
       					mergedblockTables->at(it->first) = prevValue + it->second;
       				} else {
       					mergedblockTables->insert(pair<unsigned int, unsigned int>(it->first, it->second));
       				}
       			}
       		}
       		*overflowedArr[currBlock] = *mergedblockTables;
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

   	    int itemIdx = idToIDx.at(item);

        populateIncTable(unsigned int itemIdx);
    }

    /* Last frame in the current block */
    /* NEW BLOCK */
    if ((frameItems % blockSize) == 0) {
        endBlock(blockNumber);
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
    }

    rssEstimation = (int) rssEstimation % (int) this->threshold; //TODO
    return (this->threshold * (minOverFlows + 2 ) + rssEstimation);
}

double ACC::intervalQuery(unsigned int item, int b1, int b2)
{
    int first = ((int)(lastBlock - i) % (int) (blocksNumber+1));
    int last = ((int)(lastBlock - j) % (int) (blocksNumber+1));

	int minOverFlows;
	int itemIdx;

	if (idToIDx.find(item) == idToIDx.end()) // item has no overflows
		minOverFlows = 0;
	else {
		int overTillSecond = 0;
		int overTillFirst = 0;

		itemIdx = idToIDx.at(item);

		if (overflowedArr[last]->find(itemIdx) != overflowedArr[last]->end()) {
			overTillSecond = overflowedArr[last]->at(itemIdx);
		}

		if (overflowedArr[first]->find(itemIdx) != overflowedArr[first]->end()) {
			overTillFirst = overflowedArr[first]->at(itemIdx);
		}
		minOverFlows = overTillSecond - overTillFirst;
	}

	return threshold * (minOverFlows + 1);
}
