//
//  HIT.cpp
//
//
//  Created by Rana Shahout on 3/15/17.
//
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "HIT.hpp"

using namespace std;

static inline int computeBlockLevels(unsigned int n)
{
	if (!n) return 1;
	int clz = __builtin_ctz(n);
	return clz + 1;
}

void HIT::populateIncTable(unsigned int itemIdx)
{

    unordered_map<unsigned int, unsigned int>* foundedTable = incTable;

    if(foundedTable->find(itemIdx) == foundedTable->end())
        foundedTable->insert(pair<int, int> (itemIdx, 1));
    else
        foundedTable->at(itemIdx) = foundedTable->at(itemIdx) + 1;
    
    incTable = foundedTable;
}


void HIT::populateSkipListLevels(unsigned int blockNumber)
{
     /* Level >=1 */
    for (int level = 1; level < computeBlockLevels(blockNumber); ++level) {

    	unsigned int block = levelToidx[level] + ((blockNumber) >> level);
    	unsigned int prevLevelBlock = levelToidx[level -1] + ((blockNumber) >>  (level -1));
    	unsigned int prevBlock = levelToidx[level -1] + (((blockNumber) - (1 << (level - 1))) >> (level -1));
    	unordered_map<unsigned int, unsigned int>* prevLevelTable = skiplistMap[prevLevelBlock];
    	unordered_map<unsigned int, unsigned int>* prevBlockTable = skiplistMap[prevBlock];

    	unordered_map<unsigned int, unsigned int>* mergedblockTables = new unordered_map<unsigned int, unsigned int> (maxOverflows);

    	if ( prevBlockTable->empty()) {
    		if (prevLevelTable->empty()) {
    			continue;
    		}

    		for (auto it = prevLevelTable->begin(); it != prevLevelTable->end(); ++it) {
    			mergedblockTables->insert(pair<unsigned int, unsigned int>(it->first, it->second));
    		}

    	} else {


    		for (auto it = prevBlockTable->begin(); it != prevBlockTable->end(); ++it) {
    			mergedblockTables->insert(pair<unsigned int, unsigned int>(it->first, it->second));
    		}


    		if (!prevLevelTable->empty()) {
        		for (auto it = prevLevelTable->begin(); it != prevLevelTable->end(); ++it) {
        			if (mergedblockTables->find(it->first) != mergedblockTables->end()) {
        				int prevValue = mergedblockTables->find(it->first)->second;
        				mergedblockTables->at(it->first) = prevValue + it->second;
        			} else {
        				mergedblockTables->insert(pair<unsigned int, unsigned int>(it->first, it->second));
        			}
        		}
    		}

    	}
        skiplistMap[block] = mergedblockTables;
    }
}

HIT::HIT(unsigned int windowSize, float gamma, unsigned int m, float epsilon)
{
    frameItems = 0;
    overflowsNumber = 0;
    tail = 0;
    idx = 0;
    indexTail = 0;
    this->blockSize = ceil((windowSize * epsilon)/6.f); // W/k; k= 6/epsilon
    this->windowSize = windowSize;
    this->blocksNumber = windowSize / blockSize;
    this->m = m;
    this->epsilon = epsilon;
    lastBlock = 1;

    maxOverflows = min(blocksNumber * 2, windowSize);
    indexSize = maxOverflows + blocksNumber;
    head = maxOverflows - 1;
    indexHead = blocksNumber - 1;
    index = new vector<int> (indexSize); // 0 means end of block
    overflowsElements = new unsigned int[maxOverflows];
    rss = new RSS_CPP(epsilon/6., m, gamma);//y
    totalOverflows = new unordered_map<unsigned int, unsigned int> (maxOverflows);//B TODO: allocate static?

    /* Query Interval Section */
    skiplistSize = 2*blocksNumber -1;

    skiplistMap = new unordered_map<unsigned int, unsigned int>*[skiplistSize];

    for(int i = 0 ; i < skiplistSize ; i++)
    	skiplistMap[i] =  new unordered_map<unsigned int, unsigned int> (maxOverflows);

    incTable = new unordered_map<unsigned int, unsigned int> (maxOverflows);

    int levelToidxSize = ceil(log2(blocksNumber));
    levelToidx = new unsigned int[levelToidxSize + 1];

    for(int i = 0 ; i < levelToidxSize ; i++) {
    	levelToidx[i] =  (2 * blocksNumber) - (blocksNumber * pow(2, 1-i));
    }
}

HIT::~HIT()
{
	delete[] levelToidx;
    for (int i = 0; i < skiplistSize; i++)
    	delete(skiplistMap[i]);
    delete[] skiplistMap;

    delete(totalOverflows);
    delete(rss);
    delete [] overflowsElements;
    delete(index);
}

int HIT::getIndexInskiplist(int blockNumber, int level)
{
    return levelToidx[level] +  ((blockNumber -1) >> level);

}
void HIT::endBlock(int blockNumber)
{

    lastBlock = 1 + (lastBlock % (int) blocksNumber);

    indexTail = (indexTail + 1) % indexSize;
    indexHead = (indexHead + 1) % indexSize;
    index->at(indexHead) = 0;

    for (auto it = incTable->begin(); it != incTable->end(); ++it)
        skiplistMap[blockNumber]->insert(pair<unsigned int, unsigned int>(it->first, it->second));

    incTable->clear();

    populateSkipListLevels(blockNumber);

}


void HIT::update(unsigned int item, int wieght)
{
	++frameItems;
	int blockNumber = ceil(((double)frameItems / (double)blockSize));



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

    // overflow
    //	printf("IS OVERFLOW? : %d\n", this->rss->query(item));
    if ((this->rss->query(item)%blockSize) == 0) {
        head = (head + 1) % maxOverflows;
        overflowsElements[head] = item;
        if (idToIDx.find(item) == idToIDx.end()) {
        	idToIDx.insert(pair<unsigned int, unsigned int> (item, idx));
        	++idx;
        }
        ++overflowsNumber;
        assert(overflowsNumber <= maxOverflows);
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 1;
        if (totalOverflows->find(item) == totalOverflows->end())
            totalOverflows->insert(make_pair<unsigned int, unsigned int>((unsigned int)item,1));
        else
            totalOverflows->at(item) = totalOverflows->at(item) + 1;
        populateIncTable(idToIDx.at(item));
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

double HIT::query(unsigned int item)
{
    int minOverFlows;
    double rssEstimation = this->rss->query(item);
    unordered_map<unsigned int,unsigned int>::const_iterator foundedItem = totalOverflows->find(item);
    if (foundedItem == totalOverflows->end()) // item has no oveflows
        minOverFlows = 0;
    else {
        minOverFlows = totalOverflows->at(item);
        //rssEstimation = (int) rssEstimation % (int) this->blockSize;//TODO
    }

    rssEstimation = (int) rssEstimation % (int) this->blockSize; //TODO
    return (this->blockSize * (minOverFlows + 2 ) + rssEstimation);
}

unsigned int HIT::partialIntervalQuery(unsigned int itemIdx, unsigned int b2, unsigned int b1)
{
	int b = b1;
    unsigned int count = 0;
    int d = 1 + (b1 - b2) % (blocksNumber + 1);

    while (d > 0) {
    	int level = min((computeBlockLevels(b) - 1), (int)floor(log2(d)));

    unsigned int block;
    if (b <= 0)
        block = levelToidx[level];
    else
        block = levelToidx[level] +  ((b) >> level);


        unordered_map<unsigned int, unsigned int>* foundedTable = skiplistMap[block];

		if(skiplistMap[block]->find(itemIdx) != skiplistMap[block]->end()) {
			count += skiplistMap[block]->at(itemIdx);
		}



    	int power_level = 1 << level;
    	d = d - power_level;
    	b = b - power_level;

    	if (b <= 0)
    		b = blocksNumber;
     }
     return count;
}

double HIT::intervalQuery(unsigned int item, int i, int j)
{
    int first = ((int)(lastBlock - i) % (int) (blocksNumber+1));
    int last = ((int)(lastBlock - j) % (int) (blocksNumber+1));

	unsigned int minOverFlows = 0;
	int itemIdx;
	unordered_map<unsigned int,unsigned int>::const_iterator foundedItem = idToIDx.find(item);
	if (foundedItem == idToIDx.end()) // item has no overflows
		minOverFlows = 0;
	else {
		itemIdx = idToIDx.at(item);
		minOverFlows = partialIntervalQuery(itemIdx, last, first);
	}
	return minOverFlows;
	//return block_size * (minOverFlows + 2);
}

double HIT::intervalFrequencyQuery(unsigned int item, int i, int j)
{
    return blockSize *(intervalQuery(item, ceil((double)i/(double)blockSize), floor((double)j/(double)blockSize)) + 2);
}

double HIT::intervalQueryTest(unsigned int item, int i, int j)
{
    unsigned int minOverFlows = 0;
    int itemIdx;
    unordered_map<unsigned int,unsigned int>::const_iterator foundedItem = idToIDx.find(item);
    if (foundedItem == idToIDx.end()) // item has no overflows
        minOverFlows = 0;
    else {
        itemIdx = idToIDx.at(item);
        minOverFlows = partialIntervalQuery(itemIdx, i, j);
    }
    return blockSize * (minOverFlows + 2);
}