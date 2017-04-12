//
//  WRSS.cpp
//
//
//  Created by Rana Shahout on 3/15/17.
//
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "WRSS.hpp"

using namespace std;

static inline int computeBlockLevels(int n)
{
	int trailingZeros;

	for (trailingZeros = 0; !(n&1); n >>=1) ++trailingZeros;

	return trailingZeros + 1;
}

void WRSS::populateSkipList(int blockNumber, int itemIdx)
{

	/* Level 0*/
	Block block_0 = Block(blockNumber,0);
    unordered_map<Block, unordered_map<int, int> >::const_iterator foundedTable = skiplistMap.find(block_0);

    if (foundedTable == skiplistMap.end()) {
    	unordered_map<int, int> idxToCount;
    	idxToCount.insert(pair<int, int> (itemIdx, 1));
    	skiplistMap.insert(pair<Block, unordered_map<int, int> > (block_0, idxToCount));
    } else {
    	unordered_map<int, int> blockItemMap = foundedTable->second;
    	unordered_map<int, int>::const_iterator foundedItem = blockItemMap.find(itemIdx);

    	if(foundedItem == blockItemMap.end())
    		blockItemMap.insert(pair<int, int> (itemIdx, 1));
    	else
    		blockItemMap.at(itemIdx) = blockItemMap.at(itemIdx) + 1;
    }

    /* Level >=1 */

    for (int level = 1; level < computeBlockLevels(blockNumber); ++level) {
    	Block block = Block(blockNumber, level);
    	Block prevLevelBlock = Block(blockNumber, level -1);
    	Block prevBlock = Block(blockNumber -1, level -1);
        unordered_map<Block, unordered_map<int, int> >::const_iterator prevLevelTableItr = skiplistMap.find(prevLevelBlock);
        unordered_map<Block, unordered_map<int, int> >::const_iterator prevBlockTableItr = skiplistMap.find(prevBlock);

        unordered_map<int, int> prevBlockTable = prevBlockTableItr->second;

    	unordered_map<int, int> mergedblockTables = prevLevelTableItr->second;
    	mergedblockTables.insert(prevBlockTable.begin(), prevBlockTable.end());
    	skiplistMap.insert(pair<Block, unordered_map<int, int> > (block, mergedblockTables)); // insert do nothing if the key already exist
    }


}

WRSS::WRSS(int windowSize, double gamma, int m, double epsilon)
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

    /* Query Interval Section */
    skiplistSize = ceil(maxOverflows * log (maxOverflows));

}

WRSS::~WRSS()
{
    delete(totalOverflows);
    delete(rss);
    delete [] overflowsElements;
    delete(index);
}

double WRSS::computeOverflowCount(unsigned int item)
{
    double k = 4. / this->epsilon;
    return (rss->query(item) / (this->m * (this->windowSize / k)));
}

void WRSS::update(unsigned int item, int wieght)
{

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

    // Add item to RSS_CPP
    int prevOverflowCount = computeOverflowCount(item);
    this->rss->update(item, wieght);
    int currOverflowCount = computeOverflowCount(item);

    // overflow
    if (currOverflowCount > prevOverflowCount) {
    	printf("OVERFLOW!! %u %d  \n", item, wieght);
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
        populateSkipList(frameItems/blockSize, idToIDx.at(item));
    }

    // New frame
    if (frameItems == windowSize) {
        frameItems = 0;
        rss->clear();
    }
}

double WRSS::query(unsigned int item)
{
    int minOverFlows;
    double rssEstimation = this->rss->query(item);

    printf("WRSS item: %u\n", item);
    printf("rssEstimation: %f \n", rssEstimation);

    unordered_map<int,int>::const_iterator foundedItem = totalOverflows->find(item);
    if (foundedItem == totalOverflows->end()) // item has no oveflows
        minOverFlows = 0;
    else {
        minOverFlows = totalOverflows->at(item);
        //rssEstimation = (int) rssEstimation % (int) this->threshold;//TODO
    }

    rssEstimation = (int) rssEstimation % (int) this->threshold; //TODO
    return (this->threshold * (minOverFlows + 2 ) + rssEstimation);
}
