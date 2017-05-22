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

#ifdef DEBUGGING
int testArr[8] = {1,2,1,2,1,0,1,1};
unsigned int firstOverflowed;
unsigned int secondtOverflowed;
unsigned int thirdOverflowed;
unsigned int fourthOverflowed;
unsigned int fifthOverflowed;
unsigned int sixthOverflowed;
#endif

static inline int computeBlockLevels(unsigned int n)
{
/*
	int trailingZeros;
	if (!n) return 1;
	for (trailingZeros = 0; !(n&1); n >>=1)
		++trailingZeros;
	return trailingZeros + 1;
*/

	if (!n) return 1;
	int clz = __builtin_ctz(n);
	return clz + 1;
}

void HIT::populateSkipListLevel_0(unsigned int blockNumber, unsigned int itemIdx)
{
     /* Level 0*/
     Block block_0 = Block(blockNumber,0);
    unordered_map<Block, unordered_map<unsigned int, unsigned int> >::const_iterator foundedTable = skiplistMap.find(block_0);
#ifdef DEBUGGING
    cout << "Populate skiplist level 0 for block: " << blockNumber << " item index: " << itemIdx << endl;
#endif
    if (foundedTable == skiplistMap.end()) {
     unordered_map<unsigned int, unsigned int> idxToCount;
     idxToCount.insert(pair<unsigned int, unsigned int> (itemIdx, 1));
     skiplistMap.insert(pair<Block, unordered_map<unsigned int, unsigned int> > (block_0, idxToCount));
    } else {
     unordered_map<unsigned int, unsigned int> blockItemMap = foundedTable->second;

     if(blockItemMap.find(itemIdx) == blockItemMap.end())
             blockItemMap.insert(pair<int, int> (itemIdx, 1));
     else
             blockItemMap.at(itemIdx) = blockItemMap.at(itemIdx) + 1;

     skiplistMap[block_0] = blockItemMap;
#ifdef DEBUGGING
     printf("block: %d level: 0 ", blockNumber);
     cout << "blockItemMap contains: ";
     for (auto it = blockItemMap.begin(); it != blockItemMap.end(); ++it )
       cout << " " << it->first << ":" << it->second;
     cout << endl;
#endif
    }
}

void HIT::populateSkipListLevels(unsigned int blockNumber)
{
     /* Level >=1 */
    for (int level = 1; level < computeBlockLevels(blockNumber); ++level) {

    	Block block = Block(blockNumber, level);
    	Block prevLevelBlock = Block(blockNumber, level -1);
    	Block prevBlock = Block(blockNumber - (1 << (level - 1)), level -1);
#ifdef DEBUGGING
    	cout << "populateSkipListLevels block: " << blockNumber << "level: " << level;
    	cout << " = " << "(" << prevLevelBlock.blockNumber << "," << prevLevelBlock.blockLevel << ")" << " +( " << prevBlock.blockNumber << "," << prevBlock.blockLevel << ")" << endl;
#endif
        unordered_map<Block, unordered_map<unsigned int, unsigned int> >::const_iterator prevLevelTableItr = skiplistMap.find(prevLevelBlock);
        unordered_map<Block, unordered_map<unsigned int, unsigned int> >::const_iterator prevBlockTableItr = skiplistMap.find(prevBlock);

    	unordered_map<unsigned int, unsigned int> mergedblockTables;

    	if ( prevBlockTableItr == skiplistMap.end() ) {
    		if (prevLevelTableItr == skiplistMap.end()) {
#ifdef DEBUGGING
    			cout << "BOTH prev block and prev level are empty!!" << endl;
#endif
    			continue;
    		}
#ifdef DEBUGGING
    		cout << "prev block table is empty but prev level is Not" << endl;
#endif
    		mergedblockTables = prevLevelTableItr->second;
    	} else {

    		mergedblockTables = prevBlockTableItr->second;

#ifdef DEBUGGING
            cout << "prevBlockTable contains:";
            for ( auto it = mergedblockTables.begin(); it != mergedblockTables.end(); ++it )
            	cout << " " << it->first << ":" << it->second;
            cout << endl;
#endif
    		if (prevLevelTableItr != skiplistMap.end()) {
    			unordered_map<unsigned int, unsigned int> prevLevelTable = prevLevelTableItr->second;
        		for (auto it = prevLevelTable.begin(); it != prevLevelTable.end(); ++it) {
        			if (mergedblockTables.find(it->first) != mergedblockTables.end()) {
        				int prevValue = mergedblockTables.find(it->first)->second;
        				mergedblockTables[it->first] = prevValue + it->second;
        			} else {
        				mergedblockTables.insert(pair<unsigned int, unsigned int>(it->first, it->second));
        			}
        		}
#ifdef DEBUGGING
                cout << "prevLevelTable contains:";
                for ( auto it = prevLevelTable.begin(); it != prevLevelTable.end(); ++it )
                	cout << " " << it->first << ":" << it->second;
                cout << endl;
#endif
    		}

    	}

        skiplistMap.insert(pair<Block, unordered_map<unsigned int, unsigned int> > (block, mergedblockTables)); // insert do nothing if the key already exist
#ifdef DEBUGGING
        cout << "mergedblockTables contains:";
        for ( auto it = mergedblockTables.begin(); it != mergedblockTables.end(); ++it )
        	cout << " " << it->first << ":" << it->second;
        cout << endl;
#endif
    }
}

HIT::HIT(unsigned int windowSize, float gamma, unsigned int m, float epsilon)
{
    frameItems = 0;
    overflowsNumber = 0;
    tail = 0;
    idx = 0;
    indexTail = 0;
    this->blockSize = ceil((windowSize * epsilon)/4); // W/k; k= 4/epsilon
    this->windowSize = windowSize;
    this->blocksNumber = windowSize / blockSize;
    this->m = m;
    this->epsilon = epsilon;

    maxOverflows = min(blocksNumber * 2, windowSize);
    indexSize = maxOverflows + blocksNumber;
    head = maxOverflows - 1;
    indexHead = blocksNumber - 1;
    index = new vector<int> (indexSize); // 0 means end of block
    overflowsElements = new unsigned int[maxOverflows];
    rss = new RSS_CPP(epsilon, m, gamma);//y
    threshold = ceil(windowSize*m*epsilon / 4.f); // W*M/k
    totalOverflows = new unordered_map<unsigned int, unsigned int> (maxOverflows);//B TODO: allocate static?

    /* Query Interval Section */
    skiplistSize = ceil(maxOverflows * log (maxOverflows));
}

HIT::~HIT()
{
    delete(totalOverflows);
    delete(rss);
    delete [] overflowsElements;
    delete(index);
}

void HIT::update(unsigned int item, int wieght)
{
	++frameItems;
#ifdef DEBUGGING
	cout << "*** UPDATE for frame: " << frameItems << "***" << endl;
#endif
	int blockNumber = (ceil((double)frameItems / (double)blockSize));

	blockNumber = blockNumber % (blocksNumber + 1);

	/* Last frame in the current block */
    if ((frameItems % blockSize) == 0) {
#ifdef DEBUGGING
    	cout << "NEW BLOCK! " << blockNumber << endl;
#endif
        indexTail = (indexTail + 1) % indexSize;
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 0;
        // Populate previous block levels
        if (blockNumber > 1) {
#ifdef DEBUGGING
        	cout << "populate levels of block number: " << blockNumber << endl;
#endif
        	populateSkipListLevels(blockNumber);
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

#ifdef DEBUGGING
    if (testArr[blockNumber -1]) {

    	cout << "overflow for block number: " << blockNumber << endl;
    	--testArr[blockNumber -1];

        if (blockNumber == 1)
        	firstOverflowed = item;
        else if (blockNumber == 2 && testArr[blockNumber -1] == 1){
        	secondtOverflowed = item;
        } else if (blockNumber == 2 && testArr[blockNumber -1] == 0){
        	thirdOverflowed = item;
        } else if (blockNumber == 3) {
        	item = firstOverflowed;
        } else if (blockNumber == 4 && testArr[blockNumber -1] == 1) {
        	item = firstOverflowed;
        } else if (blockNumber == 4 && testArr[blockNumber -1] == 0) {
        	fourthOverflowed = item;
        } else if (blockNumber == 5 ) {
        	item = secondtOverflowed;
        } else if (blockNumber == 7) {
        	item = thirdOverflowed;
        }else if (blockNumber == 8) {
        	sixthOverflowed = item;
        }
    	printf("OVERFLOW!! %u %d  \n", item, wieght);
#else
    // overflow
        if ((this->rss->query(item)%threshold) == 0) {
#endif
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

        populateSkipListLevel_0(blockNumber, idToIDx.at(item));
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
#ifdef DEBUGGING
    printf("HIT item: %u\n", item);
    printf("rssEstimation: %f \n", rssEstimation);
#endif
    unordered_map<unsigned int,unsigned int>::const_iterator foundedItem = totalOverflows->find(item);
    if (foundedItem == totalOverflows->end()) // item has no oveflows
        minOverFlows = 0;
    else {
        minOverFlows = totalOverflows->at(item);
        //rssEstimation = (int) rssEstimation % (int) this->threshold;//TODO
    }

    rssEstimation = (int) rssEstimation % (int) this->threshold; //TODO
    return (this->threshold * (minOverFlows + 2 ) + rssEstimation);
}

double HIT::partialIntervalQuery(unsigned int itemIdx, unsigned int b2, unsigned int b1)
{
	unsigned int b = b1;
    double count = 0;
    unsigned int d = 1 + (b1 - b2) % (blocksNumber + 1);
    int level_number = 0;
#ifdef DEBUGGING
     cout << "d: " << d << " b: " << b << endl;
#endif
    while (d > 0) {
    	int level = min((computeBlockLevels(b) - 1), (int)floor(log2(d)));

    	Block block = Block(b,level);
#ifdef DEBUGGING
    	cout << "level: " << level << " block: " << b << endl;
#endif

    	unordered_map<Block, unordered_map<unsigned int, unsigned int> >::const_iterator foundedTable = skiplistMap.find(block);

    	if (foundedTable == skiplistMap.end()) {
#ifdef DEBUGGING
    		 cout << "The item has not overflowed!! in block: " << block.blockNumber << " level: " << block.blockLevel << endl;
#endif
    	} else {
#ifdef DEBUGGING
    		cout << "Overflowed item " << "block: " << b << "level: " << level << "d: " << d << endl;
#endif
    		unordered_map<unsigned int, unsigned int>* blockItemMap = &(skiplistMap.at(block));

    		if(blockItemMap->find(itemIdx) == blockItemMap->end()) {
    		 //TODO
    			printf("Not founded item at block: %d level: %d\n", block.blockNumber, block.blockLevel);
    		}
    		else
    			count += blockItemMap->at(itemIdx);

    	}

    	int power_level = 1 << level;
    	d = d - power_level;
    	b = b - power_level;

    	if (b == 0)
    		b = blocksNumber;
     }
     return count;
}

double HIT::intervalQuery(unsigned int item, int b2, int b1)
{

	int firstBlock = (int) ceil((double) b2 / (double) blockSize) % (int) (blocksNumber + 1);
	int secondBlock = (int) ceil((double) b1 / (double)blockSize) % (int) (blocksNumber + 1);
	double minOverFlows = 0;
	int itemIdx;
#ifdef DEBUGGING
	//printf("item: %d\n", item);
	printf("second block: %d first block: %d\n", secondBlock, firstBlock); //TODO:
#endif
	unordered_map<unsigned int,unsigned int>::const_iterator foundedItem = idToIDx.find(item);
	if (foundedItem == idToIDx.end()) // item has no overflows
		minOverFlows = 0;
	else {
		itemIdx = idToIDx.at(item);
		minOverFlows = partialIntervalQuery(itemIdx, firstBlock, secondBlock);

#ifdef DEBUGGING
		cout << "number of overflow: " << minOverFlows << endl;
#endif
	}
	return threshold * (minOverFlows + 1);
}

#ifdef DEBUGGING
double HIT::testIntervalQuery(unsigned int b2, unsigned int b1)
{
	unsigned int overflowed = firstOverflowed;
	unordered_map<int,int>::const_iterator foundedItem = idToIDx.find(overflowed);
	if (foundedItem == idToIDx.end()) { // item has no overflows
		cout << "item has no overflow: "<< overflowed << endl;
		return 0;
	}
	else {
		int itemIdx = idToIDx.at(overflowed);
		cout << "item: " << overflowed << "overflowed at index: " << itemIdx << endl;
		return intervalQuery(overflowed, b2, b1);
	}
}
#endif

