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
#ifdef DEBUGGING2
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

void HIT::populateCurrSkipListLevel(unsigned int blockNumber, unsigned int itemIdx)
{
     unsigned int block_l;
#ifdef DEBUGGING2
     cout << "populateCurrSkipListLevel Levels:" << computeBlockLevels(blockNumber) << " blockNumber: " << blockNumber << endl;
#endif
     for (int level = 0; level < computeBlockLevels(blockNumber); ++level) {

    	 block_l = levelToidx[level] + ((blockNumber - 1) >> level);

#ifdef DEBUGGING
    cout << "populateCurrSkipListLevel:: Populate skiplist level "<< level << " for block: " << blockNumber << " item index: " << itemIdx << endl;
    cout << "index in skiplist: " << block_l << endl;
#endif

         unordered_map<unsigned int, unsigned int>* foundedTable = skiplistMap[block_l];

         if(foundedTable->find(itemIdx) == foundedTable->end())
        	 foundedTable->insert(pair<int, int> (itemIdx, 1));
         else
        	 foundedTable->at(itemIdx) = foundedTable->at(itemIdx) + 1;

         skiplistMap[block_l] = foundedTable;


#ifdef DEBUGGING
     printf("block: %d level: %d indexofSkiplist: %d", blockNumber, level, block_l);
     cout << "foundedTable contains: ";
     for (auto it = foundedTable->begin(); it != foundedTable->end(); ++it )
       cout << " " << it->first << ":" << it->second;
     cout << endl;
#endif
     }
}


void HIT::populateSkipListLevels(unsigned int blockNumber)
{
     /* Level >=1 */
    for (int level = 1; level < computeBlockLevels(blockNumber); ++level) {

    	unsigned int block = levelToidx[level] + ((blockNumber - 1) >> level);
    	unsigned int prevLevelBlock = levelToidx[level -1] + ((blockNumber - 1) >>  (level -1));
    	unsigned int prevBlock = levelToidx[level -1] + (((blockNumber - 1) - (1 << (level - 1))) >> (level -1));

#ifdef DEBUGGING
    	cout << "populateSkipListLevels block: " << blockNumber << "level: " << level;
    	cout << " = " << "(" << blockNumber << "," << level -1 << ")" << " +( " << (blockNumber - (1 << (level - 1))) << "," << level -1 << ")" << endl;
    	cout << "index of block: " << block << " index of prev level block: " << prevLevelBlock << " index of prev Block: " << prevBlock << endl;
#endif
    	unordered_map<unsigned int, unsigned int>* prevLevelTable = skiplistMap[prevLevelBlock];
    	unordered_map<unsigned int, unsigned int>* prevBlockTable = skiplistMap[prevBlock];

    	unordered_map<unsigned int, unsigned int>* mergedblockTables = new unordered_map<unsigned int, unsigned int> (maxOverflows);

    	if ( prevBlockTable->empty()) {
    		if (prevLevelTable->empty()) {
#ifdef DEBUGGING
    			cout << "BOTH prev block and prev level are empty!!" << endl;
#endif
    			continue;
    		}
#ifdef DEBUGGING
    		cout << "prev block table is empty but prev level is Not" << endl;
#endif

    		for (auto it = prevLevelTable->begin(); it != prevLevelTable->end(); ++it) {
    			mergedblockTables->insert(pair<unsigned int, unsigned int>(it->first, it->second));
    		}

    	} else {


    		for (auto it = prevBlockTable->begin(); it != prevBlockTable->end(); ++it) {
    			mergedblockTables->insert(pair<unsigned int, unsigned int>(it->first, it->second));
    		}


#ifdef DEBUGGING
            cout << "prevBlockTable contains:";
            for ( auto it = mergedblockTables->begin(); it != mergedblockTables->end(); ++it )
            	cout << " " << it->first << ":" << it->second;
            cout << endl;
#endif
    		if (!prevLevelTable->empty()) {
        		for (auto it = prevLevelTable->begin(); it != prevLevelTable->end(); ++it) {
        			if (mergedblockTables->find(it->first) != mergedblockTables->end()) {
        				int prevValue = mergedblockTables->find(it->first)->second;
        				mergedblockTables->at(it->first) = prevValue + it->second;
        			} else {
        				mergedblockTables->insert(pair<unsigned int, unsigned int>(it->first, it->second));
        			}
        		}
#ifdef DEBUGGING
                cout << "prevLevelTable contains:";
                for ( auto it = prevLevelTable->begin(); it != prevLevelTable->end(); ++it )
                	cout << " " << it->first << ":" << it->second;
                cout << endl;
#endif
    		}

    	}
        skiplistMap[block] = mergedblockTables;
#ifdef DEBUGGING
        cout << "mergedblockTables contains:";
        for ( auto it = mergedblockTables->begin(); it != mergedblockTables->end(); ++it )
        	cout << " " << it->first << ":" << it->second;
        cout << endl;
#endif

#ifdef DEBUGGING
        unordered_map<unsigned int, unsigned int>* tmpfoundedTable = skiplistMap[1938];

		printf("!!!!! block: 64 level: 5 ");
     cout << "!!!!! contains: ";
     for (auto it = tmpfoundedTable->begin(); it != tmpfoundedTable->end(); ++it )
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
    this->blockSize = ceil((windowSize * epsilon)/6); // W/k; k= 6/epsilon
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
    rss = new RSS_CPP(epsilon/6., m, gamma);//y
    threshold = ceil(windowSize*m*epsilon / 6.f); // W*M/k
    totalOverflows = new unordered_map<unsigned int, unsigned int> (maxOverflows);//B TODO: allocate static?

    /* Query Interval Section */
    skiplistSize = 2*blocksNumber -1;

    skiplistMap = new unordered_map<unsigned int, unsigned int>*[skiplistSize];

    for(int i = 0 ; i < skiplistSize ; i++)
    	skiplistMap[i] =  new unordered_map<unsigned int, unsigned int> (maxOverflows);

    int levelToidxSize = ceil(log2(blocksNumber));
    levelToidx = new unsigned int[levelToidxSize + 1];

    for(int i = 0 ; i < levelToidxSize ; i++) {
    	levelToidx[i] =  (2 * blocksNumber) - (blocksNumber * pow(2, 1-i));
#ifdef DEBUGGING
        cout << "blocksNumber: " << blocksNumber << " blockSize: " << blockSize << " window: " << windowSize << endl;
    	cout << "levelToidx["<<i <<"]: " << levelToidx[i] << endl;
#endif
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

void HIT::update(unsigned int item, int wieght)
{
	++frameItems;
#ifdef DEBUGGING
	cout << "*** UPDATE for frame: " << frameItems << "***" << endl;
#endif
	int blockNumber = ceil(((double)frameItems / (double)blockSize));

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

#ifdef DEBUGGING2
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
    	++idx;
#else
    // overflow
    //	printf("IS OVERFLOW? : %d\n", this->rss->query(item));
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
#ifdef DEBUGGING2
        populateCurrSkipListLevel(blockNumber, idx);
#else
        populateCurrSkipListLevel(blockNumber, idToIDx.at(item));
#endif
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

unsigned int HIT::partialIntervalQuery(unsigned int itemIdx, unsigned int b2, unsigned int b1)
{
	int b = b1;
    unsigned int count = 0;
    int d = 1 + (b1 - b2) % (blocksNumber + 1);
#ifdef DEBUGGING
     cout << "d: " << d << " b: " << b << endl;
#endif

    while (d > 0) {
    	int level = min((computeBlockLevels(b) - 1), (int)floor(log2(d)));
#ifdef DEBUGGING
    	cout << "level: " << level <<" d: "<< d << endl;
#endif

    	unsigned int block = levelToidx[level] +  ((b -1) >> level);

#ifdef DEBUGGING
    	cout << "level: " << level << " block: " << b << " value: " << skiplistMap[block]->at(itemIdx) << " index of skiplist: " << block << endl;
#endif


        unordered_map<unsigned int, unsigned int>* foundedTable = skiplistMap[block];


#ifdef DEBUGGING
		printf("### block: %d level: %d ", b, level);
        cout << "#### contains: ";
        for (auto it = foundedTable->begin(); it != foundedTable->end(); ++it )
        	cout << " " << it->first << ":" << it->second;
        cout << endl;
#endif

		if(skiplistMap[block]->find(itemIdx) != skiplistMap[block]->end()) {
			count += skiplistMap[block]->at(itemIdx);
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
	int firstBlock = ((int) ceil((double) b2 / (double) blockSize) % (int) blocksNumber) + 1;
	int secondBlock = ((int) ceil((double) b1 / (double)blockSize) % (int) blocksNumber);
	unsigned int minOverFlows = 0;
	int itemIdx;
#ifdef DEBUGGING
	//printf("item: %d\n", item);
	printf("first block: %d, second block: %d\n", firstBlock, secondBlock); //TODO:
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
	//return threshold * (minOverFlows + 1); EMP_ERROR TEST
	return minOverFlows;
}

#ifdef DEBUGGING2
double HIT::testIntervalQuery(unsigned int b2, unsigned int b1)
{
	unsigned int overflowed = firstOverflowed;
	unordered_map<unsigned int,unsigned int>::const_iterator foundedItem = idToIDx.find(overflowed);
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

