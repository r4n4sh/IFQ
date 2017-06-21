/*
 * ACC_K.cpp
 *
 *  Created on: May 9, 2017
 *      Author: ranashahout
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "ACC_K.hpp"

using namespace std;

#ifdef ACC_K_DEBUGGING
int testArr[24] = {1,2,3,0,0,1,1,0,1,2,3,0,0,1,1,0,1,2,3,0,0,1,1,0};
unsigned int firstOverflow;
unsigned int thirdOverflow;
#endif

ACC_K::ACC_K(unsigned int windowSize, float gamma, unsigned int m, float epsilon, unsigned int k)
{
	this->k = k;
	idx = 0;
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
    this->step = floor(pow(blocksNumber, (1.0/float(k))));
    maxOverflows = min(blocksNumber * 2, windowSize);
    indexSize = maxOverflows + blocksNumber;
    head = maxOverflows - 1;
    indexHead = blocksNumber - 1;
    index = new vector<int> (indexSize); // 0 means end of block
    overflowsElements = new unsigned int[maxOverflows];
    rss = new RSS_CPP(epsilon, m, gamma);//y
    threshold = ceil(windowSize*m*epsilon / 6.0); // W*M/k
    totalOverflows = new unordered_map<int, int> (maxOverflows);//B TODO: allocate static?
    overflowedArrLevels = new unordered_map<unsigned int, unsigned int>**[k];

    unsigned int level_size = blocksNumber;
    for(int level = 0; level < k; level++) {
#ifdef ACC_K_DEBUGGING
        cout << " window_size: " << windowSize << " step: " << step <<" blockNumber: " << blocksNumber << " blockSize: "<< blockSize << " level: " << level << " size: " << level_size<< endl;
#endif
        overflowedArrLevels[level] = new unordered_map<unsigned int, unsigned int>* [level_size];
        for(int i = 0 ; i < level_size ; i++)
        	overflowedArrLevels[level][i] = new unordered_map<unsigned int, unsigned int> (maxOverflows);

        level_size = ceil(double(level_size / step));
    }
}

ACC_K::~ACC_K()
{
    int level_size = blocksNumber;

    for(int level = 0; level < k; level++) {
        for(int i = 0 ; i < level_size ; i++)
        	delete (overflowedArrLevels[level][i]);
        delete(overflowedArrLevels[level]);
        level_size = ceil(double(level_size / step));
    }

    delete(overflowedArrLevels);
    delete(totalOverflows);
    delete(rss);
    delete [] overflowsElements;
    delete(index);
}

void ACC_K::update(unsigned int item, int wieght)
{
	++frameItems;
#ifdef ACC_K_DEBUGGING
	cout << "*** UPDATE for frame: " << frameItems << "***" << endl;
#endif

	int currBlock = (floor((double)frameItems / (double)blockSize));

	currBlock = (currBlock % blocksNumber) + 1;

	/* New Block */
    if (((frameItems) % blockSize) == 0) {
#ifdef ACC_K_DEBUGGING
    	cout << "NEW BLOCK! " << currBlock << endl;
#endif

    	indexTail = (indexTail + 1) % indexSize;
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 0;

        if (currBlock >= 1 && currBlock <= blocksNumber) {
        	for (int level = 1; level < k - 2; ++level) {
        		if (currBlock % (int)pow(step, level + 1) != 0 && (currBlock - pow(step, level)) % (int)pow(step, level + 1) != 0 && ((currBlock) / pow(step, pow(step, level)) ) != 1) {
#ifdef ACC_K_DEBUGGING
        			cout << "update level: " << k -1 << " for currBlock: " << currBlock << endl;
        			cout << "copy level: " << k -1 << " table: " << (int)((currBlock)/ pow(step, k -1) ) - 2 << " to table: " << (int)((currBlock) / pow(step, k -1) ) - 1 << endl;
#endif
					unordered_map<unsigned int, unsigned int>* prevL1Table = overflowedArrLevels[k -1][(int)((currBlock)/ pow(step, level) ) - 2];
					unordered_map<unsigned int, unsigned int>* mergedblockTables = overflowedArrLevels[k -1][(int)((currBlock) / pow(step, level) ) - 1];

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
					*overflowedArrLevels[level][((int) (currBlock/pow(step, level))) - 1] = *mergedblockTables;
        		}
        	}
        	// Always copy tables from previous block in level k -1

       		if((currBlock) % (int)pow(step, k -1) == 0 && ((currBlock) / pow(step, k -1) ) != 1) {
#ifdef ACC_K_DEBUGGING
        	cout << "update level: " << k -1 << " for currBlock: " << currBlock << endl;
   			cout << "copy level: " << k -1 << " table: " << (int)((currBlock)/ pow(step, k -1) ) - 2 << " to table: " << (int)((currBlock) / pow(step, k -1) ) - 1 << endl;
#endif
           		unordered_map<unsigned int, unsigned int>* prevL1Table = overflowedArrLevels[k -1][(int)((currBlock)/ pow(step, k -1) ) - 2];
           		unordered_map<unsigned int, unsigned int>* mergedblockTables = overflowedArrLevels[k -1][(int)((currBlock) / pow(step, k -1) ) - 1];

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
           		*overflowedArrLevels[k -1][((int) (currBlock/pow(step, k -1))) - 1] = *mergedblockTables;
       		}
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

#ifdef ACC_K_DEBUGGING
    if (testArr[currBlock - 1]) { //TODO for testing only
    	--testArr[currBlock - 1];
    	if(currBlock == 1) {
    		cout << "block 1 overflow" << endl;
    		item = firstOverflow;
    	}
/*
    	if(currBlock == 6 && testArr[currBlock - 1] == 0) {
    		cout << "block 3 overflow" << endl;
    		item = thirdOverflow;
    	}
*/
    	if(currBlock == 3) {
    		cout << "block 3 overflow" << endl;
    		item = firstOverflow;
    	}

    	if(currBlock == 7) {
    		cout << "block 7 overflow" << endl;
    		item = thirdOverflow;
    	}

    	printf("OVERFLOW!! %u %d  \n", item, wieght);


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
        assert(overflowsNumber <= maxOverflows);
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 1;
        unordered_map<int,int>::const_iterator foundedItem = totalOverflows->find(item);
        if (foundedItem == totalOverflows->end())
            totalOverflows->insert(make_pair<int,int>(item,1));
        else
            totalOverflows->at(item) = totalOverflows->at(item) + 1;
#ifdef ACC_K_DEBUGGING
        printf("before overflowArrs currBlock: %d blocksNumber:%d itemIdx: %d\n", currBlock, blocksNumber, itemIdx);// TODO: testing
    	//cout << "update level 0 for: " << itemIdx << " from: " << currBlock << " To: " << nextMulOf4(currBlock) << endl;
#endif

    	for (int level = 0; level < k; level++) {
    		int b = ceil(double(currBlock) / pow(step, level));
            unordered_map<unsigned int, unsigned int>* blockMapL = overflowedArrLevels[level][b -1];

            if(blockMapL->find(itemIdx) == blockMapL->end())
            	blockMapL->insert(pair<int, int> (itemIdx, 1));
            else
            	blockMapL->at(itemIdx) = blockMapL->at(itemIdx) + 1;

            //overflowedArrLevels[level][b -1] = blockMapL;
    	}

#ifdef ACC_K_DEBUGGING
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

double ACC_K::query(unsigned int item)
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


unsigned int ACC_K::withinFrameFrequency(unsigned int required_block, int itemIdx)
{
       int l_min = 1;
       int l_max;
       int sum = 0;

       while (required_block % (int)pow(step, l_min) == 0){
               ++l_min;
       }
       l_min -= 1;
       if (l_min == k -1)
               l_max = k -1;
       else
               l_max = ceil(log(required_block) / log(step)) - 1;
       int block = required_block;
       int l = l_min;
#ifdef ACC_K_DEBUGGING
       cout << "second Block: l_min: " << l_min << " l_max: " << l_max << endl;
#endif
       while (l <= l_max) {
#ifdef ACC_K_DEBUGGING
               cout << "lmin: " << l << " l_max: " << l_max << " block: " << block << endl;
               cout << "table [" << l << "]" << "[" <<  (int)((block - 1)/pow(step, l)) << "]" << endl;
#endif
               unordered_map<unsigned int, unsigned int>* table = overflowedArrLevels[l][(int)((block - 1)/pow(step, l))];
               if (table->find(itemIdx) != table->end()) {
                       sum += table->at(itemIdx);
               }
               block -= pow(step, l);
               l++;
       }
       return sum;
}

double ACC_K::intervalQuery(unsigned int item, int b1, int b2)
{
	int firstBlock = ((int) ceil((double) b2 / (double) blockSize) % (int) blocksNumber) + 1;
	int secondBlock = ((int) floor((double) b1 / (double)blockSize) % (int) blocksNumber) + 1;
	int sum = 0;
	int sum_till_second = 0;
	int sum_till_first = 0;
	int itemIdx;

#ifdef ACC_K_DEBUGGING
	printf("item: %d\n", item);
	printf("second block: %d first block: %d\n", secondBlock, firstBlock);
#endif
	unordered_map<unsigned int,unsigned int>::const_iterator foundedItem = idToIDx.find(item);
	if (foundedItem == idToIDx.end()) // item has no overflows
		sum = 0;
	else {
		itemIdx = idToIDx.at(item);

		if (firstBlock <= secondBlock) {
			sum_till_second = withinFrameFrequency(secondBlock, itemIdx);
			sum_till_first = withinFrameFrequency(firstBlock, itemIdx);
			sum = sum_till_second - sum_till_first;
		}
	}
	return threshold * (sum + 1);
}

#ifdef ACC_K_DEBUGGING
void ACC_K::printHashMaps()
{
	unsigned int level_size = blocksNumber;

    for(int level = 0 ; level < k ; level++) {
        for(int i = 0 ; i < level_size ; i++) {
    		cout << "overflowedArrLevels["<< level<< "]["<< i <<"] contains: ";
    		for ( auto it = overflowedArrLevels[level][i]->begin(); it != overflowedArrLevels[level][i]->end(); ++it )
    			cout << " " << it->first << ":" << it->second;
    		cout << endl;

        }
        level_size = ceil(double(level_size / step));
    }

    intervalQuery(firstOverflow, 41, 11);
}

#endif




