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

static inline int nextMulOf4(double number) {
	if (number == 1) return 4;
	return 4 * (int)ceil(number/4);

}

static inline int prevL1Block(int number) {
	return (int)floor(double (number)/4.0);
}

ACC1::ACC1(unsigned int windowSize, float gamma, unsigned int m, float epsilon)
{
	idx = 0;
    frameItems = 0;
    overflowsNumber = 0;
    tail = 0;
    indexTail = 0;
    this->blockSize = ceil((windowSize * epsilon)/8); // W/k; k= 8/epsilon
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
    threshold = ceil(windowSize*m*epsilon / 8.f); // W*M/k
    totalOverflows = new unordered_map<int, int> (maxOverflows);//B TODO: allocate static?

    overflowedArrL0 = new unordered_map<unsigned int, unsigned int>*[blocksNumber];
    int l1Size = ceil(blocksNumber/4); //TODO: define the magic number 4
#ifdef ACC1_DEBUGGING
    cout << "L0 blocks: " << blocksNumber << " L1 blocks: " << l1Size << endl;
#endif

    for(int i =0 ; i < blocksNumber ; i++)
    	overflowedArrL0[i] =  new unordered_map<unsigned int, unsigned int> (maxOverflows);

    overflowedArrL1 = new unordered_map<unsigned int, unsigned int>*[l1Size];

    for(int i = 0 ; i < l1Size ; i++)
    	overflowedArrL1[i] = new unordered_map<unsigned int, unsigned int> (maxOverflows);
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
	++frameItems;
#ifdef ACC1_DEBUGGING
	cout << "*** UPDATE for frame: " << frameItems << "***" << endl;
#endif

	int currBlock = (ceil((double)frameItems / (double)blockSize));

	currBlock = currBlock % (blocksNumber + 1);


	/* New Block */
    if (((frameItems) % blockSize) == 0) {
#ifdef ACC1_DEBUGGING
    	cout << "NEW BLOCK! " << currBlock << endl;
#endif

    	indexTail = (indexTail + 1) % indexSize;
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 0;
        if (currBlock >= 1 && currBlock < blocksNumber) {
        	//TODO

        	if (((currBlock) % 4) != 0)
        		*(overflowedArrL0[currBlock]) = *(overflowedArrL0[currBlock - 1]);

        	if (((currBlock + 1) % 4) == 0 && ((currBlock + 1) / 4) != 1) {
        		unordered_map<unsigned int, unsigned int>* prevL1Table = overflowedArrL1[((currBlock + 1)/ 4) - 2];
        		unordered_map<unsigned int, unsigned int>* mergedblockTables = overflowedArrL1[((currBlock + 1) / 4) - 1];

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

        		*overflowedArrL1[((currBlock + 1) / 4) - 1] = *mergedblockTables;
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

#ifdef ACC1_DEBUGGING
    if (testArr[currBlock - 1]) { //TODO for testing only
    	--testArr[currBlock - 1];
    	if(currBlock == 6) {
    		cout << "block 6 overflow" << endl;
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
#ifdef ACC1_DEBUGGING
        printf("before overflowArrs currBlock: %d blocksNumber:%d itemIdx: %d\n", currBlock, blocksNumber, itemIdx);// TODO: testing
    	//cout << "update level 0 for: " << itemIdx << " from: " << currBlock << " To: " << nextMulOf4(currBlock) << endl;
#endif
        /* Update Level 0 */
/*
    	cout << "Updating level 0, form: " << currBlock << " to: " << nextMulOf4(currBlock) << endl;
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
*/
#ifdef ACC1_DEBUGGING
    	cout << "update level 0 for: " << itemIdx << " index in L0 array: " << currBlock -1 << endl;
#endif

        unordered_map<unsigned int, unsigned int>* blockMapL0 = overflowedArrL0[currBlock -1];

        if(blockMapL0->find(itemIdx) == blockMapL0->end())
        	blockMapL0->insert(pair<int, int> (itemIdx, 1));
        else
        	blockMapL0->at(itemIdx) = blockMapL0->at(itemIdx) + 1;

        overflowedArrL0[currBlock - 1] = blockMapL0;

/*
 *
         Update Level 1


#ifdef ACC1_DEBUGGING
    	cout << "update level 1 for: " << itemIdx << " from: " << floor(log(currBlock)/log(4)) << " To: " << ceil(blocksNumber/4) - 1<< endl;
#endif

    	cout << "Updating level 1, levels number from: " <<  floor(log(currBlock)/log(4)) << " to: " << ceil(blocksNumber/4) << endl;

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
*/

    	int l1Block = (int)ceil(double(currBlock) / 4.0) - 1;
#ifdef ACC1_DEBUGGING
    	cout << "update level 1 for: " << itemIdx << " index in L1 array: " << l1Block << endl;
#endif

        unordered_map<unsigned int, unsigned int>* blockMapL1 = overflowedArrL1[l1Block];

        if(blockMapL1->find(itemIdx) == blockMapL1->end())
        	blockMapL1->insert(pair<int, int> (itemIdx, 1));
        else
        	blockMapL1->at(itemIdx) = blockMapL1->at(itemIdx) + 1;


        overflowedArrL1[l1Block] = blockMapL1;


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
	int secondBlock = (int) floor((double) b2 / (double)blockSize) % (int) (blocksNumber + 1);
	int minOverFlows;
	int itemIdx;
	int tables = 0;

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
#ifdef ACC1_DEBUGGING
			cout << "overflowedArrL1["<< prevL1Block(secondBlock) - 1 <<"] contains: ";

			for ( auto it = overflowedArrL1[prevL1Block(secondBlock) - 1]->begin(); it != overflowedArrL1[prevL1Block(secondBlock) - 1]->end(); ++it )
				cout << " " << it->first << ":" << it->second;
			cout << endl;
#endif


			if (overflowedArrL1[prevL1Block(secondBlock) - 1]->find(itemIdx) != overflowedArrL1[prevL1Block(secondBlock) - 1]->end()) {
				overTillSecond = overflowedArrL1[prevL1Block(secondBlock) - 1]->at(itemIdx);
				//++tables;
			}
		}

		if (4*prevL1Block(secondBlock) != secondBlock) {
#ifdef ACC1_DEBUGGING

			cout << "overflowedArrL0["<< secondBlock - 1 <<"] contains: ";

			for ( auto it = overflowedArrL0[secondBlock - 1]->begin(); it != overflowedArrL0[secondBlock - 1]->end(); ++it )
				cout << " " << it->first << ":" << it->second;
			cout << endl;
#endif


			if (overflowedArrL0[secondBlock - 1]->find(itemIdx) != overflowedArrL0[secondBlock - 1]->end()) {
				overTillSecond += overflowedArrL0[secondBlock - 1]->at(itemIdx);
				//++tables;
			}
		}

		if (prevL1Block(firstBlock)) {
#ifdef ACC1_DEBUGGING

			cout << "overflowedArrL1["<< prevL1Block(firstBlock) - 1 <<"] contains: ";

			for ( auto it = overflowedArrL1[prevL1Block(firstBlock) - 1]->begin(); it != overflowedArrL1[prevL1Block(firstBlock) - 1]->end(); ++it )
				cout << " " << it->first << ":" << it->second;
			cout << endl;
#endif

			if (overflowedArrL1[prevL1Block(firstBlock) - 1]->find(itemIdx) != overflowedArrL1[prevL1Block(firstBlock) - 1]->end()) {
				overTillFirst = overflowedArrL1[prevL1Block(firstBlock) - 1]->at(itemIdx);
				//++tables;
			}
		}

		if (4*prevL1Block(firstBlock)!= firstBlock) {
#ifdef ACC1_DEBUGGING

			cout << "overflowedArrL0["<< firstBlock - 1 <<"] contains: ";

			for ( auto it = overflowedArrL0[firstBlock - 1]->begin(); it != overflowedArrL0[firstBlock - 1]->end(); ++it )
				cout << " " << it->first << ":" << it->second;
			cout << endl;
#endif


			if (overflowedArrL0[firstBlock - 1]->find(itemIdx) != overflowedArrL0[firstBlock - 1]->end()) {
				overTillFirst += overflowedArrL0[firstBlock - 1]->at(itemIdx);
				//++tables;
			}
		}

#ifdef ACC1_DEBUGGING
		cout << "overTillFirst: " << overTillFirst << " overTillSecond: " << overTillSecond << endl;
#endif
		minOverFlows = overTillSecond - overTillFirst;
	}
	// printf("minoverflowed: %d \n", minOverFlows); //TODO: testing
	//cout << "number of tables: " << tables << endl;
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
