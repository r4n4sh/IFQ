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
int testArr[30] = {1,2,2,0,0,1,1,0,1,2,2,0,0,1,1,0,1,2,2,0,0,1,1,0,1,2,2,1,1,1};
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
    this->windowSize = windowSize;
    this->blockSize = ceil((windowSize * epsilon)/6.0); // W/k; k= 6/epsilon
    this->blocksNumber = windowSize / blockSize;
    this->m = m;
    this->epsilon = epsilon;
    this->gamma = gamma;
    this->step = floor(pow(blocksNumber, (1.0/float(k))));
    this->blocksNumber = ceil(blocksNumber / pow(step, k-1)) * (pow(step, k - 1));
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
    ghost_tables = new unordered_map<unsigned int, unsigned int>*[k];


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

    for (int i = 0; i < k; ++i) {
    	ghost_tables[i] = new unordered_map<unsigned int, unsigned int> (maxOverflows);
    }
}

ACC_K::~ACC_K()
{

    for (int i = 0; i < k -1; ++i) {
    	delete (ghost_tables[i]);
    }
    delete[] ghost_tables;
    int level_size = blocksNumber;

    for(int level = 0; level < k; level++) {
        for(int i = 0 ; i < level_size ; i++)
        	delete (overflowedArrLevels[level][i]);
        delete(overflowedArrLevels[level]);
        level_size = ceil(double(level_size / step));
    }
    delete [] overflowedArrLevels;
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
        	for (int level = 1; level < k - 1; ++level) {
        		if (currBlock % (int)pow(step, level + 1) != 0 && (currBlock - (int)pow(step, level)) > 0 && currBlock % (int)pow(step, level) == 0  && (currBlock - (int)pow(step, level)) % (int)(pow(step, level + 1)) != 0 && ((currBlock) / pow(step, level)) > 1) {
                    int cell = (int)((currBlock)/ pow(step, level) ) - 2;
                    if (cell < 0)
                 	   cell = 0;
#ifdef ACC_K_DEBUGGING
        			cout << "update level: " << level << " for currBlock: " << currBlock << endl;
        			cout << "copy level: " << level << " table: " << cell << " to table: " << (int)((currBlock) / pow(step, level) ) - 1 << endl;
#endif

        			unordered_map<unsigned int, unsigned int>* prevL1Table = overflowedArrLevels[level][cell];

                    cell = (int)((currBlock) / pow(step, level) ) - 1;
                    if (cell < 0)
                 	   cell = 0;
        			unordered_map<unsigned int, unsigned int>* mergedblockTables = overflowedArrLevels[level][cell];

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
					*overflowedArrLevels[level][cell] = *mergedblockTables;
        		}
        	}
        	// Always copy tables from previous block in level k -1

       		if((currBlock) % (int)pow(step, k -1) == 0 && ((currBlock) / pow(step, k -1) ) != 1) {

   				int cell = (int)((currBlock)/ pow(step, k -1) ) - 2;
   				if (cell < 0)
   					cell = 0;
#ifdef ACC_K_DEBUGGING
        	cout << "update level: " << k -1 << " for currBlock: " << currBlock << endl;
   			cout << "copy level: " << k -1 << " table: " << cell << " to table: " << (int)((currBlock) / pow(step, k -1) ) - 1 << endl;
#endif

   				unordered_map<unsigned int, unsigned int>* prevL1Table = overflowedArrLevels[k -1][cell];
                cell = (int)((currBlock) / pow(step, k -1) ) - 1;
                if (cell < 0)
             	   cell = 0;

           		unordered_map<unsigned int, unsigned int>* mergedblockTables = overflowedArrLevels[k -1][cell];

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
           		*overflowedArrLevels[k -1][cell] = *mergedblockTables;
       		}
        }


        /*
         * Ghost tables
         */

        for (int level = 1; level < k; level++) {
        	if(currBlock % (int)(pow(step, level)) == 0) {
        		ghost_tables[level] = overflowedArrLevels[level][(int)((currBlock)/ pow(step,level) ) - 1];
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
    if (frameItems == blocksNumber * blockSize) {
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

       while (required_block % (int)pow(step, l_min) == 0 && l_min < k){
               ++l_min;
       }

       l_min -= 1;
       int l = l_min;

       if (l_min == k -1)
               l_max = k -1;
       else
               l_max = ceil(log(required_block) / log(step)) - 1;
       int block = required_block;
#ifdef ACC_K_DEBUGGING
       cout << "Block: " << required_block << " l_min: " << l_min << " l_max: " << l_max << " step: " << step << endl;
#endif
       while (l <= l_max) {
               int cell = (int)floor((block)/pow(step, l)) - 1;
               if (cell < 0)
            	   cell = 0;
#ifdef ACC_K_DEBUGGING
               cout << "lmin: " << l << " l_max: " << l_max << " block: " << block << endl;
               cout << "table [" << l << "]" << "[" <<  cell << "]" << endl;
#endif

               unordered_map<unsigned int, unsigned int>* table = overflowedArrLevels[l][cell];
               if (table->find(itemIdx) != table->end()) {
                       sum += table->at(itemIdx);
               }
               block -= pow(step, l);
               l++;
       }
#ifdef ACC_K_DEBUGGING
       cout << "sum till: " << required_block << " is: " << sum << endl;
#endif
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
		} else{

			sum_till_second = withinFrameFrequency(secondBlock, itemIdx);

            int cell = (int)((blocksNumber)/pow(step, k - 1)) - 1;
            if (cell < 0)
         	   cell = 0;

			unordered_map<unsigned int, unsigned int>* table = overflowedArrLevels[k -1][cell];
			if (table->find(itemIdx) != table->end()) {
				sum_till_second += table->at(itemIdx);
			}

			int offset = ((int)(floor((double)frameItems / (double)blockSize)) % (int)(blocksNumber)) + 1;
			int location = firstBlock - 1;
			int l = 1;

			while (location % (int)pow(step, l) == 0 && l < k){
				++l;
			}

			l -= 1;

			while (l <= k -1) {
				while (location >= offset +1) {
#ifdef ACC_K_DEBUGGING
					cout << "l: " << l << " location: " << location << " offset: " << offset << endl;
#endif
		               int cell = (int)((location)/pow(step, l)) - 1;
		               if (cell < 0)
		            	   cell = 0;

					unordered_map<unsigned int, unsigned int>* table = overflowedArrLevels[l][cell];
					if (table->find(itemIdx) != table->end()) {
						sum_till_first += table->at(itemIdx);
					}
					++l;
					location = floor(location / pow(step, l)) * pow(step, l);
				}
				//read ghost from l till k -1

				while (l <= k -1) {
#ifdef ACC_K_DEBUGGING
					cout << "reading ghost_tables["<<l<<"]" << endl;
#endif
					unordered_map<unsigned int, unsigned int>* table = ghost_tables[l];
					if (table->find(itemIdx) != table->end()) {
						sum_till_first += table->at(itemIdx);
					}
					++l;
				}
			}
		}

		sum = sum_till_second - sum_till_first;
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

    intervalQuery(firstOverflow, 52, 66);
}

#endif




