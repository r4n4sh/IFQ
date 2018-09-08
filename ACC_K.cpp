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
    lastBlock = 1;
    index = new vector<int> (indexSize); // 0 means end of block
    overflowsElements = new unsigned int[maxOverflows];
    rss = new RSS_CPP(epsilon, m, gamma);//y
    totalOverflows = new unordered_map<int, int> (maxOverflows);//B TODO: allocate static?
    overflowedArrLevels = new unordered_map<unsigned int, unsigned int>**[k];
    ghost_tables = new unordered_map<unsigned int, unsigned int>*[k];


    unsigned int level_size = blocksNumber + 1;
    for(int level = 0; level < k; level++) {
        overflowedArrLevels[level] = new unordered_map<unsigned int, unsigned int>* [level_size];
        for(int i = 0 ; i < level_size ; i++)
        	overflowedArrLevels[level][i] = new unordered_map<unsigned int, unsigned int> (maxOverflows);

      //  level_size = ceil(double(level_size / step)) + 1;
    }

    for (int i = 0; i < k; ++i) {
    	ghost_tables[i] = new unordered_map<unsigned int, unsigned int> (maxOverflows);
    }

    incTable = new unordered_map<unsigned int, unsigned int>*[blocksNumber];
    for(int i =0 ; i < blocksNumber ; i++)
        incTable[i] =  new unordered_map<unsigned int, unsigned int> (maxOverflows);

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

void ACC_K::populateIncTable(unsigned int itemIdx)
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

void ACC_K::endBlock(int blockNumber)
{
    indexTail = (indexTail + 1) % indexSize;
    indexHead = (indexHead + 1) % indexSize;
    index->at(indexHead) = 0;

    int l = 0;
    while ((l < k-1) && (lastBlock % (int)(pow(step, l+1)) == 0)) {
            incTable[l]->clear();
            ghost_tables[l]->clear();
            ++l;
    }

    for (auto it = overflowedArrLevels[l][lastBlock]->begin(); it != overflowedArrLevels[l][lastBlock]->end(); ++it)
        ghost_tables[l]->insert(pair<unsigned int, unsigned int>(it->first, it->second));

    for (auto it = incTable[l]->begin(); it != incTable[l]->end(); ++it)
        overflowedArrLevels[l][lastBlock]->insert(pair<unsigned int, unsigned int>(it->first, it->second));

    if ((lastBlock % blocksNumber) == 0) {
      incTable[k-1]->clear();
      ghost_tables[k-1]->clear();
    }

    lastBlock = 1 + (lastBlock % (int) (blocksNumber + 1));
}

void ACC_K::update(unsigned int item, int wieght)
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
   if ((this->rss->query(item) %blockSize) == 0) {
    	head = (head + 1) % maxOverflows;
        overflowsElements[head] = item;
        if (idToIDx.find(item) == idToIDx.end()) {
        	idToIDx.insert(pair<unsigned int, unsigned int> (item, idx));
        	++idx;
        }
   	    int itemIdx = idToIDx.at(item);
        ++overflowsNumber;
        //assert(overflowsNumber <= maxOverflows);
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 1;
        unordered_map<int,int>::const_iterator foundedItem = totalOverflows->find(item);
        if (foundedItem == totalOverflows->end())
            totalOverflows->insert(make_pair<int,int>(item,1));
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
        //rssEstimation = (int) rssEstimation % (int) this->blockSize;//TODO
    }

    rssEstimation = (int) rssEstimation % (int) this->blockSize; //TODO
    return (this->blockSize * (minOverFlows + 2 ) + rssEstimation);
	// return intervalQuery(lastOverflowed, 2, 9); //TODO: testing
}


double ACC_K::winQuery(unsigned int itemIdx, int w) {

  int cFreq = 0;
  unordered_map<unsigned int, unsigned int>* ctable = incTable[0];
//  unordered_map<unsigned int, unsigned int>* ctable = overflowedArrLevels[0][lastBlock];

  if (ctable->find(itemIdx) != ctable->end()) {
    cFreq = incTable[0]->at(itemIdx);
  }

  int lastLevel = min((int) (log(lastBlock) / log(step)), (int) (k - 1));
  int cblock;

  for (int l = 1; l <= lastLevel; ++l) {
    cblock = pow(step, l) * floor(lastBlock / (pow(step, l)));

    unordered_map<unsigned int, unsigned int>* table = overflowedArrLevels[l][cblock];
    if (table->find(itemIdx) != table->end()) {
      cFreq += table->at(itemIdx);
    }
  }

  //if (w <= lastBlock + 1) {
  if (w <= lastBlock) {

    int tmp_res = 0;

    //lastLevel = (int)(log(lastBlock + 1 - w) / log(step));
    lastLevel = min(((int)(log(lastBlock - w) / log(step))), (int)(k - 1));

    for (int l = 0; l <= lastLevel; ++l) {
      //cblock = pow(step, l) * floor((lastBlock + 1 -w) / (pow(step, l)));
      cblock = pow(step, l) * floor((lastBlock -w) / (pow(step, l)));

      unordered_map<unsigned int, unsigned int>* table = overflowedArrLevels[l][cblock];
      if (table->find(itemIdx) != table->end()) {
        tmp_res += table->at(itemIdx);
      }
    }

    return abs(cFreq - tmp_res);
  }
  else {
    int b = blocksNumber + lastBlock + 1 - w;

    int last_l = 0;
    while (pow(step, last_l) * floor(lastBlock + 1 -w / (pow(step, last_l))))
      ++last_l;
    
    int pre_w = 0;
    
    for (int l = 0; l < last_l; ++l) {
      cblock = pow(step, l) * floor(b/ (pow(step, l)));
      unordered_map<unsigned int, unsigned int>* table = overflowedArrLevels[l][cblock];
      if (table->find(itemIdx) != table->end()) {
        pre_w += table->at(itemIdx);
      }
    }

    for (int l = last_l + 1; l < k - 1; ++l) {
      pre_w += ghost_tables[l]->at(itemIdx);
    }


    unordered_map<unsigned int, unsigned int>* table = overflowedArrLevels[k-1][blocksNumber];
    if (table->find(itemIdx) != table->end()) {
      cFreq += table->at(itemIdx);
    }

    return cFreq - pre_w;
  }
}

double ACC_K::intervalQuery(unsigned int item, int i, int j)
{
	int sum = 0;
	int itemIdx;

	unordered_map<unsigned int,unsigned int>::const_iterator foundedItem = idToIDx.find(item);
	if (foundedItem == idToIDx.end()) // item has no overflows
		sum = 0;
	else {
		itemIdx = idToIDx.at(item);// founded->second
    if (i == 0)
      sum = winQuery(itemIdx, j);
    else
      sum = abs((int)(winQuery(itemIdx, j) - winQuery(itemIdx, i)));
	}

	//return blockSize * (sum + 2);
  return blockSize * (sum );

}

double ACC_K::intervalFrequencyQuery(unsigned int item, int i, int j)
{
    return blockSize *(intervalQuery(item, ceil((double)i/(double)blockSize), floor((double)j/(double)blockSize)) + 2);
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




