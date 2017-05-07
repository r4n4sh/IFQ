//
//  BaseWRSS.cpp
//
//
//  Created by Rana Shahout on 3/15/17.
//
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "BaseWRSS.hpp"

using namespace std;

BaseWRSS::BaseWRSS(unsigned int windowSize, float gamma, unsigned int m, float epsilon)
{
    frameItems = 0;
    overflowsNumber = 0;
    tail = 0;
    indexTail = 0;
    this->blockSize = ceil((windowSize * epsilon)/4); // W/k; k= 4/epsilon
    this->windowSize = windowSize;
    unsigned int blocksNumber = windowSize / blockSize;
    this->m = m;
    this->epsilon = epsilon;

    maxOverflows = min(blocksNumber * 2, windowSize);
    indexSize = maxOverflows + blocksNumber;
    head = maxOverflows - 1;
    indexHead = blocksNumber - 1;
    index = new vector<int> (indexSize); // 0 means end of block
    overflowsElements = new unsigned int[maxOverflows];
    rss = new RSS_CPP(epsilon, m, gamma);//y
    threshold = windowSize*m*epsilon / 4.f; // W*M/k
    totalOverflows = new unordered_map<unsigned int, unsigned int> (maxOverflows);//B TODO: allocate static?
}

BaseWRSS::~BaseWRSS()
{
    delete(totalOverflows);
    delete(rss);
    delete [] overflowsElements;
    delete(index);
}

void BaseWRSS::update(unsigned int item, int wieght)
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
//            --overflowsNumber;
            indexTail = (indexTail + 1) % indexSize;
        }
    } catch (const out_of_range) {
    }

    // Add item to RSS_CPP
    int prevQuery = this->rss->query(item);
    this->rss->update(item, wieght);

    // overflow
    if ((prevQuery%threshold) + wieght > threshold) {
        head = (head + 1) % maxOverflows;
        overflowsElements[head] = item;
  //      ++overflowsNumber;
  //      assert(overflowsNumber < maxOverflows);
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 1;
        if (totalOverflows->find(item) == totalOverflows->end())
            totalOverflows->insert(make_pair<unsigned int,unsigned int>((unsigned int)item,1));
        else
            totalOverflows->at(item) = totalOverflows->at(item) + 1;
    }

    // New frame
    if (frameItems == windowSize) {
        frameItems = 0;
        rss->clear();
    }
}

double BaseWRSS::query(unsigned int item)
{
    int minOverFlows;
    double rssEstimation = this->rss->query(item);

    unordered_map<unsigned int,unsigned int>::const_iterator foundedItem = totalOverflows->find(item);
    if (foundedItem == totalOverflows->end()) // item has no oveflows
        minOverFlows = 0;
    else {
        minOverFlows = totalOverflows->at(item);
    }

    rssEstimation = (int) rssEstimation % (int) this->threshold; //TODO
    return (this->threshold * (minOverFlows + 2 ) + rssEstimation);
}
