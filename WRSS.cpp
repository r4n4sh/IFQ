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
    maxOverflows = std::min(blocksNumber * 2, windowSize);
    indexSize = maxOverflows + blocksNumber;
    head = maxOverflows - 1;
    indexHead = blocksNumber - 1;
    index = new std::vector<int> (indexSize);
    overflowsElements = new unsigned int[maxOverflows]; //b
    rss = new RSS_CPP(epsilon, m, gamma);//y
    threshold = windowSize*m*epsilon / 4.; // W*M/k
    totalOverflows = new std::unordered_map<int, int> (maxOverflows);//B
}

WRSS::~WRSS()
{
    free(totalOverflows);
    free(rss);
    free(overflowsElements);
    free(index);
}

double WRSS::computeOverflowCount(unsigned int item)
{
    double k = 4. / this->epsilon;
    return (rss->query(item) / (this->m * this->windowSize / k));
}

void WRSS::update(unsigned int item, int wieght)
{
    printf("WRSS update \n");
    // End of a block
    if (((++frameItems) % blockSize) == 0)
        indexTail = (indexTail + 1) % indexSize;
    
    // New block is started
    if (!(frameItems % blockSize)) {
        indexHead = (indexHead + 1) % indexSize;
        index->at(indexHead) = 0;
    }

    // Remove oldest element in oldest block 

    try {
        if (index->at(indexTail)) {
            printf("inside index->at[indexTail] \n");

            int oldId = overflowsElements[tail];
            if (!(totalOverflows->at(oldId) - 1))
                totalOverflows->erase(oldId);
            else
                totalOverflows->insert(std::make_pair(oldId,totalOverflows->at(oldId) - 1));
            tail = (tail + 1) % maxOverflows;
            --overflowsNumber;
            indexTail = (indexTail + 1) % indexSize;
        }
    } catch (const std::out_of_range) {
    }

    // Add item to RSS_CPP
    int prevOverflowCount = computeOverflowCount(item);
    this->rss->update(item, wieght);
    int currOverflowCount = computeOverflowCount(item);

    // overflow
    if (currOverflowCount > prevOverflowCount) {
        head = (head + 1) % maxOverflows;
        ++overflowsNumber;
        overflowsElements[head] = item;
        indexHead = (indexHead + 1) % indexSize;
        index->at(head) = 1;
        std::unordered_map<int,int>::const_iterator foundedItem = totalOverflows->find(item);
        if (foundedItem == totalOverflows->end())
            totalOverflows->insert(std::make_pair<int,int>(item,1));
        else
            totalOverflows->at(item) = totalOverflows->at(item) + 1;
    }

    // New frame
    if (frameItems == windowSize) {
        frameItems = 0;
        free(rss);
        rss = new RSS_CPP(epsilon, m, this->gamma);
    }
}

double WRSS::query(unsigned int item)
{
    int minOverFlows;
    double rssEstimation = this->rss->query(item);

    std::unordered_map<int,int>::const_iterator foundedItem = totalOverflows->find(item);
    if (foundedItem == totalOverflows->end()) // item has no oveflows
        minOverFlows = 0;
    else {
        minOverFlows = totalOverflows->at(item);
        rssEstimation = (int) rssEstimation % (int) this->threshold;
    }

    return (this->threshold * (minOverFlows + 2 ) + rssEstimation);
}
