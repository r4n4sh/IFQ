/*
 * RAW.cpp
 *
 *  Created on: May 8, 2017
 *      Author: ranashahout
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "RAW.hpp"

RAW::RAW(unsigned int windowSize, float gamma, unsigned int m, float epsilon)
{
    frameItems = 0;
    this->blockSize = ceil((windowSize * epsilon)/4); // W/k; k= 4/epsilon
    this->windowSize = windowSize;
    this->blocksNumber = windowSize / blockSize;
    this->m = m;
    this->epsilon = epsilon;
    maxOverflows = min(blocksNumber * 2, windowSize);
    threshold = ceil(windowSize*m*epsilon / 4.f); // W*M/k

    windows = new BaseWRSS*[blocksNumber];
    for (int i = 0; i < blocksNumber; i++)
    	windows[i] = new BaseWRSS(windowSize - (i*blockSize), gamma, m, epsilon);
}

RAW::~RAW()
{
    for (int i = 0; i < blocksNumber; i++)
    	delete(windows[i]);

    delete[](windows);
}

void RAW::update(unsigned int item, int wieght)
{
    for (int i = 0; i < blocksNumber; i++)
    	windows[i]->update(item, wieght);

}
double RAW::query(unsigned int item)
{
	return windows[0]->query(item);
}
double RAW::intervalQuery(unsigned int item, int b1, int b2)
{
	int firstBlock = (int) floor(b1 / blockSize) % (int) blocksNumber;
	int secondBlock = (int) floor(b2 / blockSize) % (int) blocksNumber;
	int windowidx1 = ceil((windowSize - firstBlock)/ blockSize);
	int windowidx2 = ceil((windowSize - secondBlock)/ blockSize);

	if (windowidx2 > windowidx1)
		return windows[windowidx2]->query(item) - windows[windowidx1]->query(item);
	else
		return windows[windowidx1]->query(item) - windows[windowidx2]->query(item);
}


