/*
 * RAW.hpp
 *
 *  Created on: May 8, 2017
 *      Author: ranashahout
 */

#ifndef RAW_HPP_
#define RAW_HPP_

#include "BaseWRSS.hpp"

using namespace std;

class RAW {
public:
	RAW(unsigned int windowSize, float gamma, unsigned int m, float epsilon);
    ~RAW();
    void update(unsigned int item, int wieght);
    double query(unsigned int item);
    double intervalQuery(unsigned int item, int b1, int b2);

private:
    int frameItems;
    int blockSize;
    unsigned int windowSize;
    unsigned int blocksNumber;
    int maxOverflows;
    float epsilon;
    unsigned int threshold;
    int m; //maximal value of an element in the stream
    BaseWRSS** windows;
};

#endif /* RAW_HPP_ */
