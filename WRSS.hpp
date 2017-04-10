//
//  WRSS.hpp
//  
//
//  Created by Rana Shahout on 3/15/17.
//
//

#ifndef WRSS_hpp
#define WRSS_hpp

#include <stdio.h>
#include <unordered_map>
#include <vector>

#include "RSS_CPP.hpp"

class WRSS {
public:
    WRSS(int windowSize, double gamma, int m, double epsilon);
    ~WRSS();
    void update(unsigned int item, int wieght);
    double query(unsigned int item);

private:
    int frameItems;
    int blockSize;
    int windowSize;
    int blocksNumber;
    std::vector<int> *index;
    int tail;
    int overflowsNumber;
    int indexTail;
    int head;
    int indexHead;
    std::unordered_map<int, int> *totalOverflows; //B
    unsigned int *overflowsElements;// b
    RSS_CPP *rss;
    int indexSize;
    int maxOverflows;
    double epsilon;
    double threshold;
    int m; //maximal value of an element in the stream
    double alpha = 0.2;
    double gamma;
    double computeOverflowCount(unsigned int item);
};

#endif /* WRSS_hpp */
