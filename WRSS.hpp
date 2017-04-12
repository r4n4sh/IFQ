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

using namespace std;

typedef struct Block {
	int blockNumber;
	int blockLevel;

	Block(int n, int l){
		blockNumber = n;
		blockLevel = l;
	}

	bool operator== (const Block &b) const {
		return blockNumber == b.blockNumber && blockLevel == b.blockLevel;
	}

	bool operator< (Block &b) {
		return (blockNumber < b.blockNumber || (blockNumber == b.blockNumber && b.blockLevel < b.blockLevel));
	}
} Block;

namespace std {
template <>
struct hash<Block>
{
  std::size_t operator()(const Block& b) const
  {
    using std::size_t;
    using std::hash;

    // Compute individual hash values for blockNumber,
    // and blockLevel and combine them using XOR
    // and bit shifting:

    return ((hash<int>()(b.blockNumber)
             ^ (hash<int>()(b.blockLevel) << 1)) >> 1);
  }
};
}

class WRSS {
public:
    WRSS(int windowSize, double gamma, int m, double epsilon);
    ~WRSS();
    void update(unsigned int item, int wieght);
    double query(unsigned int item);
    /* Query Interval Section */
    double intervalQuery(unsigned int item, int b1, int b2);

private:
    int frameItems;
    int blockSize;
    int windowSize;
    int blocksNumber;
    vector<int> *index;
    int tail;
    int overflowsNumber;
    int indexTail;
    int head;
    int indexHead;
    unordered_map<int, int> *totalOverflows; //B
    unsigned int *overflowsElements;// b
    RSS_CPP *rss;
    int indexSize;
    int maxOverflows;
    double epsilon;
    double threshold;
    int m; //maximal value of an element in the stream
    //double alpha = 0.2;
    double gamma;
    double computeOverflowCount(unsigned int item);

    /* Query Interval Section */
    void populateSkipList(int blockNumber, int itemIdx);

    int skiplistSize;
    unordered_map<int, int> idToIDx;
    unordered_map<Block, unordered_map<int, int> > skiplistMap;
};

#endif /* WRSS_hpp */
