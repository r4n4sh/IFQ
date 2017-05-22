//
//  HIT.hpp
//  
//
//  Created by Rana Shahout on 3/15/17.
//
//

#ifndef HIT_hpp
#define HIT_hpp

#include <stdio.h>
#include <unordered_map>
#include <vector>

#include "RSS_CPP.hpp"

//#define DEBUGGING 1

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

class HIT {
public:
    HIT(unsigned int windowSize, float gamma, unsigned int m, float epsilon);
    ~HIT();
    void update(unsigned int item, int wieght);
    double query(unsigned int item);
    /* Query Interval Section */
    double intervalQuery(unsigned int item, int b1, int b2);

#ifdef DEBUGGING
    double testIntervalQuery(unsigned int b2, unsigned int b1);
#endif
private:
    unsigned int frameItems;
    unsigned int blockSize;
    unsigned int windowSize;
    unsigned int blocksNumber;
    vector<int> *index;
    unsigned int tail;
    unsigned int overflowsNumber;
    unsigned int idx;
    unsigned int indexTail;
    unsigned int head;
    unsigned int indexHead;
    unordered_map<unsigned int, unsigned int> *totalOverflows; //B
    unsigned int *overflowsElements;// b
    RSS_CPP *rss;
    unsigned int indexSize;
    unsigned int maxOverflows;
    double epsilon;
    unsigned int threshold;
    int m; //maximal value of an element in the stream

    double computeOverflowCount(unsigned int item);

    /* Query Interval Section */
    void populateSkipListLevel_0(unsigned int blockNumber, unsigned int itemIdx);
    void populateSkipListLevels(unsigned int blockNumber);
    double partialIntervalQuery(unsigned int itemIdx, unsigned int secondBlock, unsigned int firstBlock);

    unsigned int skiplistSize;
    unordered_map<unsigned int, unsigned int> idToIDx;
    unordered_map<Block, unordered_map<unsigned int, unsigned int> > skiplistMap;
};

#endif /* HIT_hpp */
