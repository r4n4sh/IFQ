#pragma once
#pragma once
#ifndef RSS_CPP_h
#define RSS_CPP_h
#include "lossycount.hpp"
#define RSS_HASHMULT 3
#define UNDER_ESTIMATOR 0

typedef struct rss_item RSSITEM;
typedef struct rss_group RSSGROUP;

typedef long long LCUWT;

struct rss_group
{
	LCUWT count;
	RSSITEM *items;
	RSSGROUP *previousg, *nextg;
};

struct rss_item
{
	unsigned int item;
	int hash;
	LCUWT delta;
	int remainder;
	RSSGROUP *parentg;
	RSSITEM *previousi, *nexti;
	RSSITEM *nexting, *previousing;
};

class RSS_CPP {
public:
	RSS_CPP(float fPhi, int M, float gamma);
	~RSS_CPP();
	void update(int item, int weight);
	int const size();
	void clear();
	unsigned int const query(unsigned int item);
private:
	void recycleGroup(RSSGROUP& oldgroup);
	inline void connectToGroup(RSSITEM &newi, RSSGROUP &tmpg);
	void const showGroups();
	void insertIntoHashtable(RSSITEM &newi, int i, unsigned int newitem);
	RSSITEM & takoverMinimal(unsigned int item, int h);
	RSSITEM * getCounter(unsigned int item, int h);
	RSSITEM * getNewCounter(unsigned int nrIncs);
	inline int const getHash(unsigned int item);
	void advanceCounter(RSSITEM &il, unsigned int nrIncs);
	void removeFromHash(RSSITEM & il);
	RSSGROUP & getLastGroup(RSSGROUP &oldgroup, unsigned int goalVal);
	void putInNewGroup(RSSITEM &newi, RSSGROUP & tmpg);
	void moveGroup(RSSITEM &il, RSSGROUP &group, unsigned int goalVal);
	void addToGroup(RSSGROUP & group, RSSITEM & il);
	void newGroup(RSSITEM &il, RSSGROUP &group, unsigned int goalVal);
	void disconnectCounter(RSSITEM & il);
	void computeNrIncsAndAdvance(int newVal, RSSITEM &il);
	//LCUWT m_n;
	int m_gpt;
	int m_tblsz;
	int m_maxRootVictim;
	long long m_a, m_b;
	RSSGROUP * m_root;
	int m_stepSize;
	int m_stepSizeMinusOne;
	int m_M;
	float m_gamma;
	int m_nrCounters;
	unsigned int m_nrFree;
#ifdef LCU_SIZE
	LCUITEM items[LCU_SIZE];
	LCUGROUP groups[LCU_SIZE];
	LCUGROUP *freegroups[LCU_SIZE];
	LCUITEM* hashtable[LCU_TBLSIZE];
#else
	RSSITEM  *m_items;
	RSSGROUP *m_groups;
	RSSGROUP **m_freegroups;
	RSSITEM  **m_hashtable;

#endif
};

#endif
