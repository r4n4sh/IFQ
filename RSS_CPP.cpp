#include "RSS_CPP.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>

#include "prng.hpp"
/**********************************************************************************************
Implementation of Rounded Space Saving algorithm to Find Frequent Items on weighted streams
Based on the paper of Ben Basat, Einziger and Friedman, 2016
Implementation by Ben Basat and Einziger, 2016

Original Code: 2016-06
This version: 2016-06

This work is licensed under the Creative Commons
Attribution-NonCommercial License. To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc/1.0/ or send a letter
to Creative Commons, 559 Nathan Abbott Way, Stanford, California
94305, USA.
***********************************************************************************************/

RSS_CPP::RSS_CPP(float fPhi, int M, float gamma)
{
	int i;
	//int k = 1 + (int) 1.0 / fPhi;
	m_nrCounters = ceil((int) 1.0 / fPhi);

	//RSS_type* result = (RSS_type*)calloc(1, sizeof(RSS_type));
	//std::cout << "M = " << M << std::endl;
	m_M = M;
	m_gamma = gamma;
	m_stepSize =  1 + ceil(M/ 2.* gamma) ;
	m_stepSizeMinusOne = m_stepSize - 1;
	if (m_stepSize <= 0)
		//throw std::exception::exception("bad config, make sure gamma >= 0");
	assert(m_stepSize);
	assert(gamma > 0);
	m_nrCounters *= (1 + gamma);

	m_a = (long long)698124007;
	m_b = (long long)5125833;
	assert(m_nrCounters > 0);
	//m_n = 0;
#if UNDER_ESTIMATOR
	m_maxRootVictim = 0;
#else
	m_maxRootVictim = m_stepSizeMinusOne;
#endif

	m_tblsz = RSS_HASHMULT*m_nrCounters;
	m_hashtable = (RSSITEM**)calloc(m_tblsz, sizeof(RSSITEM *));
	m_groups = (RSSGROUP *)calloc(m_nrCounters, sizeof(RSSGROUP));
	m_items = (RSSITEM *)calloc(m_nrCounters, sizeof(RSSITEM));
	m_freegroups = (RSSGROUP **)calloc(m_nrCounters, sizeof(RSSGROUP*));

	for (i = 0; i<m_tblsz; i++)
		m_hashtable[i] = NULL;

	m_root = m_groups;
	m_groups->count = 0;
	m_groups->nextg = NULL;
	m_groups->previousg = NULL;

	m_groups->items = m_items;
	for (i = 0; i<m_nrCounters; i++)
		m_freegroups[i] = &m_groups[i];
	m_gpt = 1; // initialize list of free groups

	for (i = 0; i<m_nrCounters; i++)
	{
		m_items[i].item = 0;
		m_items[i].delta = -1;
		m_items[i].hash = 0;
		m_items[i].nexti = NULL;
		m_items[i].previousi = NULL;  // initialize values

		m_items[i].parentg = m_groups;
		m_items[i].nexting = &(m_items[i + 1]);
		m_items[i].previousing = &(m_items[i - 1]);
		// create doubly linked list
	}
	m_items[0].previousing = &(m_items[m_nrCounters - 1]);
	m_items[m_nrCounters - 1].nexting = &(m_items[0]);
	// fix start and end of linked list

	m_nrFree = m_nrCounters;
}

RSS_CPP::~RSS_CPP()
{
	free(m_items);
	free(m_groups);
	free(m_freegroups);
	free(m_hashtable);
}



void RSS_CPP::update(int newitem, int weight) {
	int h = getHash(newitem);
	RSSITEM *til = getCounter(newitem, h);
	RSSITEM &il = (til ? *til : takoverMinimal(newitem, h));
#if UNDER_ESTIMATOR
	if (!(til || m_nrFree))
			il.delta = m_root->count;
	//il.parentg->items = il.parentg->items->nexting;

#endif
	int newVal = weight + (til ? il.remainder : m_maxRootVictim);
	unsigned int nrIncs = newVal / m_stepSize;

	il.remainder = newVal % m_stepSize; //TODO: consider changing % to shifts
	if (nrIncs && (!m_nrFree || til)) {
		advanceCounter(il, nrIncs);
		return;
	}
	if (m_nrFree && !til)
	{ // We deliberately use duplicated code as this `else' is rarely accessed
		newVal = weight;
		nrIncs = newVal / m_stepSize;
		--m_nrFree;
		il.parentg->items = il.parentg->items->nexting;
		il.remainder = newVal % m_stepSize; //TODO: consider changing % to shifts
		if (nrIncs)
			advanceCounter(il, nrIncs);
	}
}

inline int const RSS_CPP::getHash(unsigned int item) {
	return hash31(m_a, m_b, item) % m_tblsz;
}

inline RSSITEM * RSS_CPP::getCounter(unsigned int item, int h) {
	RSSITEM * il = m_hashtable[h];
	while (il && (il->item != item))
		il = il->nexti;
	return il;
}

inline RSSITEM & RSS_CPP::takoverMinimal(unsigned int item, int h) {
	RSSITEM & il = *(m_root->items);
	removeFromHash(il);
	insertIntoHashtable(il, h, item);
#if UNDER_ESTIMATOR
		if (m_maxRootVictim < il.remainder)
		{
			m_maxRootVictim = il.remainder;
		}
#endif
	return il;
}

void RSS_CPP::advanceCounter(RSSITEM &il, unsigned int nrIncs) {
	RSSGROUP &oldgroup = *(il.parentg);
	unsigned int goalVal = oldgroup.count + nrIncs;
	RSSGROUP &group = getLastGroup(oldgroup, goalVal);

	if ((il.nexting) == &il) {      // if the counter is the only one in its group

		if (&group == &oldgroup) {    // if we can simply increase the oldgroup's count
			oldgroup.count = goalVal;

			return;
		}
		if (group.count == goalVal) { // if there exists a group with count = goalVal
			putInNewGroup(il, group);
			return;
		}
		moveGroup(il, group, goalVal);
		return;
	}
	disconnectCounter(il);
	assert(group.count <= goalVal);

	if (group.count == goalVal) {
		putInNewGroup(il, group); // if the goal group exists
		return;
	}
	newGroup(il, group, goalVal);
}

inline void RSS_CPP::disconnectCounter(RSSITEM & il) { //disconnect only if there exists other counters in the same group
	if (il.parentg->items == &il) {
		il.parentg->items = il.nexting;
	}
	//il.parentg->items = (il.parentg->items == &il) ? il.nexting: il.parentg->items);

	assert(il.previousing);
	assert(il.nexting);
	il.previousing->nexting = il.nexting;
	il.nexting->previousing = il.previousing;
	il.nexting = &il;
	il.previousing = &il;
}

inline RSSGROUP & RSS_CPP::getLastGroup(RSSGROUP & oldgroup, unsigned int goalVal) {
	RSSGROUP *group = &oldgroup;
	while (group->nextg && (group->nextg->count <= goalVal)) {
		group = group->nextg;
	}
	assert((group->items->nexting == group->items) == (group->items->previousing == group->items));
	return *group;
}

inline void RSS_CPP::removeFromHash(RSSITEM & il) {
	// if il was first - the chain will point on the next item.
	if (m_hashtable[il.hash] == &il)
		m_hashtable[il.hash] = il.nexti;
	// if there is another node with the same hash.
	if (il.nexti)
	{
		// next item - previous - points to my previous.
		il.nexti->previousi = il.previousi;

	}
	if (il.previousi) {
		//prev item next - points to my next.
		il.previousi->nexti = il.nexti;
	}
}

inline void RSS_CPP::addToGroup(RSSGROUP & group, RSSITEM & il) {
	il.previousing = group.items->previousing;
	il.nexting = group.items;
	il.parentg = &group;
	group.items->previousing->nexti = &il;
	group.items->previousing = &il;
}



inline void const RSS_CPP::showGroups() {
	RSSGROUP *g;
	RSSITEM *i, *first;
	int n, wt;

	g = m_groups;
	wt = 0;
	n = 0;
	while (g != NULL)
	{
		printf("Group %lld :", g->count);
		first = g->items;
		i = first;
		if (i != NULL)
			do
			{
				printf("%d -> ", i->item);
				i = i->nexting;
				wt += g->count;
				n++;
			} while (i != first);
		else printf(" empty");
		printf(")");
		g = g->nextg;
		if ((g != NULL) && (g->previousg->nextg != g))
			printf("Badly linked");
		printf("\n");
	}
	printf("In total, %d items, with a total count of %d\n", n, wt);
}

void RSS_CPP::insertIntoHashtable(RSSITEM &newi, int i, unsigned int newitem) {
	newi.nexti = m_hashtable[i];
	newi.item = newitem; // overwrite the old item
	newi.hash = i;
	newi.previousi = NULL;
	// insert item into the hashtable
	if (m_hashtable[i])
		m_hashtable[i]->previousi = &newi;
	m_hashtable[i] = &newi;
}

unsigned int const RSS_CPP::query(unsigned int x) {
	int h;
	RSSITEM *il;
	h = hash31(m_a, m_b, x) % m_tblsz;
	il = m_hashtable[h];
	// if x has a counter.
	while (il && il->item != x)
		il = il->nexti;
	if (il){
		if (il->delta == -1)
			return il->parentg->count * m_stepSize + il->remainder;
		return il->parentg->count * m_stepSize + il->remainder -(il->delta + 1)*m_stepSize - 1;
	}
	if (UNDER_ESTIMATOR)
		return 0;
	int minCount = m_root->count;
	// if all counters are used.
	if (minCount)
		return minCount * m_stepSize - 1;
	// if there are unused counters and x does not have a counter, we are certain that x never arrived.
	return 0;
}


inline void RSS_CPP::connectToGroup(RSSITEM &newi, RSSGROUP &tmpg) {
	newi.nexting = tmpg.items;
	newi.previousing = tmpg.items->previousing;
	assert((newi.previousing != &newi) && (newi.nexting != &newi));
	newi.previousing->nexting = &newi;
	newi.nexting->previousing = &newi;
}

inline void RSS_CPP::putInNewGroup(RSSITEM &newi, RSSGROUP & tmpg) {
	RSSGROUP * oldgroup = newi.parentg;
	assert(oldgroup != &tmpg);
	// put item in the tmpg group
	newi.parentg = &tmpg;
	assert((newi.previousing == &newi) && (newi.nexting == &newi));
	if (oldgroup->items == &newi) { // group will be empty
		recycleGroup(*oldgroup);
	}
	connectToGroup(newi, tmpg);
}

inline void RSS_CPP::moveGroup(RSSITEM &il, RSSGROUP &group, unsigned int goalVal) {
	assert(il.parentg->nextg && il.parentg->nextg->count < goalVal);
	RSSGROUP &tmpg = *(il.parentg);

	if (&tmpg == m_root) // might fail if we allocate only a single counter
	{
		m_root = tmpg.nextg;
		m_root->previousg = NULL;
#if UNDER_ESTIMATOR
		m_maxRootVictim = 0;
#endif
	}
	assert(!(m_root->previousg));
	if (tmpg.previousg)
		tmpg.previousg->nextg = tmpg.nextg;
	assert(tmpg.nextg);
	tmpg.nextg->previousg = tmpg.previousg;
	if (group.nextg) {
		group.nextg->previousg = &tmpg;
	}
	tmpg.nextg = group.nextg;
	tmpg.previousg = &group;
	tmpg.count = goalVal;
	group.nextg = &tmpg;
	assert(m_root->nextg != m_root);
}

inline void RSS_CPP::newGroup(RSSITEM &il, RSSGROUP &group, unsigned int goalVal) {
	assert(group.count < goalVal);
	//std::cout << m_gpt << ',' << m_n << std::endl;
	assert(il.parentg->items != &il);
	// Create a new group for newi
	// RSS_AddNewGroupAfter(rss, newi,group);
	RSSGROUP &newgroup = *(m_freegroups[m_gpt++]); //get new group
	newgroup.count = goalVal; // set count to the actual value of the new group
	newgroup.items = &il;
	//newgroup->previousg = group;
	if (group.nextg) { // if there is another group
		group.nextg->previousg = &newgroup;
	}
	newgroup.nextg = group.nextg;
	newgroup.previousg = &group;
	group.nextg = &newgroup;
	il.parentg = &newgroup;
	assert(m_root->nextg != m_root);
}

int const RSS_CPP::size() {
	return 0;
	//return sizeof(RSS_CPP) + (m_tblsz) * sizeof(RSSITEM*) +
		//(m_nrCounters)*(sizeof(RSSITEM) + sizeof(RSSGROUP) + sizeof(RSSITEM*));
}

inline void RSS_CPP::recycleGroup(RSSGROUP & oldgroup)
{
	if (oldgroup.nextg) // there is another group
		oldgroup.nextg->previousg = oldgroup.previousg;
	if (oldgroup.previousg)
		oldgroup.previousg->nextg = oldgroup.nextg;
	if (m_root == &oldgroup) // this is the first group
	{
		assert(!(oldgroup.nextg->previousg));
		m_root = oldgroup.nextg;

	}

	//recycle oldgroup.
	m_freegroups[--m_gpt] = &oldgroup;
	// if we have created an empty group, remove it
	assert(m_root->nextg != m_root);
}

void RSS_CPP::clear()
{
	int i;
		/*for (i = 0; i<m_tblsz; i++)
			m_hashtable[i] = NULL;

		m_groups->count = 0;
		m_groups->nextg = NULL;
		m_groups->previousg = NULL;

		m_groups->items = m_items;// TODO: check if there is more counters to reset
		m_root = m_groups;
		m_nrFree = m_nrCounters;*/
		for (i = 0; i<m_tblsz; i++)
			m_hashtable[i] = NULL;

		m_root = m_groups;
		m_groups->count = 0;
		m_groups->nextg = NULL;
		m_groups->previousg = NULL;

		m_groups->items = m_items;
		for (i = 0; i<m_nrCounters; i++)
			m_freegroups[i] = &m_groups[i];
		m_gpt = 1; // initialize list of free groups

		for (i = 0; i<m_nrCounters; i++)
		{
			m_items[i].item = 0;
			//m_items[i].delta = -1;
			m_items[i].hash = 0;
			m_items[i].nexti = NULL;
			m_items[i].previousi = NULL;  // initialize values

			m_items[i].parentg = m_groups;
			m_items[i].nexting = &(m_items[i + 1]);
			m_items[i].previousing = &(m_items[i - 1]);
			// create doubly linked list
		}
		m_items[0].previousing = &(m_items[m_nrCounters - 1]);
		m_items[m_nrCounters - 1].nexting = &(m_items[0]);
		// fix start and end of linked list

		m_nrFree = m_nrCounters;
}
