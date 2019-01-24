#include "cuttingGlass.h"

extern Data* DATA;

void swapTwoElementsInList(int p, int q, std::vector<int>& list)
{
	int c = list[p];
	list[p] = list[q];
	list[q] = c;
}

void copyFirstVectorToAnother(std::vector<int>& v1, std::vector<int>& v2)
{
	loop(i, v1.size()) v2[i] = v1[i];
}


void generateRandomFeasibleSequence(std::vector<int>& SBS)
	{
	int tstack[MAX_ITEMS];
	std::vector<int> nmb_selected_items_from_stack(DATA->nb_stacks_, 0);
	int nbs,s,i = 0;
	do {
		nbs=0;
		for (s=0;s<DATA->nb_stacks_;s++)
			if (nmb_selected_items_from_stack[s] < DATA->nmb_items_in_stack[s])
				tstack[nbs++]=s;			
		if (nbs) 
			{
			s = tstack[rand() % nbs];
			SBS[i++] = DATA->stack[s][nmb_selected_items_from_stack[s]];
			nmb_selected_items_from_stack[s]++;
			}
		}
	while (nbs);
	}



void generateRandomFeasibleSequenceFromPosition(Solution* sol, std::vector<int>& SBS, int index_start)
{
	if(index_start < 0 || index_start >= DATA->nb_items_) return;

	std::vector<std::vector<int> > items(DATA->nb_stacks_);
	for(int i = index_start; i < DATA->nb_items_; i++)
	{
		int s = DATA->item_list_[SBS[i]].stack;
		items[s].push_back(DATA->item_list_[SBS[i]].id);
	}

	int loop = 0;
	int i = index_start;
	do {

		int s = sol->rand() % items.size();
		while (items[s].size() == 0)
		{
			s = sol->rand() % items.size();
			loop++;
		}

		SBS[i++] = items[s][0];
		items[s].erase(items[s].begin());
        if(items[s].size() == 0) items.erase(items.begin() + s);


	}
	while (i < DATA->nb_items_);
}

void generateRandomFeasibleSequenceFromPosition(std::vector<int>& SBS, int index_start)
{
	if(index_start < 0 || index_start >= DATA->nb_items_) return;

	std::vector<std::vector<int> > items(DATA->nb_stacks_);
	for(int i = index_start; i < DATA->nb_items_; i++)
	{
		int s = DATA->item_list_[SBS[i]].stack;
		items[s].push_back(DATA->item_list_[SBS[i]].id);
	}

	int loop = 0;
	int i = index_start;
	do {

		int s = rand() % items.size();
		while (items[s].size() == 0)
		{
			s = rand() % items.size();
			loop++;
		}

		SBS[i++] = items[s][0];
		items[s].erase(items[s].begin());
        if(items[s].size() == 0) items.erase(items.begin() + s);


	}
	while (i < DATA->nb_items_);
}

