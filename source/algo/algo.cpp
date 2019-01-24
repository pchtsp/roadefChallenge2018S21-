#include "cuttingGlass.h"
#include <thread>
using namespace std;

#define int_64 long long int

extern Data* DATA;
extern string SOLUTION_FILE_NAME;
extern int EVAL_TO_USE;
extern vector<Solution*> SOL_THREAD;
extern int NMB_THREADS;
extern long int SEED;
extern int EVAL_TO_USE;

namespace evalBatch0
{
	std::vector<int> evalPlate_to_Mark(Solution* sol, int i0);
	std::vector<int> evalPlateFromCut1(Solution* sol, int i0);
	int evalBatch(Solution* sol, bool print = false);
	void initDefectsPlate(int pId);
}

namespace gensol
{
	int evalBatch(Solution* sol);
}

namespace evalBatch2
{
	int evalBatch(Solution* sol, int timeLimit);
	std::vector<int> evalPlate(Solution* sol, int i0);
}

namespace gensol2
{
	int evalBatch(Solution* sol);
}

void swapTwoElementsInList(int p, int q, vector<int>& list);
void copyFirstVectorToAnother(vector<int>& v1, vector<int>& v2);
void generateRandomFeasibleSequenceFromPosition(std::vector<int>& SBS, int index_start);
void generateRandomFeasibleSequenceFromPosition(Solution* sol, std::vector<int>& SBS, int index_start);

class SwapMove {
public:
	int id;
	int i;
	int j;
	int index_start_next_plate;
	int plate_surface_used;
	int plate_last_cut_x;
	int_64 value;
	SwapMove(int idd, int ii, int jj, int ind, int s, int x0) : id(idd), i(ii), j(jj), index_start_next_plate(ind), plate_surface_used(s),
			plate_last_cut_x(x0) {}
};

bool USE_STACKS_IN_FITNESS_FUNCTION = false;
int WEIGHT_STACKS_IN_FITNESS_FUNCTION = 100;

int_64 calcval(int surf, int indexStartNextPlate, int xlast, Solution* sol)
{

	if(USE_STACKS_IN_FITNESS_FUNCTION)
	{
		int valueStacks = 0;
		std::vector<int> nmbRemItemsFromStack(DATA->nb_stacks_, 0);
		for(int i = indexStartNextPlate; i < DATA->nb_items_; i++)
		{
			nmbRemItemsFromStack[DATA->item_list_[sol->seq_[i]].stack]++;
		}
		for(int i = 0; i < DATA->nb_stacks_; i++)
		{
		   valueStacks += (nmbRemItemsFromStack[i] * nmbRemItemsFromStack[i]);
		}

		return  LONG_LONG_MAX / 2 - WEIGHT_STACKS_IN_FITNESS_FUNCTION * valueStacks + (surf - xlast) / 100;
	}

	return (surf * 100 - 10 * indexStartNextPlate - xlast);
}

bool greater_value(SwapMove m1, SwapMove m2)
{
	return (m1.value  > m2.value);
}

bool smallersecond(const pair<int, int>& p1, const pair<int, int>& p2)
{
	return p1.second < p2.second;
}

vector<int> eval_plate(Solution* sol, int index_start)
{
	if(EVAL_TO_USE == 1) return evalBatch0::evalPlate_to_Mark(sol, index_start);
	return evalBatch2::evalPlate(sol, index_start);
}

void randomshuffle(Solution* sol, std::vector<int>& v)
{
	if(v.size() == 0) return;
	for(int i = 0; i < 1000; i++)
	{
		int p = sol->rand() % v.size();
		int q = sol->rand() % v.size();
		int c = v[p];
		v[p] = v[q];
		v[q] = c;
	}
}

vector<SwapMove> evaluate_swap_moves_for_plate(Solution* sol, const int index_start, const int index_end, int timeLimit,
		bool shake = false, int currentPlateSurf = -1)
{

	int startTime = time(0);

	bool allowSwapsInsidePlate = (sol->rand() % 3 == 0);
	if(index_end >= DATA->nb_items_ - 1) allowSwapsInsidePlate = true;

	vector<SwapMove> list_of_moves;

	int id = 0;
	std::vector<int> seqcopy = sol->seq_;

	int maxNmbMovesToEvaluate = 50;
	int nmbEvaluated = 0;

	if(currentPlateSurf < 0)
	{
		vector<int> eval_plate_result = eval_plate(sol, index_start);
		currentPlateSurf = eval_plate_result[0];
	}

	std::vector<int> indices_i;
	loop2(ii, index_start, index_end + 1) indices_i.push_back(ii);

	//random_shuffle(indices_i.begin(), indices_i.end());
	randomshuffle(sol, indices_i);

	loop(iii, indices_i.size())
	{
		if(time(0) - startTime >= timeLimit) break;

		int ii = indices_i[iii];
		int i = sol->seq_[ii];
		int si = DATA->item_list_[i].stack;
		if(nmbEvaluated >= maxNmbMovesToEvaluate) break;

		vector<int> stack_already_scanned(DATA->nb_stacks_, 0);
		std::vector<int> indices_j;

		int startJJ = index_end + 1;
		if(allowSwapsInsidePlate) startJJ = ii + 1;

		loop2(jj, startJJ,  DATA->nb_items_)
		{
			int j = sol->seq_[jj];
			int sj = DATA->item_list_[j].stack;
			if(sj == si) continue;
			if(stack_already_scanned[sj]) continue;
			stack_already_scanned[sj] = 1;
			indices_j.push_back(jj);
		}

		//random_shuffle(indices_j.begin(), indices_j.end());
		randomshuffle(sol, indices_j);

		loop(jjj, indices_j.size())
		{

			if(time(0) - startTime >= timeLimit) break;

			if(nmbEvaluated >= maxNmbMovesToEvaluate) break;

			int jj = indices_j[jjj];
			int j = sol->seq_[jj];
			int sj = DATA->item_list_[j].stack;
			if (sol->tabu[j][ii] < sol->iter && DATA->item_list_[i].w != DATA->item_list_[j].w)
			{

				swapTwoElementsInList(ii, jj, sol->seq_);

				std::vector<int> items_to_reorder_positions_ii;
				std::vector<pair<int,int> > items_to_reorder_id_and_seq_ii;
				std::vector<int> items_to_reorder_positions_jj;
				std::vector<pair<int,int> > items_to_reorder_id_and_seq_jj;
				for(unsigned int m = index_start; m < sol->seq_.size(); m++)
				{
					if(DATA->item_list_[sol->seq_[m]].stack == si) {
					  items_to_reorder_positions_ii.push_back(m);
					  items_to_reorder_id_and_seq_ii.push_back(make_pair(sol->seq_[m], DATA->item_list_[sol->seq_[m]].seq));
					}
					if(DATA->item_list_[sol->seq_[m]].stack == sj) {
					   items_to_reorder_positions_jj.push_back(m);
					   items_to_reorder_id_and_seq_jj.push_back(make_pair(sol->seq_[m], DATA->item_list_[sol->seq_[m]].seq));
					}
					if(items_to_reorder_positions_ii.size() == DATA->stack[si].size() &&
					   items_to_reorder_positions_jj.size() == DATA->stack[sj].size()) break;
				}
				sort(items_to_reorder_id_and_seq_ii.begin(), items_to_reorder_id_and_seq_ii.end(), smallersecond);
				sort(items_to_reorder_id_and_seq_jj.begin(), items_to_reorder_id_and_seq_jj.end(), smallersecond);

				for(unsigned int k = 0; k < items_to_reorder_positions_ii.size(); k++)
					sol->seq_[items_to_reorder_positions_ii[k]] = items_to_reorder_id_and_seq_ii[k].first;

				for(unsigned int k = 0; k < items_to_reorder_positions_jj.size(); k++)
					sol->seq_[items_to_reorder_positions_jj[k]] = items_to_reorder_id_and_seq_jj[k].first;


				vector<int> eval_plate_result;
				if(EVAL_TO_USE == 1) eval_plate_result = evalBatch0::evalPlateFromCut1(sol, ii);
				if(EVAL_TO_USE == 2) eval_plate_result = evalBatch2::evalPlate(sol, index_start);

				int total_item_surface_in_plate = eval_plate_result[0];
				int i_start_nextplate = eval_plate_result[1];
				int last_cut_x = eval_plate_result[2];

				SwapMove m(id++, ii, jj, i_start_nextplate, total_item_surface_in_plate, last_cut_x);
				m.value = calcval(total_item_surface_in_plate, i_start_nextplate, last_cut_x, sol);

				list_of_moves.push_back(m);

				sol->seq_ = seqcopy;

				nmbEvaluated++;

				if(shake == true && total_item_surface_in_plate >= currentPlateSurf)
				{
				    sort(list_of_moves.begin(), list_of_moves.end(), greater_value);
					return list_of_moves;
				}

			}
		}
	}

    sort(list_of_moves.begin(), list_of_moves.end(), greater_value);

	return list_of_moves;
}

pair<vector<int>, int> LocalSearchSwapPlate(Solution* sol, int& x_last_cut1_in_curr_plate, int index_start, int maxiter, int timeLimit)
{
	vector<int> eval_plate_result;
	eval_plate_result = eval_plate(sol, index_start);

	int index_end = eval_plate_result[1] - 1;
	int index_start_next_plate = index_end + 1;
	int_64 bestval = 0;
	bestval = calcval(eval_plate_result[0], eval_plate_result[1], eval_plate_result[2], sol); //???
	x_last_cut1_in_curr_plate = eval_plate_result[2];

	vector<int> resulting_seq = sol->seq_;

	long int t0 = time(0);

	sol->reset_tabu();
	do {

		if(time(0) - t0 >= timeLimit) break;

		int remTime = timeLimit - (time(0) - t0);

		vector<SwapMove> list_of_moves = evaluate_swap_moves_for_plate(sol, index_start, index_end, remTime);

		if (list_of_moves.size())
		{
			int nc = 1;
			loop2(x, 1, list_of_moves.size())
			{
				if(list_of_moves[x].plate_surface_used == list_of_moves[0].plate_surface_used &&
						list_of_moves[x].plate_last_cut_x == list_of_moves[0].plate_last_cut_x) nc++;
				else break;
			}
			int k = 0;
			if(nc > 1) k = sol->rand() % nc;

			index_start_next_plate = list_of_moves[k].index_start_next_plate;

			swapTwoElementsInList(list_of_moves[k].i, list_of_moves[k].j, sol->seq_);

			int si = DATA->item_list_[sol->seq_[list_of_moves[k].i]].stack;
			int sj = DATA->item_list_[sol->seq_[list_of_moves[k].j]].stack;
			std::vector<int> items_to_reorder_positions_ii;
			std::vector<pair<int,int> > items_to_reorder_id_and_seq_ii;
			std::vector<int> items_to_reorder_positions_jj;
			std::vector<pair<int,int> > items_to_reorder_id_and_seq_jj;
			for(unsigned int m = index_start; m < sol->seq_.size(); m++)
			{
				if(DATA->item_list_[sol->seq_[m]].stack == si) {
				  items_to_reorder_positions_ii.push_back(m);
				  items_to_reorder_id_and_seq_ii.push_back(make_pair(sol->seq_[m], DATA->item_list_[sol->seq_[m]].seq));
				}
				if(DATA->item_list_[sol->seq_[m]].stack == sj) {
				   items_to_reorder_positions_jj.push_back(m);
				   items_to_reorder_id_and_seq_jj.push_back(make_pair(sol->seq_[m], DATA->item_list_[sol->seq_[m]].seq));
				}
				if(items_to_reorder_positions_ii.size() == DATA->stack[si].size() &&
				   items_to_reorder_positions_jj.size() == DATA->stack[sj].size()) break;
			}
			sort(items_to_reorder_id_and_seq_ii.begin(), items_to_reorder_id_and_seq_ii.end(), smallersecond);
			sort(items_to_reorder_id_and_seq_jj.begin(), items_to_reorder_id_and_seq_jj.end(), smallersecond);

			for(unsigned int k = 0; k < items_to_reorder_positions_ii.size(); k++)
				sol->seq_[items_to_reorder_positions_ii[k]] = items_to_reorder_id_and_seq_ii[k].first;

			for(unsigned int k = 0; k < items_to_reorder_positions_jj.size(); k++)
				sol->seq_[items_to_reorder_positions_jj[k]] = items_to_reorder_id_and_seq_jj[k].first;

			eval_plate_result = eval_plate(sol, index_start);
			index_end = eval_plate_result[1] - 1;
			index_start_next_plate = index_end + 1;

			sol->update_tabu(list_of_moves[k].i, list_of_moves[k].j);

			int_64 val = calcval(eval_plate_result[0], eval_plate_result[1], eval_plate_result[2], sol);
			if (val > bestval)
			{
				bestval = val;
				x_last_cut1_in_curr_plate = eval_plate_result[2];
				loop2(i, index_start, DATA->nb_items_) resulting_seq[i] = sol->seq_[i];
				sol->reset_tabu();
			}
		}
		else
		{
			break;
		}
	}
	while(sol->iter < maxiter);
		
	return pair<vector<int>, int>(resulting_seq, index_start_next_plate);
}

pair<vector<int>, int> shakePlate(Solution* sol, int& x_last_cut1_in_curr_plate, int index_start, int maxiter,
		int timeLimit,
		int nmbShakeMovesMax = 100)
{
	vector<int> eval_plate_result = eval_plate(sol, index_start);
	int index_end = eval_plate_result[1] - 1;
	int index_start_next_plate = index_end + 1;
	int bestsurf = 0;
	bestsurf = eval_plate_result[0];
	x_last_cut1_in_curr_plate = eval_plate_result[2];
	vector<int> resulting_seq = sol->seq_;
	int nmbShakeMoves = 0;
	sol->reset_tabu();
	long int startTime = time(0);
	do {

		vector<SwapMove> list_of_moves = evaluate_swap_moves_for_plate(sol, index_start, index_end, 2, true, bestsurf);

		if (list_of_moves.size() && list_of_moves[0].plate_surface_used >= bestsurf)
		{
			nmbShakeMoves++;

			int nc = 1;
			loop2(x, 1, list_of_moves.size())
			{
				if(list_of_moves[x].plate_surface_used == list_of_moves[0].plate_surface_used && list_of_moves[x].plate_last_cut_x == list_of_moves[0].plate_last_cut_x) nc++;
				else break;
			}
			int k = 0;
			if(nc > 1) k = sol->rand() % nc;

			index_start_next_plate = list_of_moves[k].index_start_next_plate;

			swapTwoElementsInList(list_of_moves[k].i, list_of_moves[k].j, sol->seq_);

			while(1)
			{
				bool swaped = false;
				for(unsigned int m = 0; m < sol->seq_.size(); m++)
				{
					for(unsigned int n = m + 1; n < sol->seq_.size(); n++)
					{
						if(DATA->item_list_[sol->seq_[m]].stack == DATA->item_list_[sol->seq_[n]].stack &&
						   DATA->item_list_[sol->seq_[m]].seq > DATA->item_list_[sol->seq_[n]].seq)
						{
							int p = sol->seq_[m];
							sol->seq_[m] = sol->seq_[n];
							sol->seq_[n] = p;
						}
					}
				}
				if(!swaped) break;
			}

			eval_plate_result = eval_plate(sol, index_start);
			index_end = eval_plate_result[1] - 1;
			index_start_next_plate = index_end + 1;
			sol->update_tabu(list_of_moves[k].i, list_of_moves[k].j);
			if (list_of_moves[k].plate_surface_used > bestsurf ||
					(list_of_moves[k].plate_surface_used == bestsurf && list_of_moves[k].plate_last_cut_x < x_last_cut1_in_curr_plate))
			{
				bestsurf = list_of_moves[k].plate_surface_used;
				x_last_cut1_in_curr_plate = list_of_moves[k].plate_last_cut_x;
				loop2(i, index_start, DATA->nb_items_) resulting_seq[i] = sol->seq_[i];
				sol->reset_tabu();
			}
		}
		else break;

		if(nmbShakeMoves >= nmbShakeMovesMax) break;

	}
	while(time(0) - startTime < timeLimit);

	return pair<vector<int>, int>(resulting_seq, index_start_next_plate);
}

vector<int> generateSeqByPerformingFewRandomMoves(Solution* sol, vector<int>& initial_seq, int nmbMoves, int timeLimit, int i0 = 0)
{
	vector<int> result(initial_seq.size());
	copyFirstVectorToAnother(initial_seq, result);
	long int t0 = time(0);
	loop(r, nmbMoves)
	{
		if(time(0) - t0 >= timeLimit) break;
		int delta = 0;
		vector<pair<int, int> > moveslist;
		loop2(ii, i0, DATA->nb_items_ - 1)
		{
			if(time(0) - t0 >= timeLimit) break;
			int i = result[ii];
			vector<bool> is(DATA->nb_stacks_, false);
			loop2(jj, ii + 1, DATA->nb_items_)
			{
				if(time(0) - t0 >= timeLimit) break;
				int j = result[jj];
				if(DATA->item_list_[i].stack == DATA->item_list_[j].stack) break;
				if(is[DATA->item_list_[j].stack] == true) continue;
				is[DATA->item_list_[j].stack] = true;
				if (jj - ii > delta)
				{
					moveslist.push_back(pair<int, int>(ii, jj));
				}
			}
		}
		if (moveslist.size())
		{
			int k = sol->rand() % moveslist.size();
			swapTwoElementsInList(moveslist[k].first, moveslist[k].second, result);
		}
		else {
			break;
		}
	}
	return result;
}

class SolutionPartial {
public:
	std::vector<int> seq;
	int last_plate_index_;
	int nmb_items_;
	int_64 value_for_last_plate;
	int surf_last_plate;
	SolutionPartial(std::vector<int> s, int nb, int ni, int_64 val, int surf) : seq(s), last_plate_index_(nb), nmb_items_(ni),
			value_for_last_plate(val), surf_last_plate(surf) {}
};

std::vector<SolutionPartial> solution_pool_global;
std::vector<std::vector<SolutionPartial> > solution_pool_for_thread(NMB_THREADS);

bool greater_value_last_plate(const SolutionPartial& p1, const SolutionPartial& p2)
{
	return p1.value_for_last_plate > p2.value_for_last_plate;
}


int MAX_NMB_LOOPS_IN_FILL_BIN_MULTITHREAD = 10000;
int MAX_NMB_LOOPS_IN_FILL_BIN = 10000;
int TIME_LIMIT_FOR_CURR_SOLUTION = 3600;

void fillBin(int threadID, int curr_plate, int i0, int timeLimit, bool print = true,
		int maxNmbLoopsInFillBin = 5000)
{
	Solution* sol = SOL_THREAD[threadID];
	vector<int> global_best_seq = sol->seq_;
	vector<int> solutionToResetFrom = sol->seq_;
	int_64 best_obj(LONG_LONG_MAX);
	vector<SolutionPartial>& solution_pool = solution_pool_for_thread[threadID];

	long int startTime = time(0);
	int loop = 0;
	while (time(0) - startTime < timeLimit)
	{

		if(loop >= MAX_NMB_LOOPS_IN_FILL_BIN) break;
		if(loop >= maxNmbLoopsInFillBin) break;

		vector<int> eval_plate_result = eval_plate(sol, i0);

		int totalsurf = eval_plate_result[0];
		int xlast = eval_plate_result[2];
		int_64 best_obj_loop = - calcval(totalsurf, eval_plate_result[1], xlast, sol);

		solution_pool.push_back(SolutionPartial(sol->seq_, curr_plate, eval_plate_result[1],
				calcval(totalsurf, eval_plate_result[1], xlast, sol), totalsurf));

		if(best_obj_loop < best_obj)
		{
			best_obj = best_obj_loop;
			global_best_seq = sol->seq_;
		}

		int maxiter = sol->params_->maxiterTSinitial;

		while (1)
		{
			int x_last_cut1_in_curr_plate;
			pair<vector<int>, int> res;

			for(int sh = 0; sh < 2; sh++)
			{
				int remTime = timeLimit - (time(0) - startTime);
				if(remTime < 0) remTime = 0;

				shakePlate(sol, x_last_cut1_in_curr_plate, i0, maxiter, remTime,  5);

				remTime = timeLimit - (time(0) - startTime);
				if(remTime < 0) remTime = 0;

				res = LocalSearchSwapPlate(sol, x_last_cut1_in_curr_plate, i0, maxiter, remTime);
			}
			sol->seq_ = res.first;

			eval_plate_result = eval_plate(sol, i0);

			int totalsurf = eval_plate_result[0];
			int xlast = eval_plate_result[2];
			int_64 obj_loop = - calcval(totalsurf, eval_plate_result[1], xlast, sol);

			solution_pool.push_back(SolutionPartial(sol->seq_, curr_plate, eval_plate_result[1],
					calcval(totalsurf, eval_plate_result[1], xlast, sol), totalsurf));

			if(obj_loop < best_obj_loop)
			{
				maxiter += sol->params_->iterTSdelta;
				best_obj_loop = obj_loop;
				if(obj_loop < best_obj)
				{
					best_obj = obj_loop;
					global_best_seq = sol->seq_;
				}
			}
			else break;

			if(time(0) - startTime >= timeLimit) break;

		}

		sort(solution_pool.begin(), solution_pool.end(), greater_value_last_plate);
		for(int i = 0; i < solution_pool.size() - 1; i++)
		{
			if(solution_pool[i].value_for_last_plate == solution_pool[i + 1].value_for_last_plate)
				solution_pool.erase(solution_pool.begin() + (i + 1));
		}

		int k = 50;
		if(TIME_LIMIT_FOR_CURR_SOLUTION < 300) k = 5;

		if(k > solution_pool.size()) k = solution_pool.size();

		solutionToResetFrom = solution_pool[sol->rand() % k].seq;

		// few random moves from solutionToResetFrom or completely random sequence
		int N = 2;
		if(loop > 5) N = 10;
		if((sol->rand() % N) > 0)
		{
			sol->seq_ = generateSeqByPerformingFewRandomMoves(sol, solutionToResetFrom, sol->rand() % 15, 3, i0);
		}
		else
		{
			generateRandomFeasibleSequenceFromPosition(sol, sol->seq_, i0);
		}

		loop++;
	}

}

long int TOTALWASTE = 0;

std::vector<std::vector<SolutionPartial> > solutions_for_bin(NMB_PLATES);


int fillBinMultiThread(Data* data, int curr_plate, int i0, int timeLimit,
		int TOTALWASTEbeforeplate = -1,
		int maxNmbLoops = 5000, int maxNmbLoopsInFillBin = 5000)
{

	if(TOTALWASTEbeforeplate < 0) TOTALWASTEbeforeplate = TOTALWASTE;

	evalBatch0::initDefectsPlate(curr_plate);
	gensol::initDefectsPlate(curr_plate);
	evalBatch2::initDefectsPlate(curr_plate);

	solution_pool_global.clear();

	int timeForThread = 5;
	long int startTime = time(0);
	int loop = 0;
	int besttotalsurf = 0;

	while (time(0) - startTime < timeLimit)
	{

		loop++;

		thread thr[NMB_THREADS];

		int remtime = timeLimit - (time(0) - startTime);
		if(remtime < timeForThread) timeForThread = remtime;

		for (int i = 0; i < NMB_THREADS; ++i)
		{
			solution_pool_for_thread[i] = solution_pool_global;
			thr[i] = thread(fillBin, i, curr_plate, i0, timeForThread, 0, maxNmbLoopsInFillBin);
		}

		for (int i = 0; i < NMB_THREADS; ++i) thr[i].join();

		for (int i = 0; i < NMB_THREADS; ++i)
		{
			for(int j = 0; j < solution_pool_for_thread[i].size() && j < 2000; j++)
			{
				solution_pool_global.push_back(solution_pool_for_thread[i][j]);
			}
		}

		sort(solution_pool_global.begin(), solution_pool_global.end(), greater_value_last_plate);

		for(int i = 0; i < solution_pool_global.size() - 1; i++)
		{
			if(solution_pool_global[i].value_for_last_plate == solution_pool_global[i + 1].value_for_last_plate)
			{
				solution_pool_global.erase(solution_pool_global.begin() + (i + 1));
				i--;
			}
		}

		if(solution_pool_global.size() > 200)
			solution_pool_global.erase(solution_pool_global.begin() + 200, solution_pool_global.end());


		int totalsurf = 0;
		if(solution_pool_global.size()) totalsurf = solution_pool_global[0].surf_last_plate;

		if(besttotalsurf < totalsurf)
		{
			int W = PLATE_WIDTH;
			if(SOL_THREAD[0]->i0nextcut >= DATA->nb_items_) W = SOL_THREAD[0]->x0nextcut;
			TOTALWASTE = TOTALWASTEbeforeplate + (PLATE_HEIGHT * W - totalsurf);
			besttotalsurf = totalsurf;
		}
		else
		{
			TOTALWASTE = TOTALWASTEbeforeplate + (PLATE_HEIGHT * PLATE_WIDTH - totalsurf);
		}


		if(SOL_THREAD[0]->i0nextcut >= DATA->nb_items_) // don't spent too much time on last bin optimization
		{
			break;
		}

		if(loop >= MAX_NMB_LOOPS_IN_FILL_BIN_MULTITHREAD) break;
		if(loop >= maxNmbLoops) break;

	}


	if(solution_pool_global.size())
	{
		for(int t = 0; t < NMB_THREADS; t++) { SOL_THREAD[t]->seq_ = solution_pool_global[0].seq; }
	}

	vector<int> eval_plate_result = eval_plate(SOL_THREAD[0], i0);
	return eval_plate_result[1];
}


vector<int> solve(int sol_index, int timeLimit, int startP = 0, int startI = 0,
		bool reoptimizeLast2Bins = true,
		int avgLimitPerBin = -1,
		int timeLimitForMainOpt = -1)
{

	TOTALWASTE = 0;

	loop(i, NMB_PLATES) solutions_for_bin[i].clear();

	int startTime = time(0);

	for(int t = 0; t < NMB_THREADS; t++)
	{
		if(SOL_THREAD[t]->seq_[0] < 0) {
			loop(i, DATA->nb_items_) SOL_THREAD[t]->seq_[i] = i;
			generateRandomFeasibleSequence(SOL_THREAD[t]->seq_);
		}
	}

	int nmbBinsEstimated = (int) (DATA->LB_nmbBins_DOUBLE_ + 0.06 * DATA->LB_nmbBins_DOUBLE_);

	if(timeLimitForMainOpt < 0) timeLimitForMainOpt = timeLimit;

	if(nmbBinsEstimated == 0) nmbBinsEstimated = 1;

	int timeLimitPerBinAvg = timeLimitForMainOpt / nmbBinsEstimated;
	if(avgLimitPerBin > 0) timeLimitPerBinAvg = avgLimitPerBin;

	int timelimitbin1 = timeLimitPerBinAvg;
	int timelimitbinN = timeLimitPerBinAvg;

	int P(startP), S(startI);
	int timeDelta = (timelimitbin1 - timelimitbinN) / nmbBinsEstimated;
	int tL = timelimitbin1;
	while(true)
	{
		vector<int> bestseq;
		int bestSnext(S), TOTALWASTEbeforeplate(TOTALWASTE),
				TOTALWASTEbest(TOTALWASTE + PLATE_SURFACE);

		int_64 bestval = 0;

		int LLL = (timeLimit < 900) ? 1 : 2;
		for(unsigned int kk = 0; kk < LLL; kk++)
		{
			srand(SEED + 1000 * sol_index + 10 * P + kk + 2222);
			for(int t = 0; t < NMB_THREADS; t++) SOL_THREAD[t]->fill_rand_list();
			for(int t = 0; t < NMB_THREADS; t++) SOL_THREAD[t]->R = rand() % SOL_THREAD[t]->rand_list.size();

			for(int t = 0; t < NMB_THREADS; t++) generateRandomFeasibleSequenceFromPosition(SOL_THREAD[t]->seq_, S);
			int timeLmt = tL / LLL;
			if(timeLmt <= 0) timeLmt = 1;

			int snext = fillBinMultiThread(DATA, P, S, timeLmt, TOTALWASTEbeforeplate);

			vector<int> eval_plate_result = eval_plate(SOL_THREAD[0], S);
			int_64 val = calcval(eval_plate_result[0], eval_plate_result[1], eval_plate_result[2], SOL_THREAD[0]);
			if(val > bestval)
			{
				bestval = val;
				bestseq = SOL_THREAD[0]->seq_;
				bestSnext = snext;
				if(TOTALWASTE < TOTALWASTEbest) TOTALWASTEbest = TOTALWASTE;
			}

		}
		S = bestSnext;
		for(int t = 0; t < NMB_THREADS; t++) SOL_THREAD[t]->seq_ = bestseq;
		TOTALWASTE = TOTALWASTEbest;

		for(unsigned int i = 0; i < solution_pool_global.size(); i++) {
			solutions_for_bin[P].push_back(solution_pool_global[i]);
		}

		if(S >= DATA->nb_items_) break;
		if(S < 0) break;
		P++;
		tL -= timeDelta;

		double currWaste = ((double) (TOTALWASTE)) / (PLATE_HEIGHT * PLATE_WIDTH * P);
		int nmbBinsEstimated = (int) (DATA->LB_nmbBins_DOUBLE_ + currWaste * DATA->LB_nmbBins_DOUBLE_);
		int remNmbBinsEstimated = nmbBinsEstimated - P;
		if(remNmbBinsEstimated <= 0) remNmbBinsEstimated = 1;
		int remTime = timeLimitForMainOpt - (time(0) - startTime);
		if(remTime <= 0) remTime = 10;
		if(tL > remTime / remNmbBinsEstimated)
		{
			tL = remTime / remNmbBinsEstimated;
			if(tL <= 0) tL = 2;
		}
	}

	std::vector<int> best_seq = SOL_THREAD[0]->seq_;

	// re-optimize last 2 bins
	if(reoptimizeLast2Bins && P >= 2
		&& solutions_for_bin[P - 1].size() > 1 && solutions_for_bin[P - 1][0].surf_last_plate > 0)
	{

		if(time(0) - startTime >= timeLimit) timeLimit = time(0) - startTime + 1;

		//simple re-optimization
		int bestObj = evalBatch0::evalBatch(SOL_THREAD[0]) * PLATE_HEIGHT - DATA->total_surface_of_items_;
		int MM = 20;
		for(int k = 0; k < solutions_for_bin[P - 1].size() && k < MM; k++)
		{
			if(time(0) - startTime >= timeLimit) break;
			for(int l = 0; l < 2; l++)
			{
				if(time(0) - startTime >= timeLimit) break;

				for(int t = 0; t < NMB_THREADS; t++) { SOL_THREAD[t]->seq_ = solutions_for_bin[P - 1][k].seq; }

				fillBinMultiThread(DATA, P,  solutions_for_bin[P - 1][k].nmb_items_, 1, TOTALWASTE, 1, 10);

				int obj = evalBatch0::evalBatch(SOL_THREAD[0]) * PLATE_HEIGHT - DATA->total_surface_of_items_;

				if(obj < bestObj)
				{
					if(gensol::evalBatch(SOL_THREAD[0]) >= 0)
					{
						bestObj = obj;
						best_seq = SOL_THREAD[0]->seq_;
					}
				}
			}
		}

		// more complicated, solve last 2 bins again
		int P0 = P - 1;
		P = P0;
		S = solutions_for_bin[P0 - 1][0].nmb_items_;
		int NN = 10;
		for(unsigned int k = 0; k < NN; k++)
		{
			if(time(0) - startTime >= timeLimit) break;
			for(int t = 0; t < NMB_THREADS; t++) { SOL_THREAD[t]->seq_ = solutions_for_bin[P0 - 1][0].seq; }
			tL = 5;
			if(time(0) - startTime >= timeLimit - 10) break;
			P = P0;
			S = solutions_for_bin[P0 - 1][0].nmb_items_;
			while(true)
			{
				int snext = fillBinMultiThread(DATA, P, S, tL);
				S = snext;
				if(S >= DATA->nb_items_) break;
				if(S < 0) break;
				P++;
				tL = 3;
			}
			int obj = evalBatch0::evalBatch(SOL_THREAD[0]) * PLATE_HEIGHT - DATA->total_surface_of_items_;
			if(obj < bestObj)
			{
				int w = gensol::evalBatch(SOL_THREAD[0]);
				if(w >= 0 && w < bestObj)
				{
					bestObj = w;
					best_seq = SOL_THREAD[0]->seq_;
				}
			}
		}

	} // end re-optimize 2 bins

	return best_seq;

}


int solveFinal(int timeLimit)
{
	int startTime = time(0);
	vector<int> bestseq;
	int bestWaste = INT_MAX;
	Solution* SOL = new Solution(DATA);
	int MAX_SOLUTIONS = 10000;

	std::vector<int> timeLimitPerSolution;

	bool B5like = false;

	if(DATA->nb_stacks_ == 2 &&
	   (DATA->stack[0].size() > 5 * DATA->stack[1].size() ||
	    DATA->stack[1].size() > 5 * DATA->stack[0].size()) ) // B5-like
	{
		B5like = true;
		MAX_NMB_LOOPS_IN_FILL_BIN = 1;
		MAX_NMB_LOOPS_IN_FILL_BIN_MULTITHREAD = 1;
		EVAL_TO_USE = 2;
		NMB_THREADS = 1;
		USE_STACKS_IN_FITNESS_FUNCTION = true;
		for(int i = 0; i < MAX_SOLUTIONS; i++) timeLimitPerSolution.push_back(5);
	}
	else
	if(DATA->nb_items_ / DATA->LB_nmbBins_INT_ < 16 || DATA->nb_stacks_ < 5)
		// not many items per bin (B1 B9 B10 11)
		// or only few stacks
	{
		int tL;
		if(timeLimit > 3000)
		{
			timeLimitPerSolution.push_back( 1 * 60);
			timeLimitPerSolution.push_back( 1 * 60);
			timeLimitPerSolution.push_back( 2 * 60);
			timeLimitPerSolution.push_back( 3 * 60);
			timeLimitPerSolution.push_back( 5 * 60);
			tL = 120;
			for(int i = 0; i < MAX_SOLUTIONS; i++) timeLimitPerSolution.push_back(tL);
			MAX_NMB_LOOPS_IN_FILL_BIN = 20;
		}
		else
		{
			for(int i = 0; i < MAX_SOLUTIONS; i++) timeLimitPerSolution.push_back(55);
		}
	}
	else
	{
		if(timeLimit > 3000)
		{
			timeLimitPerSolution.push_back( 1 * 60);
			timeLimitPerSolution.push_back( 1 * 60);
			timeLimitPerSolution.push_back( 2 * 60);
			timeLimitPerSolution.push_back( 3 * 60);
			timeLimitPerSolution.push_back(timeLimit - 7 * 60);
			for(int i = 0; i < MAX_SOLUTIONS; i++) timeLimitPerSolution.push_back(timeLimit - 7 * 60);
			MAX_NMB_LOOPS_IN_FILL_BIN = 20;
		}
		else
		{
			timeLimitPerSolution.push_back(10);
			timeLimitPerSolution.push_back(20);
			timeLimitPerSolution.push_back(30);
			timeLimitPerSolution.push_back(timeLimit - 60);
			for(int i = 0; i < MAX_SOLUTIONS; i++) timeLimitPerSolution.push_back(60);
		}
	}

	int sol_index = 0;

	while(true) // solutions one by one
	{

		if(B5like)
		{
			int elapsedtime = time(0) - startTime;
			if(elapsedtime < timeLimit * 1.5 / 10)
			{
				MAX_NMB_LOOPS_IN_FILL_BIN = 5; MAX_NMB_LOOPS_IN_FILL_BIN_MULTITHREAD = 2;
				EVAL_TO_USE = 1; NMB_THREADS = 8; USE_STACKS_IN_FITNESS_FUNCTION = false;
			}
			else if(elapsedtime < timeLimit * 3 / 10)
			{
				MAX_NMB_LOOPS_IN_FILL_BIN = 5; MAX_NMB_LOOPS_IN_FILL_BIN_MULTITHREAD = 2;
				EVAL_TO_USE = 1; NMB_THREADS = 8; USE_STACKS_IN_FITNESS_FUNCTION = true;
				WEIGHT_STACKS_IN_FITNESS_FUNCTION = 1 + rand() % 200;
			}
			else if(elapsedtime < timeLimit * 4.5 / 10)
			{
				MAX_NMB_LOOPS_IN_FILL_BIN = 1; MAX_NMB_LOOPS_IN_FILL_BIN_MULTITHREAD = 1;
				EVAL_TO_USE = 2; NMB_THREADS = 1; USE_STACKS_IN_FITNESS_FUNCTION = false;
			}
			else
			{
				MAX_NMB_LOOPS_IN_FILL_BIN = 1; MAX_NMB_LOOPS_IN_FILL_BIN_MULTITHREAD = 1;
				EVAL_TO_USE = 2; NMB_THREADS = 1; USE_STACKS_IN_FITNESS_FUNCTION = true;
				WEIGHT_STACKS_IN_FITNESS_FUNCTION = 1 + rand() % 200;
			}
		}
		else if(DATA->nb_stacks_ < 5)
		{
			int elapsedtime = time(0) - startTime;
			if(elapsedtime < timeLimit * 5 / 10) USE_STACKS_IN_FITNESS_FUNCTION = false;
			else {
				USE_STACKS_IN_FITNESS_FUNCTION = true;
				WEIGHT_STACKS_IN_FITNESS_FUNCTION = 1 + rand() % 200;
			}
		}

		srand(SEED + sol_index + 1);
		for(int t = 0; t < NMB_THREADS; t++) SOL_THREAD[t]->fill_rand_list();

		TIME_LIMIT_FOR_CURR_SOLUTION = timeLimitPerSolution[sol_index];

		int remTime = timeLimit - (time(0) - startTime);
		if(remTime < timeLimitPerSolution[sol_index]) timeLimitPerSolution[sol_index] = remTime;

		for(int t = 0; t < NMB_THREADS; t++) loop(i, DATA->nb_items_) SOL_THREAD[t]->seq_[i] = -1;

		int timeLimitForMainOpt = timeLimitPerSolution[sol_index];
		if(timeLimitPerSolution[sol_index] > 100)
		{
			timeLimitForMainOpt = (int) (0.95 * timeLimitPerSolution[sol_index]);
		}

		vector<int> seq = solve(sol_index, timeLimitPerSolution[sol_index], 0, 0, 1, -1, timeLimitForMainOpt);

		SOL->seq_ = seq;

		int waste = evalBatch0::evalBatch(SOL) * PLATE_HEIGHT - DATA->total_surface_of_items_;

		if (gensol::evalBatch(SOL) >= 0 && waste < bestWaste)
		{
			bestWaste = waste;
			bestseq = SOL->seq_;
			string command = "mv " + SOL->solfile + " " + SOLUTION_FILE_NAME;
			system(command.c_str());
		}

		if(timeLimit - (time(0) - startTime) > 10)
		{
			int wasteGuess = evalBatch2::evalBatch(SOL, 5) * PLATE_HEIGHT - DATA->total_surface_of_items_;
			if(gensol2::evalBatch(SOL) > 0) waste = wasteGuess;
			if (waste >= 0 && waste < bestWaste)
			{
				bestWaste = waste;
				bestseq = SOL->seq_;
				string command = "mv " + SOL->solfile + " " + SOLUTION_FILE_NAME;
				system(command.c_str());
			}
		}

		std::cout << "solution " << setw(5) << sol_index + 1 << setw(12) << waste << setw(12) << bestWaste
				            << setw(7) << time(0) - startTime << std::endl;

		sol_index++;
		int elapsedTime = time(0) - startTime;
		remTime = timeLimit - elapsedTime;

		if(elapsedTime / sol_index > remTime)
		{
			if(remTime > 200)
			{
				timeLimitPerSolution[sol_index] = remTime - 100;
			}
			else break; // not enough time for next solution
		}

		if(bestWaste == INT_MAX) //feasible not found yet
		{
			if(timeLimitPerSolution[sol_index] > 10) timeLimitPerSolution[sol_index] = 10;
		}

	}

	SOL->seq_ = bestseq;

	std::cout << "-------------------------------" << std::endl;
	std::cout << "final waste " << setw(15) << bestWaste << std::endl;
	std::cout << "-------------------------------" << std::endl;

	for(int t = 0; t < NMB_THREADS; t++) SOL_THREAD[t]->seq_ = bestseq;

	return bestWaste;
}
