#ifndef ALGO_CG_H_
#define ALGO_CG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <climits>
#include <algorithm>
#include <vector>
#include <ostream>
#include <iomanip>

#define PLATE_WIDTH 6000
#define PLATE_WIDTH_20 5980
#define PLATE_HEIGHT 3210
#define PLATE_HEIGHT_20 3190
#define PLATE_SURFACE 19260000
#define TWENTY 20
#define MIN_CUT_SIZE 100
#define MAX_PLATES 100
#define MAX_ITEMS 700
#define MAX_DEFECTS 1000
#define NMB_PLATES 100
#define MAXCUTSIZE 3500


class Solution;


namespace evalBatch2
{
  void initDefectsPlate(int pId);
}
namespace gensol
{
  void initDefectsPlate(int pId);
}


int nextY_Up(int x,int y,int w,int h);
int nextX_Right(int x,int y,int w,int h);
void swapTwoElementsInList(int p, int q, int list[]);
void generateRandomFeasibleSequence(std::vector<int>& SBS);

#define loop(x, n) for(int x = 0; x < n; ++ x)
#define loop2(x, x0, n) for(int x = x0; x < n; ++ x)

struct Item {
	int id;
	int length;
	int width;
	int stack;
	int seq;
	int surf;
	int w;
};

struct Defect {
	int id;
	int plateId;
	int x;
	int y;
	int width;
	int height;
};

class Data {


public:

	// main data
	int nb_items_;
	int nb_stacks_;
	int nb_defects_;

	std::vector<Item> item_list_;
	std::vector<Defect> defect_list_;
	std::vector<int> nmb_items_in_stack;
	std::vector<std::vector<int> > stack;

	// some useful static data
	int min_item_size_;
	int min_item_surf_;
	int total_surface_of_items_;

	double LB_nmbBins_DOUBLE_;
	int LB_nmbBins_INT_;
	int LB_x_;

	Data() : nb_items_(0), nb_stacks_(0), nb_defects_(0),
			 min_item_size_(100000), min_item_surf_(100000000), total_surface_of_items_(0),
			 LB_nmbBins_INT_(0), LB_nmbBins_DOUBLE_(0), LB_x_(0) {}

	void readBatchAndDefectsFromFiles(std::string filenameBatch, std::string filenameDefects)
	{

		//read batch

		char buff[1000];
		int a,b,c,d,e;
		FILE *f;
		f = fopen(filenameBatch.c_str(), "r");
		if(f == NULL)
		{
			printf("\nBatch file %s not found", filenameBatch.c_str());
			exit(0);
		}
		fgets(buff, 1000, f);
		nb_items_ = 0;
		min_item_size_ = PLATE_WIDTH;
		min_item_surf_ = PLATE_WIDTH * PLATE_HEIGHT;
		while(fscanf(f,"%d;%d;%d;%d;%d",&a,&b,&c,&d,&e) == 5)
		{
			Item it;
			it.id = a;
			it.length = (b <= c) ? b : c;
			it.width = (c < b) ? b : c;
			it.stack = d;
			it.seq = e;
			if(b < min_item_size_) min_item_size_ = b;
			if(c < min_item_size_) min_item_size_ = c;
			if (b < c) it.w = (3500 * c) + b;
			else it.w = (3500 * b) + c;
			it.surf = it.length * it.width;
			if(it.surf < min_item_surf_) min_item_surf_ = it.surf;
			nb_items_++;
			item_list_.push_back(it);
		}
		fclose(f);

		// count stacks
		nb_stacks_ = 0;
		std::vector<int> ids(nb_items_, 0);
		for (int i = 0; i < nb_items_; i++)
		{
			if(ids[item_list_[i].stack] == 0)
			{
				nb_stacks_++;
				ids[item_list_[i].stack] = 1;
			}
		}


		// calculate surface and bounds
		double r;
		total_surface_of_items_ = 0;
		for (int i = 0; i < nb_items_; i++)
		{
			item_list_[i].surf = item_list_[i].length * item_list_[i].width;
			total_surface_of_items_ += item_list_[i].surf;
		}

		r = ((double) total_surface_of_items_) / (PLATE_HEIGHT * PLATE_WIDTH);

		LB_nmbBins_DOUBLE_ = r;
		LB_nmbBins_INT_ = (int) LB_nmbBins_DOUBLE_;
		if(LB_nmbBins_DOUBLE_ > LB_nmbBins_INT_) LB_nmbBins_INT_++;

		r -= (int) LB_nmbBins_DOUBLE_;
		r *= PLATE_WIDTH;
		LB_x_ = (int) r;
		if (r > LB_x_) LB_x_++;


		// read defects
		char buff2[1000];
		int a2, b2;
		double c2, d2, e2, f2;
		FILE *g;
		g = fopen(filenameDefects.c_str(), "r");
		if(g == NULL)
		{
			printf("\nDefects file %s not found", filenameDefects.c_str());
			exit(0);
		}

		fgets(buff2, 1000, g);
		nb_defects_ = 0;
		while(fscanf(g, "%d;%d;%lf;%lf;%lf;%lf",&a2,&b2,&c2,&d2,&e2,&f2) == 6)
		{
			Defect def;
			def.id = a2;
			def.plateId = b2;
			def.x = c2;
			def.y = d2;
			def.width = e2;
			def.height = f2;
			defect_list_.push_back(def);
			nb_defects_++;
		}
		fclose(g);

		//genStacks
		nmb_items_in_stack.resize(nb_stacks_, 0);
		stack.resize(nb_stacks_, std::vector<int>(0));
		for(unsigned i = 0; i < nb_items_; i++)
		{
			stack[item_list_[i].stack].push_back(i);
		}
		for(unsigned i = 0; i < nb_stacks_; i++)
		{
			nmb_items_in_stack[i] = stack[i].size();
		}

	}

	std::string toString()
	{
		std::ostringstream oss("");
		oss << "nmb_items_ " << nb_items_ << "\n";
		oss << "min_item_size_ " << min_item_size_ << "\n";
		oss << "min_item_surf_ " << min_item_surf_ << "\n";
		oss << "list_of_items_.size " << item_list_.size() << "\n";
		oss << "nmb_stacks_ " << nb_stacks_ << "\n";
		oss << "nmb_defects_ " << nb_defects_ << "\n";
		oss << "LB_nmbBins_DOUBLE_ " << LB_nmbBins_DOUBLE_ << "\n";
		oss << "LB_nmbBins_INT_ " << LB_nmbBins_INT_ << "\n";
		oss << "LB_x_ " << LB_x_ << "\n";
		return oss.str();
	}
};

class Parameters {
 public:
   int MAXCUT3;
   int maxiterTSinitial;
   int iterTSdelta;
 public:
   Parameters() : MAXCUT3(5), maxiterTSinitial(4), iterTSdelta(0)  {}
};

#define MAX_ITER 4000

typedef struct {
	int i0cut1;
	int x0cut1,xmax;
	int surf;
} CUT1;

class Solution {

public:
	// main data
	Data* data_;
	Parameters* params_;
	std::vector<int> seq_;


	int R = 0;
	std::vector<int> rand_list;


	Solution(Data* d) : data_(d)
	{
		iter = 0;
		tabuTenure = 7;
		seq_.resize(data_->nb_items_, -1);
		params_ = new Parameters();
		params_->MAXCUT3 = 7;
		params_->maxiterTSinitial = 4;
		params_->iterTSdelta = 0;

		for(int i = 0; i < 3000000; i++) {
			rand_list.push_back(std::rand());
		}

	}

	int rand()
	{
		R = (R + 1) % rand_list.size();
		return rand_list[R];
	}

	void fill_rand_list()
	{
		for(int i = 0; i < rand_list.size(); i++) {
			rand_list[i] = std::rand();
		}
		R = 0;
	}


	std::string toString()
	{
		std::ostringstream oss("");
		for(unsigned int i = 0; i < seq_.size(); i++)
		{
			oss << std::setw(4) << seq_[i];
		}
		oss << "\n";
		return oss.str();
	}

	std::string solfile = "sol.csv";

	int iter;
	int tabu[MAX_ITEMS][MAX_ITEMS];
	int tabuTenure;

	void update_tabu(int ii,int jj)
	{
		int j;
		iter++;
		j = seq_[jj];
		tabu[j][ii]=iter+tabuTenure;

	}
	void reset_tabu(void)
	{
		iter = 0;
		for(int i = 0; i < data_->nb_items_; i++)
		 for(int j = 0; j < data_->nb_items_; j++)
		  tabu[i][j] = -1;
	}

	int xcut1[33];
	int ycut2[33][35];
	int xmaxcut1;
	int wastemin;
	int surfCut1;
	int x0cut1;
	int i0nextcut;
	int x0nextcut;
	int numplate;
	int numcut1;
	int ycut4[33];
	int nbcut4[33];
	int cut1OfIdx[MAX_ITEMS];
	int plateOfIdx[MAX_ITEMS];
	CUT1 tcut1[60];

	int ymaxcut2[33];

	int x0cut2[33][35];
	int x1cut2[33][35];
	int y0cut2[33][35];
	int y1cut2[33][35];
	int nbicut2[35];

	int nbdp;
	int Dx0[15];
	int Dy0[15];
	int Dx1[15];
	int Dy1[15];

	int from_i0;
	int to_i1;

};


#endif /* ALGO_CG_H_ */
