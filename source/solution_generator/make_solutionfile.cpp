#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;
#include "make_solutionfile.h"

#define PLATE_WIDTH 6000
#define PLATE_HEIGHT 3210
#define TYPE_RESIDUAL -3
#define TYPE_WASTE -1
#define TYPE_USEFUL -2

/*
PLATE  ITEM   XNEXTCUT1   YNEXTCUT2      X0    Y0     X1      Y1
    0     0	        758         888 	  0     0	 100	 100
*/
class Node {
public:
	int plate_id_;
	int  node_id_ ;
	int        X_;
	int        Y_ ;
	int    width_ ;
	int   height_ ;
	int     type_ ;
	int      cut_ ;
	int   parent_ ;

	Node() {	plate_id_ =   -1;
	node_id_ =    0;
	X_ =    0;
	Y_ =    0;
	width_ = PLATE_WIDTH;
	   height_ = PLATE_HEIGHT;
	     type_ =   -2;
	      cut_ =    0;
	   parent_ =   -1;
}
	Node(int pid, int nid, int x, int y, int w, int h, int t, int c, int p) : plate_id_(pid), node_id_(nid), X_(x), Y_(y), width_(w),
			height_(h), type_(t), cut_(c), parent_(p) {}

};

void print_solution(vector<Node*>& nodes_, string filename)
{
	ofstream file(filename.c_str(), ios_base::out);
	file << "PLATE_ID;NODE_ID;X;Y;WIDTH;HEIGHT;TYPE;CUT;PARENT\n";
	for(int i = 0; i < nodes_.size(); i++)
	{
		Node* N = nodes_[i];
		file << N->plate_id_ << ";" << N->node_id_ << ";" << N->X_ << ";" << N->Y_ << ";" << N->width_ << ";"
			<< N->height_ << ";" << N->type_ << ";" << N->cut_ << ";";
		if(N->parent_ >= 0) file << N->parent_;
		file << endl;
	}
	file.close();
}

#include "../algo/cuttingGlass.h"
extern Data* DATA;

bool isXOnDefect(int pId, int x)
{
	for (int d = 0; d < DATA->nb_defects_; d++)
     if(DATA->defect_list_[d].plateId == pId)
	  if(x >= DATA->defect_list_[d].x && x < DATA->defect_list_[d].x + DATA->defect_list_[d].width) return true;

	return false;
}


/*
PLATE  ITEM   XNEXTCUT1   YNEXTCUT2      X0    Y0     X1      Y1
    0     0	        758         888 	  0     0	 100	 100
*/
bool make_solution_file(vector<vector<int> >& data, string solution_file_name)
{
	int nmbPlates = 0;
	for(unsigned int i = 0; i < data.size(); i++) if(data[i][0] >= nmbPlates) nmbPlates++;
	vector<Node*> nodes_;

	//try to solve the problem with too many 4-cuts by pushing the item up
	vector<int> UB(data.size(), 100000);
	vector<int> LB(data.size(), 0);
	for(unsigned int i = 0; i < data.size(); i++)
	{
		UB[i] = data[i][3];
		for(int j = i - 1; j >= 0 && data[i][0] == data[j][0]; j--)
		{
			if( data[i][2] == data[j][2] &&   data[i][3] != data[j][3])
			{
				LB[i] = data[j][3];
				break;
			}
		}
	}

	for(unsigned int i = 0; i < data.size(); i++)
	{
		if(LB[i] < data[i][5] && data[i][7] < UB[i])
		{
			int delta = UB[i] - data[i][7];
			data[i][5] += delta;
		    data[i][7] += delta;
		}
	}


	//PLATE  ITEM   XNEXTCUT1   YNEXTCUT2      X0    Y0     X1      Y1
	int currNodeID = 0;
	for(int pl = 0; pl < nmbPlates; pl++)
	{
		// plate node
		int plate_id = pl;
		Node* PLATE = new Node(plate_id, currNodeID, 0, 0, PLATE_WIDTH, PLATE_HEIGHT, TYPE_USEFUL, 0, -1);
		int plateNodeID = currNodeID;
		currNodeID++;
		nodes_.push_back(PLATE);

		int oneCutParentNodeID = -1;
		int twoCutParentNodeID = -1;
		int currentCut1Start  = 0;
		int currentCut2Start  = 0;
		int currentCut1End  = 0;
		int currentCut2End  = 0;
		int prev1CutEndX = 0;
		int prev2CutEndY = 0;

		for(unsigned int i = 0; i < data.size(); i++)
		{
			if(data[i][0] != plate_id) continue; // not on this plate

			bool is_new_1cut = (data[i][2] != currentCut1End);

			if(is_new_1cut) // new 1-cut (i.e. new vertical piece)
			{

				//add empty cut before if necessary
				if(data[i][2] - prev1CutEndX > 3500)
				{
					int emptyCutStart = prev1CutEndX;
					int emptyCutEnd = prev1CutEndX + (data[i][2] - prev1CutEndX - 3500);

					int XX = emptyCutEnd;
					while(XX < 6000)
					{
					    if(isXOnDefect(XX, plate_id) || XX - emptyCutStart < 20) XX++;
					    else break;
					}
					emptyCutEnd = XX;


					Node* N = new Node(plate_id, currNodeID, emptyCutStart, 0, emptyCutEnd - emptyCutStart, PLATE_HEIGHT, TYPE_WASTE, 1, plateNodeID);
					nodes_.push_back(N);
					prev1CutEndX = emptyCutEnd;
					currNodeID++;
					std::cout << "empty cut-1 on plate " << plate_id << " from " << emptyCutStart << " to " << emptyCutEnd << std::endl;
				}

				int xmax = data[i][2];
				currentCut1Start = prev1CutEndX;
				currentCut1End = data[i][2];
				Node* N = new Node(plate_id, currNodeID, currentCut1Start, 0, currentCut1End - currentCut1Start, PLATE_HEIGHT, TYPE_USEFUL, 1, plateNodeID);
				nodes_.push_back(N);
				prev1CutEndX = xmax;
				oneCutParentNodeID = currNodeID;
				currNodeID++;

				currentCut2Start  = 0;
				prev2CutEndY = 0;

				currentCut2End = 0;
			}

			bool is_new_2cut = (data[i][3] != currentCut2End);
			if(is_new_2cut) // new 2-cut (i.e. new horizontal piece i.e. 1-2cut rectangle)
			{
				int ymax = data[i][3];
				currentCut2Start = prev2CutEndY;
				currentCut2End = data[i][3];
				Node* N = new Node(plate_id, currNodeID, currentCut1Start, currentCut2Start, currentCut1End - currentCut1Start,
						currentCut2End - currentCut2Start, TYPE_USEFUL, 2, oneCutParentNodeID);
				nodes_.push_back(N);
				prev2CutEndY = ymax;
				twoCutParentNodeID = currNodeID;
				currNodeID++;

				int X0 = data[i][4];
				int Y0 = data[i][5];
				int X1 = data[i][6];
				int Y1 = data[i][7];

				unsigned int J0 = i;
				unsigned int J1 = i;

				if(X1 - X0 == nodes_.back()->width_ && Y1 - Y0 == nodes_.back()->height_) // complete cut-2 is an item
				{
					nodes_.back()->type_ = data[i][1];
				}
				else // cut items using 3-cut and 4-cut
				{
					for(unsigned int j = J0; j < data.size(); j++)
					{
						if(data[j][0] == data[J0][0]  && data[j][2] == data[J0][2] && data[j][3] == data[J0][3]) J1 = j;
						else break;
					}

					// rectangle contains items J0,....J1
					int currentCut3Start = currentCut1Start;

					// cut-3 waste before first item if necessary
					if(data[J0][4] > currentCut1Start)
					{
						Node* N = new Node(plate_id, currNodeID, currentCut1Start, currentCut2Start, data[J0][4] - currentCut1Start,
								currentCut2End - currentCut2Start, TYPE_WASTE, 3, twoCutParentNodeID);
						nodes_.push_back(N);
						currentCut3Start = data[J0][4];
						currNodeID++;
					}


					// 3-cuts for items and potentially waste between them
					for(unsigned int j = J0; j <= J1; j++)
					{

						currentCut3Start = data[j][4];

						//create useful node (item or item + waste)
						Node* N = new Node(plate_id, currNodeID, currentCut3Start, currentCut2Start, data[j][6] - data[j][4],
								currentCut2End - currentCut2Start, TYPE_USEFUL, 3, twoCutParentNodeID);
						nodes_.push_back(N);
						int threeCutParentNodeID = currNodeID;
						currNodeID++;

						//only item
						if(currentCut2End - currentCut2Start == data[j][7] - data[j][5])
						{
							nodes_.back()->type_ = data[j][1];
						}
						// or (item + waste) or 2 items
						else
						{
							// cut-4 waste before item
							if(currentCut2Start < data[j][5])
							{
								Node* N = new Node(plate_id, currNodeID, currentCut3Start, currentCut2Start, data[j][6] - data[j][4],
										data[j][5] - currentCut2Start, TYPE_WASTE, 4, threeCutParentNodeID);
								nodes_.push_back(N);
								currNodeID++;
							}

							//item
							Node* N = new Node(plate_id, currNodeID, currentCut3Start, data[j][5], data[j][6] - data[j][4],
									data[j][7] - data[j][5], data[j][1], 4, threeCutParentNodeID);
							nodes_.push_back(N);
							currNodeID++;


							if(currentCut2End > data[j][7])
							{
								if(j == J1 || data[j][4] != data[j + 1][4]) // cut-4 waste after item
								{
									Node* N = new Node(plate_id, currNodeID, currentCut3Start, data[j][7], data[j][6] - data[j][4],
											currentCut2End - data[j][7], TYPE_WASTE, 4, threeCutParentNodeID);
									nodes_.push_back(N);
									currNodeID++;

									if(currentCut2Start < data[j][5] && currentCut2End > data[j][7])
									{
										cout << "\nToo many 4-cuts!\n";
										cout << "item " << data[j][1] << "\n";
										cout << "plate " << data[j][0] << "\n";
										cout << "currentCut2Start " << currentCut2Start << "\n";
										cout << "item y start " << data[j][5] << "\n";
										cout << "currentCut2End " << currentCut2End << "\n";
										cout << "item y end " << data[j][7] << "\n";
										return false;
										//exit(0);
									}

								}
								else // another item
								{
									j = j + 1;
									Node* N = new Node(plate_id, currNodeID, currentCut3Start, data[j][5], data[j][6] - data[j][4],
											data[j][7] - data[j][5], data[j][1], 4, threeCutParentNodeID);
									nodes_.push_back(N);
									currNodeID++;

									if(currentCut2Start != data[j - 1][5] || currentCut2End != data[j][7])
									{
										cout << "\nCut 4 for 2 items not ok, waste exists!\n";
										cout << "j " << j << endl;
										cout << "J0 " << J0 << endl;
										cout << "J1 " << J1 << endl;
										cout << "currentCut2Start " << currentCut2Start << endl;
										cout << "data[j - 1][5] " << data[j - 1][5] << endl;
										cout << "data[j][7] " << data[j][7] << endl;
										cout << "currentCut2End " << currentCut2End << endl;
										return false;
										//exit(0);
									}

								}


							}

						}

						//cut-3 waste between items j and j + 1
						if(j < J1 && data[j][6] < data[j + 1][4])
						{

							Node* N = new Node(plate_id, currNodeID, data[j][6], currentCut2Start, data[j + 1][4] - data[j][6],
									currentCut2End - currentCut2Start, TYPE_WASTE, 3, twoCutParentNodeID);
							nodes_.push_back(N);
							currNodeID++;
							//cout << "CHECK add cut-3 waste between items";
							//getchar();
						}

					}

					// cut-3 waste after last item if necessary
					if(data[J1][6] < currentCut1End)
					{
						Node* N = new Node(plate_id, currNodeID, data[J1][6], currentCut2Start, currentCut1End - data[J1][6],
								currentCut2End - currentCut2Start, TYPE_WASTE, 3, twoCutParentNodeID);
						nodes_.push_back(N);
						currNodeID++;
					}


				}

				//add cut2 waste from last cut-2 to PLATE_HEIGHT
				bool is_last = (J1 == data.size() - 1 || data[J1][0] != data[J1 + 1][0] || data[J1 + 1][2] != currentCut1End);
				if(is_last && currentCut2End < PLATE_HEIGHT)
				{
					Node* N = new Node(plate_id, currNodeID, currentCut1Start, currentCut2End, currentCut1End - currentCut1Start,
							PLATE_HEIGHT - currentCut2End, TYPE_WASTE, 2, oneCutParentNodeID);
					nodes_.push_back(N);
					prev2CutEndY = PLATE_HEIGHT;
					twoCutParentNodeID = currNodeID;
					currNodeID++;
				}
			}

		}

		// add last 1-cut (WASTE or RESIDUAL)
		if(prev1CutEndX < PLATE_WIDTH)
		{
			int xmax = PLATE_WIDTH;
			int nodetype = TYPE_WASTE;
			if(pl == nmbPlates -1) nodetype = TYPE_RESIDUAL;
			Node* N = new Node(plate_id, currNodeID, prev1CutEndX, 0, xmax - prev1CutEndX, PLATE_HEIGHT, nodetype, 1, plateNodeID);
			nodes_.push_back(N);
			prev1CutEndX = xmax;
			currNodeID++;
		}
	}

	ostringstream filename2("");
	filename2 << solution_file_name;
	print_solution(nodes_, filename2.str());
	return true;
}

