#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
using namespace std;

#include "../solution_generator/make_solutionfile.h"
#include "../algo/cuttingGlass.h"

int NMB_THREADS = 8;

string BENCH_FOLDER = "instances/";
string BENCH_NAME;
string SOLUTION_FILE_NAME;
string batchPath = BENCH_NAME + "_batch.csv";
string defectsPath = BENCH_NAME + "_defects.csv";
string optParamsPath = BENCH_FOLDER + "global_param.csv";
long int SEED;
long int PROGRAM_START_TIME;
int TIME_LIMIT = 3600;
Data* DATA;
vector<Solution*> SOL_THREAD(NMB_THREADS, 0);
int EVAL_TO_USE = 1; // 1 = minimum waste per cut_1, 2 = 'dynamic' programming

namespace evalBatch0
{
  extern Data* DATA;
}

namespace gensol
{
   extern Data* DATA;
   extern std::string batchPath;
   extern std::string defectsPath;
   extern std::string optParamsPath;
}

namespace gensol2
{
   extern std::string batchPath;
   extern std::string defectsPath;
   extern std::string optParamsPath;
}

int solveFinal(int timeLimit);

int main(int argc, char** argv)
{
	PROGRAM_START_TIME = time(0);
	SEED = time(0) % 1000;

	// ./your_program_name -t time_limit -p instance_name -o solution_name -s seed -name
	for(int param = 1; param < argc; param++)
	{
		if(string(argv[param]) == string("-p")) {
			BENCH_NAME = argv[param+1];
		}
		else
		if(string(argv[param]) == string("-o")) {
			SOLUTION_FILE_NAME = argv[param+1];
		}
		else
		if(string(argv[param]) == string("-s")) {
			SEED = atoi(argv[param+1]);
		}
		else
		if(string(argv[param]) == string("-t"))
			TIME_LIMIT = atoi(argv[param+1]);
		else
		if(string(argv[param]) == string("-name"))
		{
			cout << "S21" << endl;
			if(argc == 2) exit(0);
		}
    }

	batchPath = BENCH_NAME + "_batch.csv";
	defectsPath = BENCH_NAME + "_defects.csv";
	optParamsPath = BENCH_FOLDER + "global_param.csv";
	gensol::batchPath = batchPath;
	gensol::defectsPath = defectsPath;
	gensol::optParamsPath = optParamsPath;
	gensol2::batchPath = batchPath;
	gensol2::defectsPath = defectsPath;
	gensol2::optParamsPath = optParamsPath;

	std::cout << "batchPath = " << batchPath << std::endl;

	DATA = new Data();
	DATA->readBatchAndDefectsFromFiles(batchPath, defectsPath);
	gensol::DATA = DATA;
	evalBatch0::DATA = DATA;

	for(int t = 0; t < NMB_THREADS; t++)
	{
		SOL_THREAD[t] = new Solution(DATA);
		ostringstream oss(""); oss << "sol" << t << ".csv";
		SOL_THREAD[t]->solfile = oss.str();
	}

	srand(SEED);

	int realTimeLimit = TIME_LIMIT;
	cout << "TIME_LIMIT " << TIME_LIMIT << endl;
	if(realTimeLimit > 3000) realTimeLimit -= 20;
	if(realTimeLimit > 170) realTimeLimit -= 5;

	double OBJ = solveFinal(realTimeLimit);

	std::cout << "--------------------------------\n";
	std::cout << "         OBJ " << std::setw(15) << std::setprecision(15) << OBJ << std::endl;
	std::cout << "     OBJBINS " << std::setw(15) << std::setprecision(5) << (OBJ + DATA->total_surface_of_items_) / (PLATE_WIDTH * PLATE_HEIGHT) << std::endl;
	std::cout << "      LBBINS " << std::setw(15) << std::setprecision(5) << DATA->LB_nmbBins_DOUBLE_ << std::endl;
	std::cout << "--------------------------------\n";

	return 0;

}

