#ifndef GeneralClass
#define GeneralClass


#include <random>
#include "GeneralConsts.hh"
#include "ExtraStructAndFunc.hh"

using namespace std;
using namespace voro;

class General1Class
{
public:
	container* best_cont;
	sorted_points best_cont_points;
	orient_unit best_orient;
	map<int, atoms> at;
	
	General1Class();
	void writeAll();
	void loadAll();
	double ShiftDelta(double val);
	void GrainSizeGenerate(int argc, char *argv[]);
	void GrainOrientationGenerate(int argc, char *argv[]);
	void FillGrains();
	string RelaxFirst(int argc, char *argv[]);
	void RelaxSecond(int argc, char *argv[], string filename);
	void OutputForLAMMPS();
	
};

#endif