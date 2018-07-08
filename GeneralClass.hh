#ifndef GeneralClass
#define GeneralClass


#include <random>
#include "GeneralConsts.hh"

using namespace std;
using namespace voro;

class GeneralClass
{
public:
	container* best_cont;
	sorted_points best_cont_points;
	orient_unit best_orient;
	atoms* at;
	
	GeneralClass();
	void GrainSizeGenerate();
	void GrainOrientationGenerate();
	void FillGrains();
	void RelaxFirst();
	void RelaxSecond();
	void OutputForLAMMPS();
	
};

#endif