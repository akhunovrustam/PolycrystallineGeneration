#ifndef LatGen
#define LatGen

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "voro++.hh"
#include "GeneralConsts.hh"
#include "ExtraStructAndFunc.hh"

using namespace std;
using namespace voro;



class LatticeGeneratorClass
{
	static constexpr double x_min=0, x_max=half_boxside;
	static constexpr double y_min=0, y_max=half_boxside;
	static constexpr double z_min=0, z_max=half_boxside;
	container *con;
	
public:
	float ***planes_grains;
	
	LatticeGeneratorClass(container* cont);
	
	double rnd() {return double(rand())/RAND_MAX;}

	void calculateGrainPlanes(vector<int> &f_vert,vector<double> &v,int j, float** planes, float x, float y, float z, int planes_size);
	atoms generateLattice(float** grain, orient_unit angless, bool periodic = true);
	map<int, atoms> fillRandomlyAndBuildGrains(orient_unit angles);
};

#endif