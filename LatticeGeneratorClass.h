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

using namespace std;
using namespace voro;

#define PI 3.14159265
#define LATTICE 0.2

class LatticeGeneratorClass
{
	static const double x_min=-1,x_max=1;
	static const double y_min=-1,y_max=1;
	static const double z_min=-1,z_max=1;
	container *con;
	
public:
	float ***planes_grains;
	unsigned int particles;
	
	LatticeGeneratorClass(int centers);
	
	double rnd() {return double(rand())/RAND_MAX;}

	void calculateGrainPlanes(vector<int> &f_vert,vector<double> &v,int j, float** planes, float x, float y, float z, int planes_size);
	int generateLattice(float** grain);
	void fillRandomlyAndBuildGrains();
};

#endif