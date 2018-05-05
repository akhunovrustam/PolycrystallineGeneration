#define _USE_MATH_DEFINES
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <ctime>
#include <string>

using namespace std;
using namespace voro;

struct point_for_crossover {
	double x;
	double y;
	double z;
	int block;
	int particle;
	bool is_interchanged;
};

struct euler_angles {
	double alpha;
	double beta;
	double gamma;
};

typedef map<string, euler_angles> orient_unit;
double rnd() {return double(rand())/RAND_MAX;}
