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

double rnd() {return double(rand())/RAND_MAX;}
