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

struct config{
	int co; 	//crossover operator
	int md; 	//mutation distribution
	double mp; 	//mutation probability
	int ps; 	//population size
};

config parse_prefix(string prefix)
{
	istringstream f(prefix.c_str());
    string s;
	getline(f, s, '_');
	string co = s;
	getline(f, s, '_');
	string md = s;
	getline(f, s, '_');
	string mp = s;
	getline(f, s, '_');
	string ps = s;
	
	config cfg;
	if (co == "uniform") cfg.co = -1;
	else if (co == "1p") cfg.co = 1;
	else if (co == "2p") cfg.co = 2;
	
	if (md == "uniform") cfg.md = 0;
	else if (md == "normal") cfg.md = 1;
	else if (md == "reinit") cfg.md = 2;
	
	cfg.mp = stod(mp);
	
	cfg.ps = stoi(ps);
	
	return cfg;
}

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
