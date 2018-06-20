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
#include "Consts.hh"

using namespace std;
using namespace voro;

struct config{
	int co; 		//crossover operator
	int md; 		//mutation distribution
	double mp; 		//mutation probability
	int ps; 		//population size
	double mult; 	//multiplication for mutation amplitude
	double from; 	//from how many particles make crossover
	double to; 		//to how many particles make crossover
	double thr1; 	//first threshold
	double thr2; 	//second threshold
	double fac; 	//threshold decrease factor
	int ps_fac1; //threshold decrease factor for ps
	int ps_fac2; //threshold decrease factor for ps
	int ps_fac;  //commond factor
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
	getline(f, s, '_');
	string mult = s;
	getline(f, s, '_');
	string from = s;
	getline(f, s, '_');
	string to = s;
	getline(f, s, '_');
	string threshold1 = s;
	getline(f, s, '_');
	string threshold2 = s;
	getline(f, s, '_');
	string factor = s;
	getline(f, s, '_');
	string ps_fac1 = s;
	getline(f, s, '_');
	string ps_fac2 = s;
	
	config cfg;
	if (co == "uniform") cfg.co = -1;
	else if (co == "1p") cfg.co = 1;
	else if (co == "2p") cfg.co = 2;
	
	if (md == "uniform") cfg.md = 0;
	else if (md == "normal") cfg.md = 1;
	else if (md == "reinit") cfg.md = 2;
	
	cfg.mp = stod(mp);
	
	cfg.ps = stoi(ps);
	
	cfg.mult = stod(mult);
	
	cfg.from = stoi(from);
	
	cfg.to = stoi(to);
	
	cfg.thr1 = stoi(threshold1);
	
	cfg.thr2 = stoi(threshold2);
	
	cfg.fac = stoi(factor);
	
	cfg.ps_fac1 = stoi(ps_fac1);
	
	cfg.ps_fac2 = stoi(ps_fac2);
	
	cfg.ps_fac = 1;
	return cfg;
}

void parse_args(int argc, char* argv[], string* prefix, config* cfg, int* population_size)
{
	if (argc == 3)
	{
		*prefix = argv[2];
		*prefix += "_SEED_";
		*prefix += to_string(atoi(argv[1]) * time(NULL));
		*prefix += "_exp";
	}
	else if (argc == 2)
	{
		*prefix = ((*prefix + "_SEED_") + to_string(atoi(argv[1]) * time(NULL))) + "_exp";
	}
	
	// exit(0);
	*cfg = parse_prefix(*prefix);
	*population_size = population_size_const;
	if (*prefix != ""){
		*population_size = (*cfg).ps;
	}
}

void create_dir(string* filename, string prefix, string folder = "results_size/")
{
	stringstream ss;
	
	//create folder to save output data
	time_t t = time(0);
	tm* now = localtime(&t);
    ss << now->tm_mday << "-" << (now->tm_mon + 1) << "-" << (now->tm_year + 1900) << "_" << (now->tm_hour) << "." << (now->tm_min);
	
	*filename = folder + prefix + ss.str();
	system(("mkdir " + *filename).c_str());
	
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


typedef map<int, vector<int>> neighbors;
typedef map<double, map<double, map<double, int>>> sorted_points;
typedef map<string, euler_angles> orient_unit;
double rnd() {return double(rand())/RAND_MAX;}

struct con_and_points{
	container* con;
	sorted_points pt;
};
