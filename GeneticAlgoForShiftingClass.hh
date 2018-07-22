#ifndef GenAlgoShift
#define GenAlgoShift


#include <random>
#include "GeneralConsts.hh"

using namespace std;
using namespace voro;

class GeneticAlgoForShiftingClass
{
public:
	static constexpr double x_min=-half_boxside,x_max=half_boxside;
	static constexpr double y_min=-half_boxside,y_max=half_boxside;
	static constexpr double z_min=-half_boxside,z_max=half_boxside;
	static constexpr double a_min=0,a_max=M_PI / 2;
	int population_size;
	default_random_engine generator;
	normal_distribution<double> distribution;
	

	GeneticAlgoForShiftingClass(int pop_size);
	double writeAtoms(map<int, atoms> atm, string filename = "dump.lammpstrj");
	double pe(string filename, map<int, atoms> atm);
	void write_penalty_step(string filename, int penalties, double penalty);
	static double original_distribution(double point);
	double mutate_dist(double mutation_max_applitude_e, int dist_num = 0);
	double size_penalty(map<int, atoms> atm);
	double EnergyCalc(map<int, atoms> atm);
	double reinit(double min, double max);
	double reinit_angles();
	void mutation(map<int, atoms>** con1, config cfg);
	int tournament_selection(double* penalty, int selected = -1);
	point_for_crossover* get_ind_points(container* ind1, bool to_print = false);
	map<int, point_for_crossover> ind_points_map(sorted_points sortp);
	void select_interchange_regions(map<int, atoms> con, map<int, point_for_crossover> *id1_to_coords, int id, 
		int from, int to, neighbors neg);
	double ShiftDelta(double val);
	void select_interchange_randomly(map<int, point_for_crossover> *id1_to_coords);
	map<int, atoms>* crossover_by_mapping(double* penalty, map<int, atoms>* atms, config cfg, 
		neighbors neg, sorted_points best_cont_points);
	double normal_dist(double x);
	void output_data(string filename, orient_unit size_dist);
	orient_unit* relative_euler(container* con, orient_unit *abs_euler);
};

#endif