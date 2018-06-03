#ifndef LatGen
#define LatGen

#include "Consts.hh"
#include <random>

using namespace std;
using namespace voro;

class GeneticAlgoForSizesClass
{
public:
	static constexpr double x_min=-half_boxside;
	static constexpr double x_max=half_boxside;
	static constexpr double y_min=-half_boxside;
	static constexpr double y_max=half_boxside;
	static constexpr double z_min=-half_boxside;
	static constexpr double z_max=half_boxside;
	int population_size;
	default_random_engine generator;
	normal_distribution<double> distribution;

	
	GeneticAlgoForSizesClass(int pop_size);
	static double original_distribution(double point);
	double mutate_dist(double mutation_max_applitude_e, int dist_num = 0);
	double fitness_penalty(int points_number, double (*original_distribution)(double), map<double, double> current_distribution);
	double size_penalty(double* size_dist, int iter, int offspring_amount);
	double mutate(double val, double min, double max, bool adopted_shift = true, int dist_num = 0, double mult = 1.0);
	double reinit(double min, double max);
	void mutation(container ***con1, int iteration, bool reinit_flag = false, 
		double mutation_probability = mutation_probability_default, int dist_num = 0, sorted_points** exp = nullptr, double mult = 1.0);
	int tournament_selection(double* penalty, int iter, int offspring_amount, int selected = -1);
	point_for_crossover* get_ind_points(container* ind1, bool to_print = false);
	map<int, int> ind_to_ind(int ind1, int ind2, sorted_points** exp);
	map<int, point_for_crossover> ind_points_map(sorted_points sortp);
	void select_interchange_regions(container* ind1, map<int, point_for_crossover> *id1_to_coords, 
		int id, int from, int to, neighbors neg);
	void select_interchange_randomly(map<int, point_for_crossover> *id1_to_coords);
	container** crossover_by_mapping(container** con, double* penalty, int iter, int crossover_points = 1, 
		sorted_points** exp = nullptr, sorted_points** exp_off = nullptr, int from = 1, int to = particles, neighbors* neg = nullptr);
	double normal_dist(double x);
	double** compute_cell_sizes(container** con, neighbors** neg);
	void output_data(string filename, double* size_dist);
	void write_penalty_step(string filename, double penalties, double penalty);
	
};

#endif