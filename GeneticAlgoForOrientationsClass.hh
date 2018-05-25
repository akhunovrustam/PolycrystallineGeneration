#ifndef LatGen
#define LatGen

#define half_boxside 0.5
#define penalty_steps 20
#define penalty_step M_PI/penalty_steps/2.0
#define max_iterations 1000
#define max_allowed_penalty 0.0001
#define surviving_size 0
#define population_size_const 32
#define particles 1000
#define mutation_probability_default 0.10
#define mutation_max_applitude 0.05
#define crossover_probability 0.05

#include <random>

using namespace std;
using namespace voro;

class GeneticAlgoForOrientationsClass
{
public:
	static constexpr double x_min=-half_boxside,x_max=half_boxside;
	static constexpr double y_min=-half_boxside,y_max=half_boxside;
	static constexpr double z_min=-half_boxside,z_max=half_boxside;
	static constexpr double a_min=0,a_max=M_PI / 2;
	int population_size;
	default_random_engine generator;
	normal_distribution<double> distribution;
	
	GeneticAlgoForOrientationsClass(int pop_size);
	static double original_distribution(string point);
	double mutate_dist(double mutation_max_applitude_e, int dist_num = 0);
	double fitness_penalty(int points_number, double (*original_distribution)(string), 
		map<string, double> current_distribution, string output = "");
	double size_penalty(orient_unit parent_rel, string output = "");
	euler_angles mutate(euler_angles angles, double min, double max, bool adopted_shift = true, int dist_num = 0);
	double reinit(double min, double max);
	double reinit_angles();
	void mutation(orient_unit **con1, bool reinit_flag = false, 
		double mutation_probability = mutation_probability_default, int dist_num = 0);
	int tournament_selection(orient_unit* parents, int selected = -1);
	point_for_crossover* get_ind_points(container* ind1, bool to_print = false);
	map<int, int> ind_to_ind(container* ind1, container* ind2);
	map<int, point_for_crossover> ind_points_map(container* ind);
	void select_interchange_regions(container* con, map<int, point_for_crossover> *id1_to_coords, int id);

	void select_interchange_randomly(map<int, point_for_crossover> *id1_to_coords);
	orient_unit* crossover_by_mapping(container* con, orient_unit* parents, orient_unit* parents_rel, int crossover_points = 1);
	double normal_dist(double x);
	double** compute_cell_sizes(container** con);
	void output_data(string filename, double* size_dist);
	orient_unit* relative_euler(container* con, orient_unit *abs_euler);
};

#endif