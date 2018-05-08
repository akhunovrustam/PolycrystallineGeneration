#ifndef LatGen
#define LatGen

#define half_boxside 0.5
#define penalty_steps 1000
#define penalty_step 0.025
#define max_iterations 3000
#define max_allowed_penalty 0.0001
#define surviving_size 0
#define population_size 32
#define particles 1000
#define mutation_probability 0.10
#define mutation_max_applitude 0.05
#define crossover_probability 0.05

class GeneticAlgoForSizesClass
{
public:
	static constexpr double x_min=-half_boxside;
	static constexpr double x_max=half_boxside;
	static constexpr double y_min=-half_boxside;
	static constexpr double y_max=half_boxside;
	static constexpr double z_min=-half_boxside;
	static constexpr double z_max=half_boxside;
	
	static double original_distribution(double point);
	double mutate_dist(double mutation_max_applitude_e, int dist_num = 0);
	double fitness_penalty(int points_number, double (*original_distribution)(double), std::map<double, double> current_distribution);
	double size_penalty(double* size_dist, int iter, int offspring_amount);
	double mutate(double val, double min, double max, bool adopted_shift = true);
	double reinit(double min, double max);
	void mutation(container ***con1, bool reinit_flag = false);
	int tournament_selection(double** size_dist, int iter, int offspring_amount, int selected = -1);
	point_for_crossover* get_ind_points(container* ind1, bool to_print = false);
	std::map<int, int> ind_to_ind(container* ind1, container* ind2);
	std::map<int, point_for_crossover> ind_points_map(container* ind);
	void select_interchange_regions(container* ind1, std::map<int, point_for_crossover> *id1_to_coords, int id);
	void select_interchange_randomly(std::map<int, point_for_crossover> *id1_to_coords);
	container** crossover_by_mapping(container** con, double** size_dist, int iter, int crossover_points = 1);
	double normal_dist(double x);
	double** compute_cell_sizes(container** con);
	void output_data(std::string filename, double* size_dist);
	void write_penalty_step(std::string filename, int penalties, double penalty);
	
};

#endif