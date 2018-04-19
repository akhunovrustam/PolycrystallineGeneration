// Polycrystal rebuild from the size distribution using genetics algorithm
//
// Author   : Rustam Akhunov (GUT Gdansk)
// Email    : akhunovrustam@mail.ru
// Date     : 

#include "voro++.hh"
#include "ExtraStructAndFunc.hh"
#include "GeneticAlgoForSizesClass.cc"

int main() {
	std::stringstream ss;
	
	//create folder to save output data
	std::time_t t = std::time(0);
	std::tm* now = std::localtime(&t);
    ss << now->tm_mday << "-" << (now->tm_mon + 1) << "-" << (now->tm_year + 1900) << "_" << (now->tm_hour) << "." << (now->tm_min);
	
	std::string filename = "results/" + ss.str();
	system(("mkdir " + filename).c_str());
	
	//initial variables etc
	srand (time(NULL));
	container **con;
	container **offspring;
	std::map<double, container*> common_pool;
	con = new container*[population_size];
	double x, y, z;
	double **real_sizes;
	double penalty[population_size];
	double min_penalty = 100;
	int iterations = 0;
	
	GeneticAlgoForSizesClass *algo = new GeneticAlgoForSizesClass();
	
	//initial population
	for (int i = 0; i < population_size; i++){
		con[i] = new container(algo->x_min,algo->x_max,algo->y_min,algo->y_max,algo->z_min,algo->z_max,6,6,6,true,true,true,8);
		
		for(int j = 0; j < particles; j++) {
			x=algo->reinit(algo->x_min, algo->x_max);
			y=algo->reinit(algo->y_min, algo->y_max);
			z=algo->reinit(algo->z_min, algo->z_max);
			con[i]->put(j,x,y,z);
		}
		
	}
	real_sizes = algo->compute_cell_sizes(con);

	//iterate until reach max iterations or precision
	while (true){
		int min_penalty_index = -1;
		min_penalty = 100;
		for (int i = 0; i < population_size; i++){
			penalty[i] = algo->size_penalty(real_sizes[i]);
			common_pool[penalty[i]] = con[i];
			if (penalty[i] == 0.0)
			{
				for (int m = 0; m < particles; m++)
					std::cout << "strange size " << real_sizes[i][m] << "\n";
			}
		
			if (min_penalty > penalty[i]) {
				min_penalty = penalty[i];
				min_penalty_index = i;
			}
		}
		
		if (iterations == 0){
			algo->output_data((filename + "/dist_first.txt").c_str(), real_sizes[min_penalty_index]);
			std::cout << (filename + "/dist_first.txt\n");
		}

		std::cout << "penalty " << min_penalty << " ===========================================\n";

		if (min_penalty < max_allowed_penalty) {
			algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
	
			break;
		}
		iterations++;
		std::cout << "iter " << iterations << "\n";
		
		if (iterations > max_iterations) {
			algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
	
			break;
		}
		else algo->output_data((filename + "/dist_tmp.txt").c_str(), real_sizes[min_penalty_index]);
		
		
		std::cout << "crossover\n";
		// con = crossover_by_mapping(con, real_sizes, -1);
		offspring = algo->crossover_by_mapping(con, real_sizes, 2);
		// con = crossover_by_mapping(con, real_sizes);
		
		std::cout << "mutation\n";
		algo->mutation(&offspring);
		
		double **real_sizes_offspring = algo->compute_cell_sizes(offspring);

		//common pool
		for (int i = 0; i < population_size; i++)
			common_pool[algo->size_penalty(real_sizes_offspring[i])] = offspring[i];
		
		int selected_offspring = 0;
		for (std::map<double, container*>::iterator it = common_pool.begin(); 
			it != common_pool.end() && selected_offspring < population_size; it++, selected_offspring++)
			con[selected_offspring] = it->second;
		
		real_sizes = algo->compute_cell_sizes(con);
	}
}
