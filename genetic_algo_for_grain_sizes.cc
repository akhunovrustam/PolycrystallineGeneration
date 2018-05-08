// Polycrystal rebuild from the size distribution using genetics algorithm
//
// Author   : Rustam Akhunov (GUT Gdansk)
// Email    : akhunovrustam@mail.ru
// Date     : 

#include "voro++.hh"
#include "ExtraStructAndFunc.hh"
#include "GeneticAlgoForSizesClass.cc"
#include <thread> 
#include <chrono>

using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds
	
int main(int argc, char *argv[]) {
	string prefix = "";
	if (argc > 1)
		prefix = argv[1] + prefix + "_exp";
	
	// cout << "prefix: " << prefix << "\n";
	// exit(0);
	
	std::stringstream ss;
	
	//create folder to save output data
	std::time_t t = std::time(0);
	std::tm* now = std::localtime(&t);
    ss << now->tm_mday << "-" << (now->tm_mon + 1) << "-" << (now->tm_year + 1900) << "_" << (now->tm_hour) << "." << (now->tm_min);
	
	std::string filename = "results/" + prefix + ss.str();
	system(("mkdir " + filename).c_str());
	
	//initial variables etc
	srand (time(NULL) * atoi(argv[1]));
	container **con;
	container **offspring;
	std::map<double, container*> common_pool;
	con = new container*[population_size]; // created con
	double x, y, z;
	double **real_sizes;
	double penalty[population_size];
	double min_penalty = 100;
	int iterations = 0;
	
	GeneticAlgoForSizesClass *algo = new GeneticAlgoForSizesClass(); // created algo
	
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
	real_sizes = algo->compute_cell_sizes(con); // created real_sizes

	int penalties = 0;
	//iterate until reach max iterations or precision
	while (true){
		
		int min_penalty_index = -1;
		min_penalty = 100;
		for (int i = 0; i < population_size; i++){
			penalty[i] = algo->size_penalty(real_sizes[i], 0, 0);
			common_pool[penalty[i] - rnd() * 0.0001] = con[i];
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
		
		penalties += population_size;
		algo->write_penalty_step(filename + "/penalty_steps.txt", penalties, min_penalty);
		
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
		offspring = algo->crossover_by_mapping(con, real_sizes, iterations, 2);
		// con = crossover_by_mapping(con, real_sizes);
		
		std::cout << "mutation\n";
		algo->mutation(&offspring);
		
		double **real_sizes_offspring = algo->compute_cell_sizes(offspring); // created real_sizes_offspring

		
		//common pool
		for (int i = 0; i < population_size; i++)
			common_pool[algo->size_penalty(real_sizes_offspring[i], 0, 0) - rnd() * 0.0001] = offspring[i];
		
		cout << common_pool.size();
		// exit(0);
		int selected_offspring = 0;
		for (std::map<double, container*>::iterator it = common_pool.begin(); 
			it != common_pool.end(); it++, selected_offspring++)
			if (selected_offspring < population_size)
				con[selected_offspring] = it->second;
			else {
				// cout << "delete " << selected_offspring << " \n";
				delete it->second;
				// it->second = nullptr;
				// cout << "delete end \n";
			}
			
		
		
		// cout << "offdelete 1 \n";
		// for (int i = 0; i < population_size; i++)
		// {
			// cout << "offdelete 100 \n";
			// if (offspring[i] != nullptr)
				// delete offspring[i];
			// offspring[i] = nullptr;
		// }
			// cout << "offdelete 2 \n";
		// for (int i = 0; i < population_size; i++)
			// delete offspring[i];
		delete [] offspring;
		
		common_pool.clear();
		for (int i = 0; i < population_size; i++)
			delete [] real_sizes[i];
		delete [] real_sizes; // deleted real_sizes
		
		for (int i = 0; i < population_size; i++)
			delete [] real_sizes_offspring[i];
		delete [] real_sizes_offspring; // deleted real_sizes_offspring
		
		real_sizes = algo->compute_cell_sizes(con);
	}
	
	delete [] con; // deleted con
	delete algo; // deleted algo
}
