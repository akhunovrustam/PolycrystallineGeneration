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

int main(int argc, char *argv[]) {
	string prefix = "1p_uniform_0.8_10_exp";
	if (argc == 3)
	{
		prefix = argv[2];
		prefix += "_SEED_";
		prefix += argv[1];
		prefix += "_exp";
	}
	
	// exit(0);
	config cfg = parse_prefix(prefix);
	int population_size = population_size_const;
	if (prefix != ""){
		population_size = cfg.ps;
	}
	
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
	// cout << cfg.co << " " << cfg.md << " " << cfg.mp << " " << cfg.ps << "\n";
	container **offspring;
	std::map<double, container*> common_pool;
	// exit(0);
	con = new container*[population_size]; // created con
	double x, y, z;
	double **real_sizes;
	double penalty[population_size];
	double min_penalty = 100;
	int iterations = 0;
	
	cout << "prefix: " << prefix << "\n";
	GeneticAlgoForSizesClass *algo = new GeneticAlgoForSizesClass(population_size); // created algo
	
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
		int max_penalty_index = -1;
		min_penalty = 1000000;
		double max_penalty = 0;
		double avg_penalty = 0;
		
		for (int i = 0; i < population_size; i++){
			penalty[i] = algo->size_penalty(real_sizes[i], 0, 0);
			avg_penalty += penalty[i];
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
			
			if (max_penalty < penalty[i]) {
				max_penalty = penalty[i];
				max_penalty_index = i;
			}
		}
		avg_penalty /= population_size;
		
		penalties += population_size;
		algo->write_penalty_step(filename + "/penalty_steps_best.txt", penalties, min_penalty);
		algo->write_penalty_step(filename + "/penalty_steps_worst.txt", penalties, max_penalty);
		algo->write_penalty_step(filename + "/penalty_steps_avg.txt", penalties, avg_penalty);
		
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
		// offspring = con;
		offspring = algo->crossover_by_mapping(con, real_sizes, iterations, cfg.co);
		// con = crossover_by_mapping(con, real_sizes);
		
		std::cout << "mutation\n";
		algo->mutation(&offspring, cfg.md == 2 ? true : false, cfg.mp, cfg.md == 2 ? 0 : cfg.md);
		
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
