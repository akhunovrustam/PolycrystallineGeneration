// Polycrystal rebuild from the size distribution using genetics algorithm
//
// Author   : Rustam Akhunov (GUT Gdansk)
// Email    : akhunovrustam@mail.ru
// Date     : 

#include "voro++.hh"
#include "ExtraStructAndFunc.hh"
#include "GeneticAlgoForOrientationsClass.cc"


int main(int argc, char *argv[]) {
	string prefix = "2p_uniform_0.1_20";
	if (argc == 3)
	{
		prefix = argv[2];
		prefix += "_SEED_";
		prefix += argv[1];
		prefix += "_exp";
	}
	else if (argc == 2)
	{
		prefix = ((prefix + "_SEED_") + argv[1]) + "_exp";
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
	
	std::string filename = "results_orientation/orient_" + prefix + ss.str();
	system(("mkdir " + filename).c_str());
	
	//initial variables etc
	srand (time(NULL) * atoi(argv[1]));
	container *con;
	orient_unit parents[population_size];
	orient_unit *offspring;
	orient_unit *parents_rel;
	orient_unit *offspring_rel;
	map<double, orient_unit> common_pool;
	double x, y, z;
	double penalty[population_size];
	double min_penalty = 100;
	int iterations = 0;
	
	GeneticAlgoForOrientationsClass *algo = new GeneticAlgoForOrientationsClass(population_size);
	
	con = new container(algo->x_min, algo->x_max, algo->y_min, algo->y_max, algo->z_min, algo->z_max, 6, 6, 6, true, true, true, 8);
	for(int j = 0; j < particles; j++) {
		x = algo->reinit(algo->x_min, algo->x_max);
		y = algo->reinit(algo->y_min, algo->y_max);
		z = algo->reinit(algo->z_min, algo->z_max);
		con->put(j, x, y, z);
		stringstream sign;
		sign << x << "|" << y << "|" << z;
		string signstr = sign.str();
		for (int i = 0; i < population_size; i++){
			double alpha = algo->reinit_angles();
			double beta = algo->reinit_angles();
			double gamma = algo->reinit_angles();
			// if (alpha < 0 || beta < 0 || gamma < 0)
				// cout << "bad angles: " << alpha << " " << beta << " " << gamma << "\n";
			parents[i][signstr] = {alpha, beta, gamma};
		}
	}
	// exit(0);
		
	int penalties = 0;
	parents_rel = algo->relative_euler(con, parents);
	//iterate until reach max iterations or precision
	while (true){
		cout << "begin ========================================\n";
		int min_penalty_index = -1;
		int max_penalty_index = -1;
		min_penalty = 1000000;
		double max_penalty = 0;
		double avg_penalty = 0;
		
		for (int i = 0; i < population_size; i++){
			penalty[i] = algo->size_penalty(parents_rel[i]);
			avg_penalty += penalty[i];
			cout << penalty[i] << "\n";
			// exit(0);
			common_pool[penalty[i]] = parents[i];
			if (penalty[i] == 0.0)
			{
				// for (int m = 0; m < particles; m++)
					// cout << "strange size " << real_sizes[i][m] << "\n";
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
			algo->size_penalty(parents_rel[min_penalty_index], filename + "/dist_first.txt");
			// algo->output_data((filename + "/dist_first.txt").c_str(), real_sizes[min_penalty_index]);
			cout << (filename + "/dist_first.txt\n");
			// exit(0);
		}

		cout << "penalty " << min_penalty << " !!!!!!\n";

		if (min_penalty < max_allowed_penalty) {
			// algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
	
			break;
		}
		iterations++;
		cout << "iter " << iterations << "\n";
		
		if (iterations > max_iterations) {
			// algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
			algo->size_penalty(parents_rel[min_penalty_index], filename + "/dist_end.txt");
	
			break;
		}
		else algo->size_penalty(parents_rel[min_penalty_index], filename + "/dist_tmp.txt");
			
		// else algo->output_data((filename + "/dist_tmp.txt").c_str(), real_sizes[min_penalty_index]);
		
		// if (iterations == 2)
		// {
			// cout << parents[0].size() << "\n";
			// exit(0);
			// for (orient_unit::iterator it = parents[0].begin(); it != parents[0].end(); it++)
				// cout << it->second.alpha << " " << it->second.beta << " " << it->second.gamma << "\n";
		// }
		
		cout << "crossover\n";
		
		// con = crossover_by_mapping(con, real_sizes, -1);
		offspring = algo->crossover_by_mapping(con, parents, parents_rel, cfg.co);
		// con = crossover_by_mapping(con, real_sizes);
		
			// exit(0);
		// if (iterations == 2)
		// {
			// for (orient_unit::iterator it = offspring[0].begin(); it != offspring[0].end(); it++)
				// cout << it->second.alpha << " " << it->second.beta << " " << it->second.gamma << "\n";
			// exit(0);
		// }
		cout << "mutation\n";
		algo->mutation(&offspring, cfg.md == 2 ? true : false, cfg.mp, cfg.md == 2 ? 0 : cfg.md);
		
		cout << "post genetic algo\n";
		offspring_rel = algo->relative_euler(con, offspring);

		
		cout << "post relative cacl\n";
		//common pool
		// for (map<double, orient_unit>::iterator it = common_pool.begin(); it != common_pool.end(); it++)
			// cout << it->second.size() << "\n";
		
		for (int i = 0; i < population_size; i++)
		{
			double pen = algo->size_penalty(offspring_rel[i]);
			common_pool[pen] = offspring[i];
			cout << pen << "\n";
		}
		
		// cout << parents[0].size() << "\n";
		// cout << parents[1].size() << "\n";
		// cout << offspring[0].size() << "\n";
		// cout << offspring[1].size() << "\n";
		
		cout << "post penalty recalc\n";
		int selected_offspring = 0;
		// for (map<double, orient_unit>::iterator it = common_pool.begin(); it != common_pool.end(); it++)
			// cout << it->second.size() << "\n";
		// exit(0);
		for (map<double, orient_unit>::iterator it = common_pool.begin(); 
			it != common_pool.end(); it++, selected_offspring++)
			if (selected_offspring < population_size)
				parents[selected_offspring] = it->second;
			else {}
		
		// exit(0);
		
		common_pool.clear();
		delete [] offspring;
		// delete [] parents;
		delete [] offspring_rel;
		delete [] parents_rel;
		parents_rel = algo->relative_euler(con, parents);
		
		cout << "end ========================================\n";
	}
	delete con;
	delete algo;
}
