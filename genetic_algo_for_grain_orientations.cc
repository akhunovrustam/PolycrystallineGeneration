// Polycrystal rebuild from the size distribution using genetics algorithm
//
// Author   : Rustam Akhunov (GUT Gdansk)
// Email    : akhunovrustam@mail.ru
// Date     : 

#include "voro++.hh"
#include "GeneticAlgoForOrientationsClass.cc"


int main() {
	stringstream ss;
	
	//create folder to save output data
	time_t t = time(0);
	tm* now = localtime(&t);
    ss << now->tm_mday << "-" << (now->tm_mon + 1) << "-" << (now->tm_year + 1900) << "_" << (now->tm_hour) << "." << (now->tm_min);
	
	string filename = "results/orient_" + ss.str();
	system(("mkdir " + filename).c_str());
	
	//initial variables etc
	srand (time(NULL));
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
	
	GeneticAlgoForOrientationsClass *algo = new GeneticAlgoForOrientationsClass();
	
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
			parents[i][signstr] = {alpha, beta, gamma};
		}
	}
		
	parents_rel = algo->relative_euler(con, parents);
	//iterate until reach max iterations or precision
	while (true){
		int min_penalty_index = -1;
		min_penalty = 1000000;
		for (int i = 0; i < population_size; i++){
			penalty[i] = algo->size_penalty(parents_rel[i]);
			cout << penalty[i] << "\n";
			// exit(0);
			common_pool[penalty[i]] = parents_rel[i];
			if (penalty[i] == 0.0)
			{
				// for (int m = 0; m < particles; m++)
					// cout << "strange size " << real_sizes[i][m] << "\n";
			}
		
			if (min_penalty > penalty[i]) {
				min_penalty = penalty[i];
				min_penalty_index = i;
			}
		}
		
		
		if (iterations == 0){
			algo->size_penalty(parents_rel[min_penalty_index], filename + "/dist_first.txt");
			// algo->output_data((filename + "/dist_first.txt").c_str(), real_sizes[min_penalty_index]);
			cout << (filename + "/dist_first.txt\n");
			// exit(0);
		}

		cout << "penalty " << min_penalty << " ===========================================\n";

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
		
		
		cout << "crossover\n";
		// con = crossover_by_mapping(con, real_sizes, -1);
		offspring = algo->crossover_by_mapping(con, parents, 2);
		// con = crossover_by_mapping(con, real_sizes);
		
		// cout << "mutation\n";
		// algo->mutation(&offspring);
		
		cout << "post genetic algo\n";
		offspring_rel = algo->relative_euler(con, offspring);

		//common pool
		for (int i = 0; i < population_size; i++)
			common_pool[algo->size_penalty(offspring_rel[i])] = offspring[i];
		
		int selected_offspring = 0;
		for (map<double, orient_unit>::iterator it = common_pool.begin(); 
			it != common_pool.end(); it++, selected_offspring++)
			if (selected_offspring < population_size)
				parents[selected_offspring] = it->second;
			else {}
		
		common_pool.clear();
		delete [] offspring;
		// delete [] parents;
		delete [] offspring_rel;
		delete [] parents_rel;
		parents_rel = algo->relative_euler(con, parents);
	}
	delete con;
	delete algo;
}
