// Polycrystal rebuild from the size distribution using genetics algorithm
//
// Author   : Rustam Akhunov (GUT Gdansk)
// Email    : akhunovrustam@mail.ru
// Date     : 

#include "voro++.hh"
#include "ExtraStructAndFunc.hh"
#include "GeneticAlgoForOrientationsClass.cc"
#include <thread> 
#include <chrono>

using namespace std;
using namespace this_thread;
using namespace chrono;

int main(int argc, char *argv[]) {
	time_point<system_clock> start, end, stop1, stop2;
    start = system_clock::now();
    
	//parse app arguments
	// 5th - multiplication, 6th - from, 7th - to
	string prefix = "1p_uniform_0.8_10_1.0_10_100";
	config cfg;
	int population_size;
	parse_args(argc, argv, &prefix, &cfg, &population_size);
	
	//create dir for data
	string filename;
	create_dir(&filename, prefix, "results_orientation/");
	
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
	// exit(0);
	double koef = 1 / (double)parents_rel[0].size();
	int df = parents_rel[0].size() - 3000;
	double k = 1 / pow(penalty_steps, 3);
	double sum = df*pow(4*koef - k, 2);
	int df2 = 1000 - df;
	sum += df2*pow(3*koef - k, 2);
	sum /= parents_rel[0].size();
	cout << "min pen: " << sum << endl;
	// exit(0);
	
	//iterate until reach max iterations or precision
	while (true){
		cout << "begin iter " << iterations << " ========================\n";
		
		int min_penalty_index = -1;
		int max_penalty_index = -1;
		min_penalty = 1000000;
		double max_penalty = 0;
		double avg_penalty = 0;
		
		for (int i = 0; i < population_size; i++){
			penalty[i] = algo->size_penalty(parents_rel[i]);
			avg_penalty += penalty[i];
			// cout << penalty[i] << "\n";
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

		cout << "penalty: " << min_penalty << " \n";

		// exit(0);
		if (min_penalty < max_allowed_penalty) {
			// algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
	
			break;
		}
		
		if (iterations > max_iterations) {
			// algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
			algo->size_penalty(parents_rel[min_penalty_index], filename + "/dist_end.txt");
	
			break;
		}
		else algo->size_penalty(parents_rel[min_penalty_index], filename + "/dist_tmp.txt");
			
		
		cout << "crossover\n";
		
		// con = crossover_by_mapping(con, real_sizes, -1);
		offspring = algo->crossover_by_mapping(con, parents, parents_rel, cfg.co);
	
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
		
		cout << "post penalty recalc\n";
		int selected_offspring = 0;
	
		for (map<double, orient_unit>::iterator it = common_pool.begin(); 
			it != common_pool.end(); it++, selected_offspring++)
			if (selected_offspring < population_size)
				parents[selected_offspring] = it->second;
			else {}
		
		// exit(0);
		iterations++;
		
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
