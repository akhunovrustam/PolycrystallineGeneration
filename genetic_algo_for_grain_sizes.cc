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

using namespace std;
using namespace this_thread;
using namespace chrono;

int main(int argc, char *argv[]) {
	time_point<system_clock> start, end, stop1, stop2;
    start = system_clock::now();
    
	//parse app arguments
	// 5th - multiplication for mutation amplitude, 6th - from, 7th - to, 8th - 1st threshold, 9th - 2nd threshold
	// 10th - factor for mp, 11th - 1st threshold factor for population size, 12th - 2nd threshold factor for population size
	string prefix = "2p_uniform_0.8_10_1.0_10_100_500_1000_10_5_9";
	config cfg;
	int population_size;
	parse_args(argc, argv, &prefix, &cfg, &population_size);
	
	//create dir for data
	string filename;
	create_dir(&filename, prefix);
	
	//initialize variables etc
	srand (time(NULL) * atoi(argv[1]));
	container **con;
	container **offspring;
	map<double, con_and_points> common_pool;
	con = new container*[population_size];
	double x, y, z;
	double **real_sizes;
	double penalty[population_size];
	double min_penalty = 100;
	int iterations = 0;
	
	cout << "prefix: " << prefix << "\n";
	GeneticAlgoForSizesClass *algo = new GeneticAlgoForSizesClass(population_size); // created algo
	
	sorted_points* experiment;
	experiment = new sorted_points[population_size + 1];
	
	sorted_points* experiment_off;
	experiment_off = new sorted_points[population_size + 1];
	
	neighbors* neg;
	neg = new neighbors[population_size];
	
	//initial population
	for (int i = 0; i < population_size; i++){
		con[i] = new container(algo->x_min,algo->x_max,algo->y_min,algo->y_max,algo->z_min,algo->z_max,6,6,6,true,true,true,8);
		
		for(int j = 0; j < particles; j++) {
			x=algo->reinit(algo->x_min, algo->x_max);
			y=algo->reinit(algo->y_min, algo->y_max);
			z=algo->reinit(algo->z_min, algo->z_max);
			experiment[i][x][y][z] = j;
			con[i]->put(j,x,y,z);
		}
		
	}
	real_sizes = algo->compute_cell_sizes(con, &neg); // calculate cell sizes
	
	stop1 = system_clock::now();
	int penalties = 0;
	
	int elapsed_time1 = duration_cast<milliseconds>(stop1-start).count();
	
	cout << "initialization time: " << elapsed_time1 << " ms" << endl;
	// exit(0);
	// exit(0);
	//iterate until reach max iterations or precision
	int total_time = 0;
	double mprob = cfg.mp;
	while (true){
		stop1 = system_clock::now();
		cout << "begin iter " << iterations << " ========================\n";
		
		if (iterations > cfg.thr1) {
			cfg.ps_fac = cfg.ps_fac1;
			cfg.mp = mprob/cfg.fac;
		}
		if (iterations > cfg.thr2) {
			cfg.ps_fac = cfg.ps_fac2;
			cfg.mp = mprob/cfg.fac/cfg.fac;
		}
		
		int min_penalty_index = -1;
		int max_penalty_index = -1;
		min_penalty = 1000000;
		double max_penalty = 0;
		double avg_penalty = 0;
		
		for (int i = 0; i < population_size; i++){
			penalty[i] = algo->size_penalty(real_sizes[i], 0, 0);
			avg_penalty += penalty[i];
			con_and_points cn;
			cn.con = con[i];
			cn.pt = experiment[i];
			
			common_pool[penalty[i] - rnd() * 0.0001] = cn;
			if (penalty[i] == 0.0)
			{
				for (int m = 0; m < particles; m++)
					cout << "strange size " << real_sizes[i][m] << "\n";
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
		
		cout << "penalty " << min_penalty << " ===========================================\n";

		penalties += population_size;
		algo->write_penalty_step(filename + "/penalty_steps_best.txt", penalties, min_penalty);
		algo->write_penalty_step(filename + "/penalty_steps_worst.txt", penalties, max_penalty);
		algo->write_penalty_step(filename + "/penalty_steps_avg.txt", penalties, avg_penalty);
		
		//write data into files
		if (iterations == 0){
			algo->output_data((filename + "/dist_first.txt").c_str(), real_sizes[min_penalty_index]);
			cout << (filename + "/dist_first.txt\n");
		}

		
		if (min_penalty < max_allowed_penalty) {
			algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
	
			break;
		}
		
		cout << "crossover\n";
		
		// con = crossover_by_mapping(con, real_sizes, -1);
		// offspring = con;
		stop2 = system_clock::now();
		offspring = algo->crossover_by_mapping(con, penalty, iterations, cfg.co, &experiment, &experiment_off, cfg, neg);

		
		// exit(0);
		cout << "mutation\n";
		algo->mutation(&offspring, iterations, cfg.md == 2 ? true : false, cfg.mp, cfg.md == 2 ? 0 : cfg.md, &experiment_off, cfg.mult);
		
		cout << "offspring cell size calc\n";
		double **real_sizes_offspring = algo->compute_cell_sizes(offspring, &neg); // created real_sizes_offspring

		
		cout << "offspring to common pool\n";
		//common pool
		for (int i = 0; i < population_size; i++)
		{
			con_and_points cn;
			cn.con = offspring[i];
			cn.pt = experiment_off[i];
			common_pool[algo->size_penalty(real_sizes_offspring[i], 0, 0) - rnd() * 0.00000001] = cn;
		}
		// exit(0);
		cout << "select best\n";
		int selected_offspring = 0;
		
		for (map<double, con_and_points>::iterator it = common_pool.begin(); it != common_pool.end(); it++, selected_offspring++)
			if (selected_offspring < population_size)
			{
				con[selected_offspring] = it->second.con;
				experiment[selected_offspring] = it->second.pt;
			}
			else 
			{
				// cout << "delete " << selected_offspring << " \n";
				delete it->second.con;
				// it->second = nullptr;
				// cout << "delete end \n";
			}
		
		if (iterations % 100 == 0)
		{
			stringstream image;
	
			image << "distributions_for_iteration_" << iterations << ".png";
	
			system(("gnuplot -e \"filename='" + image.str() + "'; prefix='" + filename + "'; penalty='" + to_string(min_penalty) + "'\" plot_first_tmp_auto.gp").c_str());
		}
		
		
		iterations++;
		
		if (penalties > max_iterations) {
			algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
	
			break;
		}
		else algo->output_data((filename + "/dist_tmp.txt").c_str(), real_sizes[min_penalty_index]);
		
		
		cout << "release all\n";
		//release all memory 
		delete [] offspring;
		
		common_pool.clear();
		for (int i = 0; i < population_size; i++)
			delete [] real_sizes[i];
		delete [] real_sizes; // deleted real_sizes
		
		for (int i = 0; i < population_size; i++)
			delete [] real_sizes_offspring[i];
		delete [] real_sizes_offspring; // deleted real_sizes_offspring
		
		//recompute sizes
		real_sizes = algo->compute_cell_sizes(con, &neg);

		stop2 = system_clock::now();
		elapsed_time1 = duration_cast<milliseconds>(stop2-stop1).count();
		
		cout << endl;
		cout << "iteration time: " << elapsed_time1 << " ms" << endl;
		
		total_time += elapsed_time1;
		cout << "avg iteration time: " << total_time/iterations << " ms" << endl;
		cout << endl;
		
	}
	
	//cout << "command " << ("convert -delay 300 -loop 0 " + filename + "/*.eps progress.gif").c_str() << endl;
	system(("convert -delay 100 -loop 0 " + filename + "/*.png " + filename + "/progress.gif").c_str());
	
	delete [] con; // deleted con
	delete algo; // deleted algo

	end = system_clock::now();
	elapsed_time1 = duration_cast<milliseconds>(end-start).count();
	cout << "whole time: " << elapsed_time1 << " ms" << endl;
}
