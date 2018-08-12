#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include "voro++.hh"
#include "GeneralClass.hh"
#include <thread> 
#include <chrono>
#include "lammps/src/library.h"
#include "lammps/src/lammps.h"
#include <iostream>
#include "voro++.hh"
#include "ExtraStructAndFunc.hh"
#include "GeneticAlgoForSizesClass.cc"
#include "GeneticAlgoForOrientationsClass.cc"
#include "GeneticAlgoForShiftingClass.cc"
#include "LatticeGeneratorClass.cc"

using namespace LAMMPS_NS;
using namespace std;
using namespace this_thread;
using namespace chrono;

General1Class::General1Class(){

}

void General1Class::GrainSizeGenerate(int argc, char *argv[]){
	time_point<system_clock> start, end, stop1, stop2;
    start = system_clock::now();
    
	//parse app arguments
	// 5th - multiplication for mutation amplitude, 6th - from, 7th - to, 8th - 1st threshold, 9th - 2nd threshold
	// 10th - factor for mp, 11th - 1st threshold factor for population size, 12th - 2nd threshold factor for population size
	string prefix = "2p_uniform_0.8_30_1.0_10_100_200_500_10_1_1";
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
		
		cout << "change configs: " << cfg.ps_fac << " " << cfg.mp << endl;
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
	
			best_cont = con[min_penalty_index];
			best_cont_points = experiment[min_penalty_index];
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
			cout << "image:! " << image.str() << " " << filename << endl;
			system(("gnuplot -e \"filename='" + image.str() + "'; prefix='" + filename + "'; penalty='" + to_string(min_penalty) + "'\" plot_first_tmp_auto.gp").c_str());
		}
		
		
		iterations++;
		
		if (penalties > max_iterations) {
			algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
			
			best_cont = con[min_penalty_index];
			best_cont_points = experiment[min_penalty_index];
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

void General1Class::GrainOrientationGenerate(int argc, char *argv[]){
	time_point<system_clock> start, end, stop1, stop2;
    start = system_clock::now();
    
	//parse app arguments
	// 5th - multiplication, 6th - from, 7th - to
	string prefix = "1p_uniform_0.9_10_1.0_10_100_500_1000_10_5_9";
	config cfg;
	int population_size;
	parse_args(argc, argv, &prefix, &cfg, &population_size);
	
	//create dir for data
	string filename;
	create_dir(&filename, prefix, "results_orientation/");
	
	//initial variables etc
	srand (time(NULL) * atoi(argv[1]));
	container* con = best_cont;
	orient_unit parents[population_size];
	orient_unit *offspring;
	
	orient_unit *offspring_rel;
	map<double, orient_unit> common_pool;
	double x, y, z;
	double penalty[population_size];
	double min_penalty = 100;
	int iterations = 0;
	
	GeneticAlgoForOrientationsClass *algo = new GeneticAlgoForOrientationsClass(population_size);
	
	for (sorted_points::iterator it = best_cont_points.begin(); it != best_cont_points.end(); it++)
	{
		for (map<double, map<double, int>>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
			for (map<double, int>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++)
			{
				double x = it->first;
				double y = it2->first;
				double z = it3->first;
				stringstream sign;
				sign << x << "|" << y << "|" << z;
				string signstr = sign.str();
				for (int i = 0; i < population_size; i++){
					// double alpha = rnd()*0.15;
					double alpha = algo->reinit_angles();
					// double beta = rnd()*0.15;
					double beta = algo->reinit_angles();
					// double gamma = rnd()*0.15;
					double gamma = algo->reinit_angles();
					// if (alpha < 0 || beta < 0 || gamma < 0)
						// cout << "bad angles: " << alpha << " " << beta << " " << gamma << "\n";
					parents[i][signstr] = {alpha, beta, gamma};
				}
			}
	}
	
	// exit(0);
		
	int penalties = 0;
	orient_unit *parents_rel = algo->relative_euler(con, parents);
	
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
			algo->output_data((filename + "/dist_first.txt").c_str(), parents_rel[min_penalty_index]);
			cout << (filename + "/dist_first.txt\n");
			// exit(0);
		}

		cout << "penalty: " << min_penalty << " \n";

		// exit(0);
		if (min_penalty < max_allowed_penalty_orient) {
			// algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
	
			best_orient = parents[min_penalty_index];
			break;
		}
		
		if (penalties > max_iterations_orient) {
			// algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
			algo->size_penalty(parents_rel[min_penalty_index], filename + "/dist_end.txt");
			
			best_orient = parents[min_penalty_index];
			break;
		}
		else 
		{
			// algo->size_penalty(parents_rel[min_penalty_index], filename + "/dist_tmp.txt");
			algo->output_data((filename + "/dist_tmp.txt").c_str(), parents_rel[min_penalty_index]);
			
		}
			
		
		cout << "crossover\n";
		
		// con = crossover_by_mapping(con, real_sizes, -1);
		offspring = algo->crossover_by_mapping(con, parents, parents_rel, cfg.co);
		
		cout << "pr angle: " << parents[0].begin()->second.alpha << endl;
		cout << "fr angle: " << offspring[0].begin()->second.alpha << endl;
		
		cout << "mutation\n";
		algo->mutation(&offspring, cfg.md == 2 ? true : false, cfg.mp, cfg.md == 2 ? 0 : cfg.md);
		cout << "fr angle2: " << offspring[0].begin()->second.alpha << endl;
		
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
			// cout << pen << "\n";
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
	// delete con;
	delete algo;
}

void General1Class::FillGrains(){
	LatticeGeneratorClass gen(best_cont);
	at = gen.fillRandomlyAndBuildGrains(best_orient);
}

string General1Class::RelaxFirst(int argc, char *argv[]){
	//test data
	// at[0][0] = {0, 0, 0};
	// at[0][1] = {0, 0, 4.095};
	// at[0][2] = {0, 4.095, 0};
	// at[0][3] = {4.095, 0, 0};
	
	// at[1][0] = {4.095, 4.095, 4.095};
	// at[1][1] = {4.095, 4.095, 8.19};
	// at[1][2] = {4.095, 8.19, 4.095};
	// at[1][3] = {8.19, 4.095, 4.095};
	
	// best_cont = new container(0 , half_boxside, 0, half_boxside, 0, half_boxside, 6, 6, 6, true, true, true, 8);
	// best_cont->put(0, 2.0, 2.0, 2.0);
	// best_cont->put(1, 6.0, 6.0, 6.0);
	
	//parse app arguments
	// 5th - multiplication for mutation amplitude, 6th - from, 7th - to, 8th - 1st threshold, 9th - 2nd threshold
	// 10th - factor for mp, 11th - 1st threshold factor for population size, 12th - 2nd threshold factor for population size
	string prefix = "2p_uniform_0.8_20_1.0_10_100_500_1000_10_5_9";
	config cfg;
	int population_size;
	parse_args(argc, argv, &prefix, &cfg, &population_size);
	
	//create dir for data
	string filename;
	create_dir(&filename, prefix, "results_relax1/");
	
	//initialize variables etc
	srand (time(NULL) * atoi(argv[1]));
	
	map<int, atoms>* atms = new map<int, atoms>[population_size_relax];
	map<int, atoms>* offspring;
	
	atms[0] = at;
	
	GeneticAlgoForShiftingClass *algo = new GeneticAlgoForShiftingClass(population_size);
	
	algo->pe(filename + "/..", at);
	
	for (int i = 1; i < population_size_relax; i++)
	{
		for (map<int, atoms>::iterator j = at.begin(); j != at.end(); j++)
		{
			double dx = algo->ShiftDelta(0);
			double dy = algo->ShiftDelta(0);
			double dz = algo->ShiftDelta(0);
			for (atoms::iterator it = j->second.begin(); it != j->second.end(); it++)
			{
				atom_coords ac = {it->second.x + dx, it->second.y + dy, it->second.z + dz};
				atms[i][j->first][it->first] = ac;
			}
		}
		
	}

	
	map<double, map<int, atoms>> common_pool;
	double penalty[population_size];
	double min_penalty = 100;
	int iterations = 0;
	
	cout << "prefix: " << prefix << "\n";
	
	
	
	voronoicell_neighbor c;
	vector<double> v;
	neighbors neg;
	vector<int> neigh;
	double x,y,z;
	int id;
	
	c_loop_all cl(*(best_cont));
	int loop_counter = 0;
	if(cl.start()) do if(best_cont->compute_cell(c,cl)) {
		// printf("computed individual %i\n", i);
		cl.pos(x,y,z);
		id=cl.pid();
	
		// Gather information about the computed Voronoi cell
		c.vertices(x,y,z,v);
		c.neighbors(neigh);
		
		neg[id] = neigh;
		
	} while (cl.inc());
	
	double mprob = cfg.mp;
	int penalties = 0;
	while (true){
		cout << "begin iter " << iterations << " ========================\n";
		
		// if (iterations > cfg.thr1) {
			// cfg.ps_fac = cfg.ps_fac1;
			// cfg.mp = mprob/cfg.fac;
		// }
		// if (iterations > cfg.thr2) {
			// cfg.ps_fac = cfg.ps_fac2;
			// cfg.mp = mprob/cfg.fac/cfg.fac;
		// }
		
		int min_penalty_index = -1;
		int max_penalty_index = -1;
		min_penalty = 10000000000;
		double max_penalty = 0;
		double avg_penalty = 0;
		
		for (int i = 0; i < population_size; i++){
			penalty[i] = algo->size_penalty(atms[i]);
			avg_penalty += penalty[i];
			common_pool[penalty[i] - rnd() * 0.000000000001] = atms[i];
			
			if (penalty[i] == 0.0)
			{
				cout << "zero penalty" << endl;
				exit(0);
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
			// algo->writeAtoms(atms[min_penalty_index]);
			// algo->output_data((filename + "/dist_first.txt").c_str(), real_sizes[min_penalty_index]);
			cout << (filename + "/dist_first.txt\n");
		}

		ofstream myfile;
		myfile.open ("debug.etotal.txt", ofstream::out | ofstream::app);
		myfile << min_penalty << "\n";
		myfile.close();
		
		/*if (min_penalty < max_allowed_penalty) {
			// algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
	
			break;
		}*/
		
		cout << "crossover\n";
		
		offspring = algo->crossover_by_mapping(penalty, atms, cfg, neg, best_cont_points);

		
		cout << "mutation\n";
		algo->mutation(&offspring, cfg);
		
		cout << "offspring cell size calc\n";
		
		
		cout << "offspring to common pool\n";
		//common pool
		for (int i = 0; i < population_size; i++)
		{
			common_pool[algo->size_penalty(offspring[i]) - rnd() * 0.00000001] = offspring[i];
		}
		// exit(0);
		cout << "select best\n";
		int selected_offspring = 0;
		
		for (map<double, map<int, atoms>>::iterator it = common_pool.begin(); it != common_pool.end(); it++, selected_offspring++)
			if (selected_offspring < population_size)
			{
				atms[selected_offspring] = it->second;
			}
			else 
			{
				it->second.clear();
			}
		
		// if (iterations % 100 == 0)
		// {
			// stringstream image;
	
			// image << "distributions_for_iteration_" << iterations << ".png";
	
			// system(("gnuplot -e \"filename='" + image.str() + "'; prefix='" + filename + "'; penalty='" + to_string(min_penalty) + "'\" plot_first_tmp_auto.gp").c_str());
		// }
		
		
		iterations++;
		
		if (penalties > max_iterations_relax) {
			// algo->output_data((filename + "/dist_end.txt").c_str(), real_sizes[min_penalty_index]);
	
			break;
		}
		// else algo->output_data((filename + "/dist_tmp.txt").c_str(), real_sizes[min_penalty_index]);
		
		
		cout << "release all\n";
		//release all memory 
		delete [] offspring;
		
		common_pool.clear();
		
	}
	cout << "peeeeee" << endl;
	algo->pe(filename, atms[0]);
	cout << "after peeeeee" << endl;
	//cout << "command " << ("convert -delay 300 -loop 0 " + filename + "/*.eps progress.gif").c_str() << endl;
	// system(("convert -delay 100 -loop 0 " + filename + "/*.png " + filename + "/progress.gif").c_str());
	
	delete [] atms; // deleted con
	delete algo; // deleted algo
	return filename;
}

void General1Class::RelaxSecond(int argc, char *argv[], string filename){
	
	cout << "relax second 1" << endl;
	char **lmparg = new char*[8];
	lmparg[0] = NULL;
	lmparg[1] = (char *) "-e";
	lmparg[2] = (char *) "none";
	lmparg[3] = (char *) "-screen";
	lmparg[4] = (char *) "none";
	lmparg[5] = (char *) "-v";
	lmparg[6] = (char *) "f";
	lmparg[7] = (char *) filename.c_str();
	
	
	LAMMPS *lmp;
	lammps_open_no_mpi(8, lmparg, (void **) &lmp);
	lammps_file(lmp, (char *) "etot.lmp");
	double tot;
	double tot1 = 0;
	tot = lammps_get_thermo(lmp, (char *) "etotal");
	
	lammps_close(lmp);
	
	cout << "relax second 2" << endl;
	
	double petot = 0;
	int clear = 0;
	system(("cp " + filename + "/coords_cn.dump " + filename + "/coords_cn_old.dump").c_str());
	int steps = 0;
	bool full_stop = false;
	while (true)
	{
		steps++;
		ifstream dumpfile;
		if (tot > tot1) {
			tot = tot1;
			ofstream myfile;
			myfile.open (filename + "/debug.relax2.txt", ofstream::out | ofstream::app);
			myfile << "success" << endl;
			myfile.close();
			system(("cp " + filename + "/coords_cn_new.dump " + filename + "/coords_cn_old.dump").c_str());
			clear = 0;
		}
		else
		{
			ofstream myfile;
			myfile.open (filename + "/debug.relax2.txt", ofstream::out | ofstream::app);
			myfile << "fail" << endl;
			myfile.close();
			
		}
		dumpfile.open(filename + "/coords_cn_old.dump");
		
		if (dumpfile.is_open())
		{
			string line, original;
			int linenum = 0;
			bool removed = false;
			ofstream myfile, dumpfile_new;
			myfile.open(filename + "/fcc_lattice_final.data");
			// dumpfile_new.open("coords_cn_new.dump");
			int atom_amount = 0;
			while ( getline (dumpfile, line) )
			{
				original = line;
				if (linenum == 3)
				{
					atom_amount = (stoi(line)-1);
					myfile << "Oriented crystal file" << endl;
					myfile << endl;
					myfile << atom_amount << " atoms" << endl;
					myfile << "0 bonds" << endl;
					myfile << "0 angles" << endl;
					myfile << "0 dihedrals" << endl;
					myfile << "0 impropers" << endl;
					myfile << endl;
					// myfile << atm.size() << " atom types" << endl;
					myfile << "1 atom types" << endl;
					myfile << endl;
					myfile << "0 " << half_boxside << " xlo xhi" << endl;
					myfile << "0 " << half_boxside << " ylo yhi" << endl;
					myfile << "0 " << half_boxside << " zlo zhi" << endl;
					myfile << endl;
					myfile << "Atoms" << endl;
					myfile << endl;
				
					// dumpfile_new << (stoi(line)-1) << endl;
					linenum++;
					continue;
				}
				
				if (linenum > 8)
				{
				
					string delimiter = " ";

					if (line == "") continue;
					size_t pos = 0;
					string token;
					string data[8];
					int cnt = 0;
					
					while ((pos = line.find(delimiter)) != string::npos) {
						token = line.substr(0, pos);
						// if (token == " " || token == "") continue;
						data[cnt++] = token;
						line.erase(0, pos + delimiter.length());
					}
					data[cnt] = line;
					
					if (stod(data[5]) > relax_threshold && !removed && clear < linenum)
					{
						ofstream myfile1;
						myfile1.open (filename + "/debug.relax2.txt", ofstream::out | ofstream::app);
						myfile1 << "linenum: " << linenum << endl;
						myfile1.close();
						clear = linenum;
						removed = true;
						continue;
					}
					if (!removed && linenum - 9 == atom_amount) full_stop = true;
					myfile << data[0] << " " << data[1] << " " << data[2] << " " << data[3] << " " << data[4] << endl;
					
				}
				// dumpfile_new << original << endl;
				linenum++;
			}
			if (full_stop) break;
			dumpfile.close();
			// dumpfile_new.close();
			cout << tot << " " << tot1 << endl;
	// exit(0);
		}
		cout << "bufer check 2222!!" << endl;
		char **lmparg = new char*[8];
		lmparg[0] = NULL;           
		lmparg[1] = (char *) "-e";
		lmparg[2] = (char *) "none";
		lmparg[3] = (char *) "-screen";
		lmparg[4] = (char *) "none";
		lmparg[5] = (char *) "-v";
		lmparg[6] = (char *) "f";
		lmparg[7] = (char *) filename.c_str();
		
		
		LAMMPS *lmp;
		lammps_open_no_mpi(8, lmparg, (void **) &lmp);
		lammps_file(lmp, (char *) "etot_final.lmp");
		
		tot1 = lammps_get_thermo(lmp, (char *) "etotal");
		cout << tot << " " << tot1 << endl;
		// exit(0);
		lammps_close(lmp);
		ofstream myfile;
		myfile.open (filename + "/debug.relax2.txt", ofstream::out | ofstream::app);
		myfile << tot1 << endl;
		myfile.close();
	}
}


void General1Class::writeAll(){
	ofstream myfile;
	myfile.open ("all_points.txt", ofstream::out | ofstream::app);
	
	for (sorted_points::iterator it = best_cont_points.begin(); it != best_cont_points.end(); it++)
		for (map<double, map<double, int>>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
			for (map<double, int>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++)
			{
				double x = it->first;
				double y = it2->first;
				double z = it3->first;
				int id = it3->second;
				myfile << id << " " << x << " " << y << " " << z << endl;
			}
	
	myfile.close();
	
	myfile.open ("all_angles.txt", ofstream::out | ofstream::app);
	
	for (orient_unit::iterator it = best_orient.begin(); it != best_orient.end(); it++)
	{
		string str = it->first;
		double alpha = it->second.alpha;
		double beta = it->second.beta;
		double gamma = it->second.gamma;
		myfile << str << " " << alpha << " " << beta << " " << gamma << endl;
	}
	
	myfile.close();
}

void General1Class::loadAll(){
	GeneticAlgoForSizesClass *algo = new GeneticAlgoForSizesClass(10);
	best_cont = new container(algo->x_min,algo->x_max,algo->y_min,algo->y_max,algo->z_min,algo->z_max,6,6,6,true,true,true,8);
	
	string line;
	ifstream myfile("all_points.txt");
	if (myfile.is_open())
	{
		while ( getline (myfile, line) )
		{
			string delimiter = " ";

			if (line == "") continue;
			size_t pos = 0;
			string token;
			string data[8];
			int cnt = 0;
			while ((pos = line.find(delimiter)) != string::npos) {
				token = line.substr(0, pos);
				// if (token == " " || token == "") continue;
				data[cnt++] = token;
				line.erase(0, pos + delimiter.length());
			}
			data[cnt] = line;
			best_cont_points[stod(data[1])][stod(data[2])][stod(data[3])] = stoi(data[0]);
			best_cont->put(stoi(data[0]), stod(data[1]), stod(data[2]), stod(data[3]));
		}
		myfile.close();
	}
	
	myfile.open("all_angles.txt");
	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			string delimiter = " ";

			size_t pos = 0;
			string token;
			string data[4];
			int cnt = 0;
			while ((pos = line.find(delimiter)) != string::npos) {
				token = line.substr(0, pos);
				data[cnt++] = token;
				line.erase(0, pos + delimiter.length());
			}
			data[cnt] = line;
			best_orient[data[0]] = {stod(data[1]), stod(data[2]), stod(data[3])};
		}
		myfile.close();
	}
}


void General1Class::OutputForLAMMPS(){
	
}


