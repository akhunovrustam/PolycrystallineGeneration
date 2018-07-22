#include "GeneticAlgoForShiftingClass.hh"
#include <random>
#include <thread> 
#include <chrono>
#include <iterator>
#include "lammps/src/library.h"
#include "lammps/src/lammps.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


using namespace std;
using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds
using namespace LAMMPS_NS;


GeneticAlgoForShiftingClass::GeneticAlgoForShiftingClass(int pop_size)
{
	double mutation_max_applitude_e = LATTICE / 2;
	
	population_size = pop_size;
	distribution = normal_distribution<double>(mutation_max_applitude_e / 2, mutation_max_applitude_e / 6);
}

double GeneticAlgoForShiftingClass::original_distribution(double point)
{
	double sigma = 0.35;
	double mju = 0.0;
	double res = 1/(point*sigma*sqrt(2*M_PI))*exp(-pow(log(point) - mju, 2)/(2*pow(sigma, 2)));
	return res;
}


double GeneticAlgoForShiftingClass::EnergyCalc(map<int, atoms> atm){
	
	int cnt = 0;
	for (map<int, atoms>::iterator it = atm.begin(); it != atm.end(); it++)
		cnt += it->second.size();
	
	stringstream buffer;
	ofstream myfile;
	myfile.open("fcc_lattice.data");
	
	myfile << "Oriented crystal file" << endl;
	myfile << endl;
	myfile << cnt << " atoms" << endl;
	myfile << "0 bonds" << endl;
	myfile << "0 angles" << endl;
	myfile << "0 dihedrals" << endl;
	myfile << "0 impropers" << endl;
	myfile << endl;
	myfile << "1 atom types" << endl;
	myfile << endl;
	myfile << "0 " << half_boxside*2 << " xlo xhi" << endl;
	myfile << "0 " << half_boxside*2 << " ylo yhi" << endl;
	myfile << "0 " << half_boxside*2 << " zlo zhi" << endl;
	myfile << endl;
	myfile << "Atoms" << endl;
	myfile << endl;

	cnt = 0;
	for(map<int, atoms>::iterator it = atm.begin(); it != atm.end(); it++)
	{
		for (atoms::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
		{
			cnt++;
			myfile << cnt << " 1 " << it2->second.x << " " << it2->second.y << " " << it2->second.z << endl;
			
		}
	}

	myfile.close();
	
	char **lmparg = new char*[8];
	lmparg[0] = NULL;                 // required placeholder for program name
	lmparg[1] = (char *) "-screen";
	lmparg[2] = (char *) "none";
	lmparg[3] = (char *) "-echo";
	lmparg[4] = (char *) "none";
	
	LAMMPS *lmp;
	lammps_open_no_mpi(5, lmparg, (void **) &lmp);
	lammps_file(lmp, (char *) "lammps_script.lmp");
	double etot;
	etot = lammps_get_thermo(lmp, (char *) "etotal");
	
	delete [] lmparg;
	cout << "total eng: " << etot << endl;
	lammps_close(lmp);
	return etot;
}

double GeneticAlgoForShiftingClass::pe(string filename, map<int, atoms> atm){
	
	int cnt_atoms = 0;
	for (map<int, atoms>::iterator it = atm.begin(); it != atm.end(); it++)
		cnt_atoms += it->second.size();
	
	stringstream buffer;
	ofstream myfile;
	myfile.open(filename + "/fcc_lattice.data");
	
	myfile << "Oriented crystal file" << endl;
	myfile << endl;
	myfile << cnt_atoms << " atoms" << endl;
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

	int cnt = 0;
	int type = 1;
	for(map<int, atoms>::iterator it = atm.begin(); it != atm.end(); it++)
	{
		for (atoms::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
		{
			cnt++;
			myfile << cnt << " " << type << " " << it2->second.x << " " << it2->second.y << " " << it2->second.z << endl;
			
		}
		// type++;
	}

	myfile.close();
	
	
	char **lmparg = new char*[8];
	lmparg[0] = NULL;                 // required placeholder for program name
	lmparg[1] = (char *) "-v";
	lmparg[2] = (char *) "f";
	char str3[1024];
	sprintf(str3, "%s", filename.c_str());
	lmparg[3] = str3;
	// lmparg[3] = (char *) "aasasfsdf";
	lmparg[4] = (char *) "-e";
	lmparg[5] = (char *) "none";
	lmparg[6] = (char *) "-screen";
	lmparg[7] = (char *) "none";
	

  
	cout << "peee0" << endl;
	double bla = EnergyCalc(atm);
	cout << "peee1" << endl;
	LAMMPS *lmp;
	lammps_open_no_mpi(8, lmparg, (void **) &lmp);
	lammps_file(lmp, (char *) "per_atom.lmp");
	double pe;
	pe = lammps_get_thermo(lmp, (char *) "etotal");
	lammps_close(lmp);
	
	cout << "peee2" << endl;
	
	double radius = LATTICE*(sqrt(2) + 2)/4;
	ifstream dumpfile(filename + "/test.dump");
	double petot = 0;
	
	
	map<double, int>energyHist;
	
	if (dumpfile.is_open())
	{
		string line, original;
		int linenum = 0;
		ofstream outfile;
		outfile.open(filename + "/coords_cn.dump");
	
	
		while ( getline (dumpfile, line) )
		{
			// outfile << "asd;lfjasdlfjslj" << endl;
			// continue;
			if (linenum < 8)
				outfile << line << endl;
			else if (linenum == 8)
				outfile << line << "cn" << endl;
			else {
			
				string delimiter = " ";

				if (line == "") continue;
				size_t pos = 0;
				string token;
				string data[8];
				int cnt = 0;
				original = line;
				while ((pos = line.find(delimiter)) != string::npos) {
					token = line.substr(0, pos);
					// if (token == " " || token == "") continue;
					data[cnt++] = token;
					line.erase(0, pos + delimiter.length());
				}
				data[cnt] = line;
				
				energyHist[ceil(stod(data[5])/0.1)*0.1] = energyHist[ceil(stod(data[5])/0.1)*0.1] + 1;
				petot += stod(data[5]);
				double x = stod(data[2]);
				double y = stod(data[3]);
				double z = stod(data[4]);
				double x_mn = x - radius;
				double x_mx = x + radius;
				double y_mn = y - radius;
				double y_mx = y + radius;
				double z_mn = z - radius;
				double z_mx = z + radius;
				
				int cn = 0;
				for(map<int, atoms>::iterator it = atm.begin(); it != atm.end(); it++)
				{
					for (atoms::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
					{
						if (it2->second.x < x_mx && it2->second.x > x_mn && 
							it2->second.y < y_mx && it2->second.y > y_mn && 
							it2->second.z < z_mx && it2->second.z > z_mn)
							cn++;
					}
				}
				outfile << original << cn << endl;
			}
			linenum++;
		}
		outfile.close();
		dumpfile.close();
	}
	
	cout << "peee3" << endl;
	cout << "new dump" << endl;
	
	ofstream outfile;
	outfile.open(filename + "/energy_hist.txt");
	outfile << "Energy Atoms" << endl;
	for (map<double, int>::iterator it = energyHist.begin(); it != energyHist.end(); it++)
		outfile << it->first << " " << ((double)it->second/cnt_atoms) << endl;
	outfile.close();
	
	cout << "per atom energy: " << pe << " " << petot << endl;
	
	// exit(0);
	
	cout << "end of script" << endl;
	// exit(0);
}

double GeneticAlgoForShiftingClass::size_penalty(map<int, atoms> atm)
{
	return EnergyCalc(atm);
}

double GeneticAlgoForShiftingClass::mutate_dist(double mutation_max_applitude_e, int dist_num)
{
	switch(dist_num){
		//uniform dist
		case 0:
		return rnd()*mutation_max_applitude_e*(rnd() > 0.5 ? -1 : 1);

		//normal dist
		case 1:
		
		double number = distribution(generator);
		// cout << "normal number " << number << "\n";
		
		return number;
	}
}

double GeneticAlgoForShiftingClass::reinit(double min, double max)
{
	return min+rnd()*(max-min);
}

double GeneticAlgoForShiftingClass::writeAtoms(map<int, atoms> atm, string filename){
	
	int cnt = 0;
	for (map<int, atoms>::iterator it = atm.begin(); it != atm.end(); it++)
		cnt += it->second.size();
	
	stringstream buffer;
	ofstream myfile;
	myfile.open(filename);
	
	myfile << "ITEM: TIMESTEP" << endl;
	myfile << "0" << endl;
	myfile << "ITEM: NUMBER OF ATOMS" << endl;
	myfile << cnt << endl;
	myfile << "ITEM: BOX BOUNDS pp pp pp" << endl;
	myfile << "0.0000000000000000e+00 " << setprecision(16) << scientific << half_boxside << endl;
	myfile << "0.0000000000000000e+00 " << setprecision(16) << scientific << half_boxside << endl;
	myfile << "0.0000000000000000e+00 " << setprecision(16) << scientific << half_boxside << endl;
	myfile << "ITEM: ATOMS id type x y z" << fixed << endl;
	cnt = 0;
	for (map<int, atoms>::iterator it = atm.begin(); it != atm.end(); it++)
		for (atoms::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
		{
			cnt++;
			myfile << cnt << " 1 " << it2->second.x << " " << it2->second.y << " " << it2->second.z << endl;
		}
	// myfile << endl;
	myfile.close();
}

double GeneticAlgoForShiftingClass::ShiftDelta(double val){
	return val + rnd()*LATTICE/100*(rnd() < 0.5 ? -1 : 1);
}

void GeneticAlgoForShiftingClass::mutation(map<int, atoms>** con1, config cfg)
{
	map<int, atoms>* con = *con1;
	for (int i = surviving_size; i < population_size; i++){
					
		for (int j = surviving_size; j < particles; j++){
			double rnd_num = rnd();
			if (rnd_num < cfg.mp){
				double dx = ShiftDelta(0);
				double dy = ShiftDelta(0);
				double dz = ShiftDelta(0);
				atoms tmp;
				for (atoms::iterator it = con[i][j].begin(); it != con[i][j].end(); it++)
				{
					atom_coords ac = {it->second.x + dx, it->second.y + dy, it->second.z + dz};
					tmp[it->first] = ac;
				}
				con[i][j] = tmp;
			}
		}
	}
}

int GeneticAlgoForShiftingClass::tournament_selection(double* penalty, int selected)
{
	int index1 = floor(rnd()*population_size);
	int index2 = floor(rnd()*population_size);
	if (index1 == index2 && index2 == population_size - 1) index1--;
	if (index1 == index2) index1++;
	
	double penalty1 = penalty[index1];
	double penalty2 = penalty[index2];
	
	if ((penalty1 < penalty2 && selected == -1) || (selected != -1 && selected == index2)){
		return index1; 
	}
	else {
		return index2;
	}
}


point_for_crossover* GeneticAlgoForShiftingClass::get_ind_points(container* ind1, bool to_print)
{
	point_for_crossover points[particles];
	int particles_count = 0;
	
	for (int h = 0; h < ind1->nx*ind1->ny*ind1->nz; h++)
		for (int k = 0; k < ind1->co[h]; k++)
		{
			double *pp1=ind1->p[h]+3*k;
			double x1 = *(pp1++);
			double y1 = *(pp1++);
			double z1 = *(pp1++);
			point_for_crossover pnt = {x1, y1, z1, h, k, false};
			points[particles_count++] = pnt;
			
			std::cout << "failed ind points: " << " block:" << h << " index:" << k  << " coords:" << x1 << " " << y1 << " " << z1 << "\n";
		}
}


std::map<int, point_for_crossover> GeneticAlgoForShiftingClass::ind_points_map(sorted_points sortp)
{	
	std::map<int, point_for_crossover> id2_to_coords;
	for (sorted_points::iterator it = sortp.begin(); it != sortp.end(); it++)
		for (map<double, map<double, int>>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
			for (map<double, int>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++)
			{
				point_for_crossover tmp_point = {it->first, it2->first, it3->first, 0, 0, false};
				id2_to_coords[it3->second] = tmp_point;
			}

	if (id2_to_coords.size() == 0)
	{
		cout << "strange: " << sortp.size() << " " << sortp.begin()->second.size() << " " << sortp.begin()->second.begin()->second.size() << endl;
	}
	return id2_to_coords;
}

void GeneticAlgoForShiftingClass::select_interchange_regions(map<int, atoms> ind1, std::map<int, point_for_crossover> *id1_to_coords, 
	int id, int from, int to, neighbors neg)
{
	int first_particles = from + ceil(rnd()*(to - from));
	int selected_particles = 0;
	
	vector<int> neigh;
	
	int *points;
	points = new int[first_particles];
	points[selected_particles] = id;
	
	(*id1_to_coords)[id].is_interchanged = true;
	
	selected_particles++;
	int checked_particles = 0;
	
	// selecting region for a crossover
	while (selected_particles < first_particles){
		neigh = neg[points[checked_particles++]];
		for (std::vector<int>::iterator it = neigh.begin() ; it != neigh.end(); it++)
			if (selected_particles < first_particles && *it > -1) 
			{
				// printf("selected_particels %i\n", selected_particles);
				points[selected_particles++] = *it;
				(*id1_to_coords)[*it].is_interchanged = true;
			}
			else break;
	}
	
	delete [] points;
}

void GeneticAlgoForShiftingClass::select_interchange_randomly(std::map<int, point_for_crossover> *id1_to_coords)
{
	for (std::map<int, point_for_crossover>::iterator it = (*id1_to_coords).begin() ; it != (*id1_to_coords).end(); ++it)
		if (rnd() < 0.5)
			(*id1_to_coords)[it->first].is_interchanged = true;
}

map<int, atoms>* GeneticAlgoForShiftingClass::crossover_by_mapping(double* penalty, map<int, atoms>* atms, config cfg, 
	neighbors neg, sorted_points best_cont_points)
{
	
	map<int, atoms>* offspring_containers = new map<int, atoms>[population_size];
	
	int offspring_amount = 0;
	while (offspring_amount < population_size * cfg.ps_fac)
	{
		
		//tournament select of 2 individuals for a crossover
		int index1 = tournament_selection(penalty);
		
		map<int, atoms> ind1 = atms[index1];
		int index2 = tournament_selection(penalty, index1);
		map<int, atoms> ind2 = atms[index2];
		
		//making mapping ID => xyz_bool where bool is if point goes to ind1 or not
		std::map<int, point_for_crossover> id1_to_coords = ind_points_map(best_cont_points);
		
		// cout << "pre regions select" << endl;
		if (cfg.co != -1)
		{
			auto item = id1_to_coords.begin();
			advance( item, rnd() * id1_to_coords.size() );
			
			// cout << "first select" << endl;
			select_interchange_regions(ind1, &id1_to_coords, item->first, cfg.from, cfg.to, neg);
			if (cfg.co == 2)
			{
				auto item = id1_to_coords.begin();
				advance( item, rnd() * id1_to_coords.size() );
			
				// cout << "second select" << endl;
				select_interchange_regions(ind1, &id1_to_coords, item->first, cfg.from, cfg.to, neg);
			}
		} else {
			select_interchange_randomly(&id1_to_coords);
		}
		
		// cout << "post regions select" << endl;
		map<int, atoms> offspring1;
		map<int, atoms> offspring2;
		
		for (std::map<int, point_for_crossover>::iterator it=id1_to_coords.begin(); it!=id1_to_coords.end(); ++it){
			if (it->second.is_interchanged)
			{
				offspring1[it->first] = ind2[it->first];
				offspring2[it->first] = ind1[it->first];
			}
			else {
				offspring1[it->first] = ind1[it->first];
				offspring2[it->first] = ind2[it->first];
			}
		}
		
		offspring_containers[offspring_amount++] = offspring1;
		offspring_containers[offspring_amount++] = offspring2;
		
	}
	
	// std::cout << "cross end - " << offspring_amount << "\n";
	
	return offspring_containers;
}

double GeneticAlgoForShiftingClass::normal_dist(double x)
{
	double sigma = 10;
	double mju = 0;
	
	return 1/(sigma*sqrt(2*M_PI))*exp(-pow(x - mju, 2)/2/pow(sigma, 2));
}

/*
void GeneticAlgoForShiftingClass::output_data(std::string filename, double* size_dist)
{
	double avg = 0;
	std::map<double, double> current_distribution;
	
	for (int i = 0; i < particles; i++)
		avg += size_dist[i];
	
	avg = avg / particles;
	
	
	//normalization factor
	double coef = 1 / penalty_step / particles;
	for (int i = 0; i < penalty_steps; i++)
	{
		current_distribution[i*penalty_step + penalty_step/2] = 0;
		for (int j = 0; j < particles; j++){
			// std::cout << "size: " << size_dist[j] << " i: " << i*penalty_step << "\n";
			if ((size_dist[j] / avg) > i*penalty_step && (size_dist[j] / avg) < (i+1)*penalty_step)
				current_distribution[i*penalty_step + penalty_step/2] = current_distribution[i*penalty_step + penalty_step/2] + coef;
		}
	}
	
	ofstream myfile;
	myfile.open (filename.c_str());
	myfile << "Size(A)  Probability_real    Probability_expected\n";
	double sum = 0;
	int cnt = 0;
	for (std::map<double, double>::iterator it=current_distribution.begin(); it!=current_distribution.end(); ++it)
	{
		cnt++;
		double orig = original_distribution(it->first);
		double diff = abs(orig - it->second);
		sum += pow(diff, 2);
		myfile << it->first << "  " << it->second << "  " << orig << " " << diff << " " << endl;
	}
	myfile.close();
}*/

void GeneticAlgoForShiftingClass::write_penalty_step(std::string filename, int penalties, double penalty)
{	
	ofstream myfile;
	myfile.open (filename.c_str(), ofstream::out | ofstream::app);
	if (penalties == population_size)
		myfile << "Step  Penalty\n";
	myfile << penalties << "  " << penalty << "\n";
	myfile.close();
}