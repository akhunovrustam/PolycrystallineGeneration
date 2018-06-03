#include "GeneticAlgoForOrientationsClass.hh"
#include <random>
#include <sstream>
#include <iostream>


using namespace std;
using namespace voro;

GeneticAlgoForOrientationsClass::GeneticAlgoForOrientationsClass(int pop_size)
{
	double mutation_max_applitude_e = pow(pow(half_boxside*2, 3)/particles, 1.0/3.0) / 2;
	
	population_size = pop_size;
	distribution = normal_distribution<double>(mutation_max_applitude_e / 2, mutation_max_applitude_e / 6);
}

double GeneticAlgoForOrientationsClass::original_distribution(string point)
{
	double k = 1 / pow(penalty_step, 3) / particles;

	istringstream f(point.c_str());
    string s;
	getline(f, s, '|');
	double alpha = stod(s);
	getline(f, s, '|');
	double beta = stod(s);
	getline(f, s, '|');
	double gamma = stod(s);
	return k;
}

orient_unit* GeneticAlgoForOrientationsClass::relative_euler(container* con, orient_unit *abs_euler)
{
	orient_unit *output = new orient_unit[population_size];
	map<int, point_for_crossover> id1_to_coords = ind_points_map(con);
	voronoicell_neighbor c;
	vector<int> neigh;
	double x, y, z;
	int id;

		
	c_loop_all cl(*(con));
	if(cl.start()) do if(con->compute_cell(c, cl)) {
		cl.pos(x,y,z);
		id=cl.pid();
	
		c.neighbors(neigh);
		stringstream sign;
		sign << x << "|" << y << "|" << z;
		string signstr = sign.str();
		for (vector<int>::iterator it = neigh.begin() ; it != neigh.end(); ++it)
		{
			stringstream sign_i0;
			sign_i0 << id1_to_coords[*it].x << "|" << id1_to_coords[*it].y << "|" << id1_to_coords[*it].z;
			stringstream sign_i;
			sign_i << signstr << ":" << sign_i0.str();
			string signstr10 = sign_i0.str();
			string signstr1 = sign_i.str();
			for (int i = 0; i < population_size; i++)
			{
				// cout << i << " " << signstr1 << " bla\n";
				output[i][signstr1] 
					= {	abs(abs_euler[i][signstr].alpha - abs_euler[i][signstr10].alpha), 
											abs(abs_euler[i][signstr].beta - abs_euler[i][signstr10].beta),
											abs(abs_euler[i][signstr].gamma - abs_euler[i][signstr10].gamma)	};
	
			}
		}
	} while (cl.inc());
	
	cout << "size: " << output[0].size() << "\n";
	return output;
}

double GeneticAlgoForOrientationsClass::fitness_penalty(int points_number, double (*original_distribution)(string), 
	map<string, double> current_distribution, string output)
{
	double sum = 0;
	ofstream myfile;
	if (output != "")
	{
		myfile.open (output.c_str());
		myfile << "x  y  z   error\n";
	}
	
	for (map<string, double>::iterator it=current_distribution.begin(); it!=current_distribution.end(); ++it){
		// cout << "original: " << it->first << " => " << (*original_distribution)(it->first) << "\n";
		double dif = (*original_distribution)(it->first) - it->second;
		sum += pow(dif, 2);
		
		istringstream f(it->first.c_str());
		string s;
		getline(f, s, '|');
		double alpha = stod(s);
		getline(f, s, '|');
		double beta = stod(s);
		getline(f, s, '|');
		double gamma = stod(s);
		
		if (output != "")
			myfile << alpha << "  " << beta << "  " << gamma << "  " << dif << "\n";
	}
	
	if (output != "")
		myfile.close();
	return sum / points_number;
}

double GeneticAlgoForOrientationsClass::size_penalty(orient_unit parent_rel, string output)
{
	map<string, double> current_distribution;
	
	double koef = 1 / pow(penalty_step, 3) / particles;

	for (int i = 0; i < penalty_steps; i++)
		for (int j = 0; j < penalty_steps; j++)
			for (int k = 0; k < penalty_steps; k++)
			{
				stringstream index;
				index << (i*penalty_step + penalty_step/2) << "|" << (j*penalty_step + penalty_step/2) 
					<< "|" << (k*penalty_step + penalty_step/2);
				
				string index_real = index.str();
				current_distribution[index_real] = -0.0000001;
			}
	
	// cout << "size current dist: " << current_distribution.size() << "\n";
	// cout << "check: " << current_distribution["0.075|0.675|0.225"] << "\n";
	
	for (orient_unit::iterator it = parent_rel.begin(); it != parent_rel.end(); it++){
		// cout << it->second.alpha << " " << it->second.beta << " " << it->second.gamma << "\n";
		int i = floor(it->second.alpha / (penalty_step));
		int j = floor(it->second.beta / (penalty_step));
		int k = floor(it->second.gamma / (penalty_step));
		
			// cout << "step: " << it->second.alpha << " " << penalty_step << " " << it->second.alpha / penalty_step << "\n";
			// cout << "out of range angles: " << it->second.alpha << " " << it->second.beta << " " << it->second.gamma << "\n";
			// cout << "out of range indexes: " << i << " " << j << " " << k << "\n";
			// exit(0);
		if (i == 10) i = 9;
		if (j == 10) j = 9;
		if (k == 10) k = 9;
		// if (i < 0 || i > penalty_steps || j < 0 || j > penalty_steps || k < 0 || k > penalty_steps)
		// {
			// continue;
		// }
		
		double alpha = i*penalty_step + penalty_step/2;
		double beta = j*penalty_step + penalty_step/2;
		double gamma = k*penalty_step + penalty_step/2;
		
		stringstream index;
		index << (alpha) << "|" << (beta) << "|" << (gamma);
		string index_real = index.str();
		if (current_distribution[index_real] == 0)
		{
			cout << "double check: " << current_distribution[index_real] << "\n";
			cout << "index: " << index_real << "\n";
			cout << "out of range angles: " << it->second.alpha << " " << it->second.beta << " " << it->second.gamma << "\n";
			cout << "out of range indexes: " << i << " " << j << " " << k << "\n";
			exit(0);
		}
	
		current_distribution[index_real] = current_distribution[index_real] + koef;
	}
	
	// cout << "size current dist: " << current_distribution.size() << "\n";
	// stringstream index;
	// index << (penalty_step/2) << "|" << (penalty_step/2) << "|" << (penalty_step/2);
	
	// string index_real = index.str();
	// string index_real = "0.825|0.975|1.425";
	// cout << "cube: " << index_real << " - " << current_distribution[index_real] << "\n";
	
	// int amnt = 0;
	// for (map<string, double>::iterator it = current_distribution.begin(); it != current_distribution.end(); it++)
		// if (it->second != 0)
			// amnt++;
			// cout << it->first << " " << it->second << "\n";
	// cout << "how much not null: " << amnt << "\n"; 
	
	// exit(0);
	double ret = fitness_penalty(pow(penalty_steps, 3), this->original_distribution, current_distribution, output);
	if (ret > 1000)
	{
		for (map<string, double>::iterator it = current_distribution.begin(); it != current_distribution.end(); it++)
			if (it->second != 0)
				cout << it->first << " " << it->second << "\n";
		exit(0);
	}
	// std:cout << "penalty " << ret << "\n";
	return ret;
}

double GeneticAlgoForOrientationsClass::mutate_dist(double mutation_max_applitude_e, int dist_num)
{
	switch(dist_num){
		//uniform dist
		case 0:
		return rnd()*mutation_max_applitude_e;

		//normal dist
		case 1:
		default_random_engine generator;
		normal_distribution<double> distribution(mutation_max_applitude_e / 2, mutation_max_applitude_e / 6);

		double number = distribution(generator);
		return number;
	}
}

euler_angles GeneticAlgoForOrientationsClass::mutate(euler_angles angles, double min, double max, bool adopted_shift, int dist_num)
{
	double mutation_max_applitude_e = M_PI / 2 / 100;
	
	double val = angles.alpha;
	int max_iter = 100;
	int iter = 0;
	do {
		iter++;
		if (iter > max_iter) {
			val = (max + min)/2;
			break;
		}
		if (adopted_shift)val += mutate_dist(mutation_max_applitude_e, dist_num);
		else val += rnd()*mutation_max_applitude;
	} while (val < min || val > max);
	angles.alpha = val;
	
	val = angles.beta;
	max_iter = 100;
	iter = 0;
	do {
		iter++;
		if (iter > max_iter) {
			val = (max + min)/2;
			break;
		}
		if (adopted_shift)val += mutate_dist(mutation_max_applitude_e);
		else val += rnd()*mutation_max_applitude*(rnd() > 0.5 ? -1 : 1);
	} while (val < min || val > max);
	angles.beta = val;
	
	val = angles.gamma;
	max_iter = 100;
	iter = 0;
	do {
		iter++;
		if (iter > max_iter) {
			val = (max + min)/2;
			break;
		}
		if (adopted_shift)val += mutate_dist(mutation_max_applitude_e);
		else val += rnd()*mutation_max_applitude*(rnd() > 0.5 ? -1 : 1);
	} while (val < min || val > max);
	angles.gamma = val;
	
	return angles;
}

double GeneticAlgoForOrientationsClass::reinit(double min, double max)
{
	return min+rnd()*(max-min);
}

double GeneticAlgoForOrientationsClass::reinit_angles()
{
	return rnd() * M_PI / 2;
}

void GeneticAlgoForOrientationsClass::mutation(orient_unit **con1, bool reinit_flag, double mutation_probability, int dist_num)
{
	orient_unit *con = *con1;
	for (int i = surviving_size; i < population_size; i++){
		
		orient_unit ind1 = con[i];
		for (orient_unit::iterator it = ind1.begin(); it != ind1.end(); it++)
		{
			double num = rnd();
			if (num < mutation_probability){
				if (reinit_flag == false){
					ind1[it->first] = mutate(it->second, a_min, a_max, true, dist_num);
				} else {
					ind1[it->first] = {reinit_angles(), reinit_angles(), reinit_angles()};
				}
			}
		}
	}
}

int GeneticAlgoForOrientationsClass::tournament_selection(orient_unit* parents, int selected)
{
	int index1 = rnd()*population_size;
	int index2 = rnd()*population_size;
	if (index1 == index2 && index2 == population_size - 1) index1--;
	if (index1 == index2) index1++;
	
	double penalty1 = size_penalty(parents[index1]);
	double penalty2 = size_penalty(parents[index2]);
	
	if ((penalty1 < penalty2 && selected == -1) || (selected != -1 && selected == index2)){
		return index1; 
	}
	else {
		return index2;
	}
}


point_for_crossover* GeneticAlgoForOrientationsClass::get_ind_points(container* ind1, bool to_print)
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
			
			cout << "failed ind points: " << " block:" << h << " index:" << k  << " coords:" << x1 << " " << y1 << " " << z1 << "\n";
		}
}

map<int, int> GeneticAlgoForOrientationsClass::ind_to_ind(container* ind1, container* ind2)
{
	map<int, int> id_to_id;
	for (int i = 0; i < ind1->nx*ind1->ny*ind1->nz; i++)
		for (int j = 0; j < ind1->co[i]; j++)
		{
			double *pp=ind1->p[i]+3*j;
			double x = *(pp++);
			double y = *(pp++);
			double z = *(pp++);
			
			map<double, int> id_to_dist;
			for (int h = 0; h < ind2->nx*ind2->ny*ind2->nz; h++)
				for (int k = 0; k < ind2->co[h]; k++)
				{
					double *pp1=ind2->p[h]+3*k;
					double x1 = *(pp1++);
					double y1 = *(pp1++);
					double z1 = *(pp1++);
					
					double cur_dist = sqrt(pow(x-x1, 2) + pow(y-y1, 2) + pow(z-z1, 2));
					id_to_dist[cur_dist] = ind2->id[h][k];
				}
			
			for (map<double, int>::iterator it=id_to_dist.begin(); it!=id_to_dist.end(); ++it)
			{
				bool dist_exist = false;
				for (map<int, int>::iterator it2=id_to_id.begin(); it2!=id_to_id.end(); ++it2)
					if (it2->second == it->second)
					{
						dist_exist = true;
						break;
					}
				
				if (dist_exist) continue;
				
				id_to_id[ind1->id[i][j]] = it->second;
				break;
			}
			
		}
	return id_to_id;
}

map<int, point_for_crossover> GeneticAlgoForOrientationsClass::ind_points_map(container* ind)
{	
	map<int, point_for_crossover> id2_to_coords;
	for (int i = 0; i < ind->nx*ind->ny*ind->nz; i++)
		for (int j = 0; j < ind->co[i]; j++)
		{
			double *pp=ind->p[i]+3*j;
			point_for_crossover tmp_point = {*(pp++), *(pp++), *(pp++), i, j, false};
			id2_to_coords[ind->id[i][j]] = tmp_point;
		}
	return id2_to_coords;
}

void GeneticAlgoForOrientationsClass::select_interchange_regions(container* con, map<int, point_for_crossover> *id1_to_coords, int id)
{
	int first_particles = 1 + ceil(rnd()*particles);
	int selected_particles = 0;
	
	voronoicell_neighbor c;
	vector<int> neigh;
	
	int *points;
	points = new int[first_particles];
	points[selected_particles] = id;
	(*id1_to_coords)[id].is_interchanged = true;
	
	selected_particles++;
	int checked_particles = 0;
	
	// selecting region for a crossover
	while (selected_particles < first_particles){
		if(con->compute_cell(c, (*id1_to_coords)[points[checked_particles]].block, 
			(*id1_to_coords)[points[checked_particles]].particle)) {
			
			
			checked_particles++;
			
			c.neighbors(neigh);
			
			for (vector<int>::iterator it = neigh.begin() ; it != neigh.end(); ++it)
				if (selected_particles < first_particles) 
				{
					// printf("selected_particels %i\n", selected_particles);
					points[selected_particles++] = *it;
					(*id1_to_coords)[*it].is_interchanged = true;
				}
				else break;
				
		}
		else {
			// printf("FAIL FOR individual - point_id: %i, partic_per_block: %i, block: %i, partic_index: %i, xyz: %f %f %f\n", points[checked_particles], ind1->co[(*id1_to_coords)[points[checked_particles]].block], 
				// (*id1_to_coords)[points[checked_particles]].block, (*id1_to_coords)[points[checked_particles]].particle,
				// (*id1_to_coords)[points[checked_particles]].x, (*id1_to_coords)[points[checked_particles]].y,
				// (*id1_to_coords)[points[checked_particles]].z);
			
			// get_ind_points(ind1, true);
			exit(0);
			break;
		}
	}
	delete [] points;
}

void GeneticAlgoForOrientationsClass::select_interchange_randomly(map<int, point_for_crossover> *id1_to_coords)
{
	for (map<int, point_for_crossover>::iterator it = (*id1_to_coords).begin() ; it != (*id1_to_coords).end(); ++it)
		if (rnd() < 0.5)
			(*id1_to_coords)[it->first].is_interchanged = true;
}

orient_unit* GeneticAlgoForOrientationsClass::crossover_by_mapping(container* con, orient_unit* parents, orient_unit* parents_rel,
	int crossover_points)
{
	int offspring_amount = surviving_size;
	voronoicell_neighbor c1;
	orient_unit *offspring_containers = new orient_unit[population_size];
	
	map<double, int> penalty_rating;
	if (surviving_size != 0)
		for (int i = 0; i < population_size; i++)
			penalty_rating[size_penalty(parents[i])] = i;
	
	int survived = 0;
	for (map<double, int>::iterator it=penalty_rating.begin(); it!=penalty_rating.end(); ++it)
		if (survived < surviving_size)
			offspring_containers[survived++] = parents[it->second];
	
	while (offspring_amount < population_size)
	{
		cout << "off amount " << offspring_amount << "\n";
		//tournament select of 2 individuals for a crossover
		int indx = tournament_selection(parents_rel);
		orient_unit ind1 = parents[indx];
		orient_unit ind2 = parents[tournament_selection(parents_rel, indx)];
		
		//making mapping ID => xyz_bool where bool is if point goes to ind1 or not
		map<int, point_for_crossover> id1_to_coords = ind_points_map(con);
		
		
		if (crossover_points != -1)
		{
			select_interchange_regions(con, &id1_to_coords, id1_to_coords.begin()->first);
			if (crossover_points == 2)
			{
				map<int, point_for_crossover>::iterator it=id1_to_coords.begin();
				select_interchange_regions(con, &id1_to_coords, (++it)->first);
			}
		} else {
			select_interchange_randomly(&id1_to_coords);
		}
		orient_unit offspring1;
		orient_unit offspring2;
		
		for (map<int, point_for_crossover>::iterator it=id1_to_coords.begin(); it!=id1_to_coords.end(); ++it){
			stringstream ss;
			ss << it->second.x << "|" << it->second.y << "|" << it->second.z;
			string indx = ss.str();
			if (it->second.is_interchanged)
			{
				offspring1[indx.c_str()] = ind2[indx.c_str()];
				offspring2[indx.c_str()] = ind1[indx.c_str()];
			} else {
				offspring1[indx.c_str()] = ind1[indx.c_str()];
				offspring2[indx.c_str()] = ind2[indx.c_str()];
			}
		}
		// cout << "cross6 - " << offspring_amount << "\n";
		
		offspring_containers[offspring_amount++] = offspring1;
		offspring_containers[offspring_amount++] = offspring2;
		// cout << "cross7 - " << offspring_amount << "\n";
	}
	
	// cout << "cross end - " << offspring_amount << "\n";
	
	return offspring_containers;
}

double GeneticAlgoForOrientationsClass::normal_dist(double x)
{
	double sigma = 10;
	double mju = 0;
	
	return 1/(sigma*sqrt(2*M_PI))*exp(-pow(x - mju, 2)/2/pow(sigma, 2));
}

double** GeneticAlgoForOrientationsClass::compute_cell_sizes(container** con)
{
	double x, y, z;
	int id;
	
	double **real_sizes;
	real_sizes = new double*[population_size];
	for (int i = 0; i < population_size; i++)
		real_sizes[i] = new double[particles];
	for (int i = 0; i < population_size; i++){
		voronoicell_neighbor c;
		vector<double> v;
		
		c_loop_all cl(*(con[i]));
		int loop_counter = 0;
		if(cl.start()) do if(con[i]->compute_cell(c,cl)) {
			// printf("computed individual %i\n", i);
			cl.pos(x,y,z);
			id=cl.pid();
		
			// Gather information about the computed Voronoi cell
			c.vertices(x,y,z,v);

			int planes_size = 1;
			double max = 0;
			for (vector<double>::iterator it = v.begin(); it != v.end(); ++it){
				double x = *it;
				double y = *(++it);
				double z = *(++it);
				
				for (vector<double>::iterator it2 = v.begin(); it2 != v.end(); ++it2){
					double x1 = *it2;
					double y1 = *(++it2);
					double z1 = *(++it2);
				
					double dist = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1) + (z - z1)*(z - z1));
						
					// cout << "tes vertex " << dist << " = " << x << " - " << x1 << " " << y << " - " << y1 << " " << z << " - " << z1 << "\n"; 
					// cout << dist << " " << "" << "\n"; 
					if (dist > max) max = dist;
				}
			}
			// cout << "real size before\n";
			real_sizes[i][loop_counter++] = max;
			// cout << "real size after\n";
		} while (cl.inc());
	}
	
	return real_sizes;
}

void GeneticAlgoForOrientationsClass::write_penalty_step(std::string filename, int penalties, double penalty)
{	
	ofstream myfile;
	myfile.open (filename.c_str(), ofstream::out | ofstream::app);
	if (penalties == population_size)
		myfile << "Step  Penalty\n";
	myfile << penalties << "  " << penalty << "\n";
	myfile.close();
}

void GeneticAlgoForOrientationsClass::output_data(string filename, double* size_dist)
{
	/*double avg = 0;
	map<double, double> current_distribution;
	
	for (int i = 0; i < particles; i++)
		avg += size_dist[i];
	
	avg = avg / particles;
	
	//normalization factor
	double coef = 1 / penalty_step / particles;
	for (int i = 0; i < penalty_steps; i++)
	{
		current_distribution[i*penalty_step + penalty_step/2] = 0;
		for (int j = 0; j < particles; j++){
			// cout << "size: " << size_dist[j] << " i: " << i*penalty_step << "\n";
			if ((size_dist[j] / avg) > i*penalty_step && (size_dist[j] / avg) < (i+1)*penalty_step)
				current_distribution[i*penalty_step + penalty_step/2] = current_distribution[i*penalty_step + penalty_step/2] + coef;
		}
	}
	
	ofstream myfile;
	myfile.open (filename.c_str());
	myfile << "Size(A)  Probability_real    Probability_expected\n";
	for (map<double, double>::iterator it=current_distribution.begin(); it!=current_distribution.end(); ++it)
	{
		myfile << it->first << "  " << it->second << "  " << original_distribution(it->first) << "\n";
	}
	myfile.close();*/
}