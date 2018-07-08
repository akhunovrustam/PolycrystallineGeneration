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

double GeneticAlgoForOrientationsClass::original_distribution(double point)
{
	double orig;
	if (point < M_PI/4)
		orig = 2.0/15.0*(1-cos(point));
	else orig = 2.0/15.0*(3*(sqrt(2)-1)*sin(point) - 2*(1-cos(point)));
	return orig;
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
				double f1 = abs_euler[i][signstr].alpha;
				double F = abs_euler[i][signstr].beta;
				double f2 = abs_euler[i][signstr].gamma;
				double ff1 = abs_euler[i][signstr10].alpha;
				double FF = abs_euler[i][signstr10].beta;
				double ff2 = abs_euler[i][signstr10].gamma;
				
				double m1[3][3];
				double m2[3][3];
				
				m1[0][0] = cos(f1)*cos(f2) - sin(f1)*sin(f2)*cos(F);
				m1[0][1] = sin(f1)*cos(f2) + cos(f1)*sin(f2)*cos(F);
				m1[0][2] = sin(f2)*sin(F);
				m1[1][0] = -cos(f1)*sin(f2) - sin(f1)*cos(f2)*cos(F);
				m1[1][1] = -sin(f1)*sin(f2) + cos(f1)*cos(f2)*cos(F);
				m1[1][2] = cos(f2)*sin(F);
				m1[2][0] = sin(f1)*sin(F);
				m1[2][1] = -cos(f1)*sin(F);
				m1[2][2] = cos(F);
				
				m2[0][0] = cos(ff1)*cos(ff2) - sin(ff1)*sin(ff2)*cos(FF);
				m2[0][1] = sin(ff1)*cos(ff2) + cos(ff1)*sin(ff2)*cos(FF);
				m2[0][2] = sin(ff2)*sin(FF);
				m2[1][0] = -cos(ff1)*sin(ff2) - sin(ff1)*cos(ff2)*cos(FF);
				m2[1][1] = -sin(ff1)*sin(ff2) + cos(ff1)*cos(ff2)*cos(FF);
				m2[1][2] = cos(ff2)*sin(FF);
				m2[2][0] = sin(ff1)*sin(FF);
				m2[2][1] = -cos(ff1)*sin(FF);
				m2[2][2] = cos(FF);
				
				double det = m2[0][0]*m2[1][1]*m2[2][2] + m2[0][1]*m2[1][2]*m2[2][0]
					+ m2[1][0]*m2[2][1]*m2[0][2] - m2[2][0]*m2[1][1]*m2[0][2] 
					- m2[0][0]*m2[2][1]*m2[1][2] - m2[1][0]*m2[0][1]*m2[2][2];
					
				double inv[3][3];
				inv[0][0] = (m2[1][1]*m2[2][2] - m2[1][2]*m2[2][1])/det;
				inv[1][0] = -(m2[1][0]*m2[2][2] - m2[2][0]*m2[1][2])/det;
				inv[2][0] = (m2[1][0]*m2[2][1] - m2[2][0]*m2[1][1])/det;
				inv[0][1] = -(m2[0][1]*m2[2][2] - m2[2][1]*m2[0][2])/det;
				inv[1][1] = (m2[0][0]*m2[2][2] - m2[2][0]*m2[0][2])/det;
				inv[2][1] = -(m2[0][0]*m2[2][1] - m2[2][0]*m2[0][1])/det;
				inv[0][2] = (m2[0][1]*m2[1][2] - m2[1][1]*m2[0][2])/det;
				inv[1][2] = -(m2[0][0]*m2[1][2] - m2[1][0]*m2[0][2])/det;
				inv[2][2] = (m2[0][0]*m2[1][1] - m2[1][0]*m2[0][1])/det;
				// cout << i << " " << signstr1 << " bla\n";

				double minangle = 100;
				for (int j = 0; j < 24; j++)
				{
					double tt[3][3];
					tt[0][0] = m1[0][0]*rot[j][0][0] + m1[0][1]*rot[j][1][0] + m1[0][2]*rot[j][2][0];
					tt[0][1] = m1[0][0]*rot[j][0][1] + m1[0][1]*rot[j][1][1] + m1[0][2]*rot[j][2][1];
					tt[0][2] = m1[0][0]*rot[j][0][2] + m1[0][1]*rot[j][1][2] + m1[0][2]*rot[j][2][2];
					tt[1][0] = m1[1][0]*rot[j][0][0] + m1[1][1]*rot[j][1][0] + m1[1][2]*rot[j][2][0];
					tt[1][1] = m1[1][0]*rot[j][0][1] + m1[1][1]*rot[j][1][1] + m1[1][2]*rot[j][2][1];
					tt[1][2] = m1[1][0]*rot[j][0][2] + m1[1][1]*rot[j][1][2] + m1[1][2]*rot[j][2][2];
					tt[2][0] = m1[2][0]*rot[j][0][0] + m1[2][1]*rot[j][1][0] + m1[2][2]*rot[j][2][0];
					tt[2][1] = m1[2][0]*rot[j][0][1] + m1[2][1]*rot[j][1][1] + m1[2][2]*rot[j][2][1];
					tt[2][2] = m1[2][0]*rot[j][0][2] + m1[2][1]*rot[j][1][2] + m1[2][2]*rot[j][2][2];
					
					double g11 = tt[0][0]*inv[0][0] + tt[0][1]*inv[1][0] + tt[0][2]*inv[2][0];
					double g22 = tt[1][0]*inv[0][1] + tt[1][1]*inv[1][1] + tt[1][2]*inv[2][1];
					double g33 = tt[2][0]*inv[0][2] + tt[2][1]*inv[1][2] + tt[2][2]*inv[2][2];
					
					double trace = (g11 + g22 + g33 - 1)/2;
					double tetta = acos(trace);
					
					// if (isnan(tetta))
					// {
						// cout << tetta << " " << trace << endl;
						// cout << g11 << " " << g22 << " " << g33 << endl;
						// exit(0);
					// }
					if (minangle > tetta) minangle = tetta;
				}
				output[i][signstr1] 
					= {	minangle, 0, 0 };
	
			}
		}
	} while (cl.inc());
	
	cout << "size: " << output[0].size() << "\n";
	return output;
}

double GeneticAlgoForOrientationsClass::fitness_penalty(int points_number, double (*original_distribution)(double), 
	map<double, double> current_distribution, string output)
{
	double sum = 0;
	for (map<double, double>::iterator it = current_distribution.begin(); it != current_distribution.end(); it++)
	{
		double orig;
		double par = it->first/180*M_PI;
		orig = original_distribution(par);
		sum += pow(it->second - orig, 2);
		// cout << it->first << " " << it->second << " " << orig << "\n";
	}
	
	return sum / points_number;
}

double GeneticAlgoForOrientationsClass::size_penalty(orient_unit parent_rel, string output)
{
	map<double, double> current_distribution;
	
	double koef = 1.0 / parent_rel.size() / penalty_step_orient * 10000;
	
	for (orient_unit::iterator it = parent_rel.begin(); it != parent_rel.end(); it++){
		int i = floor((it->second.alpha * 180 / M_PI) / (penalty_step_orient));
		
		double alpha = i*penalty_step_orient + penalty_step_orient/2;
		
		current_distribution[alpha] = current_distribution[alpha] + koef;
	}
	
	double ret = fitness_penalty(pow(penalty_steps_orient, 3), this->original_distribution, current_distribution, output);
	
	return ret;
}

double GeneticAlgoForOrientationsClass::mutate_dist(double mutation_max_applitude_e, int dist_num)
{
	switch(dist_num){
		//uniform dist
		case 0:
		return rnd()*mutation_max_applitude_e*(rnd() > 0.5 ? -1 : 1);

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
	double mutation_max_applitude_e = M_PI / 2 / 10;
	
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
		
		// orient_unit ind1 = con[i];
		
		for (orient_unit::iterator it = con[i].begin(); it != con[i].end(); it++)
		{
			double num = rnd();
			if (num < mutation_probability){
				if (reinit_flag == false){
					euler_angles tmp = mutate(it->second, a_min, a_max, true, dist_num);
					// cout << "mutated: " << ind1[it->first].alpha << " " << tmp.alpha << endl;
					con[i][it->first] = tmp;
				} else {
					con[i][it->first] = {reinit_angles(), reinit_angles(), reinit_angles()};
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
		// cout << "off amount " << offspring_amount << "\n";
		//tournament select of 2 individuals for a crossover
		int indx = tournament_selection(parents_rel);
		orient_unit ind1 = parents[indx];
		orient_unit ind2 = parents[tournament_selection(parents_rel, indx)];
		
		//making mapping ID => xyz_bool where bool is if point goes to ind1 or not
		map<int, point_for_crossover> id1_to_coords = ind_points_map(con);
		
		
		if (crossover_points != -1)
		{
			auto item = id1_to_coords.begin();
			advance( item, rnd() * id1_to_coords.size() );
			
			select_interchange_regions(con, &id1_to_coords, item->first);
			if (crossover_points == 2)
			{
				auto item = id1_to_coords.begin();
				advance( item, rnd() * id1_to_coords.size() );
			
				select_interchange_regions(con, &id1_to_coords, item->first);
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

void GeneticAlgoForOrientationsClass::output_data(string filename, orient_unit parent_rel)
{
	map<double, double> current_distribution;
	
	double koef = 1.0 / parent_rel.size() / penalty_step_orient * 10000;
	
	for (orient_unit::iterator it = parent_rel.begin(); it != parent_rel.end(); it++){
		int i = floor((it->second.alpha * 180 / M_PI) / (penalty_step_orient));
		
		double alpha = i*penalty_step_orient + penalty_step_orient/2;
		
		current_distribution[alpha] = current_distribution[alpha] + koef;
	}
	
	
	ofstream myfile;
	myfile.open (filename.c_str());
	myfile << "Angle(A)  Probability_real    Probability_expected\n";
	for (map<double, double>::iterator it=current_distribution.begin(); it!=current_distribution.end(); ++it)
	{
		myfile << it->first << "  " << it->second << "  " << original_distribution(it->first/180*M_PI) << "\n";
	}
	myfile.close();
}