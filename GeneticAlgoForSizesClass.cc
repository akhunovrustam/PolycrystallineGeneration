#include "GeneticAlgoForSizesClass.hh"
#include <random>
#include <thread> 
#include <chrono>

using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds


GeneticAlgoForSizesClass::GeneticAlgoForSizesClass(int pop_size)
{
	double mutation_max_applitude_e = pow(pow(half_boxside*2, 3)/particles, 1.0/3.0) / 2;
	
	population_size = pop_size;
	distribution = normal_distribution<double>(mutation_max_applitude_e / 2, mutation_max_applitude_e / 6);
}

double GeneticAlgoForSizesClass::original_distribution(double point)
{
	double sigma = 0.35;
	double mju = 0.0;
	double res = 1/(point*sigma*sqrt(2*M_PI))*exp(-pow(log(point) - mju, 2)/(2*pow(sigma, 2)));
	return res;
}

double GeneticAlgoForSizesClass::fitness_penalty(int points_number, double (*original_distribution)(double), 
	std::map<double, double> current_distribution)
{
	double sum = 0;
	for (std::map<double, double>::iterator it=current_distribution.begin(); it!=current_distribution.end(); ++it){
		// std::cout << "original: " << it->first << " => " << (*original_distribution)(it->first) << "\n";
		sum += pow((*original_distribution)(it->first) - it->second, 2);
	}
	
	return sum / points_number;
}

double GeneticAlgoForSizesClass::size_penalty(double* size_dist, int iter, int offspring_amount)
{
	std::map<double, double> current_distribution;
	double max = 0, min = 100;
	
	double avg = 0;
	
	for (int i = 0; i < particles; i++)
		avg += size_dist[i];
	
	avg = avg / particles;
	
	double k = 1 / penalty_step / particles;

	for (int i = 0; i < penalty_steps; i++)
	{
		current_distribution[i*penalty_step + penalty_step/2] = 0;
		for (int j = 0; j < particles; j++){
			// std::cout << "size: " << size_dist[j] << " i: " << i*penalty_step << "\n";
			if (size_dist[j] > i*penalty_step && size_dist[j] < (i+1)*penalty_step)
				current_distribution[i*penalty_step + penalty_step/2] = current_distribution[i*penalty_step + penalty_step/2] + k;
		}
	}
	
	double ret = fitness_penalty(penalty_steps, this->original_distribution, current_distribution);
	current_distribution.clear();
	// std:cout << "penalty " << ret << "\n";
	return ret;
}

double GeneticAlgoForSizesClass::mutate_dist(double mutation_max_applitude_e, int dist_num)
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

double GeneticAlgoForSizesClass::mutate(double val, double min, double max, bool adopted_shift, int dist_num)
{
	double mutation_max_applitude_e = pow(pow(half_boxside*2, 3)/particles, 1.0/3.0) / 2;
	
	int max_iter = 100;
	int iter = 0;
	do {
		iter++;
		if (iter > max_iter) {
			val = (max + min)/2;
			break;
		}
		if (adopted_shift)val += mutate_dist(mutation_max_applitude_e, dist_num);
		else val += rnd()*mutation_max_applitude*(rnd() > 0.5 ? -1 : 1);
	} while (val < min || val > max);

	return val;
}

double GeneticAlgoForSizesClass::reinit(double min, double max)
{
	return min+rnd()*(max-min);
}

void GeneticAlgoForSizesClass::mutation(container ***con1, bool reinit_flag, double mutation_probability, int dist_num)
{
	container **con = *con1;
	for (int i = surviving_size; i < population_size; i++){
		double new_points[particles][3];
		bool point_changed = false;
		
		container* ind1 = con[i];
		int new_points_cntr = 0;
		for (int k = 0; k < ind1->nx*ind1->ny*ind1->nz; k++)
		{
			// printf("points per box - %i - %i\n", i, ind1->co[i]);
			for (int j = 0; j < ind1->co[k]; j++)
			{
				double *pp=ind1->p[k]+3*j;
				double x = *(pp++);
				double y = *(pp++);
				double z = *(pp++);
				
				int rnd_index1 = ceil(rnd()*ceil(sqrt(1/mutation_probability)));
				int rnd_index2 = ceil(rnd()*ceil(sqrt(1/mutation_probability)));
				
				// std::cout << "rnd mutat " << (rnd_index1 == rnd_index2) << "\n";
				if (rnd_index1 == rnd_index2){
					point_changed = true;
					if (reinit_flag == false){
						new_points[new_points_cntr][0] = mutate(x, x_min, x_max, true, dist_num);
						new_points[new_points_cntr][1] = mutate(y, y_min, y_max, true, dist_num);
						new_points[new_points_cntr][2] = mutate(z, z_min, z_max, true, dist_num);
					} else {
						new_points[new_points_cntr][0] = reinit(x_min, x_max);
						new_points[new_points_cntr][1] = reinit(y_min, y_max);
						new_points[new_points_cntr][2] = reinit(z_min, z_max);
					}
				}
				else {
					new_points[new_points_cntr][0] = x;
					new_points[new_points_cntr][1] = y;
					new_points[new_points_cntr][2] = z;
				}
				new_points_cntr++;
				
			}
		}
		
		if (point_changed){
			container *con_subst;
			con_subst = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
			
			for(int j = 0; j < particles; j++) {
				con_subst->put(j,new_points[j][0],new_points[j][1],new_points[j][2]);
			}
			delete con[i];
			con[i] = con_subst;
		}
	}
}

int GeneticAlgoForSizesClass::tournament_selection(double** size_dist, int iter, int offspring_amount, int selected)
{
	int index1 = rnd()*population_size;
	int index2 = rnd()*population_size;
	if (index1 == index2 && index2 == population_size - 1) index1--;
	if (index1 == index2) index1++;
	
	double penalty1 = size_penalty(size_dist[index1], iter, offspring_amount);
	double penalty2 = size_penalty(size_dist[index2], iter, offspring_amount);
	
	if ((penalty1 < penalty2 && selected == -1) || (selected != -1 && selected == index2)){
		return index1; 
	}
	else {
		return index2;
	}
}


point_for_crossover* GeneticAlgoForSizesClass::get_ind_points(container* ind1, bool to_print)
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

std::map<int, int> GeneticAlgoForSizesClass::ind_to_ind(container* ind1, container* ind2)
{
	std::map<int, int> id_to_id;
	for (int i = 0; i < ind1->nx*ind1->ny*ind1->nz; i++)
		for (int j = 0; j < ind1->co[i]; j++)
		{
			double *pp=ind1->p[i]+3*j;
			double x = *(pp++);
			double y = *(pp++);
			double z = *(pp++);
			
			std::map<double, int> id_to_dist;
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
			
			for (std::map<double, int>::iterator it=id_to_dist.begin(); it!=id_to_dist.end(); ++it)
			{
				bool dist_exist = false;
				for (std::map<int, int>::iterator it2=id_to_id.begin(); it2!=id_to_id.end(); ++it2)
					if (it2->second == it->second)
					{
						dist_exist = true;
						break;
					}
				
				if (dist_exist) continue;
				
				id_to_id[ind1->id[i][j]] = it->second;
				break;
			}
			id_to_dist.clear();
		}
	return id_to_id;
}

std::map<int, point_for_crossover> GeneticAlgoForSizesClass::ind_points_map(container* ind)
{	
	std::map<int, point_for_crossover> id2_to_coords;
	for (int i = 0; i < ind->nx*ind->ny*ind->nz; i++)
		for (int j = 0; j < ind->co[i]; j++)
		{
			double *pp=ind->p[i]+3*j;
			point_for_crossover tmp_point = {*(pp++), *(pp++), *(pp++), i, j, false};
			id2_to_coords[ind->id[i][j]] = tmp_point;
		}
	return id2_to_coords;
}

void GeneticAlgoForSizesClass::select_interchange_regions(container* ind1, std::map<int, point_for_crossover> *id1_to_coords, int id)
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
		if(ind1->compute_cell(c, (*id1_to_coords)[points[checked_particles]].block, 
			(*id1_to_coords)[points[checked_particles]].particle)) {
			
			
			checked_particles++;
			
			c.neighbors(neigh);
			
			for (std::vector<int>::iterator it = neigh.begin() ; it != neigh.end(); ++it)
				if (selected_particles < first_particles) 
				{
					// printf("selected_particels %i\n", selected_particles);
					points[selected_particles++] = *it;
					(*id1_to_coords)[*it].is_interchanged = true;
				}
				else break;
				
		}
		else {
			printf("FAIL FOR individual - point_id: %i, partic_per_block: %i, block: %i, partic_index: %i, xyz: %f %f %f\n", points[checked_particles], ind1->co[(*id1_to_coords)[points[checked_particles]].block], 
				(*id1_to_coords)[points[checked_particles]].block, (*id1_to_coords)[points[checked_particles]].particle,
				(*id1_to_coords)[points[checked_particles]].x, (*id1_to_coords)[points[checked_particles]].y,
				(*id1_to_coords)[points[checked_particles]].z);
			
			get_ind_points(ind1, true);
			exit(0);
			break;
		}
	}
	
	delete [] points;
}

void GeneticAlgoForSizesClass::select_interchange_randomly(std::map<int, point_for_crossover> *id1_to_coords)
{
	for (std::map<int, point_for_crossover>::iterator it = (*id1_to_coords).begin() ; it != (*id1_to_coords).end(); ++it)
		if (rnd() < 0.5)
			(*id1_to_coords)[it->first].is_interchanged = true;
}

container** GeneticAlgoForSizesClass::crossover_by_mapping(container** con, double** size_dist, int iter, int crossover_points)
{
	
	int offspring_amount = surviving_size;
	voronoicell_neighbor c1;
	container** offspring_containers = new container*[population_size + 1];
	
	std::map<double, int> penalty_rating;
	for (int i = 0; i < population_size; i++)
		penalty_rating[size_penalty(size_dist[i], 0, 0)] = i;
	
	
	int survived = 0;
	for (std::map<double, int>::iterator it=penalty_rating.begin(); it!=penalty_rating.end(); ++it)
		if (survived < surviving_size)
			offspring_containers[survived++] = con[it->second];
	
	penalty_rating.clear();
	
	while (offspring_amount < population_size)
	{
		
		//tournament select of 2 individuals for a crossover
		int index1 = tournament_selection(size_dist, iter, offspring_amount);
		
		container* ind1 = con[index1];
		int index2 = tournament_selection(size_dist, iter, offspring_amount, index1);
		container* ind2 = con[index2];
		
		//making mapping ID => xyz_bool where bool is if point goes to ind1 or not
		std::map<int, point_for_crossover> id1_to_coords = ind_points_map(ind1);
		std::map<int, point_for_crossover> id2_to_coords = ind_points_map(ind2);
		
		//making mapping from particels of ind1 to ind2
		std::map<int, int> id_to_id = ind_to_ind(ind1, ind2);
		
		if (crossover_points != -1)
		{
			select_interchange_regions(ind1, &id1_to_coords, id1_to_coords.begin()->first);
			if (crossover_points == 2)
			{
				std::map<int, point_for_crossover>::iterator it=id1_to_coords.begin();
				select_interchange_regions(ind1, &id1_to_coords, (++it)->first);
			}
		} else {
			select_interchange_randomly(&id1_to_coords);
		}
		container *offspring1 = NULL;
		container *offspring2 = NULL;
		
		offspring1 = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
		// std::cout << "cross5.1" << "\n";
		offspring2 = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
		// std::cout << "cross5.2" << "\n";
		
		std::map<double, std::map<double, double> > points_map1;
		std::map<double, std::map<double, double> > points_map2;
		// std::cout << "before copy of points\n";
		for (std::map<int, point_for_crossover>::iterator it=id1_to_coords.begin(); it!=id1_to_coords.end(); ++it){
				
			//here might be problems because of =
			double x = it->second.x, y = it->second.y, z = it->second.z; 
			double x2 = id2_to_coords[id_to_id[it->first]].x, y2 = id2_to_coords[id_to_id[it->first]].y, z2 = id2_to_coords[id_to_id[it->first]].z; 
			if (it->second.is_interchanged)
			{
				// std::cout << "new point1-1 " << x << " " << y << " " << z << " points map: " << points_map1[x][y] << "\n";
				// std::cout << "new point1-2 " << x2 << " " << y2 << " " << z2 << "\n";
				if (points_map1[x][y] != z)
				{
					points_map1[x][y] = z;
					offspring1->put(it->first, x, y, z);
				} else {
					double shifted_x = mutate(x, x_min, x_max);
					points_map1[shifted_x][y] = z;
					offspring1->put(it->first + 10000, shifted_x, y, z);
				}
				
				if (points_map2[x2][y2] != z2)
				{
					points_map2[x2][y2] = z2;
					offspring2->put(it->first, x2, y2, z2);
				} else {
					double shifted_x = mutate(x2, x_min, x_max);
					points_map2[shifted_x][y2] = z2;
					offspring2->put(it->first, shifted_x, y2, z2);
				}
				
			}
			else {
				// std::cout << "new point2-2 " << x << " " << y << " " << z << "\n";
				// std::cout << "new point2-1 " << x2 << " " << y2 << " " << z2 << " points map: " << points_map1[x2][y2] << "\n";
				if (points_map1[x2][y2] != z2)
				{
					// std::cout << "NON mutated x " << x2 << "\n";
					points_map1[x2][y2] = z2;
					offspring1->put(it->first, x2, y2, z2);
				} else {
					double shifted_x = mutate(x2, x_min, x_max);
					// std::cout << "mutated x " << shifted_x << "\n";
					points_map1[shifted_x][y2] = z2;
					offspring1->put(it->first + 10000, shifted_x, y2, z2);
				}
				
				if (points_map2[x][y] != z)
				{
					points_map2[x][y] = z;
					offspring2->put(it->first, x, y, z);
				} else {
					double shifted_x = mutate(x, x_min, x_max);
					points_map2[shifted_x][y] = z;
					offspring2->put(it->first, shifted_x, y, z);
				}
				// offspring2->put(it->first, x, y, z);
				// offspring1->put(it->first, x2, y2, z2);
			}
		}
		
		// std::cout << "cross6 - " << offspring_amount << "\n";
		
		offspring_containers[offspring_amount++] = offspring1;
		offspring_containers[offspring_amount++] = offspring2;
		// std::cout << "cross7 - " << offspring_amount << "\n";
		points_map1.clear();
		points_map2.clear();
		id1_to_coords.clear();
		id2_to_coords.clear();
		id_to_id.clear();
	}
	
	// std::cout << "cross end - " << offspring_amount << "\n";
	
	return offspring_containers;
}

double GeneticAlgoForSizesClass::normal_dist(double x)
{
	double sigma = 10;
	double mju = 0;
	
	return 1/(sigma*sqrt(2*M_PI))*exp(-pow(x - mju, 2)/2/pow(sigma, 2));
}

double** GeneticAlgoForSizesClass::compute_cell_sizes(container** con)
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
			for (std::vector<double>::iterator it = v.begin(); it != v.end(); ++it){
				double x = *it;
				double y = *(++it);
				double z = *(++it);
				
				for (std::vector<double>::iterator it2 = v.begin(); it2 != v.end(); ++it2){
					double x1 = *it2;
					double y1 = *(++it2);
					double z1 = *(++it2);
				
					double dist = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1) + (z - z1)*(z - z1));
						
					// std::cout << "tes vertex " << dist << " = " << x << " - " << x1 << " " << y << " - " << y1 << " " << z << " - " << z1 << "\n"; 
					// std::cout << dist << " " << "" << "\n"; 
					if (dist > max) max = dist;
				}
			}
			// std::cout << "real size before\n";
			real_sizes[i][loop_counter++] = max;
			// std::cout << "real size after\n";
		} while (cl.inc());
	}
	
	return real_sizes;
}


void GeneticAlgoForSizesClass::output_data(std::string filename, double* size_dist)
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
	for (std::map<double, double>::iterator it=current_distribution.begin(); it!=current_distribution.end(); ++it)
	{
		myfile << it->first << "  " << it->second << "  " << original_distribution(it->first) << "\n";
	}
	myfile.close();
}

void GeneticAlgoForSizesClass::write_penalty_step(std::string filename, int penalties, double penalty)
{	
	ofstream myfile;
	myfile.open (filename.c_str(), ofstream::out | ofstream::app);
	if (penalties == population_size)
		myfile << "Step  Penalty\n";
	myfile << penalties << "  " << penalty << "\n";
	myfile.close();
}