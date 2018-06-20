#include "GeneticAlgoForSizesClass.hh"
#include <random>
#include <thread> 
#include <chrono>
#include <iterator>

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
		// sum += abs((*original_distribution)(it->first) - it->second);
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
			if ((size_dist[j] / avg) > i*penalty_step && (size_dist[j] / avg) < (i+1)*penalty_step)
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

double GeneticAlgoForSizesClass::mutate(double val, double min, double max, bool adopted_shift, int dist_num, double mult)
{
	double mutation_max_applitude_e = 1*pow(pow(half_boxside*2, 3)/particles, 1.0/3.0) / 2;
	
	int max_iter = 10;
	int iter = 0;
	do {
		iter++;
		if (iter > max_iter) {
			val = reinit(min, max);
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

void GeneticAlgoForSizesClass::mutation(container ***con1, int iteration, bool reinit_flag, double mutation_probability, 
	int dist_num, sorted_points** exp, double mult)
{
	container **con = *con1;
	sorted_points* experiment = *exp;
	for (int i = surviving_size; i < population_size; i++){
		double new_points[particles][3];
		bool point_changed = false;
		
		sorted_points pt = experiment[i];
		int new_points_cntr = 0;
		for (sorted_points::iterator it = pt.begin(); it != pt.end(); it++)
			for (map<double, map<double, int>>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
				for (map<double, int>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++)
				{
					double x = it->first;
					double y = it2->first;
					double z = it3->first;
					
					double rnd_num = rnd();
					
					if (rnd_num < mutation_probability){
						point_changed = true;
						if (reinit_flag == false){
							new_points[new_points_cntr][0] = mutate(x, x_min, x_max, true, dist_num, mult);
							new_points[new_points_cntr][1] = mutate(y, y_min, y_max, true, dist_num, mult);
							new_points[new_points_cntr][2] = mutate(z, z_min, z_max, true, dist_num, mult);
							
							if (iteration < 20)
							{
								double null_rnd = rnd();
								if (null_rnd < 0.1)
									if (new_points[new_points_cntr][0] != 0 && new_points[new_points_cntr][1] != 0 
										&& new_points[new_points_cntr][2] != 0)
										new_points[new_points_cntr][(int)floor(rnd()*3)] = 0;
							}
						
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
		
		
		if (point_changed){
			container *con_subst;
			con_subst = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
			
			experiment[i].clear();
			for(int j = 0; j < particles; j++) {
				con_subst->put(j,new_points[j][0],new_points[j][1],new_points[j][2]);
				experiment[i][new_points[j][0]][new_points[j][1]][new_points[j][2]] = j;
			}
			delete con[i];
			con[i] = con_subst;
		}
	}
}

int GeneticAlgoForSizesClass::tournament_selection(double* penalty, int iter, int offspring_amount, int selected)
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

std::map<int, int> GeneticAlgoForSizesClass::ind_to_ind(int ind1, int ind2, sorted_points** exp)
{
	sorted_points* experiment = *exp;
	// cout << "size ind to ind: " << experiment[ind1].size() << endl;
	// cout << "size ind to ind: " << experiment[ind2].size() << endl;
	int expcnt = 0;
	double amp = 20 * pow(pow(half_boxside*2, 3)/particles, 1.0/3.0) / 2;
	
	// cout << "amp: " << amp << endl;
	
	map<int, int> pairs;
	map<int, int> busy;
	int c1 = 0, c2 = 0, c3 = 0;
	for (map<double, map<double, map<double, int>>>::iterator it = experiment[ind1].begin(); it != experiment[ind1].end(); it++)
	{
		double x = it->first;
		double xmin = x - amp;
		double xmax = x + amp;
		for (map<double, map<double, int>>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
		{
			double y = it2->first;
			double ymin = y - amp;
			double ymax = y + amp;
		
			for (map<double, int>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++)
			{
				double z = it3->first;
				double zmin = z - amp;
				double zmax = z + amp;
				int id1 = it3->second;
				
				double mindist = 10000;
				int lastid;
				for (map<double, map<double, map<double, int>>>::iterator itt = experiment[ind2].begin(); itt != experiment[ind2].end(); itt++)
				{
					double xt = itt->first;
					// cout << xt << " " << xmax << " " << xmin << endl;
					if (xt > xmax) break;
					if (xt > xmin)
					{
						c1++;
						for (map<double, map<double, int>>::iterator itt2 = itt->second.begin(); itt2 != itt->second.end(); itt2++)
						{
							double yt = itt2->first;
							if (yt > ymax) break;
							if (yt > ymin)
							{
								c2++;
								for (map<double, int>::iterator itt3 = itt2->second.begin(); itt3 != itt2->second.end(); itt3++)
								{
									double zt = itt3->first;
									int id2 = itt3->second;
									// cout << "id2: " << id2 << " id1: " << id1 << " " << busy[id2] << endl;
									if (pairs[busy[id2 + 12345]] == id2 + 12345) continue;
									if (zt > zmax) break;
									if (zt > zmin)
									{
										c3++;
										expcnt++;
										double dist = abs(x - xt) + abs(y - yt) + abs(z - zt);
										if (mindist > dist)
										{
											mindist = dist;
											lastid = id2;
											pairs[id1] = lastid + 12345;
											busy[lastid + 12345] = id1;
											
										}
										
										// cout << x << " " << y << " " << z << endl;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	for (map<int, int>::iterator it = pairs.begin(); it != pairs.end(); it++)
		pairs[it->first] = pairs[it->first] - 12345;
	
	map<int, int> bla;
	for (map<int, int>::iterator it = pairs.begin(); it != pairs.end(); it++)
		if (bla[it->second] > 1)
		{
			cout << it->first << " - " << it->second << endl;
			exit(0);
		}
		else
		{
			bla[it->second]++;
		}
	
	// cout << c1 << " " << c2 << " " << c3 << " " << pairs.size() << endl;
	// exit(0);
	return pairs;
}

std::map<int, point_for_crossover> GeneticAlgoForSizesClass::ind_points_map(sorted_points sortp)
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

void GeneticAlgoForSizesClass::select_interchange_regions(container* ind1, std::map<int, point_for_crossover> *id1_to_coords, 
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

void GeneticAlgoForSizesClass::select_interchange_randomly(std::map<int, point_for_crossover> *id1_to_coords)
{
	for (std::map<int, point_for_crossover>::iterator it = (*id1_to_coords).begin() ; it != (*id1_to_coords).end(); ++it)
		if (rnd() < 0.5)
			(*id1_to_coords)[it->first].is_interchanged = true;
}

container** GeneticAlgoForSizesClass::crossover_by_mapping(container** con, double* penalty, int iter, int crossover_points, 
	sorted_points** exp, sorted_points** exp_off, config cfg, neighbors* neg)
{
	sorted_points* experiment = *exp;
	sorted_points* experiment_tmp = *exp_off;
	
	time_point<system_clock> stop1, stop2;
	
	int offspring_amount = surviving_size;
	voronoicell_neighbor c1;
	container** offspring_containers = new container*[population_size + 1];
	
	std::map<double, int> penalty_rating;
	for (int i = 0; i < population_size; i++)
		penalty_rating[penalty[i]] = i;
	
	
	int survived = 0;
	for (std::map<double, int>::iterator it=penalty_rating.begin(); it!=penalty_rating.end(); ++it)
		if (survived < surviving_size)
			offspring_containers[survived++] = con[it->second];
	
	penalty_rating.clear();
	
	while (offspring_amount < population_size * cfg.ps_fac)
	{
		
		//tournament select of 2 individuals for a crossover
		int index1 = tournament_selection(penalty, iter, offspring_amount);
		
		container* ind1 = con[index1];
		int index2 = tournament_selection(penalty, iter, offspring_amount, index1);
		container* ind2 = con[index2];
		
		//making mapping ID => xyz_bool where bool is if point goes to ind1 or not
		std::map<int, point_for_crossover> id1_to_coords = ind_points_map(experiment[index1]);
		std::map<int, point_for_crossover> id2_to_coords = ind_points_map(experiment[index2]);
		
		//making mapping from particels of ind1 to ind2
		std::map<int, int> id_to_id = ind_to_ind(index1, index2, exp);
		
		// cout << "pre regions select" << endl;
		if (crossover_points != -1)
		{
			auto item = id1_to_coords.begin();
			advance( item, rnd() * id1_to_coords.size() );
			
			// cout << "first select" << endl;
			select_interchange_regions(ind1, &id1_to_coords, item->first, cfg.from, cfg.to, neg[index1]);
			if (crossover_points == 2)
			{
				auto item = id1_to_coords.begin();
				advance( item, rnd() * id1_to_coords.size() );
			
				// cout << "second select" << endl;
				select_interchange_regions(ind1, &id1_to_coords, item->first, cfg.from, cfg.to, neg[index1]);
			}
		} else {
			select_interchange_randomly(&id1_to_coords);
		}
		
		// cout << "post regions select" << endl;
		container *offspring1 = NULL;
		container *offspring2 = NULL;
		
		offspring1 = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
		// std::cout << "cross5.1" << "\n";
		offspring2 = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
		// std::cout << "cross5.2" << "\n";
		
		std::map<double, std::map<double, double> > points_map1;
		std::map<double, std::map<double, double> > points_map2;
		experiment_tmp[offspring_amount].clear();
		experiment_tmp[offspring_amount + 1].clear();
		// std::cout << "before copy of points\n";
		for (std::map<int, point_for_crossover>::iterator it=id1_to_coords.begin(); it!=id1_to_coords.end(); ++it){
				
			//here might be problems because of =
			double x = it->second.x, y = it->second.y, z = it->second.z; 
			double x2 = id2_to_coords[id_to_id[it->first]].x, y2 = id2_to_coords[id_to_id[it->first]].y, z2 = id2_to_coords[id_to_id[it->first]].z; 
			// if (x == 0 || y == 0 || z == 0 || x2 == 0 || y2 == 0 || z2 == 0)
			// {
				// cout << it->first << " " << id_to_id[it->first] << " indices: " << index1 << " " << index2 << endl;
				// cout << x << " " << y << " " << z << " " << x2 << " " << y2 << " " << z2 << endl;
				// exit(0);
				
			// }
			if (it->second.is_interchanged)
			{
				// std::cout << "new point1-1 " << x << " " << y << " " << z << " points map: " << points_map1[x][y] << "\n";
				// std::cout << "new point1-2 " << x2 << " " << y2 << " " << z2 << "\n";
				if (points_map1[x][y] != z)
				{
					points_map1[x][y] = z;
					offspring1->put(it->first, x, y, z);
					experiment_tmp[offspring_amount][x][y][z] = it->first;
				} else {
					double shifted_x = mutate(x, x_min, x_max);
					points_map1[shifted_x][y] = z;
					offspring1->put(it->first + 10000, shifted_x, y, z);
					experiment_tmp[offspring_amount][shifted_x][y][z] = it->first + 10000;
					// cout << "bad id: " << it->first << ": " << x << " " << y << " " << z << endl;
					// cout << index1 << " " << index2 << endl;
				}
				
				if (points_map2[x2][y2] != z2)
				{
					points_map2[x2][y2] = z2;
					offspring2->put(it->first, x2, y2, z2);
					experiment_tmp[offspring_amount + 1][x2][y2][z2] = it->first;
				} else {
					double shifted_x = mutate(x2, x_min, x_max);
					points_map2[shifted_x][y2] = z2;
					offspring2->put(it->first, shifted_x, y2, z2);
					experiment_tmp[offspring_amount + 1][shifted_x][y2][z2] = it->first;
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
					experiment_tmp[offspring_amount][x2][y2][z2] = it->first;
				} else {
					double shifted_x = mutate(x2, x_min, x_max);
					// std::cout << "mutated x " << shifted_x << "\n";
					points_map1[shifted_x][y2] = z2;
					offspring1->put(it->first + 10000, shifted_x, y2, z2);
					experiment_tmp[offspring_amount][shifted_x][y2][z2] = it->first + 10000;
					// cout << "bad id2: " << it->first << ": " << x2 << " " << y2 << " " << z2 << endl;
					
					// cout << index1 << " " << index2 << endl;
				}
				
				if (points_map2[x][y] != z)
				{
					points_map2[x][y] = z;
					offspring2->put(it->first, x, y, z);
					experiment_tmp[offspring_amount + 1][x][y][z] = it->first;
				} else {
					double shifted_x = mutate(x, x_min, x_max);
					points_map2[shifted_x][y] = z;
					offspring2->put(it->first, shifted_x, y, z);
					experiment_tmp[offspring_amount + 1][shifted_x][y][z] = it->first;
				}
				// offspring2->put(it->first, x, y, z);
				// offspring1->put(it->first, x2, y2, z2);
			}
		}
		// cout << "post points copy" << endl;
		
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

double** GeneticAlgoForSizesClass::compute_cell_sizes(container** con, neighbors** neg)
{
	neighbors* neg_t = *neg;
	double x, y, z;
	int id;
	
	double **real_sizes;
	real_sizes = new double*[population_size];
	for (int i = 0; i < population_size; i++)
		real_sizes[i] = new double[particles];
	for (int i = 0; i < population_size; i++){
		voronoicell_neighbor c;
		vector<double> v;
		vector<int> neigh;
		
		c_loop_all cl(*(con[i]));
		int loop_counter = 0;
		if(cl.start()) do if(con[i]->compute_cell(c,cl)) {
			// printf("computed individual %i\n", i);
			cl.pos(x,y,z);
			id=cl.pid();
		
			// Gather information about the computed Voronoi cell
			c.vertices(x,y,z,v);
			c.neighbors(neigh);
			
			neg_t[i][id] = neigh;
			
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
}

void GeneticAlgoForSizesClass::write_penalty_step(std::string filename, double penalties, double penalty)
{	
	ofstream myfile;
	myfile.open (filename.c_str(), ofstream::out | ofstream::app);
	if (penalties == population_size)
		myfile << "Step  Penalty\n";
	myfile << penalties << "  " << penalty << "\n";
	myfile.close();
}