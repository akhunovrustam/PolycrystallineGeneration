// Polycrystal rebuild from the size distribution using genetics algorithm
//
// Author   : Rustam Akhunov (GUT Gdansk)
// Email    : akhunovrustam@mail.ru
// Date     : 
#define _USE_MATH_DEFINES
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include "voro++.hh"

	
using namespace std;
using namespace voro;
// using namespace LatticeGeneratorClass;

#define penalty_steps 100
#define penalty_step 0.01
#define max_iterations 100
#define max_allowed_penalty 0.01
#define population_size 32
#define particles 100
#define mutation_probability 0.91
#define mutation_max_applitude 0.05
#define crossover_probability 0.05

static const double x_min=-1,x_max=1;
static const double y_min=-1,y_max=1;
static const double z_min=-1,z_max=1;
static const double shift_value = 0.0001;
	
double rnd() {return double(rand())/RAND_MAX;}

double original_distribution(double point)
{
	double sigma = 0.35;
	double mju = 0.0;
	double res = 1/(point*sigma*sqrt(2*M_PI))*exp(-pow(log(point) - mju, 2)/(2*pow(sigma, 2)));
	return res;
}



double fitness_penalty(int points_number, double (*original_distribution)(double), 
	std::map<double, double> current_distribution)
{
	double sum = 0;
	for (std::map<double, double>::iterator it=current_distribution.begin(); it!=current_distribution.end(); ++it){
		// std::cout << "original: " << it->first << " => " << (*original_distribution)(it->first) << "\n";
		sum += pow((*original_distribution)(it->first) - it->second, 2);
	}
	
	return sum / points_number;
}

double size_penalty(double* size_dist)
{
	std::map<double, double> current_distribution;
	double max = 0, min = 100;
	
	double avg = 0;
	
	for (int i = 0; i < particles; i++)
		avg += size_dist[i];
	
	avg = avg / particles;
	
	double k = 1 / penalty_step / particles;
	// for (int i = 0; i < particles; i++)
	// {
		// if (size_dist[i] > max) max = size_dist[i];
		// if (size_dist[i] < min) min = size_dist[i];
	// }
	// double step = (max - min)/penalty_points;
	
	// for (int i = 0; i < particles; i++)
	// {
		// for (int j = 0; j < penalty_points; j++)
			// if (size_dist[i] > min + j*step && size_dist[i] < min + (j+1)*step)
				// if (current_distribution.find(min + (j + 0.5) * step) == current_distribution.end())
					// current_distribution[min + (j + 0.5) * step] = 1;
				// else current_distribution[min + (j + 0.5) * step] = current_distribution[min + (j + 0.5) * step] + 1;
	// }
	for (int i = 0; i < penalty_steps; i++)
	{
		current_distribution[i*penalty_step + penalty_step/2] = 0;
		for (int j = 0; j < particles; j++){
			// std::cout << "size: " << size_dist[j] << " i: " << i*penalty_step << "\n";
			if (size_dist[j] > i*penalty_step && size_dist[j] < (i+1)*penalty_step)
				current_distribution[i*penalty_step + penalty_step/2] = current_distribution[i*penalty_step + penalty_step/2] + k;
		}
	}
	
	// for (std::map<double, double>::iterator it=current_distribution.begin(); it!=current_distribution.end(); ++it)
		// std::cout << "dist: " << it->first << " => " << it->second << "\n";
	
	double ret = fitness_penalty(penalty_steps, original_distribution, current_distribution);
	// std:cout << "penalty " << ret << "\n";
	return ret;
}

void mutation(container ***con1)
{
	container **con = *con1;
	for (int i = 0; i < population_size; i++){
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
					new_points[new_points_cntr][0] = x + rnd()*mutation_max_applitude*(rnd() > 0.5 ? -1 : 1);
					new_points[new_points_cntr][1] = y + rnd()*mutation_max_applitude*(rnd() > 0.5 ? -1 : 1);
					new_points[new_points_cntr][2] = z + rnd()*mutation_max_applitude*(rnd() > 0.5 ? -1 : 1);
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
			con[i] = con_subst;
		}
	}
}

int tournament_selection(double** size_dist, int selected = -1)
{
	int index1 = rnd()*population_size;
	int index2 = rnd()*population_size;
	if (index1 == index2 && index2 == population_size - 1) index1--;
	if (index1 == index2) index1++;
	
	double penalty1 = size_penalty(size_dist[index1]);
	double penalty2 = size_penalty(size_dist[index2]);
	
	if ((penalty1 < penalty2 && selected == -1) || (selected != -1 && selected == index2)){
		// std::cout << "return1 " << index1 << " - " << index2 << "\n";
		return index1; 
	}
	else {
		// std::cout << "return2 " << index2 << "\n";
		return index2;
	}
}
/* 
container** interchage_halves(container* ind1, container* ind2, k, b)
{
	container *offspring1 = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
	container *offspring2 = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
	
	for(int j = 0; j < particles; j++) {
		//here might be problems because of =
		if (k*ind1->p[j][0] + b - ind1->p[j][1] >= 0)
			offspring1[i]->put(j,ind1->p[j][0],ind1->p[j][1],ind1->p[j][2]);
		else offspring2[i]->put(j,ind1->p[j][0],ind1->p[j][1],ind1->p[j][2]);
		
		//here might be problems because of =
		if (k*ind2->p[j][0] + b - ind2->p[j][1] >= 0)
			offspring1[i]->put(j,ind2->p[j][0],ind2->p[j][1],ind2->p[j][2]);
		else offspring2[i]->put(j,ind2->p[j][0],ind2->p[j][1],ind2->p[j][2]);
	}
	
	return [offspring1, offspring2];
}
 */
struct point_for_crossover {
	double x;
	double y;
	double z;
	int block;
	int particle;
	bool is_interchanged;
} ;

container** crossover_by_mapping(container** con, double** size_dist)
{
	int offspring_amount = 0;
	voronoicell_neighbor c1;
	container** offspring_containers = new container*[population_size + 1];
		
	while (offspring_amount < population_size)
	{
		// std::cout << "cross1" << "\n";
		bool cannot_crossover = false;
		int index1 = tournament_selection(size_dist);
		// std::cout << "cross2" << "\n";
		container* ind1 = con[index1];
		int index2 = tournament_selection(size_dist, index1);
		container* ind2 = con[index2];
		// std::cout << "cross2" << "\n";
		//making mapping ID->xyz_bool where bool is if point goes to ind1 or not
		std::map<int, point_for_crossover> id1_to_coords;
		for (int i = 0; i < ind1->nx*ind1->ny*ind1->nz; i++)
		{
			// printf("points per box - %i - %i\n", i, ind1->co[i]);
			for (int j = 0; j < ind1->co[i]; j++)
			{
				double *pp=ind1->p[i]+3*j;
				point_for_crossover tmp_point = {*(pp++), *(pp++), *(pp++), i, j, false};
				id1_to_coords[ind1->id[i][j]] = tmp_point;
			}
		}
		// for (std::map<int, point_for_crossover>::iterator it=id1_to_coords.begin(); it!=id1_to_coords.end(); ++it)
			// std::cout << it->first << " => " << it->second.block << " => " << it->second.particle << '\n';
				
		std::map<int, point_for_crossover> id2_to_coords;
		for (int i = 0; i < ind2->nx*ind2->ny*ind2->nz; i++)
			for (int j = 0; j < ind2->co[i]; j++)
			{
				double *pp=ind2->p[i]+3*j;
				point_for_crossover tmp_point = {*(pp++), *(pp++), *(pp++), i, j, false};
				id2_to_coords[ind2->id[i][j]] = tmp_point;
			}
				
		
		int mapping[particles];
		double dist = 100;
		//making mapping from particels of ind1 to ind2
		// std::cout << "cross3" << "\n";
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
						// std::cout << "cur dist " << cur_dist << " size: " << id_to_dist.size() << "\n";
					}
				// std::cout << "dist map size " << id_to_dist.size() << "\n";
				
				
				for (std::map<double, int>::iterator it=id_to_dist.begin(); it!=id_to_dist.end(); ++it)
				{
					bool dist_exist = false;
					for (std::map<int, int>::iterator it2=id_to_id.begin(); it2!=id_to_id.end(); ++it2)
						if (it2->second == it->second)
						{
							// nearest_id = it->second;
							dist_exist = true;
							break;
						}
					
					if (dist_exist) continue;
					
					// std::cout << "dist exist: " << dist_exist << "\n";
					
					id_to_id[ind1->id[i][j]] = it->second;
					break;
				}
				
				// if (dist > cur_dist)
				// {
					// dist = cur_dist;
					// if (id_to_id.find(ind1->id[i][j]) != id_to_id.end()) break;
					// id_to_id[ind1->id[i][j]] = ind2->id[h][k];
				// }
		
			}
		// std::cout << "ID to ID size: " << id_to_id.size() << '\n';
		// for (std::map<int, int>::iterator it=id_to_id.begin(); it!=id_to_id.end(); ++it)
			// std::cout << "ID to ID " << it->first << " => " << it->second << " " << id_to_id.size() << '\n';
		
		int first_particles = ceil(rnd()*particles);
		first_particles = first_particles == 0 ? 1 : first_particles;
		int selected_particles = 0;
		
		
		voronoicell_neighbor c;
		vector<int> neigh;
		int id = id1_to_coords.begin()->first;
		
		int *points;
		points = new int[first_particles];
		points[selected_particles] = id;
		id1_to_coords[id].is_interchanged = true;
		
		selected_particles++;
		int checked_particles = 0;
		// std::cout << "first particle " << first_particles << "\n";
		while (selected_particles < first_particles){
			if(ind1->compute_cell(c, id1_to_coords[points[checked_particles]].block, 
				id1_to_coords[points[checked_particles]].particle)) {
				
				
				checked_particles++;
				
				c.neighbors(neigh);
				
				for (std::vector<int>::iterator it = neigh.begin() ; it != neigh.end(); ++it)
					if (selected_particles < first_particles) 
					{
						// printf("selected_particels %i\n", selected_particles);
						points[selected_particles++] = *it;
						id1_to_coords[*it].is_interchanged = true;
					}
					else break;
					
			}
			else {
				printf("FAIL - point_id: %i, partic_per_block: %i, block: %i, partic_index: %i, xyz: %f %f %f\n", points[checked_particles], ind1->co[id1_to_coords[points[checked_particles]].block], 
					id1_to_coords[points[checked_particles]].block, id1_to_coords[points[checked_particles]].particle,
					id1_to_coords[points[checked_particles]].x, id1_to_coords[points[checked_particles]].y,
					id1_to_coords[points[checked_particles]].z);
				break;
			}
		}
		
		container *offspring1 = NULL;
		container *offspring2 = NULL;
		
		offspring1 = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
		// std::cout << "cross5.1" << "\n";
		offspring2 = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
		// std::cout << "cross5.2" << "\n";
		
		for (std::map<int, point_for_crossover>::iterator it=id1_to_coords.begin(); it!=id1_to_coords.end(); ++it){
				
			//here might be problems because of =
			if (it->second.is_interchanged)
			{
				// std::cout << "new point1-1 " << it->second.x << " " << it->second.y << " " << it->second.z << "\n";
				// std::cout << "new point1-2 " << id2_to_coords[id_to_id[j]].x << " " << id2_to_coords[id_to_id[j]].y << " " << id2_to_coords[id_to_id[j]].z << "\n";
				offspring1->put(it->first, it->second.x, it->second.y, it->second.z);
				offspring2->put(it->first, id2_to_coords[id_to_id[it->first]].x, id2_to_coords[id_to_id[it->first]].y, 
					id2_to_coords[id_to_id[it->first]].z);
			}
			else {
				// std::cout << "new point2-1 " << id2_to_coords[id_to_id[j]].x << " " << id2_to_coords[id_to_id[j]].y << " " 
					// << id2_to_coords[id_to_id[j]].z << "\n";
				// std::cout << "new point2-2 " << it->second.x << " " << it->second.y << " " << it->second.z << "\n";
				offspring2->put(it->first, it->second.x, it->second.y, it->second.z);
				offspring1->put(it->first, id2_to_coords[id_to_id[it->first]].x, id2_to_coords[id_to_id[it->first]].y, 
					id2_to_coords[id_to_id[it->first]].z);
			}
		}
		// std::cout << "cross6 - " << offspring_amount << "\n";
		
		offspring_containers[offspring_amount++] = offspring1;
		offspring_containers[offspring_amount++] = offspring2;
		// std::cout << "cross7 - " << offspring_amount << "\n";
	}
	
	// std::cout << "cross end - " << offspring_amount << "\n";
	
	return offspring_containers;
}

double normal_dist(double x){
	double sigma = 10;
	double mju = 0;
	
	return 1/(sigma*sqrt(2*M_PI))*exp(-pow(x - mju, 2)/2/pow(sigma, 2));
}

double** compute_cell_sizes(container** con)
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
			// printf("computed cell %i\n", con[i]->compute_cell(c,cl));
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
			
			real_sizes[i][loop_counter++] = max;
		} while (cl.inc());
	}
	
	return real_sizes;
}


void output_data(std::string filename, double* size_dist)
{
	double avg = 0;
	std::map<double, double> current_distribution;
	
	for (int i = 0; i < particles; i++)
		avg += size_dist[i];
	
	avg = avg / particles;
	
	double coef = 1 / penalty_step / particles;
	for (int i = 0; i < penalty_steps; i++)
	{
		current_distribution[i*penalty_step + penalty_step/2] = 0;
		for (int j = 0; j < particles; j++){
			// std::cout << "size: " << size_dist[j] << " i: " << i*penalty_step << "\n";
			if (size_dist[j] > i*penalty_step && size_dist[j] < (i+1)*penalty_step)
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

int main() {
	srand (time(NULL));
	container **con;
	con = new container*[population_size];
	double x, y, z;
	double **real_sizes;
	double penalty[population_size];
	double min_penalty = 100;
	int iterations = 0;
	
	
	//initial population
	for (int i = 0; i < population_size; i++){
		con[i] = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
		
		for(int j = 0; j < particles; j++) {
			x=x_min+rnd()*(x_max-x_min);
			y=y_min+rnd()*(y_max-y_min);
			z=z_min+rnd()*(z_max-z_min);
			con[i]->put(j,x,y,z);
			// std::cout << "x - " << x << " y - " << y << " z - " << z << "\n";
		}
		
	}
	
	real_sizes = compute_cell_sizes(con);
	
	// for (int i = 0; i < population_size; i++){
		// for (int j = 0; j < particles; j++){
			// printf("size of tessellation no %i: %f\n", i, real_sizes[i][j]);
		// }
	// }

	//iterate until reach max iterations or precision
	while (true){
		int min_penalty_index = -1;
		min_penalty = 100;
		for (int i = 0; i < population_size; i++){
			penalty[i] = size_penalty(real_sizes[i]);
			if (penalty[i] == 0.0)
			{
				for (int m = 0; m < particles; m++)
					std::cout << "strange size " << real_sizes[i][m] << "\n";
			}
		
			if (min_penalty > penalty[i]) {
				min_penalty = penalty[i];
				min_penalty_index = i;
			}
		}
		
		if (iterations == 0){
			output_data("dist_first.txt", real_sizes[min_penalty_index]);
			
			// return 0;
		}
		std::cout << "penalty " << min_penalty << " ===========================================\n";
		// printf("penalty %f  -- max pen %f\n", min_penalty, max_allowed_penalty);
		if (min_penalty < max_allowed_penalty) break;
		
		iterations++;
		std::cout << "iter " << iterations << "\n";
		
		if (iterations > max_iterations) {
			// std::cout << "end of iter " << "\n";
			output_data("dist_end.txt", real_sizes[min_penalty_index]);
	
			break;
		}
		
		
		std::cout << "crossover\n";
		con = crossover_by_mapping(con, real_sizes);
		
		std::cout << "mutation\n";
		mutation(&con);
		
		real_sizes = compute_cell_sizes(con);
		
		// for (int i = 0; i < con[0]->nx*con[0]->ny*con[0]->nz; i++)
			// for (int j = 0; j < con[0]->co[i]; j++)
			// {
				// double *pp=con[0]->p[i]+3*j;
				// std::cout << "points: " << *(pp++) << " " << *(pp++) << " " << *(pp++) << "\n";
			// }
	}
}
