// Polycrystal rebuild from the size distribution using genetics algorithm
//
// Author   : Rustam Akhunov (GUT Gdansk)
// Email    : akhunovrustam@mail.ru
// Date     : 

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
#include "LatticeGeneratorClass.cc"

using namespace std;
using namespace voro;
// using namespace LatticeGeneratorClass;

#define population_size 2
#define particles 3
#define mutation_probability 0.01
#define mutation_max_applitude 0.05
#define crossover_probability 0.05

static const double x_min=-1,x_max=1;
static const double y_min=-1,y_max=1;
static const double z_min=-1,z_max=1;
static const double shift_value = 0.0001;
	
double rnd() {return double(rand())/RAND_MAX;}

double fitness_penalty_wrap(container* cont)
{
	return 1.0;
}

double fitness_penalty(int points_number, double (*original_distribution)(double), 
	std::map<double, double> current_distribution)
{
	double sum = 0;
	for (std::map<double, double>::iterator it=current_distribution.begin(); it!=current_distribution.end(); ++it)
		sum += pow((*original_distribution)(it->first) - it->second, 2);
	
	return sum / points_number;
}

void mutation(container** con)
{
	for (int i = 0; i < population_size; i++){
		double new_points[particles][3];
		bool point_changed = false;
		
		for (int j = 0; j < particles; j++){
			srand (time(NULL));
			int rnd_index1 = ceil(rnd()*ceil(sqrt(1/mutation_probability)));
			int rnd_index2 = ceil(rnd()*ceil(sqrt(1/mutation_probability)));
			
			if (rnd_index1 == rnd_index2){
				point_changed = true;
				new_points[j][0] = con[i]->p[j][0] + rnd()*mutation_max_applitude*(rnd() > 0.5 ? -1 : 1);
				new_points[j][1] = con[i]->p[j][1] + rnd()*mutation_max_applitude*(rnd() > 0.5 ? -1 : 1);
				new_points[j][2] = con[i]->p[j][2] + rnd()*mutation_max_applitude*(rnd() > 0.5 ? -1 : 1);
			}
			else {
				new_points[j][0] = con[i]->p[j][0];
				new_points[j][1] = con[i]->p[j][1];
				new_points[j][2] = con[i]->p[j][2];
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

container* tournament_selection(container** con)
{
	int index1 = rnd()*population_size;
	int index2 = rnd()*population_size;
	double penalty1 = fitness_penalty_wrap(con[index1]);
	double penalty2 = fitness_penalty_wrap(con[index2]);
	
	if (penalty1 < penalty2) return con[index1];
	return con[index2];
}

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

container** crossover(container** con)
{
	// std::map<int,bool> parent_pool_index;
	int offspring_amount = 0;
	container** offspring_containers = new container*[population_size];
	
	while (offspring_amount < population_size)
	{
		bool cannot_crossover = false;
		container* ind1 = tournament_selection(con);
		container* ind2 = tournament_selection(con);
		// while (parent_pool_index.size() < 10) int i = 1;
		// mymap['b'] = 100;
		// mymap['a'] = 200;
		// mymap['c'] = 300;

		double sum1;
		for (int i = 0; i < particles; i++)
			sum1 += ind1->p[i][0];
		
		//get avg line for ox axis for first individual
		double avg1 = sum1 / particles;
		
		double sum2;
		for (int i = 0; i < particles; i++)
			sum2 += ind2->p[i][0];
		
		double avg2 = sum2 / particles;
		
		//count particles for one-half of avg (see above) line for the second individual
		int particles_upper = 0;
		for (int i = 0; i < particles; i++)
			if (ind2->p[i][0] < avg1) particles_upper++;
		
		int half_particles = particles / 2;
		double min_dist = 100;
		double adj_point[3];
		int adj_index = -1;
		
		//if amount of particles for second individual differs to next steps
		if (half_particles > particles_upper)
		{
			for (int i = 0; i < particles; i++)
				if (avg1 - ind1->p[i][0] < min_dist) 
					if (avg1 - ind1->p[i][0] < min_dist)
					{
						min_dist = avg1 - ind1->p[i][0];
						adj_point = ind1->p[i];
						adj_index = i;
					}
			
			//shift avg line for the max distance still keeping the same division for the individual 
			avg1 += min_dist - shift_value;
			
			particles_upper = 0;
			for (int i = 0; i < particles; i++)
				if (ind2->p[i][0] < avg1) particles_upper++;
			
			double pos_tan = 100;
			double neg_tan = -100;
			
			//if amount of particles still differs in the same direction
			if (half_particles > particles_upper)
			{
				for (int i = 0; i < particles; i++)
				{
					if (i != adj_index)
					{
						double k = (adj_point[1] - ind1->p[i][1])/(adj_point[0] - ind1->p[i][0]);
						if (k > 0)
						{
							if (pos_tan > k) pos_tan = k;
						}
						else if (neg_tan < k) neg_tan = k;
					}
				}
				
				//max allowed positive and negative inclination
				pos_tan -= shift_value;
				neg_tan += shift_value;
				
				//lean the avg line
				double b = adj_point[1] - pos_tan*adj_point[0];
				int up_sign = pos_tan*(adj_point[0]-shift_value) + b - adj_point[1];
				
				particles_upper = 0;
				for (int i = 0; i < particles; i++)
				{
					if ((pos_tan*ind2->p[i][0] + b - adj_point[1])*up_sign > 0) particles_upper++;
				}
				
				if (half_particles > particles_upper)
				{
					double b = adj_point[1] - neg_tan*adj_point[0];
					int up_sign = neg_tan*(adj_point[0]-shift_value) + b - adj_point[1];
					
					particles_upper = 0;
					for (int i = 0; i < particles; i++)
					{
						if ((neg_tan*ind2->p[i][0] + b - adj_point[1])*up_sign > 0) particles_upper++;
					}
					
					if (half_particles > particles_upper) cannot_crossover = true;
					else {
						if (half_particles == particles_upper)
						{
							//interchange halves of the choosen individuals
							container** offspring = interchage_halves(ind1, ind2);
							offspring_containers[offspring_amount++] = offspring[0];
							offspring_containers[offspring_amount++] = offspring[1];
							
						}
					}
				}
			}
		}
	}
	
	
	  // show content:
	  // for (std::map<char,int>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
		// std::cout << it->first << " => " << it->second << '\n';

	  // return 0;


}

int main() {
	// srand (time(NULL));
	container **con;
	con = new container*[population_size];
	double x, y, z;
	double real_sizes[population_size][particles];
	int id;
	
	for (int i = 0; i < population_size; i++){
		con[i] = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
		voronoicell_neighbor c;
		vector<int> neigh,f_vert;
		vector<double> v;
		
		for(int j = 0; j < particles; j++) {
			x=x_min+rnd()*(x_max-x_min);
			y=y_min+rnd()*(y_max-y_min);
			z=z_min+rnd()*(z_max-z_min);
			con[i]->put(j,x,y,z);
		}
		
		c_loop_all cl(*(con[i]));
		int loop_counter = 0;
		if(cl.start()) do if(con[i]->compute_cell(c,cl)) {
			cl.pos(x,y,z);
			id=cl.pid();

			// Gather information about the computed Voronoi cell
			c.vertices(x,y,z,v);
	
			int planes_size = 1;
			double max = 0;
			for(int k = 0; k < c.current_vertices; k++) {
				for(int j = 0; j < c.current_vertices; j++) {
					double dist = sqrt((v[3*k] - v[3*j])*(v[3*k] - v[3*j]) + (v[3*k + 1] - v[3*j + 1])*(v[3*k + 1] - v[3*j + 1])
						+ (v[3*k + 2] - v[3*j + 2])*(v[3*k + 2] - v[3*j + 2]));
						
					if (dist > max) max = dist;
				}
			}
			
			real_sizes[i][loop_counter++] = max;
		} while (cl.inc());
		
	}
	
	for (int i = 0; i < population_size; i++){
		for (int j = 0; j < particles; j++){
			printf("size of tessellation no %i: %f\n", i, real_sizes[i][j]);
		}
	}
}