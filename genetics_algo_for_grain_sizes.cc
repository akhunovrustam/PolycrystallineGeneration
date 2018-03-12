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

static const double x_min=-1,x_max=1;
static const double y_min=-1,y_max=1;
static const double z_min=-1,z_max=1;
	
double rnd() {return double(rand())/RAND_MAX;}

double fitness()
{
	
}

void mutation(container** con)
{
	for (int i = 0; i < population_size; i++){
		double new_points[particles][3];
		bool point_changed = false;
		
		for (int j = 0; j < particles; j++){
			srand (time(NULL));
			int rnd_index1 = ceil(rnd()*ceil(sqrt(population_size)));
			int rnd_index2 = ceil(rnd()*ceil(sqrt(population_size)));
			
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