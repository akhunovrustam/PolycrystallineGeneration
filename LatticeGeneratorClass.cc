#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "voro++.hh"
#include "LatticeGeneratorClass.h"

using namespace std;
using namespace voro;

LatticeGeneratorClass::LatticeGeneratorClass(int centers){
	con = new container(x_min,x_max,y_min,y_max,z_min,z_max,6,6,6,true,true,true,8);
	particles = centers;
}

void LatticeGeneratorClass::fillRandomlyAndBuildGrains(){
	unsigned int i,j;
	double x,y,z;
	int id,nx,ny,nz;
	voronoicell_neighbor c;
	vector<int> neigh,f_vert;
	vector<double> v;
	planes_grains = new float**[particles];
	
	
	
	// srand (time(NULL));
	
	for(i=0;i<particles;i++) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		con->put(i,x,y,z);
	}
	
	// Loop over all particles in the container and compute each Voronoi cell
	c_loop_all cl(*con);
	int loop_counter = 0;
	if(cl.start()) do if(con->compute_cell(c,cl)) {
		cl.pos(x,y,z);
		id=cl.pid();

		// Gather information about the computed Voronoi cell
		c.neighbors(neigh);
		c.face_vertices(f_vert);
		c.vertices(x,y,z,v);

		// Loop over all faces of the Voronoi cell
		float **planes = new float*[neigh.size()+1];
		planes[0] = new float[1];
		planes[0][0] = neigh.size();
		
		int planes_size = 1;
		for(i=0,j=0;i<neigh.size();i++) {
			calculateGrainPlanes(f_vert,v,j, planes, x, y, z, planes_size);
			planes_size++;
			
			// Skip to the next entry in the face vertex list
			j+=f_vert[j]+1;
		}
		
		planes_grains[loop_counter++] = planes;
		
	} while (cl.inc());
	
	
	// Output the particle positions in gnuplot format
	con->draw_particles("random_points_p.gnu");
	// Output the Voronoi cells in gnuplot format
	con->draw_cells_gnuplot("random_points_v.gnu");
	
	//generate lattice for polycrystal grains
	for (i = 0; i < particles; i++)
	{
		printf("draw_particle inside: i=%i\n", i);
		generateLattice(planes_grains[i]);
	}
}

void LatticeGeneratorClass::calculateGrainPlanes(vector<int> &f_vert,vector<double> &v,int j, float** planes, float x, float y, float z, int planes_size) {
	static char s[6][128];
	int k,l,n=f_vert[j];
	ofstream myfile;
	
	double pnts[3][3];
	for(k=0;k<3;k++) {
		l=3*f_vert[j+k+1];
		pnts[k][0] = v[l];
		pnts[k][1] = v[l+1];
		pnts[k][2] = v[l+2];
	}
	
	// max size from particle to grain vertex
	myfile.open("all_vertex.xyz", std::fstream::in | std::fstream::out | std::fstream::app);
	
	float max = 0;
	for(k=0;k<n;k++) {
		l=3*f_vert[j+k+1];
		myfile << "Al " << v[l] << " " << v[l+1] << " " << v[l+2] << "\n";
	
		float dist = sqrt((v[l]-x)*(v[l]-x) + (v[l+1]-y)*(v[l+1]-y) + (v[l+2]-z)*(v[l+2]-z));
		if (max < dist) max = dist;
	}
	myfile.close();
	
	double a, b, c, d;
	double a1, b1, c1, d1;
	double v1[3], v2[3];
	
	//plane vectors calculation
	v1[0] = pnts[0][0] - pnts[1][0];
	v1[1] = pnts[0][1] - pnts[1][1];
	v1[2] = pnts[0][2] - pnts[1][2];
	
	v2[0] = pnts[0][0] - pnts[2][0];
	v2[1] = pnts[0][1] - pnts[2][1];
	v2[2] = pnts[0][2] - pnts[2][2];
	
	//plane constants calculation
	a = v1[1]*v2[2] - v1[2]*v2[1];
	b = v1[2]*v2[0] - v1[0]*v2[2];
	c = v1[0]*v2[1] - v1[1]*v2[0];
	d = a*pnts[0][0] + b*pnts[0][1] + c*pnts[0][2];

	
	planes[planes_size] = new float[9];
	planes[planes_size][0] = a;
	planes[planes_size][1] = b;
	planes[planes_size][2] = c;
	planes[planes_size][3] = d;
	
	planes[planes_size][4] = a*x + b*y + c*z - d > 0 ? 1 : -1;
	planes[planes_size][5] = x;
	planes[planes_size][6] = y;
	planes[planes_size][7] = z;
	planes[planes_size][8] = max;
}

int LatticeGeneratorClass::generateLattice(float** grain) {
	
	float alpha, beta, gamma;
	int box_dimension = 10;
	// srand (time(NULL));
	
	// float basis[4][3] = { {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0, 0.5}, {0, 0.5, 0.5} };
	float basis[1][3] = { {0.0, 0.0, 0.0}};
	
	int ii = 0;
	box_dimension = 0;
	
	//atoms generation box dimension determination
	for(int ii = 1; ii <= grain[0][0]; ii++)
	{
		if (box_dimension < grain[ii][8])
		{
			box_dimension = ceil(grain[ii][8]);
		}
	}
	
	box_dimension = ceil(box_dimension / LATTICE);
	
	float xc = grain[1][5];
	float yc = grain[1][6];
	float zc = grain[1][7];
	
	//rotation angles
	alpha = rnd() * PI / 2;
	beta = rnd() * PI / 2;
	gamma = rnd() * PI / 2;
	
	
	stringstream buffer;
	ofstream myfile;
	myfile.open("fcc_lattice.xyz", std::fstream::in | std::fstream::out | std::fstream::app);
	
	//generation and cutting off of atoms
	int atoms_quantity = 0;
	for (int i = -box_dimension*2; i < box_dimension*2; i++)
		for (int j = -box_dimension*2; j < box_dimension*2; j++)
			for (int k = -box_dimension*2; k < box_dimension*2; k++){
				int basis_length = sizeof(basis)/sizeof(basis[0]);
				
				for (int m = 0; m < basis_length; m++){
				
					float x, y, z;
				
					//basis atom coords in norotated system
					x = xc + LATTICE * i + basis[m][0] * LATTICE;
					y = yc + LATTICE * j + basis[m][1] * LATTICE;
					z = zc + LATTICE * k + basis[m][2] * LATTICE;

					//rotation OX
					float tmpy, tmpz, tmpx;
					tmpy = y;
					tmpz = z;
					y = cos(alpha) * tmpy - sin(alpha) * tmpz;
					z = sin(alpha) * tmpy + cos(alpha) * tmpz;
					
					//rotation OY
					tmpx = x;
					tmpz = z;
					x = cos(beta) * tmpx + sin(beta) * tmpz;
					z = -sin(beta) * tmpx + cos(beta) * tmpz;
					
					//rotation OZ
					tmpx = x;
					tmpy = y;
					x = cos(gamma) * tmpx - sin(gamma) * tmpy;
					y = sin(gamma) * tmpx + cos(gamma) * tmpy;
					
					//cutting off excess atoms
					bool inside_grain = true;
					for(int ii = 1; ii <= grain[0][0]; ii++)
					{
						int planepos = grain[ii][0]*x + grain[ii][1]*y + grain[ii][2]*z - grain[ii][3] > 0 ? 1 : -1;
						if (planepos != grain[ii][4])
						{
							inside_grain = false;
							break;
						}
					}
					
					
					if (inside_grain)
					{
						//put external atoms inside emulation box
						// for x
						if (x < x_min) x = x_max - (x_min - x);
						if (x > x_max) x = x_min - (x_max - x);
						//for y
						if (y < y_min) y = y_max - (y_min - y);
						if (y > y_max) y = y_min - (y_max - y);
						//for z
						if (z < z_min) z = z_max - (z_min - z);
						if (z > z_max) z = z_min - (z_max - z);
						
						atoms_quantity++;
						buffer << "Al " << x << " " << y << " " << z << "\n";
					}
				}
			}
	
	myfile << atoms_quantity << "\n";
	myfile << "FCC \n";
	myfile << buffer.str();
	
	myfile.close();
}