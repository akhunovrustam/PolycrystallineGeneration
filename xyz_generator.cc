#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>

#define PI 3.14159265
#define LATTICE 1.1
#define SIZE 2

using namespace std;

double rnd() {return double(rand())/RAND_MAX;}

int main(int argc, char* argv[]) {
	int box_dimension;
	float alpha, beta, gamma;
	
	srand (time(NULL));
	
	float basis[4][3] = { {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0, 0.5}, {0, 0.5, 0.5} };
	
	if (argc > 1)
		box_dimension = atoi(argv[1]);
	else box_dimension = SIZE;
	
	//rotation angles
	// alpha = 0;
	// beta = 0;
	// gamma = 0;
	alpha = rnd() * PI / 2;
	beta = rnd() * PI / 2;
	gamma = rnd() * PI / 2;
	
	ofstream myfile;
	myfile.open ("fcc_lattice.xyz");
	myfile << box_dimension * box_dimension * box_dimension * 4 << "\n";
	myfile << "FCC \n";
	
	for (int i = 0; i < box_dimension; i++)
		for (int j = 0; j < box_dimension; j++)
			for (int k = 0; k < box_dimension; k++){
				int basis_length = sizeof(basis)/sizeof(basis[0]);
				
				for (int m = 0; m < basis_length; m++){
				
					float x, y, z;
				
					//basis atom coords in norotated system
					x = LATTICE * i + basis[m][0] * LATTICE;
					y = LATTICE * j + basis[m][1] * LATTICE;
					z = LATTICE * k + basis[m][2] * LATTICE;

					//rotation OX
					y = cos(alpha) * y - sin(alpha) * z;
					z = sin(alpha) * y + cos(alpha) * z;
					
					//rotation OY
					x = cos(alpha) * x + sin(alpha) * z;
					z = -sin(alpha) * x + cos(alpha) * z;
					
					//rotation OZ
					x = cos(alpha) * x - sin(alpha) * y;
					y = sin(alpha) * x + cos(alpha) * y;
					
					myfile << "Al " << x << " " << y << " " << z << "\n";
				}
			}
	
	myfile.close();
	return 0;
}