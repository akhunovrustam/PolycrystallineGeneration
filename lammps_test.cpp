#include <stdio.h>
#include <stdlib.h>
#include "lammps/src/library.h"
#include "lammps/src/lammps.h"
#include <iostream>

using namespace LAMMPS_NS;

int main()
{

	LAMMPS *lmp;
	lammps_open_no_mpi(1,NULL,(void **) &lmp);
	lammps_file(lmp, "lammps_script.lmp");
	double bla = lammps_get_thermo(lmp, "etotal");
	std::cout << "energy!!!!!!!!!!!! " << bla << "\n";
	return 0;
}
