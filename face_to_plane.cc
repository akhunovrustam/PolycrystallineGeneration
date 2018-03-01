// Polycrystal lattice generation
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



int main() {
	LatticeGeneratorClass gen(3);
	gen.fillRandomlyAndBuildGrains();
}