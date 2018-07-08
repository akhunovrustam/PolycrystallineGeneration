// Main programm calling all subclasses
//
// Author   : Rustam Akhunov (GUT Gdansk)
// Email    : akhunovrustam@mail.ru
// Date     : 

#include "voro++.hh"
#include "LatticeGeneratorClass.cc"
#include "ExtraStructAndFunc.hh"
#include "GeneticAlgoForSizesClass.cc"
#include "GeneticAlgoForOrientationsClass.cc"

using namespace std;
using namespace this_thread;
using namespace chrono;


int main(int argc, char *argv[]) {
	GeneralClass generator = new GeneralClass();
	
	generator->GrainSizeGenerate();
	
	generator->GrainOrientationGenerate();
	
	generator->FillGrains();
	
	generator->RelaxFirst();
	
	generator->RelaxSecond();
	
	generator->OutputForLAMMPS();
	
	delete generator;
}