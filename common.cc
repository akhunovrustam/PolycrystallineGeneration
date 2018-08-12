// Main programm calling all subclasses
//
// Author   : Rustam Akhunov (GUT Gdansk)
// Email    : akhunovrustam@mail.ru
// Date     : 


#include "GeneralClass.cc"

using namespace std;
using namespace this_thread;
using namespace chrono;


int main(int argc, char *argv[]) {
	General1Class* gn = new General1Class();
	
	cout << "size generation grains" << endl;
	gn->GrainSizeGenerate(argc, argv);
	exit(0);
	
	cout << "orient generation grains" << endl;
	gn->GrainOrientationGenerate(argc, argv);
	
	gn->writeAll();
	
	// gn->loadAll();
	cout << "fill grains" << endl;
	gn->FillGrains();
	
	cout << "relax 1" << endl;
	string str = gn->RelaxFirst(argc, argv);
	
	// cout << "relax 2" << endl;
	// gn->RelaxSecond(argc, argv, str);
	
	// gn->OutputForLAMMPS();
	
	delete gn;
}