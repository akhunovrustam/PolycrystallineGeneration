// Polycrystal rebuild from the size distribution using genetics algorithm
//
// Author   : Rustam Akhunov (GUT Gdansk)
// Email    : akhunovrustam@mail.ru
// Date     : 

#include "voro++.hh"
#include "ExtraStructAndFunc.hh"
#include "GeneticAlgoForSizesClass.cc"
#include <thread> 
#include <chrono>

using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds

int main(int argc, char *argv[]) {
	
	GeneticAlgoForSizesClass *algo = new GeneticAlgoForSizesClass(20); // created algo
	
	//initial population
	for (int i = 0; i < 50; i++){
		double point = (float)i * 0.1;
		algo->write_penalty_step("./test/original_dist.txt", point, algo->original_distribution(point));
		
	}
}
