#include <thread>
#include <future>
#include <iostream>
#include <chrono>
#include <thread>
#include <time.h>
#include <math.h>

using namespace std;

int simplefunc(std::string a)
{
	double sum = 0;
    for (int i = 0; i < 100000000; i++)
		sum += sqrt(i);
    return sum;
}

int main()
{
    auto start = chrono::steady_clock::now();

	// double simple = simplefunc("hello world");
	// double simple1 = simplefunc("hello world");
	// double simple2 = simplefunc("hello world");
	
	auto future = std::async(simplefunc, "hello world");
	auto future1 = std::async(simplefunc, "hello worl!");
	auto future2 = std::async(simplefunc, "hello world!!");
	
	cout << "futures \n";
	double simple = future.get();
	double simple1 = future1.get();
	double simple2 = future2.get();

	auto end = chrono::steady_clock::now();
	auto diff = end - start;

	cout << chrono::duration <double, milli> (diff).count() << " ms" << endl;
	cout << simple << " " << simple1 << " " << simple2 << "\n";
	return simple;
}