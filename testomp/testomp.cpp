#include <iostream>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {
	using namespace std;
	ostringstream oss;
	oss << "Test OpenMP\n";
#ifdef _OPENMP
	int threadId;
#pragma omp parallel private(threadId) 
    {
	threadId = omp_get_thread_num();
	if (0 == threadId) {
		int maxThreads = omp_get_max_threads();
		int numThreads = omp_get_num_threads();
		oss << "\tOpenMP thread counts\n";
		oss << "\tMax: " << maxThreads << "\n";
		oss << "\tActual: " << numThreads << "\n";
	}
    }
#endif
	oss << "End of test.\n";
	cout << oss.str() << endl;
}

