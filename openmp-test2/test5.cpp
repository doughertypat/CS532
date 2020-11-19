#include <iostream>
#include <sstream>
#include <omp.h>

int main(int argc, char *argv[]) {
    const int N = 5;
    int tmp = 0;

#pragma omp parallel for lastprivate(tmp)
    for (int i = 0; i < N; i++) {
        tmp = i;
        std::ostringstream oss;
        oss << "inner tmp: " << tmp << std::endl;
        std::cout << oss.str();
    }

    std::cout << "tmp: " << tmp << std::endl;
}

