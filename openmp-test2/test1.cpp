#include <iostream>
#include <sstream>
#include <omp.h>

int main(int argc, char *argv[]) {
    int x = 5;

    std::cout << omp_get_wtick() << std::endl;

#pragma omp parallel
    {
        x = x+1;
        std::ostringstream oss;
        oss << "shared: x is " << x << std::endl;
        std::cout << oss.str();
    }
}

