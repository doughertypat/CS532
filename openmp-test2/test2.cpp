#include <iostream>
#include <sstream>
#include <omp.h>

int main(int argc, char *argv[]) {
    int x = 5;
        
#pragma omp parallel
    {
        int x;
        x = 3;
        std::ostringstream oss;
        oss << "local: x is " << x << std::endl;
        std::cout << oss.str();
    }
    std::cout << "the global x: " << x << std::endl;
}

