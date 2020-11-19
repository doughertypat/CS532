#include <iostream>
#include <sstream>
#include <omp.h>

int main(int argc, char *argv[]) {
    int x = 5;

#pragma omp parallel private(x)
    {
        x = x + 1;
        std::ostringstream oss;
        oss << "private x: " << x << std::endl;
        std::cout << oss.str();
    }
    std::cout << "after x: " << x << std::endl;
}

