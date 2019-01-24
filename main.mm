
#include <ctime>
#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "solver.hpp"


int main(int argc, const char * argv[]) {
    
    unsigned long start_s = clock(); // display time execution
    
    Pde_solver my_solver;
    my_solver.pricing(true); // display result

    std::cout << std::endl << "time execution: " << (clock()-start_s)/double(CLOCKS_PER_SEC) << "sec" << std::endl;

    return 0;
}
