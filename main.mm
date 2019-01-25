
#include <ctime>
#include <iostream>
#include <vector>
#include "matrix.hpp"
#include "solver.hpp"
#include "closed_form.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>


int main(int argc, const char * argv[]) {
    
    unsigned long start_s = clock(); // display time execution
    
    double theorical_price = dauphine::bs_price(100, 100, 0.2, 1, true);
    Pde_solver my_solver;
    double price = my_solver.display_price();
    std::vector<double> greeks = my_solver.display_greeks();
    
    std::cout << "Contract price today : " << " " << price << std::endl;
    std::cout << "Error pricing : " << " " << abs(theorical_price-price) << std::endl << std::endl << std::endl;
    std::cout << "Delta : " << " " << greeks[0] << std::endl;
    std::cout << "Gamma : " << " " << greeks[1] << std::endl;
    std::cout << "Theta : " << " " << greeks[2] << std::endl;

    std::cout << std::endl << "time execution: " << (clock()-start_s)/double(CLOCKS_PER_SEC) << "sec" << std::endl;

    return 0;
}
