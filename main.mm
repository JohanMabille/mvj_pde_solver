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
    
    std::string res;
    while(res.compare("y")!=0 && res.compare("n")!=0)
    {
        std::cout << "Would you like to use predefined parameter for pricing a call ? (y/n)" << std::endl;
        std::cin >> res;
    }
    
    unsigned long start_s = clock(); // display time execution
    
    if(res=="y")
    {
        double theorical_price = dauphine::bs_price(100, 100, 0.2, 1, true);
        Pde_solver my_solver;
        my_solver.pricing();
        double price = my_solver.display_price();
        std::vector<double> greeks = my_solver.display_greeks();
        
        std::cout << "Contract price today : " << " " << price << std::endl;
        std::cout << "Error pricing : " << " " << abs(theorical_price-price) << std::endl << std::endl << std::endl;
        std::cout << "Delta : " << " " << greeks[0] << std::endl;
        std::cout << "Gamma : " << " " << greeks[1] << std::endl;
        std::cout << "Theta : " << " " << greeks[2] << std::endl;
        
        std::cout << std::endl << "time execution: " << (clock()-start_s)/double(CLOCKS_PER_SEC) << "sec" << std::endl;
    }
    
    else
    {
        double S0;
        std::cout << "Forward ?" << std::endl;
        std::cin >> S0;
        
        double T;
        std::cout << "Maturity ? (1 for 1 year)" << std::endl;
        std::cin >> T;
        
        double sigma;
        std::cout << "Volatility ? (0.2 for 20%)" << std::endl;
        std::cin >> sigma;
        
        double r;
        std::cout << "Interest rate ? (0.2 for 20%)" << std::endl;
        std::cin >> r;
        
        double theta=2;
        while(theta>1 || theta<0)
        {
            std::cout << "Theta ? (between 0 and 1)" << std::endl;
            std::cin >> theta;
        }
        
        size_t Nx;
        std::cout << "Space size ?" << std::endl;
        std::cin >> Nx;
        
        size_t Nt;
        std::cout << "Time size ?" << std::endl;
        std::cin >> Nt;
        
        std::string type;
        while(type.compare("dirichlet")!=0 && type.compare("neumann")!=0)
        {
            std::cout << "Type of boundary condition ? ('dirichlet' or 'neumann')" << std::endl;
            std::cin >> type;
        }
        
        std::vector<double> value_boundary(2, 0);
        std::cout << "Initial value of " << type << " condition" << std::endl;
        std::cin >> value_boundary[0];
        std::cout << "Terminal value of " << type << " condition" << std::endl << std::endl;
        std::cin >> value_boundary[1];
        std::cout << "If you want to change the payoff, please change it in the solver.cpp file" << std::endl;
        
        Pde_solver my_solver(S0, T, sigma, r, theta, Nx, Nt, call, type, value_boundary);
        my_solver.pricing();
        double price = my_solver.display_price();
        std::vector<double> greeks = my_solver.display_greeks();
        
        std::cout << "Contract price today : " << " " << price << std::endl << std::endl;
        std::cout << "Delta : " << " " << greeks[0] << std::endl;
        std::cout << "Gamma : " << " " << greeks[1] << std::endl;
        std::cout << "Theta : " << " " << greeks[2] << std::endl;
        
        std::cout << std::endl << "time execution: " << (clock()-start_s)/double(CLOCKS_PER_SEC) << "sec" << std::endl;
    }
    
    
    
    return 0;
}
