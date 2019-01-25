
#ifndef tridiagonal_solver_system_hpp
#define tridiagonal_solver_system_hpp

#include <string>
#include "matrix.hpp"
#include "volatility.hpp"

class Pde_solver
{
    
public:
    
    Pde_solver();
    Pde_solver(double S0, double T, double sigma, double r, double theta, size_t Nx, size_t Nt, double (*payoff)(double), std::string boundary, std::vector<double> value_boundary);
    void define_matrixes(); // method to compute _A, _Aprime, _u (will be called inside the constructor)
    std::vector<double> vector_system(const std::vector<double> &f) const; // the vector of the right member of : A(θ)f(n) = A(θ-1)f(n+1)+u. Need to first solve f(n+1) to then obtain the vector.
    void pricing();
    void compute_greeks(const std::vector <std::vector <double>> &f);
    double display_price() const;
    std::vector <double> display_greeks() const;
    
private:
    
    const double _S0; // spot at t=0
    const double _T;
    Volatility _sigma;
    const double _r;
    const double _theta;
    double _price; // price of the contract today
    std::vector <double> _greeks;
    double (*_payoff)(double); // pointer to function (need to define it outside the class)
    
    const size_t _Nx; // number of space points
    const size_t _Nt; // number of time points (T*365)
    const double _dx; // space step
    const double _dt; // time step (1/365 by default)
    std::vector<double> _space_mesh;
    std::vector<double> _time_mesh;
    
    std::string _boundary; // type of the boundary : dirichlet or neumann
    std::vector<double> _value_boundary; // values of f0 and fN if _boundary == 'dirichlet'; values of f'0 and f'N if _boundary == 'neumann'
    Tridiagonal_matrix _A; // matrix A(θ) of the linear system to solve : A(θ)f(n) = A(θ-1)f(n+1)+u
    Tridiagonal_matrix _Aprime; // matrix A(θ-1) of the linear system to solve : A(θ)f(n) = A(θ-1)f(n+1)+u
    std::vector<double> _u; // vector b of the linear system to solve : A(θ)f(n) = A(θ-1)f(n+1)+u
    
};

double call(double S); // payoff by default
std::vector<double> thomas_algorithm(const Tridiagonal_matrix &A, const std::vector<double> &y); // find x such that Ax=y when A is tridiagonal (and diagonally dominante to be preferred for stability)


#endif /* solver_hpp */
