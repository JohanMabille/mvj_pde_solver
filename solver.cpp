//
//  tridiagonal_solver_system.cpp
//  
//


#include <iostream>
#include "solver.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>


// METHODS


Pde_solver::Pde_solver(double S0, double T, double sigma, double r, double theta, size_t Nx, size_t Nt, double dx, double dt,  double (*payoff)(double), std::string boundary, std::vector<double> value_boundary): _S0(S0), _T(T), _sigma(sigma), _r(r), _theta(theta), _Nx(Nx), _Nt(Nt), _dx(5*_sigma*sqrt(_T)), _dt(_T/_Nt), _space_mesh(_Nx+1, 10*_sigma*sqrt(_T)/_Nx), _time_mesh(_Nt+1, _dt), _payoff(payoff), _boundary(boundary), _value_boundary(value_boundary), _A(_Nx-1), _Aprime(_Nx-1), _u(_Nx-1, 0)
{
    define_matrixes();

    _time_mesh[0] = 0;
    std::partial_sum(_time_mesh.begin(), _time_mesh.end(), _time_mesh.begin(), std::plus<double>()); // build time mesh : (0, dt, 2dt, ... _Nt*dt)
    
    _space_mesh[0] = 0;
    std::partial_sum(_space_mesh.begin(), _space_mesh.end(), _space_mesh.begin(), std::plus<double>());
    std::transform(_space_mesh.begin(), _space_mesh.end(), _space_mesh.begin(), bind2nd(std::plus<double>(), log(_S0)-5*_sigma*sqrt(_T))); // build space mesh : centered in log(_S0)
}

Pde_solver::Pde_solver(): _S0(100), _T(1), _sigma(0.2), _r(0.03), _theta(0.5), _payoff(call), _Nx(100), _Nt(_T*365), _dx(5*_sigma*sqrt(_T)), _dt((double)1/365), _space_mesh(_Nx+1, 10*_sigma*sqrt(_T)/_Nx), _time_mesh(_Nt+1, _dt), _boundary("dirichlet"), _value_boundary({0, exp(_space_mesh[_Nx])-1500}), _A(_Nx-1), _Aprime(_Nx-1), _u(_Nx-1, 0)
{
    define_matrixes();

    _time_mesh[0] = 0;
    std::partial_sum(_time_mesh.begin(), _time_mesh.end(), _time_mesh.begin(), std::plus<double>());// build time mesh : (0, dt, 2dt, ... _Nt*dt)
    
    _space_mesh[0] = 0;
    std::partial_sum(_space_mesh.begin(), _space_mesh.end(), _space_mesh.begin(), std::plus<double>());
    std::transform(_space_mesh.begin(), _space_mesh.end(), _space_mesh.begin(), bind2nd(std::plus<double>(), log(_S0)-5*_sigma*sqrt(_T))); // build space mesh : centered in log(_S0)
    
}

void Pde_solver::define_matrixes()
{
    // a(θ), b(θ), c(θ)
    auto a = [&](double theta) { return -_dt*_theta*(pow(_sigma,2)/(2*pow(_dx,2))) + (pow(_sigma,2)-_r)/(4*_dx); };
    auto b = [&](double theta) { return 1 + pow(_sigma,2)*_theta*_dt/pow(_dx,2) + _r*_theta*_dt; };
    auto c = [&](double theta) { return _dt* theta*(-pow(_sigma,2)/(2*pow(_dx,2)) + (pow(_sigma,2)-_r)/(4*_dx)); };
    
    if(_boundary=="dirichlet")
    {
        std::vector <std::vector <double>> matrix_boundary({{b(_theta), c(_theta)}, {a(_theta), b(_theta)}});
        std::vector<double> matrix_interior({a(_theta), b(_theta), c(_theta)});
        _A.fill_matrix(matrix_boundary, matrix_interior);
        
        std::vector <std::vector <double>> matrix_boundary_prime({{b(_theta-1), c(_theta-1)}, {a(_theta-1), b(_theta-1)}});
        std::vector<double> matrix_interior_prime({a(_theta-1), b(_theta-1), c(_theta-1)});
        _Aprime.fill_matrix(matrix_boundary_prime, matrix_interior_prime);
        
        _u[0] = _value_boundary[0]*(a(_theta-1)-a(_theta));
        _u[_Nx-2] = _value_boundary[1]*(c(_theta-1)-c(_theta));
    }
    else
    {
        std::vector <std::vector <double>> matrix_boundary({{a(_theta)+b(_theta), c(_theta)}, {a(_theta), b(_theta)+c(_theta)}});
        std::vector<double> matrix_interior({a(_theta), b(_theta), c(_theta)});
        _A.fill_matrix(matrix_boundary, matrix_interior);
        
        std::vector <std::vector <double>> matrix_boundary_prime({{a(_theta-1)+b(_theta-1), c(_theta-1)}, {a(_theta-1), b(_theta-1)+c(_theta-1)}});
        std::vector<double> matrix_interior_prime({a(_theta-1), b(_theta-1), c(_theta-1)});
        _Aprime.fill_matrix(matrix_boundary_prime, matrix_interior_prime);
        
        _u[0] = _value_boundary[0]*_dx*(a(_theta)-a(_theta-1));
        _u[_Nx-2] = _value_boundary[1]*_dx*(c(_theta-1)-c(_theta));
    }
}

std::vector<double> Pde_solver::vector_system(const std::vector<double> &f) const
{
    std::vector<double> res = _Aprime*f;
    std::transform(res.begin(), res.end(), _u.begin(), res.begin(), std::plus<double>());
    
    return res;
}

std::vector <double> Pde_solver::pricing(bool display) const
{
    std::cout << "Is A dominant ? (0,1) : " << bool(_A.is_dominant()) << std::endl; // look if A is dominant to apply thomas algorithm
    
    std::vector<double> f(_Nx-1, 0); // compute f(Nt) which os the terminal payoff
    std::transform(_space_mesh.begin(), _space_mesh.end(), f.begin(), [&](double arg1) { return (*_payoff)(exp(arg1)); });
    
    for(std::size_t i=0; i<_Nt; i++)
    {
        f = thomas_algorithm(_A, vector_system(f)); // compute f(n) knowing f(n+1)
    }
    
    if(display)
    {
        dispaly_price(f); // display f(0) : price of the contract at t=0
    }
    
    return f;
}

void Pde_solver::dispaly_price(const std::vector <double> &f) const
{
    for(std::size_t i=0; i<f.size(); i++)
    {
        std::cout << "Space index : " << i+1 << ", " << "Contract Price : " << f[i] << ", " << "Spot : " << exp(_space_mesh[i+1]) << std::endl;
    }
}


// FUNCTIONS




double call(double S)
{
    return std::max(S-70, 0.); // call strike K=70
}

std::vector<double> thomas_algorithm(const Tridiagonal_matrix &A, const std::vector<double> &y)
{
    size_t N = A.dim();
    std::vector<double> alpha(N, 0);
    std::vector<double> beta(N, 0);
    std::vector<double> x(N, 0);
    
    alpha[1] = -A(0,1)/A(0,0);
    beta[1] = y[0]/A(0,0);
    for (std::size_t i=1; i<N-1; i++)
    {
        alpha[i+1] = -A(i,i+1)/(A(i,i-1)*alpha[i]+A(i,i));
        beta[i+1] = (y[i]-beta[i]*A(i,i-1)) / (A(i,i-1)*alpha[i]+A(i,i));
    }
    
    x[N-1] = (y[N-1]-beta[N-1]*A(N-1,N-2)) / (A(N-1,N-2)*alpha[N-1]+A(N-1,N-1));
    for (std::size_t i=N-1; i>=1; i--)
    {
        x[i-1] = alpha[i]*x[i] + beta[i];
    }
    
    return x;
}
