
#include <iostream>
#include "solver.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

// METHODS


Pde_solver::Pde_solver(double S0, double T, double sigma, double r, double theta, size_t Nx, size_t Nt, double (*payoff)(double), std::string boundary, std::vector<double> value_boundary): _S0(S0), _T(T), _sigma(sigma), _r(r), _theta(theta), _price(0), _greeks(3, 0.), _Nx(Nx), _Nt(Nt), _dx(10*sigma*sqrt(_T)/_Nx), _dt(_T/_Nt), _space_mesh(_Nx+1, 10*_sigma.vol()*sqrt(_T)/_Nx), _time_mesh(_Nt+1, _dt), _payoff(payoff), _boundary(boundary), _value_boundary(value_boundary), _A(_Nx-1), _Aprime(_Nx-1), _u(_Nx-1, 0)
{
    define_matrixes();
    
    _time_mesh[0] = 0;
    std::partial_sum(_time_mesh.begin(), _time_mesh.end(), _time_mesh.begin(), std::plus<double>()); // build time mesh : (0, dt, 2dt, ... _Nt*dt)
    
    _space_mesh[0] = 0;
    std::partial_sum(_space_mesh.begin(), _space_mesh.end(), _space_mesh.begin(), std::plus<double>());
    std::transform(_space_mesh.begin(), _space_mesh.end(), _space_mesh.begin(), bind2nd(std::plus<double>(), log(_S0)-5*_sigma.vol()*sqrt(_T))); // build space mesh : centered in log(_S0) : first we created the discrete interval [0, 2c] of N+1 point by using the partial_sum function (= cummulative sum). So we have the set {0, 2c/N, 4c/N, ..., 2c} with c=5*_sigma*sqrt(_T). Then we translate the set by adding log(spot)-c to each point. So we have the set {log(spot)-c, ..., log(spot), ..., log(spot)+c}
}

Pde_solver::Pde_solver(): _S0(100), _T(1), _sigma(0.2), _r(0.0), _theta(0.5), _price(0.), _greeks(3, 0.), _payoff(call), _Nx(100), _Nt(_T*365), _dx(10*_sigma.vol()*sqrt(_T)/_Nx), _dt((double)1/365), _space_mesh(_Nx+1, 10*_sigma.vol()*sqrt(_T)/_Nx), _time_mesh(_Nt+1, _dt), _boundary("dirichlet"), _value_boundary({0, _S0*exp(5*_sigma.vol()*sqrt(_T))-100}), _A(_Nx-1), _Aprime(_Nx-1), _u(_Nx-1, 0)
{
    define_matrixes();
    
    _time_mesh[0] = 0;
    std::partial_sum(_time_mesh.begin(), _time_mesh.end(), _time_mesh.begin(), std::plus<double>());
    
    _space_mesh[0] = 0;
    std::partial_sum(_space_mesh.begin(), _space_mesh.end(), _space_mesh.begin(), std::plus<double>());
    std::transform(_space_mesh.begin(), _space_mesh.end(), _space_mesh.begin(), bind2nd(std::plus<double>(), log(_S0)-5*_sigma.vol()*sqrt(_T)));
    
}

void Pde_solver::define_matrixes()
{
    // a(θ), b(θ), c(θ)
    auto a = [&](double theta) { return -_dt*theta*(pow(_sigma.vol(),2)/(2*pow(_dx,2)) + (pow(_sigma.vol(),2)-_r)/(4*_dx)); };
    auto b = [&](double theta) { return 1 + (pow(_sigma.vol(),2)*theta*_dt/pow(_dx,2)) + _r*theta*_dt; };
    auto c = [&](double theta) { return _dt*theta*(-pow(_sigma.vol(),2)/(2*pow(_dx,2)) + (pow(_sigma.vol(),2)-_r)/(4*_dx)); };
    
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

void Pde_solver::pricing()
{
    std::cout << std::endl << "Is A dominant ? (0,1) : " << bool(_A.is_dominant()) << std::endl << std::endl; // look if A is dominant to apply thomas algorithm
    
    std::vector<double> f(_Nx-1, 0); // compute f(Nt) which os the terminal payoff
    std::transform(_space_mesh.begin(), _space_mesh.end(), f.begin(), [&](double arg1) { return (*_payoff)(exp(arg1)); });
    
    for(std::size_t i=0; i<_Nt-1; i++)
    {
        f = thomas_algorithm(_A, vector_system(f)); // compute f(n) knowing f(n+1)
    }
    
    std::vector <std::vector <double>> res(2); // store 2 price vector in time to compute the theta later
    res[1] = {f[(_Nx/2)-2], f[(_Nx/2)-1], f[(_Nx/2)]}; // prices around the forward at day 1
    f = thomas_algorithm(_A, vector_system(f));
    res[0] = {f[(_Nx/2)-2], f[(_Nx/2)-1], f[(_Nx/2)]}; // prices around the forward at day 0 (today)
    
    _price = f[(_Nx/2)-1];
    compute_greeks(res);
}

void Pde_solver::compute_greeks(const std::vector <std::vector <double>> &f)
{
    double h_1 = exp(_space_mesh[(_Nx/2)+2])-exp(_space_mesh[(_Nx/2)+1]);
    double h_2 = exp(_space_mesh[(_Nx/2)+1])-exp(_space_mesh[_Nx/2]);
    
    _greeks[0] = (f[0][2]-f[0][0])/(h_1+h_2); // delta
    _greeks[1] = 2*(f[0][2]+f[0][0]-2*f[0][1]-_greeks[0]*(h_1-h_2))/(pow(h_1,2)+pow(h_2,2)); // gamma
    _greeks[2] = f[1][1]-f[0][1]; // theta
}

double Pde_solver::display_price() const
{
    return _price;
}

std::vector <double> Pde_solver::display_greeks() const
{
    return _greeks;
}


// FUNCTIONS




double call(double S)
{
    return std::max(S-100, 0.);
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
