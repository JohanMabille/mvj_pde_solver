#include <iostream>
#include "matrix.hpp"
#include <vector>
#include <numeric>
#include <string>
// Missing include
#include <algorithm>
 
// METHODS

Tridiagonal_matrix::Tridiagonal_matrix(size_t N): _N(N)
{
    // With vector of vector, this can be put in the initializer list (before
    // the body of the constructor)
    // _tab(_N, std::vector<double>(_N, 0.));
    // With strided index scheme, this is even simpler:
    // _tab(_N * _N, 0.)
    _tab = new double *[_N]; // dynamic allocation of the pointer
    for (std::size_t i=0; i<_N; i++) // fill tab (pointer of pointer = matrix)
    {
        _tab[i] = new double[_N]; // dynamic allocation of the pointer of pointer
        std::fill(&_tab[i][0], &_tab[i][_N], 0); // all the coefficient to zero
    }
}

Tridiagonal_matrix::Tridiagonal_matrix(Tridiagonal_matrix const& temp): _N(temp._N)
{
    // WIth vector of vector or 1D buffer, no need for that
    _tab = new double *[_N];
    for (std::size_t i=0; i<_N; i++)
    {
        _tab[i] = new double[_N];
        for (std::size_t j=0; j<_N; j++)
        {
            *(_tab[i]+j) = temp(i,j); // copy coef by coef
        }
    }
}

Tridiagonal_matrix::~Tridiagonal_matrix()
{
    for (std::size_t i=0; i<_N; i++) // first, delete the pointer of pointer (second layer)
    {
        delete _tab[i];
    }
    
    delete [] _tab; // then the first layer of pointer
}

std::size_t Tridiagonal_matrix::dim() const
{
    return _N;
}

void Tridiagonal_matrix::fill_matrix(std::vector <std::vector <double>> boundary, std::vector<double> interior)
{
    for (std::size_t i=0; i<_N; i++) // fill tab (pointer of pointer = matrix)
    {
        if(i==0 || i==_N-1) // if in the boundary of the matrix
        {
            if(i==0)
            {
                *(_tab[i]) = boundary[0][0];
                *(_tab[i]+1) = boundary[0][1];
            }
            if(i==_N-1)
            {
                *(_tab[i]+_N-2) = boundary[1][0];
                *(_tab[i]+_N-1) = boundary[1][1];
            }
        }
        else // if not in the boundary of the matrix
        {
            *(_tab[i]+i-1) = interior[0];
            *(_tab[i]+i) = interior[1];
            *(_tab[i]+i+1) = interior[2];
        }
    }
}

const double& Tridiagonal_matrix::operator()(std::size_t i, std::size_t j) const
{
    return *(_tab[i]+j); // A(i,j)
}

bool Tridiagonal_matrix::is_dominant() const
{
    if(*(_tab[0])<=*(_tab[0]+1)) // condition for the first line
    {
        return false;
    }
    
    for(std::size_t i=1; i<_N-1; i++) // condition for the interior line
    {
        if(*(_tab[i]+i)<=*(_tab[i]+i-1)+*(_tab[i]+i+1))
        {
            return false;
        }
    }
    
    if(*(_tab[_N-1])<=*(_tab[_N-1]-1)) // condition for the second line
    {
        return false;
    }
    
    return true;
}

std::vector<double> Tridiagonal_matrix::col(std::size_t j) const
{
    // Why storing array of row if you operate on columns
    // (and not on rows)?
    std::vector<double> res(_N);
    for(std::size_t i=0; i<_N; i++)
    {
        res[i] = this->operator()(i,j);
    }
    return res;
}

double Tridiagonal_matrix::det(size_t N)
{
    // First, compute the two cases N=1 and N=2
    if(N==1) // 1*1 dterminant is itself
    {
        return this->operator()(0,0);
    }
    
    if(N==2) // 2*2 dterminant can be easily computed
    {
        return this->operator()(0,0)*this->operator()(1,1) - this->operator()(0,1)*this->operator()(1,0);
    }
    
    // Then compute the determinant according to the type of the matrix (tridiagonal or not)
    return this->operator()(N-1,N-1)*Tridiagonal_matrix::det(N-1) - this->operator()(N-2,N-1)*this->operator()(N-1,N-2)*Tridiagonal_matrix::det(N-2);
    // complexity in O(N)
    // Remark : initial matrix is tridiagonal ==> its minor are tridiagonal

}



// FUNCTIONS



std::ostream& operator<<(std::ostream& out, const Tridiagonal_matrix& A)
{
    for(std::size_t i=0; i<A.dim(); i++)
    {
        for(std::size_t j=0; j<A.dim(); j++)
        {
            out << A(i, j) << ", ";
        }
        out << std::endl;
    }
    return out;
}

const std::vector<double> operator* (const Tridiagonal_matrix &A, const std::vector<double> &b)
{
    size_t N = A.dim();
    std::vector<double> res(N, 0);

    // Product of a tridiagonal matrix with a vector should be a O(3N) algorithm, not O(N²)
    for (std::size_t i=0; i<N; i++)
    {
        for (std::size_t j=0; j<N; j++)
        {
            std::vector<double> prod(N); // store the products a(ik)*b(kj), k=0...N-1
            std::transform(A._tab[i], A._tab[i]+N, b.begin(), prod.begin(), [](double arg1, double arg2) { return arg1 * arg2; }); // compute the products
            res[i] = std::accumulate(prod.begin(), prod.end(), 0., std::plus<double>()); // sum of the products
        }
    }
    return res;
}
