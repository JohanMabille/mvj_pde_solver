
#ifndef matrix_hpp
#define matrix_hpp

#include <vector>

class Tridiagonal_matrix
{

public:
    
    Tridiagonal_matrix(size_t N);
    Tridiagonal_matrix(Tridiagonal_matrix const& temp); // copy constructor : have to redefine it since the attribute _tab is a pointer (by default, pointer _tab to be create and pointer _tab copy will be equal ==> big issue). But in this programm, copy constructor is never called.
    ~Tridiagonal_matrix();  // destructor (need to destroy pointers)
    std::size_t dim() const; // dimension of the Tridiagonal matrix
    void fill_matrix(std::vector <std::vector <double>> boundary, std::vector<double> interior); // Specific construction of a tridiagonal matrix for the pde solver (specific coefficient at the boundary) with 
    const double& operator()(std::size_t i, std::size_t j) const; // return A(i,j)
    bool is_dominant() const; // is the matrix diagonally dominant ? Helpfull for the pde solver to use thomas algorithm after (need this hypothesis for stability)
    std::vector<double> col(std::size_t i) const; // ith column of A
    double det(size_t N); // compute the determinant of A (maybe helpfull for the pde solver to the LU decomposition). N : dimension of the current matrix computed (since the function is called by reccurence by taking into account the minor of the initial matrix, the dimension N change)
    
private:
    
    const size_t _N; // dimension
    double **_tab; // data matrix
    
    friend const std::vector<double> operator* (const Tridiagonal_matrix&, const std::vector<double>&); // multiplication of a matrix with a vector
    
};

std::ostream& operator<<(std::ostream& out, const Tridiagonal_matrix& m); // display the matrix
const std::vector<double> operator* (const Tridiagonal_matrix &A, const std::vector<double> &b); // multiplication of a matrix with a vector


#endif
