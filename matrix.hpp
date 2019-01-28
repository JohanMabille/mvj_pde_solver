
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
    
    // You forgot to define assignment operator, thus the following will crash:
    // Trigiagonal matrix a(10), b(10);
    // ....
    // a = b;
    
    // You prevent the generation of move constructor and move assign operator
    // and don't implement them ...

private:
    
    const size_t _N; // dimension
    // Why double** instead of std::vector<std::vector<double>>?
    // Besides, vector of vector (or array of array) is not efficient,
    // Prefer strided scheme (1D buffer storing the rows one after the other
    // or the columns one after the other);
    // Storing n² numbers while the matrix is tridiagonal is pretty inefficient
    double **_tab; // data matrix
    
    friend const std::vector<double> operator* (const Tridiagonal_matrix&, const std::vector<double>&); // multiplication of a matrix with a vector
    
};

std::ostream& operator<<(std::ostream& out, const Tridiagonal_matrix& m); // display the matrix
const std::vector<double> operator* (const Tridiagonal_matrix &A, const std::vector<double> &b); // multiplication of a matrix with a vector


#endif
