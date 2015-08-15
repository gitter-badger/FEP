#ifndef SPARSE_SYMMETRIC_MATRIX
#define SPARSE_SYMMETRIC_MATRIX 1

#include <vector>
#include <stdint.h>

class BandedSymmetricMatrix
{
/* store the symmetric matrix as a Lower Triangular matrix using
* a series of row vectors as the element storage. These row vectors
* have their first element at the main diagonal, and expand from
* right to left along the row, so that the internal representation
* row vector's last element is the first non zero element in the 
* actual matrix. The offset vector keeps track of the column locatio
* of the first non zero element tracked by the matrix.
* Every other element of the matrix not specified is zero.
* Diagram of matrix
* offset  matrix
*   0       [0] 
*   1           [0]
*   2               [0]
*   2               [1   0]
*   1           [3   2   1   0] <--- this is internal vector for row 4
*   5                           [0]
*
*/
public:
    BandedSymmetricMatrix(uint32_t rows, uint32_t columns);
    ~BandedSymmetricMatrix();

    /*read only operator*/
    double operator()(std::size_t ii, std::size_t jj) const;

    double& operator()(std::size_t ii, std::size_t jj);

    /* output format corresponds to matlab's sparse matrix import*
    * format for ease of checking */
    friend std::ostream& operator<< (std::ostream& s, BandedSymmetricMatrix const& M)
    {
        for(std::size_t ii = 0; ii < M.rows; ++ii) {
            /*compute the column index max*/
            std::size_t col_max = ii - M.offsets[ii];
            /*iterate over each row vector backwards to keep ouput with 
            * ordered as if we were reading each row left to right*/
            for(std::size_t jj = 0; jj <= col_max; ++jj) {
                s << ii << '\t' << (M.offsets[ii] + jj) << '\t';
                s << M.elements[ii][col_max - jj] << '\n';
                /*ouput the symmetric element*/
                if( ii != (M.offsets[ii] + jj)) {
                    s << (M.offsets[ii] + jj) << '\t' << ii << '\t';
                    s << M.elements[ii][col_max - jj] << '\n';
                }
            }
        }
        return s;
    }

private:
    uint32_t rows;
    uint32_t columns;
    std::vector<std::size_t> offsets;
    std::vector <std::vector<double> > elements;
};

#endif
