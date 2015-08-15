#ifndef SPARSE_SYMMETRIC_MATRIX
#define SPARSE_SYMMETRIC_MATRIX 1

#include <vector>
#include <stdint.h>

class BandedSymmetricMatrix
{
/* store

*/
public:
	BandedSymmetricMatrix(uint32_t rows, uint32_t columns);
	~BandedSymmetricMatrix();

	/*read only operator*/
	double operator()(std::size_t ii, std::size_t jj) const;

	double& operator()(std::size_t ii, std::size_t jj);

	friend std::ostream& operator<< (std::ostream& s, BandedSymmetricMatrix const& M)
	{
		for(std::size_t ii = 0; ii < M.rows; ++ii) {
			/*compute the column index max*/
			std::size_t col_max = ii - M.offsets[ii];
			/*iterate over each row vector backwards to keep the 
			* output order lexographically ordered */
			for(std::size_t jj = 0; jj <= col_max; ++jj) {
				s << ii << '\t' << (M.offsets[ii] + jj) << '\t';
				s << M.elements[ii][col_max - jj] << '\n';
			}
		}
		return s;
	}

private:
	uint32_t rows;
	uint32_t columns;
	std::vector<std::size_t> offsets;
	std::vector <std::vector<double> > elements;


	/* data */
};


#endif
