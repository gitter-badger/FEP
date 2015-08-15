#include <exception>
#include <iostream>

#include "BandedSymmetricMatrix.h"

BandedSymmetricMatrix::BandedSymmetricMatrix(uint32_t r) :
	rows(r)
{
	this->elements.reserve(r);
	this->elements.clear();
	this->offsets.reserve(r);
	/*initialize the row of offsets*/
	for(std::size_t ii = 0; ii < r; ++ii){
		this->offsets[ii] = ii;
		
		/*initialize the main diagonal*/
		this->elements.push_back(std::vector< double >());
		this->elements[ii].clear();
		this->elements[ii].push_back(0.0);
	}

}

BandedSymmetricMatrix::~BandedSymmetricMatrix()
{

}

double BandedSymmetricMatrix::operator()(std::size_t ii, std::size_t jj) const
{
	/*internally the storage is lower triangular, so swap indicies if
	* requested component is upper triangular*/
	if( jj > ii ) {
		/*swap rows an columns*/
		jj ^= ii;
		ii ^= jj;
		jj ^= ii;
	}
	/*if the requested element is no tracked, it is zero*/
	if(jj < this->offsets[ii]) {
		return 0.0;
	}
	return this->elements[ii][(ii - jj)];
}

double& BandedSymmetricMatrix::operator()(std::size_t ii, std::size_t jj)
{
	/*internally the storage is lower triangular, so swap indicies if
	* requested component is upper triangular*/
	if( jj > ii ) {
		/*swap rows an columns*/
		jj ^= ii;
		ii ^= jj;
		jj ^= ii;
	}
	/*check if that row has the column initialized, otherwise
	* we must allocate new mememory up to that point*/
	while(jj < this->offsets[ii]) {
		this->elements[ii].push_back(0.0);
		this->offsets[ii] -= 1;
	}
	return this->elements[ii][(ii-jj)];
}

