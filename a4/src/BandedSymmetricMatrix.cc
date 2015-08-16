#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <assert.h>

#include "BandedSymmetricMatrix.h"

BandedSymmetricMatrix::BandedSymmetricMatrix()
{
	this->rows = 0;
}

BandedSymmetricMatrix::~BandedSymmetricMatrix()
{
}

double BandedSymmetricMatrix::operator()(std::size_t ii, std::size_t jj) const
{
	/*sanity checks*/
	if(ii >= this->rows) {
		throw std::out_of_range("row index out of range");
	}
	if( jj >= this->rows) {
		throw std::out_of_range("column index out of range");
	}
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
	if(jj < this->offsets[ii]) {
		this->elements[ii].resize((ii-jj+1), 0.0);
		this->offsets[ii] = jj;
	}
	return this->elements[ii][(ii-jj)];
}

void BandedSymmetricMatrix::setSize(std::size_t r)
{
	/*nothing to do if it is same size*/
	if(r == this->rows) {
		return;
	} else if(r < this->rows) {
		/*shrink the vector by truncating*/
		this->elements.resize(r);
		this->offsets.resize(r);
		this->rows = r;
	} else {
		/*expand the vector, adding in new once cell vectors*/
		this->elements.resize(r, std::vector <double>(1, 0.0));
		this->offsets.resize(r);
		for(std::size_t ii = this->rows; ii < r; ++ii) {
			this->offsets[ii] = ii;
		}
		this->rows = r;
	}
	assert(this->rows == this->elements.size());
	assert(this->rows == this->offsets.size());
}

void BandedSymmetricMatrix::zero()
{
	for(std::size_t ii = 0; ii < this->rows; ++ii){
		std::fill(this->elements[ii].begin(), this->elements[ii].end(), 0.0);
	}
}
