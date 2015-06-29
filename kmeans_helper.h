#pragma once

#include "base.h"

#if DISTRIBUTED

#include "dmlc/data.h"

#else

#include "local/data.h"
#include <iostream>

#endif

#include <string>
#include <cstring>
#include <random>
#include <vector>
#include <cassert>

namespace dddml{
using namespace dmlc;
/*
* Class for possible types of centers
*
*/

//typedef std::size_t size_t //TODO

//inline real_t square(real_t x){ return x*x;}


template <typename I>
inline real_t squareDistBetweenRows(const Row<I> &r1, const Row<I> &r2)
{
	size_t i,j;
	real_t sqdist = 0.0;
	for (i = 0, j = 0; (i < r1.length && j < r2.length); )
	{
		if (r1.index[i] == r2.index[j])
		{
			sqdist += (r1.value[i] - r2.value[j])*(r1.value[i] - r2.value[j]);
			++i; ++j;
		}
		else if (r1.index[i] > r2.index[j])
		{
			sqdist += (r2.value[j])*(r2.value[j]);
			++j;
		}
		else
		{
			sqdist += (r1.value[i])*(r1.value[i]);
			++i;
		}
	}
	return sqdist;
}

template<typename I>
class center_t
{
	public:
	real_t *arr;
	size_t dim;
	const Row<I> *row_ptr;
	// type: denotes type of center
	// 1: array; 2: vector; 3:Row<I> 
	int type; 	
	
	center_t(int type, int dim) {constructor(type, dim);}
	
	~center_t() {destructor();}
	
	center_t(const Row<I> &row, int dim)
	{
		this->type = 3;
		this->row_ptr = &(row);
		this->dim = dim;
	}
	
	void print()
	{
		using namespace std;
		switch (this->type)
		{
			case 1:
				for (int i = 0; i < this->dim; ++i)
				{
					cout << arr[i] << " ";
				}
				cout << endl;
				break;
			case 3:
				for (int i = 0; i < row_ptr->length; ++i)
				{
					cout << row_ptr->index[i] << ':' << row_ptr->value[i] << ' ';
				}
				cout << endl;
		}
	}

	real_t squareDist(const Row<I> &r1)
	{
		real_t r = -1.0;
		switch(this->type)
		{
			case 1: 
			{
				r = squareDist12(r1, this->arr);
				break;
			}
			case 3:
				r = squareDistBetweenRows(r1, *row_ptr);
				break;
		}
		return r;	
	}
	
	real_t squareDist(const real_t *array)
	{
		//ensure sizes are same before calling
		switch(this->type)
		{
			case 1:
			{
				real_t sqdist  = 0;
				for (int i = 0; i < dim; ++i)
				{
					sqdist += (arr[i] - array[i]) * (arr[i] - array[i]);
				}
				return sqdist;
				break;
			}
			case 3:
				return squareDist12(*(this->row_ptr), array);
		}
		return -1;
	}
	real_t squareDist( center_t &center)
	{
		CHECK(this->dim == center.dim);
		switch(this->type)
		{
			case 1:
			{
				return center.squareDist(this->arr);
				break;
			}
			case 3:
				return center.squareDist(*(this->row_ptr));
		}
		return -1;
	}
	
	/* resets *this to zeros */
	inline void reset()
	{
		CHECK(type == 1);
		std::memset(arr, 0, sizeof(real_t)*dim);
	}
	void reset(int new_type)
	{
		CHECK(new_type == 1);
		if (this->type != new_type)
		{
			size_t dim1 = this->dim;
			destructor();
			constructor(dim1);
		}
		else
			this->reset();
	}
	
	/* add a row (r1) to *this */
	void add_into(const Row<I> &r1)
	{	
		CHECK(type == 1);
		for(size_t i = 0; i < r1.length; ++i)
		{ 
			arr[r1.index[i]] += r1.weight * r1.value[i];
			
		}
	}
	
	/* divide *this by a scalar div */
	void divide_by(const real_t div)
	{
		CHECK(type == 1);
		for(size_t i = 0; i < dim; ++i)
		{
			arr[i] = (div == 0) ? 0: arr[i] / div;
		}
	}
	
	private:
	inline void constructor(int dim)
	{
		//CHECK(type == 1);
		this->type = 1;
		this->dim = dim;
		arr = new real_t[dim];
		std::memset(arr, 0, sizeof(real_t)*dim);
				
	}
	inline void destructor()
	{
		
		if ((this->type == 1) && (this->arr != NULL)) delete[] arr;
		
	}
	
	
	inline real_t squareDist12(const Row<I> &r1, const real_t *r2)
	{
		#if 0
		real_t sqdist = 0;
		for (int i = 0; i < this->dim; ++i) sqdist += r2[i] * (r2[i]);
		for (int i = 0; i < r1.length; ++i) sqdist += (r1.value[i])*(r1.value[i]);
		sqdist -= 2 * r1.SDot(r2, this->dim);
		return sqdist;
		#endif
		//#if 0
		size_t i,j;
		real_t sqdist = 0.;
		for (i = 0, j = 0; (i < r1.length && j < this->dim); )
		{
			if (r1.index[i] == j)
			{
				sqdist += (r1.value[i] - r2[j])*(r1.value[i] - r2[j]);
				++i; ++j;
			}
			else if (r1.index[i] > j)
			{
				sqdist += (r2[j] * r2[j]);
				++j;
			}
			else
			{
				CHECK(-2 == 0);
				//sqdist += (r1.value[i] * r1.value[i]);
				//++i;
			}
		}
		return sqdist;
		//#endif
	}
	
};


}//namespace dddml
