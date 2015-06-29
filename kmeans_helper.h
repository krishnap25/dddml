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


template<typename IndexType>
class MyRow {
 public:
  /*! \brief label of the instance */
  real_t label;
  /*! \brief weight of the instance */
  real_t weight;
  /*! \brief length of the sparse vector */
  size_t length;
  /*!
   * \brief index of each instance    
   */
   IndexType *index;
  /*!
   * \brief array value of each instance, this can be NULL
   *  indicating every value is set to be 1
   */
  real_t *value;
  /*! \return i-th feature index */
  inline IndexType get_index(size_t i) const {
    return index[i];
  }
  /*!    
   * \return i-th feature value, this function is always
   *  safe even when value == NULL
   */
  inline real_t get_value(size_t i) const {
    return value == NULL ? 1.0f : value[i];
  }


  MyRow(const Row<IndexType> &row)
  {
    label = row.label;
    weight = row.weight;
    length = row.length;
    index = new IndexType[length];
    for (int i = 0; i < length; ++i)
    {
      index[i] = row.index[i];
    }
    value = new real_t[length];
    for (int i = 0; i < length; ++i)
    {
      value[i] = row.value[i];
    }
  }
  MyRow(real_t label1, real_t weight1, IndexType *index1, real_t *value1)
  {
  	label = label1;
  	weight = weight1;
  	index = index1;
  	value = value1;
  }

};

inline real_t square(real_t x){ return x*x;}

template<typename I>
class center_t
{
	public:
	real_t *arr;
	size_t dim;
	std::vector<real_t>* vec;
	MyRow<I> *row_ptr;
	// type: denotes type of center
	// 1: array; 2: vector; 3:Row<I> 
	int type; 
	//real_t sqnorm;
	//bool norm_updated;
	
	
	center_t(int type, int dim) {constructor(type, dim);}
	
	~center_t() {destructor();}
	
	center_t(const Row<I> &row, int dim)
	{
		this->type = 3;
		this->row_ptr = new MyRow<I>(row);
		this->dim = dim;
		//updateNorm();
	}
	
	void print()
	{
		using namespace std;
		//cout << this->type << " print called\n";
		switch (this->type)
		{
			case 1:
				for (int i = 0; i < this->dim; ++i)
				{
					cout << arr[i] << " ";
				}
				cout << endl;
				break;
			case 2:
				for (int i = 0; i < this->dim; ++i)
				{
					cout << (*vec)[i] << " ";
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
			case 2:
			{
				auto r2 = (this->type == 1) ? this->arr : &((*(this->vec))[0]); //&vector[0] gives pointer to underlying array
				r = squareDist12(r1, r2);
				break;
			}
			case 3:
				r = squareDist3(r1);
				break;
		}
		return r;	
	}

	real_t squareDist(const MyRow<I> &r1)
	{
		real_t r = -1.0;
		switch(this->type)
		{
			case 1: 
			case 2:
			{
				auto r2 = (this->type == 1) ? this->arr : &((*(this->vec))[0]); //&vector[0] gives pointer to underlying array
				r = squareDist12(r1, r2);
				break;
			}
			case 3:
				r = squareDist3(r1);
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
			case 2:
			{
				real_t sqdist  = 0;
				auto r2 = (this->type == 1) ? this->arr : &((*(this->vec))[0]); //&vector[0] gives pointer to underlying array
				for (int i = 0; i < dim; ++i)
				{
					sqdist += square(r2[i] - array[i]);
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
		assert(this->dim == center.dim);
		switch(this->type)
		{
			case 1:
			case 2:
			{
				auto r2 = (this->type == 1) ? this->arr : &((*(this->vec))[0]); //&vector[0] gives pointer to underlying array
				return center.squareDist(r2);
				break;
			}
			case 3:
				return center.squareDist(*(this->row_ptr));
		}
		return -1;
	}
	
	/* resets *this to zeros */
	void reset()
	{
		switch(this->type)
		{
			case 1: 
				std::memset(arr, 0, sizeof(real_t)*dim);
				break;
			case 2:
				std::fill(vec->begin(), vec->end(), 0.0f);
				break;
			case 3:
				std::memset(row_ptr->value, 0, sizeof(real_t) * this->dim);		
		}
		//this->sqnorm = 0;
	}
	void reset(int new_type)
	{
		if (this->type != new_type)
		{
			size_t dim1 = this->dim;
			destructor();
			constructor(new_type, dim1);
		}
		else
			this->reset();
		//this->sqnorm = 0;
	}
	
	/* add a row (r1) to *this; leaves norm inconsistent*/
	void add_into(const Row<I> &r1)
	{
		for(size_t i = 0; i < r1.length; ++i)
		{
			if (type == 1) 
				arr[r1.index[i]] += r1.weight * r1.value[i];
			else if (type == 2)
				vec[0][r1.index[i]] += r1.weight * r1.value[i];
			else
			{
				auto r2 = this->row_ptr;
				size_t i,j;
				for (i = 0, j = 0; (i <= r1.length && j <= r2->length); )
				{
					if (r1.index[i] == r2->index[j])
					{
						r2->value[j] += r1.weight * r1.value[i];
						++i; ++j;
					}
					else if (r1.index[i] > r2->index[j])
					{
						++j;
					}
					else
					{
						++i;
					}
				}
			}
		}
		//this-> norm_updated = false;
	}
	
	
	/* divide *this by a scalar div; leaves norm inconsistent */
	void divide_by(const real_t div)
	{
		if (type == 3)
		{	
			if (div == 0) std::memset(row_ptr->value, 0, sizeof(real_t) * row_ptr->length);
			else 
			{
				for (int i = 0; i < row_ptr-> length; ++i)
				{
					row_ptr->value[i] /= div;
				}
			}
		}
		else //types 1 and 2
		{
			for(size_t i = 0; i < dim; ++i)
			{
				if (type == 1) 
					arr[i] = (div == 0) ? 0: arr[i] / div;
				else if (type == 2)
					vec[0][i] = (div == 0) ? 0: vec[0][i] / div;
			}
		}
		//this->sqnorm = (div == 0) ? 0: this->sqnorm / div;
		//this->norm_updated = this->norm_updated && true;
	}
	
	/*
	void updateNorm()
	{
		this->sqnorm = 0;
		switch(this->type)
		{
			case 1: 
				for(int i = 0; i < dim; ++i)
				{
					this->sqnorm += square(this->arr[i]);
				}
				break;
			case 2:
				for(int i = 0; i < vec->size(); ++i)
				{
					this->sqnorm += square( this->vec[0][i]);
				}
				break;
			case 3:
				for(int i = 0; i < row_ptr->length; ++i)
				{
					this->sqnorm += square(this->row_ptr->value[i]);
				}
		}
		this->norm_updated = true;
	}
	*/
	private:
	inline void constructor(int type, int dim)
	{
		assert(type > 0 && type <= 3);
		this->type = type;
		this->dim = dim;
		switch(this->type)
		{
			case 1: 
				arr = new real_t[dim];
				std::memset(arr, 0, sizeof(real_t)*dim);
				break;
			case 2:
				vec = new std::vector<real_t>(dim, 0.0f);
				break;
			case 3:
				int label1 = 1, weight1 = 1;
				int length1 = dim;
				I *index1 = new I[dim];
				for (int i = 0; i < dim; ++i) index1[i] = i;
				real_t * value1 = new real_t[dim];
				std::memset(value1, 0, sizeof(real_t) * dim);
				
				row_ptr = new MyRow<I>(label1, weight1, index1, value1);	

				break;
		}
		//this->sqnorm = 0;
	}
	inline void destructor()
	{
		switch(this->type)
		{
			case 1: 
				if (arr != NULL) delete[] arr;
				break;
			case 2:
				if (vec != NULL) delete vec;
				break;
			case 3:
				if (row_ptr->index != NULL) delete[] row_ptr->index;
				if (row_ptr->value != NULL) delete[] row_ptr->value;
				if (row_ptr != NULL) delete row_ptr; 			
		}
	}
	
	
	
	
	inline real_t squareDist12(const Row<I> &r1, const real_t *r2)
	{
		
		size_t i,j;
		real_t sqdist = 0.;
		for (i = 0, j = 0; (i < r1.length && j < this->dim); )
		{
			if (r1.index[i] == j)
			{
				sqdist += square(r1.value[i] - r2[j]);
				++i; ++j;
			}
			else if (r1.index[i] > j)
			{
				sqdist += square(r2[j]);
				++j;
			}
			else
			{
				sqdist += square(r1.value[j]);
				++i;
			}
		}
		return sqdist;
	}
	inline real_t squareDist12(const MyRow<I> &r1, const real_t *r2)
	{
		
		size_t i,j;
		real_t sqdist = 0.;
		for (i = 0, j = 0; (i < r1.length && j < this->dim); )
		{
			if (r1.index[i] == j)
			{
				sqdist += square(r1.value[i] - r2[j]);
				++i; ++j;
			}
			else if (r1.index[i] > j)
			{
				sqdist += square(r2[j]);
				++j;
			}
			else
			{
				sqdist += square(r1.value[j]);
				++i;
			}
		}
		return sqdist;
	}
	
	inline real_t squareDist3(const Row<I> &r1)
	{
		auto r2 = this->row_ptr;
		size_t i,j;
		real_t sqdist = 0.0;
		for (i = 0, j = 0; (i < r1.length && j < r2->length); )
		{
			if (r1.index[i] == r2->index[j])
			{
				sqdist += square(r1.value[i] - r2->value[j]);
				++i; ++j;
			}
			else if (r1.index[i] > r2->index[j])
			{
				sqdist += square(r2->value[j]);
				++j;
			}
			else
			{
				sqdist += square(r1.value[j]);
				++i;
			}
		}
		return sqdist;
	}
	inline real_t squareDist3(const MyRow<I> &r1)
	{
		auto r2 = this->row_ptr;
		size_t i,j;
		real_t sqdist = 0.0;
		for (i = 0, j = 0; (i < r1.length && j < r2->length); )
		{
			if (r1.index[i] == r2->index[j])
			{
				sqdist += square(r1.value[i] - r2->value[j]);
				++i; ++j;
			}
			else if (r1.index[i] > r2->index[j])
			{
				sqdist += square(r2->value[j]);
				++j;
			}
			else
			{
				sqdist += square(r1.value[j]);
				++i;
			}
		}
		return sqdist;
	}
};



template <typename I>
real_t squareDistBetweenRows(const Row<I> &r1, const Row<I> &r3)
	{
		auto r2 = &r3;
		size_t i,j;
		real_t sqdist = 0.0;
		for (i = 0, j = 0; (i < r1.length && j < r2->length); )
		{
			if (r1.index[i] == r2->index[j])
			{
				sqdist += square(r1.value[i] - r2->value[j]);
				++i; ++j;
			}
			else if (r1.index[i] > r2->index[j])
			{
				sqdist += square(r2->value[j]);
				++j;
			}
			else
			{
				sqdist += square(r1.value[j]);
				++i;
			}
		}
		return sqdist;
	}


}//namespace dddml
