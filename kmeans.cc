#include <string>
#include <random>
#include <vector>
#include "kmeans_helper.h"
#include "sample_helper.h"
#include <memory>
#include <queue>
#include <array>
#include <algorithm>

#if DISTRIBUTED 

#include <dmlc/data.h>
#include <dmlc/logging.h>
#include <dmlc/io.h>
#include "data/row_block.h"
#include <dmlc/timer.h>

#else

#include <iostream>
#include "local/data.h"
#include "local/io.h"
#include "local/row_block.h"
#include "local/timer.h"
#include "local/libsvm_parser.h"

#endif



/*
* -----------------------------------------------------------
*	IMPLEMENTATION of k-MEANS with k-means++ INITIALIZATION
* -----------------------------------------------------------
*/

namespace dddml{
using namespace dmlc;
using namespace dmlc::data;

template<typename I>
using center_t_ptr = std::shared_ptr<center_t<I>>;

template<typename I>
using vector_center_ptr = std::shared_ptr<std::vector<center_t_ptr<I>>>;

typedef std::shared_ptr<std::vector<int>> vector_int_ptr;

typedef std::shared_ptr<std::vector<std::vector<int>>> vector_vector_int_ptr;

/*
* functions to find nearest point from a given collection
*/

template<typename I>
int find_closest(const Row<I> &row, const std::vector<center_t_ptr<I>> &block)
{
	if (block.size() == 0) return -1;
	else if (block.size() == 1) return 0;
	int min_index = 0;
	real_t min_dist = block[0]->squareDist(row);
	real_t cur_dist;
	for (size_t i = 1; i < block.size(); ++i)
	{
		cur_dist = block[i]->squareDist(row);
		if (cur_dist < min_dist)
		{
			min_dist = cur_dist;
			min_index = i;
		} 
	}
	return min_index;
}


/*
* Find closest point from vector ignoring index i
*/
template<typename I>
int find_closest(const center_t_ptr<I> center, const vector_center_ptr<I> all_centers, int ignore_index)
{
	if (all_centers->size() == 0) return -1;
	else if (all_centers->size() == 1) return 0;
	int min_index = (ignore_index != 0) ? 0 : 1;
	real_t min_dist = (*all_centers)[min_index]->squareDist(*center);
	real_t cur_dist;
	for (size_t i = 1 + min_index; i < all_centers->size(); ++i)
	{
		if (i == ignore_index)
		{
			continue;
		}
		else
		{
			cur_dist = (*all_centers)[i]->squareDist(*center);
			if (cur_dist < min_dist)
			{
				min_dist = cur_dist;
				min_index = i;
			} 
		}	
	}
	return min_index;
}

template<typename I>
int *find_p_closest(int p, const Row<I> &row, const std::vector<center_t_ptr<I>> &centers)
{
	if (centers.size() < p)
	{
	//should not occur:
		exit(0);
	}
	else
	{
		real_t sqdist;
		//heap select to efficiently find p things
		std::priority_queue<std::pair<real_t, int>> pq;
		for (int i = 0; i < p; ++i) //insert for p elements into the heap
		{
			pq.push(std::pair<real_t, int>(centers[i]->squareDist(row), i));
		}
		//insert rest of the distances while maintaing p smallest in the heap
		for (int i = p; i < centers.size(); ++i)
		{
			sqdist = centers[i]->squareDist(row);
			if (pq.top().first > sqdist)
			{
				pq.pop();
				pq.push(std::pair<real_t, int>(sqdist, i));
			}
		}
		//read off assignments in reverse order
		int *assignments = new int[p];
		std::memset(assignments, 0, sizeof(int) * p);
		for (int i = p-1; i >= 0; --i)
		{
			assignments[i] = pq.top().second;
			pq.pop();
		}
		return assignments;
	}
}

/*
* Update assignment step of k-means based on Voronoi partitions
* /param centers: vector of centers
* /param data: datapoints that are to be reassigned
* /param assignments: vector of assignments from 1 to k. It will be modified
* /return true if any assignments have changed. False otherwise
*/
template<typename I>
bool update_assignments(const std::vector<center_t_ptr<I>> &centers, const RowBlock<I> &data, std::vector<int> &assignments)
{
	size_t k = centers.size(), n = data.size; 
	assert(data.size == assignments.size());
	bool changed = false;
	int assignment;
	for (size_t i = 0; i < n; ++i)
	{
		assignment = find_closest(data[i], centers);
		if (assignment != assignments[i])
		{
			assignments[i] = assignment;
			changed = true;
		}
	}
	return changed;
}

template<typename I>
bool update_assignments(const std::vector<center_t_ptr<I>> &centers, const RowBlock<I> &data, std::vector<std::vector<int>> &assignments)
{
	size_t k = centers.size(), p = assignments[0].size(), n = data.size; 
	assert(data.size == assignments.size());
	bool changed = false;
	std::vector<int> assignment;
	for (size_t i = 0; i < n; ++i)
	{
		assignment = find_p_closest(p, data[i], centers);
		if (assignment != assignments[i])
		{
			assignments[i] = assignment;
			changed = true;
		}
	}
	return changed;
}

/*
* Update centers step of k-means 
* /param centers: vector of centers (will be changed)
* /param data: datapoints that are to be reassigned
* /param assignments: vector of assignments from 1 to k. It will be modified
*/

template<typename I>
void update_centers(std::vector<center_t_ptr<I>> &centers, const RowBlock<I> &data, const std::vector<int> &assignments, int center_type)
{
	size_t k = centers.size(), n = data.size; 
	assert(data.size == assignments.size());
	for (auto c : centers)
	{
		c->reset(center_type);
	}
	int counts[k];
	for (int i = 0; i < k; ++i) counts[i] = 0;
	//std::memset(counts, 0, sizeof(int) * k);
	// aggregate all points in a cluster
	for (int i = 0; i < assignments.size(); ++i)
	{
		if (assignments[i] >= 0)
		{ //valid assignment
			centers[assignments[i]]->add_into(data[i]);
			counts[assignments[i]] += (data.weight == NULL) ? 1 : data.weight[i];
		}
	}
	//average
	for (int i = 0; i < k; ++i)
	{
		centers[i]->divide_by(counts[i]);
	}
}


template<typename I>
void update_centers(std::vector<center_t_ptr<I>> &centers, const RowBlock<I> &data, const std::vector<std::vector<int>> &assignments, int center_type)
{
	size_t k = centers.size(), p = assignments[0].size(), n = data.size; 
	assert(data.size == assignments.size());
	for (auto c : centers)
	{
		c->reset(center_type);
	}
	int counts[k];
	for (int i = 0; i < k; ++i) counts[i] = 0;
	//std::memset(counts, 0, sizeof(int) * k);
	// aggregate all points in a cluster
	for (int i = 0; i < assignments.size(); ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			if (assignments[i][j] >= 0)
			{ //valid assignment
				centers[assignments[i][j]]->add_into(data[i]);
				counts[assignments[i][j]] += (data.weight == NULL) ? 1 : data.weight[i];
			}
		}
	}
	//average
	for (int i = 0; i < k; ++i)
	{
		centers[i]->divide_by(counts[i]);
		//centers[i]->updateNorm();
	}
}

/*
* Update a single center: for balancing heuristics 
* /param centers: vector of centers (will be changed)
* /param data: datapoints that are to be reassigned
* /param assignments: vector of assignments from 1 to k. It will be modified
* /param center_id: ID of center to recompute
*/

template<typename I>
void update_one_center(std::vector<center_t_ptr<I>> &centers, const RowBlock<I> &data, const std::vector<int> &assignments, int center_id)
{
	size_t k = centers.size(), n = data.size; 
	assert(data.size == assignments.size());
	
	centers[center_id]->reset();
	
	int count = 0;
	// aggregate all points in a cluster
	for (int i = 0; i < assignments.size(); ++i)
	{
		if (assignments[i] == center_id)
		{ 
			centers[center_id]->add_into(data[i]);
			count += (data.weight == NULL) ? 1 : data.weight[i];
		}
	}
	//average
	centers[center_id]->divide_by(count);
}


template<typename I>
void update_one_center(std::vector<center_t_ptr<I>> &centers, const RowBlock<I> &data, const std::vector<std::vector<int>> &assignments, int center_id)
{
	size_t k = centers.size(), p = assignments[0].size(), n = data.size; 
	assert(data.size == assignments.size());
	centers[center_id]->reset();
	int count = 0;
	// aggregate all points in a cluster
	for (int i = 0; i < assignments.size(); ++i)
	{
		for (int j = 0; j < p; ++j)
		{
			if (assignments[i][j] == center_id)
			{ //valid assignment
				centers[center_id]->add_into(data[i]);
				count += (data.weight == NULL) ? 1 : data.weight[i];
			}
		}
	}
	//average
	centers[center_id]->divide_by(count);
}


/*
* Pick up random centers to initialize k-means (code: 0)
*/
template<typename I>
vector_center_ptr<I> random_init(const RowBlock<I> &data, int k, int dim, std::mt19937_64 &rng)
{
	size_t numData = data.size;
	int *sample = SampleWithoutReplacement(k, numData, rng);
	//alternative 1:
	vector_center_ptr<I> centers = std::make_shared<std::vector<center_t_ptr<I>>> ();
	centers->reserve(k);
	for (int i = 0; i < k; ++i)
	{
		(centers)->push_back(std::make_shared<center_t<I>>(data[sample[i]], dim));
	}
	//alternative 2:
	//vector_center_ptr<I> centers = std::make_shared<std::vector<center_t_ptr<I>>> (new std::vector<center_t_ptr<I>>[k]);
	// for (int i = 0; i < k; ++i)
	// {
	// 	(*centers)[i] = std::make_shared<center_t<I>>(data[sample[i]], dim);
	// }

	return centers;
}

/* 
*	k-means++ initialization (code: 1)
*/
template<typename I>
vector_center_ptr<I> kmpp_init(const RowBlock<I> &data, int k, int dim, std::mt19937_64 &rng)
{
	real_t weight; 
	//alterantive 2:
	//vector_center_ptr<I> centers = std::make_shared<std::vector<center_t_ptr<I>>> (new std::vector<center_t_ptr<I>>[k]);
	//vector_center_ptr<I> centers = std::make_shared<std::vector<center_t_ptr<I>>> (std::array<center_t_ptr<I>, k>());

	//alternative 1:
	vector_center_ptr<I> centers = std::make_shared<std::vector<center_t_ptr<I>>> ();
	centers->reserve(k);

	size_t numData = data.size;
	real_t *sqdists = new real_t[numData]; //distances
	for (int i = 0; i < numData; ++i) sqdists[i] = 0;
	// initialize first center
	std::uniform_int_distribution<> dis(0, numData-1);
		//(*centers)[0] = std::make_shared<center_t<I>>(data[dis(rng)], dim);
	(*centers).push_back ( std::make_shared<center_t<I>>(data[dis(rng)], dim));

	//initialize distances
	for (int j = 0; j < numData; ++j)
	{
		weight = (data.weight == NULL) ? 1 : data.weight[j];
		sqdists[j] = (weight) * (*centers)[0]->squareDist(data[j]);
	}
	//loop for the next (k-1) centers
	for (int i = 1; i < k; ++i)
	{
		//sample next center
		int next_center = weightedSample(sqdists, numData, rng);
			//(*centers)[i] = std::make_shared<center_t<I>>(data[next_center], dim);
		(*centers).push_back ( std::make_shared<center_t<I>>(data[next_center], dim));
		//update distances
		for (int j = 0; j < numData; ++j)
		{
			weight = (data.weight == NULL) ? 1 : data.weight[j];
			sqdists[j] = std::min(sqdists[j], (*centers)[i]->squareDist(data[j]));
		}
	}
	delete[] sqdists;
	return centers;	
}

template<typename I>
real_t kmeans_objective(const RowBlock<I> &data, vector_int_ptr assignments, vector_center_ptr<I> centers)
{
	size_t n = data.size;
	real_t obj = 0;
	int count = 0;
	for (int i = 0; i < n; ++i)
	{
		if ((*assignments)[i] != -1)
		{
			++ count;
			real_t weight = (data.weight == NULL) ? 1 : data.weight[i];
			obj += weight * (*centers)[(*assignments)[i]]->squareDist(data[i]);
		}
	}
	if (count != n) std::cerr << "ERROR@!\n";
	return obj/count;
}

template<typename I>
real_t kmeans_objective(const RowBlock<I> &data, vector_vector_int_ptr assignments, vector_center_ptr<I> centers)
{
	size_t n = data.size;
	int p = (*assignments)[0].size();
	int count = 0;
	real_t obj = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < p; ++j)
		{	
			if ((*assignments)[i][j] != -1)
			{
				++count;
				real_t weight = (data.weight == NULL) ? 1 : data.weight[i];
				obj += weight * (*centers)[(*assignments)[i][j]]->squareDist(data[i]);
			}
		}
	}
	return obj / count;
}

/*
* -----------------------------------------------------------------------
*	k-means
* ----------------------------------------------------------------------
*/
template<typename I>
std::pair<vector_int_ptr, vector_center_ptr<I>> kmeans(const RowBlock<I> &data, int k, int dim, int center_type, std::mt19937_64 &rng, int init = 1)
{
	double start, time;
	
	#if DISTRIBUTED
	LOG(INFO) << "Starting k-means";
	#else
	std::cout << "Starting k-means" << std::endl;
	#endif
	#if DISTRIBUTED
	LOG(INFO) << "Initializing centers.."; 
	#else
	std::cout << "Initializing centers.." << std::endl;
	#endif

	vector_center_ptr<I> centers;
	start = GetTime();
	if (init == 0) 
		centers = random_init(data, k, dim, rng);
	else //kmpp
		centers = kmpp_init(data, k, dim, rng);
	time = GetTime() - start;

	#if DISTRIBUTED
	LOG(INFO) << "Initialization done in " << time << " sec";
	#else
	std::cout << "Initialization done in " << time << " sec" << std::endl;
	#endif


	#if DISTRIBUTED
	LOG(INFO) << "Starting Lloyd's iterations";
	#else
	std::cout << "Starting Lloyd's iterations" << std::endl;
	#endif

	bool changed = true;
	vector_int_ptr assignments = std::make_shared<std::vector<int>>(data.size);
	start = GetTime();
	int iter_count = 0;
	while(changed)
	{
		changed = update_assignments(*centers, data, *assignments);
		update_centers(*centers, data, *assignments, center_type);
		++iter_count;
		std::cout << iter_count << ": objective: " << kmeans_objective(data, assignments, centers) << std::endl;
	}
	time = GetTime() - start;

	#if DISTRIBUTED
	LOG(INFO) << "Finished k-means: " << iter_count << " iterations in " << time << " sec";
	#else
	std::cout << "Finished k-means: " << iter_count << " iterations in " << time << " sec" << std::endl;
	#endif

	return std::pair<vector_int_ptr, vector_center_ptr<I>>(assignments, centers);
}

template<typename I>
std::pair<vector_vector_int_ptr, vector_center_ptr<I>> kmeans(const RowBlock<I> &data, int k, int p, int dim, int center_type, std::mt19937_64 &rng, int init = 1)
{
	double start, time;
	
	#if DISTRIBUTED
	LOG(INFO) << "Starting k-means";
	#else
	std::cout << "Starting k-means" << std::endl;
	#endif
	#if DISTRIBUTED
	LOG(INFO) << "Initializing centers.."; 
	#else
	std::cout << "Initializing centers.." << std::endl;
	#endif

	vector_center_ptr<I> centers;
	start = GetTime();
	if (init == 0) 
		centers = random_init(data, k, dim, rng);
	else //kmpp
		centers = kmpp_init(data, k, dim, rng);
	time = GetTime() - start;

	#if DISTRIBUTED
	LOG(INFO) << "Initialization done in " << time << " sec";
	#else
	std::cout << "Initialization done in " << time << " sec" << std::endl;
	#endif


	#if DISTRIBUTED
	LOG(INFO) << "Starting Lloyd's iterations";
	#else
	std::cout << "Starting Lloyd's iterations" << std::endl;
	#endif

	bool changed = true;
	vector_vector_int_ptr assignments = std::make_shared<std::vector<std::vector<int>>>(data.size, std::vector<int>(p));
	start = GetTime();
	int iter_count = 0;
	while(changed)
	{
		changed = update_assignments(*centers, data, *assignments);
		update_centers(*centers, data, *assignments, center_type);
		++iter_count;
	}
	time = GetTime() - start;

	#if DISTRIBUTED
	LOG(INFO) << "Finished k-means: " << iter_count << " iterations in " << time << " sec";
	#else
	std::cout << "Finished k-means: " << iter_count << " iterations in " << time << " sec" << std::endl;
	#endif

	return std::pair<vector_vector_int_ptr, vector_center_ptr<I>>(assignments, centers);
}


template<typename I>
void save_data_to_file(const char * filename, const RowBlock<I> &data, vector_int_ptr assignments, Stream *fo = NULL)
{
	libsvmwrite<I>(filename, data, *assignments);
}

template<typename I>
void save_data_to_file(const char *filename, const RowBlock<I> &data, vector_vector_int_ptr assignments, Stream *fo = NULL)
{
	libsvmwrite(filename, data, *assignments);
}


template<typename I>
void save_centers_to_file(Stream *fo, vector_center_ptr<I> centers)
{
	//TODO
}


/*
*	MERGE AND SPLIT HEURISTICS:
*	 Heuristics to merge clusters that are too small and split clusters that are too big 
*
*/

std::vector<int> _ones(int len)
{
	std::vector<int> ones(len);
	for (int i = 0; i < len; ++i) ones[i] = 1;
	return ones;
}

template <typename I>
int merge_and_split(const RowBlock<I> &data, vector_int_ptr assignments, vector_center_ptr<I> centers, real_t lower_bound, real_t upper_bound, int k, int dim, std::mt19937_64 &rng)
{
	//preprocess
	int current_k = k;
	real_t current_lower = lower_bound,
			current_upper = upper_bound;
	int counts[k] ; 
	for (int i = 0; i < k; ++i) counts[i] = 0;
	int nmerge = 0,  nsplit = 0;
	for (int i = 0; i < assignments->size(); ++i)
	{
		++counts[(*assignments)[i]];
	}
	std::vector<int> permutation(k);
	std::vector<int> deleted_indices;
	std::cout << current_k << std::endl;
	//merge
	for (int i = 0; i < k; ++i) permutation[i] = i;
	std::shuffle(permutation.begin(), permutation.end(), rng);
	for (int i: permutation)
	{
		if (counts[i] >= 0 && counts[i] < current_lower)
		{
			//cluster too small. Need to merge it with nearest cluster.
			deleted_indices.push_back(i);
			int new_index = find_closest((*centers)[i], centers, i);
			++nmerge;
			counts[new_index] += counts[i];
			counts[i] = -1; //destroyed cluster
			--current_k;
			std::cout << "Merged cluster " << i << " with cluster " << new_index << std::endl;
			for (int j = 0; j < assignments->size(); ++j)
			{
				if ((*assignments)[j] == i) (*assignments)[j] = new_index;
			}
			//update center:
			update_one_center(*centers, data, (*assignments), new_index);
		}
	}
	//modify upper bound (it gets looser) and leave lower bound intact
	current_upper *= k / current_k;
	
	
	
	//split
	for (int i = 0; i < k; ++i) permutation[i] = i;
	std::shuffle(permutation.begin(), permutation.end(), rng);
	int current_index = k; //starting index for new clusters formed
	for (int i: permutation)
	{ 
		if (counts[i] >= 0 && counts[i] > current_upper)
		{
			//cluster too big. Need to split.
			nsplit++;
			int num_new_clusters_for_this_cluster = counts[i] / current_upper + (counts[i] != current_upper); //ceil
			std::cout << "split cluster " << i << " to clusters " << current_index  << " ... " << (current_index + num_new_clusters_for_this_cluster - 1) << std::endl;
			deleted_indices.push_back(i); //delete this index
			std::vector<int> ones = _ones(num_new_clusters_for_this_cluster);
			std::discrete_distribution<int> dis (ones.begin(), ones.end());
			for (int j = 0; j < assignments->size(); ++j)
			{
				if ((*assignments)[j] == i)  //ramdomly assign to one of the split clusters
				{
					int temp = dis(rng);
					//std::cout << '*' << temp << std::endl;
					(*assignments)[j] = temp + current_index;
				}
			}
			current_index += num_new_clusters_for_this_cluster; //next cluster starts from this index
			/*
			if we split into t clusters, t-1 new clusters are formed
			*/
			current_k += (num_new_clusters_for_this_cluster - 1);	
		}
	}
	
	//pre-clean-up
	int counts1[current_index];
	for (int i = 0; i < current_index; ++i) {counts1[i] = 0; }

	for (int i = 0; i < assignments->size(); ++i)
	{
		int cl = (*assignments)[i];
		++counts1[cl];
	}
	for (int i = 0; i < current_index; ++i) {std::cout << i << ": " << counts1[i] << std::endl; }
	
	
	//clean-up
	//re-map indices >=current_k to deleted indices
	for (int j = 0; j < assignments->size(); ++j)
	{
		int temp = (*assignments)[j] - current_k;
		if (temp >= 0) (*assignments)[j] = deleted_indices[temp];
	}
//	for (int j = 0; j < assignments->size(); ++j)
//		std::cout << (*assignments)[j] << ' ' ;
	#if DISTRIBUTED
	LOG(INFO) << "Merged: " << nmerge << " and split " << nsplit << ".\n Old k: " << k << "; Current k: " << current_k;
	#else
	std::cout << "Merged: " << nmerge << " and split " << nsplit << ".\n Old k: " << k << "; Current k: " << current_k << std::endl;
	#endif
	return current_k;
}



} //namespace dddml



int main()
{
	using namespace dddml;
	using namespace std;
	std::random_device rd; 
	std::mt19937_64 rng(rd());
	dmlc::data::RowBlockContainer<int> rbc = dmlc::data::libsvmread("./rcv.txt");
	RowBlock<int> block = rbc.GetBlock();
	int n = block.size;
	auto rb = block;
	int k = 10;
	int dim = 48000;


	#if 0
	for (int i = 0; i < block.size; ++i)
	{
		real_t dist = 1e230;
		for (int j = 0; j < block.size; ++j)
			if (i != j) dist = std::min(dist, squareDistBetweenRows(block[i], block[j]));
		cout << dist << endl;
	}
	cout << block.size << endl << endl;
	#endif
	
	//#if 0
	auto output = kmeans(block,  k , /*dim */ dim, /* type */ 1, rng, 1);
	std::cout << "-----------------------\n";
	vector_int_ptr assignments = output.first;
	auto centers = output.second;
	//for (auto i : *assignments)
	//	std::cout << i << std::endl;
	int counts[k];
	real_t distances[k];
	for (int i = 0; i < k; ++i) {counts[i] = 0; distances[i] = 0;}

	for (int i = 0; i < assignments->size(); ++i)
	{
		int cl = (*assignments)[i];
		++counts[cl];
		distances[cl] += (*centers)[cl]->squareDist(block[i]);
	}
	for (int i = 0; i < k; ++i)
	{
		std::cout << i << ": " << counts[i] << "; " << std::endl;
	}
	std::cout << "-----------------------\n";
	real_t lb = 0.5 * n / k , ub = 2 * n / k;
	std::cout << "Bounds: " << lb << ' ' << ub << std::endl;
	int new_k = merge_and_split(block, assignments, centers, lb, ub, k, dim, rng);
	
	int counts1[new_k];
	for (int i = 0; i < new_k; ++i) {counts1[i] = 0; }

	for (int i = 0; i < assignments->size(); ++i)
	{
		int cl = (*assignments)[i];
		++counts1[cl];
	}
	for (int i = 0; i < new_k; ++i)
	{
		std::cout << i << ": " << counts1[i] << "; " << std::endl;
	}
	save_data_to_file("./sample_out", block, assignments);
	
}



