#include <algorithm>
#include <iostream>
#include <random>

#include "RPTS.h"
#include "local/io.h"
#include "local/libsvm_parser.h"
#include "kmeans_helper.h"

using namespace std;
using namespace std::placeholders;
using namespace dddml;
using namespace dddml::rpt_impl;
using namespace dmlc;
using namespace dmlc::data;

auto query_dist(function<int(Row<int>)> f, RowBlockContainer<int> data,
                Row<int> row) -> real_t {
  int nnidx = f(row);
  return sqrt(squareDistBetweenRows(data.GetBlock()[nnidx], row));
}

auto mean_dist(function<int(Row<int>)> f, const RowBlockContainer<int> &data,
               int start) -> real_t {
  real_t total_dist = 0.0;
  for (int i = start; i < data.Size(); ++i) {
    total_dist += query_dist(f, data, data.GetBlock()[i]);
  }
  return total_dist / (data.Size() - start);
}

int main() {
  mt19937_64 rng(1);
  auto data = libsvmread("./mnist.txt");

  vector<int> idxs(1000);
  for (size_t i = 0; i < idxs.size(); ++i)
    idxs[i] = i;

  auto random_idx = [&rng, idxs](Row<int> row) {
    return uniform_int_distribution<int>(0, idxs.size())(rng);
  };

  RandomPartitionTree<int> rpt_exhaustive(rng, 784, idxs.size(), data, idxs);
  auto exhaustive_nn =
      bind(&RandomPartitionTree<int>::find_nn, &rpt_exhaustive, _1);

  RandomPartitionTree<int> rpt(rng, 784, 50, data, idxs);
  auto small_nn = bind(&RandomPartitionTree<int>::find_nn, &rpt, _1);

  cout << "Random dist =        " << mean_dist(random_idx, data, 1000) << endl;
  flush(cout);
  cout << "approx nn dist =     " << mean_dist(small_nn, data, 1000) << endl;
  flush(cout);
  cout << "exhaustive nn dist = " << mean_dist(exhaustive_nn, data, 1000)
       << endl;
}
