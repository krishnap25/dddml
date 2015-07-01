#ifndef RPTS_H
#define RPTS_H

#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "local/row_block.h"
#include "local/data.h"

namespace dddml {

namespace rpt_impl { // Namespace for the implementation

using namespace std;
using namespace dmlc;
using namespace dmlc::data;

enum RouteDirection { LEFT, RIGHT };

struct RPTSplit {
  vector<real_t> d; // The direction for the split.
  real_t t;         // The split point or threshold.

  RPTSplit() : d(0), t(0) {}

  template <typename IndexType>
  RPTSplit(mt19937_64 &rng, size_t dimension, const RowBlock<IndexType> data,
           const vector<IndexType> &idxs)
      : d(dimension) {
    // sample a random direction d
    auto std_normal = normal_distribution<real_t>();
    for (size_t i = 0; i < dimension; i++) {
      d[i] = std_normal(rng);
    }
    // project each sample onto d
    vector<real_t> ps(idxs.size());
    for (int i = 0; i < idxs.size(); ++i) {
      ps[i] = data[idxs[i]].SDot(d.data(), dimension);
    }
    // Pick a random fractile between 1/4 and 3/4
    real_t fractile = std::uniform_real_distribution<real_t>(0.25, 0.75)(rng);
    // Use std::nth_element to find the split threshold
    int findex = int(ps.size() * fractile);
    nth_element(ps.begin(), ps.begin() + findex, ps.end());
    t = ps[findex];
  };

  template <typename IndexType>
  inline auto route(Row<IndexType> row) -> RouteDirection {
    real_t p = row.SDot(d.data(), d.size());
    return (p > t) ? RIGHT : LEFT;
  }
};

template <typename IndexType> struct RPTNode {
  unique_ptr<RPTNode> left;
  unique_ptr<RPTNode> right;
  RPTSplit split;
  vector<IndexType> idxs;

  inline auto depth() -> size_t {
    if (!(left && right)) {
      return 1;
    } else {
      return 1 + max(left->depth(), right->depth());
    }
  }

  inline auto search(Row<IndexType> row) -> const vector<IndexType> * {
    if (!(left && right)) {
      return &idxs;
    } else {
      if (split.route(row) == LEFT) {
        return left->search(row);
      } else {
        return right->search(row);
      }
    }
  };
};

template <typename IndexType>
unique_ptr<RPTNode<IndexType>> make_rptree(mt19937_64 &rng, int dimension,
                                           int n0, RowBlock<IndexType> data,
                                           vector<IndexType> &idxs) {
  unique_ptr<RPTNode<IndexType>> node(new RPTNode<IndexType>());
  if (idxs.size() < n0) {
    node->idxs = idxs;
  } else {
    node->split = RPTSplit(rng, dimension, data, idxs);
    vector<IndexType> left_idxs(0);
    vector<IndexType> right_idxs(0);
    for (auto idx : idxs) {
      if (node->split.route(data[idx]) == LEFT)
        left_idxs.push_back(idx);
      else
        right_idxs.push_back(idx);
    }
    node->left = make_rptree(rng, dimension, n0, data, left_idxs);
    node->right = make_rptree(rng, dimension, n0, data, right_idxs);
  }
  return node;
}

template <typename IndexType>
unique_ptr<RPTNode<IndexType>> make_rptree(mt19937_64 &rng, int dimension,
                                           int n0, RowBlock<IndexType> data) {
  vector<IndexType> idxs(data.size);
  for (size_t i = 0; i < data.size; ++i) {
    idxs[i] = static_cast<IndexType>(i);
  }
  return make_rptree(rng, dimension, n0, data, idxs);
}
};

template <typename IndexType> class RandomPartitionTree {
public:
  RandomPartitionTree(std::mt19937_64 &rng, int dimension, int n0,
                      dmlc::data::RowBlockContainer<IndexType> data)
      : data(data),
        tree(rpt_impl::make_rptree(rng, dimension, n0, data.GetBlock())) {}

  auto depth() -> size_t { return tree->depth(); }

private:
  dmlc::data::RowBlockContainer<IndexType> data;
  std::unique_ptr<rpt_impl::RPTNode<IndexType>> tree;
};
};
#endif
