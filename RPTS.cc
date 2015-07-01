#include <algorithm>
#include <iostream>
#include <random>

#include "RPTS.h"
#include "local/io.h"
#include "local/libsvm_parser.h"

using namespace std;
using namespace dddml;
using namespace dddml::rpt_impl;
using namespace dmlc;
using namespace dmlc::data;

int main() {
  mt19937_64 rng(1);
  RandomPartitionTree<int> rpt(rng, 784, 100, libsvmread("./mnist.txt"));
  cout << rpt.depth() << endl;
}
