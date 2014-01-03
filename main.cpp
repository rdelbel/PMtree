#include "pmtree.h"
#include <iostream>

int main() {
  vector<double> times {306, 455, 1010, 210, 883, 1022, 310, 361, 218, 166, 170, 654, 728, 71, 567, 144, 613, 707, 61, 88 };
  vector<int> status {1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  vector<vector<double> > adjustments;
  vector<vector<double> > treatments{{0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,0}};
  vector<vector<double> > splits{    {1,1,0,0,0,1,1,1,0,0,0,1,1,0,0,1,1,0,0,0}};

  tree pmtree(times, status, adjustments, treatments, splits);
  pmtree.build_tree(10, 10);
  return 0;
}
