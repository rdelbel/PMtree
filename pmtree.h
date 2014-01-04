//#include <Rcpp.h>
#include <stdio.h>
#include <set>
#include <list>
#include <vector>
#include <cmath>
#include <numeric>

using namespace std;
//using namespace Rcpp;

struct node {
  int node_num;
  int is_split;         // has the node been split? 0=no 1=yes
  double score;        // Node's score
  int split;           // Which feature this node splits
  int num_events; 
  int direction;       // direction currently is set to 1 for the side where 
                       //   trt2 is worse than trt1
  set<int> rem_splits; // Remaining splits
  set<int> rem_obs;    // remaining observations

  node * left;
  node * right;

  node(int i, set<int> sp, set<int> ro, int e) : 
    node_num(i),is_split(0), score(-1), split(-1), 
    num_events(e), direction(-1), rem_splits(sp),
    rem_obs(ro), left(NULL), right(NULL) 
  { }
 // ~node() { delete left; delete right; }
};

class tree {
  public:
    int num_splits;

    vector<double> times;
    vector<int> status;
    vector<vector<double> > adjustment;
    vector<vector<double> > treatments;
    vector<vector<double> > splits;

    node * root;

    tree(vector<double> tm, vector<int> st,
        vector<vector<double> > ad,
        vector<vector<double> > tr,
        vector<vector<double> > sp) :
      num_splits(0), times(tm), status(st),
      adjustment(ad), treatments(tr), splits(sp), root(NULL)
    { }
    ~tree() { delete root; }

    void build_tree(int, int);
    void print_tree(bool);

  private:
    void find_best_split_node(node *);
    list<node *>::iterator find_best_split_tree(list<node *> &);
    void make_split(node *, int, bool);
};

