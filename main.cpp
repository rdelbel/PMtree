//#include <Rcpp.h>
#include <stdio.h>
#include <set>
#include <list>
#include <vector>
#include <cmath>

using namespace std;
//using namespace Rcpp;

struct node {
  int node_num;
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
    node_num(i), score(-1), split(-1), 
    num_events(e), direction(-1), rem_splits(sp),
    rem_obs(ro), left(NULL), right(NULL) 
  { }
  ~node() { delete left; delete right; }
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

  private:
    void find_best_split_node(node *, double **);
    vector<vector<double> > make_big_matrix();
};

vector<vector<double> > tree::make_big_matrix(){
  // initiate big matrix
  int nrow = times.size();
  int ncol = adjustment[0].size() + 2*treatments[0].size() + 1;
  int i, j;
  // Make matrix to input into coxph. Every call to coxph will change the 
  // last columns to the appropriate information apparently this works in
  // C++11
  vector<vector<double> > big_mat(nrow,vector<double>(ncol));

  for(i = 0; i < nrow; i++){
    //big_mat[i][0]=times[i];
    //big_mat[i][1]=status[i];
    for(j = 0; j < adjustment[0].size(); j++){
      big_mat[i][j] = adjustment[i][j];
    }
    for(j = 0; j < treatments[0].size(); j++){
      big_mat[i][j+adjustment[0].size()] = treatments[i][j];
    }
    //dont need to set the reamaining columns. We will change them later.
  }
  return big_mat;
}

/*
// build the bigmatrix. We will use this in find_best_split_node and throw it
// into the coxph function
vector<vector<double> > 
make_big_matrix(vector<vector<double> > adj, vector<vector<double> > tmts, vector<vector<double> > splts) {
  // initiate big matrix
  int nrow = times.size();
  int ncol = adj[0].size() + 2*tmts[0].size() + 1;

  // Make matrix to input into coxph. Every call to coxph will change the 
  // last columns to the appropriate information apparently this works in
  // C++11
  vector<vector<double> > big_mat[nrow][ncol];
  for(int j=0;i<nrow;i++){
    big_mat[i][0]=times[i];
    big_mat[i][1]=status[i];
    for(int j=0; j<adj[0].size();j++){
      big_mat[i][j+2]=adj[i,j];
    }
    for(int j=0; j<tmts[0].size();j++){
      big_mat[i][j+2+adj[0].size()]=tmts[i,j];
    }
    //dont need to set the reamaining columns. We will change them later.
  }

  return big_mat;
}
*/

// find the best split of a single node. assumes all splitting covariates are
// binary. We will call this function to split the root. Later calls will be
// made to the children of a split automatically in the make_split function
void tree::find_best_split_node(node * nd, double ** big_mat) {
  //if no split is viable, cur_max=-1 at end of function call
  double cur_max = -1;
  int cur_max_index = -1;
  int cur_direction = -1;
  double times[nd->rem_obs.size()];
  double status[nd->rem_obs.size()];
  double * mat[nd->rem_obs.size()];
  int pref_size = adjustment.size() + treatments.size();
  set<int>::iterator col, row;
  int i, j;

  // For each remaining split..
  for (col = nd->rem_splits.begin(); col != nd->rem_splits.end(); col++) {
    // Make matrix to input into coxph
    for (row = nd->rem_obs.begin(), i = 0; row != nd->rem_obs.end(); row++, i++) {
      times[i]  = times[*row];
      status[i] = status[*row];
      mat[i] = big_mat[*row];
      mat[i][pref_size] = splits[*row][*col];
      for (j = 0; j < treatments.size(); j++) {
        mat[i][pref_size + 1 +j] = splits[*row][*col] * treatments[*row][j];
      }
    }

    vector<double> base_info; // = coxph(times,status,mat,pref_size+1);
    //  if base model does notconverge, go to next splitting covariate.
    if (base_info[0] = 100) {
      continue;
    }

    vector<double> split_info; // = coxph(times, status, mat,
      //  pref_size+1+treatment.size());

    //if split is viable, get score
    if(!(split_info[0] == 10 || abs(split_info[1]) > 10 || split_info[1] == 0)){
      //split_score is LRT of the interactive terms == 0
      double split_score = 2*(split_info[2]-base_info[2]);
        if (split_score > cur_max) {
          cur_max = split_score;
          cur_max_index = *col;
          //level of covariate that has worse response to treatment goes
          //to has direction 1. Direction 1 will go to right of tree.
          cur_direction = split_info[1] >= 0 ? 1 : 0;
        }
    }
  }

  //  cur_max=-1 if no viable split
  nd->score = cur_max;
  nd->split = cur_max_index;
  nd->direction = cur_direction;
}    

// find best split in the tree. Simply searches for the splitting statistic of
// each node    
list<node *>::iterator find_best_split_tree(list<node *> leaves) {
    list<node *>::iterator it, max_node = leaves.begin();

    for (it = leaves.begin(); it != leaves.end(); it++) {
      if ((*it)->score > (*max_node)->score) max_node = it;
    }

    // if the max score is -1 then there are no more possible splits
    // otherwise return the node with the best split    
    if ((*max_node)->score == -1) return leaves.end();
    else return max_node;
}


// make the split, update the leaves with the two new children nodes, and find
// the best split for them
void make_split(node * leaf, int min_num_events) {
    set<int>::iterator it;
    set<int> r_rem_obs, l_rem_obs, rem_splits = leaf->rem_splits;
    int l_events = 0;
    int r_events = 0;

    // dont try and split by previously split variables. This assumes that all
    // splitting variables are binary
    rem_splits.erase(leaf->split);

    // find the remaining observations of each child node as well as the number
    // of evets. We assume status 0 = censor 1 = event
    for (it = leaf->rem_obs.begin(); it != leaf->rem_obs.end(); it++) {
        if (splits[*it][leaf->split] == leaf->direction) {
            r_rem_obs.insert(* it);
            if (status[*it] == 1) l_events++;
        } else {
            l_rem_obs.insert(* it);
            if (status[*it] == 1) r_events++;
        }
    }

    leaf->left  = new node(leaf->node_num * 2, rem_splits, l_rem_obs, l_events);
    leaf->right = new node((leaf->node_num * 2) + 1, rem_splits, r_rem_obs, r_events);

    // we dont want to split the children if they dont have enough events. in 
    // this case child.score will be -1
    if (l_events >= min_num_events) find_best_split_node(leaf->left);
    if (r_events >= min_num_events) find_best_split_node(leaf->right);

}   

// actually build the tree.
void tree::build_tree(int min_num_splits, int min_num_events) {
    // initiate list of leaves
    list<node *> leaves;
    set<int> splits, obs;
    int i, num_events;

    // initiate root node. Each node will contain the information of the sample
    // before the split as well as information from the split 
    for (i = 0; i <= splits[0].size(); i++) splits.insert(i);
    for (i = 0; i <= status.size();    i++) obs.insert(i);
    // number of events is sum of status vector. Assumes 0=no event, 1=event.
    num_events = accumulate(status.begin(), status.end(), 0);
    root = new node(1, rem_splits, rem_obs, num_events);
    // get the first split information.
    find_best_split_node(root);

    leaves.push_back(root);

    // the best splits have already been found for the new leaves
    // the two possible conditions to stop building the tree are: after the
    // tree has max_num_splits splits, or when each leaf has too few elements.
    // We will contorl for the first condition here, while the second condition
    // is accoutned for in find_best_split_tree 
    // remove the parent node from the list of leaves and add the two children
    // nodes from the leaves
    while (leaves.size() < max_num_splits) {          
        node_to_split = find_best_split_tree(leaves)
        // node_to_split will be null if no more splits (each leaf has too few
        // ovservations or models do not fit properly)
        if (node_to_split == leaves.end()) break;
        else {
          // this will update leaves and find the best split for the new leaves
          make_split(node_to_split);
          leaves.erase(node_to_split);
          leaves.push_back((*node_to_split)->left);
          leaves.push_back((*node_to_split)->right);
        }
    }
    // tree is now built
    num_splits = leaves.size();
}

/*
//todo iterface coxph with our function
// interface our function with r
// build the bigmatrix
// need to sort the rows of our input data in ascending order of time
// Extract information out of a C++ node and make an Rcpp column.
Rcpp::NumericMatrix::Row extract_info(node * nd){
    return Rcpp::NumericMatrix::Row(
            nd->node_num, nd->rem_obs.size(), nd->num_events, nd->split,
            nd->score, nd->direction);
}

//Convert C++ tree to R tree. Use BFS algorithm to traverse through the tree
//and write node information as a line in a matrix
Rcpp::NumericMatrix makeRtree(node * ctree){
    Rcpp::NumericMatrix Rtree(num_splits, 6);
    Q=queue<node *>;
    Q.push(ctree);
    int i=0
    while(!Q.empty()){
        t=Q.front();
        Q.pop();
        Rtree(i,_)=extract_info(t);
        Q.push(t->left);
        Q.push(t->right);
        i++;
    }
    return Rtree;
}

//interface with R. Take R inputs, covert to C++, build tree in CPP, convert
//tree to R, return it.
Rcpp::NumericMatrix survtree(Rcpp::NumericVector tim, Rcpp::IntegerVector sta, 
        Rcpp::NumericMatrix trt, Rcpp::NumericMatrix adj,
        Rcpp::NumericMatrix splt, int mns, int mne){
            // turn Rcpp objects into C++ objects.    
            vector<double> ctim=as<vector<double> >(tim);
            vector<int> csta=as<vector<int> >(sta);
            vector<vector<double> > ctrt(trt.nrow(),vector<double>(trt.ncol()));
            for(int i=0; i<trt.nrow(); i++){
              for(int j=0; j<trt.ncol();j++){
                ctrt[i][j]=trt(i,j);
              }       
            }

           vector<vector<double> > cadj(adj.nrow(),vector<double>(adj.ncol()));
            for(int i=0; i<adj.nrow(); i++){
              for(int j=0; j<adj.ncol();j++){
                cadj[i][j]=adj(i,j);
              }       
            }
             vector<vector<double> > csplt(splt.nrow(),vector<double>(splt.ncol()));
           for(int i=0; i<splt.nrow(); i++){
              for(int j=0; j<splt.ncol();j++){
                csplt[i][j]=splt(i,j);
              }       
            }
    //load the data into the tree
           // load data into the tree
    tree  t(ctim, csta, cadj, ctrt,csplt);
    //build tree
    ...;
    // convert tree into R list of list
    Rcpp::NumericMatrix myRtree=makeRtree(mytree);
    // return list back to R
    return myRtree;
}
*/

//dont need main when combining with R.

int main() {
  return 0;
}
