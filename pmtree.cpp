#include "pmtree.h"
#include "coxph.cpp"
#include <iostream>

// find the best split of a single node. assumes all splitting covariates are
// binary. We will call this function to split the root. Later calls will be
// made to the children of a split automatically in the make_split function
void tree::find_best_split_node(node * nd) {
  //if no split is viable, cur_max=-1 at end of function call
  double cur_max = -1;
  int cur_max_index = -1;
  int cur_direction = -1;
  double ctime[nd->rem_obs.size()];
  int cstat[nd->rem_obs.size()];
  double ** mat = new double *[adjustment.size() + treatments.size()*2 + 1];
  for (int j = 0; j < adjustment.size() + treatments.size()*2 + 1; j++) {
    mat[j] = new double[nd->rem_obs.size()];
  }
  int pref_size = adjustment.size() + treatments.size();
  set<int>::iterator splt, person;
  int i, j;

  // For each remaining split..
  for (splt = nd->rem_splits.begin(); splt != nd->rem_splits.end(); splt++) {
    // Make matrix to input into coxph
    for (person = nd->rem_obs.begin(), i = 0; person != nd->rem_obs.end(); person++, i++) {
      ctime[i] = times[*person];
      cstat[i] = status[*person];
      for (j = 0; j < adjustment.size(); j++) {
        mat[j][i] = adjustment[j][*person];
      }
      for (j = 0; j < treatments.size(); j++) {
        mat[j + adjustment.size()][i] = treatments[j][*person];
      }
      mat[adjustment.size() + treatments.size()][i] = splits[*splt][*person];
      for (j = 0; j < treatments.size(); j++) {
        mat[adjustment.size() + treatments.size() + 1 + j][i]
          = splits[*splt][*person] * treatments[j][*person];
      }
    }

    vector<double> base_info = coxph(ctime, cstat, mat, 
                                     nd->rem_obs.size(), 
                                     adjustment.size() + treatments.size() + 1);
    //  if base model does notconverge, go to next splitting covariate.
    if (base_info[0] == 1000) continue;

    vector<double> split_info = coxph(ctime, cstat, mat,
                                      nd->rem_obs.size(),
                                      adjustment.size() + treatments.size()*2 + 1);

    //if split is viable, get score
    bool flag = false;
    for (int i = 0; i < treatments.size(); i++) {
      flag = flag || (abs(split_info[adjustment.size() + treatments.size() + 1 + i]) > 10)
                  || (split_info[adjustment.size() + treatments.size() + 1 + i] == 0);
    }
    if ((split_info[0] != 1000) && !flag) {
      //split_score is LRT of the interactive terms == 0
      double split_score = 2*(split_info[1]-base_info[1]);
        if (split_score > cur_max) {
          cur_max = split_score;
          cur_max_index = *splt;
          //level of covariate that has worse response to treatment goes
          //to has direction 1. Direction 1 will go to right of tree.
          cur_direction = split_info[2] >= 0 ? 1 : 0;
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
list<node *>::iterator tree::find_best_split_tree(list<node *> & leaves) {
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
void tree::make_split(node * leaf, int min_num_events) {
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
        if (splits[leaf->split][*it] == leaf->direction) {
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
void tree::build_tree(int max_num_splits, int min_num_events) {
    // initiate list of leaves
    list<node *> leaves;
    list<node *>::iterator node_to_split;
    set<int> splts, obs;
    int i, num_events;

    // initiate root node. Each node will contain the information of the sample
    // before the split as well as information from the split 
    for (i = 0; i < splits.size(); i++) splts.insert(i);
    for (i = 0; i < status.size(); i++) obs.insert(i);
    // number of events is sum of status vector. Assumes 0=no event, 1=event.
    num_events = accumulate(status.begin(), status.end(), 0);

    root = new node(1, splts, obs, num_events);
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
        node_to_split = find_best_split_tree(leaves);
        // node_to_split will be null if no more splits (each leaf has too few
        // ovservations or models do not fit properly)
        if (node_to_split == leaves.end()) break;
        else {
          // this will update leaves and find the best split for the new leaves
          make_split(*node_to_split, min_num_events);
          leaves.push_back((*node_to_split)->left);
          leaves.push_back((*node_to_split)->right);
          node_to_split = leaves.erase(node_to_split);
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
