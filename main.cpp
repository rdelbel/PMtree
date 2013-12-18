cclass tree{
    struct node{
        double score;
        int split;
        int num_events;
        int node_num;
        // direction currently is set to 1 for the side where trt2 is worse
        // than trt1
        int direction;
        set<int> rem_splts;
        set<int> rem_obs;
        node * left;
        node * right;
    }
    int num_splits;
    struct node * root;
    vector<double> times;
    vector<int> status;

    vector<vector<double> > adjustment;
    vector<vector<double> > treatments;
    vector<vector<double> > splits;

    tree(vector<double>, vector<int>,
            vector<vector<double> >,
            vector<vector<double> >,
            vector<vector<double> >);

    void build_tree(int, int);

// find the best split of a single node. assumes all splitting covariates are
// binary. We will call this function to split the root. Later calls will be
// made to the children of a split automatically in the make_split function
void find_best_split_node(node * nd, double ** big_mat){
    //if no split is viable, cur_max=-1 at end of function call
    double cur_max=-1;
    int cur_max_index=-1;
    int cur_direction=-1;
    double times[rem_obs.size()];
    double status[rem_obs.size()];
    double ** mat[rem_obs.size()];
    pref_size = adjustment.size() + treatments.size();
    // For each remaining split..
    for (col = rem_splts.begin(); col != rem_splts.end(); col++) {
       // Make matrix to input into coxph
       for (row = rem_obs.begin(), i=0; row != rem_obs.end(); row++, i++) {
           times[i]  = times[*rem_obs];
           status[i] = status[*rem_obs];
           mat[i] = big_mat[*row];
           mat[i][pref_size] = splits[*row][*col];
           for (j = 0; j < treatments.size(); j++) {
               mat[i][pref_size + 1 +j] = splits[*row][*col]*treatment[*row][j];
           }
        }
        vector<double> base_info = coxph(times,status,mat,pref_size+1);
        //  if base model does notconverge, go to next splitting covariate.
        if(base_info[0]=100){
            continue;
        }
        vector<double> split_info = coxph(times, status, mat,
                pref_size+1+treatment.size());
        //if split is viable, get score
        if(!(split_info[0]==10 || abs(split_info[1])>10 || split_info[1]==0)){
            //split_score is LRT of the interactive terms == 0
            split_score=2*(split_info[2]-base_info[2])
            if(split_score>cur_max){
                cur_max=split_score;
                cur_max_index=* col;
                //level of covariate that has worse response to treatment goes
                //to has direction 1. Direction 1 will go to right of tree.
                cur_direction=split_info[1]>=0?1:0;
            }
        }
    }
    //  cur_max=-1 if no viable split
    nd->score=cur_max;
    nd->split=cur_max_index;
    nd->direction=cur_direction;
}    
    
// find best split in the tree. Simply searches for the splitting statistic of
// each node    
list<node *>::iterator find_best_split_tree(list<node *> leaves){
    max_node= node *;
    list<node *>::iterator it;
    for(it=leaves.begin(); it!=leaves.end(); it++){
        if (it == leaves.begin() || it->score>max_node){
            max_node = * it;
        }
    }
    // if the max score is -1 then there are no more possible splits
    if(max_node->score==-1){
        return NULL;
    // otherwise return the node with the best split    
    }else{
        return max_node;
    }
}

// make the split, update the leaves with the two new children nodes, and find
// the best split for them
void make_split(list<node *>::iterator leaf_it, node * leaves){
  // initiate variables
    node * leaf = *leaf_it;
    templ = malloc(sizeof(node));
    // dont try and split by previously split variables. This assumes that all
    // splitting variables are binary
    templ->rem_splits=leaf->rem_splits;
    templ->rem_splits=templ->rem_splits.erase(leaf->split);
    templ->rem_obs=set<int>();
    templ->score=-1;
    templ->split=-1;
    templ->direction=-1;
    templ->node_num=leaf->node_num*2
    templ->left=NULL;
    templ->right=NULL;
    tempr = malloc(sizeof(node));
    // dont try and split by previously split variables. This assumes that all
    // splitting variables are binary
    tempr->rem_splits=leaf->rem_splits;
    tempr->rem_splits=tempr->rem_splits.erase(leaf->split);
    tempr->rem_obs=set<int>();
    tempr->score=-1;
    tempr->split=-1;
    tempr->direction=-1;
    tempr->node_num=leaf->node_num*2+1
    tempr->left=NULL;
    tempr->right=NULL;
    // find the remaining observations of each child node as well as the number
    // of evets. We assume status 0 = censor 1 = event
    set<int>::iterator it;
    for(it=leaf->rem_obs.begin(), int left_events=0, int right_events=0;
            it!=leaf->rem_obs.end(); it++){
        if(splits[* it][leaf->split]==leaf->direction){
            tempr->rem_obs.insert(* it);
            if(status[* it]==1){
                left_events++;
            }
        }else{
            templ->rem_obs.insert(* it);
            if(status[* it]==1){
                right_events++;
            }
        }
    }
    templ->num_events=left_events;
    tempr->num_events=right_events;
    // we dont want to split the children if they dont have enough events. in 
    // this case child.score will be -1
    if(left_events>=min_num_events){
        find_best_split_node(templ);
    }
    if(right_events>=min_num_events){
        find_best_split_node(tempr);
    }
    leaf->left=templ;
    leaf->right=tempr;
    // remove the parent node from the list of leaves and add the two children
    // nodes from the leaves
    leaves.erase(leaf_it);
    leaves.push_back(leaf->left);
    leaves.push_back(leaf->right);

}   

tree::tree(vector<double> tm, vector<int> st,
            vector<vector<double> > ad,
            vector<vector<double> > tr,
            vector<vector<double> > sp) {
    times = tm;
    status = st;
    adjustment = ad;
    treatments = tr;
    splits = sp;
    root = NULL;
}

// build the bigmatrix. We will use this in find_best_split_node and throw it
// into the coxph function
vector<vector<double> > make_big_matrix(vector<vector<double> > adjustment,
vector<vector<double> > treatments, vector<vector<double> > splits){
    // initiate big matrix
    nrow=times.size()
    ncol = adjustment[0].size() + 2*treatments[0].size() + 1
    // Make matrix to input into coxph. Every call to coxph will change the 
    // last columns to the appropriate information apparently this works in
    // C++11
    vector<vector<double> > big_mat[nrow][ncol];
    for(int j=0;i<nrow;i++){
        big_mat[i][0]=times[i];
        big_mat[i][1]=status[i];
        for(int j=0; j<adjustment[0].size();j++){
            big_mat[i][j+2]=adjustment[i,j];
        }
        for(int j=0; j<treatments[0].size();j++){
            big_mat[i][j+2+adjustment[0].size()]=treatments[i,j];
        }
        //dont need to set the reamaining columns. We will change them later.
    }
}


// actually build the tree. Return the pointer to the node.
node * tree::build_tree(min_num_splits, max_num_events){
    // initiate list of leaves
    list <node *> leaves= NULL;
    // initiate root node. Each node will contin the information of the sample
    // before the split as well as information from the split 
    root=malloc(sizeof(node));
    root->rem_splits=set<int>;
    for(int i=0; i<=splits[0].size();i++){
        root->rem_splits.insert(i);
    }
    root->rem_obs=set<int>;
    for(int i=0; i<=status.size();i++){
        root->rem_splits.insert(i);
    }
    // number of events is sum of status vector. Assumes 0=no event, 1=event.
    root->num_events=accumulate(status.begin(),status.end(),0);
    root->score=-1;
    root->split=-1;
    root->direction=-1;
    root->node_num=1;
    // get the first split information.
    find_best_split_node(root);
    // add the root to the list of leaves. After this first split, the list of 
    // leaves will be managed by the make_split function
    leaves.push_back(root);
    //make first split
    make_split(leaves)
    // the children of root are now in leaves
    // the best splits have already been found for the new leaves
    // the two possible conditions to stop building the tree are: after the
    // tree has max_num_splits splits, or when each leaf has too few elements.
    // We will contorl for the first condition here, while the second condition
    // is accoutned for in find_best_split_tree 
    while(leaves.size()<max_num_splits){          
        node_to_split=find_best_split_tree(leaves)
        // node_to_split will be null if no more splits (each leaf has too few
        // ovservations or models do not fit properly)
        if(node_to_split==NULL){
          break;
        }else{
          // this will update leaves and find the best split for the new leaves
          make_split(node_to_split);
        }
    }
    // tree is now built. Return root which is a pointer to the root node of 
    // the tree.
    num_splits=leaves.size();
    return root;
}
    
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
    vector<vector<double> > ctrt[trt.nrow][trt.ncol];
    for(int i=1; i<trt.nrow(); i++){
       ctrt[i]=Rcpp::as<vector<double> >(trt(i_))
    }
   vector<vector<double> > cadj[adj.nrow][adj.ncol];
    for(int i=1; i<adj.nrow(); i++){
       cadj[i]=Rcpp::as<vector<double> >(adj(i_))
    }
     vector<vector<double> > csplt[splt.nrow][splt.ncol];
    for(int i=1; i<splt.nrow(); i++){
       csplt[i]=Rcpp::as<vector<double> >(splt(i_))
    }
    //load the data into the tree
    tree(Rcpp::as<vector<double> > tim, Rcpp::as<vector<int> > sta,
            cadj, ctrt, csplt);
    //build tree
    node * tree mytree=build_tree(mns, mne);
    // convert tree into R list of list
    Rcpp::NumericMatrix myRtree=makeRtree(mytree);
    // return list back to R
    return myRtree;
}

int main(){
    return 0;
}
