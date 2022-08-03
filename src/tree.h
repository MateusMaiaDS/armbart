#include<iostream>
#include <vector>
using namespace std;

// Creating one sequence of values
vector<int> seq_along(int n){
  vector<int> vec_seq;
  for(int i =0; i<n;i++){
    vec_seq.push_back(i);
  }
  return vec_seq;
}

class node{

  // Storing parameters
public:
  int index; // Storing the node index
  vector<int> obs_train; // Storing train observations index that belong to that node
  vector<int> obs_test; // Storing test observations index that belong to that node
  int left; // Index of the left node
  int right; // Index of the right node
  int parent; // Index of the parent
  int depth; //  Depth of the node

  int var; // Variable which the node was split
  double var_split; // Variable which the node was split

  double mu; // Mu parameter from the node
  double loglike; // Loglikelihood of the residuals in that terminal node
public:

  // Getting the constructor
  node(int index_n, vector<int> obs_train_numb, vector<int> obs_test_numb,
       int left_i, int right_i, int depth_i, int var_i, double var_split_i,
       double mu_i){

    index = index_n;
    obs_train = obs_train_numb;
    obs_test = obs_test_numb;
    left = left_i;
    right = right_i;
    depth = depth_i;
    var = var_i;
    var_split = var_split_i;
    mu = mu_i;
  }

  void DisplayNode(){

    cout << "Node Number: " << index << endl;
    cout << "Decision Rule -> Var:  " << var << " & Rule: " << var_split << endl;
    cout << "Left <-  " << left << " & Right -> " << right << endl;

    if(true){
      cout << "Observations train: " ;
      for(int i = 0; i<obs_train.size(); i++){
        cout << obs_train[i] << " ";
      }
      cout << endl;
    }

    if(true){
      cout << "Observations test: " ;
      for(int i = 0; i<obs_test.size(); i++){
        cout << obs_test[i] << " ";
      }
      cout << endl;
    }

  }

  bool isTerminal(){
    return ((left == -1) && (right == -1) );
  }

};



class Tree{

  public:
    // Defining the main element of the tree structure
    vector<node> list_node;

    // Getting the vector of nodes
    Tree(int n_obs_train,int n_obs_test){
      // Creating a root node
      list_node.push_back(node(0,
                               seq_along(n_obs_train),
                               seq_along(n_obs_test),
                               -1, // left
                               -1, // right
                               0, //depth
                               -1, // var
                               -1.1, // var_split
                               0 )); // loglike
    }

    // void DisplayNodesNumber(){
    //   cout << "The tree has " << list_node.size() << " nodes" << endl;
    // }

    void DisplayNodes(){
      for(int i = 0; i<list_node.size(); i++){
        list_node[i].DisplayNode();
      }
      cout << "# ====== #" << endl;
    }


    // Getting terminal nodes
    vector<node> getTerminals(){

      // Defining terminal nodes
      vector<node> terminalNodes;

      for(int i = 0; i<list_node.size(); i++){
        if(list_node[i].isTerminal()==1){
          terminalNodes.push_back(list_node[i]); // Adding the terminals to the list
        }
      }
      return terminalNodes;
    }

    // Getting terminal nodes
    vector<node> getNonTerminals(){

      // Defining terminal nodes
      vector<node> NonTerminalNodes;

      for(int i = 0; i<list_node.size(); i++){
        if(list_node[i].isTerminal()==0){
          NonTerminalNodes.push_back(list_node[i]); // Adding the terminals to the list
        }
      }
      return NonTerminalNodes;
    }


    // Getting the number of n_terminals
    int n_terminal(){

      // Defining the sum value
      int terminal_sum = 0;
      for(int i = 0; i<list_node.size(); i++){
        if(list_node[i].isTerminal()==1){
          terminal_sum++;
        }
      }

      return terminal_sum;
    }

    // Getting the number of non-terminals
    int n_internal(){

      // Defining the sum value
      int internal_sum = 0;
      for(int i = 0; i<list_node.size(); i++){
        if(list_node[i].isTerminal()==0){
          internal_sum++;
        }
      }

      return internal_sum;
    }

    // Get the number of NOG (branches parents of terminal nodes)
    int n_nog(){
        int nog_counter = 0;
          for(int i=0;i<list_node.size();i++){
            if(list_node[i].isTerminal()==0){
              if(list_node[list_node[i].left].isTerminal()==1 && list_node[list_node[i].right].isTerminal()==1){
                nog_counter++;
              }
            }
        }
        return nog_counter;
      }
};

RCPP_EXPOSED_CLASS(node)
RCPP_EXPOSED_CLASS(Tree)

