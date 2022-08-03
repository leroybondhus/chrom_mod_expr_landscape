
library(stringr)

convert_id <- function(hpo_id){
  id_only <- as.numeric(str_match(hpo_id, pattern = ".*:(\\d+)")[[2]])
  return(id_only)
}

### to ensure efficient ordering, make sure to add nodes to tree randomly
set.seed(123)
rand_seq <- sample(1:length(hpo_list))
hpo_tree <- list()

### initialize hpo_tree with first node and add left and right children
hpo_tree <- hpo_list[[(rand_seq[1])]]
hpo_tree$left <- NA
hpo_tree$right <- NA
###now loop through adding each node
for(i in 1:length(rand_seq)){
  hpo_tree <- add_to_tree(hpo_list[[rand_seq[i]]], hpo_tree) #hpo_tree? 
}
#### add_to_tree function to use
add_to_tree <- function(node, tree){
  if(node$Id < tree$Id){ #PROBLEM: if we replace tree$Id with HPO Id, what is "tree" in calling the add_tree, because there is a tree$left and a tree$right, you need it to be some sort of a node
    if(all(is.na(tree$left))){
      node$left <- NA
      node$right <- NA
      tree$left <- node
      return(tree)
    } else {
      tree$left <- add_to_tree(node, tree$left)
      return(tree) 
    }
  } else {
    if(all(is.na(tree$right))){
      node$left <- NA
      node$right <- NA
      tree$right <- node
      return(tree)
    } else {
      tree$right <- add_to_tree(node, tree$right)
      return(tree)
    }
  }
}

get_node <- function(node_id,tree){ 
  #gets node of HPO term from the hpo_tree
  #check if tree_id == node_id 
  #(tree left = NULL, tree right = NULL)
  #if equal --> create a new node called return node, which we set as equal to the tree 
  #return (tree)
  #if not the case
  #if hpo_id < tree_id
  #get_node(node_id, tree$left)
  #else get_node(node_id, tree$right)
  
  if(node_id == tree$Id){
    tree$left <- NULL
    tree$right <- NULL
    return(tree)
  } else {
    if(node_id < tree$Id){
      get_node(node_id, tree$left)
    } else {
      get_node(node_id, tree$right)
    }
  }
  
}

### testing get_node
get_node("HP:5000044",hpo_tree)
