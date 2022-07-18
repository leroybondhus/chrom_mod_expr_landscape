#load libraries

library(dplyr)
library(stringr)
library(data.table)
omim_dataset <-  tempfile()
download.file("https://www.omim.org/static/omim/data/mim2gene.txt", omim_dataset)

omim <- read.table(omim_dataset, fill = TRUE, header = TRUE, sep= "\t")
#Renaming
colnames(omim) <- c("MIM","MIM entry","Entrez Gene ID", "Approved Gene Symbol", "Gene Ensemble ID")

### Reading in file with OMIM and HP
hp_and_omim_temp <- tempfile()
download.file("https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/genes_to_phenotype.txt", hp_and_omim_temp)
hp_and_omim_id <- read.table(hp_and_omim_temp, sep = "\t")


### extract the number for respective OMIM ID for each row
omim_id_respective <- as.character(str_extract(hp_and_omim_id[,9],"\\d+"))

### creates a new column without the "OMIM" text in front of digits for OMIM ID

hp_and_omim_id <- data.frame(cbind(hp_and_omim_id), omim_id_respective)
#renaming dataframe
colnames(hp_and_omim_id) <- c("Entrez-gene-id","Entrez-gene-symbol","HPO Term ID","HPO Term Name","Frequency-Raw","Frequency-HPO","Additional Info from GD Source","GD Source","Disease ID for link","OMIM ID Numerical only")

final_df <- merge(hp_and_omim_id, omim, by.x = "Entrez-gene-id",by.y = "Entrez Gene ID")

## BUILDING ACTUAL ONTOLOGY TREE

##reading in line by line, readlines function cannot read in file directly from connection 
#download to computer


hpo_file <- readLines("hpo_obo.txt")

hpo_list <- list()

relationship <- 0 #initializing
current_item_in_list <- 0 # initializing
temp_is_a <- 0 # initializing


for(i in 1:length(hpo_file)){
 
  temp_line <- hpo_file[i]
  
  if(temp_line == "[Term]" ){
    #cannot put temp_list_structure here because you could be hitting a new term and need to replace
    if(current_item_in_list>0){
      hpo_list[[current_item_in_list]] <- temp_list_structure
    }
    current_item_in_list <- current_item_in_list + 1
    temp_list_structure <- list() #resetting to be empty
    temp_list_structure$parent <- list()
    temp_is_a <- 1 #so the first is_a parent will be stored in first element of list "relationship"
  }
  if(str_detect(temp_line,"^id")){
    temp_list_structure$Id <- str_match(temp_line, "id: (.+)")[2]
  }
  if(str_detect(temp_line,"name")){
    temp_list_structure$name <- temp_line
  }
  if(str_detect(temp_line, "is_a")){
    temp_list_structure$parent[[temp_is_a]] <- str_match(temp_line, "is_a: (HP:\\d+).+")[2]
    temp_is_a <- temp_is_a + 1
  }
}

###creating get_node function for more efficient ordering

library(stringr)

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


### Here we write code to compute the information content


#helper function that finds parents of a node

## Added a vector_to_store argument to keep cocatenating, if we create a new vector inside the function it will keep 'emptying' the parent nodes everytime the function is called recursively
find_ancestors <- function(node, tree, vector_to_store = c()){
  ### ASSUMPTION: node is not an ancestor of itself
  if(length(node$parent)==0){
    #condition where root node is the desired node
    return(vector_to_store)
  }
  #find all parents of node
  parents_of_current_node <- as.vector(unlist(node$parent))
  #adding parents to return vector
  vector_to_store <- c(vector_to_store, parents_of_current_node)
  
  for(i in 1:length(parents_of_current_node)){
    parent_node <- get_node(parents_of_current_node[i], tree)
    vector_to_store <- c(find_ancestors(parent_node,tree), vector_to_store)
  
  }
  return(unique(vector_to_store))
}

test_node <- get_node("HP:5000044",hpo_tree)
find_ancestors(test_node, hpo_tree)



##### LOG FREQUENCY

## Extract all ID's in vector form from hpo_list for use later on

#extract all hpo_list ID's
all_hpo_list_id <- c()

for(i in 1:length(hpo_list)){
  all_hpo_list_id[i] <- hpo_list[[i]]$Id
}

phenotypes <- unique(hp_and_omim_id$`HPO Term ID`) ## phenotypes from hp_and_omim_id df & all_hpo_list from HPO database
diseases <- unique(hp_and_omim_id$`OMIM ID Numerical only`)

#finds overlap between two dataframes
all(phenotypes %in% all_hpo_list_id)
#hence phenotypes is a subset of HPO_LIST

hpo_list_overlaps <- hpo_list[(which(is.element(all_hpo_list_id, phenotypes)))]

#extract HPO ID's in order (order matters in subsetting)
hpo_list_overlap_ids <- c()
for(i in 1:length(hpo_list_overlaps)){
  hpo_list_overlap_ids[i] <- hpo_list_overlaps[[i]]$Id
}


#disease_phenotype_df = count number of diseases each phenotype appears in (table form)
disease_phenotype_df <- data.frame(cbind("phenotype" = all_hpo_list_id, "number_of_diseases_appeared" = rep(0,length(all_hpo_list_id))))


for(i in 1:length(diseases)){
  #For each diseases, finding the corresponding phenotypes
  #Find all indexes in the hp_and_omim_id that have the corresponding diseases ID
  indexes_disease <- which(hp_and_omim_id$`OMIM ID Numerical only`==diseases[i])
  #For each indexes, find (1) the corresponding HP Terms and (2) its associated parents
  current_hpo_ids <- c()
  for(j in indexes_disease){
    current_hpo_ids <- c(current_hpo_ids,hp_and_omim_id$`HPO Term ID`[j]) #find all HPO terms
  }
  #Find all ancestors of current_hpo_ids
  ancestors <- c()
  for(k in current_hpo_ids){
    #Find corresponding index in HPO list 
    current_node <- get_node(k,hpo_tree)
    #for all current_hpo_ids, find ancestors than cocatenante everything in one big vector
    temp <- find_ancestors(current_node,hpo_tree)
    ancestors <- c(temp, ancestors)
  }
  #cocatenating ancestors + child and applying unique function to prevent overlaps
  ancestors_and_child <- c(current_hpo_ids, ancestors)
  ancestors_and_child <- unique(ancestors_and_child)
  
  for(m in 1:length(ancestors_and_child)){
    item_id <- which(disease_phenotype_df$phenotype==ancestors_and_child[m])
    new_count <- as.numeric(disease_phenotype_df[item_id,"number_of_diseases_appeared"])+1
    disease_phenotype_df[item_id, "number_of_diseases_appeared"] <- new_count
  }
}

disease_phenotype_df$number_of_diseases_appeared <- as.numeric(disease_phenotype_df$number_of_diseases_appeared)

### Creating log frequency table
log_freq_table <- data.frame(cbind("phenotype" = all_hpo_list_id, "freq" = as.numeric(rep(0,length(all_hpo_list_id)))))

## function log_freq that calculates log frequency

log_freq <- function(x){
   frequency <- (x+1)/length(diseases)
   -1*log(frequency) ###** ASK: what if we hit the case where HP:00001, so log freq of 1 = 0?
}

for(j in 1:nrow(disease_phenotype_df)){
  log_freq_table[j,2] <- log_freq(as.numeric(disease_phenotype_df[j,2]))
}

log_freq_table$freq <- as.numeric(log_freq_table$freq) #technically, it is log freq + 1 (to account for log (0) case)


### write a function called least_common_ancestor that finds least common parent amongst two HPO phenotypes

least_common_ancestor <- function(hpo_1, hpo_2, log_freq_object = log_freq_table){
  #log_freq_table uses colnames: freq and phenotype
  node_one <- get_node(hpo_1,hpo_tree)
  node_two <- get_node(hpo_2,hpo_tree)
  
  parents_of_hpo_one <- find_ancestors(node_one, hpo_tree)
  parents_of_hpo_two <- find_ancestors(node_two, hpo_tree)
  common_parents <- intersect(parents_of_hpo_one, parents_of_hpo_two)
  
  #find the indices of common_parents so we can find their log frequency
  
  current_maximum_frequency<- 0 #sets index
  current_least_common_ancestor <- c()
  for(k in 1:length(common_parents)){
    temp_log_frequency <- log_freq_object$freq[which(log_freq_object$phenotype==common_parents[k])]
    if(temp_log_frequency >= current_maximum_frequency){ ###** ASK: equal because what is we hit HP:0001, log freq = 0?
      current_maximum_frequency <- temp_log_frequency
      current_least_common_ancestor <- common_parents[k]
    }
  }
  
##### **** NOTE: basing least common ancestor out of maximum log frequency, is this accurate? ****
  return(current_least_common_ancestor)
  #if there is a tie, return one element
}

#testing this function
least_common_ancestor(phenotypes[43],phenotypes[45], log_freq_table)
least_common_ancestor(phenotypes[56],phenotypes[77], log_freq_table)
least_common_ancestor(phenotypes[5],phenotypes[7], log_freq_table) # HP: 0000001
least_common_ancestor(phenotypes[5500],phenotypes[67], log_freq_table)
least_common_ancestor(phenotypes[5690],phenotypes[7670], log_freq_table) #All HP: 000002


find_hpo_similarity <- function(hpo_1, hpo_2){
  #specify log freq table
  least_common_parent <- least_common_ancestor(hpo_1,hpo_2, log_freq_table)
  hpo_similarity <-log_freq_table$freq[which(log_freq_table$phenotype==least_common_parent)]
  ### find log frequency of the least common parent
  #### **** I have been calling this log freq, just double checking its information content. **
  
  as.numeric(hpo_similarity)
  
}

find_hpo_similarity(phenotypes[2],phenotypes[4])
find_hpo_similarity(phenotypes[45],phenotypes[56])

### Find OMIM similarity

find_omim_similarity <- function(omim_1, omim_2){
  
  omim_1_hpos <- hp_and_omim_id$`HPO Term ID`[which(hp_and_omim_id$`OMIM ID Numerical only`==omim_1)]
  omim_2_hpos <- hp_and_omim_id$`HPO Term ID`[which(hp_and_omim_id$`OMIM ID Numerical only`==omim_2)]
  sum_1 <- 0 
  
  for(i in 1:length(omim_1_hpos)){
    max_1 <- 0 
    for(j in 1:length(omim_2_hpos)){
      temp_1 <- as.numeric(find_hpo_similarity(omim_1_hpos[i], omim_2_hpos[j]))
      if(temp_1 > max_1){
        max_1 <- temp_1
      }
    }
    sum_1 <- sum_1 + max_1
  }
  
  sum_2 <- 0 
  
  for(k in 1:length(omim_2_hpos)){
    max_2 <- 0 #goes to zero every time you hit loop
    for(m in 1:length(omim_1_hpos)){
      temp_2 <- as.numeric(find_hpo_similarity(omim_2_hpos[k], omim_1_hpos[m]))
      if(temp_2 > max_2){
        max_2 <- temp_2
      }
    }
    sum_2 <- sum_2 + max_2
  }
  
  omim_similarity <- (sum_1 + sum_2)/2
  omim_similarity
  ### CHECK WITH LEROY, normalization
}

#testing the function
find_omim_similarity(diseases[23],diseases[45])
find_omim_similarity(diseases[34], diseases[56])
find_omim_similarity(diseases[3442], diseases[5000])

#reading in chromatin modifier file b/c we are only interested in diseases related to chromatin modifier

chromatin_modifier_disease <- read.csv("chromatin_modifier_disease.csv")
library(stringr)
chromatin_omim_ids <- c()
for(i in 1:nrow(chromatin_modifier_disease)){
  chromatin_omim_ids <- c(chromatin_omim_ids,str_extract_all(chromatin_modifier_disease$OMIM.ID[i], pattern = "(\\d+)")[[1]])
}

chromatin_omim_ids <- chromatin_omim_ids[!is.na(chromatin_omim_ids)]
chromatin_omim_ids <- intersect(chromatin_omim_ids,hp_and_omim_id$`OMIM ID Numerical only`)
n <- length(chromatin_omim_ids)

#extracting disease names from chromatin ID's
gene_names <- c()

for(i in 1:length(chromatin_omim_ids)){
  for(j in 1:nrow(chromatin_modifier_disease)){
    if(!is.na(chromatin_modifier_disease$OMIM.ID[j])){
      if(str_detect(as.character(as.character(chromatin_modifier_disease$OMIM.ID[j])),chromatin_omim_ids[i])){
        #if more than one match gene, select one for display of data
        gene_names[i] <- chromatin_modifier_disease$Gene.Name[j]
      }
    }
  }
}

##### Creating a Table with these results:

omim_similarity_table <- matrix(rep(0,n*n), nrow = n, ncol = n)


for(i in 1:n){
  
  for(j in 1:n){
   
  temp <- (find_omim_similarity(chromatin_omim_ids[i], chromatin_omim_ids[j]))/(max(find_omim_similarity(chromatin_omim_ids[i],chromatin_omim_ids[i]), find_omim_similarity(chromatin_omim_ids[j], chromatin_omim_ids[j])))
  omim_similarity_table[i,j] <- temp
   
  }
}

rownames(omim_similarity_table) <- gene_names
colnames(omim_similarity_table) <- gene_names
diag(omim_similarity_table) <- NA
heatmap(omim_similarity_table, scale = "none")
