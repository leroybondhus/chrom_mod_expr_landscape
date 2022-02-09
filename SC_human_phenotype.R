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

{##!## Start: cannot access file on my system ##!##
##!## should add lines to download file locally if possible ##!##
hpo_file <- readLines("hpo_obo.txt", n= -10)
}##!## Stop  ##!##

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

### Here we write code to compute the information content


#helper function that finds parents of a node

## Added a vector_to_store argument to keep cocatenating, if we create a new vector inside the function it will keep 'emptying' the parent nodes everytime the function is called recursively
find_ancestors <- function(node, node_list, vector_to_store = c()){
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
    for(j in 1:length(node_list)){
      
      if(node_list[[j]]$Id == parents_of_current_node[i]){
        #find the node in hpo list that is the node of the parent of the current node
        if(node_list[[j]]$Id=="HP:0000001"){
          #if we hit root node than return the root node
          return(vector_to_store)
        }
        vector_to_store <- c(find_ancestors(node_list[[j]],hpo_list), vector_to_store)
      }
    }
  }
  return(unique(vector_to_store))
}

find_ancestors(hpo_list[[1187]],hpo_list)



##### LOG FREQUENCY

## Extract all ID's in vector form from hpo_list for use later on

#extract all hpo_list ID's
all_hpo_list_id <- c()

for(i in 1:length(hpo_list)){
  all_hpo_list_id[i] <- hpo_list[[i]]$Id
}

phenotypes <- unique(hp_and_omim_id$`HPO Term ID`)
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
    index_in_hpo_list <- which(all_hpo_list_id == current_hpo_ids[k])
    index_in_hpo_list <- index_in_hpo_list[1] # in case there is more than one match
    #for all current_hpo_ids, find ancestors than cocatenante everything in one big vector
    temp <- find_ancestors(hpo_list[index_in_hpo_list],hpo_list)
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
  x <- x + 1 #in case x is 0 
  result <- -1*log((x/length(diseases)))
  result
}

for(j in 1:nrow(disease_phenotype_df)){
  log_freq_table[j,2] <- log_freq(as.numeric(disease_phenotype_df[j,2]))
}

log_freq_table$freq <- as.numeric(log_freq_table$freq) #technically, it is log freq + 1 (to account for log (0) case)


### write a function called least_common_ancestor that finds least common parent amongst two HPO phenotypes

least_common_ancestor <- function(hpo_1, hpo_2){
  
  #the way the find_ancestors works is that the input needs to correspond to respective position in hpo_list
  #so need to find the "index" that HPO terms 1 and 2 appear in 
  n1 <- which(all_hpo_list_id==hpo_1)
  n2 <- which(all_hpo_list_id==hpo_2)
  
  parents_of_hpo_one <- find_ancestors(hpo_list[[n1]], hpo_list)
  parents_of_hpo_two <- find_ancestors(hpo_list[[n2]], hpo_list)
  common_parents <- intersect(parents_of_hpo_one, parents_of_hpo_two)
  
  max_num_parents <- -1 
  max_num_parents_identity <- common_parents[1]
  
  for(i in 1:length(common_parents)){
    index <- which(all_hpo_list_id==common_parents[i])
    temp_num_parents <- length(find_ancestors(hpo_list[[index]],hpo_list))
    if(temp_num_parents > max_num_parents){
      max_num_parents <- temp_num_parents
      max_num_parents_identity <- common_parents[i]
    }
  }
  return(max_num_parents_identity)
}



find_hpo_similarity <- function(hpo_1, hpo_2){
  least_common_parent <- least_common_ancestor(hpo_1,hpo_2)
  hpo_similarity <-log_freq_table$freq[which(log_freq_table$phenotype==least_common_parent)]
  as.numeric(hpo_similarity)
  
}

### Find OMIM similarity

find_omim_similarity <- function(omim_1, omim_2){
  
  omim_1_hpos <- hp_and_omim_id$`HPO Term ID`[which(hp_and_omim_id$`OMIM ID Numerical only`==omim_1)]
  omim_2_hpos <- hp_and_omim_id$`HPO Term ID`[which(hp_and_omim_id$`OMIM ID Numerical only`==omim_2)]
  sum_1 <- 0 
  
  for(i in 1:length(omim_1_hpos)){
    max_1 <- 0 
    for(j in 1:length(omim_2_hpos)){
      temp_1 <- find_hpo_similarity(omim_1_hpos[i], omim_2_hpos[j])
      if(temp_1 > max_1){
        max_1 <- temp_1
      }
      sum_1 <- sum_1 + max_1
    }
  }
  
  sum_2 <- 0 
  
  for(k in 1:length(omim_2_hpos)){
    max_2 <- 0 #goes to zero every time you hit loop
    for(m in 1:length(omim_1_hpos)){
      temp_2 <- find_hpo_similarity(omim_2_hpos[k], omim_1_hpos[m])
      if(temp_2 > max_2){
        max_2 <- temp_2
      }
      sum_2 <- sum_2 + max_2
    }
  }
  
  omim_similarity <- (sum_1 + sum_2)/2
  omim_similarity
}

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

##### Creating a Table with these results:

omim_similarity_table <- matrix(rep(0,n*n), nrow = n, ncol = n)
rownames(omim_similarity_table) <- chromatin_omim_ids
colnames(omim_similarity_table) <- chromatin_omim_ids
for(i in 1:n){
  for(j in 1:n){
    omim_similarity_table[i,j] <- find_omim_similarity(chromatin_omim_ids[i], chromatin_omim_ids[j])
  }
}

