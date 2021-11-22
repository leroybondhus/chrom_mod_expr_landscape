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

hpo_file <- readLines("hpo_obo.txt", n= -10)

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

#construct dataframe as storage site
hp_terms <- unique(hp_and_omim_id$`HPO Term ID`)
hp_freq <- numeric((length(hp_terms)))
hp_information_content <- data.frame(cbind("hp_terms" = hp_terms, "hp_freq" = hp_freq))

for(i in 1:length(unique(hp_and_omim_id$`OMIM ID Numerical only`))){
  #for row which contains OMIM ID
  #Then finding the hpo term ID
  # Go back to hp_freq vector and increment it by one
  # For each unique OMIM ID
  omim_id_num <- unique(hp_and_omim_id$`OMIM ID Numerical only`)[i]
  # Find the rows that have the OMIM ID and find the HPO Terms
  hpo_terms <- hp_and_omim_id$`HPO Term ID`[which(hp_and_omim_id$`OMIM ID Numerical only`==omim_id_num)]
  for(j in 1:length(hpo_terms)){
    target_hpo <- hpo_terms[i]
    #increment by one
    hp_information_content[which(hp_information_content$hp_terms==target_hpo), "hp_freq"] <- hp_information_content[which(hp_information_content$hp_terms==target_hpo), "hp_freq"] + 1
    # find all parents
  }
}

#helper function that finds parents of a node

## Added a vector_to_store argument to keep cocatenating, if we create a new vector inside the function it will keep 'emptying' the parent nodes everytime the function is called recursively
find_parent <- function(node, node_list, vector_to_store = c()){
  #find all parents of node
  parents_of_current_node <- as.vector(unlist(node$parent))
  #adding parents to return vector
  vector_to_store <- c(vector_to_store, parents_of_current_node)
  if(length(node$parent)==0){
    #condition where root node is the desired node
    return(vector_to_store)
  }
  for(i in 1:length(parents_of_current_node)){
    for(j in 1:length(node_list)){
      if(node_list[[j]]$Id == parents_of_current_node[i]){
        #find the node in hpo list that is the node of the parent of the current node
        if(node_list[[j]]$Id=="HP:0000001"){
          #if we hit root node than return the vector_to_store result
          print(vector_to_store) # IF you change this to return(vector_to_store) doesn't return anything
        }
        find_parent(node_list[[j]],hpo_list, vector_to_store)
        
      }
    }
  }
}

find_parent(hpo_list[[113]],hpo_list)



