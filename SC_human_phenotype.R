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


for(i in 1:1000){
 
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
    temp_list_structure$Id <- temp_line
  }
  if(str_detect(temp_line,"name")){
    temp_list_structure$name <- temp_line
  }
  if(str_detect(temp_line, "is_a")){
    temp_list_structure$parent[[temp_is_a]] <- temp_line
    temp_is_a <- temp_is_a + 1
  }
}

#getting rid of first element of hpo_list
hpo_list <- hpo_list[-1]


# l = list(left = list(), right = list())

#for in in length of list, check if ID is A, if not go to next node
# if that node has a child, check the list of all the child nodes
# find (id no node)

#To find ID of node
# function: find -id = function(node = full tree)
# check if its ID = ID you are searching for and node you are searching for again
# if it does, return the node and you are done, if not, check if node has any children
# check if the node is present in the children (for (i in length of children nodes, find ID of child node))
#if none children match, return -1 (to have some control)
# if thing = -1, proceed in a loop, if its a list, return that list
#on working hpo_list, call function on every single node that exists, should return ID on every node that exists

findId <- function(ID_name, tree){
  for(i in 1:length(tree)){
    if(ID_name==tree[i]){
      return(tree[i])
    } 
    else if(ID_name!=tree[i]){
      findId(ID_name, tree$children)
    }
  }
  return(i)
}