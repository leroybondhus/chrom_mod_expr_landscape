#for each unique OMIM ID 
# find all the HPO Terms
# find the parents of the HPO Terms
# Increment the frequency by one

#create data frame of all unique HPO Terms and column of frequency
#Find unique HPO Terms
hpo_terms_unique <- c()
for(i in 1:length(hpo_list)){
  hpo_terms_unique <- c(hpo_list[[i]]$Id, hpo_terms_unique)
}


information_content <- data.frame(cbind("HPO_name" = hpo_terms_unique, "freq" = rep(0,length(hpo_terms_unique))))

#Extract all HP Id's from hpo_list - helpful for index searching later
hp_ids_extracted <- c()
for(i in 1:length(hpo_list)){
  hp_ids_extracted <- c(hpo_list[[i]]$Id, hp_ids_extracted)
}
  
unique_omim_id <- unique(hp_and_omim_id$`Disease ID for link`)

for(i in 1:length(unique_omim_id)){
  #Find all HP terms associated with an omim id or disease
  hp_unique_to_omim <- hp_and_omim_id$`HPO Term ID`[which(hp_and_omim_id$`Disease ID for link`==unique_omim_id[i])]
  # Then we have to find these HP terms in the hpo_list in order to find parental relationships
  for(j in 1:length(hp_unique_to_omim)){
    #Find index in hpo_list which has same name as each HPO Term. We use hp_ids_extracted
    index <- which(hp_ids_extracted==hp_unique_to_omim[j])
    #search for parents
    parents_of_node <- find_parent(hpo_list[[index]],hpo_list)
    #cocatenating parent with child 
    parents_of_node_and_child <- c(hpo_list[[index]]$Id, parents_of_node)
    for(k in 1:length(parents_of_node_and_child)){
      #increment all HP terms in dataframe
      index_2 <- which(information_content$HPO_name==parents_of_node[k])
      information_content[index_2,2] <- as.numeric(information_content[index_2,2]) + 1
    }
  }
}
