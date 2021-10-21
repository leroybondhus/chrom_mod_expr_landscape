#load libraries

library(dplyr)
omim_dataset <-  tempfile()
download.file("https://www.omim.org/static/omim/data/mim2gene.txt", omim_dataset)

omim <- read.csv(omim_dataset, fill = TRUE, header = TRUE)

#extracts omim ids with corresponding gene ID's
omim_ensg <-  str_extract(omim[,1],"ENSG\\d+")
omim_id <- str_extract(omim[,1],"\\d+")

omim_cleaned <- data.frame(cbind("OMIM_ID" = omim_id, "OMIM_ENSG" = omim_ensg))
# remove NA's

omim_cleaned <- omim_cleaned[!is.na(omim_cleaned$OMIM_ID) & !is.na(omim_cleaned$OMIM_ENSG), ]


gene_and_omim_temp <- tempfile()
download.file("http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt", gene_and_omim_temp)
gene_and_omim_id <- read.table(gene_and_omim_temp, fill = TRUE, row.names = NULL, sep = "\t")
#eliminate first column because that is entrez
gene_and_omim_id <- gene_and_omim_id[,-1]

#combine the first three columns into one

gene_and_omim_id <- unite(gene_and_omim_id, "text", c(1,2,3))[,1]

#for each row extract HP and OMIM term

hp_id <- str_extract(gene_and_omim_id,"HP:\\d+")
#extracts only digits portion
hp_id <- str_extract(hp_id, "\\d+")
omim_id <- str_extract(gene_and_omim_id, "OMIM:\\d+")
#extract only digits portion
omim_id <- str_extract(omim_id, "\\d+")

hp_and_omim_combined <- data.frame(cbind("HP" = hp_id, "OMIM" = omim_id))
hp_and_omim_combined <- hp_and_omim_combined[(!is.na(hp_and_omim_combined$HP) & !is.na(hp_and_omim_combined$OMIM)),]
hp_and_omim_combined <- unique(hp_and_omim_combined)
colnames(hp_and_omim_combined) <- c("HP","OMIM")

# create a new data frame called gene_hp_omim to add corresponding HP terms to each OMIM ID

gene_hp_omim <- data.frame(cbind(omim_cleaned,"HP" = rep("0",nrow(omim_cleaned)))) #last column for matching HP terms

## creating an object called temp_stored to store the HP terms

for(i in 1:nrow(gene_hp_omim)){
  temp_stored <- "0"
  if(any(hp_and_omim_combined$OMIM ==gene_hp_omim$OMIM_ID[i])){
    index <- which(hp_and_omim_combined$OMIM ==gene_hp_omim$OMIM_ID[i])
    temp_stored <- hp_and_omim_combined[(index[1]),"HP"]
    if(length(index) > 1){
      for(j in 2:length(index)){
        temp_stored <- paste(temp_stored, hp_and_omim_combined[(index[j]),"HP"])
      }
    }
  }
  gene_hp_omim[i,"HP"] <- temp_stored
}