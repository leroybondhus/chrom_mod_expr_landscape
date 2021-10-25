#load libraries

library(dplyr)
library(stringr)
library(data.table)
omim_dataset <-  tempfile()
download.file("https://www.omim.org/static/omim/data/mim2gene.txt", omim_dataset)

omim <- read.table(omim_dataset, fill = TRUE, header = TRUE, sep= "\t")
#extracting column with GENE ENSEMBLE and OMIM ID
omim <- omim[,c(1,5)]
#removing rows without corresponding GENE ENSEMBL ID
row_no_ensembl_id <- omim$X.2!=""
omim <- omim[row_no_ensembl_id, ]




### Reading in file with OMIM and HP
hp_and_omim_temp <- tempfile()
download.file("https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/genes_to_phenotype.txt", hp_and_omim_temp)
hp_and_omim_id <- fread(hp_and_omim_temp, sep = "\t")

### extract the number for respective OMIM ID for each row
omim_id_respective <- as.character(str_extract(hp_and_omim_id[,9],"\\d+"))


### creates a new column to attach respective ENSG ID's to hp_and_omim_id dataframe

hp_and_omim_id <- data.frame(cbind(hp_and_omim_id, "ENSG" = rep(0,nrow(hp_and_omim_id))))

#### for each omim id find corresponding ENSG ID

for(i in 1:length(omim_id_respective)){
  if(any(as.character(omim$X100050)==omim_id_respective[i])){
    similar_index <- which(as.character(omim$X100050) == omim_id_respective[i])
    if(length(similar_index)>1){
      temp <- paste(omim$X.2[similar_index], sep = ",")
      correspondent_ensg <- c(correspondent_ensg, temp)
      hp_and_omim_id[i,"ENSG"] <- correspondent_ensg
    } else {
      correspondent_ensg <- omim$X.2[similar_index]
      hp_and_omim_id[i,"ENSG"] <- correspondent_ensg
    }
  }
}

#renaming dataframe
colnames(hp_and_omim_id) <- c("Entrez-gene-id","Entrez-gene-symbol","HPO Term ID","HPO Term Name","Frequency-Raw","Frequency-HPO","Additional Info from GD Source","GD Source","Disease ID for link","Gene Ensemble ID")
