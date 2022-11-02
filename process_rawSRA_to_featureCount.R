## author: Leroy Bondhus
## date: 2022 11 2

#### MANUAL PREP STEP ####
### SraExpPack_{SRPID}.xml files downloaded manually
### Search for SRP associated with relevant GEO or SRA in NCBI SRA
### Click 'Sent to:' button on top right of result page
### Select 'Full XML' and click 'Create File' button
### Move downloaded 'SraExperimentPackage.xml' to './data/sra_metadata/SraExpPack_{SRPID}.xml'

### read in xml of data
library(XML)

dirs <- list(data="./data/",
             sra_metadata="./data/sra/metadata/",
             scratch="$SCRATCH/")
 

##### 
if(file.exists("./data/sra_metadata/SraExpPack_SRP018525.xml")){
  root_sras <- xmlRoot(xmlParse("./data/sra_metadata/SraExpPack_SRP018525.xml"))
  # temp <- xmlToList(root_sras[[1]])
  sras <- xmlToList("./data/sra_metadata/SraExpPack_SRP018525.xml")
  temp_cols <- c("srp_id","srr_id","taxon_id","file_size","sample_title", "std_dev_time","source","library_info")
  sra_metadata <- data.frame(matrix(nrow=length(sras),ncol=length(temp_cols))) 
  colnames(sra_metadata) <- temp_cols
  for(i in 1:length(sras)){
    temp_df <- cbind(xmlToDataFrame(xmlElementsByTagName(root_sras[[i]],"TAG", recursive=TRUE)),
                     xmlToDataFrame( xmlElementsByTagName(root_sras[[i]],"VALUE", recursive=TRUE)))
    colnames(temp_df) <- c("TAG","VALUE")
    
    sra_metadata$srp_id[i] <- sras[[i]]$EXPERIMENT$STUDY_REF$IDENTIFIERS$PRIMARY_ID
    sra_metadata$srr_id[i] <- sras[[i]]$RUN_SET$RUN$IDENTIFIERS$PRIMARY_ID
    sra_metadata$taxon_id[i] <- sras[[i]]$SAMPLE$SAMPLE_NAME$TAXON_ID
    sra_metadata$file_size[i] <- sras[[i]]$RUN_SET$RUN$.attrs[["size"]] ## file size
    sra_metadata$sample_title[i] <- sras[[i]]$SAMPLE$TITLE
    sra_metadata$source[i] <- temp_df$VALUE[which(temp_df$TAG == "source_name")]
    
    sra_metadata$library_info[i] <- paste0(c(as.character(sras[[1]]$EXPERIMENT$DESIGN$LIBRARY_DESCRIPTOR[c(1,3)]),
                                        names(sras[[1]]$EXPERIMENT$DESIGN$LIBRARY_DESCRIPTOR$LIBRARY_LAYOUT)),
                                      collapse = "__")
  }
  utils:::format.object_size(sum(as.numeric(sra_metadata$file_size) ), "Gb")
}


## fasterq-dump -> fastq -> bam -> featureCount


# -> loop through 
# -> grab particular chunk size 
# -> download to scratch 
# -> process from fastq->bam->featureCount
# repeat until done print available space and num remaining as goes

