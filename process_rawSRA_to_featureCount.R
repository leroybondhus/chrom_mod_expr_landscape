### read in xml of data
library(XML)

## make loop IF way of accessing relevant data pieces is standard ELSE keep SRP by SRP approach
# xml_files <- list.files("./data/sra_metadata",full.names = T)
# xml_files <- xml_files[grep("\\.xml",xml_files)]
# for(f in xml_files){
#   temp <- xmlToList(f)
# }

if(file.exists("./data/sra_metadata/SraExpPack_SRP018525.xml")){
  root_sras <- xmlRoot(xmlParse("./data/sra_metadata/SraExpPack_SRP018525.xml"))
  # temp <- xmlToList(root_sras[[1]])
  sras <- xmlToList("./data/sra_metadata/SraExpPack_SRP018525.xml")
  temp_cols <- c("srp_id","srr_id","taxon_id","file_size","sample_title", "sample_age","source")
  sra_metadata <- data.frame(matrix(nrow=length(sras),ncol=length(temp_cols))) 
  colnames(sra_metadata) <- temp_cols
  for(i in 1:length(sras)){
    sra_metadata$srp_id[i] <- sras[[i]]$EXPERIMENT$STUDY_REF$IDENTIFIERS$PRIMARY_ID
    sra_metadata$srr_id[i] <- sras[[i]]$RUN_SET$RUN$IDENTIFIERS$PRIMARY_ID
    sra_metadata$taxon_id[i] <- sras[[i]]$SAMPLE$SAMPLE_NAME$TAXON_ID
    sra_metadata$file_size[i] <- sras[[i]]$RUN_SET$RUN$.attrs[["size"]] ## file size
    # utils:::format.object_size(as.numeric( sras[[i]]$RUN_SET$RUN$.attrs["size"]),"Mb")
    sra_metadata$sample_title[i] <- sras[[i]]$SAMPLE$TITLE
    
    temp_df <- cbind(xmlToDataFrame(xmlElementsByTagName(root_sras[[i]],"TAG", recursive=TRUE)),
                     xmlToDataFrame( xmlElementsByTagName(root_sras[[i]],"VALUE", recursive=TRUE)))
    colnames(temp_df) <- c("TAG","VALUE")
    sra_metadata$source[i] <- temp_df$VALUE[which(temp_df$TAG == "source_name")]
  }
  utils:::format.object_size(sum(as.numeric(sra_metadata$file_size) ), "Gb")
}



# -> loop through 
# -> grab particular chunk size 
# -> download to scratch 
# -> process from fastq->bam->featureCount
# repeat until done print available space and num remaining as goes

