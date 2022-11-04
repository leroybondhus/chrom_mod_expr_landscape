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
library(stringr)

file_expansion_factor <- 12 ## from SRA documentation, fasterq-dump requires up to 8x space available
dirs <- list(data="./data/",
             sra="./data/sra/",
             sra_metadata="./data/sra/metadata/",
             scratch="/u/scratch/l/leroybon/")
date <- format(Sys.time(),"%Y%m%d")
files <- list(sra_run_log = paste0(dirs$sra,date,"_sra.log"))
for(d in dirs){
  if(!dir.exists(d)){dir.create(d)}
}
file.remove(files$sra_run_log)
for(f in files){
  system(paste0("touch ",f))
}
to_log <- paste("time", "srr_id","status", sep="\t")
write(to_log,files$sra_run_log, sep="\t", append = T)

##### 
if(file.exists("./data/sra_metadata/SraExpPack_SRP018525.xml")){
  root_sras <- xmlRoot(xmlParse("./data/sra_metadata/SraExpPack_SRP018525.xml"))
  # temp <- xmlToList(root_sras[[1]])
  sras <- xmlToList("./data/sra_metadata/SraExpPack_SRP018525.xml")
  temp_cols <- c("srp_id","srr_id","library_info","barcode_file","taxon_id","sample_title", "std_dev_time","source","file_size")
  sra_metadata <- data.frame(matrix(nrow=length(sras),ncol=length(temp_cols))) 
  colnames(sra_metadata) <- temp_cols
  for(i in 1:length(sras)){
    temp_df <- cbind(xmlToDataFrame(xmlElementsByTagName(root_sras[[i]],"TAG", recursive=TRUE)),
                     xmlToDataFrame( xmlElementsByTagName(root_sras[[i]],"VALUE", recursive=TRUE)))
    colnames(temp_df) <- c("TAG","VALUE")
    
    sra_metadata$srp_id[i] <- sras[[i]]$EXPERIMENT$STUDY_REF$IDENTIFIERS$PRIMARY_ID
    sra_metadata$srr_id[i] <- sras[[i]]$RUN_SET$RUN$IDENTIFIERS$PRIMARY_ID
    sra_metadata$taxon_id[i] <- sras[[i]]$SAMPLE$SAMPLE_NAME$TAXON_ID
    sra_metadata$file_size[i] <- as.numeric(sras[[i]]$RUN_SET$RUN$.attrs[["size"]]) ## file size
    sra_metadata$sample_title[i] <- sras[[i]]$SAMPLE$TITLE
    sra_metadata$source[i] <- temp_df$VALUE[which(temp_df$TAG == "source_name")]
    
    sra_metadata$library_info[i] <- paste0(c(as.character(sras[[1]]$EXPERIMENT$DESIGN$LIBRARY_DESCRIPTOR[c(1,3)]),
                                        names(sras[[1]]$EXPERIMENT$DESIGN$LIBRARY_DESCRIPTOR$LIBRARY_LAYOUT)),
                                      collapse = "__")
  }
  utils:::format.object_size(sum(as.numeric(sra_metadata$file_size) ), "Gb")
}

sra_metadata <- sra_metadata[which(sra_metadata$taxon_id == "10090"),]

## fasterq-dump -> fastq -> bam -> featureCount

## check how much space is available on scratch dir for processing

sra_metadata$status <- factor("-",levels=c("-","running","success","fail"))
space_total <- system(paste0("df ", dirs$scratch),intern = T)
space_total <- as.numeric(strsplit(space_total[2], " +")[[1]][4]) * 1024
space_total <- space_total * 0.8 ## leave space for margin of error
print(paste0(utils:::format.object_size(space_total, "Gb"),
             " available at '", dirs$scratch, "' "))

space_total

not_yet_run <- "-"
## write a log file to keep track of file status
while(any(sra_metadata$status==not_yet_run )){
  
  which_not_run <- which(sra_metadata$status==not_yet_run)
  which_running <- which(sra_metadata$status=="running")
  space_unavail <- sum(as.numeric(sra_metadata$file_size[which_running])) * file_expansion_factor
  space_avail <- space_total - space_unavail
  
  which_next <- which_not_run[1]
  if(space_avail >= sra_metadata$file_size[which_next]){
    sra_metadata$status[which_next] <- "running" ## success and fail status should be assigned at end of script
    
    ## UPDATE SYSTEM LOG FILE
    to_log <- paste(Sys.time(),sra_metadata$srr_id[which_next],as.character(sra_metadata$status[which_next]), sep="\t")
    write(to_log,files$sra_run_log, sep="\t", append = T)
    ## UPDATE SCRIPT LOG (CURRENTLY aka sra_metadata)
    current_log <- fread(files$sra_run_log)
    which_running <- which(is.element(sra_metadata$srr_id, current_log$srr_id[which(current_log$status=="running")]))
    which_success <- which(is.element(sra_metadata$srr_id, current_log$srr_id[which(current_log$status=="success")]))
    which_fail <- which(is.element(sra_metadata$srr_id, current_log$srr_id[which(current_log$status=="fail")]))
    which_running <- setdiff(which_running, union(which_success,which_fail))
    sra_metadata
    
    SUBMIT SCRIPT TO RUN sra_metadata$srr_id[which_next]
  } else {
    READ LOG FILE
    IF ALL JOBS COMPLETED AND NO SPACE EXIT
  }
  
}

# 
# while(jobs_remain){
#   check space_avail from space_unavail from jobs table object
#   if(no_space_avail){
#     if(no_jobs_running){stop "NO SPACE AVAIL - CLEAR SCRATCH"}
#     sleep 1 min, then check again
#   } else{
#     queu_next_job
#     update jobs table object
#     SUBMIT NEXT JOB TO QUEUE3
#     
#   }
  

# run processes and store pids in array

##### https://stackoverflow.com/questions/356100/how-to-wait-in-bash-for-several-subprocesses-to-finish-and-return-exit-code-0
# for i in $n_procs; do
# ./procs[${i}] &
#   pids[${i}]=$!
#   done
# 
# # wait for all pids
# for pid in ${pids[*]}; do
# wait $pid
# done

#### 


# -> loop through 
# -> grab particular chunk size 
# -> download to scratch 
# -> process from fastq->bam->featureCount
# repeat until done print available space and num remaining as goes

