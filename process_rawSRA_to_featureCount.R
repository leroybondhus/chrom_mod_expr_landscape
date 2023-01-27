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
library(data.table)

file_expansion_factor <- 12 ## from SRA documentation, fasterq-dump requires up to 8x space available
dirs <- list(data="./data/",
             sra="/u/scratch/l/leroybon/chrom_mod_proj_data/sra/",
             sra_metadata="./data/sra_metadata/",
             scratch="/u/scratch/l/leroybon/")
date <- format(Sys.time(),"%Y%m%d")
## Note may want to read in old log if planning to only run failed runs
## NOTE: LEROY code this up to read in previous run log
files <- list(sra_run_log = paste0(dirs$sra,date,"_sra.log"))
for(d in dirs){
  if(!dir.exists(d)){dir.create(d)}
}
if(file.exists(files$sra_run_log)){
  file.remove(files$sra_run_log) 
}
for(f in files){
  system(paste0("touch ",f))
}
to_log <- paste("time", "srr_id","status", sep="\t")
write(to_log,files$sra_run_log, sep="\t", append = T)



temp_cols <- c("srp_id","srr_id",
               "library_info","barcode_file",
               "taxon_id","sample_title",
               "std_dev_time","source",
               "file_size","file_space_req")
sra_metadata_full <- data.frame(matrix(nrow=0,ncol=length(temp_cols)))
colnames(sra_metadata_full) <- temp_cols

#####
files$srps <- list(SRP018525 = list( xml = paste0(dirs$sra_metadata,"SraExpPack_SRP018525.xml")),
                   SRP110669 = list( xml = paste0(dirs$sra_metadata,"SraExpPack_SRP110669.xml")),
                   SRP190004 = list( xml = paste0(dirs$sra_metadata,"SraExpPack_SRP190004.xml")),
                   SRP161714 = list( xml = paste0(dirs$sra_metadata,"SraExpPack_SRP161714.xml")),
                   SRP126776 = list( xml = paste0(dirs$sra_metadata,"SraExpPack_SRP126776.xml"))
                   )
for(srp in files$srps ){
  if(!file.exists(srp$xml)){ stop(paste0(srp$xml, ": file does not exist")) }
  root_sras <- xmlRoot(xmlParse(srp$xml))
  # temp <- xmlToList(root_sras[[1]])
  sras <- xmlToList(srp$xml)
  sra_metadata <- data.frame(matrix(nrow=length(sras),ncol=ncol(sra_metadata_full))) 
  colnames(sra_metadata) <- colnames(sra_metadata_full)
  for(i in 1:length(sras)){
    temp_df <- cbind(xmlToDataFrame(xmlElementsByTagName(root_sras[[i]],"TAG", recursive=TRUE)),
                     xmlToDataFrame( xmlElementsByTagName(root_sras[[i]],"VALUE", recursive=TRUE)))
    colnames(temp_df) <- c("TAG","VALUE")
    
    sra_metadata$srp_id[i] <- sras[[i]]$EXPERIMENT$STUDY_REF$IDENTIFIERS$PRIMARY_ID
    sra_metadata$srr_id[i] <- sras[[i]]$RUN_SET$RUN$IDENTIFIERS$PRIMARY_ID
    sra_metadata$taxon_id[i] <- sras[[i]]$SAMPLE$SAMPLE_NAME$TAXON_ID
    sra_metadata$file_size[i] <- as.numeric(sras[[i]]$RUN_SET$RUN$.attrs[["size"]]) ## file size
    sra_metadata$file_space_req[i] <- as.numeric(sras[[i]]$RUN_SET$RUN$.attrs[["size"]])*file_expansion_factor
    sra_metadata$sample_title[i] <- sras[[i]]$SAMPLE$TITLE
    sra_metadata$source[i] <- paste0(temp_df$VALUE,collapse="---") #temp_df$VALUE[which(temp_df$TAG == "source_name")]
    
    sra_metadata$library_info[i] <- paste0(c(as.character(sras[[1]]$EXPERIMENT$DESIGN$LIBRARY_DESCRIPTOR[c(1,3)]),
                                        names(sras[[1]]$EXPERIMENT$DESIGN$LIBRARY_DESCRIPTOR$LIBRARY_LAYOUT)),
                                      collapse = "__")
  }
  sra_metadata_full <- rbind(sra_metadata_full, sra_metadata)
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
  space_unavail <- sum(sra_metadata$file_space_req[which_running])
  space_avail <- space_total - space_unavail
  
  which_next <- which_not_run[1]
  if(space_avail >= sra_metadata$file_space_req[which_next] ){
    space_avail <- space_avail - sra_metadata$file_space_req[which_next] 
    ### SUBMIT SCRIPT TO RUN
    # SCRIPT SUBMISSION sra_metadata$srr_id[which_next]
    ### 
    sra_metadata$status[which_next] <- "running" ## success and fail status must be assigned at end of external SCRIPT
    ## UPDATE SYSTEM LOG FILE
    to_log <- paste(Sys.time(),sra_metadata$srr_id[which_next],as.character(sra_metadata$status[which_next]), sep="\t")
    write(to_log,files$sra_run_log, sep="\t", append = T)
    ## UPDATE SESSION LOG (CURRENTLY aka sra_metadata)
    current_log <- fread(files$sra_run_log)
    which_running <- which(is.element(sra_metadata$srr_id, current_log$srr_id[which(current_log$status=="running")]))
    which_success <- which(is.element(sra_metadata$srr_id, current_log$srr_id[which(current_log$status=="success")]))
    which_fail <- which(is.element(sra_metadata$srr_id, current_log$srr_id[which(current_log$status=="fail")]))
    which_running <- setdiff(which_running, union(which_success,which_fail))
    sra_metadata$status[which_running] <- "running" 
    sra_metadata$status[which_success] <- "success"
    sra_metadata$status[which_fail] <- "fail"
  } else {
    current_log <- fread(files$sra_run_log)
    which_running <- which(is.element(sra_metadata$srr_id, current_log$srr_id[which(current_log$status=="running")]))
    which_success <- which(is.element(sra_metadata$srr_id, current_log$srr_id[which(current_log$status=="success")]))
    which_fail <- which(is.element(sra_metadata$srr_id, current_log$srr_id[which(current_log$status=="fail")]))
    which_running <- setdiff(which_running, union(which_success,which_fail))
    if(length(which_running)==0){
      stop("Jobs remain but no space avail and no jobs running: Need to free up space to continue")
    } else{
      ## no space but jobs running which will eventually complete and free up space.
      ## Let's just wait a little while then check if anything has completed ok?
      system("sleep 2.5m")
      system("date +'%Y-%m-%d %H:%M:%S PDT' ",intern = T) ## unix to R format from Sys.time()
    }
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

