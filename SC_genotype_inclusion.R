library(ggplot2)
library(ggpubr)
library(biomaRt)
library(stringr)
library(dendextend)
library(dplyr)
library(Hmisc)
library(doParallel)
library(foreach)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(reldist)
library(gridExtra)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(topGO) 
library(enrichplot)
library(data.table)

figures_dir <- paste(getwd(),"figures_2", sep = "/")
date <- format(Sys.time(), format="%Y%m%d") 
registerDoParallel(detectCores()-1)


## load gtex here

## import gtex medians data

temp <- tempfile()
download.file("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",temp)
gtex <- read.table( temp, skip=2, header = TRUE, sep = "\t")
gtex_rowinfo <- data.frame(Name=gtex$Name, Description=gtex$Description)
rownames(gtex) <- gtex$Name
gtex <- as.matrix(gtex[,3:ncol(gtex)])
unlink(temp)
rm(temp)

## import ensembl gene data

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genes <- getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol', 'ensembl_gene_id','gene_biotype'),
               filters = list('biotype'='protein_coding'),
               mart = ensembl, useCache = F) 
genes <- genes[which(is.element(genes$chromosome_name, c(1:22, "X", "Y", "MT")) & genes$hgnc_symbol != "" ) ,]


## clean gtex dataset here
## show mitochondrial genes drive a large part of sample similarity
## Note sum of medians in gtex not quite 1M - likely artifact of taking medians

temp <- data.frame(names=colnames(gtex), total_transcripts=colSums(gtex))
ggplot( temp ,aes(y=total_transcripts, x=names))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
## match genes between gtex and ensembl
gtex_names <- str_split_fixed(gtex_rowinfo$Name, "[.]", 2)[,1]
which <- which(genes$chromosome_name != "MT" )
gtex_cleaned <- gtex[which(is.element(gtex_names, genes$ensembl_gene_id[which])),]
which <- which(genes$chromosome_name == "MT" )
gtex_cleanedMT <- gtex[which(is.element(gtex_names, genes$ensembl_gene_id[which])),]
##non-mitochondrial TPM sum
temp <- data.frame(names=colnames(gtex_cleaned), total_transcripts=colSums(gtex_cleaned))
ggplot( temp ,aes(y=total_transcripts, x=names))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
##mitochondrial TPM sum
temp <- data.frame(names=colnames(gtex_cleanedMT), total_transcripts=colSums(gtex_cleanedMT))
ggplot( temp ,aes(y=total_transcripts, x=names))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
rm(gtex_cleanedMT)
##renormalize TPM without mitochondrial genes
for(i in 1:ncol(gtex_cleaned)){
  gtex_cleaned[,i] <- (gtex_cleaned[,i]*1e6 / sum(gtex_cleaned[,i]))
  ## set all very low counts to zero to avoid variable read depth issues
  gtex_cleaned[,i] [which(gtex_cleaned[,i]  < 1)] <- 0
  gtex_cleaned[,i] <- (gtex_cleaned[,i]*1e6 / sum(gtex_cleaned[,i]))
}

temp <- data.frame(names=colnames(gtex_cleaned), total_transcripts=colSums(gtex_cleaned))
ggplot( temp ,aes(y=total_transcripts, x=names))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gtex <- gtex_cleaned; rm(gtex_cleaned)
exp_mat <- gtex
## log10(TPM+1) transform of data
exp_mat <- log10(exp_mat+1)
## remove mitochondrial contribution
## median normalization
temp <- exp_mat
temp[which(temp==0 | is.infinite(temp))] <- NA
boxplot(temp)
median_normalize <- TRUE
if(median_normalize){
  for(i in 1:ncol(exp_mat)){
    exp_mat[,i] <- exp_mat[,i] / median(exp_mat[,i][which(exp_mat[,i] > 0)])
  }
}
temp <- exp_mat
temp[which(temp==0)] <- NA
boxplot(temp)

calc_zscore_matrix<- function(dat) {
  zscores <- dat; zscores[] <- 0 
  means <- rowMeans(dat)
  sds <- as.numeric(rep(NA,length(means)))
  which <- which(means != 0)
  sds[which] <- apply(dat[which,],1,sd)
  for(j in 1:ncol(dat)){zscores[,j] <- (dat[,j] - means)/sds  }
  return(zscores)
}

calc_dot_product_similarity_matrix <- function(dat) {
  dot_product_similarity_matrix <- matrix(0, nrow = ncol(dat), ncol = ncol(dat))
  colnames(dot_product_similarity_matrix) <- colnames(dat)
  rownames(dot_product_similarity_matrix) <- colnames(dat)
  for(i in 1:ncol(dat)){
    for(j in 1:ncol(dat)){
      which_i <- which(!is.na(dat[,i])) ## ignore NAs
      which_j <- which(!is.na(dat[,j])) ## ignore NAs
      dot_product_similarity_matrix[i,j] <- sum(dat[which_i,i] * dat[which_j,j]) / (norm(dat[which_i,i],"2")*norm(dat[which_j,j],"2"))
    }
  }
  return(dot_product_similarity_matrix)
}

### uses Equation 1. from paper 
add_dist_to_parent <- function(dend, dist_to_parent=0){
  ## note: distance to parent is fed in at the start of the function
  attributes(dend) <- c(attributes(dend), dist_to_parent=dist_to_parent)
  ## test if at leaf node
  if(!is.null(attributes(dend)$leaf) && attributes(dend)$leaf){
    return(dend)
  }
  for(i in 1:length(dend)){ ## length of dend should be number of child nodes
    ## distance to parent is simply the difference in height between parent and child
    dist_to_parent <- attributes(dend)$height - attributes(dend[[i]])$height 
    dend[[i]] <- add_dist_to_parent(dend[[i]], 
                                    dist_to_parent = dist_to_parent)
  }
  return(dend)
}

## this functions calculates and adds weights to dendrogram object using the 'dist_to_parent' attribute added previously
## weight_of_parent parameter exists only for recursion and should not be manually adjusted without understanding it's function
add_weights <- function(dend, weight_of_parent=0){
  weight <- (attributes(dend)$dist_to_parent / attributes(dend)$members) + weight_of_parent 
  attributes(dend) <- c(attributes(dend), weight=weight)
  ## test if at leaf node
  if(!is.null(attributes(dend)$leaf) && attributes(dend)$leaf){
    return(dend)
  }
  for(i in 1:length(dend)){ ## length of dend should be number of child nodes
    dend[[i]] <- add_weights(dend[[i]], weight_of_parent=weight)
  }
  return(dend)
}

## this function returns the weights from a dendrogram object that has a "weight" attribute at leaves. Also requires the order of the vector to return based on names of leaves
get_weights <- function(dend, name_order){
  weights <- setNames(get_leaves_attr(dend,"weight"),nm=get_leaves_attr(dend,"lab") )
  weights <- weights[order(factor(names(weights),levels = name_order))]
  return(weights)
}

# function to calculate weighted zscores given matrix and vector of weights. column names of the matrix and names of the weight vector must match
calc_weighted_zscore_matrix <- function(mat, weights){
  if(any( colnames(mat) != names(weights) )){stop("WARNING: mismatch in weights names and matrix colnames order")}
  weighted_mat <- mat; weighted_mat[] <- 0
  for (i in 1:length(weights)){
    weighted_mat[,i] <- weights[i]*mat[,i]
  }
  weighted_means <- numeric(length = nrow(weighted_mat))
  sum_of_weights <- sum(weights)
  for (i in 1:nrow(weighted_mat)){
    weighted_means[i] <- sum(weighted_mat[i,]) / sum_of_weights
  }
  weighted_var <- numeric(length=nrow(mat))
  for (i in 1:nrow(mat)){
    weighted_var[i] <- Hmisc::wtd.var(mat[i,],weights=weights)
  }
  weighted_sd <- sqrt(weighted_var)
  for(i in 1:ncol(mat)){
    mat[,i] <- (mat[,i]-weighted_means)/weighted_sd
  }
  weighted_zscores <- mat
  return(weighted_zscores)
}
# weighted tau
calc_weighted_tau <- function(te_matrix, weights_vector){
  xhat_matrix <- matrix(nrow=nrow(te_matrix),ncol=ncol(te_matrix))
  te_row_maxima <- apply(te_matrix, 1, max)
  for(j in 1:ncol(te_matrix)){
    xhat_matrix[,j] <- te_matrix[,j] / te_row_maxima
  }
  temp_matrix <- matrix(nrow=nrow(te_matrix),ncol=ncol(te_matrix))
  for (i in 1:nrow(te_matrix)){
    temp_matrix[i,] <- weights_vector - (xhat_matrix[i,] * weights_vector)
  }
  tau <- numeric(length = nrow(temp_matrix))
  for (i in 1:nrow(temp_matrix)){
    temp <- sum(temp_matrix[i,]) / (sum(weights_vector) - weights_vector[which.max(temp_matrix[i,])])
    tau[i] <- ifelse(length(temp)==0,NA,temp)
  }
  
  ## add normalization (believe this is a numeric instability issue from dividing small numbers)
  # tau <- tau / max(tau, na.rm=T)
  ## alternative, set all > 1 to 1 (when looking at plots for different cutoffs, normalizing true 1 values causes issue)
  tau[which(tau > 1)] <- 1
  return(tau)
}

calc_weighted_tsi <- function(te_matrix,weights_vector){
  weighted_matrix <- t(apply(te_matrix,1, "*", weights_vector))
  tsi <- numeric(length=nrow(weighted_matrix))
  for(i in 1:nrow(weighted_matrix)){
    tsi[i] <- weights_vector[which.max(te_matrix[i,])] * max(te_matrix[i,]) / sum(weighted_matrix[i,]) 
  }
  names(tsi) <- rownames(te_matrix)
  return(tsi)
}

calc_weighted_gini <- function(te_matrix, weights_vector){
  gini_values <-  c()
  for (i in 1:nrow(te_matrix)){
    temp <- as.numeric(te_matrix[i,])
    temp <- reldist::gini(temp, weights_vector)
    gini_values <- append(gini_values,temp)
  }
  return(gini_values)
}



### generalizing a list to store results in - this will make it easier to extend later if necessary.
## Note: only need the weighted version of each equation as each simplifies to flat version when all weights are 1
specificity_measures <- list(func_names=c("Zscore", "Tau", "Tsi","Gini"),
                             funcs=list(Zscore=calc_weighted_zscore_matrix,
                                        Tau=calc_weighted_tau,
                                        Tsi=calc_weighted_tsi,
                                        Gini=calc_weighted_gini),
                             out_type=list(Zscore="matrix",
                                           Tau="vector",
                                           Tsi="vector",
                                           Gini="vector")
)

## only 1 similarity function tested for now, can make as list later
similarity_func <- function(exp_mat){calc_dot_product_similarity_matrix(calc_zscore_matrix(exp_mat))}
## only 1 clustering fucntion tested for now, can make as a list later
cluster_func <- function(sim_mat){add_weights(add_dist_to_parent(as.dendrogram(hclust(as.dist(1-sim_mat), method = "average") ) ))} 


## select which to be all non-brain and one brain 
which <- c(grep("Brain", colnames(exp_mat), invert = TRUE ),sample(grep("Brain", colnames(exp_mat)),1))
flat <- rep(1,ncol(exp_mat)); names(flat) <- colnames(exp_mat)
spec_one_all_list <- list()

for(measure in specificity_measures$func_names){
  specificity_func <- specificity_measures$funcs[[measure]]
  if(!is.function(specificity_func)){print(paste(measure, "func is not a function")); next;}
  spec_one_all_list[[paste(measure, "flat_one",sep="_")]] <- specificity_func(exp_mat[,which], flat[which])
  spec_one_all_list[[paste(measure, "flat_all",sep="_")]] <- specificity_func(exp_mat, flat)
  if(measure=="Zscore"){
    spec_one_all_list[[paste(measure, "flat_all",sep="_")]] <- spec_one_all_list[[paste(measure, "flat_all",sep="_")]][,which]
  }
  dot_sim <- similarity_func(exp_mat[,which])
  sim_tree <- cluster_func(dot_sim)
  weights <- get_weights(sim_tree, colnames(exp_mat)[which])
  spec_one_all_list[[paste(measure, "weighted_one",sep="_")]] <- specificity_func(exp_mat[,which], weights)
  dot_sim <- similarity_func(exp_mat)
  sim_tree <- cluster_func(dot_sim)
  weights <- get_weights(sim_tree, colnames(exp_mat))
  spec_one_all_list[[paste(measure, "weighted_all",sep="_")]] <- specificity_func(exp_mat, weights)
  if(measure=="Zscore"){
    spec_one_all_list[[paste(measure, "weighted_all",sep="_")]] <- spec_one_all_list[[paste(measure, "weighted_all",sep="_")]][,which]
  }
}

grep("ENSG00000268173", as.character(exp_mat[,33]))

###gene expression = tau measure
## KAT6A --> less specific, low tau value (spec_one_all_list)
## more specific = high tau value (TBX5)

### tbx5 disease associated is heart related, kat6a --> different HPO terms
##what we need to do is develop a measure for a given gene disease, how specific is the phenotype to a given system
###if you have something like KAT6A (20 different phenotypes and 20 different syndromes)
### relatively high information content --> not occur in different diseases
##KAT6A not focused on a specific system, not a specific system that is affected
###TBX5 --> relatively high specificity of systems affected for TBX5
### gene expression specificity (Tau) vs phenotype specificity for a given syndrome, is the phenotype relatively specific to particular system
#for a given disease, is the phenotype relatively specific to a given system


### 