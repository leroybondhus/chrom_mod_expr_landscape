library(ggplot2)
library(ggpubr)
library(biomaRt)
library(stringr)
library(dendextend)
library(dplyr)
library(Hmisc)
library(doParallel)
library(foreach)
registerDoParallel(detectCores())
library(ComplexHeatmap)
library(circlize)

temp <- tempfile()
download.file("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",temp)
gtex <- read.table( temp, skip=2, header = TRUE, sep = "\t")
gtex_rowinfo <- data.frame(Name=gtex$Name, Description=gtex$Description)
rownames(gtex) <- gtex$Name
gtex <- as.matrix(gtex[,3:ncol(gtex)])
unlink(temp); rm(temp)

## import ensembl gene data
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genes <- getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol', 'ensembl_gene_id','gene_biotype'),
               filters = list('biotype'='protein_coding'),
               mart = ensembl, useCache = F) 
genes <- genes[which(is.element(genes$chromosome_name, c(1:22, "X", "Y", "MT")) & genes$hgnc_symbol != "" ) ,]


## show mitochondrial genes drive a large part of sample similarity
## Note sum of medians in gtex not quite 1M - likely artifact of taking medians
barplot(colSums(gtex))
## match genes between gtex and ensembl
gtex_names <- str_split_fixed(gtex_rowinfo$Name, "[.]", 2)[,1]
which <- which(genes$chromosome_name != "MT" )
gtex_cleaned <- gtex[which(is.element(gtex_names, genes$ensembl_gene_id[which])),]
which <- which(genes$chromosome_name == "MT" )
gtex_cleanedMT <- gtex[which(is.element(gtex_names, genes$ensembl_gene_id[which])),]
barplot(colSums(gtex_cleaned))     ##non-mitochondrial TPM sum
barplot(colSums(gtex_cleanedMT))   ##mitochondrial TPM sum
rm(gtex_cleanedMT)
##renormalize TPM without mitochondrial genes
for(i in 1:ncol(gtex_cleaned)){
  gtex_cleaned[,i] <- (gtex_cleaned[,i]*1e6 / sum(gtex_cleaned[,i]))
}
barplot(colSums(gtex_cleaned))
gtex <- gtex_cleaned; rm(gtex_cleaned)
exp_mat <- gtex
## log10(TPM+1) transform of data
exp_mat <- log10(exp_mat+1)
## remove mitochondrial contribution


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
### uses Equation 1. 
## modelled on Yenifer's code
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
    weighted_var[i] <- wtd.var(mat[i,],weights=weights)
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
  tau <- c()
  den <- sum(weights_vector) - 1
  for (i in 1:nrow(temp_matrix)){
    temp <- sum(temp_matrix[i,]) / den
    tau <- append(tau,temp)
  }
  ## add normalization (believe this is a numeric instability issue from dividing small numbers)
  tau <- tau / max(tau, na.rm=T)
  return(tau)
}


### generalizing a list to store results in - this will make it easier to extend later if necessary.
## Note: only need the weighted version of each equation as each simplifies to flat version when all weights are 1
specificity_measures <- list(func_names=c("Zscore", "Tau", "Tsi","Gini"),
                             funcs=list(Zscore=calc_weighted_zscore_matrix,
                                        Tau=calc_weighted_tau,
                                        Tsi=NA,
                                        Gini=NA),
                             out_type=list(Zscore="matrix",
                                           Tau="vector",
                                           Tsi="vector",
                                           Gini="vector")
)
## only 1 similarity function tested for now, can make as list later
similarity_func <- function(exp_mat){calc_dot_product_similarity_matrix(calc_zscore_matrix(exp_mat))}
## only 1 clustering fucntion tested for now, can make as a list later
cluster_func <- function(sim_mat){add_weights(add_dist_to_parent(as.dendrogram(hclust(as.dist(1-sim_mat), method = "single") ) ))}  


## dot similarity from initial Z scores
dot_sim <- similarity_func(exp_mat)
heatmap(dot_sim)
sim_tree <- cluster_func(dot_sim)
## plot similarity tree with weights
sim_tree %>% set("nodes_pch",19) %>% set("nodes_cex", 2.2*sqrt(get_nodes_attr(sim_tree,"weight"))) %>% plot
## take a look at the weights
weights <- get_weights(sim_tree, colnames(exp_mat))
ggplot( data.frame(names=names(weights), weights=weights),aes(x=names, y=weights))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
spec_mat_weighted <- calc_weighted_zscore_matrix(exp_mat, weights)   #### these are Z-scores as a measure of specificity
flat <- weights; flat[1:length(flat)] <- 1
spec_mat_flat  <- calc_weighted_zscore_matrix(exp_mat, flat)
## plot to look at global dif between flat and weighted



#### End: Code lifted from gene_specificity project: specificity_paper_main.Rmd ####

chromatin_modifier_table <- read.csv("chromatin_modifier_table_shortened.csv")

pattern <- "[a-zA-Z0-9_]+"

chromatin_modifier_vector <- str_extract(chromatin_modifier_table[,1], pattern)

index <- c()

#subset to find corresponding ENSEMBL ID 
for(i in 1:length(chromatin_modifier_vector)){
  temp <- which(str_extract(rownames(gtex), pattern) == chromatin_modifier_vector[i])
  index <- c(temp, index)
}

gtex_subset <- gtex[index,] #extract rows with the corresponding ENSEMBL ID

log_of_chromatin_modifier_vector <- log10(gtex_subset+1)

####### MODIFIED CHROMATIN MODIFIER VECTOR

## dot similarity from initial Z scores

dot_sim <- similarity_func(log_of_chromatin_modifier_vector)

dot_sim_matrix <- as.matrix(scale(dot_sim))

png("dot_sim.png", width = 1000, height = 1000)
Heatmap(dot_sim_matrix,show_row_dend = FALSE, show_column_dend = FALSE)
dev.off()

sim_tree <- cluster_func(dot_sim)

## plot similarity tree with weights

sim_tree %>% set("nodes_pch",19) %>% set("nodes_cex", 2.2*sqrt(get_nodes_attr(sim_tree,"weight"))) %>% plot

## take a look at the weights
weights <- get_weights(sim_tree, colnames(exp_mat))
ggplot( data.frame(names=names(weights), weights=weights),aes(x=names, y=weights))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
spec_mat_weighted <- calc_weighted_zscore_matrix(log_of_chromatin_modifier_vector, weights)   #### these are Z-scores as a measure of specificity
flat <- weights; flat[1:length(flat)] <- 1
spec_mat_flat  <- calc_weighted_zscore_matrix(log_of_chromatin_modifier_vector, flat)

#NaN produced
weighted_remove_na_n <- as.matrix(spec_mat_weighted[which(!is.nan(spec_mat_weighted[,1])),])

png("weighted_scaled.png", width = 1300, height = 1300)
Heatmap(weighted_remove_na_n, show_row_dend = FALSE, show_column_dend = FALSE)
dev.off()

hist(spec_mat_weighted, breaks = 100)

##Adjust the scale, which values are color legends, do not need dendograms
# unweighted - weighted
# push a file to github
# human phenotype ontology - what omen is? inheritance?

## Don't understand human phenotype ontology and difference between OMIM - online medelian inheritance in man
