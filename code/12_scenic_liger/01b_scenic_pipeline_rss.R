#!/usr/bin/env Rscript
# scenic_rss_jenson_shannon_divergence_entropy_par.R
# Do in parallel
# Run through shannon entropy alorithm from cell paper

# Requirements
library(AUCell)
library(SCENIC)
library(SingleCellExperiment)
library(foreach)
library(doParallel)
library(argparse)


### Options/Args
args <- commandArgs(trailingOnly=TRUE)
in_region <- as.character(args[1])
in_celltype <- as.character(args[2])
scenic_dir <- as.character(args[3])
num_cores <- detectCores() / 2
registerDoParallel(num_cores)
#########################################################################################
### Functions
#########################################################################################
h_func <- function(p){
  # calculate the JSD
  p <- p*log(p)
  p[is.na(p)] <- 0
  p <- -sum(p)
  return(p)
}

#########################################################################################
### Procedure
#########################################################################################
setwd(scenic_dir)
regulon_auc <- readRDS("int/3.4_regulonAUC.Rds")
scenicOptions <- readRDS("int/scenicOptions.Rds")
cellInfo <- readRDS("int/cellInfo.Rds")

# formating 
auc_scores <- getAUC(regulon_auc)
cellInfo$cluster <- factor(cellInfo$cluster)
levels <- levels(factor(cellInfo$cluster))
# remove repeated extended regulons
auc_scores <- auc_scores[onlyNonDuplicatedExtended(rownames(auc_scores)), ]
# making normalized auc matrix 
normalized_regulons <- auc_scores/rowSums(auc_scores)

## Proc
# making an "identity" matrix and then normalizing it
identity <- c()
for(i in 1:length(levels(cellInfo$cluster))){
  identity <- rbind(identity, as.character(cellInfo$cluster))
}

for(i in 1:nrow(identity)){
  for(n in 1:ncol(identity)){
    if(identity[i,n] == levels[i]){
      identity[i,n] = 1
    }
    else{
      identity[i,n] = 0
    }
  }
}

identity <- as.data.frame(identity)
identity <- data.frame(sapply(identity, function(x) as.numeric(as.character(x))))
normalized_identity <- identity/rowSums(identity)
rownames(normalized_identity) <- levels
colnames(normalized_identity) <- colnames(auc_scores)

# Parallel version
returned_matrix <- c()
returned_matrix <- foreach(i = 1:length(levels), .combine = "cbind") %:%
  foreach(n = 1:nrow(normalized_regulons), .combine = "c") %dopar% { 
    h_rc <- h_func((normalized_identity[i,] + normalized_regulons[n,])/2)
    h_r <- h_func(normalized_regulons[n,])
    h_c <- h_func(normalized_identity[i,])
    rss <- 1-sqrt(h_rc - ((h_c + h_r)/2))
    return(rss)
  }


colnames(returned_matrix) <- levels
rownames(returned_matrix) <- rownames(normalized_regulons)

saveRDS(
  returned_matrix,
  file = paste0("int/regulon_specificity_score_",in_celltype,"_", in_region,".rds")
)

celltypes <- colnames(returned_matrix)
for (i in 1:length(celltypes)) {
  ct <- celltypes[i]
  returned_matrix <- as.data.frame(cbind(returned_matrix, rank(-returned_matrix[,ct])))
  colnames(returned_matrix)[length(celltypes) + i] <- paste0(ct,"_rank")
}

write.csv(
  returned_matrix, 
  file =paste0("int/regulon_specificity_score_",in_celltype,"_", in_region,".csv")
)
