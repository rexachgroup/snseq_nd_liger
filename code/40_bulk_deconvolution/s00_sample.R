source('./COMMON_DEFS.R')

base_dir <- getBaseDir()
input_dir <- getInputDir(base_dir)

out_path <- file.path(base_dir, 'Step01_PCA_OutlierRemoval')

load(file.path(input_dir, 'Salmon_count.matrix_Batch1_All.rda'))
load(file.path(input_dir, 'Salmon_count.matrix_Batch2_MergeAll.rda'))
load(file.path(input_dir, 'Salmon_count.matrix_Batch3_MergeAll.rda'))

logTPM_Batch1 <- as.data.frame(log2(count.matrix_Batch1_Merged+1))
logTPM_Batch2 <- as.data.frame(log2(count.matrix_Batch2_Merged+1))
logTPM_Batch3 <- as.data.frame(log2(count.matrix_Batch3_Merged+1))

rawExpr <- as.data.frame(cbind(logTPM_Batch1, logTPM_Batch2, logTPM_Batch3))

# (1) get whitelist of columns to include in metadata / sequencing --------
targets <- as.data.frame(rbind(targets_Batch1, targets_Batch2, targets_Batch3))
colnames(targets)
#  [1] "Autopsy.ID"                   "Source"                      
#  [3] "Clinical.Dx"                  "Primary.Neuropath.Dx"        
#  [5] "Age"                          "Sex"                         
#  [7] "PMI"                          "Batch"                       
#  [9] "SampleID"                     "Block.."                     
# [11] "Region"                       "Date_of_RNA_extraction"      
# [13] "RNA_conc"                     "Total_RNA_yield"             
# [15] "RIN"                          "FileName"                    
# [17] "PF_BASES"                     "PF_ALIGNED_BASES"            
# [19] "RIBOSOMAL_BASES"              "CODING_BASES"                
# [21] "UTR_BASES"                    "INTRONIC_BASES"              
# [23] "INTERGENIC_BASES"             "CORRECT_STRAND_READS"        
# [25] "INCORRECT_STRAND_READS"       "PCT_RIBOSOMAL_BASES"         
# [27] "PCT_CODING_BASES"             "PCT_UTR_BASES"               
# [29] "PCT_INTRONIC_BASES"           "PCT_INTERGENIC_BASES"        
# [31] "PCT_MRNA_BASES"               "PCT_USABLE_BASES"            
# [33] "PCT_CORRECT_STRAND_READS"     "MEDIAN_CV_COVERAGE"          
# [35] "MEDIAN_5PRIME_BIAS"           "MEDIAN_3PRIME_BIAS"          
# [37] "MEDIAN_5PRIME_TO_3PRIME_BIAS" "TOTAL_CLUSTERS"              
# [39] "ALIGNED_READS"                "AT_DROPOUT"                  
# [41] "GC_DROPOUT"                   "UNPAIRED_READS_EXAMINED"     
# [43] "UNMAPPED_READS"               "READ_PAIR_DUPLICATES"        
# [45] "READ_PAIR_OPTICAL_DUPLICATES" "PERCENT_DUPLICATION"         
# [47] "ESTIMATED_LIBRARY_SIZE"       "TOTAL_READS"                 
# [49] "PF_READS"                     "PF_READS_ALIGNED"            
# [51] "PCT_PF_READS_ALIGNED"         "PF_ALIGNED_BASES.1"          
# [53] "PF_HQ_ALIGNED_READS"          "PF_HQ_ALIGNED_BASES"         
# [55] "PF_HQ_ALIGNED_Q20_BASES"      "PF_MISMATCH_RATE"            
# [57] "PF_HQ_ERROR_RATE"             "PF_INDEL_RATE"               
# [59] "READS_ALIGNED_IN_PAIRS"       "STRAND_BALANCE"              
# [61] "PCT_CHIMERAS"                 "PCT_ADAPTER"                 
# [63] "Seq.PC1"                      "Seq.PC2"                     
# [65] "Seq.PC3"                      "Seq.PC4"                     
# [67] "Seq.PC5"

# (1a.) blacklist sample 19-2 from datExpr and targets. sample 19-2 has
# NAs in dx and other metadata. -------
rownames(targets) <- targets$SampleID
print(which(is.na(targets), arr.ind=TRUE))
targets <- subset(targets, targets$SampleID != "19-2")
rawExpr <- subset(rawExpr, select=(colnames(rawExpr) != "19-2"))

# (1b.) drop whitespace from Region - mid insula -> mid-insula
targets[["Region"]] <- gsub("\\s", "-", targets[["Region"]])

# (1c.) Reorder Dx such that Control comes first
targets[["Primary.Neuropath.Dx"]] <- relevel(targets[["Primary.Neuropath.Dx"]], "Control")

# (1d.) Write out unscaled targets
saveRDS(targets, file.path(out_path, "targets_unscaled.rds"))

# (2) Separate into metadata / sequencing and scale / convert to factors ---------
meta_cols <- c("Autopsy.ID","Source","Primary.Neuropath.Dx","Age",
			   "Sex","PMI","Batch","Region", "RIN")
seq_cols <- c("RNA_conc","Total_RNA_yield",
			  "PF_BASES","PF_ALIGNED_BASES","RIBOSOMAL_BASES","CODING_BASES","UTR_BASES","INTRONIC_BASES","INTERGENIC_BASES",
			  "CORRECT_STRAND_READS","INCORRECT_STRAND_READS","PCT_RIBOSOMAL_BASES","PCT_CODING_BASES","PCT_UTR_BASES",
			  "PCT_INTRONIC_BASES","PCT_INTERGENIC_BASES","PCT_MRNA_BASES","PCT_USABLE_BASES","PCT_CORRECT_STRAND_READS",
			  "MEDIAN_CV_COVERAGE","MEDIAN_5PRIME_BIAS","MEDIAN_3PRIME_BIAS","MEDIAN_5PRIME_TO_3PRIME_BIAS","TOTAL_CLUSTERS",
			  "ALIGNED_READS","AT_DROPOUT","GC_DROPOUT","UNPAIRED_READS_EXAMINED","UNMAPPED_READS","READ_PAIR_DUPLICATES",
			  "READ_PAIR_OPTICAL_DUPLICATES","PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE","TOTAL_READS","PF_READS",
			  "PF_READS_ALIGNED","PCT_PF_READS_ALIGNED","PF_ALIGNED_BASES.1","PF_HQ_ALIGNED_READS","PF_HQ_ALIGNED_BASES",
			  "PF_HQ_ALIGNED_Q20_BASES","PF_MISMATCH_RATE","PF_HQ_ERROR_RATE","PF_INDEL_RATE","READS_ALIGNED_IN_PAIRS",
			  "STRAND_BALANCE","PCT_CHIMERAS","PCT_ADAPTER")
targets_meta <- targets[, meta_cols]
targets_seq <- targets[, seq_cols]

targets_meta[, "Age"] <- scale(as.numeric(targets_meta[, "Age"]))
targets_meta[, "PMI"] <- scale(as.numeric(targets_meta[, "PMI"]))
targets_meta[, "RIN"] <- scale(as.numeric(targets_meta[, "RIN"]))
targets_meta[, "Autopsy.ID"] <- droplevels(as.factor(targets_meta[, "Autopsy.ID"]))
targets_meta[, "Batch"] <- droplevels(as.factor(targets_meta[, "Batch"]))
targets_meta[, "Source"] <- droplevels(as.factor(targets_meta[, "Source"]))
targets_meta[, "Primary.Neuropath.Dx"] <- droplevels(as.factor(targets_meta[, "Primary.Neuropath.Dx"]))
targets_meta[, "Region"] <- droplevels(as.factor(targets_meta[, "Region"]))
targets_meta[, "Sex"] <- droplevels(as.factor(targets_meta[, "Sex"]))

targets_seq <- apply(targets_seq, 2, as.numeric)
targets_seq <- apply(targets_seq, 2, scale)
targets_seq <- as.data.frame(targets_seq)

rownames(targets_meta) <- targets$SampleID
rownames(targets_seq) <- targets$SampleID

# (2a). Pre-quantile filter genes such that the gene set is consistent across regions.
qnt <- function(x) { return( x[apply(x,1,quantile,0.8) > 5, ] ) }
normExpr <- qnt(rawExpr)

targets_meta_region <- split(targets_meta, targets$Region)
targets_seq_region <- split(targets_seq, targets$Region)
normExpr_region_cols <- split(colnames(normExpr), targets$Region)
normExpr_region <- lapply(normExpr_region_cols, function(x) { return (normExpr[, x])} )

targets_region <- list()
targets_region$meta <- targets_meta_region
targets_region$seq <- targets_seq_region
saveRDS(targets_region, file.path(out_path, "targets_regions.rds"))
saveRDS(normExpr_region, file.path(out_path, "normExpr_regions.rds"))

saveRDS(targets_meta, file.path(out_path, "targets_meta.rds"))
saveRDS(targets_seq, file.path(out_path, "targets_seq.rds"))
saveRDS(normExpr, file.path(out_path, "normExpr.rds"))
saveRDS(rawExpr, file.path(out_path, "rawExpr.rds"))
