source('./COMMON_DEFS.R')
suppressPackageStartupMessages(library(WGCNA, quietly=TRUE))
options(stringAsFactors = FALSE)

base_dir <- getBaseDir()
input_dir <- file.path(base_dir, 'Step01_PCA_OutlierRemoval')
out_path <- file.path(base_dir, 'Step01_PCA_OutlierRemoval')

# (1) load in files. ---------

targets <- readRDS(file.path(input_dir, 'targets_regions.rds'))
normExpr <- readRDS(file.path(input_dir, 'normExpr_regions.rds'))
targets_meta_all <- readRDS(file.path(input_dir, "targets_meta.rds"))
targets_seq_all <- readRDS(file.path(input_dir, "targets_seq.rds"))
normExpr_all <- readRDS(file.path(input_dir, "normExpr.rds"))

# (2.) filter sample connectivity / pca z-score. --------

#' Filter datExpr based on biweight midcorrelation sample connectivity.
#' 
#' @param datExpr: a gene expression dataframe with rows = genes, columns = samples.
#' @param threshold: a threshold; samples above this number of standard deviations
#' from mean connectivity are marked TRUE in the return value.
#' @param verbose: print some output about the outliers.
#' @return outliers: a binary vector indicating whether a given sample is an outlier.
gene_bicor_conn <- function(datExpr, threshold=2, verbose=FALSE) {
	normadj <- (0.5+0.5*bicor(datExpr))

	## Calculate connectivity
	netsummary <- fundamentalNetworkConcepts(normadj)
	ku <- netsummary$Connectivity
	z.ku <- ku-(mean(ku))/sqrt(var(ku))
	## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
	outliers <- (z.ku > mean(z.ku)+threshold*sd(z.ku))|(z.ku < mean(z.ku)-threshold*sd(z.ku))
	if (verbose) {
		print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep=""))
		print(colnames(datExpr)[outliers])
		print(table(outliers))
	}
	return(outliers)
}


#' Filter datExpr based on z-score from mean gene pca.
#'
#' @param datExpr: a gene expression dataframe with rows = genes, columns = samples.
#' @param threshold: a threshold; samples above this number of standard deviations
#' from the mean principal component are marked TRUE in the return .
#' @param verbose: print some output about the outliers.
#' @return outliers: a binary vector indicating whether a given sample is an outlier.
gene_pca_zscore <- function(datExpr, threshold=3, verbose=FALSE) {
	outliers <- vector(length=ncol(datExpr)) 
	exprPCs <- getTopPCs(scale(datExpr, scale=F), 5)
	for (i in 1:ncol(exprPCs)) {
		zscore <- (exprPCs[, i] - mean(exprPCs[, i])) / sqrt(var(exprPCs[, i]))
		outl <- which(abs(zscore) > threshold)
		outliers[outl] <- TRUE
	}
	return(outliers)
}


regions <- length(targets$meta)
targets_rmsampl <- targets
datExpr_rmsampl <- normExpr

# datExpr_rmsampl_new <- lapply(normExpr, outliers,
#     FUN = function(expr, out) {
#         expr_rm <- expr[,!(colnames(expr) %in% out)]
#         return(expr_rm)
#     })
# 
# targets_rmsampl_new <- list(
#     meta = lapply(targets$meta, outliers, 
#         FUN = function(t_meta, out) {
#             include <- which(!(rownames(t_meta) %in% out))
#             meta_rm <- as.data.frame(t_meta[include, ])
#             return(meta_rm)
#         }),
#     seq = lapply(targets$seq, outliers,
#         FUN = function(t_seq, out) {
#             include <- which(!(rownames(t_seq) %in% out))
#             seq_rm <- as.data.frame(t_seq[include, ])
#             return(seq_rm)
#         })
#     )
for (i in 1:regions) {
	outliers <- gene_bicor_conn(normExpr[[i]]) | gene_pca_zscore(normExpr[[i]])
	outliers <- names(which(outliers))
    print(outliers)
	datExpr_rmsampl[[i]] <- normExpr[[i]][,!(colnames(normExpr[[i]]) %in% outliers)]
	datExpr_rmsampl[[i]] <- datExpr_rmsampl[[i]][unique(rownames(datExpr_rmsampl[[i]])), ]

	include <- which(!rownames(targets_rmsampl$meta[[i]]) %in% outliers)
	targets_rmsampl$meta[[i]] <- as.data.frame(targets$meta[[i]][include,])
	targets_rmsampl$seq[[i]] <- as.data.frame(targets$seq[[i]][include,])
}

saveRDS(datExpr_rmsampl, file.path(out_path, 'datExpr_regions_rmsampl.rds'))
saveRDS(targets_rmsampl, file.path(out_path, 'targets_regions_rmsampl.rds'))

targets <- do.call(rbind, targets_rmsampl$meta)
rownames(targets) <- gsub('^.+\\.(.+)', '\\1', rownames(targets))
saveRDS(targets, file.path(out_path, "targets.rds"))

# (3) filter sequencing variables based on collinearity --------

collinearity <- function(targets_seq) {
	compare_seq <- list()

	for (i in colnames(targets_seq)) {
		testSeq1 <- targets_seq[, i]
		compare_seq[[i]] <- list()
		for (j in colnames(targets_seq)) {
			testSeq2 <- targets_seq[, j]
			if (i != j) {
				r2 <- summary(lm(testSeq1 ~ testSeq2))$adj.r.squared
				if (r2 >= 0.95) {
					compare_seq[[i]] <- c(compare_seq[[i]], j)
				}
			}
		}
	}


	seqPCs <- prcomp(t(scale(targets_seq, scale=F)), center=FALSE)$rotation[, 1:5]

	# (3b) for vars that are similar to each other: --------
	# loop through and compare confounds to the sequencing pcs
	# the var with the highest rsq is selected and all others are excluded
	# from target_seq_rmout_var.rds
	exclude <- list()
	for (ex_var in names(compare_seq)) {
		if (!ex_var %in% exclude && length(compare_seq[[ex_var]]) != 0) {
			keep0 = which(!compare_seq[[ex_var]] %in% exclude)
			ex_vars_test = c(ex_var,compare_seq[[ex_var]][keep0])
			r2_mat = matrix(NA,nrow=5,ncol=length(ex_vars_test))
			colnames(r2_mat)= ex_vars_test
			for(i in 1:5){
				for(j in 1:length(ex_vars_test)){
					ex_var <- which(colnames(targets_seq) == ex_vars_test[j])
					mod = summary(lm(seqPCs[,i] ~ targets_seq[,ex_var]))
					r2_mat[i,j]=mod$adj.r.squared
				}
			}
			r2_sum = apply(r2_mat,2,sum)
			max_val = max(r2_sum)
			max_ex_var = names(r2_sum)[which(r2_sum==max_val)]
			if(length(max_ex_var) > 1){
				keep = sample(max_ex_var,1)
			}else{
				keep = max_ex_var
			}
			compare_seq[[keep]] = compare_seq[[keep]][-c(1:length(compare_seq[[keep]]))]
			exclude = c(exclude, ex_vars_test[-which(ex_vars_test %in% keep)])
			print(length(exclude))
		}
	}
	return(exclude)
}

# (3) writeout rmvar targets --------

targets_rmvar <- targets_rmsampl
for (i in 1:regions) {
	exclude <- collinearity(targets_rmsampl$seq[[i]])
	targets_rmvar$seq[[i]] <- targets_rmsampl$seq[[i]][, which(!colnames(targets_rmsampl$seq[[i]]) %in% exclude)]
}
saveRDS(targets_rmvar, file.path(out_path, 'targets_regions_rmvar.rds'))

# (2lmer) filter on all regions to create input to s04_lmer.R. ------

outliers <- gene_bicor_conn(normExpr_all) | gene_pca_zscore(normExpr_all)
datExpr_all <- normExpr_all[,!outliers]
targets_meta_rmsampl_all <- targets_meta_all[!outliers,]
targets_seq_rmsampl_all <- targets_seq_all[!outliers,]

saveRDS(datExpr_all, file.path(out_path, "datExpr_all.rds"))
saveRDS(targets_meta_rmsampl_all, file.path(out_path, "targets_meta_rmsampl_all.rds"))
saveRDS(targets_seq_rmsampl_all, file.path(out_path, "targets_seq_rmsampl_all.rds"))

# (4) export outlier samples

gene_bicor_outliers_region <- mapply(normExpr, names(normExpr), FUN = function(expr, region) {
    df = data.frame(sample = names(which(gene_bicor_conn(expr))))
    if (nrow(df) > 0) {
        df$Autopsy.ID <- targets_meta_all[as.character(df$sample), "Autopsy.ID"]
        df$region <- region
        df$filter <- "connectivity"
    }
    return(df)
}, SIMPLIFY = FALSE)
gene_pca_outliers_region <- mapply(normExpr, names(normExpr), FUN = function(expr, region) {
    df <- data.frame(sample = names(which(gene_pca_zscore(expr))))
    if (nrow(df) > 0) {
        df$Autopsy.ID <- targets_meta_all[as.character(df$sample), "Autopsy.ID"]
        df$region <- region
        df$filter <- "pca"
    }
    return(df)
})
bor_df <- Reduce(rbind, gene_bicor_outliers_region)
por_df <- Reduce(rbind, gene_pca_outliers_region)
region_outliers <- rbind(bor_df, por_df)
write.csv(bor_df, file.path(out_path, "gene_bicor_outliers.csv"))
write.csv(por_df, file.path(out_path, "gene_pca_outliers.csv"))

gene_bicor_outliers <- names(which(gene_bicor_conn(normExpr_all)))
gene_pca_outliers <- names(which(gene_pca_zscore(normExpr_all)))

b_out_df <- data.frame(sample = gene_bicor_outliers, filter = "connectivity")
b_out_df$Autopsy.ID = targets_meta_all[as.character(b_out_df$sample), "Autopsy.ID"]
p_out_df <- data.frame(sample = gene_pca_outliers)
p_out_df$Autopsy.ID = targets_meta_all[p_out_df$sample, "Autopsy.ID"]
write.csv(bor_df, file.path(out_path, "gene_bicor_outliers_allregion.csv"))
write.csv(por_df, file.path(out_path, "gene_pca_outliers_allregion.csv"))
