source("./COMMON_DEFS.R")
library(parallel)
library(abind)
library(progress)
library(lme4)

base_dir <- getBaseDir()
out_dir <- file.path(base_dir, "Step04_Regression")

ncores <- detectCores()
cl <- makeCluster(ncores, type = "FORK")
# memory usage is largely determined by the number of processes spawned; reduce ncores if necessary

datExpr <- readRDS(file.path(base_dir, "Step01_PCA_OutlierRemoval", "datExpr_regions_rmsampl.rds"))
targets <- readRDS(file.path(base_dir, "Step01_PCA_OutlierRemoval", "targets_regions_rmsampl.rds"))
targets_seqpc <- list()
seq_names <- paste0("seqpc", 1:5)
for (region in names(targets$meta)) {
    targets_seqpc[[region]] <- as.data.frame(getTopPCs(t(scale(targets$seq[[region]])), 5))
    colnames(targets_seqpc[[region]]) <- seq_names
}

# construct per-gene linear mixed effects model --------


# Fit model. --------

## export datExpr, targets, targets_seqpc to cluster.
## hopefully this saves memory copies when doing parLapply.
clusterExport(cl = cl, list("datExpr", "targets", "targets_seqpc"))
clusterEvalQ(cl = cl, gc())

coef_list <- list()
system.time(
    for (region in names(targets$meta)) {
        crange <- 1:nrow(datExpr[[region]])
        coef_list[[region]] <- parLapply(cl, crange, function(i, region) {
            targets_seqpc <- as.data.frame(getTopPCs(t(scale(targets$seq[[region]])), 5))
            colnames(targets_seqpc) <- paste0("seqpc", 1:5)

            targets_meta <- targets$meta[[region]]
            targets_region <- cbind(targets_meta, targets_seqpc)

            datExpr_col <- as.numeric(datExpr[[region]][i, ])
            lm_dat <- cbind(datExpr_col, targets_region)

            form <- as.formula(paste0("datExpr_col ~ Primary.Neuropath.Dx + Batch + RIN + PMI + Sex + Age + Source + ", paste0("seqpc", 1:5, collapse = " + ")))
            lm1 <- lm(form, data = lm_dat)
            return(coef(summary(lm1)))
        }, region)
    }
)
saveRDS(coef_list, file.path(out_dir, "coef_list.rds"))

stopCluster(cl)

# coef_list <- readRDS(file.path(out_dir, "coef_list.rds"))
# Correct for sequencing effects. --------
# Stop + remake cluster; regression is single-core, so we want
# all cores to be in use.
cl <- makeCluster(detectCores(), type = "FORK")
clusterEvalQ(cl = cl, gc())

datExpr_region <- list()
system.time(
    for (region in names(targets$meta)) {
        crange <- 1:nrow(datExpr[[region]])
        batch_binary <- model.matrix(~ 0 + targets$meta[[region]]$Batch)
        clusterExport(cl = cl, list("datExpr", "targets", "targets_seqpc", "coef_list", "batch_binary", "seq_names"))
        coef_corr_mat <- parSapply(cl, crange, function(i, region) {
            targets_meta <- targets$meta[[region]]
            targets_region <- cbind(targets$meta[[region]], targets_seqpc[[region]])

            datExpr_col <- as.numeric(datExpr[[region]][i, ])
            age_effect <- as.matrix(targets_region[, "Age"]) %*% as.matrix(coef_list[[region]][[i]]["Age", "Estimate"])
            batch_effect <- as.matrix(batch_binary[, -1]) %*% as.matrix(coef_list[[region]][[i]][c("Batch2", "Batch3"), "Estimate"])
            pmi_effect <- as.matrix(targets_region[, "PMI"]) %*% as.matrix(coef_list[[region]][[i]]["PMI", "Estimate"])
            rin_effect <- as.matrix(targets_region[, "RIN"]) %*% as.matrix(coef_list[[region]][[i]]["RIN", "Estimate"])
            sex_effect <- as.numeric(targets_region[, "Sex"]) %*% as.matrix(coef_list[[region]][[i]]["SexM", "Estimate"])
            source_effect <- as.numeric(targets_region[, "Source"]) %*% as.matrix(coef_list[[region]][[i]]["SourceUPENN", "Estimate"])
            seq_effect <- as.matrix(targets_region[, seq_names]) %*% as.matrix(coef_list[[region]][[i]][seq_names, "Estimate"])

            return(datExpr_col - age_effect - batch_effect - pmi_effect - rin_effect - sex_effect - source_effect - seq_effect)
        }, region)
        datExpr_region[[region]] <- t(coef_corr_mat)
        datExpr_region[[region]] <- apply(datExpr_region[[region]], c(1, 2), as.numeric)
        rownames(datExpr_region[[region]]) <- rownames(datExpr[[region]])[crange]
		colnames(datExpr_region[[region]]) <- colnames(datExpr[[region]])
    }
)
datExpr_reg <- abind(datExpr_region, along=2)
saveRDS(datExpr_region, file.path(out_dir, "datExpr_region.rds"))
saveRDS(datExpr_reg, file.path(out_dir, "datExpr_reg.rds"))

stopCluster(cl)

print(gc())
