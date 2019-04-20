# Damon Polioudakis
# 2018-10-18
# Plots of seurat analysis

# Must load modules:
#  module load gcc/7.2.0
#  module load R/3.5.1
################################################################################

rm(list = ls())
set.seed(27)

require(Seurat)
require(tidyverse)
require(SCENIC)
source("function_library.R")
source("ggplot_theme.R")
options(stringsAsFactors = FALSE)
sessionInfo()

## inputs
in_seurat <- paste0(
  "/u/flashscratch/d/dpolioud/p1_p2_seurat_subcluster/20190221"
  , "/p1_p2_seurat_subcluster_astro.rdat")

  in_seurat <- paste0(
    "/u/flashscratch/d/dpolioud/p1_p2_seurat_subcluster/20190227/cca_lib"
    , "/p1_p2_seurat_subcluster_astro_controls.rdat")

## Variables
date <- format(Sys.Date(), "%Y%m%d")
script_name <- "scenic.R"
graph_subtitle <- "Astrocyte control samples"

## outputs
out_graph <- paste0("../analysis/scenic/", date
  , "/graphs/subcluster_astro_remp27/scenic_")
out_table <- paste0("../analysis/scenic/", date
  , "/tables/subcluster_astro_remp27/scenic_")

# make directories
# dir.create(dirname(out_graph), recursive = TRUE)
# dir.create(dirname(out_table), recursive = TRUE)
################################################################################

load(in_seurat)
# exprMat <- p1_p2_astro_so@raw.data %>% as.matrix()

seurat_obj <- SubsetData(seurat_obj, subset.raw = TRUE, cells.use = sample(seurat_obj@cell.names)[1:200])

exprMat <- seurat_obj@raw.data %>% as.matrix()
metadata <- seurat_obj@meta.data

save(exprMat, metadata, file = "../tmp/scenic_test_data.rdat")
load("../tmp/scenic_test_data.rdat")

# Initialize SCENIC settings

org="hgnc" # or hgnc, or dmel
dbDir="../resources/scenic_cistarget_databases" # RcisTarget databases location
myDatasetTitle="SCENIC on P1 P2 astrocyte samples" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=10)
## Motif databases selected:
##   mm9-500bp-upstream-7species.mc9nr.feather
##   mm9-tss-centered-10kb-7species.mc9nr.feather

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- metadata
# scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
# saveRDS(scenicOptions, file="int/scenicOptions.Rds")
################################################################################


# (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)


runCorrelation <- function(exprMat_filtered,scenicOptions)
{
  corrMat <- cor(t(exprMat_filtered), method="spearman")
  saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
}
runCorrelation(exprMat_filtered, scenicOptions)


# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1)

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)




scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

# For a very quick run:
# coexMethod=c("top5perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
# save...

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
runSCENIC_3_scoreCells(scenicOptions, logMat)
