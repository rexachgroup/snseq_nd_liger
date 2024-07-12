# Introduction
This repo contains code for the RNA analysis for Rexach et al. (2024). For the scATAC analysis see [rexachgroup/nd_atac_archr](https://github.com/rexachgroup/nd_atac_archr)
# Installation / Running
Dependencies for this code are specified using a conda environment; create and activate the environment by running
```
conda create -p conda/ -f environment.yml
conda activate conda/
```

Some steps require R packages that are not available from `conda-forge`  or `bioconda`. Install the dependencies by running `install_doubletfinder.sh` , `install_scenic.sh` after the conda environment is created.
# analysis/ folder organization
Subcomponents of the analysis are organized in run order under `code/`:
## 00_cellranger
initial cellranger count and aggregate runs per region.
## 01_doubletfinder
doublet detection and filtering.
## 10_seurat_cluster
Initial Seurat preprocessing, Seurat clustering, and celltype assignment.
## 11_liger_subcluster
Subclustering of Seurat celltypes using LIGER.
## 12_scenic
SCENIC runs on Seurat celltypes astrocyte, excitatory, microglia, oligodendrocyte per region.
## 13_seurat_metadata_prep
Miscellaneous metadata operations.
- generate filter list with subclusters above 10umi
- testing of metadata across dx
- Import neurodegeneration scores from bulk metadata file
## 14_seurat_azimuth
Celltype annotation of Seurat clusters using the [Azimuth human motor cortex reference](https://azimuth.hubmapconsortium.org/references/#Human - Motor Cortex).
## 20_seurat_cluster_dge
DGE analyses at Seurat cluster level.
## 21_liger_subcluster_dge
DGE analyses at Liger subcluster level.
## 23_seurat_excitatory_layers_dge
Re-calling excitatory cells by neuronal layer and differential expression across layers.
## 30_subcluster_counts_composition
Subcluster cell type composition testing.
## 31_liger_subcluster_hier
Hierarchical clustering of subclusters based on DGE beta.
## 40_bulk_deconvolution
Bulk RNA sample filtering, analysis, and Bisque / CIBERSORTx deconvolution.
