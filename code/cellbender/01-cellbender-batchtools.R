set.seed(0)
liblist <- c("Seurat", "tidyverse", "batchtools")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf)

in_seurat_meta <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_cellbender_tb <- "../../analysis/seurat_lchen/cellbender/cellranger_h5.csv"
batchtools <- "../../analysis/seurat_lchen/cellbender/batchtools"
cellbender_env <- "../../ext/cellbender/conda/"
OUT_DIR <- "../../analysis/seurat_lchen/cellbender/runs"

RESOURCES <- list(
    ncpus = 1,
    memory = 16,
    walltime = 6 * 60 * 60,
    gpu = "V100"
)

CELLBENDER_BIN <- str_glue("conda run --no-capture-output -p {cellbender_env} cellbender remove-background")

main <- function() {
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    }
    dir.create(OUT_DIR)
    cellbender_tb <- read_csv(in_cellbender_tb)
    meta <- readRDS(in_seurat_meta)

    cellbender_tb <- cellbender_tb %>%
        mutate(
            counts_h5 = str_glue("{dirname(molecule_h5)}/raw_feature_bc_matrix.h5"),
            out_dir = str_glue("{OUT_DIR}/{library_id}"),
            out_h5 = str_glue("{out_dir}/{library_id}_cellbender.h5")

        )
    cellbender_f <- cellbender_tb %>%
        filter(!file.exists(out_h5))
    print(duplicated(cellbender_f$counts_h5))
    cellbender_f <- cellbender_f %>% 
        mutate(
            cmd = pmap(list(counts_h5, out_h5), fmt_cellbender_cmd),
            cmd_path = str_glue("{out_dir}/cmd.txt")
        )
    map(cellbender_f$out_dir, dir.create)
    pwalk(cellbender_f, function(...) {
        cr <- list(...)
        writeLines(cr$cmd, cr$cmd_path)
    })

    clearRegistry()
    ids <- batchMap(sys_run, args = list(cellbender_f$cmd_path))
    submitJobs(findNotSubmitted(), RESOURCES)
    waitForJobs()
}

sys_run <- function(cmd_path) {
    cmd <- readLines(cmd_path)
    system(cmd)
}

fmt_cellbender_cmd <- function(h5_path, out_h5) {
    str_glue(
        CELLBENDER_BIN, 
        "--input",
        h5_path,
        "--output",
        out_h5,
        "--expected-cells 10000",
        "--total-droplets-included 50000",
        "--epochs 150",
        "--learning-rate 5e-05",
        "--cuda",
        "--posterior-batch-size 10",
        .sep = " "
    )
}

if (!interactive()) main()
