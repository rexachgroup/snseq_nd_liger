import logging
import os

import subprocess
import pandas as pd
import anndata as ad

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from dask.diagnostics import ProgressBar
from dask.distributed import LocalCluster, Client

IN_EXPR = "../../analysis/12_liger_scenic/insula/oligodendrocyte/expr.loom"
IN_GENES = "../../analysis/12_liger_scenic/insula/oligodendrocyte/expr_genes.tsv"
OUT_MAT = "../../analysis/12_liger_scenic/insula/oligodendrocyte/grn_network.txt"
GRNBOOST_SEED = 123
LOG = logging.getLogger()
LOG.setLevel(logging.DEBUG)

def main(): 
    LOG.info("dask init")
    local_cluster = LocalCluster(n_workers = 4, threads_per_worker = 8, memory_limit = "50GiB")
    local_client = Client(local_cluster)

    LOG.info("read ad")
    expr_obj = load_expr(IN_EXPR)
    LOG.info("tf load")
    genes = load_tf_names(IN_GENES)

    LOG.info("grnboost2")
    with ProgressBar():
        adjacencies = grnboost2(expr_obj.to_df(), tf_names = genes, verbose = True, seed = GRNBOOST_SEED, client_or_address = local_client)

    LOG.info("write adjacency network")
    adjacencies.to_csv(OUT_MAT) 

def load_expr(path):
    return(ad.read_loom(path, sparse = True))

if __name__ == "__main__":
    main()
