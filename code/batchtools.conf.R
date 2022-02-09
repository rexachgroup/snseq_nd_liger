library(rprojroot)
root <- rprojroot::has_file(".git/index")
ORION_SLURM_TMPL = root$find_file("code/slurm.tmpl")
#HOFFMAN_SGE_TMPL = "sge.tmpl"

is_orion <- function() dir.exists("/geschwindlabshares")
is_nessie <- function() system("hostname", intern = TRUE) == "nessie.bmap.ucla.edu"

if (is_orion()) {
    compress = FALSE
    cluster.functions <- makeClusterFunctionsSlurm(template = ORION_SLURM_TMPL)
}

if (is_nessie()) {
    compress = FALSE
    cluster.functions = makeClusterFunctionsMulticore()
}
