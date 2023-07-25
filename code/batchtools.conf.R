library(rprojroot)
root <- rprojroot::has_file(".git/index")
ORION_SLURM_TMPL = root$find_file("code/slurm.tmpl")
HOFFMAN_SGE_TMPL = root$find_file("code/sge.tmpl")

is_orion <- function() dir.exists("/geschwindlabshares")
is_hoffman <- function() dir.exists("/u/project/geschwind")
is_nessie <- function() system("hostname", intern = TRUE) == "nessie.bmap.ucla.edu"
compress = FALSE

if (is_orion()) {
    cluster.functions <- makeClusterFunctionsSlurm(template = ORION_SLURM_TMPL, array.jobs = TRUE)
}

if (is_hoffman()) {
    cluster.functions <- makeClusterFunctionsSGE(template = HOFFMAN_SGE_TMPL) 
}

if (is_nessie()) {
    cluster.functions = makeClusterFunctionsMulticore()
}
