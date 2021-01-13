ORION_SLURM_TMPL = "slurm.tmpl"
HOFFMAN_SGE_TMPL = "sge.tmpl"

is_orion <- function() dir.exists("/geschwindlabshares")

if (is_orion()) {
    cluster.functions <- makeClusterFunctionsSlurm(template = "slurm.tmpl")
}
