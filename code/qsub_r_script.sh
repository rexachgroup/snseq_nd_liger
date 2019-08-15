#!/bin/bash

# Damon Polioudakis
# 2018-11-15
# job submission script for R scripts

# to submit this script:
# qsub qsub_r_script.sh -p [path to R script] options: -r [path to Rscript]
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N Rscript
#$ -o logs/qsub_r_script_seurat_plots_ad_ctrl_$JOB_NAME_$JOB_ID_$TASK_ID.log
#$ -e logs/qsub_r_script_seurat_plots_ad_ctrl_$JOB_NAME_$JOB_ID_$TASK_ID.error
#$ -l h_data=128G,h_rt=12:00:00
# #$ -pe shared 8
#$ -t 1
# #$ -m bea
# #$ -hold_jid

# reminder: make /logs directory in code directory
################################################################################
echo ""
echo "Starting qsub_r_script.sh ${SGE_TASK_ID}... "$(date)
echo ""
################################################################################

### functions

usage () {
  message=$1
  if [ -n "${message}" ]
    then
      echo "${message}"
  fi
  echo "usage: $0"
  echo "  [-p path to the R script (script to run)]"
  echo "  [-r path to Rscript (for using different versions of R)]"
}
################################################################################

### main

## variables

# path to Rscript
rscript_path=/u/local/apps/R/3.6.0/gcc-4.9.3_MKL-2018/bin/Rscript
# arguments
# required
# path to script to run
in_r_script=

## handle arguments

while [ $# -gt 0 ]; do
  case "$1" in
  	-p) in_r_script="$2"; shift;;
    -r) rscript_path="$2"; shift;;
  	--)	shift; break;;
    # usage message and terminates if an unknown command line flag starting with
    # a dash was specified
  	-*)
      echo >&2 usage "Invalid command line flag"
      exit 1;;
    *) break;;	# terminate while loop
  esac
  shift
done
# all command line switches are processed,
# "$@" contains all file names

# check arguments
# -z switch will test if the expansion of "$1" is a null string or not. If it is
# a null string then the body is executed.
if  [ -z "${in_r_script}" ]
  then
    usage "Missing required arguments"
    exit 1
fi

echo "Arguments:"
echo "Path to the R script to run: ${in_r_script}"
echo "Optional arguments:"
echo "Path to Rscript: ${rscript_path}"

## run Rscript
${rscript_path} ${in_r_script} ${SGE_TASK_ID}
##########################################################################

echo ""
echo "End of qsub_r_script.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
