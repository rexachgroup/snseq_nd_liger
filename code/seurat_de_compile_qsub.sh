#!/bin/bash

# Damon Polioudakis
# 2019-06-03
# Script to run seurat_de.R with compile flag

################################################################################
echo ""
echo "Starting seurat_de_compile_qsub.sh ${SGE_TASK_ID}... "$(date)
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
  echo "  [-p path to seurat object]"
  echo "  [-t path to output DE csv table]"
  echo "  [-s name of seurat object]"
  echo "  [-c name of column in seurat object metadata with cluster IDs]"
}
################################################################################

### main

## variables

# arguments
# required
# path to seurat object
# /u/flashscratch/d/dpolioud/seurat/20190427/p1_p2_p3_p4_filtered_rmp27_test.rdat
path_seurat=
# output path for compiled csv table
# /u/flashscratch/d/dpolioud/seurat/20190427/p1_p2_p3_p4_filtered_rmp27_cluster_enriched_de.csv
out_table=
# "p1_p2_p3_p4_so"
seurat_obj_name=
# "RNA_snn_res.0.6"
cluster_col_name=

# directory for temporary output
path_tmp=/u/flashscratch/d/dpolioud/tmp/seurat/$(date +%Y%m%d)/cluster_${SGE_TASK_ID}.csv
# path to Rscript
rscript_path=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript

## handle arguments

while [ $# -gt 0 ]; do
  case "$1" in
  	-p) path_seurat="$2"; shift;;
    -t) out_table="$2"; shift;;
    -s) seurat_obj_name="$2"; shift;;
    -c) cluster_col_name="$2"; shift;;
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
if [ -z "${path_seurat}" ] | [ -z "${out_table}" ] | [ -z "${seurat_obj_name}" ] | [ -z "${cluster_col_name}" ]
  then
    usage "Missing required arguments"
    exit 1
fi

echo "Arguments:"
echo "Path to seurat object: ${path_seurat}"
echo "Seurat object name: ${seurat_obj_name}"
echo "Cluster column name: ${cluster_col_name}"
# echo "Directory for temporary output: ${path_tmp}"
echo "Output path for compiled csv table: ${out_table}"

## submit jobs

## run Rscript
${rscript_path} seurat_de.R "NA" ${path_seurat} ${seurat_obj_name} ${cluster_col_name} ${path_tmp} ${out_table} "TRUE"
##########################################################################

echo ""
echo "End of seurat_de_compile_qsub.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
