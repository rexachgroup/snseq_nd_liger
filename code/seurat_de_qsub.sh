#!/bin/bash

# Damon Polioudakis
# 2019-06-03
# job submission script for R scripts

# to submit this script:
# qsub seurat_de_qsub.sh -p [input seurat object] -t [output csv]
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N srtde
#$ -o logs/seurat_de_qsub_$JOB_NAME_$JOB_ID_$TASK_ID.log
#$ -e logs/seurat_de_qsub_$JOB_NAME_$JOB_ID_$TASK_ID.error
#$ -l h_data=128G,h_rt=12:00:00
# #$ -pe shared 8
#$ -t 1-27
#$ -tc 4
# #$ -m bea
# #$ -hold_jid

# reminder: make /logs directory in code directory
################################################################################
echo ""
echo "Starting seurat_de_qsub.sh ${SGE_TASK_ID}... "$(date)
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
rscript_path=/u/local/apps/R/3.6.0/gcc-4.9.3_MKL-2018/bin/Rscript

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
${rscript_path} seurat_de.R ${SGE_TASK_ID} ${path_seurat} ${seurat_obj_name} ${cluster_col_name} ${path_tmp}

# compile seurat_de output - hold until array jobs are complete
if [ $SGE_TASK_ID == $SGE_TASK_FIRST ]; then
  qsub -cwd -S /bin/bash -V -N srt_de_c -o logs/seurat_de_compile_qsub_$JOB_ID.log -e logs/seurat_de_compile_qsub_$JOB_ID.error -l h_data=128G,h_rt=1:00:00 -hold_jid srtde seurat_de_compile_qsub.sh -p ${path_seurat} -s ${seurat_obj_name} -c ${cluster_col_name} -t ${out_table}
fi
##########################################################################

echo ""
echo "End of seurat_de_qsub.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
