#!/bin/bash

# Damon Polioudakis
# 2018-09-07
# run cell ranger reanalyze command to normalize and cluster data

# to submit this script:
# qsub cell_ranger_reanalyze_data.sh [path to csv args file] [outdir]
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N CellR_Ag
#$ -o logs/cell_ranger_reanalyze_data_$JOB_ID.log
#$ -e logs/cell_ranger_reanalyze_data_$JOB_ID.error
#$ -l h_data=64G,h_rt=24:00:00
#$ -m bea

# reminder: make /logs directory in code directory
################################################################################
echo ""
echo "Starting cell_ranger_reanalyze_data.sh ... "$(date)
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
  echo "  [-i run ID]"
  echo "  [-m path to matrices_h5.h5]"
  echo "optional arguments:"
  echo "  [-p path to csv file of parameters]"
}
################################################################################

### run cell ranger reanalyze

## variables

# arguments
# required
id=
matrix=
# optional
parameters_csv=

## handle arguments

while [ $# -gt 0 ]
do
  case "$1" in
  	-i) id="$2"; shift;;
    -m) matrix="$2"; shift;;
    -p) parameters_csv="$2"; shift;;
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
if [ -z "${id}" ] || [ -z "${matrix}" ]
  then
    usage "Missing required arguments: -i -m"
    exit 1
fi

echo "Arguments:"
echo "run ID: ${id}"
echo "path to matrices_h5.h5: ${matrix}"
echo "(optional) path to a CSV file containing parameters: ${parameters_csv}"


## run cell ranger reanalyze

# the cellranger reanalyze command reruns secondary analysis performed on the
# gene-barcode matrix (dimensionality reduction, clustering and visualization)
# using different parameter settings.

cellranger reanalyze \
  --id=${id} \
  --matrix=${matrix} \
  --params=${parameters_csv}
  # --id=P1_combined \
  # --matrix=P1/outs/filtered_gene_bc_matrices_h5.h5 \
  # --params=P1_reanalysis.csv
################################################################################

echo ""
echo "End of cell_ranger_reanalyze_data.sh ... "$(date)
################################################################################
