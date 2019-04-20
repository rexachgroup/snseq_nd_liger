#!/bin/bash

# Damon Polioudakis
# 2018-09-07
# run cell ranger aggregate to combine aligned samples into one counts matrix

# to submit this script:
# qsub cell_ranger_aggregate_data.sh [-f path to csv args file] [-o out_dir] [-i run id]
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N CellR_Ag
#$ -o logs/cell_ranger_aggregate_data_$JOB_ID.log
#$ -e logs/cell_ranger_aggregate_data_$JOB_ID.error
#$ -l h_data=256G,h_rt=24:00:00
#$ -m bea

# csv args file saved as excel default CSV UTF-8
# Reminder: make /logs directory in code directory
################################################################################
echo ""
echo "Starting cell_ranger_aggregate_data.sh ... "$(date)
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
  echo "  [-f cell ranger csv of arguments]"
  echo "  [-o output directory]"
}
################################################################################

### run cell ranger aggregate

## variables

project_dir="/u/home/d/dpolioud/project-geschwind/nucseq_nd"

## handle arguments

in_args=
out_dir=
while [ $# -gt 0 ]
do
  case "$1" in
    -i) id="$2"; shift;;
  	-f) in_args="$2"; shift;;
    -o) out_dir="$2"; shift;;
  	--)	shift; break;;
    # usage message and terminates if an unknown command line flag starting with a
    # dash was specified
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
if [ -z "${id}" ] || [ -z "${in_args}" ] || [ -z "${out_dir}" ]
  then
    usage "Missing required arguments: -i -f -o"
    exit 1
fi

echo "Arguments:"
echo "run ID: ${id}"
echo "cell ranger csv arguments: ${in_args}"
echo "output directory: ${out_dir}"

# csv file arguments

# create a CSV file with a header line containing the following columns:

# library_id: Unique identifier for this input library. This will be used for
# labeling purposes only; it doesn't need to match any previous ID you've
# assigned to the library.

# molecule_h5: Path to the molecule_info.h5 file produced by cellranger count.
# So if you processed your library by calling cellranger count --id=ID in some
# directory /DIR, this path would be /DIR/ID/outs/molecule_info.h5.

# remove <U+FEFF> character from csv files
mkdir -p ../tmp
sed 's/\xEF\xBB\xBF//' < ${in_args} > ../tmp/tmp.csv


## run cell ranger aggregate

# The cellranger aggr command takes a CSV file specifying a list of cellranger
# count output files (specifically the molecule_info.h5 from each run), and
# produces a single gene-barcode matrix containing all the data.

# When combining multiple libraries, the barcode sequences for each library are
# distinguished by their GEM group (see Gem Groups).

# By default, each library's reads are subsampled such that all libraries have
# the same effective sequencing depth, measured in terms of reads per cell (see
# Depth Normalization).

# cell ranger outputs in dir it is executed from
# output directory
mkdir -p ${out_dir}
cd ${out_dir}

cellranger aggr \
  --id=${id} \
  --csv=${project_dir}/tmp/tmp.csv \
  --normalize=mapped
  # --normalize=none

# cleanup
rm /u/home/d/dpolioud/project-geschwind/nucseq_nd/tmp/tmp.csv
################################################################################

echo ""
echo "End of cell_ranger_aggregate_data.sh ... "$(date)
################################################################################
