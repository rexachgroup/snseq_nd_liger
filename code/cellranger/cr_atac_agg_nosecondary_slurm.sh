#!/bin/bash

# Damon Polioudakis
# 2018-09-07
# run cell ranger aggregate to combine aligned samples into one counts matrix
# modified by Misty OCT 2020
# Lawrence Chen: convert to orion / SLURM, Dec 2020

# to submit this script:
# sbatch cell_ranger_atac_aggregate_data.sh [-f path to csv args file] [-o out_dir] [-i run id]
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N CRat_ag
#$ -o logs/cell_ranger_atac_aggregate_data_$JOB_ID.log
#$ -e logs/cell_ranger_atac_aggregate_data_$JOB_ID.error
#$ -l h_data=32G,h_rt=180:00:00,highp,highmem
#$ -pe shared 8
#$ -m bea
#SBATCH --nodes 1
#SBATCH --cpus-per-task 8
#SBATCH -t 160:00:00
#SBATCH -o logs/cell_ranger_atac_aggregate_data_$SLURM_JOBID.log
#SBATCH --export=ALL
#SBATCH --mem 128G

################################################################################
echo ""
echo "Starting cell_ranger_atac_aggregate_data_slurm.sh ... "$(date)
echo ""
################################################################################

### functions

#usage () {
#  message=$1
#  if [ -n "${message}" ]
#    then
#      echo "${message}"
#  fi
#  echo "usage: $0"
#  echo "  [-i run ID]"
#  echo "  [-f cell ranger csv of arguments]"
#  echo "  [-o output directory]"
#}
################################################################################

### run cell ranger aggregate

## variables

#project_dir="/u/project/geschwind/mdknight/ATAC"
id=ATAC_nosecondary
in_args=/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/code/cell_ranger_agg_args_four.csv
out_dir=/geschwindlabshares/lchenprj01/nucseq_nd_dpolioud/analysis/
CELLRANGER_ATAC=/geschwindlabshares/lchenprj01/software/seq/cellranger-atac-1.2.0/cellranger-atac
CELLRANGER_REFERENCE=/geschwindlabshares/lchenprj01/software/seqdata/cellranger-atac/refdata-cellranger-atac-GRCh38-1.2.0

## handle arguments

#in_args=
#out_dir=
#while [ $# -gt 0 ]
#do
#  case "$1" in
#    -i) id="$2"; shift;;
#  	-f) in_args="$2"; shift;;
#    -o) out_dir="$2"; shift;;
#  	--)	shift; break;;
    # usage message and terminates if an unknown command line flag starting with a
    # dash was specified
#  	-*)
#      echo >&2 usage "Invalid command line flag"
#      exit 1;;
#    *) break;;	# terminate while loop
#  esac
#  shift
#done
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

# fragments: Path to the fragments.tsv.gz file produced by cellranger-atac count. For example, if you processed your GEM well by calling cellranger-atac count --id=ID in some directory /DIR, the fragments would be /DIR/ID/outs/fragments.tsv.gz.

# cells: Path to the singlecell.csv file produced by cellranger-atac count.

# (Optional) peaks: Path to the peaks.bed file produced by cellranger-atac count.

# (Optional) Additional custom columns containing library meta-data (e.g., lab or sample origin). These custom library annotations do not affect the analysis pipeline but can be visualized downstream in the Loupe Cell Browser. Note that unlike other CSV inputs to Cell Ranger ATAC, these custom columns may contain characters outside the ASCII range (e.g., non-Latin characters).

# remove <U+FEFF> character from csv files
#TMPFILE=$PWD/tmp/tmp.csv
#mkdir -p $PWD/tmp
#sed 's/\xEF\xBB\xBF//' < ${in_args} > $TMPFILE
# mkdir -p ../tmp
# sed 's/\xEF\xBB\xBF//' < ${in_args} > ../tmp/tmp.csv


## run cell ranger atac aggregate

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

${CELLRANGER_ATAC} aggr \
  --id=${id} \
  --csv=${in_args} \
  --reference=${CELLRANGER_REFERENCE} \
  --normalize=none \
  --localmem=128 \
  --localcores=8 \
  --nosecondary

# cleanup
# rm /u/home/d/dpolioud/project-geschwind/nucseq_nd/tmp/tmp.csv
################################################################################

echo ""
echo "End of cell_ranger_atac_aggregate_data_slurm.sh ... "$(date)
################################################################################
