#!/bin/bash

# Damon Polioudakis
# 2018-09-03
# run cell ranger count to process fastqs

# to submit this script:
# qsub cell_ranger_process_data.sh -p [path to tsv args file]
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N CellR_PD
#$ -o logs/cell_ranger_process_data_$JOB_ID_$TASK_ID.log
#$ -e logs/cell_ranger_process_data_$JOB_ID_$TASK_ID.error
#$ -l h_data=192G,h_rt=12:00:00,highp,highmem
# #$ -pe shared 8
#$ -t 1-4
#$ -tc 6
# #$ -m bea

# tsv args file saved as excel default tsv UTF-8
# run as array job, set -t option to number of samples
# reminder: make /logs directory in code directory

# to copy html summary files to one folder and label with sample id
# for dir in ../data/20190306/P*; do id=$(basename ${dir}); cp ../data/20190306/$id/outs/web_summary.html ../data/20190306/web_summaries/"$id"_web_summary.html; done
################################################################################
echo ""
echo "Starting cell_ranger_process_data.sh ${SGE_TASK_ID}... "$(date)
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
  echo "  [-p path to tsv file of cell ranger parameters]"
  echo "tsv required arguments:"
  echo "  [--id]"
  echo "  [--transcriptome [path to reference]]"
  echo "  [--fastqs [path to fastqs]]"
  echo "  [--sample]"
  echo "  [path to out directory]"
  echo "tsv optional arguments:"
  echo "  [--expect-cells [number of expected cells]]"
  echo "  [--force-cells [number of cells to force]]"
}
################################################################################

### main

## variables

cellranger_bin="/u/project/geschwind/chenlo/Software/SeqBuilds/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger"
project_dir="/u/project/geschwind/chenlo/nucseq-atac"
# arguments
# required
in_args_tsv=
# required from tsv file of paramaters
id=
reference=
fastqs=
sample=
outdir=
# optional
expect_cells=
force_cells=


## handle arguments

while [ $# -gt 0 ]; do
  case "$1" in
  	-p) in_args_tsv="$2"; shift;;
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

# make tmp directory
mkdir -p ../tmp

# remove <U+FEFF> character from tsv files
sed 's/\xEF\xBB\xBF//' < ${in_args_tsv} > ../tmp/tmp_${SGE_TASK_ID}.tsv

# read arguments from tsv
# RS='\r\n' to clean windows line ends
read -r -a argsA <<< $(awk -v taskID="${SGE_TASK_ID}" 'BEGIN {FS="\t"; OFS=" ";}
	NR == taskID {print $1; print $2; print $3; print $4; print $5; print $6; print $7}' RS='\r\n' ../tmp/tmp_${SGE_TASK_ID}.tsv)

# assign args from tsv to variables
id=${argsA[0]}
reference=${argsA[1]}
fastqs=${argsA[2]}
sample=${argsA[3]}
expect_cells=${argsA[4]}
outdir=${argsA[5]}
force_cells=${argsA[6]}
# remove quotes (awk adds quotes to string with comma for some reason)
fastqs=$(echo ${fastqs} | sed 's/"//g')

# check arguments
# -z switch will test if the expansion of "$1" is a null string or not. If it is
# a null string then the body is executed.
if  [ -z "${in_args_tsv}" ] || \
		[ -z "${reference}" ] || \
		[ -z "${fastqs}" ] || \
		[ -z "${sample}" ] || \
		[ -z "${expect_cells}" ] || \
		[ -z "${outdir}" ]
  then
    usage "Missing required arguments"
    exit 1
fi

echo "Arguments:"
echo "tsv of cell ranger parameters: ${in_args_tsv}"
echo "ID: ${id}"
echo "path to reference genome: ${reference}"
echo "path to fastqs directory: ${fastqs}"
echo "sample ID: ${sample}"
echo "number of expected: ${expect_cells}"
echo "output directory: ${outdir}"
echo "number of cells to force: ${force_cells}"

## run cell ranger count

# output directory (cell ranger outputs in dir it is executed from)
mkdir -p ${outdir}
cd ${outdir}

# run cell ranger count
${cellranger_bin} count \
  --id=${id} \
  --transcriptome=${reference} \
  --fastqs=${fastqs} \
  --sample=${sample} \
  --force-cells=${force_cells}
  # --expect-cells=${expect_cells} \

  # ${cellranger_bin} count \
  #   --id=P2_7B \
  #   --transcriptome=/u/nobackup/dhg/common_resources/ref_genomes/cellranger/refdata-cellranger-GRCh38-3.0.0_premrna \
  #   --fastqs=/u/home/d/dpolioud/project-geschwind/nucseq_nd/raw_data/fastq/20190214/HYJCKBBXX/JR001 \
  #   --sample=JR001 \
  #   --force-cells=7000

# cleanup
rm ${project_dir}/tmp/tmp_${SGE_TASK_ID}.tsv
##########################################################################

echo ""
echo "End of cell_ranger_process_data.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
