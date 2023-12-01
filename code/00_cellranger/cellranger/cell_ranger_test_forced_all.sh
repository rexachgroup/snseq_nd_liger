#!/bin/bash

# Damon Polioudakis
# 2019-12-03
# run cell ranger atac count to process fastqs

# to submit this script:
# qsub cell_ranger_atac_process_data.sh -p [path to tsv args file]
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -N CRatac_tst
#$ -o logs/cell_ranger_atac_process_data_$JOB_ID_$TASK_ID.log
#$ -e logs/cell_ranger_atac_process_data_$JOB_ID_$TASK_ID.error
#$ -l h_data=24G,h_rt=36:00:00,highp
#$ -pe shared 8
#$ -t 1-52
#$ -tc 6
# #$ -m bea

# tsv args file saved as excel default tsv UTF-8
# run as array job, set qsub -t option to number of samples
# reminder: make /logs directory in code directory

# to copy html summary files to one folder and label with sample id
# for dir in /u/project/geschwind/drewse/g_singlecell/cellranger/data/*/*; do
#   id=$(basename ${dir})
#   cp /u/project/geschwind/drewse/g_singlecell/cellranger/data/*/$id/outs/web_summary.html ../analysis/cell_ranger/20190815_web_summaries/"$id"_web_summary.html
# done

# # to compile csv statistics files
# for dir in /u/project/geschwind/drewse/g_singlecell/cellranger/data/*/*; do
#   id=$(basename ${dir})
#   in_file=/u/project/geschwind/drewse/g_singlecell/cellranger/data/*/${id}/outs/metrics_summary.csv
#   out_file="../analysis/cell_ranger/20190815_metrics_summary.csv"
#   if [ ! -f $out_file ]; then
#     awk 'FNR == 1 {print "ID",$0}' OFS=, ${in_file} > ${out_file};
#     awk -v id=$id 'FNR == 2 {print id,$0}' OFS=, ${in_file} >> ${out_file}
#   else
#     awk -v id=$id 'FNR == 2 {print id,$0}' OFS=, ${in_file} >> ${out_file}
#   fi
# done
################################################################################
echo ""
echo "Starting cell_ranger_atac_process_data.sh ${SGE_TASK_ID}... "$(date)
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
  echo "  [--force-cells cell#/N]"
  echo "  [path to out directory]"
}
################################################################################

### main

## variables

#project_dir="/u/scratch/m/mdknight/cell_ranger_atac/202008"
# arguments
# required
in_args_tsv=
# required from tsv file of paramaters
id=
reference=
fastqs=
outdir=
# optional
force_cells=
# SGE_TASK_ID=1

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
	NR == taskID {print $1; print $2; print $3; print $4; print $5; print $6}' RS='\r\n' ../tmp/tmp_${SGE_TASK_ID}.tsv)

# assign args from tsv to variables
id=${argsA[0]}
reference=${argsA[1]}
fastqs=${argsA[2]}
sample=${argsA[3]}
force_cells=${argsA[4]}
outdir=${argsA[5]}
# remove quotes (awk adds quotes to string with comma for some reason)
fastqs=$(echo ${fastqs} | sed 's/"//g')

echo "Arguments:"
echo "tsv of cell ranger parameters: ${in_args_tsv}"
echo "ID: ${id}"
echo "path to reference genome: ${reference}"
echo "path to fastqs directory: ${fastqs}"
echo "sample: ${sample}"
echo "output directory: ${outdir}"
echo "number of cells to force: ${force_cells}"
echo ""

# check arguments
# -z switch will test if the expansion of "$1" is a null string or not. If it is
# a null string then the body is executed.
if  [ -z "${in_args_tsv}" ] || \
		[ -z "${reference}" ] || \
		[ -z "${fastqs}" ] || \
    [ -z "${sample}" ] || \
		[ -z "${force_cells}" ] || \
		[ -z "${outdir}" ]
  then
    usage "Missing required arguments"
    exit 1
fi

## run cell ranger count

# output directory (cell ranger outputs in dir it is executed from)
mkdir -p ${outdir}
cd ${outdir}

# run cell ranger count
# Arguments:
#     id      A unique run id, used to name output folder [a-zA-Z0-9_-]+.
#     fastqs  Path of folder created by mkfastq or bcl2fastq.
#     sample  Prefix of the filenames of FASTQs to select.
# Options:
# --reference=PATH    Path of folder containing a 10x-compatible reference.
#                       Required.
# --force-cells=N     Define the top N barcodes with the most reads as
#                       cells. N must be a positive integer <= 20,000. Please
#                       consult the documentation before using this option.
if  [[ "${force_cells}" =~ ^[0-9]+$ ]]
  then
    echo "Calling cellranger-atac count"
    echo "Using --force-cells parameter"
    echo ""
    /u/home/d/dpolioud/bin/cellranger-atac-1.2.0/cellranger-atac count \
      --id=${id} \
      --reference=${reference} \
      --fastqs=${fastqs} \
      --sample=${sample} \
      --force-cells=${force_cells} \
      --localcores=8 \
      --localmem=192
  else
    echo "Calling cellranger-atac count"
    echo ""
    /u/home/d/dpolioud/bin/cellranger-atac-1.2.0/cellranger-atac count \
      --id=${id} \
      --reference=${reference} \
      --fastqs=${fastqs} \
      --sample=${sample} \
      --localcores=8 \
      --localmem=192
fi

#
# /u/home/d/dpolioud/bin/cellranger-atac-1.2.0/cellranger-atac count \
#   --id=PFC1_1_at_ \
#   --reference=/u/nobackup/dhg/common_resources/ref_genomes/cellranger/refdata-cellranger-atac-GRCh38-1.2.0 \
#   --fastqs=/u/project/geschwind/dpolioud/nucseq_nd/raw_data/fastq/20191004/PFC1_1_at/ \
#   --sample=PFC1_1_at
#   --force-cells=7000
#
# /u/home/d/dpolioud/bin/cellranger-atac-1.2.0/cellranger-atac count --id=AT2_ --reference=/u/project/geschwind/bwamsley/programs/refdata-cellranger-atac-GRCh38-1.1.0 --fastqs=/u/project/geschwind/bwamsley/ATAC/AT2 --sample=AT2

# /u/home/d/dpolioud/bin/cellranger-atac-1.2.0/cellranger-atac upload damonp@g.ucla.edu PFC1_1_at/PFC1_1_at.mri.tgz

# cleanup
#rm ${project_dir}/tmp/tmp_${SGE_TASK_ID}.tsv
##########################################################################

echo ""
echo "End of cell_ranger_atac_process_data.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
