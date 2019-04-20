#!/bin/bash

# cellranger count \
#   --id=P1_1 \
#   --transcriptome=/u/nobackup/dhg/common_resources/ref_genomes/cellranger/refdata-cellranger-GRCh38-3.0.0_premrna \
#   --fastqs=/u/home/d/dpolioud/project-geschwind/nucseq_nd/raw_data/fastq/20180811/P1/P1_1,/u/home/d/dpolioud/project-geschwind/nucseq_nd/raw_data/fastq/20181115/P1_1 \
#   --sample=P1_1 \
#   --force-cells=6000

/u/home/d/dpolioud/bin/cellranger-2.2.0/cellranger count \
  --id=P1_2 \
  --transcriptome=/u/nobackup/dhg/common_resources/ref_genomes/cellranger/refdata-cellranger-GRCh38-3.0.0_premrna \
  --fastqs=/u/home/d/dpolioud/project-geschwind/nucseq_nd/raw_data/fastq/20180811/P1/P1_2,/u/home/d/dpolioud/project-geschwind/nucseq_nd/raw_data/fastq/20181115/P1_2 \
  --sample=P1_2 \
  --force-cells=6000





# cellranger count \
#   --id=P1_4 \
#   --transcriptome=/u/nobackup/dhg/common_resources/ref_genomes/cellranger/refdata-cellranger-GRCh38-3.0.0_premrna \
#   --fastqs=/u/home/d/dpolioud/project-geschwind/nucseq_nd/raw_data/fastq/20180811/P1/P1_4,/u/home/d/dpolioud/project-geschwind/nucseq_nd/raw_data/fastq/20181115/P1_4 \
#   --sample=P1_4 \
#   --force-cells=6000
