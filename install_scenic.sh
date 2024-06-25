#source ext/scenic/bin/activate
Rscript -e 'devtools::install_github("aertslab/SCENIC", ref = "v1.1.2", upgrade = "never")'
#pip install pyscenic
curl -L 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather' -o resources/cisTopic_databases/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather
curl -L 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather' -o resources/cisTopic_databases/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather
