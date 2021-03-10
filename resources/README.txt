################################################################################

disease_gene_markers_jessica_20200429.csv

Disease gene marker list curated by Jessica.

https://mail.google.com/mail/u/0/#label/Follow-up/FMfcgxwHNCspDbXjrnDnvgckGFpBTChM
################################################################################

grubman_2019_st5_AD_GWAS_genes.xlsx

Alexandra Grubman biorxiv 2019

Supplementary Table 4 Genes associated with AD and AD-related traits by GWAS, grouped by GWAS trait (AD, LOAD, biomarkers, neuropathologic change)
################################################################################

chen_2019_WGCNA_modules.csv

WGCNA gene modules from Chen 2019 biorxiv paper on the AD mouse model.

"Attachment is Table S3 for WGCNA modules. PIGs are the purple module."

https://www.biorxiv.org/content/10.1101/719930v1?ct=
https://mail.google.com/mail/u/0/#inbox/FMfcgxwDqfNXGBFhcWJFqCVkdsDVGPsv
################################################################################

AstrocyteMarkers_Brie_20190110.xlsx

Astrocyte subtype markers from Jessica Rexach.
################################################################################

cluster_markers_mixed_20181019.csv

Cell type markers compiled by Jessica Rexach.
################################################################################

excitatory_markers_20191023.csv

Excitatory markers derived from Lake et al and Hodge et al.
################################################################################

inhibitory_markers_20200711.csv

Inhibitory markers compiled from Inma's list:

https://mail.google.com/mail/u/0/#label/Hold/FMfcgxwHNqDbsqmmkbBpDPChcjWwLflf
################################################################################

RNASEQ_SELECTIVE_GENES.csv

From Jessica:

https://mail.google.com/mail/u/0/#inbox/FFNDWLHwpLLMlnqtfMhtHhsJLrndjWpr
################################################################################

scenic_cistarget_databases

Downloaded with R and SCENIC package using these commands:

# Download the *.feather* files (and *.descr*, if available) for the relevant organism. [Aprox file size: 1GB]
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
"https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")

# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
  descrURL <- gsub(".feather$", ".descr", featherURL)
  if(file.exists(descrURL)) download.file(descrURL, destfile=basename(descrURL))
}
################################################################################