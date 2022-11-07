## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(jjtools)
library(dplyr)
library(readr)
library(data.table)

## -----------------------------------------------------------------------------
addinfo(data_RNASeq_raw_counts,hg38) #Selected hg38_RefSeq_exon dataset (in Rsubread folder)

## -----------------------------------------------------------------------------
# You can use the function extract() to obtain gene data from Rsubread package and obtain the gene_lenght_file"
datastart <- data.frame("GeneID"=hg38$GeneID,"Start"=hg38$Start)
dataend <- data.frame("GeneID"=hg38$GeneID,"End"=hg38$End)
setDT(datastart)[ , .SD[which.min(Start)], by = "GeneID"]  -> min_start
setDT(dataend)[ , .SD[which.max(End)], by = "GeneID"]  -> max_start

left_join(max_start,min_start, by="GeneID")->joined_max_min_gene

data.frame("GeneID"=joined_max_min_gene$GeneID, "Gene_lenght"=joined_max_min_gene$End-joined_max_min_gene$Start) -> gene_lentgh
gene_lentgh$Gene_lenght/1000 -> gene_lentgh$kb_Gene_lenght
gene_lentgh -> gene_lenght_file

transform(gene_lenght_file, Gene_lenght = as.character(Gene_lenght),
          kb_Gene_lenght = as.character(kb_Gene_lenght))->gene_lenght_file2
left_join(raw_counts_addinfo,gene_lenght_file2, 
          by="GeneID")->raw_counts_addinfo2
raw_counts_addinfo2[,c(1:6,11:12)] #See that the last two columns are 'chr'

## -----------------------------------------------------------------------------
tpm_method(raw_counts_addinfo,gene_lenght_file)
# Testing that the sum of each column is = 1000000, for example with column S1:
sum(tpm_sample$S1)

## -----------------------------------------------------------------------------
calc(tpm_sample)

## -----------------------------------------------------------------------------
norm_counts_added_stats[,c(1:3,7,11:13)]

## -----------------------------------------------------------------------------
surface(norm_counts_added_stats) #We use the result in the previous chunk -calc()-
counts_filtered_by_surfaceome

## -----------------------------------------------------------------------------
# First of all you MUST use only the matrix with the dataset, so you must remove columns without data:
tpm_sample[,2:6] -> new_sample
# Use a cut-off value to filter the dataset. In this case we will use 10000:
filter_value=10000
# Apply the function:
inform(new_sample,filter_value)

## -----------------------------------------------------------------------------
data_cells <- read_table2("https://www.ebi.ac.uk/gxa/experiments-content/E-PROT-25/resources/ExperimentDownloadSupplier.Proteomics/tsv", 
    col_names = FALSE, skip = 5)


## -----------------------------------------------------------------------------
replace(data_cells, is.na(data_cells), 0) -> file_line_cells_NA
WantedData=file_line_cells_NA[file_line_cells_NA$X2 %in% surfy$`Gene symbol UniProt`, ]
WantedData
#Now you can apply tpm_method() function in jjtools or other. Or use simply these values, as in this example to ilustrate the group() function in jjtools

## -----------------------------------------------------------------------------
group(WantedData,1000,250,10,337,3)
names_group_cells

