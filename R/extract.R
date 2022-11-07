#' @import readr
#' @import dplyr
#' @import data.table
#' @title Extract gene length
#' @description Function that extracts the genes and the length of the transcript.
#' @details It extracts from the "RefSeq" file of the Rsubread package the genes and length of the transcript and
#' returns it in an R object (gene_lenght_file).
#' @examples
#' hg38 # When you use the function extract(), you can load a gene data file from Rsubread package
#' @author Jose Jordan
#' @export
extract <- function() {
  message('Select the file to extract each gene length')
  system.file("annot", package = "Rsubread") -> pathRsubread
  choose.files(default=paste0(pathRsubread, "/hg38_RefSeq_exon.txt")) -> fileref_seq
  hg38_RefSeq_exon <-read_delim(fileref_seq,
                                "\t", escape_double = FALSE, trim_ws = TRUE)
  message('CALCULATING...(O.o o.O)')
  datastart <- data.frame("GeneID"=hg38_RefSeq_exon$GeneID,"Start"=hg38_RefSeq_exon$Start)
  dataend <- data.frame("GeneID"=hg38_RefSeq_exon$GeneID,"End"=hg38_RefSeq_exon$End)
  setDT(datastart)[ , .SD[which.min(Start)], by = "GeneID"]  -> min_start
  setDT(dataend)[ , .SD[which.max(End)], by = "GeneID"]  -> max_start

  left_join(max_start,min_start, by="GeneID")->joined_max_min_gene

  data.frame("GeneID"=joined_max_min_gene$GeneID, "Gene_lenght"=joined_max_min_gene$End-joined_max_min_gene$Start) -> gene_lentgh
  gene_lentgh$Gene_lenght/1000 -> gene_lentgh$kb_Gene_lenght
  gene_lentgh ->> gene_lenght_file
  print(as_tibble(gene_lenght_file))
  rm(joined_max_min_gene,max_start,min_start,dataend,datastart,hg38_RefSeq_exon,fileref_seq,pathRsubread,gene_lentgh)
  message('In the object **gene_lenght_file** you have for each "GeneID" the length of the transcript (Gene_lenght columns)')
}
