#' @import dplyr
#' @title PKM method
#' @description Normalize RNASeq counts by PKM
#' @details Computes the RPKM / FPKM that normalizes RNASeq data (RPKM for single reads and FPKM for paired reads)
#' and it returns in a R object (pkm_sample).
#' @details
#' 1. Count up the total reads in a sample and divide that number by 1,000,000. This is your 'per million' scaling factor.
#'
#' 2. Divide the read counts by the 'per million' scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM).
#'
#' 3. Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM/FPKM.
#' @param raw_counts Object containing raw_counts table - with addinfo () -
#' @param genes Object containing gene lenght table -with  extract () -
#' @references
#' Misuse of RPKM or TPM normalization when comparing across samples and sequencing protocols. Zhao, S.; Stanton, R. RNA (2020), Vol. 26, No. 8; DOI: \href{https://doi.org/10.1261/rna.074922.120 }{10.1261/rna.074922.120 }
#'
#' \href{https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/}{https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/}
#'
#' \href{https://www.youtube.com/watch?v=TTUrtCY2k-w}{https://www.youtube.com/watch?v=TTUrtCY2k-w}
#' @examples
#' # We prepare the sample data:
#' c(1,2,10,55,100) -> genes
#' c(0,5,10,20,50) -> sampleA
#' c(0,0,15,30,50) -> sampleB
#' c(0,1,2,3,100) -> sampleC
#' data.frame("GeneID"= genes,"A_counts"= sampleA,"B_counts"=sampleB,"C_counts"=sampleC) -> test_counts
#' hg38 # When you use the function extract(), you can load a gene data file from Rsubread package
#' pkm_method(test_counts, gene_lenght_file)
#' ####################################################################
#' addinfo(data_RNASeq_raw_counts) # Internal dataset
#' hg38 # When you use the function extract(), you can load a gene data file from Rsubread package
#' pkm_method(raw_counts_addinfo,gene_lenght_file)
#' ####################################################################
#' # Using hg38 dataset from jjtools instead of extract() function:
#' addinfo(data_RNASeq_raw_counts) # Internal dataset
#' as.data.frame(hg38)->hg38dt
#' datastart <- data.frame("GeneID"=hg38dt$GeneID,"Start"=hg38dt$Start)
#' dataend <- data.frame("GeneID"=hg38dt$GeneID,"End"=hg38dt$End)
#' setDT(datastart)[ , .SD[which.min(Start)], by = "GeneID"]  -> min_start
#' setDT(dataend)[ , .SD[which.max(End)], by = "GeneID"]  -> max_start
#' left_join(max_start,min_start, by="GeneID")->joined
#' data.frame("GeneID"=joined$GeneID, "Gene_lenght"=joined$End-joined$Start) -> gene_length
#' gene_length$Gene_lenght/1000 -> gene_length$kb_Gene_lenght
#' pkm_method(raw_counts_addinfo,gene_length)
#' @author Jose Jordan
#' @export
pkm_method <- function(raw_counts,genes) {
  left_join(raw_counts, genes, by="GeneID")->joi
  joi %>% dplyr::select(where(is.numeric)) ->joi2
  joi %>% dplyr::select(where(is.character)) ->joi3
  n = ncol(joi2)
  colSums(joi2[,2:(n-2)])/1000000 -> sca
  t(joi2[,2:(n-2)])/sca -> rpm
  t(rpm)/joi2$kb_Gene_lenght -> pkm1
  cbind("GeneID"=c(joi$GeneID, as.character()), pkm1, joi3) -> pkm_sampl
  as_tibble(pkm_sampl) ->> pkm_sample
  print(pkm_sample)
  rm(joi,joi2,joi3,sca,rpm,pkm1,pkm_sampl)
  message('Calculated RPKM / FPKM value of each gene per sample')
  message('Available in **pkm_sample** object')
}
