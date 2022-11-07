#' @import dplyr
#' @title TPM method
#' @description Normalize RNASeq counts by TPM
#' @details Computes the TPM that normalizes RNASeq data and it returns in a R object (tpm_sample).
#' @details
#' 1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
#'
#' 2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your 'per million' scaling factor.
#'
#' 3. Divide the RPK values by the 'per million' scaling factor. This gives you TPM.
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
#' test_counts
#' extract() # When you use the function extract(), you can load a gene data file from Rsubread package (in this case we will use hg38 data)
#' tpm_method(test_counts, gene_lenght_file)
#' ####################################################################
#' addinfo(data_RNASeq_raw_counts) # Internal dataset
#' extract() # When you use the function extract(), you can load a gene data file from Rsubread package (in this case we will use hg38 data)
#' tpm_method(raw_counts_addinfo,gene_lenght_file)
#' # Testing that the sum of each column is = 1000000, for example with column S1:
#' # sum(tpm_sample$S1)
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
#' tpm_method(raw_counts_addinfo,gene_length)
#' # Testing that the sum of each column is = 1000000, for example with column S1:
#' # sum(tpm_sample$S1)
#' @author Jose Jordan
#' @export
tpm_method <- function(raw_counts,genes) {
  left_join(raw_counts, genes, by="GeneID")->joi
  joi %>% dplyr::select(where(is.numeric)) ->joi2
  joi %>% dplyr::select(where(is.character)) ->joi3
  n = ncol(joi2)
  joi2[,2:(n-2)]/joi2$kb_Gene_lenght -> rpk
  colSums(rpk)/1000000->sca
  t(rpk)/sca->ttpm
  t(ttpm)->tpm1
  cbind("GeneID"=c(joi$GeneID, as.character()), tpm1, joi3) -> tpm_sampl
  as_tibble(tpm_sampl) ->> tpm_sample
  print(tpm_sample)
  rm(joi,joi2,joi3,sca,rpk,ttpm,tpm1,tpm_sampl)
  message('Calculated TPM value of each gene per sample')
  message('Available in **tpm_sample** object')
}
