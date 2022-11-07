#' @import dplyr
#' @title RPM method
#' @description Normalize RNASeq counts by RPM.
#' @details Computes the RPM that normalizes RNASeq data (also called CPM)
#' and it returns in a R object (cpm_sample). RPM does not consider the transcript length normalization.
#' @details
#' 1. Multiply number of reads mapped to gene by 1,000,000.
#'
#' 2. Divide this value by total number of mapped reads.
#' @param raw_counts Object containing raw_counts table - with addinfo () -
#' @references
#' \href{https://www.reneshbedre.com/blog/expression_units.html}{https://www.reneshbedre.com/blog/expression_units.html}
#' @examples
#' # We prepare the sample data:
#' c(1,2,10,55,100) -> genes
#' c(0,5,10,20,50) -> sampleA
#' c(0,0,15,30,50) -> sampleB
#' c(0,1,2,3,100) -> sampleC
#' data.frame("GeneID"= genes,"A_counts"= sampleA,"B_counts"=sampleB,"C_counts"=sampleC) -> test_counts
#' rpm_method(test_counts)
#' ####################################################################
#' addinfo(data_RNASeq_raw_counts) # Internal dataset
#' rpm_method(raw_counts_addinfo)
#' @author Jose Jordan
#' @export
rpm_method <- function(raw_counts) {
  raw_counts$GeneID <- as.character(raw_counts$GeneID)
  raw_counts %>% dplyr::select(where(is.numeric)) ->joi2
  raw_counts %>% dplyr::select(where(is.character)) ->joi3
  joi2*1000000 -> pre_rpm
  colSums(joi2)-> pre_rpm_col_sums
  t(pre_rpm)/(pre_rpm_col_sums)->jj
  t(jj) -> cpm
  cbind(joi3,cpm) -> pre_cpm_sample
  as_tibble(pre_cpm_sample) ->> cpm_sample
  print(cpm_sample)
  rm(joi2,joi3,pre_rpm,pre_rpm_col_sums,jj,cpm,pre_cpm_sample)
  message('Calculated RPM (= CPM) value of each gene per sample')
  message('Available in **cpm_sample** object')
}
