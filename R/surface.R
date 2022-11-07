#' @import dplyr
#' @title Surfaceome filter
#' @description Surfaceome filter function using \link{surfy} list
#' @details Filter the list counts by a list of 2911 surface proteins and viral receptors and it returns in a R object.
#' @param counts Object containing RNASeq counts. GeneID column is required to apply the function.
#' @references
#' Original surfaceome list of 2886 proteins available in: \href{http://wlab.ethz.ch/surfaceome/}{http://wlab.ethz.ch/surfaceome/}
#'
#' The in silico human surfaceome. Bausch-Fluck,D.; Goldmann, .U; Muller,S.; van Oostrum, M.; Muller, M.; Schubert,O.T.; Wollscheid, B. PNAS. Nov 2018, 115(46) E10988-E10997; DOI: \href{https://doi.org/10.1073/pnas.1808790115}{10.1073/pnas.1808790115}
#'
#' 25 viral receptors added to the original surfaceome list, founded in: \href{https://viralzone.expasy.org/5356}{https://viralzone.expasy.org/5356}
#' @examples
#' surface(data_RNASeq_raw_counts)
#' ####################################################################
#' addinfo(data_RNASeq_raw_counts) # Internal dataset
#' hg38 # When you use the function extract(), you can load a gene data file from Rsubread package (in this case we will use hg38 data)
#' tpm_method(raw_counts_addinfo,gene_lenght_file)
#' surface(tpm_sample) # Result: TPM table filtered by 'surfy' dataset
#' @author Jose Jordan
#' @export
surface <- function(counts) {
  counts$GeneID <- as.numeric(counts$GeneID)
  counts %>% inner_join(surfy) ->> counts_filtered_by_surfaceome
  counts_filtered_by_surfaceome
  message('Result available in **counts_filtered_by_surfaceome** object')
  }
