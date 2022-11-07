#' @import dplyr
#' @import Matrix
#' @import matrixStats
#' @title Statistics
#' @description Calcule PKM & TPM statistics
#' @details Computes statistics from normalized counts and it returns in a R object (norm_counts_added_stats).
#' @details
#' |-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-|
#'
#' Mean: Calculated as the arithmetic mean, the sum of all of the numbers divided by the number of numbers (n).
#'
#' Standard Deviation (StDev): Calculated as a quantity expressing by how much the members of a group
#' differ from the mean value for the group. It is the square root of the sum of the difference between each
#' value with respect to the mean divided by n-1.
#'
#' Mean Standard Error (SEM): Calculated as the value that quantifies how much the values deviate from
#' the population mean. In other words, it quantifies the oscillations of the sample mean. It results from dividing
#' Standard Deviation by the square root of n.
#' @param x Object containing raw_counts table normalized by TPM or PKM - tpm_method() or pkm_method()-
#' @references
#'
#' @examples
#' # We prepare the sample data:
#' c(1,2,10,55,100) -> genes
#' c(0,5,10,20,50) -> sampleA
#' c(0,0,15,30,50) -> sampleB
#' c(0,1,2,3,100) -> sampleC
#' data.frame("GeneID"= genes,"A_counts"= sampleA,"B_counts"=sampleB,"C_counts"=sampleC) -> test_counts
#' hg38 # When you use the function extract(), you can load a gene data file from Rsubread package (in this case we will use hg38 data)
#' tpm_method(test_counts, gene_lenght_file)
#' calc(tpm_sample)
#' @author Jose Jordan
#' @export
calc <- function(x) {
  x$GeneID <- as.character(x$GeneID)
  x %>% mutate(Mean = rowMeans(across(where(is.numeric)))) -> x_mean
  x %>% mutate(StDev = rowSds(as.matrix(across(where(is.numeric))))) -> x_devest
  x$Mean <- x_mean$Mean
  x$StDev <- x_devest$StDev
  x$SEM <- x_devest$StDev/sqrt(nrow(x_devest))
  x ->> norm_counts_added_stats
  print(norm_counts_added_stats)
  message('Calculated "Mean", "Standard deviation" and "Mean Standard Error" by GeneID')
  rm(x_mean,x_devest)
  message('Available in **norm_counts_added_stats** object')
  }
