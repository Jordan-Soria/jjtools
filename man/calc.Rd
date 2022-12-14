% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc.R
\name{calc}
\alias{calc}
\title{Statistics}
\usage{
calc(x)
}
\arguments{
\item{x}{Object containing raw_counts table normalized by TPM or PKM - tpm_method() or pkm_method()-}
}
\description{
Calcule PKM & TPM statistics
}
\details{
Computes statistics from normalized counts and it returns in a R object (norm_counts_added_stats).

|-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-||-|

Mean: Calculated as the arithmetic mean, the sum of all of the numbers divided by the number of numbers (n).

Standard Deviation (StDev): Calculated as a quantity expressing by how much the members of a group
differ from the mean value for the group. It is the square root of the sum of the difference between each
value with respect to the mean divided by n-1.

Mean Standard Error (SEM): Calculated as the value that quantifies how much the values deviate from
the population mean. In other words, it quantifies the oscillations of the sample mean. It results from dividing
Standard Deviation by the square root of n.
}
\examples{
# We prepare the sample data:
c(1,2,10,55,100) -> genes
c(0,5,10,20,50) -> sampleA
c(0,0,15,30,50) -> sampleB
c(0,1,2,3,100) -> sampleC
data.frame("GeneID"= genes,"A_counts"= sampleA,"B_counts"=sampleB,"C_counts"=sampleC) -> test_counts
hg38 # When you use the function extract(), you can load a gene data file from Rsubread package (in this case we will use hg38 data)
tpm_method(test_counts, gene_lenght_file)
calc(tpm_sample)
}
\references{

}
\author{
Jose Jordan
}
