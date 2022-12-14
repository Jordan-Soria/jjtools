% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rpm_method.R
\name{rpm_method}
\alias{rpm_method}
\title{RPM method}
\usage{
rpm_method(raw_counts)
}
\arguments{
\item{raw_counts}{Object containing raw_counts table - with addinfo () -}
}
\description{
Normalize RNASeq counts by RPM.
}
\details{
Computes the RPM that normalizes RNASeq data (also called CPM)
and it returns in a R object (cpm_sample). RPM does not consider the transcript length normalization.

1. Multiply number of reads mapped to gene by 1,000,000.

2. Divide this value by total number of mapped reads.
}
\examples{
# We prepare the sample data:
c(1,2,10,55,100) -> genes
c(0,5,10,20,50) -> sampleA
c(0,0,15,30,50) -> sampleB
c(0,1,2,3,100) -> sampleC
data.frame("GeneID"= genes,"A_counts"= sampleA,"B_counts"=sampleB,"C_counts"=sampleC) -> test_counts
rpm_method(test_counts)
####################################################################
addinfo(data_RNASeq_raw_counts) # Internal dataset
rpm_method(raw_counts_addinfo)
}
\references{
\href{https://www.reneshbedre.com/blog/expression_units.html}{https://www.reneshbedre.com/blog/expression_units.html}
}
\author{
Jose Jordan
}
