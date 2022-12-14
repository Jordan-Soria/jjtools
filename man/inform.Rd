% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inform.R
\name{inform}
\alias{inform}
\title{More informatives}
\usage{
inform(x, Filter_value)
}
\arguments{
\item{x}{Object in R that contains the table (raw counts, tpm counts, pkm counts...) to apply the filter.}

\item{Filter_value}{Cut-off value. Optional. By default = 1.}
}
\description{
More informatives in a dataset
}
\details{
Computes a binary table from a table counts using a cut-off value (= 1 BY DEFAULT
if you don't specify a value).
If Filter_value is = or > will be computed as 1. If is < will be computed as 0.

Then, the function computes the 20 most informative columns / samples / cells of the dataset,
showing you on the screen the percentage increase of information and the list
of the 20 most informative cells/groups (ordered from highest to lowest).

\strong{IMPORTANT!: YOU MUST USE A MATRIX ONLY WITH DATA}
}
\examples{
addinfo(data_RNASeq_raw_counts) # Internal dataset
hg38 # When you use the function extract(), you can load a gene data file from Rsubread package
tpm_method(raw_counts_addinfo,gene_lenght_file)
# To use inform() you MUST use only the matrix with the dataset
tpm_sample[,2:6] -> new_sample
inform(new_sample,10000)
# The most informative columns are S1,S2,S3 and S5 (S4 is not informative).
# S1.1, ... appears because we have less than 20 columns in the dataset.
# Don't worry, ignore.
}
\author{
Jose Jordan
}
