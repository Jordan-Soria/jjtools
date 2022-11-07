#' Surfaceome data
#'
#' Surfaceome data (modified) from surfaceome web page (see source).
#' The modifications consist of include the 'GeneID' of each protein, eliminate columns with low-relevance and added
#' some genes (25) not indicated in the original dataset, but they were identified in virus receptors (indicated in the 'Surfaceome Label Source' column with the target 'Viral list').
#'
#' @docType data
#'
#' @usage surfy
#'
#' @format A tibble: 2,911 x 9
#'
#' @keywords datasets
#'
#' @references Bausch-Fluck et al. (2018) PNAS. 115 (46) E10988-E10997, doi:
#' (\href{https://doi.org/10.1073/pnas.1808790115}{10.1073/pnas.1808790115})
#'
#' @source \href{http://wlab.ethz.ch/surfaceome/}{http://wlab.ethz.ch/surfaceome/} ........Surfaceome list (Bausch-Fluck et al. 2018)
#'
#' \href{https://viralzone.expasy.org/5356}{https://viralzone.expasy.org/5356} .....Viral receptors table
#'
#' @examples
#'
"surfy"
