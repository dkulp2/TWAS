#' Sample transcriptome data for package \code{twas}.
#'
#' \code{twas.tx} is a list containing two data.frames as follows:
#'
#' @format transcriptome - A data.frame of 922 rows with one ID column and 4 columns of gene expression:
#' \describe{
#'   \item{ID}{subject identifier}
#'   \item{COL5A2}{COL5A2 RNASeq transcript counts}
#'   \item{LRRC16A}{LRRC16A RNASeq counts}
#'   \item{SLC9A8}{SLC9A8 RNASeq counts}
#'   \item{PDIA5}{PDIA5 RNASeq counts}
#' }
#'
#' genotypes - A data.frame of 922 rows with one ID column and genotypes for 30 markers
#' \describe{
#'   \item{ID}{subject identifier}
#'   \item{rs13392501}{genotypes for marker rs13392501}
#'   \item{etc...}
#' }
#' 
"twas.tx"

#' Sample trait data for package \code{twas}.
#'
#' \code{twas.traits} is a list containing two data.frames as follows:
#'
#' @format traits - A data.frame of 174 rows with one ID column and 2 traits
#' \describe{
#'   \item{ID}{subject identifier}
#'   \item{logCRP}{log CRP}
#'   \item{logIL6}{log IL6}
#' }
#'
#' genotypes - A data.frame of 174 rows with one ID column and genotypes for 30 markers
#' \describe{
#'   \item{ID}{subject identifier}
#'   \item{rs13392501}{genotypes for marker rs13392501}
#'   \item{etc...}
#' }
#' 
"twas.traits"

#' Sample sets of markers proximal to genes
#'
#' \code{twas.gene.markers} is a list of character vectors of marker
#' names for four sample gene names. The gene symbol is the name of
#' each character vector.
"twas.gene.markers"

