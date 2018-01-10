#' twas: A package for transcriptome-wide association studies
#'
#' The twas package provides a single function, \code{twas}, for
#' jointly modeling two independent data sets, trait data (genome-wide
#' association, GWA) and expression data (transcriptome-wide
#' association, TWA). The analytical goal is to leverage the
#' transcriptome data in analysis of the GWA data to inform the
#' mechanisms underlying gene-trait associations.
#' 
#' @docType package
#' @name twas
NULL

#' An EM model fit jointly to transcriptome and trait data
#'
#' \code{twas} returns a model fit from two data sets from different
#' subjects, one set of transcriptome data and the other of a trait
#' data. \code{twas} iteratively fits the trait data to the subjects
#' in the transcriptome data and vice versa, until convergence.
#'
#' @param transcriptome a data frame of expression levels for one or
#'   more genes. Each row is a subject. There is one ID column and the
#'   remaining columns are genes. Expression values are non-negative
#'   integers.
#' @param traits a data frame of one or more traits. Each row is a
#'   subject. There is one ID column and the remaining columns are
#'   traits. Trait values are numeric.
#' @param tx_genotypes a data frame of genotypes for the transcriptome
#'   subjects. Each row is a subject. There is one ID column with the
#'   same name as used in \code{transcriptome} and the remaining
#'   columns are markers. Genotype values are non-negative integers.
#' @param trait_genotypes a data frame of genotypes for the trait
#'   subjects. Each row is a subject. There is one ID column with the
#'   same name as used in \code{traits} and the remaining columns are
#'   markers. Genotype values are non-negative integers.
#' @param markers a character vector of the names of the markers to
#'   use for analysis. Defaults to the shared marker names in
#'   \code{tx_genotypes} and \code{trait_genotypes}.
#' @param LD_reduction a logical scalar or a function to perform LD
#'   reduction. If true, then the ZZZ package will be used on the
#'   combined genotypes of the two data sets for the specified
#'   \code{markers} to generate a subset of \code{markers}, which will
#'   be used in the analysis.
#' @param tx_id the index or name of the column containing the subject
#'   ID in the \code{transcriptome} and \code{tx_genotype} data
#'   frames.
#' @param trait_id the index or name of the column containing the
#'   subject ID in the \code{traits} and \code{trait_genotype} data
#'   frames.
#' @param max_tx an integer representing the maximum possible
#'   transcript count value to model. The default is the maximum
#'   observed value in \code{transcriptome}. Large values increase the
#'   modeling time and some observed counts are believed to be
#'   outliers. If no set and the observed count is greater than 1000,
#'   then a warning is issues.
#' @param rounds an integer for the maximum number of EM
#'   steps. Default is 1000.
#' @param epsilon a numeric for the EM stopping condition. If the
#'   absolute difference is less than \code{epsilon} then modeling
#'   halts. Default is 1e-6.
#' @param tiny_weights a numeric. Any weights that are less than or
#'   equal to \code{tiny_weights} are excluded from modeling. Default
#'   is 1e-16.
twas <- function(transcriptome, traits, tx_genotypes,
                 trait_genotypes,
                 markers=intersect(colnames(transcriptome),colnames(traits)),
                 LD_reduction=FALSE,
                 tx_id=1,
                 trait_id=1,
                 max_tx,
                 rounds=1000,
                 epsilon=1e-6,
                 tiny_weights=1e-16) {
}

