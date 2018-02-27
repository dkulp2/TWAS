# calculating posterior probabilities - for trait study only
#
# For each person with trait, create tx.max (e.g. 1000) posterior probs corresponding to prob that RNAseq = 1 ... 1000
# See equation 6: Ri is ri in equation. y1[i] is yi in equation.

posterior.prob <- function(geno, y1, Ri, beta.hat, gamma.hat, sige.hat, phi.hat, eta.hat) {
  pii <- matrix(nrow=length(Ri), ncol=length(y1))
  for (i in 1:ncol(pii)) {
    pii[,i] <-
      dnorm(y1[i],
            mean=beta.hat + Ri*gamma.hat,
            sd=sige.hat,log=FALSE) *
      dnbinom(Ri,
              size=phi.hat,
              mu=exp(eta.hat[1]+sum(geno[i,]*eta.hat[-1])), log = FALSE)
    pii[,i] <- pii[,i] / sum(pii[,i])
  }
  return(pii)
}

#  calculating negative binomial likelihood
logLik3 <- function(phi.c, r.wAll, mu, weightsAll) {
  dens <- log(dnbinom(r.wAll, size=phi.c, mu=mu, log = FALSE))
  return(sum(weightsAll*dens,na.rm=T))
}

# perform LD reduction on matrix m using SNPRelate, returning SNP IDs.
# Assumes chromosome positions are sequential in matrix.
def.LDprune <- function(m) {
  f <- tempfile(fileext=".gds")
  SNPRelate::snpgdsCreateGeno(f, genmat = m,
    sample.id = rownames(m), snp.id = colnames(m),
    snp.chromosome = rep(0, ncol(m)),
    snp.position = seq(1, ncol(m)), snpfirstdim=FALSE)

  genofile <- snpgdsOpen(f)
  unlist(SNPRelate::snpgdsLDpruning(genofile))
}

# returns either a 'progress_bar' object, a character object, or NULL
pbar.create <- function(show.progress=TRUE) {
  if (show.progress) {
    if (requireNamespace("progress", quietly = TRUE)) {
      pbar <- progress::progress_bar$new(
        format = "k=:k eps=:eps [:bar] :elapsed"
      )
    } else {
      pbar <- 'text'
    }
  } else {
    pbar <- NULL
  }
}

# update the progress bar with the current iteration and epsilon
pbar.update <- function(pbar, k, abs.diff, epsilon) {
  if (!is.null(pbar)) {
    if (is(pbar, 'progress_bar')) {
      est.frac.done <- max(0, min(1,log(abs.diff)/log(epsilon)))  # bar shows exponent
      pbar$update(est.frac.done, tokens=list(k=k, eps=sprintf("%.03e",abs.diff)))
    } else {
      cat("k=",k," of ",rounds,"\tepsilon: ",abs.diff," > ",epsilon,"\n")
    }
  }
  invisible()
}

# finalize the progress bar
pbar.done <- function(pbar, k, abs.diff) {
  if (!is.null(pbar)) {
    if (is(pbar, 'progress_bar')) {
      pbar$update(1, tokens=list(k=k, eps=sprintf("%.03e",abs.diff)))
    } else {
      cat("Done.\n")
    }
  }
  invisible()
}

# always warn immediately without code ref
warning <- function(...) {
  base::warning(..., call.=FALSE, immediate.=TRUE)
}

#
#' An EM model fit jointly to transcriptome and trait data
#'
#' \code{twas} returns a model fit from two data sets from different
#' subjects: one set of transcriptome data and the other of trait
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
#' @param tx.genotypes a data frame of genotypes for the transcriptome
#'   subjects. Each row is a subject. There is one ID column with the
#'   same name as used in \code{transcriptome} and the remaining
#'   columns are markers. Genotype values are non-negative integers.
#'   Subjects with any NA genotypes are ignored.
#' @param trait.genotypes a data frame of genotypes for the trait
#'   subjects. Each row is a subject. There is one ID column with the
#'   same name as used in \code{traits} and the remaining columns are
#'   markers. Genotype values are non-negative integers. Subjects with
#'   any NA genotypes are ignored.
#' @param gene a character identifying the column in
#'   \code{transcriptome} of the gene of interest. Default is the
#'   first column in \code{transcriptome} excluding the subject ID.
#' @param trait.names a character vector of one or more column names
#'   in \code{traits} to be used in analysis. Defaults to all columns
#'   in \code{traits} except the subject ID.
#' @param markers a character vector of the names of the markers to
#'   use for analysis. Defaults to the shared marker names in
#'   \code{tx.genotypes} and \code{trait.genotypes}.
#' @param LD.reduction a logical scalar or a function to perform LD
#'   reduction. If true, then the SNPRelate package will be used on the
#'   genotypes of the two data sets for the specified
#'   \code{markers} to generate a subset of \code{markers}, which will
#'   be used in the analysis. If a function is supplied, then it must
#'   accept a matrix of columns of genotypes and return a vector of
#'   column names of the SNPs to keep. Default is false.
#' @param tx.id.col the index or name of the column containing the
#'   subject ID in the \code{transcriptome} and \code{tx_genotype}
#'   data frames. Defaults to the first column in
#'   \code{transcriptome}.
#' @param trait.id.col the index or name of the column containing the
#'   subject ID in the \code{traits} and \code{trait_genotype} data
#'   frames. Defaults to the first column in \code{traits}.
#' @param tx.max an integer representing the maximum possible
#'   transcript count value to model. The default is the maximum
#'   observed value in \code{transcriptome}. Large values increase the
#'   modeling time and some observed counts are believed to be
#'   outliers. If not set and the observed count is greater than 1000,
#'   then a warning is issues.
#' @param rounds an integer for the maximum number of EM
#'   steps. Default is 1000.
#' @param epsilon a numeric for the EM stopping condition. If the
#'   absolute difference is less than \code{epsilon} then modeling
#'   halts. Default is 1e-6.
#' @param tiny.weights a numeric. Any weights that are less than or
#'   equal to \code{tiny.weights} are excluded from modeling. Default
#'   is 1e-16.
#' @param show.progress a logical indicating whether to display progress
#'   during EM iterations. If the \code{progress} package is installed
#'   the it will be used to display a progress bar. Otherwise, text
#'   will be displayed. Default is true.
#' @return A list object containing the fit expression levels and
#'   trait values for all subjects in a single data frame (data), each
#'   model parameter (gamma, beta, eta), and the variance and p-value
#'   of the test statistic for each parameter.
#'
#' @examples
#' Not run:
#' twas(twas.tx$transcriptome, twas.traits$traits,
#'      twas.tx$genotypes, twas.traits$genotypes,
#'      gene="LRRC16A", trait.names="logIL6",
#'      markers=twas.gene.markers$LRRC16A)
#'
#' @seealso \code{\link{twas.tx}}, \code{\link{twas.traits}} and
#'   \code{\link{twas.gene.markers}} for examples of the necessary
#'   data formats.
#'
#' @importFrom assertthat assert_that
#' @export
twas <- function(transcriptome,
                 traits,
                 tx.genotypes,
                 trait.genotypes,
                 gene,
                 trait.names,
                 markers,
                 LD.reduction=FALSE,
                 tx.id.col=1,
                 trait.id.col=1,
                 tx.max,
                 rounds=1000,
                 epsilon=1e-6,
                 tiny.weights=1e-16,
                 show.progress=TRUE)
{
  # Minimal argument checking
  assert_that(assertthat::is.flag(show.progress))
  assert_that(is.data.frame(transcriptome))
  assert_that(is.data.frame(traits))
  assert_that(is.data.frame(tx.genotypes))
  assert_that(is.data.frame(trait.genotypes))
  assert_that(assertthat::is.flag(LD.reduction) || is.function(LD.reduction))
  rounds <- as.integer(rounds)
  assert_that(rounds > 0)
  assert_that(is.numeric(epsilon) && epsilon < 1)
  assert_that(is.numeric(tiny.weights) && tiny.weights < 1)

  # Use column names, not indices
  if (is.numeric(tx.id.col)) tx.id.col <- colnames(transcriptome)[tx.id.col]
  if (is.numeric(trait.id.col)) trait.id.col <- colnames(traits)[trait.id.col]
  assert_that(tx.id.col %in% colnames(transcriptome))
  assert_that(trait.id.col %in% colnames(traits))

  # Housekeeping --------------------

  # Set *.ids to the identifier column, arrange the genotypes to be the same
  # order and nrow as transcriptome/trait, remove the id
  # column leaving just expression, trait and genotype
  tx.ids <- dplyr::select(transcriptome, dplyr::one_of(tx.id.col))
  transcriptome <- dplyr::select(transcriptome, -dplyr::one_of(tx.id.col))
  tx.genotypes <- dplyr::select(dplyr::left_join(tx.ids, tx.genotypes,by=tx.id.col), -dplyr::one_of(tx.id.col))

  fx.ids <- dplyr::select(traits, dplyr::one_of(trait.id.col))
  traits <- dplyr::select(traits, -dplyr::one_of(trait.id.col))
  trait.genotypes <- dplyr::select(dplyr::left_join(fx.ids, trait.genotypes,by=tx.id.col), -dplyr::one_of(trait.id.col))

  # set markers to common markers in two files
  if (missing(markers)) {
    markers <- intersect(colnames(tx.genotypes), colnames(trait.genotypes))
  }
  assert_that(is.character(markers)) # 1 or more

  # Optionally perform LD reduction
  if (assertthat::is.flag(LD.reduction) && LD.reduction) {
    # use default LD function
    if (!requireNamespace("SNPRelate", quietly = TRUE)) {
      stop("SNPRelate needed for LD reduction. Please install it with source('http://bioconductor.org/biocLite.R'); biocLite('SNPRelate')",
           call. = FALSE)
    }
    LD.reduction <- def.LDprune
  }
  if (is.function(LD.reduction)) {
      markers <- LD.reduction(as.matrix(rbind(tx.genotypes[,markers], trait.genotypes[,markers])))
  }

  # use only selected marker columns from _genotypes
  tx.genotypes <- dplyr::select(tx.genotypes, dplyr::one_of(markers))
  trait.genotypes <- dplyr::select(trait.genotypes, dplyr::one_of(markers))

  # filter to just selected gene
  if (missing(gene)) {
    gene <- colnames(transcriptome)[1]
  }
  assert_that(assertthat::is.string(gene))  # only 1

  transcriptome <- dplyr::as_tibble(dplyr::select(transcriptome, dplyr::one_of(gene)))

  # filter to just select traits
  if (missing(trait.names)) {
    trait.names <- colnames(traits)
  }
  traits <- dplyr::as_tibble(dplyr::select(traits, dplyr::one_of(trait.names)))

  # Check that extraction of ID columns and filtering to specified
  # genes, markers and traits results in valid data.frames
  assert_that(ncol(tx.ids)==1)
  assert_that(ncol(fx.ids)==1)
  assert_that(ncol(transcriptome)>=1)
  assert_that(ncol(traits)>=1)
  assert_that(ncol(tx.genotypes)>=1)
  assert_that(ncol(trait.genotypes)>=1)
  assert_that(nrow(transcriptome)==nrow(tx.genotypes))
  assert_that(nrow(traits)==nrow(trait.genotypes))

  # RNA values must be non-negative integers
  assert_that(all(sapply(transcriptome, function(col) is.numeric(col) && all(col >= 0))))
  if (any(sapply(transcriptome, function(col) !is.integer(col)))) {
    warning("Coercing real valued transcriptome to integers")
    transcriptome <- data.frame(lapply(transcriptome, function(col) as.integer(col)))
  }

  # trait values can be any numeric
  assert_that(all(sapply(traits, function(col) is.numeric(col))))

  if (missing(tx.max)) {
    tx.max <- max(transcriptome[[gene]])
    if (tx.max > 1000) { warning(paste0("maximum expression, ", tx.max, ", is greater than 1000. May be outlier. Consider setting tx.max parameter.")) }
  }
  assert_that(is.numeric(tx.max))

  # eliminate subjects with NA genotype
  tx.genotypes <- dplyr::as_tibble(na.omit(tx.genotypes))
  if (!is.null(attributes(tx.genotypes)$na.action)) {  # na.omit stores dropped row indices as attribute
    transcriptome <- transcriptome[-attributes(tx.genotypes)$na.action,]
    tx.omit <- tx.ids[attributes(tx.genotypes)$na.action,1]
  } else {
    tx.omit <- c()
  }


  trait.genotypes <- dplyr::as_tibble(na.omit(trait.genotypes))
  if (!is.null(attributes(trait.genotypes)$na.action)) {
    traits <- traits[-attributes(trait.genotypes)$na.action,]
    fx.omit <- fx.ids[attributes(trait.genotypes)$na.action,1]
  } else {
    fx.omit <- c()
  }

  # report any dropped samples as a warning.
  if (length(tx.omit) > 0) warning(paste0(length(tx.omit)," subjects removed from transcriptome set due to NA genotypes values. (", vec2str(as.character(tx.omit)), ")"))
  if (length(fx.omit) > 0) warning(paste0(length(fx.omit)," subjects removed from trait set due to NA genotype values. (", vec2str(as.character(fx.omit)), ")"))

  ## This is the end of the data clean.
  ## From here on, we just deal with tx and fx as source data

  # tx: transcript X genotypes
  tx=dplyr::bind_cols(transcriptome, tx.genotypes)

  # fx: features (aka traits) X genotypes
  fx=dplyr::bind_cols(traits, trait.genotypes)

  tx.formula <- stats::reformulate(response=gene, termlabels='.')  # gene~.
  mod.tx.nb <- MASS::glm.nb(tx.formula, data=tx)

  fx.formula <- as.formula(paste(paste(trait.names, collapse='+'),'~1')) # trait.names~1
  mod.fx <- stats::lm(fx.formula, data=fx)

  # FIXME: Original does not allow multiple traits. Can y1 be a matrix here?
  y1 <- fx[[trait.names]]

  r2 <- tx[[gene]]

  n.tx <- nrow(tx)  # number of transcriptome subjects
  n.fx <- nrow(fx)  # number of trait subjects
  n.snps <- length(markers)  # number of SNPs

  # FIXME
  trait.geno <- as.matrix(dplyr::select(fx, -dplyr::one_of(trait.names)))
  tx.geno <- as.matrix(dplyr::select(tx, -dplyr::one_of(gene)))

  # model 1 (eq 2) (just fx data to start)
  # mod.fx, above, is initial model that just fits the interecept
  updated.beta.hat <- mod.fx$coefficients[1] # intercept parameter. FIXME: multiple traits
  updated.gamma.hat <- 0
  updated.sige.hat <- summary(mod.fx)$sigma  # std dev from initial estimate

  # model 2 (eq 3) (just tx data to start)
  updated.phi.hat <- mod.tx.nb$theta
  updated.eta.hat <- as.matrix(mod.tx.nb$coefficients)

  # ri in Eq 6
  Ri <- seq(0, tx.max-1)

  # expanding y and r
  y.w <- kronecker(y1,rep(1, tx.max)) # expanded version of fx
  r.w <- rep(Ri, n.fx)
  r.wAll <- c(r.w, r2)                # combine tx and fx

  # same idea expanding genotypes
  geno.w <- kronecker(trait.geno, rep(1, tx.max))
  geno.wAll <- rbind(geno.w, tx.geno) # hmmm...

  k=1
  abs.diff=1

  # display a progress bar, if available
  pbar <- pbar.create(show.progress)
  pbar.update(pbar, k, abs.diff, epsilon)

  while (k<=rounds && abs.diff>epsilon){

    updated.beta.hat.orig <- updated.beta.hat
    updated.gamma.hat.orig <- updated.gamma.hat
    updated.sige.hat.orig <- updated.sige.hat
    updated.phi.hat.orig <- updated.phi.hat
    updated.eta.hat.orig <- updated.eta.hat

    # estimating mean components
    # model 1?
    pii <- posterior.prob(trait.geno, y1, Ri,
                          updated.beta.hat, updated.gamma.hat, updated.sige.hat,
                          updated.phi.hat, updated.eta.hat)
    weights.1 <- as.vector(pii)             # fx only
    weightsAll <- c(weights.1,rep(1,n.tx))  # because we know the RNAseq value for tx, then set weights to 1 for all of those individuals

    # see practical implications. eta needs to be estimated with a modified model
    significant.weights <- weights.1>tiny.weights
    mod2w <- lm(y.w[significant.weights]~r.w[significant.weights],weights=weights.1[significant.weights])
    updated.beta.hat <- mod2w$coefficients[1]
    updated.gamma.hat <- mod2w$coefficients[2]

    significant.weightsAll <- weightsAll>tiny.weights
    mod1w <- MASS::glm.nb(r.wAll[significant.weightsAll]~geno.wAll[significant.weightsAll,],weights=weightsAll[significant.weightsAll])
    updated.eta.hat <- mod1w$coefficients

    # estimating variance components conditional on mean components
    pii <- posterior.prob(trait.geno, y1, Ri,
                          updated.beta.hat, updated.gamma.hat, updated.sige.hat,
                          updated.phi.hat, updated.eta.hat)
    weights.2 <- as.vector(pii)
    weightsAll <- c(weights.2, rep(1,n.tx))

    updated.sige.hat <- sqrt(sum(weights.2[weights.1>1e-16]*(mod2w$residuals)^2)/sum(weights.2[weights.1>1e-16]))*sqrt(n.fx/(n.fx-5))
    mu <- exp(updated.eta.hat[1]+(geno.wAll%*%updated.eta.hat[-1]))

    updated.phi.hat <- optimize(logLik3, interval=c(0,8), maximum=TRUE, r.wAll, mu, weightsAll)$maximum

    # determining if convergence is met
    diff <- c(updated.beta.hat.orig-updated.beta.hat,
              updated.gamma.hat.orig-updated.gamma.hat,
              updated.sige.hat.orig-updated.sige.hat,
              updated.phi.hat.orig-updated.phi.hat,
              updated.eta.hat.orig-updated.eta.hat)
    abs.diff <- max(abs(diff))

    pbar.update(pbar, k, abs.diff, epsilon)
    k <- k+1
  }
  pbar.done(pbar, k, abs.diff)

  im <- info.matrix(y.w=y.w, r.w=r.w, r.wAll=r.wAll,
                    weights=weights.2, weights.all=weightsAll,
                    geno.w=geno.w, geno.wAll=geno.wAll,
                    n.fx=n.fx, n.tx=n.tx, n.snps=n.snps, tx.max=tx.max,
                    beta.hat=updated.beta.hat, gamma.hat=updated.gamma.hat,
                    sige.hat=updated.sige.hat, phi.hat=updated.phi.hat, eta.hat=updated.eta.hat)

  sqrt(diag(solve(im))) # solve(): inverse of Info
  eigen(solve(im))$values

  ### testing parameters
  Var <- diag(solve(im))
  gamma.var <- Var[1]
  beta.var <- Var[2]
  eta.var <- Var[4:(4+n.snps)]

  # analyze each parameter, output each of these as parameter and its corresponding Variance, and test statistic (p-value of test that gamma is zero)
  gamma.pval <- (1-pnorm(abs(updated.gamma.hat/sqrt(gamma.var))))*2
  beta.pval <- (1-pnorm(abs(updated.beta.hat/sqrt(beta.var))))*2

  # eta is a vector, so will be different sizes. this is interactive. choose one.
  eta.pval <- (1-pnorm(abs(updated.eta.hat/sqrt(eta.var))))*2

  return(list(beta=unname(updated.beta.hat), beta.var=beta.var, beta.pval=unname(beta.pval),
              gamma=unname(updated.gamma.hat), gamma.var=gamma.var, gamma.pval=unname(gamma.pval),
              sige=updated.sige.hat, phi=updated.phi.hat,
              eta=unname(updated.eta.hat), eta.var=eta.var, eta.pval=unname(eta.pval)))


}

