\name{fisher.method.perm}
\alias{fisher.method.perm}

\title{
Derive a p-value for a summary statistic of p-values by permutation
}
\description{
  Given a set of p-values and a summary statistic S: \deqn{S = -2 \sum
    \log p,} a p-value for this statistic can be derived by randomly
  generating summary statistics [Rhodes,2002]. Therefore, a p-value is
  randomly sampled from each contributing study and a random statistic is 
  calculated. The fraction of random statistics greater or equal to S 
  then gives the p-value.
}
\usage{
  fisher.method.perm(pvals, p.corr = c("bonferroni", "BH", "none"), 
  zero.sub = 1e-05, B = 10000, mc.cores = NULL, blinker = 1000)
}
\arguments{
 \item{pvals}{
  A matrix or data.frame containing the p-values from the single tests
}
  \item{p.corr}{
  Method for correcting the summary p-values. BH: Benjamini-Hochberg
  (default); Bonferroni's method or no ('none') correction are currently
  supported.
}
  \item{zero.sub}{
  Replacement for p-values of 0
}
  \item{B}{
Number of random statistics
}
  \item{mc.cores}{
Number of cores used for calculating the permutations. If not
\code{NULL} the \code{multicore} package is used for parallelization
with the given number of cores.
}
\item{blinker}{
  An indicator that prints "\code{=}" after each \code{blinker} rows of
  \code{pvals} in order to follow the progress.
}
}
\details{
  At the moment this function only supports situations were all passed 
  p-values are  not \code{NA}. We plan on extending this functionality 
  in upcoming versions.
  
  For large data sets and/or large \code{B} we strongly recommend using
  the \code{mc.cores} option as the calculation will otherwise be
  computationally demanding. This will call the \code{mclapply}
  function from the \cite{multicore} package, which you will have to
  install in that case.

  By default a blinker (a small string "=") is shown after each 1000
  rows that were computed. This function allows you to assess the
  progress of the analysis. If you don't want to see the blinker set 
  it to \code{NA}.

  As \code{log(0)} results in \code{Inf} we replace p-values of 0 by 
  default with a small float. If you want to keep them as 0 you have 
  to provide 0 as a parameter in \code{zero.sub}.

  Note that only p-values between 0 and 1 are allowed to be passed to 
  this method. 
}
\value{
  This method returns a data.frame containing the following columns
 \item{S }{The statistic}
 \item{num.p }{The number of p-values used to calculate S}
 \item{p.value }{The overall p-value}
 \item{p.adj}{The adjusted p-value}
}
\note{
  This function was copied from the CRAN package MADAM which is no longer
  maintained. Recognition goes to the original author(s) below.
}
\references{
Rhodes, D. R., (2002). Meta-analysis of microarrays: interstudy 
alidation of gene expression profiles reveals pathway dysregulation 
in prostate cancer. Cancer research, 62(15), 4427-33.
}
\author{
Karl Kugler <karl@eigenlab.net>
}
\seealso{
\code{\link{fisher.sum}}, \code{\link{fisher.method}}
}
\examples{
  set.seed(123)
  pp <- matrix(c(runif(20),c(0.001,0.02,0.03,0.001)), ncol=4)
  fisher.method.perm(pp, B=10, blinker=1)
  \dontrun{
    fisher.method.perm(pp, B=10000, mc.cores=3, blinker=1) #use multicore
  }
}
