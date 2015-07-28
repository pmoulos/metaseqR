\name{fisher.method}
\alias{fisher.method}
\title{
Perform Fisher's Method for combining p-values
}
\description{
  Function for combining p-values by performing Fisher's method.
  The approach as described by Fisher's combines p-values to a statistic
  \deqn{S = -2 \sum^k \log p}, which follows a \eqn{\chi^2} distribution
  with \eqn{2k} degrees of freedom.
}
\usage{
fisher.method(pvals, method = c("fisher"), p.corr = c("bonferroni",
"BH", "none"), zero.sub = 1e-05, na.rm = FALSE, mc.cores=NULL)
}
\arguments{
  \item{pvals}{
    A matrix or data.frame containing the p-values from the single tests
}
  \item{method}{
    A string indicating how to combine the p-values for deriving a
    sumary p-value. Currently only the classical approach described by
    Fisher is implemented.
}
\item{p.corr}{
  Method for correcting the summary p-values. BH: Benjamini-Hochberg
  (default); Bonferroni's method or no ('none') correction are currently
  supported.
}
  \item{zero.sub}{
    Replacement for p-values of 0
}
  \item{na.rm}{
    A flag indicating whether NA values should be removed from the analysis.
}
  \item{mc.cores}{
    Currently ignored
}
}
\details{
  As \code{log(0)} results in \code{Inf} we replace p-values of 0 by default 
  with a small float. If you want to keep them as 0 you have to provide 0
  as a parameter in \code{zero.sub}.

  Note that only p-values between 0 and 1 are allowed to be passed to this 
  method.
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
  Fisher, R.A. (1925). Statistical Methods for Research Workers. Oliver and 
  Boyd (Edinburgh).
  Moreau, Y.et al. (2003). Comparison and meta-analysis of microarray data: 
  from the bench to the computer desk. Trends in Genetics, 19(10), 570-577.
}
\author{
Karl Kugler <karl@eigenlab.net>
}
\seealso{
\code{\link{fisher.method.perm}}
}
\examples{
set.seed(123)
pp <- matrix(c(runif(20),c(0.001,0.02,0.03,0.001)), ncol=4)
pp[2,3] <- NA
fisher.method(pp) #returns one NA row
fisher.method(pp, na.rm=TRUE) #ignore NA entry in that row
}
