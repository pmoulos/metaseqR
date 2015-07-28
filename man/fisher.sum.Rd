\name{fisher.sum}
\alias{fisher.sum}
\title{
A function to calculate Fisher's sum for a set of p-values
}
\description{
  This method combines a set of p-values using Fisher's method:
  \deqn{S = -2 \sum \log p}
}
\usage{
fisher.sum(p, zero.sub=0.00001, na.rm=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
    A vector of p-values
  }
  \item{zero.sub}{
    Replacement for 0 values.
  }
  \item{na.rm}{
    Should NA values be removed before calculating the sum
  }
  
}
\details{
  As \code{log(0)} results in \code{Inf} we replace p-values of 0 by 
  default with a small float. If you want to keep them as 0 you have 
  to provide 0 as a parameter in \code{zero.sub}.

  Note that only p-values between 0 and 1 are allowed to be passed to 
  this method.
}
\value{
  Fisher's sum as described above.
}
\note{
  This function was copied from the CRAN package MADAM which is no longer
  maintained. Recognition goes to the original author(s) below.
}
\references{
  Fisher, R.A. (1925). Statistical Methods for Research Workers. Oliver 
  and Boyd (Edinburgh).
}
\author{
 Karl Kugler <karl@eigenlab.net>
}
\seealso{
\code{\link{fisher.method}}
}
\examples{
fisher.sum(c(0.2,0.05,0.05))
fisher.sum(c(0.2,0.05,0.05, NA), na.rm=TRUE)
}
