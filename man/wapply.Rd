\name{wapply}
\alias{wapply}
\title{List apply helper}
\usage{
    wapply(m, ...)
}
\arguments{
    \item{m}{a logical indicating whether to execute in
    parallel or not.}

    \item{...}{the rest arguments to \code{\link{lapply}} (or
    \code{mclapply})}
}
\description{
    A wrapper around normal and parallel apply
    (\code{\link{mclapply}} or multicore package) to avoid
    excessive coding for control of single or parallel code
    execution. Internal use.
}
\examples{
\donttest{
multic <- check.parallel(0.8)
# Test meaningful only in machines where parallel computation 
# supported
if (multic) {
    system.time(r<-wapply(TRUE,1:10,function(x) runif(1e+6)))
    system.time(r<-wapply(FALSE,1:10,function(x) runif(1e+6)))
}
}
}
\author{
    Panagiotis Moulos
}

