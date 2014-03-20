\name{combine.bonferroni}
\alias{combine.bonferroni}
\title{Combine p-values with Bonferroni's method}
\usage{
    combine.bonferroni(p)
}
\arguments{
    \item{p}{a p-value matrix (rows are genes, 
    columns are statistical tests).}
}
\value{
    A vector of combined p-values. 
}
\description{
    This function combines p-values from the 
    various statistical tests supported by 
    metaseqR using the Bonferroni's method (see 
    reference in the main \code{\link{metaseqr}}
    help page or in the vignette).
}
\examples{
p <- matrix(runif(300),100,3)
pc <- combine.bonferroni(p)
}
\author{
    Panagiotis Moulos
}

