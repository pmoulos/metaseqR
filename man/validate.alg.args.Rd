\name{validate.alg.args}
\alias{validate.alg.args}
\title{Validate normalization and statistical 
    algorithm arguments}
\usage{
    validate.alg.args(normalization, statistics, 
        norm.args, stat.args)
}
\arguments{
    \item{normalization}{a keyword determining the 
    normalization strategy to be performed by 
    metaseqR. See \code{\link{metaseqr}} main help 
    page for details.}

    \item{statistics}{the statistical tests to be 
    performed by metaseqR. See \code{\link{metaseqr}} 
    main help page for details.}
    
    \item{norm.args}{the user input list of 
    normalization arguments. See 
    \code{\link{metaseqr}} main help page for 
    details.}

    \item{stat.args}{the user input list of 
    statistical test arguments. See 
    \code{\link{metaseqr}} main help page for 
    details.}
}
\value{
    A list with two members (\code{norm.args}, 
    \code{stat.args}) with valid arguments to 
    be used as user input for the algorithms 
    supported by metaseqR.
}
\description{
    This function checks and validates the arguments 
    passed by the user to the normalization and 
    statistics algorithms supported by metaseqR. As 
    these are given into lists and passed to the 
    algorithms, the list members must be checked
    for \code{NULL}, valid names etc. This function 
    performs these checks and ignores any invalid 
    arguments.
}
\examples{
\donttest{
normalization <- "edaseq"
statistics <- "edger"
norm.args <- get.defaults("normalization","edaseq")
stat.args <- get.defaults("statistics","deseq")
# Will return as is
val <- validate.alg.args(normalization,statistics,norm.args,stat.args)
val$norm.args
val$stat.args
# but...
stat.args <- c(stat.args,my.irrelevant.arg=999)
val <- validate.alg.args(normalization,statistics,norm.args,stat.args)
# irrelevant argument will be removed
val$norm.args
val$stat.args
}
}
\author{
    Panagiotis Moulos
}

