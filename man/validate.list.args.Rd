\name{validate.list.args}
\alias{validate.list.args}
\title{Validate list parameters for several metaseqR functions}
\usage{
    validate.list.args(what, method = NULL, arg.list)
}
\arguments{
    \item{what}{what a keyword determining the procedure 
    for which to validate arguments. It can be one of 
    \code{"normalization"}, \code{"statistics"}, 
    \code{"gene.filter"}, \code{"exon.filter"} or 
    \code{"biotype.filter"}.}

    \item{method}{the normalization/statistics/filtering
    algorithm included in metaseqR for which to validate 
    user input. When \code{what} is \code{"normalization"}, 
    method is one of \code{"edaseq"}, \code{"deseq"}, 
    \code{"edger"}, \code{"noiseq"} or \code{"nbpseq"}. 
    When \code{what} is \code{"statistics"}, method is 
    one of \code{"deseq"}, \code{"edger"}, \code{"noiseq"}, 
    \code{"bayseq"}, \code{"limma"} or \code{"nbpseq"}. 
    When \code{method} is \code{"biotype.filter"}, 
    \code{what} is the input organism (see the main 
    \code{\link{metaseqr}} help page for a list of 
    supported organisms).}
    
    \item{arg.list}{the user input list of arguments.}
}
\value{
    A list with valid arguments to be used as user 
    input in the algorithms supported by metaseqR.
}
\description{
    This function validates the arguments passed 
    by the user to the normalization, statistics 
    and filtering algorithms supported by metaseqR. 
    As these are given into lists and passed to 
    the algorithms, the list member names must 
    be valid algorithm arguments for the pipeline 
    not to crash. This function performs these 
    checks and ignores any invalid arguments.
}
\examples{
\donttest{
norm.args.edger <- list(method="TMM",refColumn=NULL,
    logratioTrim=0.3,sumTrim=0.05,doWeighting=TRUE,
    Bcutoff=-1e10,p=0.75)
# Bcutoff does not exist, will throw a warning and ignore it.
norm.args.edger <- validate.list.args("normalization",
    "edger",norm.args.edger)
}
}
\author{
    Panagiotis Moulos
}
