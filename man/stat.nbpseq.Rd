\name{stat.nbpseq}
\alias{stat.nbpseq}
\title{Statistical testing with NBPSeq}
\usage{
    stat.nbpseq(object, sample.list, contrast.list = NULL,
        stat.args = NULL, libsize.list = NULL)
}
\arguments{
    \item{object}{a matrix or an object specific to each
    normalization algorithm supported by metaseqR, containing
    normalized counts. Apart from matrix (also for NOISeq),
    the object can be a SeqExpressionSet (EDASeq),
    CountDataSet (DESeq), DGEList (edgeR) or list (NBPSeq).}

    \item{sample.list}{the list containing condition names
    and the samples under each condition.}

    \item{contrast.list}{a named structured list of contrasts
    as returned by \code{\link{make.contrast.list}} or just
    the vector of contrasts as defined in the main help page
    of \code{\link{metaseqr}}.}

    \item{stat.args}{a list of NBPSeq statistical algorithm
    parameters. See the result of
    \code{get.defaults("statistics",} \code{"nbpseq")} for an
    example and how you can modify it. It is not required
    when the input object is already a list from NBPSeq
    normalization as the dispersions are already estimated.}

    \item{libsize.list}{an optional named list where names
    represent samples (MUST be the same as the samples
    \code{in sample.list}) and members are the library sizes
    (the sequencing depth) for each sample. If not provided,
    the default is the column sums of the \code{gene.counts}
    matrix.}
}
\value{
    A named list of p-values, whose names are the names of
    the contrasts.
}
\description{
    This function is a wrapper over NBPSeq statistical
    testing. It accepts a matrix of normalized gene counts or
    an S4 object specific to each normalization algorithm
    supported by metaseqR.
}
\note{
    There is currently a problem with the NBPSeq package and
    the workflow that is specific to the NBPSeq package. The
    problem has to do with function exporting as there are
    certain functions which are not recognized from the
    package internally. For this reason and until it is
    fixed, only the Smyth workflow will be available with the
    NBPSeq package.
}
\examples{
\donttest{
require(DESeq)
data.matrix <- counts(makeExampleCountDataSet())
sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
norm.data.matrix <- normalize.nbpseq(data.matrix,sample.list)
p <- stat.nbpseq(norm.data.matrix,sample.list,contrast)
}
}
\author{
    Panagiotis Moulos
}

