\name{make.permutation}
\alias{make.permutation}
\title{Create counts matrix permutations}
\usage{
    make.permutation(counts, sample.list, contrast,
        repl = FALSE)
}
\arguments{
    \item{counts}{the gene read counts matrix.}

    \item{sample.list}{the list containing condition names
    and the samples under each condition.}

    \item{contrast}{the contrasts vector. See the main
    \code{\link{metaseqr}} help page.}

    \item{repl}{the same as the replace argument in
    \code{\link{sample}} function.}
}
\value{
    A list with three members: the matrix of permuted per
    sample read counts, the virtual sample list and the
    virtual contrast to be used with the \code{stat.*}
    functions.
}
\description{
    This function creates a permuted read counts matrix based
    on the \code{contrast} argument (to define new virtual
    contrasts of the same number) and on the
    \code{sample.list} to derive the number of samples for
    each virtual condition.It is a helper for the
    \code{\link{meta.perm}} function.
}
\examples{
\donttest{
data("mm9.gene.data",package="metaseqR")
per <- make.permutation(mm9.gene.counts,sample.list.mm9,
    "e14.5_vs_adult_8_weeks")
}
}
\author{
    Panagiotis Moulos
}

