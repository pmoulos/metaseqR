\name{diagplot.noiseq}
\alias{diagplot.noiseq}
\title{Diagnostic plots based on the NOISeq package}
\usage{
    diagplot.noiseq(x, sample.list, covars,
        which.plot = c("biodetection", "countsbio", "saturation", "rnacomp", 
                "biodist"),
        output = "x11",
        biodist.opts = list(p = NULL, pcut = NULL, name = NULL),
        path = NULL, is.norm = FALSE, ...)
}
\arguments{
    \item{x}{the count data matrix.}

    \item{sample.list}{the list containing condition names
    and the samples under each condition.}

    \item{covars}{a list (whose annotation elements are
    ideally a subset of an annotation data frame produced by
    \code{\link{get.annotation}}) with the following members:
    data (the data matrix), length (gene length), gc (the
    gene gc_content), chromosome (a data frame with
    chromosome name and co-ordinates), factors (a factor with
    the experimental condition names replicated by the number
    of samples in each experimental condition) and biotype
    (each gene's biotype as depicted in Ensembl-like
    annotations).}

    \item{which.plot}{the NOISeq package plot to generate. It
    can be one or more of \code{"biodetection"},
    \code{"countsbio"}, \code{"saturation"},
    \code{"rnacomp"}, \code{"readnoise"} or \code{"biodist"}.
    Please refer to the documentation of the EDASeq package
    for details on the use of these plots. The
    \code{which.plot="saturation"} case is modified to be
    more informative by producing two kinds of plots. See
    \code{\link{diagplot.noiseq.saturation}}.}

    \item{biodist.opts}{a list with the following members: p
    (a vector of p-values, e.g. the p-values of a contrast),
    pcut (a unique number depicting a p-value cutoff,
    required for the \code{"biodist"} case), name (a name for
    the \code{"biodist"} plot, e.g. the name of the
    contrast.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"x11"}
    (default), \code{"png"}, \code{"jpg"}, \code{"bmp"},
    \code{"pdf"} or \code{"ps"}.}

    \item{path}{the path to create output files.}

    \item{is.norm}{a logical indicating whether object
    contains raw or normalized data. It is not essential and
    it serves only plot annotation purposes.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    The filenames of the plots produced in a named list with
    names the \code{which.plot} argument. If
    \code{output="x11"}, no output filenames are produced.
}
\description{
    A wrapper around the plotting functions availale in the
    NOISeq Bioconductor package. For analytical explanation
    of each plot please see the vignette of the NOISeq 
    package. It is best to use this function through the 
    main plotting function \code{\link{diagplot.metaseqr}}.
}
\note{
    Please note that in case of \code{"biodist"} plots, the
    behavior of the function is unstable, mostly due to the
    very specific inputs this plotting function accepts in
    the NOISeq package. We have tried to predict unstable
    behavior and avoid exceptions through the use of tryCatch
    but it's still possible that you might run onto an error.
}
\examples{
\donttest{
require(DESeq)
data.matrix <- counts(makeExampleCountDataSet())
sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
lengths <- round(1000*runif(nrow(data.matrix)))
starts <- round(1000*runif(nrow(data.matrix)))
ends <- starts + lengths
covars <- list(
    data=data.matrix,
    length=lengths,
    gc=runif(nrow(data.matrix)),
    chromosome=data.frame(
        chromosome=c(rep("chr1",nrow(data.matrix)/2),
            rep("chr2",nrow(data.matrix)/2)),
        start=starts,
        end=ends
    ),
    factors=data.frame(class=as.class.vector(sample.list)),
    biotype=c(rep("protein_coding",nrow(data.matrix)/2),rep("ncRNA",
            nrow(data.matrix)/2))
)
p <- runif(nrow(data.matrix))
diagplot.noiseq(data.matrix,sample.list,covars=covars,
        biodist.opts=list(p=p,pcut=0.1,name="A_vs_B"))
}
}
\author{
    Panagiotis Moulos
}

