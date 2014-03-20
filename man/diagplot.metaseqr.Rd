\name{diagplot.metaseqr}
\alias{diagplot.metaseqr}
\title{Diagnostic plots for the metaseqR package}
\usage{
    diagplot.metaseqr(object, sample.list, annotation = NULL,
        contrast.list = NULL, p.list = NULL,
        thresholds = list(p = 0.05, f = 1),
        diagplot.type = c("mds", "biodetection", "countsbio", "saturation", 
        "readnoise", "rnacomp", "correl", "pairs", "boxplot", "gcbias", 
        "lengthbias", "meandiff", "meanvar", "deheatmap", "volcano", 
        "biodist", "filtered", "venn"),
        is.norm = FALSE, output = "x11", path = NULL, ...)
}
\arguments{
    \item{object}{a matrix or a data frame containing count
    data derived before or after the normalization procedure,
    filtered or not by the metaseqR's filters and/or p-value.
    The object can be fed to any of the
    \code{diagplot.metaseqr} plotting systems but not every
    plot is meaningful. For example, it's meaningless to
    create a \code{"biodist"} plot for a count matrix before
    normalization or statistical testing.}

    \item{sample.list}{the list containing condition names
    and the samples under each condition.}

    \item{annotation}{a data frame containing annotation
    elements for each row in object. Usually, a subset of the
    annotation obtained by \code{\link{get.annotation}} or a
    subset of possibly embedded annotation with the input
    counts table. This parameter is optional and required
    only when diagplot.type is any of \code{"biodetection"},
    \code{"countsbio"}, \code{"saturation"},
    \code{"rnacomp"}, \code{"readnoise"}, \code{"biodist"},
    \code{"gcbias"}, \code{"lengthbias"} or
    \code{"filtered"}.}

    \item{contrast.list}{a named structured list of contrasts
    as returned by \code{\link{make.contrast.list}} or just
    the vector of contrasts as defined in the main help page
    of \code{\link{metaseqr}}. This parameter is optional and
    required only when \code{diagplot.type} is any of
    \code{"deheatmap"}, \code{"volcano"} or
    \code{"biodist"}.}

    \item{p.list}{a list of p-values for each contrast as
    obtained from any of the \code{stat.*} methods of the
    metaseqr package. This parameter is optional and required
    only when \code{diagplot.type} is any of
    \code{"deheatmap"}, \code{"volcano"} or
    \code{"biodist"}.}

    \item{thresholds}{a list with the elements \code{"p"} and
    \code{"f"} which are the p-value and the fold change
    cutoff when \code{diagplot.type="volcano"}.}

    \item{diagplot.type}{one or more of the diagnostic plots
    supported in metaseqR package. Many of these plots
    require the presence of additional package, something
    that is checked while running the main metaseqr function.
    The supported plots are \code{"mds"},
    \code{"biodetection"}, \code{"countsbio"},
    \code{"saturation"}, \code{"rnacomp"}, \code{"boxplot"},
    \code{"gcbias"}, \code{"lengthbias"}, \code{"meandiff"},
    \code{"meanvar"}, \code{"deheatmap"}, \code{"volcano"},
    \code{"biodist"}, \code{"filtered"}, \code{"readnoise"},
    \code{"venn"}, \code{"correl"}, \code{"pairwise"}. For a
    brief description of these plots please see the main
    \code{\link{metaseqr}} help page.}

    \item{is.norm}{a logical indicating whether object
    contains raw or normalized data. It is not essential and
    it serves only plot annotation purposes.}

    \item{output}{one or more R plotting device to direct the
    plot result to. Supported mechanisms: \code{"png"},
    \code{"jpg"}, \code{"bmp"}, \code{"pdf"}, \code{"ps"} or
    \code{"json"}. The latter is currently available for the
    creation of interactive volcano plots only when reporting
    the output, through the highcharts javascript library.
    The default plotting (\code{"x11"}) is not supported due
    to instability in certain devices.}

    \item{path}{the path to create output files.}

    \item{...}{further arguments to be passed to plot
    devices, such as parameter from \code{\link{par}}.}
}
\value{
    A named list containing the file names of the produced
    plots. Each list member is names according to the
    selected plotting device and is also a named list, whose
    names are the plot types. The final contents are the file
    names in case the plots are written to a physical
    location (not meaningful for \code{"x11"}).
}
\description{
    This is the main function for producing sructured quality
    control and informative graphs base on the results of the
    various steps of the metaseqR package. The graphs
    produced span a variety of issues like good sample
    reproducibility (Multi-Dimensional Scaling plot, biotype
    detection, heatmaps. diagplot.metaseqr, apart from
    implementing certain package-specific plots, is a wrapper
    around several diagnostic plots present in other RNA-Seq
    analysis packages such as EDASeq and NOISeq.
}
\note{
    In order to make the best out of this function, you
    should generally provide the annotation argument as most
    and also the most informative plots depend on this. If
    you don't know what is inside your counts table or how
    many annotation elements you can provide by embedding it,
    it's always best to set the annotation parameter of the
    main metaseqr function to \code{"download"} to use
    predefined annotations that work better with the
    functions of the whole package.
}
\examples{
\donttest{
require(DESeq)
data.matrix <- counts(makeExampleCountDataSet())
sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
contrast <- "A_vs_B"
diagplot.metaseqr(data.matrix,sample.list,diagplot.type=c("mds","boxplot"))

norm.args <- get.defaults("normalization","deseq")
object <- normalize.deseq(data.matrix,sample.list,norm.args)
diagplot.metaseqr(object,sample.list,diagplot.type="boxplot")

p <- stat.deseq(object,sample.list)
diagplot.metaseqr(object,sample.list,contrast.list=contrast,p.list=p,
        diagplot.type="volcano")
}
}
\author{
    Panagiotis Moulos
}

