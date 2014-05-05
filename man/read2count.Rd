\name{read2count}
\alias{read2count}
\title{SAM/BAM/BED file reader helper for the metaseqr pipeline}
\usage{
    read2count(targets, annotation, file.type = targets$type,
        has.all.fields = FALSE, multic = FALSE)
}
\arguments{
    \item{targets}{a named list, the output of 
    \code{\link{read.targets}}.}

    \item{annotation}{see the \code{annotation} argument in
    the main \code{\link{metaseqr}} function. The
    \code{"annotation"} parameter here is the result of the
    same parameter in the main function. See also
    \code{\link{get.annotation}}.}

	\item{file.type}{the type of raw input files. It can be
    \code{"bed"} for BED files or \code{"sam"}, \code{"bam"}
    for SAM/BAM files. See the same argument in the main
    \code{\link{metaseqr}} function for the case of
    auto-guessing.}

    \item{has.all.fields}{a logical variable indicating if
    all annotation fields used by \code{metaseqr} are
    available (that is apart from the main chromosome, start,
    end, unique id and strand columns, if also present are
    the gene name and biotype columns). The default is
    \code{FALSE}.}

    \item{multic}{a logical value indicating the presence 
    of multiple cores. Defaults to \code{FALSE}. Do not 
    change it if you are not sure whether package parallel
    has been loaded or not.}
}
\value{
    A data frame with counts for each sample, ready to be
    passed to the main \code{\link{metaseqr}} pipeline.
}
\description{
    This function is a helper for the \code{metaseqr}
    pipeline, for reading SAM/BAM or BED files when a read
    counts file is not available.
}
\examples{
\dontrun{
my.targets <- read.targets("my_mm9_study_bam_files.txt")
sample.list <- my.targets$samples
file.list <- my.targets$files
gene.data <- get.annotation("mm9","gene")
r2c <- read2count(files.list=file.list,
    file.type=my.targets$type,annotation=gene.data)
gene.counts <- r2c$counts
libsize.list <- r2s$libsize
}
}
\author{
    Panagiotis Moulos
}

