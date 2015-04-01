\name{make.sample.list}
\alias{make.sample.list}
\title{Creates sample list from file}
\usage{
    make.sample.list(input, type=c("simple","targets"))
}
\arguments{
    \item{input}{a tab-delimited file structured as follows:
    the first line of the external tab delimited file should
    contain column names (names are not important). The first
    column MUST contain UNIQUE sample names and the second
    column MUST contain the biological condition where each
    of the samples in the first column should belong to.}
    
    \item{type}{one of \code{"simple"} or \code{"targets"} to 
    indicate if the input is a simple two column text file or 
    the targets file used to launch the main analysis pipeline.}
}
\value{
    A named list whose names are the conditions of the
    experiments and its members are the samples belonging to
    each condition.
}
\description{
    Create the main sample list from an external file.
}
\examples{
\donttest{
targets <- data.frame(sample=c("C1","C2","T1","T2"),
    condition=c("Control","Control","Treatment","Treatment"))
write.table(targets,file="targets.txt",sep="\\t",row.names=FALSE,
    quote="")
sample.list <- make.sample.list("targets.txt")
}
}
\author{
    Panagiotis Moulos
}

