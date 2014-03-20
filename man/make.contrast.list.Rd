\name{make.contrast.list}
\alias{make.contrast.list}
\title{Create contrast lists from contrast vectors}
\usage{
    make.contrast.list(contrast, sample.list)
}
\arguments{
    \item{contrast}{a vector of contrasts in the form
    "ConditionA_vs_ConditionB" or "ConditionA_
    vs_ConditionB_vs_ConditionC_vs_...". In case of Control
    vs Treatment designs, the Control condition should ALWAYS
    be the first.}

    \item{sample.list}{the list of samples in the experiment.
    See also the main help page of \code{\link{metaseqr}}.}
}
\value{
    A named list whose names are the contrasts and its
    members are named vectors, where the names are the sample
    names and the actual vector members are the condition
    names. See the example.
}
\description{
    Returns a list, properly structured to be used within the
    \code{stat.*} functions of the metaseqr package. See the
    main documentation for the structure of this list and the
    example below. This function is mostly for internal use,
    as the \code{stat.*} functions can be supplied directly
    with the contrasts vector which is one of the main
    \code{\link{metaseqr}} arguments.
}
\examples{
sample.list <- list(Control=c("C1","C2"),TreatmentA=c("TA1","TA2"),TreatmentB=c("TB1","TB2"))
contrast <- c("Control_vs_TreatmentA","Control_vs_TreatmentA_vs_TreatmentB")
cl <- make.contrast.list(contrast,sample.list)
}
\author{
    Panagiotis Moulos
}

