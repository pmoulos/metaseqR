\name{check.num.args}
\alias{check.num.args}
\title{Numeric argument validator}
\usage{
    check.num.args(arg.name, arg.value, arg.type, arg.bounds,
        direction)
}
\arguments{
    \item{arg.name}{the name of the argument that is checked
    (for display purposes).}

    \item{arg.value}{the value(s) of the argument to be
    checked.}

    \item{arg.type}{either the string \code{"numeric"} to
    denote generic double-like R numerics or \code{"integer"}
    for integer values.}

    \item{arg.bounds}{a numeric or a vector with 2 elements,
    restraining \code{arg.value} to be within the bounds
    defined by the input vector or e.g. larger (smaller) than
    the numeric value. See examples.}

    \item{direction}{a string denoting to which direction the
    \code{arg.value} should be compared with
    \code{arg.bounds}. For example, \code{"both"} should be
    given with a two element vector against which,
    \code{arg.value} will be checked to see whether it is
    smaller than the low boundary or larger than the higher
    boundary. In that case, the function will throw an error.
    The direction parameter can be one of: \code{"both"}
    (described above), \code{"botheq"} (as above, but the
    \code{arg.val} is also checked for equality -closed
    intervals), \code{"gt"} or \code{"gte"} (check whether
    \code{arg.val} is smaller or smaller than or equal to the
    first value of \code{arg.bounds}), \code{"lt"} or
    \code{"lte"} (check whether \code{arg.val} is larger or
    larger than or equal to the first value of
    \code{arg.bounds}).}
}
\description{
    Checks if one or more given numeric argument(s) satisfy
    several rules concerning numeric arguments, e.g. proper
    bounds or proper format (e.g. it must be a number and not
    a character). Mostly for internal use.
}
\examples{
\donttest{
pcut <- 1.2 # A probability cannot be larger than 1! It will throw an error!
#check.num.args("pcut",pcut,"numeric",c(0,1),"botheq")
pcut <- 0.05 # Pass
check.num.args("pcut",pcut,"numeric",c(0,1),"botheq")
gc.col <- 3.4 # A column in a file cannot be real! It will throw an error!
#check.num.args("gc.col",gc.col,"integer",0,"gt")
gc.col <- 5L # Pass
check.num.args("gc.col",gc.col,"integer",0,"gt")
}
}
\author{
    Panagiotis Moulos
}
