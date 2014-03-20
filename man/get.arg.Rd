\name{get.arg}
\alias{get.arg}
\title{Argument getter}
\usage{
    get.arg(arg.list, arg.name)
}
\arguments{
    \item{arg.list}{the initial list of a method's (e.g.
    normalization) arguments. Can be created with the
    \code{\link{get.defaults}} function.}

    \item{arg.name}{the argument name inside the argument
    list to fetch its value.}
}
\value{
    The argument sub-list.
}
\description{
    Get argument(s) from a list of arguments, e.g.
    normalization arguments.
}
\examples{
\donttest{
norm.list <- get.defaults("normalization","egder")
a <- get.arg(norm.list,c("main.method","logratioTrim"))
}
}
\author{
    Panagiotis Moulos
}

