\name{set.arg}
\alias{set.arg}
\title{Argument setter}
\usage{
    set.arg(arg.list, arg.name, arg.value = NULL)
}
\arguments{
    \item{arg.list}{the initial list of a method's (e.g.
    normalization) arguments. Can be created with the
    \code{\link{get.defaults}} function.}

    \item{arg.name}{a named list with names the new arguments
    to be set, and mebers the values to be set or a vector of
    argument names. In this case, \code{arg.value} must be
    supplied.}

    \item{arg.value}{when \code{arg.name} is a vector of
    argument names, the values corresponding to these
    arguments.}
}
\value{
    the \code{arg.list} with the changed \code{arg.value} for
    \code{arg.name}.
}
\description{
    Set argument(s) to a list of arguments, e.g.
    normalization arguments.
}
\examples{
\donttest{
norm.list <- get.defaults("normalization","egder")
set.arg(norm.list,list(main.method="glm",logratioTrim=0.4))
}
}
\author{
    Panagiotis Moulos
}
