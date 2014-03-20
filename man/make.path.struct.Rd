\name{make.path.struct}
\alias{make.path.struct}
\title{Project path constructor helper}
\usage{
    make.path.struct(main.path)
}
\arguments{
    \item{main.path}{The desired project path.}
}
\value{
    A named list whose names are the conditions of the
    experiments and its members are the samples belonging to
    each condition.
}
\description{
    Helper for \code{make.project.path}. Internal use only.
}
\author{
    Panagiotis Moulos
}

