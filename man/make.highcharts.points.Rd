\name{make.highcharts.points}
\alias{make.highcharts.points}
\title{Interactive volcano plot helper}
\usage{
    make.highcharts.points(x, y, a)
}
\arguments{
    \item{x}{The x coordinates (should be a named vector!).}

    \item{y}{The y coordinates.}

    \item{a}{Alternative names for each point.}
}
\value{
    A list that is later serialized to JSON.
}
\description{
    Creates a list which contains the data series of a
    scatterplot, to be used for serialization with highcharts
    JavaScript plotting. framework. Internal use only.
}
\author{
    Panagiotis Moulos
}

