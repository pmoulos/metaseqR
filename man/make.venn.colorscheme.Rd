\name{make.venn.colorscheme}
\alias{make.venn.colorscheme}
\title{Helper for Venn diagrams}
\usage{
    make.venn.colorscheme(n)
}
\arguments{
    \item{n}{the number of the sets used for the Venn
    diagram.}
}
\value{
    A list with colors for fill and font.
}
\description{
    This function returns a list of colorschemes accroding to
    the number of sets. Internal use.
}
\examples{
\donttest{
sets <- c("apple","pear","banana")
cs <- make.venn.colorscheme(length(sets))
}
}
\author{
    Panagiotis Moulos
}

