\name{calc.otr}
\alias{calc.otr}
\title{Calculate the ratio TP/(FP+FN)}
\usage{
    calc.otr(truth, p, sig = 0.05)
}
\arguments{
    \item{truth}{the ground truth differential 
    expression vector. It should contain only 
    zero and non-zero elements, with zero denoting
    non-differentially expressed genes and non-zero, 
    differentially expressed genes. Such a vector 
    can be obtained for example by using the 
    \code{\link{make.sim.data.sd}} function, which 
    creates simulated RNA-Seq read counts based on 
    real data. It MUST be named with gene names, 
    the same as in \code{p}.}

    \item{p}{a p-value matrix whose rows correspond 
    to each element in the \code{truth} vector. If 
    the matrix has a \code{colnames} attribute, a 
    legend will be added to the plot using these 
    names, else a set of column names will be 
    auto-generated. \code{p} can also be a list or 
    a data frame. In any case, each row (or element) 
    MUST be named with gene names (the same as in 
    \code{truth}).}

    \item{sig}{a significance level (0 < \code{sig} 
    <=1).}
}
\value{
    A named list with two members. The first member 
    is a data frame with the numbers used to 
    calculate the TP/(FP+FN) ratio and the second 
    member is the ratio TP/(FP+FN) for each 
    statistical test.
}
\description{
    This function calculates the ratio of True 
    Positives to the sum of False Positives and 
    False Negatives given a matrix of p-values 
    (one for each statistical test used) and a 
    vector of ground truth (DE or non-DE). This
    function serves as a method evaluation helper.
}
\examples{
p1 <- 0.001*matrix(runif(300),100,3)
p2 <- matrix(runif(300),100,3)
p <- rbind(p1,p2)
rownames(p) <- paste("gene",1:200,sep="_")
colnames(p) <- paste("method",1:3,sep="_")
truth <- c(rep(1,40),rep(-1,40),rep(0,20),rep(1,10),
    rep(2,10),rep(0,80))
names(truth) <- rownames(p)
otr <- calc.otr(truth,p)
}
\author{
    Panagiotis Moulos
}

