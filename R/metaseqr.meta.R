#' Meta-analysis using several RNA-Seq statistics
#'
#' This function calculates the combined p-values when multiple statistical algorithms
#' are applied to the input dataset. It is a helper and it requires very specific 
#' arguments so it should not be used individually.
#'
#' @param cp.list a named list whose names are the contrasts requested from metaseqR.
#' Each member is a p-value matrix whose colnames are the names of the statistical
#' tests applied to the data. See the main \code{\link{metaseqr}} help page.
#' @param meta.p the p-value combination method to use. See the main 
#' \code{\link{metaseqr}} help page.
#' @param counts the normalized and possibly filtered read counts matrix. See the
#' main \code{\link{metaseqr}} help page.
#' @param sample.list the list containing condition names and the samples under
#' each condition. See the main \code{\link{metaseqr}} help page.
#' @param statistics the statistical algorithms used in metaseqr. See the main 
#' \code{\link{metaseqr}} help page.
#' @param stat.args the parameters for each statistical argument. See the main
#' \code{\link{metaseqr}} help page.
#' @param libsize.list a list with library sizes. See the main \code{\link{metaseqr}}
#' and the \code{stat.*} help pages.
#' @param nperm the number of permutations (Monte Carlo simulations) to perform.
#' @param weight a numeric vector of weights for each statistical algorithm.
#' @param reprod create reproducible permutations. Ideally one would want to create
#' the same set of indices for a given dataset so as to create reproducible p-values.
#' If \code{reprod=TRUE}, a fixed seed is used by \code{meta.perm} for all the
#' datasets analyzed with \code{metaseqr}. If \code{reprod=FALSE}, then the
#' p-values will not be reproducible, although statistical significance is not
#' expected to change for a large number of resambling. Finally, \code{reprod}
#' can be a numeric vector of seeds with the same length as \code{nperm} so that
#' the user can supply his/her own seeds.
#' @param multic use multiple cores to execute the premutations. This is an
#' external parameter and implies the existence of parallel package in the execution
#' environment. See the main \code{\link{metaseqr}} help page.
#' @return A named list with combined p-values. The names are the contrasts and
#' the list members are combined p-value vectors, one for each contrast.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # This function is not exported
#'}
meta.test <- function(cp.list,meta.p=c("simes","bonferroni","fisher",
    "dperm.min","dperm.max","dperm.weight","fperm","whitlock","minp","maxp",
    "weight","pandora","none"),counts,sample.list,statistics,stat.args,
    libsize.list,nperm=10000,weight=rep(1/length(statistics),
    length(statistics)),reprod=TRUE,multic=FALSE) {
    check.text.args("meta.p",meta.p,c("simes","bonferroni","fisher","dperm.min",
        "dperm.max","dperm.weight","fperm","whitlock","minp","maxp","weight",
        "pandora","none"))
    contrast <- names(cp.list)
    disp("Performing meta-analysis with ",meta.p)
    if (meta.p=="pandora")
        meta.p <- "weight"
    switch(meta.p,
        fisher = {
            sum.p.list <- wapply(multic,cp.list,function(x) {
                tmp <- fisher.method(x,p.corr="none",
                    zero.sub=.Machine$double.xmin)
                rp <- tmp$p.value
                names(rp) <- rownames(x)
                return(rp)
            })
        },
        fperm = {
            sum.p.list <- wapply(multic,cp.list,function(x) {
                if (multic)
                    tmp <- fisher.method.perm(x,p.corr="none",B=nperm,
                        mc.cores=getOption("cores"),zero.sub=1e-32)
                else
                    tmp <- fisher.method.perm(x,p.corr="none",B=nperm,
                        zero.sub=.Machine$double.xmin)
                return(tmp$p.value)
            })
        },
        whitlock = {
            sum.p.list <- wapply(multic,cp.list,function(x) 
                return(apply(x,1,combine.test,method="z.transform")))
        },
        simes = {
            sum.p.list <- wapply(multic,cp.list,function(x) {
                return(apply(x,1,combine.simes))
            })
        },
        bonferroni = {
            sum.p.list <- wapply(multic,cp.list,function(x) {
                return(apply(x,1,combine.bonferroni))
            })
        },
        minp = {
            sum.p.list <- wapply(multic,cp.list,function(x) {
                return(apply(x,1,combine.minp))
            })
        },
        maxp = {
            sum.p.list <- wapply(multic,cp.list,function(x) {
                return(apply(x,1,combine.maxp))
            })
        },
        weight = {
            sum.p.list <- wapply(multic,cp.list,function(x) {
                return(apply(x,1,combine.weight,weight))
            })
        },
        dperm.min = {
            sum.p.list <- vector("list",length(cp.list))
            names(sum.p.list) <- names(cp.list)
            conl <- as.list(contrast)
            names(conl) <- contrast
            temp.p.list <- wapply(multic,conl,meta.perm,
                counts=counts,sample.list=sample.list,
                statistics=statistics,stat.args=stat.args,
                libsize.list=libsize.list,
                nperm=nperm,weight=weight,
                select="min",reprod=reprod,
                multic=multic)
            original.p.list <- wapply(multic,cp.list,function(x,m,w=NULL) {
                x[which(is.na(x))] <- 1
                switch(m,
                    min = {
                        return(apply(x,1,min))
                    },
                    max = {
                        return(apply(x,1,max))
                    },
                    weight = {
                        return(apply(x,1,function(p,w) return(prod(p^w)),
                            w))
                    }
                )
            },"min")
            for (cc in names(original.p.list))
            {
                pc <- cbind(temp.p.list[[cc]],original.p.list[[cc]])
                ly <- ncol(pc)
                sum.p.list[[cc]] <- apply(pc,1,function(y,m) 
                    return(length(which(y[1:(m-1)]<y[m]))/(m-1)),ly)
            }
            #assign("perm.list",temp.p.list,envir=.GlobalEnv)
            #assign("o.list",original.p.list,envir=.GlobalEnv)
        },
        dperm.max = {
            sum.p.list <- vector("list",length(cp.list))
            names(sum.p.list) <- names(cp.list)
            conl <- as.list(contrast)
            names(conl) <- contrast
            temp.p.list <- wapply(multic,conl,meta.perm,
                counts=counts,sample.list=sample.list,
                statistics=statistics,stat.args=stat.args,
                libsize.list=libsize.list,
                nperm=nperm,weight=weight,
                select="max",reprod=reprod,
                multic=multic)
            original.p.list <- wapply(multic,cp.list,function(x,m,w=NULL) {
                switch(m,
                    min = {
                        return(apply(x,1,min))
                    },
                    max = {
                        return(apply(x,1,max))
                    },
                    weight = {
                        return(apply(x,1,function(p,w) return(prod(p^w)),
                            w))
                    }
                )
            },"max")
            for (cc in names(original.p.list))
            {
                pc <- cbind(temp.p.list[[cc]],original.p.list[[cc]])
                ly <- ncol(pc)
                sum.p.list[[cc]] <- apply(pc,1,function(y,m) 
                    return(length(which(y[1:(m-1)]<y[m]))/(m-1)),ly)
            }
            #assign("perm.list",temp.p.list,envir=.GlobalEnv)
            #assign("o.list",original.p.list,envir=.GlobalEnv)
        },
        dperm.weight = {
            sum.p.list <- vector("list",length(cp.list))
            names(sum.p.list) <- names(cp.list)
            conl <- as.list(contrast)
            names(conl) <- contrast
            temp.p.list <- wapply(multic,conl,meta.perm,
                counts=counts,sample.list=sample.list,
                statistics=statistics,stat.args=stat.args,
                libsize.list=libsize.list,
                nperm=nperm,weight=weight,
                select="weight",reprod=reprod,
                multic=multic)
            original.p.list <- wapply(multic,cp.list,function(x,m,w=NULL) {
                switch(m,
                    min = {
                        return(apply(x,1,min))
                    },
                    max = {
                        return(apply(x,1,max))
                    },
                    weight = {
                        return(apply(x,1,function(p,w) {return(prod(p^w))},
                            w))
                    }
                )
            },"weight",weight)
            for (cc in names(original.p.list))
            {
                pc <- cbind(temp.p.list[[cc]],original.p.list[[cc]])
                ly <- ncol(pc)
                sum.p.list[[cc]] <- apply(pc,1,function(y,m) 
                    return(length(which(y[1:(m-1)]<y[m]))/(m-1)),ly)
            }
            #assign("perm.list",temp.p.list,envir=.GlobalEnv)
            #assign("o.list",original.p.list,envir=.GlobalEnv)
        },
        none = {
            # A default value must be there to use with volcanos, we say the one
            # of the first statistic in order of input
            sum.p.list <- wapply(multic,cp.list,function(x) return(x[,1]))
        }
    )
    return(sum.p.list)
}

#' Permutation tests for meta-analysis
#'
#' This function performs permutation tests in order to derive a meta p-value by
#' combining several of the statistical algorithms of metaseqr. This is probably
#' the most accurate way of combining multiple statistical algorithms for RNA-Seq
#' data, as this issue is different from the classic interpretation of the term
#' "meta-analysis" which implies the application of the same statistical test on
#' different datasets treating the same subject/experiment. For other methods, see
#' also the main \code{\link{metaseqr}} help page. You should keep in mind that
#' the permutation procedure can take a long time, even when executed in parallel.
#'
#' @param counts a normalized read counts table, one row for each gene, one column
#' for each sample.
#' @param sample.list the list containing condition names and the samples under
#' each condition. See the main \code{\link{metaseqr}} help page.
#' @param contrast the contrasts to be tested by each statistical algorithm. See
#' the main \code{\link{metaseqr}} help page.
#' @param statistics the statistical algorithms used in metaseqr. See the main
#' \code{\link{metaseqr}} help page.
#' @param stat.args the parameters for each statistical algorithm. See the main
#' \code{\link{metaseqr}} help page.
#' @param libsize.list a list with library sizes. See the main \code{\link{metaseqr}}
#' and the \code{stat.*} help pages.
#' @param nperm the number of permutations (Monte Carlo simulations) to perform.
#' @param weight a numeric vector of weights for each statistical algorithm.
#' @param select how to select the initial vector of p-values. It can be \code{"min"}
#' to select the minimum p-value for each gene (more conservative), \code{"max"}
#' to select the maximum p-value for each gene (less conservative), \code{"weight"}
#' to apply the weights to the p-value vector for each gene and derive a weighted
#' p-value.
#' @param replace same as the \code{replace} argument in the \code{\link{sample}}
#' function. Implies bootstraping or simple resampling without replacement. It can
#' also be \code{"auto"}, to determine bootstraping or not with the following rule:
#' if \code{ncol(counts)<=6} \code{replace=FALSE else} \code{replace=TRUE}. This
#' protects from the case of having zero variability across resampled conditions.
#' In such cases, most statistical tests would crash.
#' @param multic use multiple cores to execute the premutations. This is an 
#' external parameter and implies the existence of parallel package in the
#' execution environment. See the main \code{\link{metaseqr}} help page.
#' @param reprod create reproducible permutations. Ideally one would want to
#' create the same set of indices for a given dataset so as to create reproducible
#' p-values. If \code{reprod=TRUE}, a fixed seed is used by \code{meta.perm} for
#' all the datasets analyzed with \code{metaseqr}. If \code{reprod=FALSE}, then
#' the p-values will not be reproducible, although statistical significance is not
#' expected to change for a large number of resambling. Finally, \code{reprod} can
#' be a numeric vector of seeds with the same length as \code{nperm} so that the
#' user can supply his/her own seeds.
#' @return A vector of meta p-values
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # This function is not exported
#'}
meta.perm <- function(contrast,counts,sample.list,statistics,stat.args,
    libsize.list,nperm=10000,weight=rep(1/ncol(counts),ncol(counts)),
    select=c("min","max","weight"),replace="auto",reprod=TRUE,multic=FALSE) {
    check.text.args("select",select,c("min","max","weight"))
    if (replace=="auto") {
        if (ncol(counts)<=6)
            replace=FALSE
        else
            replace=TRUE
    }
    # We will construct relist in a way so that we can assign seeds for random
    # number generation and track progress at the same time
    if (is.logical(reprod)) {
        relist <- vector("list",nperm)
        if (reprod) {
            relist <- wapply(multic,seq_along(relist),function(i) 
                {return(list(seed=i,prog=i))})
        }
        else
            relist <- wapply(multic,seq_along(relist),function(i) 
                {return(list(seed=round(1e+6*runif(1)),prog=i))})
    }
    else if (is.numeric(reprod)) {
        if (length(reprod) != nperm)
            stopwrap("When reprod is numeric, it must have the same length as ",
                "nperm!")
        relist <- wapply(multic,seq_along(reprod),function(i) 
            {return(list(seed=reprod[i],prog=i))})
    }
    else
        stopwrap("reprod must be either a logical or a numeric vector!")
    disp("  Resampling procedure started...")
    # In this case, we must not use wapply as we want to be able to track progress
    # through mc.preschedule...
    if (multic)
        pp <- mclapply(relist,meta.worker,counts,sample.list,contrast,
            statistics,replace,stat.args,libsize.list,select,weight,
            mc.preschedule=FALSE,mc.cores=getOption("cores"))
    else
        pp <- lapply(relist,meta.worker,counts,sample.list,contrast,statistics,
            replace,stat.args,libsize.list,select,weight)
    disp("  Resampling procedure ended...")
    return(do.call("cbind",pp))
}

#' Permutation tests helper
#'
#' This function performs the statistical test for each permutation. Internal use
#' only.
#'
#' @param x a virtual list with the random seed and the permutation index.
#' @param co the counts matrix.
#' @param sl the sample list.
#' @param cnt the contrast name.
#' @param s the statistical algorithms.
#' @param sa the parameters for each statistical algorithm.
#' @param ll a list with library sizes.
#' @param r same as the \code{replace} argument in the \code{\link{sample}} function.
#' @param el min, max or weight.
#' @param w the weights when \code{el="weight"}.
#' @return A matrix of p-values.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # This function is not exported
#'}
meta.worker <- function(x,co,sl,cnt,s,r,sa,ll,el,w) {
    set.seed(x$seed)
    disp("    running permutation #",x$prog)
    pl <- make.permutation(co,sl,cnt,r)
    ppmat <- matrix(NA,nrow(co),length(s))
    colnames(ppmat) <- s
    for (alg in s) {
        #disp("      running permutation tests with: ",alg)
        tcl <- make.contrast.list(pl$contrast,pl$sample.list)
        switch(alg,
            deseq = {
                p.list <- suppressMessages(stat.deseq(pl$counts,pl$sample.list,
                    tcl,sa[[alg]]))
            },
            edger = {
                p.list <- suppressMessages(stat.edger(pl$counts,pl$sample.list,
                    tcl,sa[[alg]]))
            },
            noiseq = {
                p.list <- suppressMessages(stat.noiseq(pl$counts,pl$sample.list,
                    tcl,sa[[alg]]))
            },
            bayseq = {
                p.list <- suppressMessages(stat.bayseq(pl$counts,pl$sample.list,
                    tcl,sa[[alg]],ll))
            },
            limma = {
                p.list <- suppressMessages(stat.limma(pl$counts,pl$sample.list,
                    tcl,sa[[alg]]))
            },
            nbpseq = {
                p.list <- suppressMessages(stat.nbpseq(pl$counts,pl$sample.list,
                    tcl,sa[[alg]],ll))
            }
        )
        ppmat[,alg] <- as.numeric(p.list[[1]])
    }
    ppmat[which(is.na(ppmat))] <- 1
    switch(el,
        min = {
            p.iter <- apply(ppmat,1,min)
        },
        max = {
            p.iter <- apply(ppmat,1,max)
        },
        weight = {
            p.iter <- apply(ppmat,1,function(p,w) return(prod(p^w)),w)
        }
    )
    return(p.iter)
}

#' Combine p-values with Simes' method
#'
#' This function combines p-values from the various statistical tests supported by
#' metaseqR using the Simes' method (see reference in the main \code{\link{metasqr}}
#' help page or in the vignette).
#'
#' @param p a p-value matrix (rows are genes, columns are statistical tests).
#' @return A vector of combined p-values.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p <- matrix(runif(300),100,3)
#' pc <- combine.simes(p)
#'}
combine.simes <- function(p) {
    ze <- which(p==0)
    if (length(ze)>0)
        p[ze] <- 0.1*min(p[-ze])
    m <- length(p)
    y <- sort(p)
    s <- min(m*(y/(1:m)))
    return(min(c(s,1)))
}

#' Combine p-values with Bonferroni's method
#'
#' This function combines p-values from the various statistical tests supported by
#' metaseqR using the Bonferroni's method (see reference in the main
#' \code{\link{metasqr}} help page or in the vignette).
#'
#' @param p a p-value matrix (rows are genes, columns are statistical tests).
#' @return A vector of combined p-values.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p <- matrix(runif(300),100,3)
#' pc <- combine.bonferroni(p)
#'}
combine.bonferroni <- function(p) {
    ze <- which(p==0)
    if (length(ze)>0)
        p[ze] <- 0.1*min(p[-ze])
    b <- length(p)*min(p)
    return(min(c(1,b)))
}

#' Combine p-values using weights
#'
#' This function combines p-values from the various statistical tests supported by
#' metaseqR using p-value weights.
#'
#' @param p a p-value matrix (rows are genes, columns are statistical tests).
#' @param w a weights vector, must sum to 1.
#' @return A vector of combined p-values.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p <- matrix(runif(300),100,3)
#' pc <- combine.weight(p,w=c(0.2,0.5,0.3))
#'}
combine.weight <- function(p,w) {
    ze <- which(p==0)
    if (length(ze)>0)
        p[ze] <- 0.1*min(p[-ze])
    return(prod(p^w))
}

#' Combine p-values using the minimum p-value
#'
#' This function combines p-values from the various statistical tests supported by
#' metaseqR by taking the minimum p-value.
#'
#' @param p a p-value matrix (rows are genes, columns are statistical tests).
#' @return A vector of combined p-values.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p <- matrix(runif(300),100,3)
#' pc <- combine.min(p)
#'}
combine.minp <- function(p) { return(min(p)) }

#' Combine p-values using the maximum p-value
#'
#' This function combines p-values from the various statistical tests supported by
#' metaseqR by taking the maximum p-value.
#'
#' @param p a p-value matrix (rows are genes, columns are statistical tests).
#' @return A vector of combined p-values.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p <- matrix(runif(300),100,3)
#' pc <- combine.max(p)
#'}
combine.maxp <- function(p) { return(max(p)) }

# Copied from ex-CRAN package MADAM and exported. The man pages are copied from
# the original package.
fisher.method <- function(pvals,method=c("fisher"),p.corr=c("bonferroni","BH",
    "none"),zero.sub=0.00001,na.rm=FALSE,mc.cores=NULL) {
    stopifnot(method %in% c("fisher"))
    stopifnot(p.corr %in% c("none","bonferroni","BH"))
    stopifnot(all(pvals>=0, na.rm=TRUE) & all(pvals<=1, na.rm=TRUE))
    stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
    if(is.null(dim(pvals)))
        stop("pvals must have a dim attribute")
    p.corr <- ifelse(length(p.corr)!=1, "BH", p.corr)
    ##substitute p-values of 0
    pvals[pvals == 0] <- zero.sub
    if(is.null(mc.cores)) {
        fisher.sums <- data.frame(do.call(rbind,apply(pvals,1,fisher.sum,
            zero.sub=zero.sub,na.rm=na.rm)))
    } 
    else {
        fisher.sums <- parallel::mclapply(1:nrow(pvals), function(i) {
            fisher.sum(pvals[i,],zero.sub=zero.sub,na.rm=na.rm)
        }, mc.cores=mc.cores)
        fisher.sums <- data.frame(do.call(rbind,fisher.sums))
    }
    
  rownames(fisher.sums) <- rownames(pvals)
  fisher.sums$p.value <- 1-pchisq(fisher.sums$S,df=2*fisher.sums$num.p)
  fisher.sums$p.adj <- switch(p.corr,
      bonferroni = p.adjust(fisher.sums$p.value,"bonferroni"),
      BH = p.adjust(fisher.sums$p.value,"BH"),
      none = fisher.sums$p.value
    )
    return(fisher.sums)
}

# Copied from ex-CRAN package MADAM and exported. The man pages are copied from
# the original package.
fisher.method.perm <- function(pvals,p.corr=c("bonferroni","BH","none"),
    zero.sub=0.00001,B=10000,mc.cores=NULL,blinker=1000) {
    stopifnot(is.na(blinker) || blinker>0)
    stopifnot(p.corr %in% c("none","bonferroni","BH"))
    stopifnot(all(pvals>=0,na.rm=TRUE) & all(pvals<=1,na.rm=TRUE))
    stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
    if(is.null(dim(pvals)))
        stop("pvals must have a dim attribute")
    p.corr <- ifelse(length(p.corr)!=1,"BH",p.corr)
    pvals[pvals==0] <- zero.sub
  
    res.perm <- lapply(1:nrow(pvals),function(i) {
        if(!is.na(blinker) & i%%blinker==0)
        message("=", appendLF=FALSE)
        ##which studies contribute to S (don't have a NA in row i)
        good.p <- which(!is.na(pvals[i,]))
        S.obs= fisher.sum(pvals[i,good.p], na.rm=FALSE)
        if(is.null(mc.cores)) {
            S.rand <- unlist(lapply(1:B, function(b) {
            ##get non NA p-values from studies contributing to S
            myp <- sapply(good.p, function(pc){
                sample(na.exclude(pvals[,pc]),1)
            })
            fisher.sum(myp)$S
        }))
        } else {
        S.rand <- unlist(parallel::mclapply(1:B, function(b) {
            ##get non NA p-values from studies contributing to S
            myp <- sapply(good.p, function(pc) {
                sample(na.exclude(pvals[,pc]),1)
            })
            fisher.sum(myp)$S
            }, mc.cores=mc.cores))
        }
        p.value <- sum(S.rand>=S.obs$S)/B
        data.frame(S=S.obs$S, num.p=S.obs$num.p, p.value=p.value)
    })
    res.perm <- data.frame(do.call(rbind, res.perm))
  
    if(!is.na(blinker) && blinker>0)
        message()
    ## rownames(res.perm) <- rownames(pvals)
    res.perm$p.adj <- switch(p.corr,
      bonferroni = p.adjust(res.perm$p.value,"bonferroni"),
      BH = p.adjust(res.perm$p.value,"BH"),
      none = res.perm$p.value)
    return(res.perm)
}

# Copied from ex-CRAN package MADAM and exported. The man pages are copied from
# the original package.
fisher.sum <- function(p,zero.sub=0.00001,na.rm=FALSE) {
    if(any(p>1, na.rm=TRUE)||any(p<0, na.rm=TRUE))
        stop("You provided bad p-values")
    stopifnot(zero.sub>=0 & zero.sub<=1 || length(zero.sub)!=1)
    p[p==0] <- zero.sub
    if (na.rm)
        p <- p[!is.na(p)]
    S = -2*sum(log(p))
    res <- data.frame(S=S,num.p=length(p))
    return(res)
}
