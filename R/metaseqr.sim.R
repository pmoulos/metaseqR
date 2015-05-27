#' Estimate AUFC weights
#'
#' This function automatically estimates weights for the \code{"weight"} and
#' \code{"dperm.weight"} options of metaseqR for combining p-values from multiple
#' statistical tests. It creates simulated dataset based on real data and then
#' performs statistical analysis with metaseqR several times in order to derive
#' False Discovery Curves. Then, the average areas under the false discovery curves
#' are used to construct weights for each algorithm, according to its performance
#' when using simulated data.
#'
#' @param counts the real raw counts table from which the simulation parameters
#' will be estimated. It must not be normalized and must contain only integer
#' counts, without any other annotation elements and unique gene identifiers as
#' the rownames attribute.
#' @param normalization same as \code{normalization} in \code{link{metaseqr}}.
#' @param statistics same as \code{statistics} in \code{link{metaseqr}}.
#' @param nsim the number of simulations to perform to estimate the weights. It
#' default to 10.
#' @param N the number of genes to produce. See \code{link{make.sim.data.sd}}.
#' @param samples a vector with 2 integers, which are the number of samples for
#' each condition (two conditions currently supported).
#' @param ndeg a vector with 2 integers, which are the number of differentially
#' expressed genes to be produced. The first element is the number of up-regulated
#' genes while the second is the number of down-regulated genes.
#' @param fc.basis the minimum fold-change for deregulation.
#' @param top the top \code{top} best ranked (according to p-value) to use, to
#' calculate area under the false discovery curve.
#' @param model.org the organism from which the data are derived. It must be one
#' of \code{\link{metaseqr}} supported organisms.
#' @param seed a number to be used as seed for reproducible simulations. Defaults
#' to \code{NULL} (NOT reproducible results!).
#' @param draw.fpc draw the averaged false discovery curves? Default to \code{FALSE}.
#' @param multic whether to run in parallel (if package \code{parallel} is present
#' or not.
#' @param ... Further arguments to be passed to \code{\link{estimate.sim.params}}.
#' @value A vector of weights to be used in \link{\code{metaseqr}} with the
#' \code{weights} option.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data("mm9.gene.data",package="metaseqR")
#' multic <- check.parallel(0.8)
#' weights <- estimate.aufc.weights(
#'   counts=as.matrix(mm9.gene.counts[,9:12]),
#'   normalization="edaseq",
#'   statistics=c("deseq","edger"),
#'   nsim=3,N=100,ndeg=c(10,10),top=10,model.org="mm9",
#'   seed=100,multic=multic,libsize.gt=1e+5
#' )
#'}
estimate.aufc.weights <- function(counts,normalization,statistics,nsim=10,
    N=10000,samples=c(3,3),ndeg=c(500,500),top=500,model.org="mm9",fc.basis=1.5,
    seed=NULL,draw.fpc=FALSE,multic=FALSE,...) {
    if (!require(zoo))
        stopwrap("R pacakage zoo is required in order to estimate AUFC ",
            "weights!")

    if (is.null(seed)) {
        seed.start <- round(100*runif(1))
        seed.end <- seed.start + nsim - 1
        seed <- as.list(seed.start:seed.end)
    }
    else {
        set.seed(seed)
        seed.start <- round(100*runif(1))
        seed.end <- seed.start + nsim - 1
        seed <- as.list(seed.start:seed.end)
    }
    
    if (ncol(counts)<4)
        stopwrap("Cannot estimate AUFC weights with an initial dataset with ",
            "less than 4 samples!")
    else if (ncol(counts)>=4 && ncol(counts)<10) {
        set.seed(seed.start)
        reind <- sample(1:ncol(counts),20,replace=TRUE)
        counts <- counts[,reind]
    }
    par.list <- estimate.sim.params(counts,...)

    disp("Running simulations... This procedure requires time... Please ",
        "wait...")
    sim.results <- wapply(multic,seed,function(x,normalization,statistics,N,
        par.list,samples,ndeg,fc.basis,model.org) {
        D <- make.sim.data.sd(N=N,param=par.list,samples=samples,ndeg=ndeg,
            fc.basis=fc.basis,model.org=model.org,seed=x)
        dd <- D$simdata
        
        if (!is.null(model.org)) {
            tmp <- metaseqr(
                counts=dd,
                sample.list=list(G1=paste("G1_rep",1:samples[1],sep=""),
                    G2=paste("G2_rep",1:samples[2],sep="")),
                contrast=c("G1_vs_G2"),
                annotation="embedded",
                id.col=4,
                gc.col=5,
                name.col=7,
                bt.col=8,
                org=model.org,
                count.type="gene",
                normalization=normalization,
                statistics=statistics,
                meta.p="simes",
                fig.format="png",
                preset="all.basic",
                export.where=tempdir(),
                qc.plots=NULL,
                report=FALSE,
                run.log=FALSE,
                out.list=TRUE
            )
        }
        else {
            tmp <- metaseqr(
                counts=dd,
                sample.list=list(G1=paste("G1_rep",1:samples[1],sep=""),
                    G2=paste("G2_rep",1:samples[2],sep="")),
                contrast=c("G1_vs_G2"),
                annotation="embedded",
                id.col=4,
                gc.col=5,
                name.col=7,
                bt.col=8,
                count.type="gene",
                normalization=normalization,
                statistics=statistics,
                meta.p="simes",
                fig.format="png",
                preset="all.basic",
                export.where=tempdir(),
                qc.plots=NULL,
                report=FALSE,
                run.log=FALSE,
                out.list=TRUE
            )
        }

        # Retrieve several p-values
        p.list <- vector("list",length(statistics))
        for (s in statistics) {
            field <- paste("p-value",s,sep="_")
            p.list[[s]] <- tmp$data[[1]][,field]
            names(p.list[[s]]) <- rownames(tmp$data[[1]])
        }
        p.matrix <- do.call("cbind",p.list)
        return(list(simdata=D,pvalues=p.matrix))
    },normalization,statistics,N,par.list,samples,ndeg,fc.basis,model.org)

    disp("Estimating AUFC weights... Please wait...")
    fpc.obj <- wapply(multic,sim.results,function(x) {
        true.de <- x$simdata$truedeg
        names(true.de) <- rownames(x$simdata$simdata)
        p.matrix <- x$pvalues
        true.de <- true.de[rownames(p.matrix)]
        fdc <- diagplot.ftd(true.de,p.matrix,type="fpc",draw=FALSE)
    })
    avg.fpc <- diagplot.avg.ftd(fpc.obj,draw=draw.fpc)

    x <- 1:top
    aufc <- apply(avg.fpc$avg.ftdr$means[1:top,],2,function(x,i) {
        return(sum(diff(i)*rollmean(x,2)))
    },x)
    weight.aufc <- (sum(aufc)/aufc)/sum(sum(aufc)/aufc)
    return(weight.aufc)
}

#' Create simulated counts using TCC package
#'
#' This function creates simulated RNA-Seq gene expression datasets using the
#' \code{simulateReadCounts} function from the Bioconductor
#' package TCC and it adds simulated annoation elements. For further information
#' please consult the TCC package documentation. Note that the produced data are
#' based in an Arabidopsis dataset.
#'
#' @param ... parameters to the \code{simulateReadCounts} function.
#' @return A list with the following members: \code{simdata} holding the simulated
#' dataset complying with metaseqr requirements, and \code{simparam} holding the
#' simulation parameters (see TCC documentation).
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' dd <- make.sim.data(Ngene=10000,PDEG=0.2,DEG.assign=c(0.9,0.1),
#'   DEG.foldchange=c(5,5),replicates=c(3,3))
#' head(dd$simdata)
#'}
make.sim.data.tcc <- function(...) {
    if (suppressWarnings(!require(TCC)))
        stopwrap("Bioconductor package TCC is required to create simulated data!")
    #tcc <- simulateReadCounts(Ngene=Ngene,PDEG=PDEG,DEG.assign=DEG.assign,
    #    DEG.foldchange=DEG.foldchange,replicates=replicates)
    tcc <- TCC::simulateReadCounts(...)
    n <- nrow(tcc$count)
    # Now we have to simulate annotation
    chromosome <- paste("chr",1+round(20*runif(n)),sep="")
    start <- 1 + round(1e+6*runif(n))
    end <- start + 250 + round(1e+6*runif(n))
    gene_id <- gene_name <- rownames(tcc$count)
    gc_content <- runif(n)
    strand <- sample(c("+","-"),n,replace=TRUE)
    biotype <- sample(paste("biotype",1:10),n,replace=TRUE)
    sim.data <- data.frame(
        chromosome=chromosome,
        start=start,
        end=end,
        gene_id=gene_id,
        gc_content=gc_content,
        strand=strand,
        gene_name=gene_name,
        biotype=biotype
    )
    sim.data <- cbind(sim.data,tcc$count)
    return(list(simdata=sim.data,simparam=tcc$simulation))
}

#' Create simulated counts using the Soneson-Delorenzi method
#'
#' This function creates simulated RNA-Seq gene expression datasets using the
#' method presented in (Soneson and Delorenzi, BMC Bioinformatics, 2013). For the
#' time being, it creates only simulated datasets with two conditions.
#'
#' @param N the number of genes to produce.
#' @param param a named list with negative binomial parameter sets to sample from.
#' The first member is the mean parameter to sample from (\code{mu.hat}} and the
#' second the dispersion (\code{phi.hat}). This list can be created with the
#' \code{\link{estimate.sim.params}} function.
#' @param samples a vector with 2 integers, which are the number of samples for
#' each condition (two conditions currently supported).
#' @param ndeg a vector with 2 integers, which are the number of differentially
#' expressed genes to be produced. The first element is the number of up-regulated
#' genes while the second is the number of down-regulated genes.
#' @param fc.basis the minimum fold-change for deregulation.
#' @param libsize.range a vector with 2 numbers (generally small, see the default),
#' as they are multiplied with \code{libsize.mag}.
#' These numbers control the library sized of the synthetic data to be produced.
#' @param libsize.mag a (big) number to multiply the \code{libsize.range} to
#' produce library sizes.
#' @param model.org the organism from which the real data are derived from. It
#' must be one of the supported organisms (see the main \code{\link{metaseqr}}
#' help page). It is used to sample real values for GC content.
#' @param sim.length.bias a boolean to instruct the simulator to create genes
#' whose read counts is proportional to their length. This is achieved by sorting
#' in increasing order the mean parameter of the negative binomial distribution
#' (and the dispersion according to the mean) which will cause an increasing gene
#' count length with the sampling. The sampled lengths are also sorted so that in
#' the final gene list, shorter genes have less counts as compared to the longer
#' ones. The default is FALSE.
#' @param seed a seed to use with random number generation for reproducibility.
#' @return A named list with two members. The first member (\code{simdata})
#' contains the synthetic dataset 
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # File "bottomly_read_counts.txt" from the ReCount database
#' download.file("http://bowtie-bio.sourceforge.net/recount/countTables/bottomly_count_table.txt",
#'   destfile="~/bottomly_count_table.txt")
#' N <- 10000
#' par.list <- estimate.sim.params("~/bottomly_read_counts.txt")
#' sim <- make.sim.data.sd(N,par.list)
#' synth.data <- sim$simdata
#' true.deg <- which(sim$truedeg!=0)
#'}
make.sim.data.sd <- function(N,param,samples=c(5,5),ndeg=rep(round(0.1*N),2),
    fc.basis=1.5,libsize.range=c(0.7,1.4),libsize.mag=1e+7,model.org=NULL,
    sim.length.bias=FALSE,seed=NULL) {
    if (!is.null(model.org)) {
        model.org <- tolower(model.org)
        check.text.args("model.org",model.org,c("hg18","hg19","mm9","mm10",
            "rno5","dm3","danrer7","pantro4","tair10"),multiarg=FALSE)
        ann <- get.annotation(model.org,"gene")
        real.gc <- as.numeric(ann$gc_content)
        real.start <- as.numeric(ann$start)
        real.end <- as.numeric(ann$end)
        real.strand <- as.character(ann$strand)
    }
    mu.hat <- param$mu.hat
    phi.hat <- param$phi.hat

    if (!is.null(seed)) set.seed(seed)
    if (sim.length.bias) {
        sind <- sort(mu.hat,index.return=TRUE)$ix
        mu.hat <- mu.hat[sind]
        phi.hat <- phi.hat[sind]
        if (length(mu.hat)>=N)
            ii <- sort(sample(1:length(mu.hat),N))
        else
            ii <- sort(sample(1:length(mu.hat),N,replace=TRUE))
    }
    else {
        if (length(mu.hat)>=N)
            ii <- sample(1:length(mu.hat),N)
        else
            ii <- sample(1:length(mu.hat),N,replace=TRUE)
    }

    s1 <- samples[1]
    s2 <- samples[2]
    if (!is.null(seed)) set.seed(seed)
    L1 <- round(libsize.mag*runif(s1,min=libsize.range[1],
        max=libsize.range[2]))
    if (!is.null(seed)) set.seed(2*seed)
    L2 <- round(libsize.mag*runif(s2,min=libsize.range[1],
        max=libsize.range[2]))

    lambda.1 <- do.call("cbind",rep(list(mu.hat[ii]),s1))
    mu.1 <- sweep(lambda.1,2,L1/sum(lambda.1[,1]),"*")
    sim.1 <- matrix(0,N,s1)
    for (j in 1:s1) {
        if (!is.null(seed)) set.seed(seed+j)
        sim.1[,j] <- rnbinom(N,size=1/phi.hat[ii],mu=mu.1[,j])
    }

    v <- numeric(N)
    if (sum(ndeg)>0) {
        if (!is.null(seed)) set.seed(seed)
        i.updown <- sample(1:length(v),sum(ndeg))
        reg.dir <- rep(c(1,-1),c(ndeg[1],ndeg[2]))
        v[i.updown] <- reg.dir
        if (!is.null(seed)) set.seed(seed+19051980)
        lambda.2 <- ((fc.basis + rexp(N))^v)*lambda.1
        mu.2 <- sweep(lambda.2,2,L2/sum(lambda.2[,1]),"*")
        sim.2 <- matrix(0,N,s2)
        for (j in 1:s2)
            sim.2[,j] <- rnbinom(N,size=1/phi.hat[ii],mu=mu.2[,j])
    }
    else {
        if (!is.null(seed)) set.seed(seed+19051980)
        lambda.2 <- lambda.1
        mu.2 <- sweep(lambda.2,2,L2/sum(lambda.2[,1]),"*")
        sim.2 <- matrix(0,N,s2)
        for (j in 1:s2)
            sim.2[,j] <- rnbinom(N,size=1/phi.hat[ii],mu=mu.2[,j])
    }

    # Now we have to simulate annotation
    if (!is.null(seed)) set.seed(seed)
    chromosome <- paste("chr",1+round(20*runif(N)),sep="")
    gene_id <- gene_name <- paste("gene",1:N,sep="_")
    if (!is.null(model.org)) {
        if (!is.null(seed)) set.seed(seed)
        if (length(real.gc)>=N)
            sample.ind <- sample(1:length(real.gc),N)            
        else
            sample.ind <- sample(1:length(real.gc),N,replace=TRUE)
        gc_content <- real.gc[sample.ind]
        start <- real.start[sample.ind]
        end <- real.end[sample.ind]
        strand <- real.strand[sample.ind]
        if (sim.length.bias) {
            lenix <- sort(end-start,index.return=TRUE)$ix
            start <- start[lenix]
            end <- end[lenix]
            gc_content <- gc_content[lenix]
            strand <- strand[lenix]
        }
    }
    else {
        if (!is.null(seed)) set.seed(seed)
        gc_content <- runif(N)
        if (!is.null(seed)) set.seed(seed)
        start <- 1 + round(1e+6*runif(N))
        if (!is.null(seed)) set.seed(seed)
        end <- start + 250 + round(1e+6*runif(N))
        if (!is.null(seed)) set.seed(seed)
        strand <- sample(c("+","-"),N,replace=TRUE)
        if (sim.length.bias) {
            lenix <- sort(end-start,index.return=TRUE)$ix
            start <- start[lenix]
            end <- end[lenix]
            gc_content <- gc_content[lenix]
            strand <- strand[lenix]
        }
    }
    if (!is.null(seed)) set.seed(seed)
    biotype <- sample(paste("biotype",1:10),N,replace=TRUE)
    sim.data <- data.frame(
        chromosome=chromosome,
        start=start,
        end=end,
        gene_id=gene_id,
        gc_content=gc_content,
        strand=strand,
        gene_name=gene_name,
        biotype=biotype
    )
    colnames(sim.1) <- paste("G1_rep",1:s1,sep="")
    colnames(sim.2) <- paste("G2_rep",1:s2,sep="")
    rownames(sim.1) <- rownames(sim.2) <- names(v) <- gene_id

    return(list(simdata=cbind(sim.data,sim.1,sim.2),truedeg=v))
}

#' Estimate negative binomial parameters from real data
#'
#' This function reads a read counts table containing real RNA-Seq data (preferebly
#' with more than 20 samples so as to get as much accurate as possible estimations)
#' and calculates a population of count means and dispersion parameters which can
#' be used to simulate an RNA-Seq dataset with synthetic genes by drawing from a
#' negative binomial distribution. This function works in the same way as described
#' in (Soneson and Delorenzi, BMC Bioinformatics, 2013) and (Robles et al., BMC
#' Genomics, 2012).
#'
#' @param real.counts a text tab-delimited file with real RNA-Seq data. The file
#' should strictly contain a unique gene name (e.g. Ensembl accession) in the
#' first column and all other columns should contain read counts for each gene.
#' Each column must be named with a unique sample identifier. See examples in the
#' ReCount database \link{http://bowtie-bio.sourceforge.net/recount/}.
#' @param libsize.gt a library size below which samples are excluded from parameter
#' estimation (default: 3000000).
#' @param rowmeans.gt a row means (mean counts over samples for each gene) below
#' which genes are excluded from parameter estimation (default: 5).
#' @param eps the tolerance for the convergence of \code{\link{optimize}} function.
#' Defaults to 1e-11.
#' @param restrict.cores in case of parallel optimization, the fraction of the
#' available cores to use.
#' @param seed a seed to use with random number generation for reproducibility.
#' @param draw boolean to determine whether to plot the estimated simulation
#' parameters (mean and dispersion) or not. Defaults to \code{FALSE} (do not draw
#' a mean-dispersion scatterplot).
#' @return A named list with two members: \code{mu.hat} which contains negative
#' binomial mean estimates and \code{phi.hat} which contains dispersion.
#' estimates
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # Dowload locally the file "bottomly_count_table.txt" from the ReCount datbase
#' download.file("http://bowtie-bio.sourceforge.net/recount/countTables/bottomly_count_table.txt",
#'   destfile="~/bottomly_count_table.txt")
#' # Estimate simulation parameters
#' par.list <- estimate.sim.params("~/bottomly_count_table.txt")
#'}
estimate.sim.params <- function(real.counts,libsize.gt=3e+6,rowmeans.gt=5,
    eps=1e-11,restrict.cores=0.8,seed=42,draw=FALSE) {
    multic <- check.parallel(restrict.cores)
    if (is.data.frame(real.counts))
        mat <- as.matrix(real.counts)
    else if (is.matrix(real.counts))
        mat <- real.counts
    else if (file.exists(real.counts)) {
        real.data <- read.delim(real.counts,row.names=1)
        mat <- as.matrix(real.data)
    }
    else
        stopwrap("The input count data must be either a file, a matrix or a ",
            "data frame!")
    
    low.lib <- which(apply(mat,2,sum)<libsize.gt)
    if (length(low.lib)==ncol(mat))
        stopwrap("Cannot estimate simulation parameters as the library sizes ",
            "are too small! Try lowering the value of the libsize.gt ",
            "parameter...")
    if (length(low.lib)>0)
        mat <- mat[,-low.lib]
    disp("Downsampling counts...")
    dmat <- downsample.counts(mat,seed)
    low.co <- which(apply(dmat,1,
        function(x) if (mean(x)<5) TRUE else FALSE))
    if (length(low.co)>0)
        dmat <- dmat[-low.co,]
    mu.hat <- apply(dmat,1,mean)
    disp("Estimating initial dispersion population...")
    phi.est <- apply(dmat,1,function(x) {
        m <- mean(x)
        v <- var(x)
        phi <- (v-m)/m^2
        return(phi)
    })
    phi.ind <- which(phi.est>0)
    phi.est <- phi.est[phi.ind]
    dmat <- dmat[phi.ind,]
    disp("Estimating dispersions using log-likelihood...\n")
    init <- wapply(multic,seq_along(1:nrow(dmat)),function(i,d,p) {
        list(y=d[i,],h=p[i])
    },dmat,phi.est)
    phi.hat <- unlist(wapply(multic,init,function(x,eps) {
        optimize(mlfo,c(x$h-1e-2,x$h+1e-2),y=x$y,tol=eps)$minimum
    },eps))
    if (draw) {
        dev.new()
        plot(log10(mu.hat[phi.ind]),log10(phi.hat),col="blue",pch=20,cex=0.5,
            xlab="",ylab="")
        title(xlab="mean",ylab="dispesion",font=2,cex=0.9)
        grid()
    }
    return(list(mu.hat=mu.hat[phi.ind],phi.hat=phi.hat))
}

#' Downsample read counts
#'
#' This function downsamples the library sizes of a read counts table to the lowest
#' library size, according to the methdology used in  (Soneson and Delorenzi,
#' BMC Bioinformatics, 2013).
#'
#' @param counts the read counts table which is subjected to downsampling.
#' @param seed random seed for reproducible downsampling.
#' @return The downsampled counts matrix.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # Dowload locally the file "bottomly_read_counts.txt" from
#' # the ReCount database
#' download.file(paste("http://bowtie-bio.sourceforge.net/",
#'   "recount/countTables/bottomly_count_table.txt",sep=""),
#'  destfile="~/bottomly_count_table.txt")
#' M <- as.matrix(read.delim("~/bottomly_count_table.txt",row.names=1))
#' D <- downsample.counts(M)
#'}
downsample.counts <- function(counts,seed=42) {
    lib.sizes <- apply(counts,2,sum)
    target.size <- min(lib.sizes)
    to.remove <- lib.sizes-target.size
    ii <- which(to.remove>0)
    dcounts <- counts
    for (i in ii) {
        tmp <- round(to.remove[i]*(counts[,i]/sum(counts[,i])))
        victim.size <- sum(tmp)
        if (victim.size>to.remove[i]) {
            dif <- victim.size - to.remove[i]
            #victims <- sample(1:length(tmp),dif)
            victims <- sort(tmp,decreasing=TRUE,index.return=TRUE)$ix[1:dif]
            tmp[victims] <- tmp[victims] - 1
        }
        else if (victim.size<to.remove[i]) {
            dif <- to.remove[i] - victim.size
            #victims <- sample(1:length(tmp),dif)
            victims <- sort(tmp,decreasing=TRUE,index.return=TRUE)$ix[1:dif]
            tmp[victims] <- tmp[victims] + 1
        }
        dcounts[,i] <- dcounts[,i] - tmp
    }
    return(dcounts)
}

#' MLE dispersion estimate
#'
#' MLE function used to estimate negative binomial dispersions from real RNA-Seq
#' data, as in (Soneson and Delorenzi, BMC Bioinformatics, 2013) and (Robles et al.,
#' BMC Genomics, 2012). Internal use.
#'
#' @param phi the parameter to be optimized.
#' @param y count samples used to perform the optimization.
#' @return objective function value.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # Not yet available
#'}
mlfo <- function(phi,y) {
    N <- length(y)
    mu <- mean(y)
    -(sum(lgamma(y+1/phi)) - N*lgamma(1/phi) - sum(lgamma(y+1)) + 
        sum(y*log(mu*phi/(1+mu*phi))) - (N/phi)*log(1+mu*phi))
}

#' Create counts matrix permutations
#'
#' This function creates a permuted read counts matrix based on the \code{contrast}
#' argument (to define new virtual contrasts of the same number) and on the
#' \code{sample.list} to derive the number of samples for each virtual condition.
#' It is a helper for the \code{\link{meta.perm}} function.
#'
#' @param counts the gene read counts matrix.
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param contrast the contrasts vector. See the main \code{\link{metaseqr}} help
#' page.
#' @param repl the same as the replace argument in \code{\link{sample}} function.
#' @return A list with three members: the matrix of permuted per sample read counts,
#' the virtual sample list and the virtual contrast to be used with the \code{stat.*}
#' functions.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' data("mm9.gene.data",package="metaseqR")
#' per <- make.permutation(mm9.gene.counts,sample.list.mm9,"e14.5_vs_adult_8_weeks")
#'}
make.permutation <- function(counts,sample.list,contrast,repl=FALSE) {
    cnts <- strsplit(contrast,"_vs_")[[1]]
    virtual.contrast <- paste(paste("VirtCond",1:length(cnts),sep=""),
        collapse="_vs_")
    virtual.sample.list <- vector("list",length(sample.list))
    names(virtual.sample.list) <- paste("VirtCond",1:length(sample.list),
        sep="")
    # Avoid the extreme case of returning a vector with all samples the same
    if (repl) {
        resample <- rep(1,ncol(counts))
        while(length(unique(resample))==1)
            resample <- sample(1:ncol(counts),ncol(counts),replace=repl)
    }
    else
        resample <- sample(1:ncol(counts),ncol(counts),replace=repl)
    virtual.counts <- counts[,resample]
    samples <- paste("VirtSamp",1:ncol(counts),sep="")
    colnames(virtual.counts) <- samples
    nsample <- sapply(sample.list,length)
    virtual.samples <- split(samples,rep(1:length(nsample),nsample))
    names(virtual.samples) <- names(virtual.sample.list)
    for (n in names(virtual.sample.list))
        virtual.sample.list[[n]] <- virtual.samples[[n]]
    return(list(counts=virtual.counts,sample.list=virtual.sample.list,
        contrast=virtual.contrast))
}

#' Calculate the ratio TP/(FP+FN)
#'
#' This function calculates the ratio of True Positives to the sum of False
#' Positives and False Negatives given a matrix of p-values (one for each
#' statistical test used) and a vector of ground truth (DE or non-DE). This
#' function serves as a method evaluation helper.
#'
#' @param truth the ground truth differential expression vector. It should contain
#' only zero and non-zero elements, with zero denoting non-differentially expressed
#' genes and non-zero, differentially expressed genes. Such a vector can be obtained
#' for example by using the \code{\link{make.sim.data.sd}} function, which creates
#' simulated RNA-Seq read counts based on real data. It MUST be named with gene
#' names, the same as in \code{p}.
#' @param p a p-value matrix whose rows correspond to each element in the
#' \code{truth} vector. If the matrix has a \code{colnames} attribute, a legend
#' will be added to the plot using these names, else a set of column names will
#' be auto-generated. \code{p} can also be a list or a data frame. In any case,
#' each row (or element) MUST be named with gene names (the same as in \code{truth}).
#' @param sig a significance level (0 < \code{sig} <=1).
#' @return A named list with two members. The first member is a data frame with
#' the numbers used to calculate the TP/(FP+FN) ratio and the second member is
#' the ratio TP/(FP+FN) for each statistical test.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p1 <- 0.001*matrix(runif(300),100,3)
#' p2 <- matrix(runif(300),100,3)
#' p <- rbind(p1,p2)
#' rownames(p) <- paste("gene",1:200,sep="_")
#' colnames(p) <- paste("method",1:3,sep="_")
#' truth <- c(rep(1,40),rep(-1,40),rep(0,10),rep(1,10),rep(2,10),rep(0,80))
#' names(truth) <- rownames(p)
#' otr <- calc.otr(truth,p)
#'}
calc.otr <- function(truth,p,sig=0.05) {
    if (is.list(p))
        pmat <- do.call("cbind",p)
    else if (is.data.frame(p))
        pmat <- as.matrix(p)
    else if (is.matrix(p))
        pmat <- p
    if (is.null(colnames(pmat)))
        colnames(pmat) <- paste("p",1:ncol(pmat),sep="_")

    sig.genes <- true.isects <- missed <- vector("list",ncol(pmat))
    names(sig.genes) <- names(true.isects) <- names(missed) <- colnames(pmat)
    for (n in colnames(pmat)) {
        sig.genes[[n]] <- names(which(pmat[,n]<sig))
        true.isects[[n]] <- intersect(sig.genes[[n]],names(which(truth!=0)))
        missed[[n]] <- setdiff(names(which(truth!=0)),true.isects[[n]])
    }
    result <- data.frame(
        P=sapply(sig.genes,length),
        TP=sapply(true.isects,length),
        FN=sapply(missed,length)
    )
    result$FP <- result$P - result$TP
    otr <- result$TP/(result$FP+result$FN)
    names(otr) <- rownames(result)
    return(list(result=result,otr=otr))
}

#' Calculate the F1-score
#'
#' This function calculates the F1 score (2*(precision*recall/precision+racall)
#' or 2*TP/(2*TP+FP+FN) given a matrix of p-values (one for each statistical test
#' used) and a vector of ground truth (DE or non-DE). This function serves as a
#' method evaluation helper.
#'
#' @param truth the ground truth differential expression vector. It should contain
#' only zero and non-zero elements, with zero denoting non-differentially expressed
#' genes and non-zero, differentially expressed genes. Such a vector can be obtained
#' for example by using the \code{\link{make.sim.data.sd}} function, which creates
#' simulated RNA-Seq read counts based on real data. It MUST be named with gene
#' names, the same as in \code{p}.
#' @param p a p-value matrix whose rows correspond to each element in the
#' \code{truth} vector. If the matrix has a \code{colnames} attribute, a legend
#' will be added to the plot using these names, else a set of column names will
#' be auto-generated. \code{p} can also be a list or a data frame. In any case,
#' each row (or element) MUST be named with gene names (the same as in \code{truth}).
#' @param sig a significance level (0 < \code{sig} <=1).
#' @return A named list with two members. The first member is a data frame with
#' the numbers used to calculate the F1-score and the second member is the
#' F1-score for each statistical test.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' p1 <- 0.001*matrix(runif(300),100,3)
#' p2 <- matrix(runif(300),100,3)
#' p <- rbind(p1,p2)
#' rownames(p) <- paste("gene",1:200,sep="_")
#' colnames(p) <- paste("method",1:3,sep="_")
#' truth <- c(rep(1,40),rep(-1,40),rep(0,10),rep(1,10),rep(2,10),rep(0,80))
#' names(truth) <- rownames(p)
#' f1 <- calc.f1score(truth,p)
#'}
calc.f1score <- function(truth,p,sig=0.05) {
    if (is.list(p))
        pmat <- do.call("cbind",p)
    else if (is.data.frame(p))
        pmat <- as.matrix(p)
    else if (is.matrix(p))
        pmat <- p
    if (is.null(colnames(pmat)))
        colnames(pmat) <- paste("p",1:ncol(pmat),sep="_")

    sig.genes <- true.isects <- missed <- vector("list",ncol(pmat))
    names(sig.genes) <- names(true.isects) <- names(missed) <- colnames(pmat)
    for (n in colnames(pmat)) {
        sig.genes[[n]] <- names(which(pmat[,n]<sig))
        true.isects[[n]] <- intersect(sig.genes[[n]],names(which(truth!=0)))
        missed[[n]] <- setdiff(names(which(truth!=0)),true.isects[[n]])
    }
    result <- data.frame(
        P=sapply(sig.genes,length),
        TP=sapply(true.isects,length),
        FN=sapply(missed,length)
    )
    result$FP <- result$P - result$TP
    f1 <- 2*result$TP/(2*result$TP+result$FP+result$FN)
    names(f1) <- rownames(result)
    return(list(result=result,f1=f1))
}
