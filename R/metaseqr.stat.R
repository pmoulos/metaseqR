#' Statistical testing with DESeq
#'
#' This function is a wrapper over DESeq statistical testing. It accepts a matrix
#' of normalized gene counts or an S4 object specific to each normalization
#' algorithm supported by metaseqR.
#'
#' @param object a matrix or an object specific to each normalization algorithm
#' supported by metaseqR, containing normalized counts. Apart from matrix (also
#' for NOISeq), the object can be a SeqExpressionSet (EDASeq), CountDataSet (DESeq)
#' or DGEList (edgeR).
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param contrast.list a named structured list of contrasts as returned by
#' \code{\link{make.contrast.list}} or just the vector of contrasts as defined
#' in the main help page of \code{\link{metaseqr}}.
#' @param stat.args a list of DESeq statistical algorithm parameters. See the
#' result of \code{get.defaults("statistics",} \code{"deseq")} for an example and
#' how you can modify it. It is not required when the input object is already a
#' CountDataSet from DESeq normalization
#' as the dispersions are already estimated.
#' @return A named list of p-values, whose names are the names of the contrasts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' contrast <- "A_vs_B"
#' norm.data.matrix <- normalize.deseq(data.matrix,sample.list)
#' p <- stat.deseq(norm.data.matrix,sample.list,contrast)
#'}
stat.deseq <- function(object,sample.list,contrast.list=NULL,stat.args=NULL) {
    #if (is.null(norm.args) && class(object)=="DGEList")
    #    norm.args <- get.defaults("normalization","edger")
    if (is.null(stat.args) && class(object)!="CountDataSet")
        stat.args <- get.defaults("statistics","deseq")
    if (is.null(contrast.list))
        contrast.list <- make.contrast.list(paste(names(sample.list)[1:2],
            collapse="_vs_"),sample.list)
    if (!is.list(contrast.list))
        contrast.list <- make.contrast.list(contrast.list,sample.list)
    classes <- as.class.vector(sample.list)
    the.design <- data.frame(condition=classes,row.names=colnames(object))
    p <- vector("list",length(contrast.list))
    names(p) <- names(contrast.list)
     # Check if there is no replication anywhere
    if (all(sapply(sample.list,function(x) ifelse(length(x)==1,TRUE,
        FALSE)))) {
        warnwrap("No replication detected! There is a possibility that ",
            "DESeq will fail to estimate dispersions...")
        method.disp <- "blind"
        sharingMode.disp <- "fit-only"
        fitType.disp <- "local"
    }
    else {
        method.disp <- stat.args$method
        sharingMode.disp <- stat.args$sharingMode
        fitType.disp <- stat.args$fitType
    }
    switch(class(object),
        CountDataSet = { # Has been normalized with DESeq
            cds <- object
            cds <- estimateDispersions(cds,method=method.disp,
                sharingMode=sharingMode.disp,fitType=fitType.disp)
        },
        DGEList = { # Has been normalized with edgeR
            # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
            scl <- object$samples$lib.size * object$samples$norm.factors
            cds <- newCountDataSet(round(t(t(object$counts)/scl)*mean(scl)),
                the.design$condition)
            sizeFactors(cds) <- rep(1,ncol(cds))
            cds <- estimateDispersions(cds,method=method.disp,
                sharingMode=sharingMode.disp)
        },
        matrix = { # Has been normalized with EDASeq or NOISeq
            cds <- newCountDataSet(object,the.design$condition)
            sizeFactors(cds) <- rep(1,ncol(cds))
            cds <- estimateDispersions(cds,method=method.disp,
                sharingMode=sharingMode.disp)
        },
        list = { # Has been normalized with NBPSeq and main method was "nbpseq"
            cds <- newCountDataSet(as.matrix(round(sweep(object$counts,2,
                object$norm.factors,"*"))),the.design$condition)
            sizeFactors(cds) <- rep(1,ncol(cds))
            cds <- estimateDispersions(cds,method=method.disp,
                sharingMode=sharingMode.disp)
        },
        nbp = { # Has been normalized with NBPSeq and main method was "nbsmyth"...
            cds <- newCountDataSet(as.matrix(round(object$pseudo.counts)),
                the.design$condition)
            sizeFactors(cds) <- rep(1,ncol(cds))
            cds <- estimateDispersions(cds,method=method.disp,
                sharingMode=sharingMode.disp)
        }
    )
    for (con.name in names(contrast.list)) {
        disp("  Contrast: ", con.name)
        con <- contrast.list[[con.name]]
        cons <- unique(unlist(con))
        if (length(con)==2) {
            res <- nbinomTest(cds,cons[1],cons[2])
            p[[con.name]] <- res$pval
        }
        else {
            #cind <- match(cons,the.design$condition)
            #if (any(is.na(cind)))
            #    cind <- cind[which(!is.na(cind))]
            cc <- names(unlist(con))
            cds.tmp <- cds[,cc]
            fit0 <- fitNbinomGLMs(cds.tmp,count~1)
            fit1 <- fitNbinomGLMs(cds.tmp,count~condition)
            p[[con.name]] <- nbinomGLMTest(fit1,fit0)
        }
        names(p[[con.name]]) <- rownames(object)
        p[[con.name]][which(is.na(p[[con.name]]))] <- 1
    }
    return(p)
}

#' Statistical testing with edgeR
#'
#' This function is a wrapper over edgeR statistical testing. It accepts a matrix
#' of normalized gene counts or an S4 object specific to each normalization
#' algorithm supported by metaseqR.
#'
#' @param object a matrix or an object specific to each normalization algorithm
#' supported by metaseqr, containing normalized counts. Apart from matrix (also
#' for NOISeq), the object can be a SeqExpressionSet (EDASeq), CountDataSet (DESeq)
#' or DGEList (edgeR).
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param contrast.list a named structured list of contrasts as returned by
#' \code{\link{make.contrast.list}} or just the vector of contrasts as defined in
#' the main help page of \code{\link{metaseqr}}.
#' @param stat.args a list of edgeR statistical algorithm parameters. See the
#' result of \code{get.defaults("statistics",} \code{"edger")} for an example and
#' how you can modify it.
#' @return A named list of p-values, whose names are the names of the contrasts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' contrast <- "A_vs_B"
#' norm.data.matrix <- normalize.edger(data.matrix,sample.list)
#' p <- stat.edger(norm.data.matrix,sample.list,contrast)
#'}
stat.edger <- function(object,sample.list,contrast.list=NULL,stat.args=NULL) {
    if (is.null(stat.args))
        stat.args <- get.defaults("statistics","edger")
    if (is.null(contrast.list))
        contrast.list <- make.contrast.list(paste(names(sample.list)[1:2],
            collapse="_vs_"),sample.list)
    if (!is.list(contrast.list))
        contrast.list <- make.contrast.list(contrast.list,sample.list)
    classes <- as.class.vector(sample.list)
    p <- vector("list",length(contrast.list))
    names(p) <- names(contrast.list)
    switch(class(object),
        CountDataSet = { # Has been normalized with DESeq
            dge <- DGEList(counts=counts(object,normalized=TRUE),group=classes)
        },
        DGEList = { # Has been normalized with edgeR
            dge <- object    
        },
        matrix = { # Has been normalized with EDASeq or NOISeq
            dge <- DGEList(object,group=classes)
        },
        list = { # Has been normalized with NBPSeq and main method was "nbpseq"
            dge <- DGEList(counts=as.matrix(round(sweep(object$counts,2,
                object$norm.factors,"*"))),group=classes)
        },
        nbp = { # Has been normalized with NBPSeq and main method was "nbsmyth"
            dge <- DGEList(counts=as.matrix(round(object$pseudo.counts)),
                group=classes)
        }
    )
    # Dispersion estimate step
    # Check if there is no replication anywhere
    repli = TRUE
    if (all(sapply(sample.list,function(x) ifelse(length(x)==1,TRUE,
        FALSE)))) {
        warnwrap("No replication when testing with edgeR! Consider using ",
            "another statistical test or just performing empirical analysis. ",
            "Setting to 0.2...")
        repli <- FALSE
        bcv <- 0.2
    }
    if (repli) {
        if (stat.args$main.method=="classic") {
            dge <- estimateCommonDisp(dge,rowsum.filter=stat.args$rowsum.filter)
            dge <- estimateTagwiseDisp(dge,prior.df=stat.args$prior.df,
                trend=stat.args$trend,span=stat.args$span,
                    method=stat.args$tag.method,
                grid.length=stat.args$grid.length,
                    grid.range=stat.args$grid.range)
        }
        else if (stat.args$main.method=="glm") {
            design <- model.matrix(~0+classes,data=dge$samples)
            dge <- estimateGLMCommonDisp(dge,design=design,
                offset=stat.args$offset,
                method=stat.args$glm.method,subset=stat.args$subset,
                AveLogCPM=stat.args$AveLogCPM)
            dge <- estimateGLMTrendedDisp(dge,design=design,
                offset=stat.args$offset,
                method=stat.args$trend.method,AveLogCPM=stat.args$AveLogCPM)
            dge <- estimateGLMTagwiseDisp(dge,design=design,
                offset=stat.args$offset,
                dispersion=stat.args$dispersion,prior.df=stat.args$prior.df,
                span=stat.args$span,AveLogCPM=stat.args$AveLogCPM)
        }
    }
    # Actual statistical test
    for (con.name in names(contrast.list))
    {
        disp("  Contrast: ", con.name)
        con <- contrast.list[[con.name]]
        if (length(con)==2) {
            if (repli) {
                if (stat.args$main.method=="classic") {
                    res <- exactTest(dge,pair=unique(unlist(con)))
                }
                else if (stat.args$main.method=="glm") {
                    s <- unlist(con)
                    us <- unique(s)
                    ms <- match(names(s),rownames(dge$samples))
                    if (any(is.na(ms)))
                        ms <- ms[which(!is.na(ms))]
                    design <- model.matrix(~0+s,data=dge$samples[ms,])
                    colnames(design) <- us
                    #fit <- glmFit(dge[,ms],design=design,
                    #    offset=stat.args$offset,
                    #    weights=stat.args$weights,lib.size=stat.args$lib.size,
                    #    prior.count=stat.args$prior.count,
                    #    start=stat.args$start,method=stat.args$method)
                    fit <- glmFit(dge[,ms],design=design,
                        prior.count=stat.args$prior.count,
                        start=stat.args$start,method=stat.args$method)
                    co <- makeContrasts(paste(us[2],us[1],sep="-"),
                        levels=design)
                    lrt <- glmLRT(fit,contrast=co)
                    res <- topTags(lrt,n=nrow(dge))
                }
            }
            else {
                if (stat.args$main.method=="classic") {
                    res <- exactTest(dge,pair=unique(unlist(con)),
                        dispersion=bcv^2)
                }
                else if (stat.args$main.method=="glm") {
                    s <- unlist(con)
                    us <- unique(s)
                    ms <- match(names(s),rownames(dge$samples))
                    if (any(is.na(ms)))
                        ms <- ms[which(!is.na(ms))]
                    design <- model.matrix(~0+s,data=dge$samples[ms,])
                    colnames(design) <- us
                    #fit <- glmFit(dge[,ms],design=design,
                    #    offset=stat.args$offset,weights=stat.args$weights,
                    #    lib.size=stat.args$lib.size,
                    #    prior.count=stat.args$prior.count,
                    #    start=stat.args$start,
                    #    method=stat.args$method,dispersion=bcv^2)
                    fit <- glmFit(dge[,ms],design=design,dispersion=bcv^2,
                        prior.count=stat.args$prior.count,start=stat.args$start,
                        method=stat.args$method)
                    co <- makeContrasts(paste(us[2],us[1],sep="-"),
                        levels=design)
                    lrt <- glmLRT(fit,contrast=co)
                    res <- topTags(lrt,n=nrow(dge))
                }
            }
        }
        else { # GLM only
            s <- unlist(con)
            us <- unique(s)
            #design <- model.matrix(~0+s,data=dge$samples) # Ouch!
            ms <- match(names(s),rownames(dge$samples))
            if (any(is.na(ms)))
                ms <- ms[which(!is.na(ms))]
            design <- model.matrix(~s,data=dge$samples[ms,])
            if (repli)
                #fit <- glmFit(dge[,ms],design=design,offset=stat.args$offset,
                #    weights=stat.args$weights,lib.size=stat.args$lib.size,
                #    prior.count=stat.args$prior.count,start=stat.args$start,
                #    method=stat.args$method)
                fit <- glmFit(dge[,ms],design=design,start=stat.args$start,
                    prior.count=stat.args$prior.count)
            else
                #fit <- glmFit(dge[,ms],design=design,offset=stat.args$offset,
                #    weights=stat.args$weights,lib.size=stat.args$lib.size,
                #    prior.count=stat.args$prior.count,start=stat.args$start,
                #    method=stat.args$method,dispersion=bcv^2)
                #    dispersion=bcv^2)
                fit <- glmFit(dge[,ms],design=design,dispersion=bcv^2,
                    prior.count=stat.args$prior.count,start=stat.args$start)
                #lrt <- glmLRT(fit,coef=2:ncol(fit$design))
                lrt <- glmLRT(fit,coef=2:ncol(fit$design))
                res <- topTags(lrt,n=nrow(dge))
        }
        p[[con.name]] <- res$table[,"PValue"]
        names(p[[con.name]]) <- rownames(res$table)
        p[[con.name]] <- p[[con.name]][rownames(dge)]
        p[[con.name]][which(is.na(p[[con.name]]))] <- 1
    }
    return(p)
}

#' Statistical testing with limma
#'
#' This function is a wrapper over limma statistical testing. It accepts a matrix
#' of normalized gene counts or an S4 object specific to each normalization
#' algorithm supported by metaseqR.
#'
#' @param object a matrix or an object specific to each normalization algorithm
#' supported by metaseqr, containing normalized counts. Apart from matrix (also
#' for NOISeq), the object can be a SeqExpressionSet (EDASeq), CountDataSet (DESeq)
#' or DGEList (edgeR).
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param contrast.list a named structured list of contrasts as returned by
#' \code{\link{make.contrast.list}} or just the vector of contrasts as defined in
#' the main help page of \code{\link{metaseqr}}.
#' @param stat.args a list of edgeR statistical algorithm parameters. See the
#' result of \code{get.defaults("statistics",} \code{"limma")} for an example and
#' how you can modify it.
#' @return A named list of p-values, whose names are the names of the contrasts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' contrast <- "A_vs_B"
#' norm.data.matrix <- normalize.edger(data.matrix,sample.list)
#' p <- stat.limma(norm.data.matrix,sample.list,contrast)
#'}
stat.limma <- function(object,sample.list,contrast.list=NULL,stat.args=NULL) {
    if (is.null(stat.args))
        stat.args <- get.defaults("statistics","limma")
    if (is.null(contrast.list))
        contrast.list <- make.contrast.list(paste(names(sample.list)[1:2],
            collapse="_vs_"),sample.list)
    if (!is.list(contrast.list))
        contrast.list <- make.contrast.list(contrast.list,sample.list)
    classes <- as.class.vector(sample.list)
    p <- vector("list",length(contrast.list))
    names(p) <- names(contrast.list)
    switch(class(object),
        CountDataSet = { # Has been normalized with DESeq
            dge <- DGEList(counts=counts(object,normalized=TRUE),group=classes)
        },
        DGEList = { # Has been normalized with edgeR
            dge <- object
        },
        matrix = { # Has been normalized with EDASeq or NOISeq
            dge <- DGEList(object,group=classes)
        },
        list = { # Has been normalized with NBPSeq and main method was "nbpseq"
            dge <- DGEList(counts=as.matrix(round(sweep(object$counts,2,
                object$norm.factors,"*"))),group=classes)
        },
        nbp = { # Has been normalized with NBPSeq and main method was "nbsmyth"
            dge <- DGEList(counts=as.matrix(round(object$pseudo.counts)),
                group=classes)
        }
    )
    for (con.name in names(contrast.list))
    {
        disp("  Contrast: ", con.name)
        con <- contrast.list[[con.name]]
        s <- unlist(con)
        us <- unique(s)
        ms <- match(names(s),rownames(dge$samples))
        if (any(is.na(ms)))
            ms <- ms[which(!is.na(ms))]
        if (length(con)==2) {
            design <- model.matrix(~0+s,data=dge$samples[ms,])
            colnames(design) <- us
            vom <- voom(dge[,ms],design,
                normalize.method=stat.args$normalize.method)
            fit <- lmFit(vom,design)
            fit <- eBayes(fit)
            co <- makeContrasts(contrasts=paste(us[2],us[1],sep="-"),
                levels=design)
            fit <- eBayes(contrasts.fit(fit,co))
            p[[con.name]] <- fit$p.value[,1]
        }
        else {
            design <- model.matrix(~s,data=dge$samples[ms,])
            vom <- voom(dge[,ms],design,
                normalize.method=stat.args$normalize.method)
            fit <- lmFit(vom,design)
            fit <- eBayes(fit)
            res <- topTable(fit,coef=2:ncol(fit$design),number=nrow(vom))
            p[[con.name]] <- res[,"P.Value"]
            names(p[[con.name]]) <- rownames(res)
            p[[con.name]] <- p[[con.name]][rownames(dge)]
        }
        p[[con.name]][which(is.na(p[[con.name]]))] <- 1
    }
    return(p)
}

#' Statistical testing with NOISeq
#'
#' This function is a wrapper over NOISeq statistical testing. It accepts a matrix
#' of normalized gene counts or an S4 object specific to each normalization
#' algorithm supported by metaseqR.
#'
#' @param object a matrix or an object specific to each normalization algorithm
#' supported by metaseqr, containing normalized counts. Apart from matrix (also
#' for NOISeq), the object can be a SeqExpressionSet (EDASeq), CountDataSet (DESeq)
#' or DGEList (edgeR).
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param contrast.list a named structured list of contrasts as returned by
#' \code{\link{make.contrast.list}} or just the vector of contrasts as defined in
#' the main help page of \code{\link{metaseqr}}.
#' @param stat.args a list of edgeR statistical algorithm parameters. See the
#' result of \code{get.defaults("statistics",} \code{"noiseq")} for an example
#' and how you can modify it.
#' @param gene.data an optional annotation data frame (such the ones produced by
#' \code{get.annotation} which contains the GC content for each gene and from
#' which the gene lengths can be inferred by chromosome coordinates.
#' @param log.offset a number to be added to each element of data matrix in order
#' to avoid Infinity on log type data transformations.
#' @return A named list of NOISeq q-values, whose names are the names of the
#' contrasts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' contrast <- "A_vs_B"
#' lengths <- round(1000*runif(nrow(data.matrix)))
#' starts <- round(1000*runif(nrow(data.matrix)))
#' ends <- starts + lengths
#' gc=runif(nrow(data.matrix)),
#' gene.data <- data.frame(
#'   chromosome=c(rep("chr1",nrow(data.matrix)/2),rep("chr2",nrow(data.matrix)/2)),
#'   start=starts,end=ends,gene_id=rownames(data.matrix),gc_content=gc
#' )
#' norm.data.matrix <- normalize.noiseq(data.matrix,sample.list,gene.data)
#' p <- stat.noiseq(norm.data.matrix,sample.list,contrast,gene.data=gene.data)
#'}
stat.noiseq <- function(object,sample.list,contrast.list=NULL,stat.args=NULL,
    gene.data=NULL,log.offset=1) {
    #if (is.null(norm.args) && class(object)=="DGEList")
    #    norm.args <- get.defaults("normalization","edger")
    if (is.null(stat.args))
        stat.args <- get.defaults("statistics","noiseq")
    if (is.null(contrast.list))
        contrast.list <- make.contrast.list(paste(names(sample.list)[1:2],
            collapse="_vs_"),sample.list)
    if (!is.list(contrast.list))
        contrast.list <- make.contrast.list(contrast.list,sample.list)
    if (is.null(gene.data)) {
        gc.content <- NULL
        chromosome <- NULL
        biotype <- NULL
        gene.length <- NULL
    }
    else {
        gc.content <- gene.data$gc_content
        biotype <- as.character(gene.data$biotype)
        names(gc.content) <- names(biotype) <- rownames(gene.data)
        if (is.null(attr(gene.data,"gene.length")))
            gene.length <- NULL
        else {
            gene.length <- attr(gene.data,"gene.length")
            names(gene.length) <- rownames(gene.data)
        }
    }
    classes <- as.class.vector(sample.list)
    p <- vector("list",length(contrast.list))
    names(p) <- names(contrast.list)
    switch(class(object),
        CountDataSet = { # Has been normalized with DESeq
            ns.obj <- NOISeq::readData(
                data=counts(object,normalized=TRUE),
                length=gene.length,
                gc=gc.content,
                chromosome=gene.data[,1:3],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        DGEList = { # Has been normalized with edgeR
            # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
            scl <- object$samples$lib.size * object$samples$norm.factors
            dm <- round(t(t(object$counts)/scl)*mean(scl))            
            ns.obj <- NOISeq::readData(
                data=dm,
                length=gene.length,
                gc=gc.content,
                chromosome=gene.data[,1:3],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        ExpressionSet = { # Has been normalized with NOISeq
            ns.obj <- object
        },
        matrix = { # Has been normalized with EDASeq
            ns.obj <- NOISeq::readData(
                data=object,
                length=gene.length,
                gc=gc.content,
                chromosome=gene.data[,1:3],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        list = { # Has been normalized with NBPSeq and main method was "nbpseq"
            ns.obj <- NOISeq::readData(
                data=as.matrix(round(sweep(object$counts,2,
                    object$norm.factors,"*"))),
                length=gene.length,
                gc=gc.content,
                chromosome=gene.data[,1:3],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        },
        nbp = { # Has been normalized with NBPSeq and main method was "nbsmyth"
            ns.obj <- NOISeq::readData(
                data=as.matrix(round(object$pseudo.counts)),
                length=gene.length,
                gc=gc.content,
                chromosome=gene.data[,1:3],
                factors=data.frame(class=classes),
                biotype=biotype
            )
        }
    )
    for (con.name in names(contrast.list)) {
        disp("  Contrast: ", con.name)
        con <- contrast.list[[con.name]]
        if (length(con)==2) {
            stat.args$conditions=unique(unlist(con))
            if (any(sapply(sample.list,function(x) ifelse(length(x)==1,
                TRUE,FALSE))))
                # At least one condition does not have replicates
                stat.args$replicates <- "no" 
            if (stat.args$replicates %in% c("technical","no"))
                res <- noiseq(ns.obj,k=log.offset,norm="n",
                    replicates=stat.args$replicates,factor=stat.args$factor,
                    conditions=stat.args$conditions,pnr=stat.args$pnr,
                    nss=stat.args$nss,v=stat.args$v,lc=stat.args$lc)
            else
                res <- noiseqbio(ns.obj,k=log.offset,norm="n",
                    nclust=stat.args$nclust,factor=stat.args$factor,
                        lc=stat.args$lc,conditions=stat.args$conditions,
                        r=stat.args$r,adj=stat.args$adj,a0per=stat.args$a0per,
                        cpm=stat.args$cpm,
                    filter=stat.args$filter,depth=stat.args$depth,
                    cv.cutoff=stat.args$cv.cutoff)
            # Beware! This is not the classical p-value!
            p[[con.name]] <- 1 - res@results[[1]]$prob
        }
        else {
            warnwrap(paste("NOISeq differential expression algorithm does not ",
                    "support multi-factor designs (with more than two ",
                    "conditions to be compared)! Switching to DESeq for this ",
                    "comparison:",con.name))
            M <- assayData(ns.obj)$exprs
            cds <- newCountDataSet(round(M),data.frame(condition=unlist(con),
                row.names=names(unlist(con))))
            sizeFactors(cds) <- rep(1,ncol(cds))
            cds <- estimateDispersions(cds,method="blind",
                sharingMode="fit-only")
            fit0 <- fitNbinomGLMs(cds,count~1)
            fit1 <- fitNbinomGLMs(cds,count~condition)
            p[[con.name]] <- nbinomGLMTest(fit1,fit0)
        }
        names(p[[con.name]]) <- rownames(ns.obj)
        p[[con.name]][which(is.na(p[[con.name]]))] <- 1
    }
    return(p)
}

#' Statistical testing with baySeq
#'
#' This function is a wrapper over baySeq statistical testing. It accepts a matrix
#' of normalized gene counts or an S4 object specific to each normalization
#' algorithm supported by metaseqR.
#'
#' @param object a matrix or an object specific to each normalization algorithm
#' supported by metaseqr, containing normalized counts. Apart from matrix (also
#' for NOISeq), the object can be a SeqExpressionSet (EDASeq), CountDataSet
#' (DESeq) or DGEList (edgeR).
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param contrast.list a named structured list of contrasts as returned by
#' \code{\link{make.contrast.list}} or just the vector of contrasts as defined in
#' the main help page of \code{\link{metaseqr}}.
#' @param stat.args a list of edgeR statistical algorithm parameters. See the
#' result of \code{get.defaults("statistics",} \code{"bayseq")} for an example
#' and how you can modify it.
#' @param libsize.list an optional named list where names represent samples (MUST
#' be the same as the samples in \code{sample.list}) and members are the library
#' sizes (the sequencing depth) for each sample. If not provided, they will be
#' estimated from baySeq.
#' @return A named list of the value 1-likelihood that a gene is differentially
#' expressed, whose names are the names of the contrasts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' contrast <- "A_vs_B"
#' norm.data.matrix <- normalize.edaseq(data.matrix,sample.list,gene.data)
#' p <- stat.bayseq(norm.data.matrix,sample.list,contrast)
#'}
stat.bayseq <- function(object,sample.list,contrast.list=NULL,stat.args=NULL,
    libsize.list=NULL) {
    if (is.null(stat.args))
        stat.args <- get.defaults("statistics","bayseq")
    if (is.null(contrast.list))
        contrast.list <- make.contrast.list(paste(names(sample.list)[1:2],
            collapse="_vs_"),sample.list)
    if (!is.list(contrast.list))
        contrast.list <- make.contrast.list(contrast.list,sample.list)
    classes <- as.class.vector(sample.list)
    p <- vector("list",length(contrast.list))
    names(p) <- names(contrast.list)
    switch(class(object),
        CountDataSet = { # Has been normalized with DESeq
            bayes.data <- counts(object,normalized=TRUE)
        },
        DGEList = { # Has been normalized with edgeR
            scl <- object$samples$lib.size * object$samples$norm.factors
            bayes.data <- round(t(t(object$counts)/scl)*mean(scl))
        },
        matrix = { # Has been normalized with EDASeq or NOISeq
            bayes.data <- object
        },
        list = {
            bayes.data <- as.matrix(round(sweep(object$counts,2,
                object$norm.factors,"*")))
        },
        nbp = {
            bayes.data <- as.matrix(round(object$pseudo.counts))
        }
    )
    CD <- new("countData",data=bayes.data,replicates=classes)
    if (is.null(libsize.list))
        baySeq::libsizes(CD) <- baySeq::getLibsizes(CD)
    else
        baySeq::libsizes(CD) <- unlist(libsize.list)
    for (con.name in names(contrast.list)) {
        disp("  Contrast: ", con.name)
        con <- contrast.list[[con.name]]
        #cd <- CD[,names(unlist(con))]
        cd <- CD[,match(names(unlist(con)),colnames(CD@data))]
        if (length(con)==2)
            baySeq::groups(cd) <- list(NDE=rep(1,length(unlist(con))),
                DE=c(rep(1,length(con[[1]])),rep(2,length(con[[2]]))))
        else
            baySeq::groups(cd) <- list(NDE=rep(1,length(unlist(con))),
                DE=unlist(con,use.names=FALSE)) # Maybe this will not work
        baySeq::replicates(cd) <- as.factor(classes[names(unlist(con))])
        cd <- baySeq::getPriors.NB(cd,samplesize=stat.args$samplesize,
            samplingSubset=stat.args$samplingSubset,
            equalDispersions=stat.args$equalDispersions,
            estimation=stat.args$estimation,zeroML=stat.args$zeroML,
            consensus=stat.args$consensus,cl=stat.args$cl)
        cd <- baySeq::getLikelihoods(cd,pET=stat.args$pET,
            marginalise=stat.args$marginalise,subset=stat.args$subset,
            priorSubset=stat.args$priorSubset,bootStraps=stat.args$bootStraps,
            conv=stat.args$conv,nullData=stat.args$nullData,
            returnAll=stat.args$returnAll,returnPD=stat.args$returnPD,
            discardSampling=stat.args$discardSampling,cl=stat.args$cl)
        tmp <- baySeq::topCounts(cd,group="DE",number=nrow(cd))
        p[[con.name]] <- 1 - as.numeric(tmp[,"Likelihood"])
        names(p[[con.name]]) <- rownames(tmp)
        p[[con.name]] <- p[[con.name]][rownames(CD@data)]
        p[[con.name]][which(is.na(p[[con.name]]))] <- 1
    }
    return(p)
}

#' Statistical testing with NBPSeq
#'
#' This function is a wrapper over NBPSeq statistical testing. It accepts a matrix
#' of normalized gene counts or an S4 object specific to each normalization
#' algorithm supported by metaseqR.
#'
#' @param object a matrix or an object specific to each normalization algorithm
#' supported by metaseqr, containing normalized counts. Apart from matrix (also
#' for NOISeq), the object can be a SeqExpressionSet (EDASeq), CountDataSet
#' (DESeq), DGEList (edgeR) or list (NBPSeq).
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param contrast.list a named structured list of contrasts as returned by
#' \code{\link{make.contrast.list}} or just the vector of contrasts as defined in
#' the main help page of \code{\link{metaseqr}}.
#' @param stat.args a list of NBPSeq statistical algorithm parameters. See the
#' result of \code{get.defaults("statistics",} \code{"nbpseq")}
#' for an example and how you can modify it. It is not required when the input
#' object is already a list from NBPSeq normalization as the dispersions are
#' already estimated.
#' @param libsize.list an optional named list where names represent samples
#' (MUST be the same as the samples \code{in sample.list}) and members are the
#' library sizes (the sequencing depth) for each sample. If not provided, the
#' default is the column sums of the \code{gene.counts} matrix.
#' @return A named list of p-values, whose names are the names of the contrasts.
#' @note There is currently a problem with the NBPSeq package and the workflow that
#' is specific to the NBPSeq package. The problem has to do with function exporting
#' as there are certain functions which are not recognized from the package
#' internally. For this reason and until it is fixed, only the Smyth workflow
#' will be available with the NBPSeq package.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' contrast <- "A_vs_B"
#' norm.data.matrix <- normalize.nbpseq(data.matrix,sample.list)
#' p <- stat.nbpseq(norm.data.matrix,sample.list,contrast)
#'}
stat.nbpseq <- function(object,sample.list,contrast.list=NULL,stat.args=NULL,
    libsize.list=NULL) {
    if (is.null(stat.args) && class(object)!="list")
        stat.args <- get.defaults("statistics","nbpseq")
    if (is.null(contrast.list))
        contrast.list <- make.contrast.list(paste(names(sample.list)[1:2],
            collapse="_vs_"),sample.list)
    if (!is.list(contrast.list))
        contrast.list <- make.contrast.list(contrast.list,sample.list)
    classes <- as.class.vector(sample.list)
    p <- vector("list",length(contrast.list))
    names(p) <- names(contrast.list)
    switch(class(object),
        CountDataSet = { # Has been normalized with DESeq
            counts <- round(counts(object,normalized=TRUE))
            if (is.null(libsize.list)) {
                libsize.list <- vector("list",length(classes))
                names(libsize.list) <- unlist(sample.list,use.names=FALSE)
                for (n in names(libsize.list))
                    libsize.list[[n]] <- sum(counts[,n])
            }
            lib.sizes <- unlist(libsize.list)
        },
        DGEList = { # Has been normalized with edgeR
            # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
            scl <- object$samples$lib.size * object$samples$norm.factors
            counts <- round(t(t(object$counts)/scl)*mean(scl))            
            if (is.null(libsize.list)) {
                libsize.list <- vector("list",length(classes))
                names(libsize.list) <- unlist(sample.list,use.names=FALSE)
                for (n in names(libsize.list))
                    libsize.list[[n]] <- sum(counts[,n])
            }
            lib.sizes <- unlist(libsize.list)
        },
        matrix = { # Has been normalized with EDASeq or NOISeq
            counts <- object
            if (is.null(libsize.list)) {
                libsize.list <- vector("list",length(classes))
                names(libsize.list) <- unlist(sample.list,use.names=FALSE)
                for (n in names(libsize.list))
                    libsize.list[[n]] <- sum(counts[,n])
            }
            lib.sizes <- unlist(libsize.list)
        },
        list = { # Has been normalized with NBPSeq
            object$counts <- as.matrix(object$counts)
            nb.data <- object
            if (is.null(libsize.list)) {
                libsize.list <- vector("list",length(classes))
                names(libsize.list) <- unlist(sample.list,use.names=FALSE)
                for (n in names(libsize.list))
                    libsize.list[[n]] <- sum(nb.data$counts[,n])
            }
            lib.sizes <- unlist(libsize.list)
            nb.data$pseudo.lib.sizes=rep(1e+7,dim(object$counts)[2])
        },
        nbp = { # Same...
            object$pseudo.counts <- as.matrix(object$pseudo.counts)
            nb.data <- object
            if (is.null(libsize.list)) {
                libsize.list <- vector("list",length(classes))
                names(libsize.list) <- unlist(sample.list,use.names=FALSE)
                for (n in names(libsize.list))
                    libsize.list[[n]] <- sum(nb.data$counts[,n])
            }
            lib.sizes <- unlist(libsize.list)
        }
    )
    # To avoid repeating the following chunk in the above
    if (class(object)!="list" && class(object)!="nbp") {
        #if (stat.args$main.method=="nbpseq") {
        #    nb.data <- list(
        #        counts=as.matrix(counts),
        #        lib.sizes=lib.sizes,
        #        norm.factors=rep(1,dim(counts)[2]),
        #        eff.lib.sizes=lib.sizes*rep(1,dim(counts)[2]),
        #        rel.frequencies=as.matrix(sweep(counts,2,
        #            lib.sizes*rep(1,dim(counts)[2]),"/")),
        #        tags=matrix(row.names(counts),dim(counts)[1],1)
        #    )
        #}
        #else if (stat.args$main.method=="nbsmyth") {
        #    nb.data <- new("nbp",list(
        #        counts=as.matrix(counts),
        #        lib.sizes=lib.sizes,
        #        grp.ids=classes,
        #        eff.lib.sizes=lib.sizes*rep(1,dim(counts)[2]),
        #        pseudo.counts=as.matrix(counts),
        #        pseudo.lib.sizes=colSums(as.matrix(counts))*rep(1,dim(counts)[2])
        #    ))
            nb.data <- list(
                counts=as.matrix(counts),
                lib.sizes=lib.sizes,
                grp.ids=classes,
                eff.lib.sizes=lib.sizes*rep(1,dim(counts)[2]),
                pseudo.counts=as.matrix(counts),
                #pseudo.lib.sizes=colSums(as.matrix(counts)) *
                #    rep(1,dim(counts)[2])
                pseudo.lib.sizes=rep(1e+7,dim(counts)[2])
            )
            class(nb.data) <- "nbp"
        #}
    }
    for (con.name in names(contrast.list)) {
        disp("  Contrast: ", con.name)
        con <- contrast.list[[con.name]]
        cons <- unique(unlist(con))
        if (length(con)==2) {
            #if (stat.args$main.method=="nbpseq") {
            #    dispersions <- estimate.dispersion(nb.data,
            #    model.matrix(~classes),
            #        method=stat.args$method$nbpseq)
            #    res <- test.coefficient(nb.data,dispersion=dispersions,
            #        x=model.matrix(~classes),beta0=c(NA,0),
            #        tests=stat.args$tests,
            #        alternative=stat.args$alternative,print.level=1)
            #    #res <- nb.glm.test(nb.data$counts,x=model.matrix(~classes),
            #        beta0=c(NA,0),lib.sizes=lib.sizes,
            #    #    dispersion.method=stat.args$method$nbpseq,
            #         tests=stat.args$tests)
            #    p[[con.name]] <- res[[stat.args$tests]]$p.values
            #    #p[[con.name]] <- res$test[[stat.args$tests]]$p.values
            #}
            #else if (stat.args$main.method=="nbsmyth") {
                obj <- suppressWarnings(estimate.disp(nb.data,
                    model=stat.args$model$nbsmyth,print.level=0))
                obj <- exact.nb.test(obj,cons[1],cons[2],print.level=0)
                p[[con.name]] <- obj$p.values
            #}
        }
        else {
            warnwrap(paste("NBPSeq differential expression algorithm does not ",
                "support ANOVA-like designs with more than two conditions to ",
                "be compared! Switching to DESeq for this comparison:",
                con.name))
            cds <- newCountDataSet(nb.data$counts,
                data.frame(condition=unlist(con),
                row.names=names(unlist(con))))
            sizeFactors(cds) <- rep(1,ncol(cds))
            cds <- estimateDispersions(cds,method="blind",
                sharingMode="fit-only")
            fit0 <- fitNbinomGLMs(cds,count~1)
            fit1 <- fitNbinomGLMs(cds,count~condition)
            p[[con.name]] <- nbinomGLMTest(fit1,fit0)
        }
        names(p[[con.name]]) <- rownames(nb.data$counts)
        p[[con.name]][which(is.na(p[[con.name]]))] <- 1
    }
    return(p)
}
