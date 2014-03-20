#' Normalization based on the EDASeq package
#'
#' This function is a wrapper over EDASeq normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param gene.counts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param norm.args a list of EDASeq normalization parameters. See the result of
#' \code{get.defaults("normalization",} \code{"edaseq")} for an example and how
#' you can modify it.
#' @param gene.data an optional annotation data frame (such the ones produced by
#' \code{get.annotation}) which contains the GC content for each gene and from
#' which the gene lengths can be inferred by chromosome coordinates.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the EDASeq native S4
#' object (SeqExpressionSet). In the latter case it should be handled with suitable
#' EDASeq methods.
#' @return A matrix or a SeqExpressionSet with normalized counts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplot.boxplot(data.matrix,sample.list)
#'
#' lengths <- round(1000*runif(nrow(data.matrix)))
#' starts <- round(1000*runif(nrow(data.matrix)))
#' ends <- starts + lengths
#' gc <- runif(nrow(data.matrix))
#' gene.data <- data.frame(
#'   chromosome=c(rep("chr1",nrow(data.matrix)/2),rep("chr2",nrow(data.matrix)/2)),
#'   start=starts,end=ends,gene_id=rownames(data.matrix),gc_content=gc
#' )
#' norm.data.matrix <- normalize.edaseq(data.matrix,sample.list,gene.data=gene.data)
#' diagplot.boxplot(norm.data.matrix,sample.list)
#'}
normalize.edaseq <- function(gene.counts,sample.list,norm.args=NULL,
    gene.data=NULL,output=c("matrix","native")) {
    if (is.null(norm.args))
        norm.args <- get.defaults("normalization","edaseq")
    if (!is.matrix(gene.counts))
        gene.counts <- as.matrix(gene.counts)
    if (!is.null(gene.data) && is.null(attr(gene.data,"gene.length")))
        attr(gene.data,"gene.length") <- rep(1,nrow(gene.counts))
    output <- tolower(output[1])
    check.text.args("output",output,c("matrix","native"))
    classes <- as.class.vector(sample.list)
    if (is.null(gene.data)) {
        seq.genes <- newSeqExpressionSet(
            gene.counts,
            phenoData=AnnotatedDataFrame(
                data.frame(
                    conditions=classes,
                    row.names=colnames(gene.counts)
                )
            ),
            featureData=AnnotatedDataFrame(
                data.frame(
                    length=rep(1,nrow(gene.counts)),
                    row.names=rownames(gene.counts)
                )
            )
        )
        seq.genes <- betweenLaneNormalization(withinLaneNormalization(seq.genes,
            "length",which=norm.args$within.which),
            which=norm.args$between.which)
    }
    else {
        seq.genes <- newSeqExpressionSet(
            gene.counts,
            phenoData=AnnotatedDataFrame(
                data.frame(
                    conditions=classes,
                    row.names=colnames(gene.counts)
                )
            ),
            featureData=AnnotatedDataFrame(
                data.frame(
                    gc=gene.data$gc_content,
                    length=attr(gene.data,"gene.length"),
                    row.names=rownames(gene.data)
                )
            )
        )
        seq.genes <- betweenLaneNormalization(withinLaneNormalization(seq.genes,
            "gc",which=norm.args$within.which),which=norm.args$between.which)
    }
    if (output=="matrix")
        return(exprs(seq.genes)) # Class: matrix
    else if (output=="native")
        return(seq.genes) # Class: SeqExpressionSet
}

#' Normalization based on the DESeq package
#'
#' This function is a wrapper over DESeq normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param gene.counts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param norm.args a list of DESeq normalization parameters. See the result of
#' \code{get.defaults("normalization",} \code{"deseq")} for an example and how you
#' can modify it.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the DESeq native S4
#' object (CountDataSet). In the latter case it should be handled with suitable
#' DESeq methods.
#' @return A matrix or a CountDataSet with normalized counts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplot.boxplot(data.matrix,sample.list)
#'
#' norm.data.matrix <- normalize.deseq(data.matrix,sample.list)
#' diagplot.boxplot(norm.data.matrix,sample.list)
#'}
normalize.deseq <- function(gene.counts,sample.list,norm.args=NULL,
    output=c("matrix","native")) {
    if (is.null(norm.args))
        norm.args <- get.defaults("normalization","deseq")
    output <- tolower(output[1])
    check.text.args("output",output,c("matrix","native"))
    classes <- as.class.vector(sample.list)
    cds <- newCountDataSet(gene.counts,classes)
    cds <- estimateSizeFactors(cds,locfunc=norm.args$locfunc)
    if (output=="native")
        return(cds) # Class: CountDataSet
    else if (output=="matrix")
        return(round(counts(cds,normalized=TRUE))) # Class: matrix
}

#' Normalization based on the edgeR package
#'
#' This function is a wrapper over edgeR normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param gene.counts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param norm.args a list of edgeR normalization parameters. See the result of
#' \code{get.defaults("normalization",} \code{"edger")} for an example and how
#' you can modify it.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the edgeR native S4
#' object (DGEList). In the latter case it should be handled with suitable edgeR
#' methods.
#' @return A matrix or a DGEList with normalized counts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplot.boxplot(data.matrix,sample.list)
#'
#' norm.data.matrix <- normalize.edger(data.matrix,sample.list)
#' diagplot.boxplot(norm.data.matrix,sample.list)
#'}
normalize.edger <- function(gene.counts,sample.list,norm.args=NULL,
    output=c("matrix","native")) {
    if (is.null(norm.args))
        norm.args <- get.defaults("normalization","edger")
    output <- tolower(output[1])
    check.text.args("output",output,c("matrix","native"))
    classes <- as.class.vector(sample.list)
    dge <- DGEList(counts=gene.counts,group=classes)
    dge <- calcNormFactors(dge,method=norm.args$method,
        refColumn=norm.args$refColumn,logratioTrim=norm.args$logratioTrim,
        sumTrim=norm.args$sumTrim,doWeighting=norm.args$doWeighting,
        Acutoff=norm.args$Acutoff,p=norm.args$p)
    if (output=="native")
        return(dge) # Class: DGEList
    else if (output=="matrix") {
        scl <- dge$samples$lib.size * dge$samples$norm.factors
        return(round(t(t(dge$counts)/scl)*mean(scl)))
        #return(round(dge$pseudo.counts)) # Class: matrix
    }
}

#' Normalization based on the NOISeq package
#'
#' This function is a wrapper over NOISeq normalization. It accepts a matrix of
#' gene counts (e.g. produced by importing an externally generated table of counts
#' to the main metaseqr pipeline).
#'
#' @param gene.counts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param norm.args a list of NOISeq normalization parameters. See the result of
#' \code{get.defaults("normalization",} \code{"noiseq")} for an example and how
#' you can modify it.
#' @param gene.data an optional annotation data frame (such the ones produced by
#' \code{get.annotation} which contains the GC content for each gene and from which
#' the gene lengths can be inferred by chromosome coordinates.
#' @param log.offset an offset to use to avoid infinity in logarithmic data
#' transformations.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the NOISeq native S4
#' object (SeqExpressionSet). In the latter case it should be handled with suitable
#' NOISeq methods.
#' @return A matrix with normalized counts.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplot.boxplot(data.matrix,sample.list)
#'
#' lengths <- round(1000*runif(nrow(data.matrix)))
#' starts <- round(1000*runif(nrow(data.matrix)))
#' ends <- starts + lengths
#' gc=runif(nrow(data.matrix)),
#' gene.data <- data.frame(
#'   chromosome=c(rep("chr1",nrow(data.matrix)/2),rep("chr2",nrow(data.matrix)/2)),
#'   start=starts,end=ends,gene_id=rownames(data.matrix),gc_content=gc
#' )
#' norm.data.matrix <- normalize.noiseq(data.matrix,sample.list,gene.data)
#' diagplot.boxplot(norm.data.matrix,sample.list)
#'}
normalize.noiseq <- function(gene.counts,sample.list,norm.args=NULL,
    gene.data=NULL,log.offset=1,output=c("matrix","native")) {
    if (is.null(norm.args))
        norm.args <- get.defaults("normalization","noiseq")
    output <- tolower(output[1])
    check.text.args("output",output,c("matrix","native"))
    classes <- as.class.vector(sample.list)
    if (is.null(gene.data)) {
        ns.obj <- NOISeq::readData(
            data=gene.counts,
            factors=data.frame(class=classes)
        )
    }
    else {
        gc.content <- gene.data$gc_content
        gene.length <- attr(gene.data,"gene.length")
        biotype <- as.character(gene.data$biotype)
        names(gc.content) <- names(biotype) <- names(gene.length) <- 
            rownames(gene.data)
        ns.obj <- NOISeq::readData(
            data=gene.counts,
            length=gene.length,
            gc=gc.content,
            chromosome=gene.data[,1:3],
            factors=data.frame(class=classes),
            biotype=biotype
        )
        #norm.args$long=gene.length # Set the gene length feature
    }
    norm.args$k=log.offset # Set the zero fixing constant
    switch(norm.args$method,
        rpkm = {
            #M <- NOISeq::rpkm(assayData(ns.obj)$exprs,long=norm.args$long,
            #    k=norm.args$k,lc=norm.args$lc)
            M <- rpkm(assayData(ns.obj)$exprs,gene.length=norm.args$long)
        },
        uqua = {
            M <- NOISeq::uqua(assayData(ns.obj)$exprs,long=norm.args$long,
                k=norm.args$k,lc=norm.args$lc)
        },
        tmm = {
            M <- NOISeq::tmm(assayData(ns.obj)$exprs,long=norm.args$long,
                k=norm.args$k,lc=norm.args$lc,refColumn=norm.args$refColumn,
                logratioTrim=norm.args$logratioTrim,sumTrim=norm.args$sumTrim,
                doWeighting=norm.args$doWeighting,Acutoff=norm.args$Acutoff)
        }
    )
    if (output=="native") {
        if (is.null(gene.data))
            return(NOISeq::readData(
            data=M,
            factors=data.frame(class=classes)
        )) # Class: CD
        else    
            return(NOISeq::readData(
                data=M,
                length=gene.length,
                gc=gc.content,
                chromosome=gene.data[,1:3],
                factors=data.frame(class=classes),
                biotype=biotype
            )) # Class: CD
    }
    else if (output=="matrix")
        return(as.matrix(round(M))) # Class: matrix
}

#' Normalization based on the NBPSeq package
#'
#' This function is a wrapper over DESeq normalization. It accepts a matrix of gene
#' counts (e.g. produced by importing an externally generated table of counts to
#' the main metaseqr pipeline).
#'
#' @param gene.counts a table where each row represents a gene and each column a
#' sample. Each cell contains the read counts for each gene and sample. Such a
#' table can be produced outside metaseqr and is imported during the basic metaseqr
#' workflow.
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param norm.args a list of NBPSeq normalization parameters. See the result of
#' \code{get.defaults("normalization",} \code{"nbpseq")} for an example and how
#' you can modify it.
#' @param libsize.list an optional named list where names represent samples (MUST
#' be the same as the samples in \code{sample.list}) and members are the library
#' sizes (the sequencing depth) for each sample. If not provided, the default is
#' the column sums of the \code{gene.counts} matrix.
#' @param output the class of the output object. It can be \code{"matrix"} (default)
#' for versatility with other tools or \code{"native"} for the NBPSeq native S4
#' object (a specific list). In the latter case it should be handled with suitable
#' NBPSeq methods.
#' @return A matrix with normalized counts or a list with the normalized counts
#' and other NBPSeq specific parameters.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' require(DESeq)
#' data.matrix <- counts(makeExampleCountDataSet())
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' diagplot.boxplot(data.matrix,sample.list)
#'
#' norm.data.matrix <- normalize.nbpseq(data.matrix,sample.list)
#' diagplot.boxplot(norm.data.matrix,sample.list)
#'}
normalize.nbpseq <- function(gene.counts,sample.list,norm.args=NULL,
    libsize.list=NULL,output=c("matrix","native")) {
    if (is.null(norm.args))
        norm.args <- get.defaults("normalization","nbpseq")
    output <- tolower(output[1])
    check.text.args("output",output,c("matrix","native"))
    classes <- as.class.vector(sample.list)
    if (is.null(libsize.list)) {
        libsize.list <- vector("list",length(classes))
        names(libsize.list) <- unlist(sample.list,use.names=FALSE)
        for (n in names(libsize.list))
            libsize.list[[n]] <- sum(gene.counts[,n])
    }
    lib.sizes <- unlist(libsize.list)
    norm.factors <- estimate.norm.factors(gene.counts,lib.sizes=lib.sizes,
        method=norm.args$method)
    #if (norm.args$main.method=="nbpseq")
    #    nb.data <- prepare.nb.data(gene.counts,lib.sizes=lib.sizes,
    #    norm.factors=norm.factors)
    #else if (norm.args$main.method=="nbsmyth")
        nb.data <- prepare.nbp(gene.counts,classes,lib.sizes=lib.sizes,
            norm.factors=norm.factors,thinning=norm.args$thinning)
    if (output=="native")
        return(nb.data) # Class: list or nbp
    else if (output=="matrix") {
        #if (norm.args$main.method=="nbpseq") {
        #    norm.counts <- matrix(0,nrow(gene.counts),ncol(gene.counts))
        #    for (i in 1:ncol(gene.counts))
        #        norm.counts[,i] <- norm.factors[i]*gene.counts[,i]
        #}
        #else if (norm.args$main.method=="nbsmyth") 
            norm.counts <- nb.data$pseudo.counts
        return(as.matrix(round(norm.counts))) # Class: matrix
    }
}
