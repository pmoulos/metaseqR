#' Filter gene expression based on exon counts
#'
#' This function performs the gene expression filtering based on exon read counts 
# and a set of exon filter rules. For more details see the main help pages of 
#' \code{\link{metaseqr}}.
#'
#' @param the.counts a named list created with the \code{\link{construct.gene.model}} 
#' function. See its help page for details.
#' @param gene.data an annotation data frame usually obtained with 
#' \code{\link{get.annotation}} containing the unique gene accession identifiers.
#' @param sample.list the list containing condition names and the samples under 
#' each condition.
#' @param exon.filters a named list with exon filters and their parameters. See 
#' the main help page of \code{\link{metaseqr}} for details.
#' @param restrict.cores in case of parallel execution of several subfunctions, 
#' the fraction of the available cores to use. In some cases if all available 
#' cores are used (\code{restrict.cores=1} and the system does not have sufficient 
#' RAM, the running machine might significantly slow down.
#' @return a named list with two members. The first member (\code{result} is a 
#' named list whose names are the exon filter names and its members are the filtered 
#' rownames of \code{gene.data}. The second member is a matrix of binary flags 
#' (0 for non-filtered, 1 for filtered) for each gene. The rownames of the flag 
#' matrix correspond to gene ids.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' data("hg19.exon.data",package="metaseqR")
#' exon.counts <- hg19.exon.counts
#' gene.data <- get.annotation("hg19","gene")
#' sample.list <- sample.list.hg19
#' exon.filters <- get.defaults("exon.filter")
#' the.counts <- construct.gene.model(exon.counts,sample.list,gene.data)
#' filter.results <- filter.exons(the.counts,gene.data,sample.list,exon.filters)
#'}
filter.exons <- function(the.counts,gene.data,sample.list,exon.filters,
    restrict.cores=0.8) {
    multic <- check.parallel(restrict.cores)
    exon.filter.result <- vector("list",length(exon.filters))
    names(exon.filter.result) <- names(exon.filters)
    flags <- matrix(0,nrow(gene.data),1)
    rownames(flags) <- rownames(gene.data)
    colnames(flags) <- c("MAE")
    the.genes <- as.character(gene.data$gene_id)
    if (!is.null(exon.filters)) {
        for (xf in names(exon.filters)) {
            disp("Applying exon filter ",xf,"...")
            switch(xf,
                min.active.exons = {
                    pass <- vector("list",length(unlist(sample.list)))
                    names(pass) <- names(the.counts)
                    for (n in names(pass)) {
                        disp("  Checking read presence in exons for ",n,"...")
                        pass[[n]] <- the.genes
                        names(pass[[n]]) <- the.genes
                        pass[[n]] <- wapply(multic,the.counts[[n]],
                            function(x,f) {
                            if (length(x$count) == 1)
                                if (x$count[1]!=0)
                                    return(FALSE)
                                else
                                    return(TRUE)
                            else if (length(x$count) > 1 &&
                                length(x$count) <= f$exons.per.gene)
                                if (length(which(x$count!=0)) >= f$min.exons)
                                    return(FALSE)
                                else
                                    return(TRUE)
                            else
                                if (length(which(x$count!=0)) >=
                                    ceiling(length(x$count)*f$frac))
                                    return(FALSE)
                                else
                                    return(TRUE)
                            },exon.filters$min.active.exons)
                        pass[[n]] <- do.call("c",pass[[n]])
                    }
                    pass.matrix <- do.call("cbind",pass)
                    exon.filter.result[[xf]] <- the.genes[which(apply(
                        pass.matrix,1,function(x) return(all(x))))]
                    flags[exon.filter.result[[xf]],"MAE"] <- 1
                }
                # More to come...
                # TODO: Write more rules based in exons
            )
        }
    }
    return(list(result=exon.filter.result,flags=flags))
}

#' Filter gene expression based on gene counts
#'
#' This function performs the gene expression filtering based on gene read counts 
#' and a set of gene filter rules. For more details see the main help pages of 
#' \code{\link{metaseqr}}.
#'
#' @param gene.counts a matrix of gene counts, preferably after the normalization 
#' procedure.
#' @param gene.data an annotation data frame usually obtained with 
#' \code{\link{get.annotation}} containing the unique gene accession identifiers.
#' @param gene.filters a named list with gene filters and their parameters. See 
#' the main help page of \code{\link{metaseqr}} for details.
#' @return a named list with three members. The first member (\code{result} is a 
#' named list whose names are the gene filter names and its members are the 
#' filtered rownames of \code{gene.data}. The second member (\code{cutoff} is a 
#' named list whose names are the gene filter names and its members are the cutoff 
#' values corresponding to each filter. The third member is a matrix of binary
#' flags (0 for non-filtered, 1 for filtered) for each gene. The rownames of the 
#' flag matrix correspond to gene ids.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' data("mm9.gene.data",package="metaseqR")
#' gene.counts <- mm9.gene.counts
#' sample.list <- mm9.sample.list
#' gene.counts <- normalize.edger(as.matrix(gene.counts[,9:12]),sample.list)
#' gene.data <- get.annotation("mm9","gene")
#' gene.filters <- get.defaults("gene.filter","mm9")
#' filter.results <- filter.genes(gene.counts,gene.data,gene.filters)
#'}
filter.genes <- function(gene.counts,gene.data,gene.filters,sample.list) {
    gene.filter.result <- gene.filter.cutoff <- vector("list",
        length(gene.filters))
    names(gene.filter.result) <- names(gene.filter.cutoff) <-
        names(gene.filters)
    flags <- matrix(0,nrow(gene.counts),9)
    rownames(flags) <- rownames(gene.counts)
    colnames(flags) <- c("LN","AR","MD","MN","QN","KN","CM","BT","PR")
    for (gf in names(gene.filters)) {
        disp("Applying gene filter ",gf,"...")
        switch(gf,
            length = { # This is real gene length independently of exons
                if (!is.null(gene.filters$length)) {
                    gene.filter.result$length <- rownames(gene.data)[which(
                        gene.data$end - gene.data$start <
                            gene.filters$length$length
                    )]
                    gene.filter.cutoff$length <- gene.filters$length$length
                    flags[intersect(gene.filter.result$length,
                        rownames(gene.counts)),"LN"] <- 1
                    disp("  Threshold below which ignored: ",
                        gene.filters$length$length)
                }
                else
                    gene.filter.cutoff$length <- NULL
            },
            avg.reads = {
                if (!is.null(gene.filters$avg.reads)) {
                    len <- attr(gene.data,"gene.length")
                    if (!is.null(len))
                        len <- len[rownames(gene.counts)]
                    else {
                        gg <- gene.data[rownames(gene.counts),]
                        len <- gg$end - gg$start
                    }
                    avg.mat <- sweep(gene.counts,1,
                        len/gene.filters$avg.reads$average.per.bp,"/")
                    q.t <- max(apply(avg.mat,2,quantile,
                        gene.filters$avg.reads$quantile))
                    gene.filter.result$avg.reads <- rownames(gene.data)[which(
                        apply(avg.mat,1,filter.low,q.t))]
                    gene.filter.cutoff$avg.reads <- q.t
                    flags[intersect(gene.filter.result$avg.reads,
                        rownames(gene.counts)),"AR"] <- 1
                    disp("  Threshold below which ignored: ",q.t)
                }
                else
                    gene.filter.cutoff$avg.reads <- NULL
            },
            expression = {
                if (!is.null(gene.filters$expression)) {
                    if (!is.null(gene.filters$expression$median) && 
                        gene.filters$expression$median) {
                        md <- median(gene.counts)
                        the.dead.median <- rownames(gene.counts)[which(
                            apply(gene.counts,1,filter.low,md))]
                        disp("  Threshold below which ignored: ",md)
                    }
                    else
                        the.dead.median <- md <- NULL
                    if (!is.null(gene.filters$expression$mean) &&
                        gene.filters$expression$mean) {
                        mn <- mean(gene.counts)
                        the.dead.mean <- rownames(gene.counts)[which(apply(
                            gene.counts,1,filter.low,mn))]
                        disp("  Threshold below which ignored: ",mn)
                    }
                    else
                        the.dead.mean <- mn <- NULL
                    if (!is.null(gene.filters$expression$quantile) &&
                        !is.na(gene.filters$expression$quantile)) {
                        qu <- quantile(gene.counts,
                            gene.filters$expression$quantile)
                        the.dead.quantile <- rownames(gene.counts)[which(
                            apply(gene.counts,1,filter.low,qu))]
                        disp("  Threshold below which ignored: ",qu)
                    }
                    else
                        the.dead.quantile <- qu <- NULL
                    if (!is.null(gene.filters$expression$known) &&
                        !is.na(gene.filters$expression$known)) {
                        # Think about the case of embedded
                        bio.cut <- match(gene.filters$expression$known,
                            gene.data$gene_name) 
                        bio.cut <- bio.cut[-which(is.na(bio.cut))]
                        bio.cut.counts <- as.vector(gene.counts[bio.cut,])
                        the.bio.cut <- quantile(bio.cut.counts,0.9)
                        the.dead.known <- rownames(gene.counts)[which(apply(
                            gene.counts,1,filter.low,the.bio.cut))]
                        disp("  Threshold below which ignored: ",the.bio.cut)
                    }
                    else
                        the.dead.known <- the.bio.cut <- NULL
                    if (!is.null(gene.filters$expression$custom) &&
                        !is.na(gene.filters$expression$custom)) {
                        # For future use
                        the.dead.custom <- NULL
                    }
                    else
                        the.dead.custom <- NULL
                    # Derive one common expression filter
                    gene.filter.result$expression$median <- the.dead.median
                    gene.filter.result$expression$mean <- the.dead.mean
                    gene.filter.result$expression$quantile <-
                        the.dead.quantile
                    gene.filter.result$expression$known <- the.dead.known
                    gene.filter.result$expression$custom <- the.dead.custom
                    gene.filter.cutoff$expression$median <- md
                    gene.filter.cutoff$expression$mean <- mn
                    gene.filter.cutoff$expression$quantile <- qu
                    gene.filter.cutoff$expression$known <- the.bio.cut
                    gene.filter.cutoff$expression$custom <- NULL
                    if (!is.null(the.dead.median)) flags[the.dead.median,
                        "MD"] <- 1
                    if (!is.null(the.dead.mean)) flags[the.dead.mean,"MN"] <- 1
                    if (!is.null(the.dead.quantile)) flags[the.dead.quantile,
                        "QN"] <- 1
                    if (!is.null(the.dead.known)) flags[the.dead.known,
                        "KN"] <- 1
                    if (!is.null(the.dead.custom)) flags[the.dead.custom,
                        "CM"] <- 1
                    #the.dead <- list(the.dead.median,the.dead.mean,
                    #    the.dead.quantile,the.dead.known,the.dead.custom)
                    #gene.filter.result$expression <- Reduce("union",the.dead)
                }
            },
            biotype = {
                if (!is.null(gene.filters$biotype)) {
                    filter.out <- names(which(unlist(gene.filters$biotype)))
                    # Necessary hack because of R naming system
                    if (length(grep("three_prime_overlapping_ncrna",
                        filter.out))>0) 
                        filter.out <- sub("three_prime_overlapping_ncrna",
                            "3prime_overlapping_ncrna",filter.out)
                    filter.ind <- vector("list",length(filter.out))
                    names(filter.ind) <- filter.out
                    for (bt in filter.out)
                        filter.ind[[bt]] <- rownames(gene.data)[which(
                            gene.data$biotype==bt)]
                    gene.filter.result$biotype <- Reduce("union",filter.ind)
                    gene.filter.cutoff$biotype <- paste(filter.out,
                        collapse=", ")
                    disp("  Biotypes ignored: ",paste(filter.out,
                        collapse=", "))
                }
                else
                    gene.filter.result$biotype <- NULL
                if (!is.null(gene.filter.result$biotype) && 
                    length(gene.filter.result$biotype)>0) 
                    flags[gene.filter.result$biotype,"BT"] <- 1
            },
            presence = {
				if (!is.null(gene.filters$presence)) {
                    frac <- gene.filters$presence$frac
                    minc <- gene.filters$presence$min.count
                    pc <- gene.filters$presence$per.condition
                    
                    if (pc) {
						np.tmp <- nam.tmp <- vector("list",length(sample.list))
						names(np.tmp) <- names(nam.tmp) <- names(sample.list)
						for (n in names(sample.list)) {
							tmp <- gene.counts[,sample.list[[n]]]
							nreq <- ceiling(frac*ncol(tmp))
							np.tmp[[n]] <- apply(tmp,1,function(x,m,n) {
								w <- which(x>=m)
								if (length(w)>0 && length(w)>n)
									return(FALSE)
								return(TRUE)
							},minc,nreq)
							nam.tmp[[n]] <- rownames(tmp[which(np.tmp[[n]]),,
								drop=FALSE])
						}
						the.dead.presence <- Reduce("union",nam.tmp)
					}
					else {
						nreq <- ceiling(frac*ncol(gene.counts))
						not.present <- apply(gene.counts,1,function(x,m,n) {
							w <- which(x>=m)
							if (length(w)>0 && length(w)>n)
								return(FALSE)
							return(TRUE)
						},minc,nreq)
						the.dead.presence <- 
							rownames(gene.counts[which(not.present),,
								drop=FALSE])
					}
                    
                    if (length(the.dead.presence)>0) {
						gene.filter.result$presence <- the.dead.presence
						flags[the.dead.presence,"PR"] <- 1
						gene.filter.cutoff$presence <- nreq
						disp("  Threshold below which ignored: ",nreq)
					}
					else
						gene.filter.result$presence <- 
							gene.filter.cutoff$presence <- NULL
                }
                else
                    gene.filter.result$presence <- NULL
			}
        )
    }
    return(list(result=gene.filter.result,cutoff=gene.filter.cutoff,
        flags=flags))
}
