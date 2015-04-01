#' Results export builder
#'
#' This function help build the output files of the metaseqr pipeline based on 
#' several elements produced during the pipeline execution. It is intended for 
#' internal use and not available to the users.
#'
#' @param gene.data an annotation data frame (such the ones produced by 
#' \code{\link{get.annotation}}).
#' @param raw.gene.counts a matrix of un-normalized gene counts.
#' @param norm.gene.counts a matrix of normalized gene counts.
#' @param flags a matrix of filtering flags (0,1), created by the filtering 
#' functions.
#' @param sample.list see the documentation of \code{\link{metaseqr}}.
#' @param cnt the statistical contrast for which the export builder is currently 
#' running.
#' @param statistics the statistical tests used (see the documentation of 
#' \code{\link{metaseqr}}).
#' @param raw.list a list of transformed un-normalized counts, see the documentation 
#' of \code{\link{make.transformation}}.
#' @param norm.list a list of transformed normalized counts, see the documentation 
#' of \code{\link{make.transformation}}.
#' @param p.mat a matrix of p-values, see the documentation of \code{\link{metaseqr}}.
#' @param adj.p.mat a matrix of adjusted p-values, see the documentation of 
#' \code{\link{metaseqr}}.
#' @param sum.p a vector of combined p-values, see the documentation of 
#' \code{\link{metaseqr}}.
#' @param adj.sum.p a vector of adjusted combined p-values, see the documentation 
#' of \code{\link{metaseqr}}.
#' @param export.what see the documentation of \code{\link{metaseqr}}.
#' @param export.scale see the documentation of \code{\link{metaseqr}}.
#' @param export.values see the documentation of \code{\link{metaseqr}}.
#' @param export.stats see the documentation of \code{\link{metaseqr}}.
#' @param log.offset see the documentation of \code{\link{metaseqr}}.
#' @param report see the documentation of \code{\link{metaseqr}}.
#' @return A list with three members: a data frame to be exported in a text file, 
#' a long string with the result in a html formatted table (if \code{report=TRUE}) 
#' and the column names of the output data frame.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' # Not yet available
#'}
build.export <- function(gene.data,raw.gene.counts,norm.gene.counts,flags,
    sample.list,cnt,statistics,
    raw.list,norm.list,
    p.mat=matrix(NA,nrow(gene.data),length(statistics)),
    adj.p.mat=matrix(NA,nrow(gene.data),length(statistics)),
    sum.p=rep(NA,nrow(gene.data)),
    adj.sum.p=rep(NA,nrow(gene.data)),
    export.what=c("annotation","p.value","adj.p.value","meta.p.value",
        "adj.meta.p.value","fold.change","stats","counts","flags"),
    export.scale=c("natural","log2","log10","rpgm","vst"),
    export.values=c("raw","normalized"),
    export.stats=c("mean","median","sd","mad","cv","rcv"),
    log.offset=1,report=TRUE
) {    
    if (is.null(colnames(p.mat)))
        colnames(p.mat) <- statistics
    if (is.null(adj.p.mat))
        adj.p.mat=matrix(NA,nrow(gene.data),length(statistics))
    if (is.null(colnames(adj.p.mat)))
        colnames(adj.p.mat) <- statistics

    export <- data.frame(row.names=rownames(gene.data))
    if (report) export.html <- as.matrix(export)
    the.names <- character(0)
    if ("annotation" %in% export.what) {
        disp("      binding annotation...")
        export <- cbind(export,gene.data)
        if (report)
            export.html <- cbind(export.html,make.html.cells(gene.data,
                type="text"))
        the.names <- c(the.names,colnames(gene.data))
    }
    if ("p.value" %in% export.what) {
        disp("      binding p-values...")
        export <- cbind(export,p.mat)
        if (report) 
            export.html <- cbind(export.html,make.html.cells(p.mat))
        the.names <- c(the.names,paste("p-value_",colnames(p.mat),sep=""))
    }
    if ("adj.p.value" %in% export.what) {
        disp("      binding FDRs...")
        export <- cbind(export,adj.p.mat)
        if (report) 
            export.html <- cbind(export.html,make.html.cells(adj.p.mat))
        the.names <- c(the.names,paste("FDR_",colnames(adj.p.mat),sep=""))
    }
    if ("meta.p.value" %in% export.what && length(statistics)>1) { 
        # Otherwise, it does not exist
        disp("      binding meta p-values...")
        export <- cbind(export,sum.p)
        if (report) 
            export.html <- cbind(export.html,make.html.cells(sum.p))
        the.names <- c(the.names,paste("meta_p-value_",cnt,sep=""))
    }
    if ("adj.meta.p.value" %in% export.what && length(statistics)>1) {
        disp("      binding adjusted meta p-values...")
        export <- cbind(export,adj.sum.p)
        if (report) 
            export.html <- cbind(export.html,make.html.cells(adj.sum.p))
        the.names <- c(the.names,paste("meta_FDR_",cnt,sep=""))
    }
    if ("fold.change" %in% export.what) {
        if ("normalized" %in% export.values) {
            tmp <- make.fold.change(cnt,sample.list,norm.gene.counts,log.offset)
            if ("natural" %in% export.scale) {
                disp("      binding natural normalized fold changes...")
                export <- cbind(export,tmp)
                if (report) 
                    export.html <- cbind(export.html,make.html.cells(tmp))
                the.names <- c(the.names,
                    paste("natural_normalized_fold_change_",colnames(tmp),
                        sep=""))
            }
            if ("log2" %in% export.scale) {
                disp("      binding log2 normalized fold changes...")
                export <- cbind(export,log2(tmp))
                if (report) 
                    export.html <- cbind(export.html,
                        make.html.cells(log2(tmp)))
                the.names <- c(the.names,paste("log2_normalized_fold_change_",
                    colnames(tmp),sep=""))
            }
        }
        if ("raw" %in% export.values) {
            tmp <- make.fold.change(cnt,sample.list,raw.gene.counts,log.offset)
            if ("natural" %in% export.scale) {
                disp("      binding natural raw fold changes...")
                export <- cbind(export,tmp)
                if (report) 
                    export.html <- cbind(export.html,make.html.cells(tmp))
                the.names <- c(the.names,paste("natural_raw_fold_change_",
                    colnames(tmp),sep=""))
            }
            if ("log2" %in% export.scale) {
                disp("      binding log2 raw fold changes...")
                export <- cbind(export,log2(tmp))
                if (report) 
                    export.html <- cbind(export.html,make.html.cells(log2(tmp)))
                the.names <- c(the.names,paste("log2_raw_fold_change_",
                    colnames(tmp),sep=""))
            }
        }
    }
    if ("stats" %in% export.what) {
        conds <- strsplit(cnt,"_vs_")[[1]]
        for (cond in conds) {
            if ("normalized" %in% export.values) {
                if ("mean" %in% export.stats) {
                    disp("      binding normalized mean counts...")
                    tmp <- make.stat(sample.list[[cond]],norm.list,"mean",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report) 
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_normalized_mean_counts_",cond,sep=""))
                }
                if ("median" %in% export.stats) {
                    disp("      binding normalized median counts...")
                    tmp <- make.stat(sample.list[[cond]],norm.list,"median",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report) 
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_normalized_median_counts_",cond,sep=""))
                }
                if ("sd" %in% export.stats) {
                    disp("      binding normalized count sds...")
                    tmp <- make.stat(sample.list[[cond]],norm.list,"sd",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report) 
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_normalized_sd_counts_",cond,sep=""))
                }
                if ("mad" %in% export.stats) {
                    disp("      binding normalized count MADs...")
                    tmp <- make.stat(sample.list[[cond]],norm.list,"mad",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report)
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_normalized_mad_counts_",cond,sep=""))
                }
                if ("cv" %in% export.stats) {
                    disp("      binding normalized count CVs...")
                    tmp <- make.stat(sample.list[[cond]],norm.list,"cv",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report)
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_normalized_cv_counts_",cond,sep=""))
                }
                if ("rcv" %in% export.stats) {
                    disp("      binding normalized counts RCVs...")
                    tmp <- make.stat(sample.list[[cond]],norm.list,"rcv",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report)
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_normalized_rcv_counts_",cond,sep=""))
                }
            }
            if ("raw" %in% export.values) {
                if ("mean" %in% export.stats) {
                    disp("      binding raw mean counts...")
                    tmp <- make.stat(sample.list[[cond]],raw.list,"mean",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report)
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_raw_mean_counts_",cond,sep=""))
                }
                if ("median" %in% export.stats) {
                    disp("      binding raw median counts...")
                    tmp <- make.stat(sample.list[[cond]],raw.list,"median",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report)
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_raw_median_counts_",cond,sep=""))
                }
                if ("sd" %in% export.stats) {
                    disp("      binding raw counts sds...")
                    tmp <- make.stat(sample.list[[cond]],raw.list,"sd",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report)
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_raw_sd_counts_",cond,sep=""))
                }
                if ("mad" %in% export.stats) {
                    disp("      binding raw counts MADs...")
                    tmp <- make.stat(sample.list[[cond]],raw.list,"mad",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report)
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_raw_mad_counts_",cond,sep=""))
                }
                if ("cv" %in% export.stats) {
                    disp("      binding raw counts CVs...")
                    tmp <- make.stat(sample.list[[cond]],raw.list,"cv",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report)
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_raw_cv_counts_",cond,sep=""))
                }
                if ("rcv" %in% export.stats) {
                    disp("      binding raw counts RCVs...")
                    tmp <- make.stat(sample.list[[cond]],raw.list,"rcv",
                        export.scale)
                    export <- cbind(export,tmp)
                    if (report)
                        export.html <- cbind(export.html,make.html.cells(tmp))
                    the.names <- c(the.names,paste(colnames(tmp),
                        "_raw_rcv_counts_",cond,sep=""))
                }
            }
        }
    }
    if ("counts" %in% export.what) {
        conds <- strsplit(cnt,"_vs_")[[1]]
        for (cond in conds) {
            if ("normalized" %in% export.values) {                    
                disp("      binding all normalized counts for ",cond,"...")
                tmp <- make.matrix(sample.list[[cond]],norm.list,export.scale)
                export <- cbind(export,tmp)
                if (report)
                    export.html <- cbind(export.html,make.html.cells(tmp))
                part.1 <- rep(paste(export.scale,"_normalized_counts_",sep=""),
                    each=length(sample.list[[cond]]))
                part.2 <- paste(part.1,colnames(tmp),sep="")
                the.names <- c(the.names,part.2)
            }
            if ("raw" %in% export.values) {
                disp("      binding all raw counts for ",cond,"...")
                tmp <- make.matrix(sample.list[[cond]],raw.list,export.scale)
                export <- cbind(export,tmp)
                if (report)
                    export.html <- cbind(export.html,make.html.cells(tmp))
                part.1 <- rep(paste(export.scale,"_raw_counts_",sep=""),
                    each=length(sample.list[[cond]]))
                part.2 <- paste(part.1,colnames(tmp),sep="")
                the.names <- c(the.names,part.2)
            }
        }
    }
    if ("flags" %in% export.what && !is.null(flags)) {
        disp("      binding filtering flags...")
        export <- cbind(export,as.data.frame(flags))
        if (report)
            export.html <- cbind(export.html,make.html.cells(flags))
        the.names <- c(the.names,colnames(flags))
    }
    names(export) <- the.names

    if (!report)
        export.html <- NULL

    return (list(text.table=export,html.table=export.html,headers=the.names))
}
