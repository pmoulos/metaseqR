#' Main argument validator
#'
#' Checks if the arguments passed to \code{\link{metaseqr}} are valid and throws
#' a warning about the invalid ones (which are ignored anyway because of the
#' \code{...} in \code{\link{metaseqr}}. However, for this reason this function
#' is useful as some important parameter faults might go unnoticed in the beginning
#' and cause a failure afterwards. Internal use.
#' 
#' @param main.args a list of parameters with which metaseqr is called (essentially,
#' the output of \code{\link{match.call}}.
#' @author Panagiotis Moulos
check.main.args <- function(main.args) {
    in.args <- names(main.args)[-1] # 1st member name of calling function
    valid.args <- c(
        "counts","sample.list","exclude.list","file.type","path","contrast",
        "libsize.list","id.col","gc.col","name.col","bt.col","annotation",
        "gene.file","org","trans.level","count.type","utr.flank","exon.filters",
        "gene.filters","when.apply.filter","normalization","norm.args",
        "statistics","stat.args","adjust.method","meta.p","weight","nperm",
        "reprod","pcut","log.offset","preset","qc.plots","fig.format",
        "out.list","export.where","export.what","export.scale","export.values",
        "export.stats","export.counts.table","restrict.cores","report","refdb",
        "report.top","report.template","verbose","run.log","save.gene.model"
    )
    invalid <- setdiff(in.args,valid.args)
    if (length(invalid) > 0) {
        for (i in 1:length(invalid))
            warnwrap("Unknown input argument to metaseqr pipeline: ",invalid[i],
                " ...Ignoring...",now=TRUE)
    }
}

#' Text argument validator
#'
#' Checks if one or more given textual argument(s) is/are member(s) of a list of 
#' correct arguments. It's a more package-specific function similar to 
#' \code{\link{match.arg}}. Mostly for internal use.
#' 
#' @param arg.name the name of the argument that is checked (for display purposes).
#' @param arg.value the value(s) of the argument to be checked.
#' @param arg.list a vector of valid argument values for \code{arg.value} to be 
#' matched against.
#' @param multiarg a logical scalar indicating whether \code{arg.name} accepts 
#' multiple arguments or not. In that case, all of the values in \code{arg.value} 
#' are checked against \code{arg.list}.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' check.text.args("count.type",count.type,c("gene","exon"),multiarg=FALSE)
#' check.text.args("statistics",statistics,c("deseq","edger","noiseq","bayseq",
#'  "limma"), multiarg=TRUE)
#'}
check.text.args <- function(arg.name,arg.value,arg.list,multiarg=FALSE) {
    if (multiarg) {
        arg.value <- tolower(arg.value)
        if (!all(arg.value %in% arg.list))
            stopwrap("\"",arg.name,"\""," parameter must be one or more of ",
                paste(paste("\"",arg.list,sep=""),collapse="\", "),"\"!")
    }
    else {
        arg.save <- arg.value[1]
        arg.value <- tolower(arg.value[1])
        # An exception must be added for annotation because it can be an external 
        # file too
        if (arg.name=="annotation") { 
            if (!(arg.value %in% arg.list) && !file.exists(arg.save))
                stopwrap("\"",arg.name,"\""," parameter must be one of ",
                    paste(paste("\"",arg.list,sep=""),collapse="\", "),
                        "\" or an existing file!")
        }
        else {
            if (!(arg.value %in% arg.list))
                stopwrap("\"",arg.name,"\""," parameter must be one of ",
                    paste(paste("\"",arg.list,sep=""),collapse="\", "),"\"!")
        }
    }
}

#' Numeric argument validator
#'
#' Checks if one or more given numeric argument(s) satisfy several rules concerning 
#' numeric arguments, e.g. proper bounds or proper format (e.g. it must be a number 
#' and not a character). Mostly for internal use.
#' 
#' @param arg.name the name of the argument that is checked (for display purposes).
#' @param arg.value the value(s) of the argument to be checked.
#' @param arg.type either the string \code{"numeric"} to denote generic double-like 
#' R numerics or \code{"integer"} for integer values.
#' @param arg.bounds a numeric or a vector with 2 elements, restraining 
#' \code{arg.value} to be within the bounds defined by the input vector or e.g. 
#' larger (smaller) than the numeric value. See examples.
#' @param direction a string denoting to which direction the \code{arg.value} 
#' should be compared with \code{arg.bounds}. For example, \code{"both"} should 
#' be given with a two element vector against which, \code{arg.value} will be 
#' checked to see whether it is smaller than the low boundary or larger than the 
#' higher boundary. In that case, the function will throw an error. The direction 
#' parameter can be one of: \code{"both"} (described above), \code{"botheq"} (as 
#' above, but the \code{arg.val} is also checked for equality -closed intervals), 
#' \code{"gt"} or \code{"gte"} (check whether \code{arg.val} is smaller or smaller 
#' than or equal to the first value of \code{arg.bounds}), \code{"lt"} or \code{"lte"} 
#' (check whether \code{arg.val} is larger or larger than or equal to the first 
#' value of \code{arg.bounds}).
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' pcut <- 1.2 # A probability cannot be larger than 1! It will throw an error!
#' check.num.args("pcut",pcut,"numeric",c(0,1),"botheq")
#' pcut <- 0.05 # Pass
#' check.num.args("pcut",pcut,"numeric",c(0,1),"botheq")
#' gc.col <- 3.4 # A column in a file cannot be real! It will throw an error!
#' check.num.args("gc.col",gc.col,"integer",0,"gt")
#' gc.col <- 5 # Pass
#' check.num.args("gc.col",gc.col,"integer",0,"gt")
#'}
check.num.args <- function(arg.name,arg.value,arg.type,arg.bounds,direction) {
    switch(arg.type,
        numeric = {
            if (!is.numeric(arg.value))
                stopwrap("\"",arg.name,"\"",
                    " parameter must be a numeric value!")
            if (!missing(arg.bounds)) {
                switch(direction,
                    both = {
                        if (arg.value<=arg.bounds[1] ||
                            arg.value>=arg.bounds[2])
                            stopwrap("\"",arg.name,"\""," parameter must be a ",
                                "numeric ","value larger than or equal to ",
                                arg.bounds[1]," and smaller than or equal to ",
                                arg.bounds[2],"!")
                    },
                    botheq = {
                        if (arg.value<arg.bounds[1] || arg.value>arg.bounds[2])
                            stopwrap("\"",arg.name,"\""," parameter must be a ",
                                "numeric value larger than ",arg.bounds[1],
                                " and smaller than ",arg.bounds[2],"!")
                    },
                    gt = {
                        if (arg.value<=arg.bounds[1])
                            stopwrap("\"",arg.name,"\""," parameter must be a ",
                                "numeric value greater than ",arg.bounds[1],"!")
                    },
                    lt = {
                        if (arg.value>=arg.bounds[1])
                            stopwrap("\"",arg.name,"\""," parameter must be a ",
                                "numeric value lower than ",arg.bounds[1],"!")
                    },
                    gte = {
                        if (arg.value<arg.bounds[1])
                            stopwrap("\"",arg.name,"\""," parameter must be a ",
                                "numeric value greater than or equal to ",
                                arg.bounds[1],"!")
                    },
                    lte = {
                        if (arg.value>arg.bounds[1])
                            stopwrap("\"",arg.name,"\""," parameter must be a ",
                                "numeric value lower than or equal to ",
                                arg.bounds[1],"!")
                    }
                )
            }
        },
        integer = {
            if (!is.integer(arg.value))
                stopwrap("\"",arg.name,"\""," parameter must be an integer!")
            if (!missing(arg.bounds)) {
                switch(direction,
                    both = {
                        if (arg.value<=arg.bounds[1] ||
                            arg.value>=arg.bounds[2])
                            stopwrap("\"",arg.name,"\""," parameter must be ",
                                "an integer larger than or equal to ",
                                arg.bounds[1]," and smaller than or equal to ",
                                arg.bounds[2],"!")
                    },
                    botheq = {
                        if (arg.value<arg.bounds[1] || arg.value>arg.bounds[2])
                            stopwrap("\"",arg.name,"\""," parameter must be ",
                                "an integer larger than or equal to ",
                                arg.bounds[1]," and smaller than or equal to ",
                                arg.bounds[2],"!")
                    },
                    gt = {
                        if (arg.value<=arg.bounds[1])
                            stopwrap("\"",arg.name,"\""," parameter must be ",
                                "an integer greater than ",arg.bounds[1],"!")
                    },
                    lt = {
                        if (arg.value>=arg.bounds[1])
                            stopwrap("\"",arg.name,"\""," parameter must be ",
                                "an integer lower than ",arg.bounds[1],"!")
                    },
                    gte = {
                        if (arg.value<arg.bounds[1])
                            stopwrap("\"",arg.name,"\""," parameter must be ",
                                "an integer greater than or equal to ",
                                arg.bounds[1],"!")
                    },
                    lte = {
                        if (arg.value>arg.bounds[1])
                            stopwrap("\"",arg.name,"\""," parameter must be ",
                                "an integer lower than or equal to ",
                                arg.bounds[1],"!")
                    }
                )
            }
        }
    )
}

#' File argument validator
#'
#' Checks if a file exists for specific arguments requiring a file input. Internal 
#' use only.
#'
#' @param arg.name argument name to display in a possible error.
#' @param arg.value the filename to check.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # OK
#' check.file.args("file",system.file("metaseqr_report.html",package="metaseqR"))
#' # Error!
#' check.file.args("file",system.file("metaseqr_report.htm",package="metaseqR"))
#'}
check.file.args <- function(arg.name,arg.value) {
    if (!file.exists(arg.value))
        stopwrap("\"",arg.name,"\""," parameter must be an existing file!")
}

#' Parallel run validator
#'
#' Checks existence of multiple cores and loads parallel package.
#'
#' @param rc fraction of available cores to use.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#  multic <- check.parallel(0.8)
#'}
check.parallel <- function(rc) {
    if (suppressWarnings(!require(parallel)) || .Platform$OS.type!="unix")
        multi <- FALSE
    else {
        multi <- TRUE
        ncores <- parallel::detectCores()
        if (!missing(rc) || !is.na(rc) || !is.null(rc))
            ncores <- ceiling(rc*ncores)
        options(cores=ncores)
    }
    return(multi)
}

#' Contrast validator
#'
#' Checks if the contrast vector follows the specified format. Internal use only.
#'
#' @param cnt contrasts vector.
#' @param sample.list the input sample list.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' dontrun{
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' cnt <- c("A_vs_B") # Will work
#' #cnt <- c("A_vs_C") ## Will throw error!
#' check.contrast.format(cnt,sample.list)
#}
check.contrast.format <- function(cnt,sample.list) {
    # This function will break cnt and check that all contrast counter parts are 
    # members of the names of the sample.list and contain the string "_vs_" as 
    # many times as the names of the sample.list minus 1. If satisfied return 
    # TRUE else error.
    cnts <- strsplit(cnt,"_vs_")
    #if (length(unique(unlist(cnts))) != length(names(sample.list)))
    if (!all(unique(unlist(cnts)) %in% names(sample.list)))
        stopwrap("Condition names in sample list and contrast list do not ",
            "match! Check if the contrasts follow the appropriate format (e.g.",
            " \"_vs_\" separating contrasting conditions...")
    if (length(unique(cnt))!=length(cnt))
        warnwrap("Duplicates found in the contrasts list! Duplicates will be ",
            "ignored...")
}

#' Library size validator
#'
#' Checks the names of the supplied library sizes. Internal use only.
#'
#' @param libsize.list the samples-names library size list.
#' @param sample.list the input sample list.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' sample.list <- list(A=c("A1","A2"),B=c("B1","B2","B3"))
#' libsize.list.1 <- list(A1=1e+6,A2=1.1e+6,B1=1.2e+6,B2=1.3e+6,B3=1.5e+6)
#' libsize.list.2 <- list(A1=1e+6,A2=1.1e+6,B1=1.2e+6,B2=1.3e+6)
#' check.libsize(libsize.list.1,sample.list) # Will work
#' #check.libsize(libsize.list.2,sample.list) # Will throw error!
#'}
check.libsize <- function(libsize.list,sample.list) {
    if (!is.null(libsize.list)) {
        if (length(intersect(names(libsize.list),unlist(sample.list,
            use.names=FALSE)))!=length(unlist(sample.list,
            use.names=FALSE))) {
            warnwrap("Sample names in \"libsize.list\" and \"sample.list\" do ",
                "not match! Library sizes will be estimated from count data...")
            return(NULL)
        }
        else return(libsize.list)
    }
    else
        return(NULL)
}

#' Required packages validator
#'
#' Checks if all the required packages are present according to metaseqr input 
#' options. Internal use only.
#'
#' @param m meta-analysis method
#' @param p qc plot types
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' check.packages(c("simes","whitlock"),c("gcbias","correl"))
#}
check.packages <- function(m,p) {
    # Check meta-analysis packages
    if (m=="whitlock" && !require(survcomp))
        stopwrap("Bioconductor package survcomp is required for \"whitlock\" ",
            "p-value meta analysis!")
    if ("venn" %in% p && !require(VennDiagram))
        stopwrap("R package VennDiagram is required for some of the selected ",
            "QC plots!")
}
