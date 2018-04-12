#' SAM/BAM/BED file reader helper for the metaseqr pipeline
#'
#' This function is a helper for the \code{metaseqr} pipeline, for reading SAM/BAM
#' or BED files when a read counts file is not available.
#'
#' @param targets a named list, the output of \code{\link{read.targets}}.
#' @param file.type the type of raw input files. It can be \code{"bed"} for BED
#' files or \code{"sam"}, \code{"bam"} for SAM/BAM files. See the same argument
#' in the main \code{\link{metaseqr}} function for the case of auto-guessing.
#' @param annotation see the \code{annotation} argument in the main
#' \code{\link{metaseqr}} function. The \code{"annotation"} parameter here is the
#' result of the same parameter in the main function. See also
#' \code{\link{get.annotation}}.
#' @param has.all.fields a logical variable indicating if all annotation fields
#' used by \code{metaseqr} are available (that is apart from the main chromosome,
#' start, end, unique id and strand columns, if also present are the gene name and
#' biotype columns). The default is \code{FALSE}.
#' @param multic a logical value indicating the presence of multiple cores. Defaults
#' to \code{FALSE}. Do not change it if you are not sure whether package parallel
#' has been loaded or not.
#' @return A data frame with counts for each sample, ready to be passed to the
#' main \code{\link{metaseqr}} pipeline.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' my.targets <- read.targets("my_mm9_study_bam_files.txt")
#' sample.list <- my.targets$samples
#' file.list <- my.targets$files
#' gene.data <- get.annotation("mm9","gene")
#' r2c <- read2count(files.list=file.list,file.type=my.targets$type,
#'   annotation=gene.data)
#' gene.counts <- r2c$counts
#' libsize.list <- r2s$libsize
#'}
read2count <- function(targets,annotation,file.type=targets$type,
	trans.level="gene",utr.flank=500,has.all.fields=FALSE,multic=FALSE) {
    if (missing(targets))
        stopwrap("You must provide the targets argument!")
    if (missing(annotation))
        stopwrap("You must provide an annotation data frame!")
    if (!require(GenomicRanges))
        stopwrap("The Bioconductor package GenomicRanges is required to ",
            "proceed!")
    if (file.type=="bed" && !require(rtracklayer))
        stopwrap("The Bioconductor package rtracklayer is required to process ",
            "BED files!")
    if (file.type %in% c("sam","bam")) {
        if (!require(Rsamtools))
            stopwrap("The Bioconductor package Rsamtools is required to ",
                "process BAM files!")
    }
    if (!is.list(targets) && file.exists(targets))
        targets <- read.targets(targets)
    else if (!is.list(targets) && !file.exists(targets))
        stopwrap("You must provide a targets list or a valid targets file!")
    
    # Convert annotation to GRanges
    disp("Converting annotation to GenomicRanges object...")
    if (packageVersion("GenomicRanges")<1.14) { # Classic way
        if (has.all.fields)
            annotation.gr <- GRanges(
                seqnames=Rle(annotation[,1]),
                ranges=IRanges(start=annotation[,2],end=annotation[,3]),
                strand=Rle(annotation[,6]),
                name=as.character(annotation[,4]),
                symbol=as.character(annotation[,7]),
                biotype=as.character(annotation[,8])
            )
        else
            annotation.gr <- GRanges(
                seqnames=Rle(annotation[,1]),
                ranges=IRanges(start=annotation[,2],end=annotation[,3]),
                strand=Rle(annotation[,6]),
                name=as.character(annotation[,4])
            )
    }
    else # Use native method in newer versions of GenomicRanges
        annotation.gr <- makeGRangesFromDataFrame(
            df=annotation,
            keep.extra.columns=TRUE,
            seqnames.field="chromosome"
        )
    # If the count type is "exon", we must reduce the overlapping exons belonging
    # to multiple transcripts, so as to avoid inflating the final read count when
    # summing all exons
    inter.feature <- FALSE
    if (length(grep("exon",colnames(annotation)))>0) { # count.type is exon
        if (length(grep("MEX",annotation$exon_id[1]))) # Retrieved from previous
            merged.annotation <- annotation
        else {
			if (trans.level=="gene") {
				disp("Merging exons to create unique gene models...")
				annotation.gr <- reduce.exons(annotation.gr,multic=multic)
				#merged.annotation <- as.data.frame(annotation.gr) # Bug?
			}
			#else if (trans.level=="transcript") {
			#	disp("Merging exons to create unique gene models...")
			#	annotation.gr <- reduce.exons.transcript(annotation.gr,
			#		multic=multic)
			#}
            merged.annotation <- data.frame(
                chromosome=as.character(seqnames(annotation.gr)),
                start=start(annotation.gr),
                end=end(annotation.gr),
                exon_id=if (!is.null(annotation.gr$exon_id))
                    as.character(annotation.gr$exon_id) else
                    as.character(annotation.gr$name),
                gene_id=if (!is.null(annotation.gr$gene_id))
                    as.character(annotation.gr$gene_id) else
                    as.character(annotation.gr$name),
                strand=as.character(strand(annotation.gr)),
                gene_name=if (!is.null(annotation.gr$gene_name))
                    as.character(annotation.gr$gene_name) else 
                    if (!is.null(annotation.gr$symbol))
                    as.character(annotation.gr$name) else NULL,
                biotype=if (!is.null(annotation.gr$biotype))
                    as.character(annotation.gr$biotype) else NULL
            )
            rownames(merged.annotation) <- 
                as.character(merged.annotation$exon_id)
        }
    }
    else if (length(grep("transcript",colnames(annotation)))>0) { # count.type may be utr
        if (length(grep("MET",annotation$transcript_id[1]))
			|| length(grep("MEU",annotation$transcript_id[1]))) # Retrieved from previous
            merged.annotation <- annotation
        else {
			if (trans.level=="gene") {
				disp("Merging transcript 3' UTRs to create unique ",
					"gene models...")
				annotation.gr <- 
					reduce.transcripts.utr(annotation.gr,multic=multic)
			}
			if (trans.level=="transcript") {
				disp("Merging transcript 3' UTRs to create unique ",
					"transcript models...")
				annotation.gr <- 
					reduce.transcripts.utr.transcript(annotation.gr,
						multic=multic)
			}
			if (utr.flank > 0) {
				disp("Flanking merged transcript 3' UTRs per ",utr.flank,
					"bp...")
				w <- width(annotation.gr)
				annotation.gr <- promoters(annotation.gr,upstream=utr.flank,
					downstream=0)
				annotation.gr <- resize(annotation.gr,width=w+2*utr.flank)
			}
			#merged.annotation <- as.data.frame(annotation.gr) # Bug?
			merged.annotation <- data.frame(
				chromosome=as.character(seqnames(annotation.gr)),
				start=start(annotation.gr),
				end=end(annotation.gr),
				transcript_id=if (!is.null(annotation.gr$transcript_id))
					as.character(annotation.gr$transcript_id) else
					as.character(annotation.gr$name),
				gene_id=if (!is.null(annotation.gr$gene_id))
					as.character(annotation.gr$gene_id) else
					as.character(annotation.gr$name),
				strand=as.character(strand(annotation.gr)),
				gene_name=if (!is.null(annotation.gr$gene_name))
					as.character(annotation.gr$gene_name) else 
					if (!is.null(annotation.gr$symbol))
					as.character(annotation.gr$name) else NULL,
				biotype=if (!is.null(annotation.gr$biotype))
					as.character(annotation.gr$biotype) else NULL
			)
			rownames(merged.annotation) <- 
				as.character(merged.annotation$transcript_id)
        }
        inter.feature = FALSE # Quant-Seq
    }
    else
        merged.annotation <- NULL
    # Continue
    files.list <- targets$files
    sample.names <- unlist(lapply(files.list,names),use.names=FALSE)
    sample.files <- unlist(files.list,use.names=FALSE)
    names(sample.files) <- sample.names
    if (!is.null(targets$paired)) {
        paired <- unlist(targets$paired,use.names=FALSE)
        names(paired) <- sample.names
    }
    else
        paired <- NULL
    if (!is.null(targets$stranded)) {
        stranded <- unlist(targets$stranded,use.names=FALSE)
        names(stranded) <- sample.names
    }
    else
        stranded <- NULL
    counts <- matrix(0,nrow=length(annotation.gr),ncol=length(sample.names))
    if (length(grep("exon",colnames(annotation)))>0)
        rownames(counts) <- as.character(annotation.gr$exon_id)
    else if (length(grep("transcript",colnames(annotation)))>0)
        rownames(counts) <- as.character(annotation.gr$transcript_id)
    else
        rownames(counts) <- as.character(annotation.gr$gene_id)
    colnames(counts) <- sample.names
    libsize <- vector("list",length(sample.names))
    names(libsize) <- sample.names
    if (file.type=="bed") {
        ret.val <- wapply(multic,sample.names,function(n,sample.files) {
            disp("Reading bed file ",basename(sample.files[n]),
                " for sample with name ",n,". This might take some time...")
            bed <- import.bed(sample.files[n],trackLine=FALSE)
            disp("  Checking for chromosomes not present in the annotation...")
            bed <- bed[which(!is.na(match(as(seqnames(bed),"character"),
                seqlevels(annotation.gr))))]
            libsize <- length(bed)
            if (length(bed)>0) {
                disp("  Counting reads overlapping with given annotation...")
                counts <- countOverlaps(annotation.gr,bed)
            }
            else
                warnwrap(paste("No reads left after annotation chromosome ",
                    "presence check for sample ",n,sep=""))
            gc(verbose=FALSE)
            return(list(counts=counts,libsize=libsize))
        },sample.files)
    }
    else if (file.type %in% c("sam","bam")) {
        if (file.type=="sam") {
            for (n in sample.names) {
                dest <- file.path(dirname(sample.files[n]),n)
                disp("Converting sam file ",basename(sample.files[n]),
                    " to bam file ",basename(dest),"...")
                asBam(file=sample.files[n],destination=dest,overwrite=TRUE)
                sample.files[n] <- paste(dest,"bam",sep=".")
            }
        }
        ret.val <- wapply(multic,sample.names,function(n,sample.files,paired,
            stranded) {
            disp("Reading bam file ",basename(sample.files[n])," for sample ",
                "with name ",n,". This might take some time...")
            if (!is.null(paired)) {
                p <- tolower(paired[n])
                if (p=="single") {
                    singleEnd <- TRUE
                    fragments <- FALSE
                    asMates <- FALSE
                }
                else if (p=="paired") {
                    singleEnd <- FALSE
                    fragments <- FALSE
                    asMates <- TRUE
                }
                else if (p=="mixed") {
                    singleEnd <- FALSE
                    fragments <- TRUE
                    asMates <- TRUE
                }
                else {
                    warnwrap("Information regarding single- or paired-end ",
                        "reads is not correctly provided! Assuming single...")
                    singleEnd <- TRUE
                    fragments <- FALSE
                    asMates <- FALSE
                }
            }
            else {
                singleEnd <- TRUE
                fragments <- FALSE
                asMates <- FALSE
            }
            if (!is.null(stranded)) {
                s <- tolower(stranded[n])
                if (s %in% c("forward","reverse"))
                    ignore.strand <- FALSE
                else if (s=="no")
                    ignore.strand <- TRUE
                else {
                    warnwrap("Information regarding strandedness of the reads ",
                        "is not correctly provided! Assuming unstranded...")
                    ignore.strand <- TRUE
                }
            }
            else
                ignore.strand <- TRUE
            # Check remoteness
            if (length(grep("^(http|ftp)",sample.files[n],perl=TRUE))>=1) {
                reads <- as(readGAlignments(file=sample.files[n]),"GRanges")
                libsize <- length(reads)
                is.remote <- TRUE
            }
            else {
                reads <- BamFile(sample.files[n],asMates=asMates)
                libsize <- countBam(reads,
                param=ScanBamParam(scanBamFlag(isUnmappedQuery=FALSE)))$records
                is.remote <- FALSE
            }
            if (libsize>0) {
                disp("  Counting reads overlapping with given annotation...")
                if (singleEnd & !fragments)
                    disp("    ...for single-end reads...")
                else if (!singleEnd & !fragments)
                    disp("    ...for paired-end reads...")
                else if (!singleEnd & fragments)
                    disp("    ...for mixed single- and paired-end reads...")
                if (ignore.strand)
                    disp("    ...ignoring strandedness...")
                else {
                    disp("    ...assuming ",s," sequenced reads...")
                    if (s=="reverse")
                        strand(annotation.gr) <- ifelse(strand(
                            annotation.gr)=="+","-","+")
                }
                if (is.remote)
                    disp("    ...for remote BAM file... might take longer...")
                counts <- summarizeOverlaps(annotation.gr,reads,
                    singleEnd=singleEnd,fragments=fragments,
                    ignore.strand=ignore.strand,inter.feature=inter.feature)
                counts <- assays(counts)$counts
            }
            else
                warnwrap(paste("No reads left after annotation chromosome ",
                    "presence check for sample ",n,sep=""))
            gc(verbose=FALSE)
            return(list(counts=counts,libsize=libsize))
        },sample.files,paired,stranded)
    }
    for (i in 1:length(ret.val)) {
        counts[,i] <- ret.val[[i]]$counts
        libsize[[i]] <- ret.val[[i]]$libsize
    }
    
    return(list(counts=counts,libsize=libsize,mergedann=merged.annotation))
}

#' Merges exons to create a unique set of exons for each gene
#'
#' This function uses the \code{"reduce"} function of IRanges to construct virtual
#' unique exons for each gene, so as to avoid inflating the read counts for each
#' gene because of multiple possible transcripts. If the user wants transcripts
#' instead of genes, they should be supplied to the original annotation table.
#'
#' @param gr a GRanges object created from the supplied annotation (see also the
#' \code{\link{read2count}} and \code{\link{get.annotation}} functions.
#' @param multic a logical value indicating the presence of multiple cores. Defaults
#' to \code{FALSE}. Do not change it if you are not sure whether package parallel
#' has been loaded or not.
#' @return A GRanges object with virtual merged exons for each gene/transcript.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' require(GenomicRanges)
#' ann <- get.annotation("mm9","exon")
#' gr <- makeGRangesFromDataFrame(
#'  df=ann,
#'  keep.extra.columns=TRUE,
#'  seqnames.field="chromosome"
#' )
#' re <- reduce.exons(gr)
#'}
reduce.exons <- function(gr,multic=FALSE) {
    gene <- unique(as.character(gr$gene_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    red.list <- wapply(multic,gene,function(x,a,g,b) {
        tmp <- a[a$gene_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            exon_id=paste(x,"MEX",1:n,sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt)
    return(do.call("c",red.list))
}

#' Merges 3' UTR of transcripts to create a unique set of coordinates for each
#' transcript
#'
#' This function uses the \code{"reduce"} function of IRanges to construct virtual
#' unique transcripts for each gene, so as to avoid inflating the read counts for 
#' each gene because of multiple possibly overlaping 3' UTR starts/ends when using
#' metaseqR with QuantSeq protocol.
#'
#' @param gr a GRanges object created from the supplied annotation (see also the
#' \code{\link{read2count}} and \code{\link{get.annotation}} functions.
#' @param multic a logical value indicating the presence of multiple cores. Defaults
#' to \code{FALSE}. Do not change it if you are not sure whether package parallel
#' has been loaded or not.
#' @return A GRanges object with virtual merged exons for each gene/transcript.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' require(GenomicRanges)
#' ann <- get.annotation("mm9","exon")
#' gr <- makeGRangesFromDataFrame(
#'  df=ann,
#'  keep.extra.columns=TRUE,
#'  seqnames.field="chromosome"
#' )
#' re <- reduce.exons(gr)
#'}
reduce.transcripts.utr <- function(gr,multic=FALSE) {
    gene <- unique(as.character(gr$gene_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    red.list <- wapply(multic,gene,function(x,a,g,b) {
        tmp <- a[a$gene_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            transcript_id=paste(x,"MET",1:n,sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt)
    return(do.call("c",red.list))
}

reduce.transcripts.utr.transcript <- function(gr,multic=FALSE) {
    trans <- unique(as.character(gr$transcript_id))
    if (!is.null(gr$gene_name))
        gn <- gr$gene_name
    else
        gn <- NULL
    if (!is.null(gr$biotype))
        bt <- gr$biotype   
    else
        bt <- NULL
    red.list <- wapply(multic,trans,function(x,a,g,b) {
        tmp <- a[a$transcript_id==x]
        if (!is.null(g))
            gena <- as.character(tmp$gene_name[1])
        if (!is.null(b))
            btty <- as.character(tmp$biotype[1])
        merged <- reduce(tmp)
        n <- length(merged)
        meta <- DataFrame(
            transcript_id=paste(x,"MEU",1:n,sep="_"),
            gene_id=rep(x,n)
        )
        if (!is.null(g))
            meta$gene_name <- rep(gena,n)
        if (!is.null(b))
            meta$biotype <- rep(btty,n)
        mcols(merged) <- meta
        return(merged)
    },gr,gn,bt)
    return(do.call("c",red.list))
}

#' Creates sample list and BAM/BED file list from file
#'
#' Create the main sample list and determine the BAM/BED files for each sample
#' from an external file.
#'
#' @param input a tab-delimited file structured as follows: the first line of the
#' external tab delimited file should contain column names (names are not important).
#' The first column MUST contain UNIQUE sample names. The second column MUST contain
#' the raw BAM/BED files WITH their full path. Alternatively, the \code{path}
#' argument should be provided (see below). The third column MUST contain the
#' biological condition where each of the samples in the first column should belong
#' to. There is an optional fourth column which should contain the keywords 
#' \code{"single"} for single-end reads, \code{"paired"} for paired-end reads or
#' \code{"mixed"} for BAM files that contain bith paired- and single-end reads.
#' If this column is not provided, single-end reads will be assumed. There is an
#' optional fifth column which controls stranded read assignment. It should 
#' contain the keywords \code{"forward"} for a forward (5'->3') strand library 
#' construction protocol, \code{"reverse"} for a reverse (3'->5') strand library 
#' construction protocol, or \code{"no"} for unstranded/unknown protocol. If 
#' this column is not provided, unstranded reads will be assumed.
#' @param path an optional path where all the BED/BAM files are placed, to be
#' prepended to the BAM/BED file names in the targets file.
#' @return A named list with three members. The first member is a named list whose
#' names are the conditions of the experiments and its members are the samples
#' belonging to each condition. The second member is like the first, but this time
#' the members are named vectors whose names are the sample names and the vector
#' elements are full path to BAM/BED files. The third member is the guessed type
#' of the input files (BAM or BED). It will be used if not given in the main
#' \code{\link{read2count}} function.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' targets <- data.frame(sample=c("C1","C2","T1","T2"),
#'   filename=c("C1_raw.bam","C2_raw.bam","T1_raw.bam","T2_raw.bam"),
#'   condition=c("Control","Control","Treatment","Treatment"))
#' path <- "/home/chakotay/bam"
#' write.table(targets,file="targets.txt",sep="\t",row.names=F,quote="")
#' the.list <- read.targets("targets.txt",path=path)
#' sample.list <- the.list$samples
#' bamfile.list <- the.list$files
#'}
read.targets <- function(input,path=NULL) {
    if (missing(input) || !file.exists(input))
        stopwrap("The targets file should be a valid existing text file!")
    tab <- read.delim(input,strip.white=TRUE)
    samples <- as.character(tab[,1])
    conditions <- unique(as.character(tab[,3]))
    rawfiles <- as.character(tab[,2])
    if (!is.null(path)) {
        tmp <- dirname(rawfiles) # Test if there is already a path
        if (any(tmp=="."))
            rawfiles <- file.path(path,basename(rawfiles))
    }
    # Check if they exist!!!
    for (f in rawfiles) {
        if (!file.exists(f))
            stopwrap("Raw reads input file ",f," does not exist! Please check!")
    }
    if (length(samples) != length(unique(samples)))
        stopwrap("Sample names must be unique for each sample!")
    if (length(rawfiles) != length(unique(rawfiles)))
        stopwrap("File names must be unique for each sample!")
    sample.list <- vector("list",length(conditions))
    names(sample.list) <- conditions
    for (n in conditions)
        sample.list[[n]] <- samples[which(as.character(tab[,3])==n)]
    file.list <- vector("list",length(conditions))
    names(file.list) <- conditions
    for (n in conditions) {
        file.list[[n]] <- rawfiles[which(as.character(tab[,3])==n)]
        names(file.list[[n]]) <- samples[which(as.character(tab[,3])==n)]
    }
    if (ncol(tab)>3) { # Has info about single- or paired-end reads / strand
        if (ncol(tab)==4) { # Stranded or paired
            whats <- tolower(as.character(tab[,4]))
            if (!all(whats %in% c("yes","no","forward","reverse",
                "single","paired")))
                stopwrap("Unknown options for paired-end reads and/or ",
                    "strandedness in targets file.")
            what <- whats[1]
            if (what %in% c("single","paired")) {
                has.paired.info <- TRUE
                has.stranded.info <- FALSE
            }
            else {
                if (what %in% c("yes","no")) {
                    deprecated.warning("read.targets")
                    tmp <- as.character(tab[,4])
                    tmp[tmp=="yes"] <- "forward"
                    tab[,4] <- tmp
                    has.paired.info <- FALSE
                    has.stranded.info <- TRUE
                }
                if (what %in% c("forward","reverse","no")) {
                    has.paired.info <- FALSE
                    has.stranded.info <- TRUE
                }
            }
        }
        if (ncol(tab)==5) { # Both
            whats.paired <- tolower(as.character(tab[,4]))
            if (!all(whats.paired %in% c("single","paired","mixed")))
                stopwrap("Unknown option for type of reads (single, paired, ",
                    "mixed) in targets file.")
            whats.strand <- tolower(as.character(tab[,5]))
            if (!all(whats.strand %in% c("yes","no","forward","reverse")))
                stopwrap("Unknown option for read strandedness in targets file")
            if (any(whats.strand=="yes")) {
				deprecated.warning("read.targets")
				tmp <- as.character(tab[,5])
				tmp[tmp=="yes"] <- "forward"
				tab[,5] <- tmp
			}
            has.paired.info <- TRUE
            has.stranded.info <- TRUE
        }
        if (has.paired.info && !has.stranded.info) {
            paired.list <- vector("list",length(conditions))
            names(paired.list) <- conditions
            for (n in conditions) {
                paired.list[[n]] <- character(length(sample.list[[n]]))
                names(paired.list[[n]]) <- sample.list[[n]]
                for (nn in names(paired.list[[n]]))
                    paired.list[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),4])
            }
        }
        else
            paired.list <- NULL
        if (has.stranded.info && !has.paired.info) {
            stranded.list <- vector("list",length(conditions))
            names(stranded.list) <- conditions
            for (n in conditions) {
                stranded.list[[n]] <- character(length(sample.list[[n]]))
                names(stranded.list[[n]]) <- sample.list[[n]]
                for (nn in names(stranded.list[[n]]))
                    stranded.list[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),4])
            }
        }
        else
            stranded.list <- NULL
        if (has.stranded.info && has.paired.info) {
            stranded.list <- vector("list",length(conditions))
            names(stranded.list) <- conditions
            for (n in conditions) {
                stranded.list[[n]] <- character(length(sample.list[[n]]))
                names(stranded.list[[n]]) <- sample.list[[n]]
                for (nn in names(stranded.list[[n]]))
                    stranded.list[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),5])
            }
            paired.list <- vector("list",length(conditions))
            names(paired.list) <- conditions
            for (n in conditions) {
                paired.list[[n]] <- character(length(sample.list[[n]]))
                names(paired.list[[n]]) <- sample.list[[n]]
                for (nn in names(paired.list[[n]]))
                    paired.list[[n]][nn] <- as.character(tab[which(as.character(
                        tab[,1])==nn),4])
            }
        }
    }
    else
        paired.list <- stranded.list <- NULL
    # Guess file type based on only one of them
    tmp <- file.list[[1]][1]
    if (length(grep("\\.bam$",tmp,ignore.case=TRUE,perl=TRUE))>0)
        type <- "bam"
    else if (length(grep("\\.sam$",tmp,ignore.case=TRUE,perl=TRUE))>0)
        type <- "sam"
    else if (length(grep("\\.bed$",tmp,ignore.case=TRUE,perl=TRUE)>0))
        type <- "bed"
    else
        type <- NULL
    return(list(samples=sample.list,files=file.list,paired=paired.list,
        stranded=stranded.list,type=type))
}
