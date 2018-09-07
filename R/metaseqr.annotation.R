#' Annotation downloader
#'
#' This function connects to the EBI's Biomart service using the package biomaRt
#' and downloads annotation elements (gene co-ordinates, exon co-ordinates, gene
#' identifications, biotypes etc.) for each of the supported organisms. See the
#' help page of \code{\link{metaseqr}} for a list of supported organisms. The
#' function downloads annotation for an organism genes or exons.
#'
#' @param org the organism for which to download annotation.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb the online source to use to fetch annotation. It can be
#' \code{"ensembl"} (default), \code{"ucsc"} or \code{"refseq"}. In the later two
#' cases, an SQL connection is opened with the UCSC public databases.
#' @param multic a logical value indicating the presence of multiple cores. Defaults
#' to \code{FALSE}. Do not change it if you are not sure whether package parallel
#' has been loaded or not. It is used in the case of \code{type="exon"} to process
#' the return value of the query to the UCSC Genome Browser database.
#' @return A data frame with the canonical (not isoforms!) genes or exons of the
#' requested organism. When \code{type="genes"}, the data frame has the following
#' columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' When \code{type="exon"} the data frame has the following columns: chromosome,
#' start, end, exon_id, gene_id, strand, gene_name, biotype. The gene_id and exon_id
#' correspond to Ensembl gene and exon accessions respectively. The gene_name
#' corresponds to HUGO nomenclature gene names.
#' @note The data frame that is returned contains only "canonical" chromosomes
#' for each organism. It does not contain haplotypes or random locations and does
#' not contain chromosome M.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg19.genes <- get.annotation("hg19","gene","ensembl")
#' mm9.exons <- get.annotation("mm9","exon","ucsc")
#'}
get.annotation <- function(org,type,refdb="ensembl",multic=FALSE) {
    org <- tolower(org)
    if (type %in% c("utr","transcript") && refdb %in% c("ucsc","refseq")) {
        disp("Quant-Seq (utr) and transcript analysis is not yet supported ",
            "with UCSC or RefSeq annotation. Switching to Ensembl...")
        refdb <- "ensembl"
    }
    switch(refdb,
        ensembl = { return(get.ensembl.annotation(org,type)) },
        ucsc = { return(get.ucsc.annotation(org,type,refdb,multic)) },
        refseq = { return(get.ucsc.annotation(org,type,refdb,multic)) }
    )
}

#' Ensembl annotation downloader
#'
#' This function connects to the EBI's Biomart service using the package biomaRt
#' and downloads annotation elements (gene co-ordinates, exon co-ordinates, gene
#' identifications, biotypes etc.) for each of the supported organisms. See the
#' help page of \code{\link{metaseqr}} for a list of supported organisms. The
#' function downloads annotation for an organism genes or exons.
#'
#' @param org the organism for which to download annotation.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @return A data frame with the canonical (not isoforms!) genes or exons of the
#' requested organism. When \code{type="genes"}, the data frame has the following
#' columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' When \code{type="exon"} the data frame has the following columns: chromosome,
#' start, end, exon_id, gene_id, strand, gene_name, biotype. The gene_id and exon_id
#' correspond to Ensembl gene and exon accessions respectively. The gene_name
#' corresponds to HUGO nomenclature gene names.
#' @note The data frame that is returned contains only "canonical" chromosomes
#' for each organism. It does not contain haplotypes or random locations and does
#' not contain chromosome M.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg19.genes <- get.ensembl.annotation("hg19","gene")
#' mm9.exons <- get.ensembl.annotation("mm9","exon")
#'}
get.ensembl.annotation <- function(org,type) {
    if (org=="tair10")
        dat <- "plants_mart"
    else if (org=="bmori2")
        dat <- "metazoa_mart"
    else
        dat <- "ENSEMBL_MART_ENSEMBL"
        
    mart <- tryCatch({
            useMart(biomart=dat,host=get.host(org),dataset=get.dataset(org))
        },
        error=function(e) {
            useMart(biomart=dat,host=get.alt.host(org),
            dataset=get.dataset(org))
        },
        finally={})
    chrs.exp <- paste(get.valid.chrs(org),collapse="|")
    if (type=="gene") {
        bm <- getBM(attributes=get.gene.attributes(org),mart=mart)
        ann <- data.frame(
            chromosome=paste("chr",bm$chromosome_name,sep=""),
            start=bm$start_position,
            end=bm$end_position,
            gene_id=bm$ensembl_gene_id,
            gc_content=if (org %in% c("hg18","mm9","tair10")) 
				bm$percentage_gc_content else bm$percentage_gene_gc_content,
            strand=ifelse(bm$strand==1,"+","-"),
            gene_name=if (org %in% c("hg18","mm9","tair10")) bm$external_gene_id 
                else bm$external_gene_name,
            biotype=bm$gene_biotype
        )
        rownames(ann) <- ann$gene_id
    }
    else if (type=="exon") {
        bm <- getBM(attributes=get.exon.attributes(org),mart=mart)
        if (org == "hg19") {
            disp("  Bypassing problem with hg19 Ensembl combined gene-exon ",
                "annotation... Will take slightly longer...")
            bmg <- getBM(attributes=get.gene.attributes(org),mart=mart)
            gene_name <- bmg$external_gene_name
            names(gene_name) <- bmg$ensembl_gene_id
            ann <- data.frame(
                chromosome=paste("chr",bm$chromosome_name,sep=""),
                start=bm$exon_chrom_start,
                end=bm$exon_chrom_end,
                exon_id=bm$ensembl_exon_id,
                gene_id=bm$ensembl_gene_id,
                strand=ifelse(bm$strand==1,"+","-"),
                gene_name=gene_name[bm$ensembl_gene_id],
                biotype=bm$gene_biotype
            )
            rownames(ann) <- ann$exon_id
        }
        else
            ann <- data.frame(
                chromosome=paste("chr",bm$chromosome_name,sep=""),
                start=bm$exon_chrom_start,
                end=bm$exon_chrom_end,
                exon_id=bm$ensembl_exon_id,
                gene_id=bm$ensembl_gene_id,
                strand=ifelse(bm$strand==1,"+","-"),
                gene_name=if (org %in% c("hg18","mm9","tair10")) 
                    bm$external_gene_id else bm$external_gene_name,
                biotype=bm$gene_biotype
            )
            rownames(ann) <- ann$exon_id
    }
    else if (type=="utr") {
        bm <- getBM(attributes=get.transcript.utr.attributes(org),mart=mart)
        ann <- data.frame(
            chromosome=paste("chr",bm$chromosome_name,sep=""),
            start=bm$`3_utr_start`,
            end=bm$`3_utr_end`,
            tstart=bm$transcript_start,
            tend=bm$transcript_end,
            transcript_id=bm$ensembl_transcript_id,
            gene_id=bm$ensembl_gene_id,
            strand=ifelse(bm$strand==1,"+","-"),
            gene_name=if (org %in% c("hg18","mm9","tair10")) 
                bm$external_gene_id else bm$external_gene_name,
            biotype=bm$gene_biotype
        )
        ann <- correct.transcripts(ann)
        ann <- ann[,c("chromosome","start","end","transcript_id","gene_id",
            "strand","gene_name","biotype")]
    }
    else if (type=="transcript") {
        bm <- getBM(attributes=get.transcript.attributes(org),mart=mart)
        ann <- data.frame(
            chromosome=paste("chr",bm$chromosome_name,sep=""),
            start=bm$transcript_start,
            end=bm$transcript_end,
            transcript_id=bm$ensembl_transcript_id,
            gene_id=bm$ensembl_gene_id,
            strand=ifelse(bm$strand==1,"+","-"),
            gene_name=if (org %in% c("hg18","mm9","tair10")) 
                bm$external_gene_id else bm$external_gene_name,
            biotype=bm$gene_biotype
        )
        rownames(ann) <- as.character(ann$transcript_id)
    }
    ann <- ann[order(ann$chromosome,ann$start),]
    ann <- ann[grep(chrs.exp,ann$chromosome),]
    ann$chromosome <- as.character(ann$chromosome)
    return(ann)
}

#' UCSC/RefSeq annotation downloader
#'
#' This function connects to the UCSC Genome Browser public database and downloads
#' annotation elements (gene co-ordinates, exon co-ordinates, gene identifications
#' etc.) for each of the supported organisms, but using UCSC instead of Ensembl.
#' See the help page of \code{\link{metaseqr}} for a list of supported organisms.
#' The function downloads annotation for an organism genes or exons.
#'
#' @param org the organism for which to download annotation.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb either \code{"ucsc"} or \code{"refseq"}.
#' @param multic a logical value indicating the presence of multiple cores. Defaults
#' to \code{FALSE}. Do not change it if you are not sure whether package parallel
#' has been loaded or not. It is used in the case of \code{type="exon"} to process
#' the return value of the query to the UCSC Genome Browser database.
#' @return A data frame with the canonical (not isoforms!) genes or exons of the
#' requested organism. When \code{type="genes"}, the data frame has the following
#' columns: chromosome, start, end, gene_id, gc_content, strand, gene_name, biotype.
#' When \code{type="exon"} the data frame has the following columns: chromosome,
#' start, end, exon_id, gene_id, strand, gene_name, biotype. The gene_id and exon_id
#' correspond to UCSC or RefSeq gene and exon accessions respectively. The gene_name
#' corresponds to HUGO nomenclature gene names.
#' @note The data frame that is returned contains only "canonical" chromosomes
#' for each organism. It does not contain haplotypes or random locations and does
#' not contain chromosome M. Note also that as the UCSC databases do not contain
#' biotype classification like Ensembl, this will be returned as \code{NA} and
#' as a result, some quality control plots will not be available.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg19.genes <- get.ucsc.annotation("hg19","gene","ucsc")
#' mm9.exons <- get.ucsc.annotation("mm9","exon")
#'}
get.ucsc.annotation <- function(org,type,refdb="ucsc",multic=FALSE) {
    if (org=="bmori2") {
        warnwrap("Bombyx mori (silkworm) annotation is not supported by UCSC ",
            "or RefSeq! Will use Ensembl...")
        return(get.ensembl.annotation("bmori2",type))
    }
    if (!require(RMySQL)) {
        rmysql.present <- FALSE
        warnwrap("R package RMySQL is not present! Annotation will be ",
            "retrieved by downloading temporary files from UCSC and the usage
            of a temporary SQLite database...")
    }
    else
        rmysql.present <- TRUE
    if (!require(RSQLite))
        stopwrap("R package RSQLite is required to use annotation from UCSC!")

    if (org=="tair10") {
        warnwrap("Arabidopsis thaliana genome is not supported by UCSC Genome ",
            "Browser database! Switching to Ensembl...")
        return(get.ensembl.annotation("tair10",type))
    }
    if (org=="equcab2") {
        warnwrap("Equus cabalus genome is not supported by UCSC Genome ",
            "Browser database! Switching to Ensembl...")
        return(get.ensembl.annotation("equcab2",type))
    }

    valid.chrs <- get.valid.chrs(org)
    chrs.exp <- paste("^",paste(valid.chrs,collapse="$|^"),"$",sep="")

    db.org <- get.ucsc.organism(org)
    if (rmysql.present) {
        db.creds <- get.ucsc.credentials()
        drv <- dbDriver("MySQL")
        con <- dbConnect(drv,user=db.creds[2],password=NULL,dbname=db.org,
            host=db.creds[1])
        query <- get.ucsc.query(org,type,refdb)
        raw.ann <- dbGetQuery(con,query)
        dbDisconnect(con)
    }
    else {
        # This should return the same data frame as the db query
        tmp.sqlite <- get.ucsc.dbl(org,type,refdb)
        drv <- dbDriver("SQLite")
        con <- dbConnect(drv,dbname=tmp.sqlite)
        query <- get.ucsc.query(org,type,refdb)
        raw.ann <- dbGetQuery(con,query)
        dbDisconnect(con)
    }
    if (type=="gene") {
        ann <- raw.ann
        ann <- ann[grep(chrs.exp,ann$chromosome,perl=TRUE),]
        ann$chromosome <- as.character(ann$chromosome)
        rownames(ann) <- ann$gene_id
    }
    else if (type=="exon") {
        raw.ann <- raw.ann[grep(chrs.exp,raw.ann$chromosome,perl=TRUE),]
        ex.list <- wapply(multic,as.list(1:nrow(raw.ann)),function(x,d,s) {
            r <- d[x,]
            starts <- as.numeric(strsplit(r[,"start"],",")[[1]])
            ends <- as.numeric(strsplit(r[,"end"],",")[[1]])
            nexons <- length(starts)
            ret <- data.frame(
                rep(r[,"chromosome"],nexons),
                starts,ends,
                paste(r[,"exon_id"],"_e",1:nexons,sep=""),
                rep(r[,"strand"],nexons),
                rep(r[,"gene_id"],nexons),
                rep(r[,"gene_name"],nexons),
                rep(r[,"biotype"],nexons)
            )
            names(ret) <- names(r)
            rownames(ret) <- ret$exon_id
            ret <- makeGRangesFromDataFrame(
                df=ret,
                keep.extra.columns=TRUE,
                seqnames.field="chromosome",
                seqinfo=s
            )
            return(ret)
        },raw.ann,valid.chrs)
        tmp.ann <- do.call("c",ex.list)
        ann <- data.frame(
            chromosome=as.character(seqnames(tmp.ann)),
            start=start(tmp.ann),
            end=end(tmp.ann),
            exon_id=as.character(tmp.ann$exon_id),
            gene_id=as.character(tmp.ann$gene_id),
            strand=as.character(strand(tmp.ann)),
            gene_name=as.character(tmp.ann$gene_name),
            biotype=tmp.ann$biotype
        )
        rownames(ann) <- ann$exon_id
    }
    
    gc.content <- get.gc.content(ann,org)
    ann$gc_content <- gc.content
    ann <- ann[order(ann$chromosome,ann$start),]
    return(ann)
}

#' Return a named vector of GC-content for each genomic region
#'
#' Returns a named numeric vector (names are the genomic region names, e.g. genes)
#' given a data frame which can be converted to a GRanges object (e.g. it has at
#' least chromosome, start, end fields). This function works best when the input
#' annotation data frame has been retrieved using one of the SQL queries generated
#' from \code{\link{get.ucsc.query}}, used in \code{\link{get.ucsc.annotation}}.
#'
#' @param ann a data frame which can be converted to a GRanges object, that means
#' it has at least the chromosome, start, end fields. Preferably, the output of
#' \code{link{get.ucsc.annotation}}.
#' @param org one of metaseqR supported organisms.
#' @return A named numeric vector.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' ann <- get.ucsc.annotation("mm9","gene","ucsc")
#' gc <- get.gc.content(ann,"mm9")
#'}
get.gc.content <- function(ann,org) {
    if (missing(ann))
        stopwrap("A valid annotation data frame must be provided in order to ",
            "retrieve GC-content.")
    org <- tolower(org[1])
    check.text.args("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
		"dm6","danrer7","pantro4","susscr3","tair10"),multiarg=FALSE)
    # Convert annotation to GRanges
    disp("Converting annotation to GenomicRanges object...")
    if (packageVersion("GenomicRanges")<1.14)
        ann.gr <- GRanges(
            seqnames=Rle(ann[,1]),
            ranges=IRanges(start=ann[,2],end=ann[,3]),
            strand=Rle(ann[,6]),
            name=as.character(ann[,4])
        )
    else
        ann.gr <- makeGRangesFromDataFrame(
            df=ann,
            keep.extra.columns=TRUE,
            seqnames.field="chromosome"
        )
    bsg <- load.bs.genome(org)
    disp("Getting DNA sequences...")
    seqs <- getSeq(bsg,names=ann.gr)
    disp("Getting GC content...")
    freq.matrix <- alphabetFrequency(seqs,as.prob=TRUE,baseOnly=TRUE)
    gc.content <- apply(freq.matrix,1,function(x) round(100*sum(x[2:3]),
        digits=2))
    names(gc.content) <- as.character(ann[,4])
    return(gc.content)
}

#' Return a proper formatted organism alias
#'
#' Returns the proper UCSC Genome Browser database organism alias based on what is
#' given to metaseqR. Internal use.
#'
#' @return A proper organism alias.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' org <- get.ucsc.organism("danrer7")
#'}
get.ucsc.organism <- function(org) {
    switch(org,
        hg18 = { return("hg18") },
        hg19 = { return("hg19") },
        hg38 = { return("hg38") },
        mm9 = { return("mm9") },
        mm10 = { return("mm10") },
        rn5 = { return("rn5") },
        rn6 = { return("rn6") },
        dm3 = { return("dm3") },
        dm6 = { return("dm6") },
        danrer7 = { return("danRer7") },
        pantro4 = { return("panTro4") },
        susscr3 = { return("susScr3") },
        tair10 = { return("TAIR10") },
        equcab2 = { return("equCab2") }
    )
}

#' Return a proper formatted BSgenome organism name
#'
#' Returns a properly formatted BSgenome package name according to metaseqR's
#' supported organism. Internal use.
#'
#' @return A proper BSgenome package name.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' bs.name <- get.bs.organism("hg18")
#'}
get.bs.organism <- function(org) {
    switch(org,
        hg18 = {
            return("BSgenome.Hsapiens.UCSC.hg18")
        },
        hg19 = {
            return("BSgenome.Hsapiens.UCSC.hg19")
        },
        hg38 = {
            return("BSgenome.Hsapiens.UCSC.hg38")
        },
        mm9 = {
            return("BSgenome.Mmusculus.UCSC.mm9")
        },
        mm10 = {
            return("BSgenome.Mmusculus.UCSC.mm10")
        },
        rn5 = {
            return("BSgenome.Rnorvegicus.UCSC.rn5")
        },
        rn6 = {
            return("BSgenome.Rnorvegicus.UCSC.rn6")
        },
        dm3 = {
            return("BSgenome.Dmelanogaster.UCSC.dm3")
        },
        dm6 = {
            return("BSgenome.Dmelanogaster.UCSC.dm6")
        },
        danrer7 = {
            return("BSgenome.Drerio.UCSC.danRer7")
        },
        pantro4 = {
            stopwrap("panTro4 is not yet supported by BSgenome! Please use ",
                "Ensembl as annoation source.")
        },
        susscr3 = {
            return("BSgenome.Sscrofa.UCSC.susScr3")
        },
        tair10 = {
            stopwrap("TAIR10 is not yet supported by BSgenome! Please use ",
                "Ensembl as annoation source.")
        }
    )
}

#' Loads (or downloads) the required BSGenome package
#'
#' Retrieves the required BSgenome package when the annotation source is \code{"ucsc"}
#' or \code{"refseq"}. These packages are required in order to estimate the
#' GC-content of the retrieved genes from UCSC or RefSeq.
#'
#' @param org one of \code{\link{metaseqr}} supported organisms.
#' @return The BSgenome object for the requested organism.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' bs.obj <- load.bs.genome("mm9")
#'}
load.bs.genome <- function(org) {
    if (!require(BiocManager))
        stopwrap("The Bioconductor package BiocManager is required to ",
            "proceed!")
    if (!require(BSgenome))
        stopwrap("The Bioconductor package BSgenome is required to ",
            "proceed!")
    bs.org <- get.bs.organism(org)
    if (bs.org %in% installed.genomes())
        bs.obj <- getBSgenome(get.ucsc.organism(org))
    else {
        BiocManager::install(bs.org)
        bs.obj <- getBSgenome(get.ucsc.organism(org))
    }
    return(bs.obj)
}

#' Biotype converter
#'
#' Returns biotypes as character vector. Internal use.
#'
#' @param a the annotation data frame (output of \code{\link{get.annotation}}).
#' @return A character vector of biotypes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg18.genes <- get.annotation("hg18","gene")
#' hg18.bt <- get.biotypes(hg18.genes)
#'}
get.biotypes <- function(a) {
    return(as.character(unique(a$biotype)))
}

#' Annotation downloader helper
#'
#' Returns the appropriate Ensembl host address to get different versions of
#' annotation from. Internal use.
#'
#' @param org the organism for which to return the host address.
#' @return A string with the host address.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' mm9.host <- get.host("mm9")
#'}
get.host <- function(org) {
    switch(org,
        hg18 = { return("may2009.archive.ensembl.org") },
        hg19 = { return("grch37.ensembl.org") },
        hg38 = { return("www.ensembl.org") },
        mm9 = { return("may2012.archive.ensembl.org") },
        mm10 = { return("www.ensembl.org") },
        rn5 = { return("may2012.archive.ensembl.org") },
        rn6 = { return("www.ensembl.org") },
        dm3 = { return("grch37.ensembl.org") },
        dm6 = { return("www.ensembl.org") },
        danrer7 = { return("www.ensembl.org") },
        pantro4 = { return("www.ensembl.org") },
        susscr3 = { return("www.ensembl.org") },
        tair10 = { return("plants.ensembl.org") },
        equcab2 = { return("www.ensembl.org") },
        bmori2 = { return("metazoa.ensembl.org") }
    )
}

#' Annotation downloader helper
#'
#' Returns the appropriate Ensembl host address to get different versions of
#' annotation from (alternative hosts). Internal use.
#'
#' @param org the organism for which to return the host address.
#' @return A string with the host address.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' mm9.host <- get.alt.host("mm9")
#'}
get.alt.host <- function(org) {
    switch(org,
        hg18 = { return("may2009.archive.ensembl.org") },
        hg19 = { return("grch37.ensembl.org") },
        hg38 = { return("uswest.ensembl.org") },
        mm9 = { return("may2012.archive.ensembl.org") },
        mm10 = { return("uswest.ensembl.org") },
        rn5 = { return("uswest.ensembl.org") },
        dm3 = { return("grch37.ensembl.org") },
        dm6 = { return("uswest.ensembl.org") },
        danrer7 = { return("uswest.ensembl.org") },
        pantro4 = { return("uswest.ensembl.org") },
        susscr3 = { return("uswest.ensembl.org") },
        tair10 = { return("plants.ensembl.org") },
        equcab2 = { return("uswest.ensembl.org") },
        bmori2 = { return("metazoa.ensembl.org") }
    )
}

#' Annotation downloader helper
#'
#' Returns a dataset (gene or exon) identifier for each organism recognized by
#' the Biomart service for Ensembl. Internal use.
#'
#' @param org the organism for which to return the identifier.
#' @return A string with the dataset identifier.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' dm6.id <- get.dataset("dm6")
#'}
get.dataset <- function(org) {
    switch(org,
        hg18 = { return("hsapiens_gene_ensembl") },
        hg19 = { return("hsapiens_gene_ensembl") },
        hg38 = { return("hsapiens_gene_ensembl") },
        mm9 = { return("mmusculus_gene_ensembl") },
        mm10 = { return("mmusculus_gene_ensembl") },
        rn5 = { return("rnorvegicus_gene_ensembl") },
        rn6 = { return("rnorvegicus_gene_ensembl") },
        dm3 = { return("dmelanogaster_gene_ensembl") },
        dm6 = { return("dmelanogaster_gene_ensembl") },
        danrer7 = { return("drerio_gene_ensembl") },
        pantro4 = { return("ptroglodytes_gene_ensembl") },
        susscr3 = { return("sscrofa_gene_ensembl") },
        tair10 = { return("athaliana_eg_gene") },
        equcab2 = { return("ecaballus_gene_ensembl") },
        bmori2 = { return("bmori_eg_gene") },
    )
}

#' Annotation downloader helper
#'
#' Returns a vector of chromosomes to maintain after annotation download. Internal
#' use.
#'
#' @param org the organism for which to return the chromosomes. 
#' @return A character vector of chromosomes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' hg18.chr <- get.valid.chrs("hg18")
#'}
get.valid.chrs <- function(org)
{
    switch(org,
        hg18 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg19 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        hg38 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3",
                "chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm9 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        mm10 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX","chrY"
            ))
        },
        rn5 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX"
            ))
        },
        rn6 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr3","chr4","chr5","chr6",
                "chr7","chr8","chr9","chrX"
            ))
        },
        dm3 = {
            return(c(
                "chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet",
                "chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet",
                "chrYHet"
            ))
        },
        dm6 = {
            return(c(
                "chr2L","chr2LHet","chr2R","chr2RHet","chr3L","chr3LHet",
                "chr3R","chr3RHet","chr4","chrU","chrUextra","chrX","chrXHet",
                "chrYHet"
            ))
        },
        danrer7 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr3","chr4","chr5","chr6","chr7","chr8","chr9"
            ))
        },
        pantro4 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr20","chr21","chr22","chr2A","chr2B",
                "chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        susscr3 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr2","chr3","chr4","chr5","chr6","chr7",
                "chr8","chr9","chrX","chrY"
            ))
        },
        tair10 = {
            return(c(
                "chr1","chr2","chr3","chr4","chr5"
            ))
        },
        equcab2 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr26","chr27","chr28","chr29","chr3","chr30",
                "chr31","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"
            ))
        },
        bmori2 = {
            return(c(
                "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                "chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23",
                "chr24","chr25","chr26","chr27","chr28","chr3","chr4","chr5",
                "chr6","chr7","chr8","chr9"
            ))
        }
    )
}

#' Annotation downloader helper
#'
#' Returns a vector of genomic annotation attributes which are used by the biomaRt
#' package in order to fetch the gene annotation for each organism. It has no
#' parameters. Internal use.
#'
#' @param org one of the supported organisms.
#' @return A character vector of Ensembl gene attributes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' gene.attr <- get.gene.attributes()
#'}
get.gene.attributes <- function(org) {
    if (org %in% c("hg18","mm9","tair10","bmori2"))
        return(c(
            "chromosome_name",
            "start_position",
            "end_position",
            "ensembl_gene_id",
            "percentage_gc_content",
            "strand",
            "external_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "start_position",
            "end_position",
            "ensembl_gene_id",
            "percentage_gene_gc_content",
            "strand",
            "external_gene_name",
            "gene_biotype"
        ))
}

#' Annotation downloader helper
#'
#' Returns a vector of genomic annotation attributes which are used by the biomaRt
#' package in order to fetch the exon annotation for each organism. It has no
#' parameters. Internal use.
#'
#' @param org one of the supported organisms.
#' @return A character vector of Ensembl exon attributes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' exon.attr <- get.exon.attributes()
#'}
get.exon.attributes <- function(org) {
    if (org %in% c("hg18","mm9","tair10","bmori2"))
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_id",
            "gene_biotype"
        ))
    else if (org == "hg19")
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "exon_chrom_start",
            "exon_chrom_end",
            "ensembl_exon_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_name",
            "gene_biotype"
        ))
}

#' Annotation downloader helper
#'
#' Returns a vector of genomic annotation attributes which are used by the biomaRt
#' package in order to fetch the exon annotation for each organism. It has no
#' parameters. Internal use.
#'
#' @param org one of the supported organisms.
#' @return A character vector of Ensembl transcript attributes.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' trans.attr <- get.transcript.attributes()
#'}
get.transcript.utr.attributes <- function(org) {
    if (org %in% c("hg18","mm9","tair10","bmori2"))
        return(c(
            "chromosome_name",
            "transcript_start",
            "transcript_end",
            "3_utr_start",
            "3_utr_end",
            "ensembl_transcript_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "transcript_start",
            "transcript_end",
            "3_utr_start",
            "3_utr_end",
            "ensembl_transcript_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_name",
            "gene_biotype"
        ))
}

get.transcript.attributes <- function(org) {
    if (org %in% c("hg18","mm9","tair10","bmori2"))
        return(c(
            "chromosome_name",
            "transcript_start",
            "transcript_end",
            "ensembl_transcript_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_id",
            "gene_biotype"
        ))
    else
        return(c(
            "chromosome_name",
            "transcript_start",
            "transcript_end",
            "ensembl_transcript_id",
            "strand",
            "ensembl_gene_id",
            "external_gene_name",
            "gene_biotype"
        ))
}

correct.transcripts <- function(ann) {
    rownames(ann) <- paste("T",1:nrow(ann),sep="_")
    len <- ann[,3] - ann[,2]
    len <- len[-which(is.na(len))]
    len[len==0] <- 1
    def.utr.len <- round(2^mean(log2(len)))
    nas <- which(is.na(ann$start))
    ann.na <- ann[nas,]
    ann.na$start <- ann.na$tstart
    ann.na$end <- ann.na$tend
    tmp <- makeGRangesFromDataFrame(df=ann.na)
    tmp <- flank(resize(tmp,width=1,fix="end"),start=FALSE,width=def.utr.len)
    ann[names(tmp),"start"] <- start(tmp)
    ann[names(tmp),"end"] <- end(tmp)
    return(ann)
}
