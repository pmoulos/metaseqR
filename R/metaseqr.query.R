#' Download annotation from UCSC servers, according to organism and source
#'
#' Directly downloads UCSC and RefSeq annotation files from UCSC servers to be
#' used with metaseqR. This functionality is used when the package RMySQL is not
#' available for some reason, e.g. Windows machines.
#'
#' @param org one of metaseqR supported organisms.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb one of \code{"ucsc"} or \code{"refseq"} to use the UCSC or RefSeq
#' annotation sources respectively.
#' @return A data frame with annotation elements.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.file <- get.ucsc.file("hg18","gene","ucsc")
#'}
get.ucsc.dbl <- function(org,type,refdb="ucsc") {
    type <- tolower(type[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    check.text.args("type",type,c("gene","exon"))
    check.text.args("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "dm6","danrer7","pantro4","susscr3","tair10"),multiarg=FALSE)
    check.text.args("refdb",refdb,c("ucsc","refseq"))
    
    if (!require(RSQLite))
        stopwrap("R package RSQLite is required to use annotation from UCSC!")

    http.base <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/",
        get.ucsc.organism(org),"/database/",sep="")
    table.defs <- get.ucsc.tabledef(org,type,refdb,"fields")
    file.list <- vector("list",length(table.defs))
    names(file.list) <- names(table.defs)
    for (n in names(file.list))
        file.list[[n]] <- paste(http.base,n,".txt.gz",sep="")
        
    # Fill the fields for each table
    drv <- dbDriver("SQLite")
    db.tmp <- tempfile()
    con <- dbConnect(drv,dbname=db.tmp)
    #disp("  Defining tables for temporary SQLite ",refdb," ",org," ",
    #    type," subset database")
    #for (n in names(file.list)) {
    #    disp("    Creating table ",n,"\n")
    #    dbSendQuery(con,table.defs[[n]])
    #}
    disp("  Retrieving tables for temporary SQLite ",refdb," ",org," ",type,
        " subset database")
    for (n in names(file.list)) {
        disp("    Retrieving table ",n)
        download.file(file.list[[n]],file.path(tempdir(),
            paste(n,".txt.gz",sep="")),quiet=TRUE)
        if (.Platform$OS.type == "unix")
            system(paste("gzip -df",file.path(tempdir(),
                paste(n,".txt.gz",sep=""))))
        else
            unzip(file.path(tempdir(),paste(n,".txt.gz",sep="")))
        sql.df <- read.delim(file.path(tempdir(),paste(n,".txt",sep="")),
            row.names=NULL,header=FALSE,strip.white=TRUE)
        names(sql.df) <- table.defs[[n]]
        dbWriteTable(con,n,sql.df,row.names=FALSE)
    }
    dbDisconnect(con)
    return(db.tmp)
}

#' Get SQLite UCSC table defintions, according to organism and source
#'
#' Creates a list of UCSC Genome Browser database tables and their SQLite
#' definitions with the purpose of creating a temporary SQLite database to be 
#' used used with metaseqR. This functionality is used when the package RMySQL 
#' is not available for some reason, e.g. Windows machines.
#'
#' @param org one of metaseqR supported organisms.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb one of \code{"ucsc"} or \code{"refseq"} to use the UCSC or RefSeq
#' annotation sources respectively.
#' @return A list with SQLite table definitions.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.tabledefs <- get.ucsc.tables("hg18","gene","ucsc")
#'}
get.ucsc.tabledef <- function(org,type,refdb="ucsc",what="queries") {
    type <- tolower(type[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    what <- tolower(what[1])
    check.text.args("type",type,c("gene","exon"))
    check.text.args("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "dm6","danrer7","pantro4","susscr3","tair10"),multiarg=FALSE)
    check.text.args("refdb",refdb,c("ucsc","refseq"))
    check.text.args("what",what,c("queries","fields"))
    switch(type,
        gene = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(list(
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownGene=get.ucsc.tbl.tpl("knownGene",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                refFlat=get.ucsc.tbl.tpl("refFlat",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownGene=get.ucsc.tbl.tpl("knownGene",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what),
                                refFlat=get.ucsc.tbl.tpl("refFlat",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownGene=get.ucsc.tbl.tpl("knownGene",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                refFlat=get.ucsc.tbl.tpl("refFlat",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownGene=get.ucsc.tbl.tpl("knownGene",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what),
                                refFlat=get.ucsc.tbl.tpl("refFlat",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownGene=get.ucsc.tbl.tpl("knownGene",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what),
                                refFlat=get.ucsc.tbl.tpl("refFlat",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                mgcGenes=get.ucsc.tbl.tpl("mgcGenes",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                flyBaseCanonical=
                                    get.ucsc.tbl.tpl("flyBaseCanonical",what),
                                flyBaseGene=
                                    get.ucsc.tbl.tpl("flyBaseGene",what),
                                flyBaseToRefSeq=
                                    get.ucsc.tbl.tpl("flyBaseToRefSeq",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        dm6 = {
                            # Stub until we find out what will happen with
                            # Augustus
                        },
                        danrer7 = {
                            return(list(
                                mgcGenes=get.ucsc.tbl.tpl("mgcGenes",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            warnwrap("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                now=TRUE)
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            warnwrap("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                now=TRUE)
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        tair10 = {
                            warnwrap("Arabidopsis thaliana genome is not ",
                                "supported by UCSC Genome Borwser database! ",
                                "Will automatically switch to Ensembl...",
                                now=TRUE)
                            return(FALSE)
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        dm6 = {
                            # Stub until we find out what is going on with
                            # Augustus
                        },
                        danrer7 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        tair10 = {
                            warnwrap("Arabidopsis thaliana genome is not ",
                                "supported by UCSC Genome Borwser database! ",
                                "Will automatically switch to Ensembl...",
                                now=TRUE)
                            return(FALSE)
                        }
                    )
                }
            )
        },
        exon = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(list(
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownGene=get.ucsc.tbl.tpl("knownGene",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                refFlat=get.ucsc.tbl.tpl("refFlat",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownGene=get.ucsc.tbl.tpl("knownGene",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what),
                                refFlat=get.ucsc.tbl.tpl("refFlat",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownGene=get.ucsc.tbl.tpl("knownGene",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                refFlat=get.ucsc.tbl.tpl("refFlat",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownGene=get.ucsc.tbl.tpl("knownGene",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what),
                                refFlat=get.ucsc.tbl.tpl("refFlat",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownGene=get.ucsc.tbl.tpl("knownGene",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what),
                                refFlat=get.ucsc.tbl.tpl("refFlat",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                mgcGenes=get.ucsc.tbl.tpl("mgcGenes",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                flyBaseCanonical=
                                    get.ucsc.tbl.tpl("flyBaseCanonical",what),
                                flyBaseGene=
                                    get.ucsc.tbl.tpl("flyBaseGene",what),
                                flyBaseToRefSeq=
                                    get.ucsc.tbl.tpl("flyBaseToRefSeq",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        dm6 = {
                            # Stub until we find out what is going on with
                            # Augustus
                        },
                        danrer7 = {
                            return(list(
                                mgcGenes=get.ucsc.tbl.tpl("mgcGenes",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            warnwrap("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                now=TRUE)
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            warnwrap("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                now=TRUE)
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        tair10 = {
                            warnwrap("Arabidopsis thaliana genome is not ",
                                "supported by UCSC Genome Borwser database! ",
                                "Will automatically switch to Ensembl...",
                                now=TRUE)
                            return(FALSE)
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what)
                            ))
                        },
                        hg19 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        hg38 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what)
                            ))
                        },
                        mm9 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        mm10 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                knownToRefSeq=
                                    get.ucsc.tbl.tpl("knownToRefSeq",what),
                                knownCanonical=
                                    get.ucsc.tbl.tpl("knownCanonical",what),
                                knownToEnsembl=
                                    get.ucsc.tbl.tpl("knownToEnsembl",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        rn5 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        dm3 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        dm6 = {
                            # Stub until we find out what is going on with
                            # Augustus
                        },
                        danrer7 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        pantro4 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        susscr3 = {
                            return(list(
                                refFlat=get.ucsc.tbl.tpl("refFlat",what),
                                ensemblToGeneName=
                                    get.ucsc.tbl.tpl("ensemblToGeneName",what),
                                ensemblSource=
                                    get.ucsc.tbl.tpl("ensemblSource",what)
                            ))
                        },
                        tair10 = {
                            warnwrap("Arabidopsis thaliana genome is not ",
                                "supported by UCSC Genome Borwser database! ",
                                "Will automatically switch to Ensembl...",
                                now=TRUE)
                            return(FALSE)
                        }
                    )
                }
            )
        }
    )
}

#' Create SQLite UCSC table template defintions
#'
#' Returns an SQLIte table template defintion, according to  UCSC Genome Browser 
#' database table schemas. This functionality is used when the package RMySQL 
#' is not available for some reason, e.g. Windows machines. Internal use only.
#'
#' @param tab name of UCSC database table.
#' @param what \code{"queries"} for SQLite table definitions or \code{"fields"}
#' for table column names.
#' @return An SQLite table definition.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.table.tmpl <- get.ucsc.tbl.tpl("knownCanonical")
#'}
get.ucsc.tbl.tpl <- function(tab,what="queries") {
    if (what=="queries") {
        switch(tab,
            knownCanonical = {
                return(paste(
                    "CREATE TABLE",
                    "`knownCanonical` (",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`chromStart` INTEGER NOT NULL DEFAULT '0',",
                    "`chromEnd` INTEGER NOT NULL DEFAULT '0',",
                    "`clusterId` INTEGER NOT NULL DEFAULT '0',",
                    "`transcript` TEXT NOT NULL DEFAULT '',",
                    "`protein` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            knownGene = {
                return(paste(
                    "CREATE TABLE",
                    "`knownGene` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`strand` TEXT NOT NULL DEFAULT '',",
                    "`txStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`txEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonCount` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL,",
                    "`proteinID` TEXT NOT NULL DEFAULT '',",
                    "`alignID` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            knownToRefSeq = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToRefSeq` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            refFlat = {
                return(paste("CREATE TABLE",
                    "`refFlat` (",
                    "`geneName` TEXT NOT NULL,",
                    "`name` TEXT NOT NULL,",
                    "`chrom` TEXT NOT NULL,",
                    "`strand` TEXT NOT NULL,",
                    "`txStart` UNSIGNED INTEGER NOT NULL,",
                    "`txEnd` UNSIGNED INTEGER NOT NULL,",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL,",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL,",
                    "`exonCount` UNSIGNED INTEGER NOT NULL,",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            knownToEnsembl = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToEnsembl` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            ensemblSource = {
                return(paste(
                    "CREATE TABLE",
                    "`ensemblSource` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`source` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            mgcGenes = {
                return(paste(
                    "CREATE TABLE `mgcGenes` (",
                    "`bin` UNSIGNED INTEGER NOT NULL,",
                    "`name` TEXT NOT NULL,",
                    "`chrom` TEXT NOT NULL,",
                    "`strand` TEXT NOT NULL,",
                    "`txStart` UNSIGNED INTEGER NOT NULL,",
                    "`txEnd` UNSIGNED INTEGER NOT NULL,",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL,",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL,",
                    "`exonCount` UNSIGNED INTEGER NOT NULL,",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL,",
                    "`score` INTEGER DEFAULT NULL,",
                    "`name2` TEXT NOT NULL,",
                    "`cdsStartStat` TEXT NOT NULL,",
                    "`cdsEndStat` TEXT NOT NULL,",
                    "`exonFrames` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            ensemblToGeneName = {
                return(paste(
                    "CREATE TABLE",
                    "`knownToGeneName` (",
                    "`name` TEXT NOT NULL,",
                    "`value` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            flyBaseCanonical = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseCanonical` (",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`chromStart` INTEGER NOT NULL DEFAULT '0',",
                    "`chromEnd` INTEGER NOT NULL DEFAULT '0',",
                    "`clusterId` INTEGER unsigned NOT NULL DEFAULT '0',",
                    "`transcript` TEXT NOT NULL DEFAULT '',",
                    "`protein` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            },
            flyBaseGene = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseGene` (",
                    "`bin` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`chrom` TEXT NOT NULL DEFAULT '',",
                    "`strand` TEXT NOT NULL DEFAULT '',",
                    "`txStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`txEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsStart` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`cdsEnd` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonCount` UNSIGNED INTEGER NOT NULL DEFAULT '0',",
                    "`exonStarts` TEXT NOT NULL,",
                    "`exonEnds` TEXT NOT NULL",
                    ")",collapse=" "
                ))
            },
            flyBaseToRefSeq = {
                return(paste(
                    "CREATE TABLE",
                    "`flyBaseToRefSeq` (",
                    "`name` TEXT NOT NULL DEFAULT '',",
                    "`value` TEXT NOT NULL DEFAULT ''",
                    ")",collapse=" "
                ))
            }
        )
    }
    else if (what=="fields") {
        switch(tab,
            knownCanonical = {
                return(c("chrom","chromStart","chromEnd","clusterId",
                "transcript","protein"))
            },
            knownGene = {
                return(c("name","chrom","strand","txStart","txEnd","cdsStart",
                    "cdsEnd","exonCount","exonStarts","exonEnds","proteinID",
                    "alignID"))
            },
            knownToRefSeq = {
                return(c("name","value"))
            },
            refFlat = {
                return(c("geneName","name","chrom","strand","txStart","txEnd",
                    "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds"))
            },
            knownToEnsembl = {
                return(c("name","value"))
            },
            ensemblSource = {
                return(c("name","source"))
            },
            mgcGenes = {
                return(c("name","chrom","strand","txStart","txEnd","cdsStart",
                    "cdsEnd","exonCount","exonStarts","exonEnds","score",
                    "name2","cdsStartStat","cdsEndStat","exonFrames"
                ))
            },
            ensemblToGeneName = {
                return(c("name","value"))
            },
            flyBaseCanonical = {
                return(c("chrom","chromStart","chromEnd","clusterId",
                    "transcript","protein"))
            },
            flyBaseGene = {
                return(c("bin","name","chrom","strand","txStart","txEnd",
                    "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds"))
            },
            flyBaseToRefSeq = {
                return(c("name","value"))
            }
        )
    }
}

#' Return queries for the UCSC Genome Browser database, according to organism and
#' source
#'
#' Returns an SQL query to be used with a connection to the UCSC Genome Browser
#' database and fetch metaseqR supported organism annotations. This query is
#' constructed based on the data source and data type to be returned.
#'
#' @param org one of metaseqR supported organisms.
#' @param type either \code{"gene"} or \code{"exon"}.
#' @param refdb one of \code{"ucsc"} or \code{"refseq"} to use the UCSC or RefSeq
#' annotation sources respectively.
#' @return A valid SQL query.
#' @export
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.query <- get.ucsc.query("hg18","gene","ucsc")
#'}
get.ucsc.query <- function(org,type,refdb="ucsc") {
    type <- tolower(type[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    check.text.args("type",type,c("gene","exon"))
    check.text.args("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "dm6","danrer7","pantro4","susscr3","tair10"),multiarg=FALSE)
    check.text.args("refdb",refdb,c("ucsc","refseq"))
    switch(type,
        gene = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,knownGene.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `knownCanonical` INNER ",
                                "JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",                                
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "knownGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownCanonical` INNER JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,knownGene.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `knownCanonical` INNER ",
                                "JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",                                
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                            # Should be the same as hg19 but is like hg18
                        },
                        mm9 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "knownGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownCanonical` INNER JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT knownCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "knownGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownCanonical` INNER JOIN `knownGene` ON ",
                                "knownCanonical.transcript=knownGene.name ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "`gene_id` ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`txStart` AS `start`,`txEnd` ",
                                "AS `end`,mgcGenes.name AS `gene_id`,0 AS ",
                                "`gc_content`,mgcGenes.strand AS `strand`,",
                                "`name2` AS `gene_name`,`source` AS `biotype` ",
                                "FROM `mgcGenes` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT flyBaseCanonical.chrom AS ",
                                "`chromosome`,`chromStart` AS `start`,",
                                "`chromEnd` AS `end`,`transcript` AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "flyBaseGene.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`flyBaseCanonical` INNER JOIN `flyBaseGene` ",
                                "ON flyBaseCanonical.transcript=",
                                "flyBaseGene.name INNER JOIN ",
                                "`flyBaseToRefSeq` ON ",
                                "flyBaseCanonical.transcript=",
                                "flyBaseToRefSeq.name INNER JOIN `refFlat` ON ",
                                "flyBaseToRefSeq.value=refFlat.name INNER ",
                                "JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm6 = {
                            # Stub until we find out what is going on with
                            # Augustus
                        },
                        danrer7 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`txStart` AS `start`,`txEnd` ",
                                "AS `end`,mgcGenes.name AS `gene_id`,0 AS ",
                                "`gc_content`,mgcGenes.strand AS `strand`,",
                                "`name2` AS `gene_name`,`source` AS `biotype` ",
                                "FROM `mgcGenes` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            warnwrap("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                now=TRUE)
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            warnwrap("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                now=TRUE)
                            return(paste(
                                "SELECT refFlat.chrom AS `chromosome`,",
                                "refFlat.txStart AS `start`, refFlat.txEnd AS ",
                                "`end`, refFlat.name AS `gene_id`, 0 AS ",
                                "`gc_content`, refFlat.strand AS `strand`,",
                                "`geneName` AS `gene_name`, `source` AS ",
                                "`biotype` FROM `refFlat` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`,",
                                "`start`",
                                sep=""
                            ))
                        },
                        tair10 = {
                            warnwrap("Arabidopsis thaliana genome is not ",
                                "supported by UCSC Genome Browser database! ",
                                "Will automatically switch to Ensembl...",
                                now=TRUE)
                            return(FALSE)
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,refFlat.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `refFlat` INNER JOIN ",
                                "`knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER ",
                                "JOIN `knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,refFlat.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,'NA' ",
                                "AS `biotype` FROM `refFlat` INNER JOIN ",
                                "`knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                            # Should be the same as hg19 but is as hg18
                        },
                        mm9 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER ",
                                "JOIN `knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT  refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER ",
                                "JOIN `knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY refFlat.name ORDER BY `chromosome`,",
                                " `start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,refFlat.strand ",
                                "AS `strand`,`geneName` AS `gene_name`,",
                                "`source` AS `biotype` FROM `refFlat` INNER ",
                                "JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm6 = {
                            # Stub until we find out what is going on with
                            # Augustus
                        },
                        danrer7 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.txStart AS `start`,",
                                "refFlat.txEnd AS `end`,refFlat.name AS ",
                                "`gene_id`,0 AS `gc_content`,",
                                "refFlat.strand AS `strand`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            return(paste(
                                "SELECT refFlat.chrom AS `chromosome`,",
                                "refFlat.txStart AS `start`, refFlat.txEnd AS ",
                                "`end`, refFlat.name AS `gene_id`, 0 AS ",
                                "`gc_content`, refFlat.strand AS `strand`,",
                                "`geneName` AS `gene_name`, `source` AS ",
                                "`biotype` FROM `refFlat` INNER JOIN ",
                                "`ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`,",
                                "`start`",
                                sep=""
                            ))
                        },
                        tair10 = {
                            warnwrap("Arabidopsis thaliana genome is not ",
                                "supported by UCSC Genome Browser database! ",
                                "Will automatically switch to Ensembl...",
                                now=TRUE)
                            return(FALSE)
                        }
                    )
                }
            )
        },
        exon = {
            switch(refdb,
                ucsc = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                            # Should be the same as hg19 but is as hg18
                        },
                        mm9 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT knownGene.chrom AS ",
                                "`chromosome`,knownGene.exonStarts AS `start`,",
                                "knownGene.exonEnds AS `end`,knownGene.name ",
                                "AS `exon_id`,knownGene.strand AS `strand`,",
                                "`transcript` AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`knownGene` INNER JOIN `knownCanonical` ON ",
                                "knownGene.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "knownCanonical.transcript=knownToRefSeq.name ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "INNER JOIN `refFlat` ON ",
                                "knownToRefSeq.value=refFlat.name GROUP BY ",
                                "knownGene.name ORDER BY `chromosome`, `start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`exonStarts` AS `start`,",
                                "`exonEnds` AS `end`,mgcGenes.name AS ",
                                "`exon_id`,mgcGenes.strand AS `strand`,",
                                "mgcGenes.name AS `gene_id`,`name2` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`mgcGenes` INNER JOIN `ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT flyBaseCanonical.chrom AS ",
                                "`chromosome`,flyBaseGene.exonStarts AS ",
                                "`start`,flyBaseGene.exonEnds AS `end`,",
                                "`transcript` AS `exon_id`,flyBaseGene.strand ",
                                "AS `strand`,`transcript` AS `gene_id`,",
                                "`geneName` AS `gene_name`,`source` AS ",
                                "`biotype` FROM `flyBaseCanonical` INNER JOIN ",
                                "`flyBaseGene` ON ",
                                "flyBaseCanonical.transcript=flyBaseGene.name ",
                                "INNER JOIN `flyBaseToRefSeq` ON ",
                                "flyBaseCanonical.transcript=",
                                "flyBaseToRefSeq.name INNER JOIN `refFlat` ON ",
                                "flyBaseToRefSeq.value=refFlat.name ",
                                "INNER JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm6 = {
                            # Stub until we find out what is going on with
                            # Augustus
                        },
                        danrer7 = {
                            return(paste("SELECT mgcGenes.chrom AS ",
                                "`chromosome`,`exonStarts` AS `start`,",
                                "`exonEnds` AS `end`,mgcGenes.name AS ",
                                "`exon_id`,mgcGenes.strand AS `strand`,",
                                "mgcGenes.name AS `gene_id`,`name2` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`mgcGenes` INNER JOIN `ensemblToGeneName` ON ",
                                "mgcGenes.name2=ensemblToGeneName.value INNER ",
                                "JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            warnwrap("No UCSC Genome annotation for Pan ",
                                "troglodytes! Will use RefSeq instead...",
                                now=TRUE)
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            warnwrap("No UCSC Genome annotation for Sus ",
                                "scrofa! Will use RefSeq instead...",
                                now=TRUE)
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        tair10 = {
                            warnwrap("Arabidopsis thaliana genome is not ",
                                "supported by UCSC Genome Browser database! ",
                                "Will automatically switch to Ensembl...",
                                now=TRUE)
                            return(FALSE)
                        }
                    )
                },
                refseq = {
                    switch(org,
                        hg18 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM `refFlat` ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        hg19 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        hg38 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,'NA' AS `biotype` FROM `refFlat` ",
                                "INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                            # Should be the same as hg19 but is as hg18
                        },
                        mm9 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        mm10 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `knownToRefSeq` ON ",
                                "refFlat.name=knownToRefSeq.value INNER JOIN ",
                                "`knownCanonical` ON ",
                                "knownToRefSeq.name=knownCanonical.transcript ",
                                "INNER JOIN `knownToEnsembl` ON ",
                                "knownCanonical.transcript=knownToEnsembl.name",
                                " INNER JOIN `ensemblSource` ON ",
                                "knownToEnsembl.value=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        rn5 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm3 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "ensemblToGeneName.value=refFlat.geneName ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `gene_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        dm6 = {
                            # Stub until we find out what is going on with
                            # Augustus
                        },
                        danrer7 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        pantro4 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        susscr3 = {
                            return(paste("SELECT refFlat.chrom AS ",
                                "`chromosome`,refFlat.exonStarts AS `start`,",
                                "refFlat.exonEnds AS `end`,refFlat.name AS ",
                                "`exon_id`,refFlat.strand AS `strand`,",
                                "refFlat.name AS `gene_id`,`geneName` AS ",
                                "`gene_name`,`source` AS `biotype` FROM ",
                                "`refFlat` INNER JOIN `ensemblToGeneName` ON ",
                                "refFlat.geneName=ensemblToGeneName.value ",
                                "INNER JOIN `ensemblSource` ON ",
                                "ensemblToGeneName.name=ensemblSource.name ",
                                "GROUP BY `exon_id` ORDER BY `chromosome`, ",
                                "`start`",
                                sep=""))
                        },
                        tair10 = {
                            warnwrap("Arabidopsis thaliana genome is not ",
                                "supported by UCSC Genome Browser database! ",
                                "Will automatically switch to Ensembl...",
                                now=TRUE)
                            return(FALSE)
                        }
                    )
                }
            )
        }
    )
}

#' Return host, username and password for UCSC Genome Browser database
#'
#' Returns a character vector with a hostname, username and password to connect
#' to the UCSC Genome Browser database to retrieve annotation. Internal use.
#'
#' @return A named character vector.
#' @author Panagiotis Moulos
#' @examples
#' \dontrun{
#' db.creds <- get.ucsc.credentials()
#'}
get.ucsc.credentials <- function() {
    return(c(
        host="genome-mysql.cse.ucsc.edu",
        user="genome",
        password=""
    ))
}
