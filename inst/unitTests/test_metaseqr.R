test_metaseqr <- function() {
    data("mm9.gene.data",package="metaseqR")
    ex.dir <- tempdir()

    result.1 <- metaseqr(
        counts=mm9.gene.counts,
        sample.list=sample.list.mm9,
        contrast=c("e14.5_vs_adult_8_weeks"),
        libsize.list=libsize.list.mm9,
        annotation="download",
        org="mm9",
        count.type="gene",
        normalization="edger",
        statistics=c("edger","limma"),
        meta.p="simes",
        preset="medium.basic",
        qc.plots="mds",
        fig.format="png",
        export.where=ex.dir,
        out.list=TRUE
    )
    checkTrue(file.exists(file.path(ex.dir,"index.html")))
    checkTrue(file.exists(file.path(ex.dir,"plots","qc","mds.png")))
    checkTrue(file.exists(file.path(ex.dir,"lists")))
    checkTrue(nrow(result.1[[1]][[1]])>0)
    checkEqualsNumeric(ncol(result.1[[1]][[1]]),16)
}
