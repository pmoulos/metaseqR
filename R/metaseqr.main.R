#' The main metaseqr pipeline
#'
#' This function is the main metaseqr workhorse and implements the main metaseqr 
#' workflow which performs data read, filtering, normalization and statistical 
#' selection, creates diagnostic plots and exports the results and a report if 
#' requested. The metaseqr function is responsible for assembling all the steps 
#' of the metaseqr pipeline which i) reads the input gene or exon read count table 
#' ii) performs prelimininary filtering of data by removing chrM and other 
#' non-essential information for a typical differential gene expression analysis 
#' as well as a preliminary expression filtering based on the exon counts, if an 
#' exon read count file is provided. iii) performs data normalization with one of 
#' currently widely used algorithms, including EDASeq (Risso et al., 2011), DESeq re
#' (Anders and Huber, 2010), edgeR (Robinson et al., 2010), NOISeq (Tarazona et 
#' al., 2012) or no normalization iv) performs a second stage of filtering based 
#' on the normalized gene expression according to several gene filters v) performs
#' statistical testing with one or more of currently widely used algorithms, 
#' including DESeq (Anders and Huber, 2010), edgeR (Robinson et al., 2010), NOISeq 
#' (Tarazona et al., 2012), limma (Smyth et al., 2005) for RNA-Seq data, baySeq 
#' (Hardcastle et al., 2012) vi) in the case of multiple statistical testing 
#' algorithms, performs meta-analysis using one of five available methods (see the
#' meta.p argument) vii) exports the resulting differentially expressed gene list 
#' in text tab-delimited format viii) creates a set of diagnostic plots either 
#' available in the aforementioned packages or metaseqr specific ones and ix) 
#' creates a comprehensive HTML report which summarizes the run information, the 
#' results and the diagnostic plots. Certain diagnostic plots (e.g. the volcano 
#' plot) can be interactive with the use of the external Highcharts 
#' (http://www.highcharts.com) JavaScript library for interactive graphs. Although 
#' the inputs to the metaseqr workflow are many, in practice, setting only very 
#' few of them and accepting the defaults as the rest can result in quite 
#' comprehensible results for mainstream organisms like mouse, human, fruitfly and 
#' rat.
#'
#' @aliases metaseqr.main
#' @param counts a text tab-delimited file containing gene or exon counts in one 
#' of the following formats: i) the first column contains unique gene or exon 
#' identifiers and the rest of the columns contain the read counts for each sample. 
#' Thus the first cell of each row is a gene or exon accession and the rest are 
#' integers representing the counts for that accession. In that case, the 
#' \code{annotation} parameter should strictly be \code{"download"} or an external 
#' file in proper format. ii) The first n columns should contain gene or exon 
#' annotation elements like chromosomal locations, gene accessions, exon accessions, 
#' GC content etc. In that case, the \code{annotation} parameter can also be 
#' \code{"embedded"}. The ideal embedded annotation contains 8 columns, chromosome, 
#' gene or exon start, gene or exon end, gene or exon accession, GC-content 
#' (fraction or percentage), strand, HUGO gene symbol and gene biotype (e.g.
#' "protein_coding" or "ncRNA"). When the \code{annotation} parameter is "embedded", 
#' certain of these features are mandatory (co-ordinates and accessions). If they 
#' are not present, the pipeline will not run. If additional elements are not 
#' present (e.g. GC content or biotypes), certain features of metaseqr will not 
#' be available. For example, EDASeq normalization will not be performed based on
#' a GC content covariate but based on gene length which is not what the authors 
#' of EDASeq suggest. If biotypes are not present, a lot of diagnostic plots will 
#' not be available. If the HUGO gene symbols are missing, the final annotation 
#' will contain only gene accessions and thus be less comprehensible. Generally, 
#' it's best to set the \code{annotation} parameter to \code{"download"} to ensure 
#' the most comprehensible results. Counts can be a data frame satisfying the 
#' above conditions. It is a data frame by default when \code{read2count} is used.
#' counts can also be an .RData file (output of \code{\link{save}} function
#' which contains static input elements (list containing the gene model (exon 
#' counts for each gene constructed by the \code{\link{construct.gene.model}} 
#' function, gene and exon annotation to avoid re-downloading and/or gene counts
#' depending on \code{count.type}). This kind of input facilitates the 
#' re-analysis of the same experiment, using different filtering, normalization 
#' and statistical algorithms. Finally, counts can be a list representing the 
#' gene model (exon counts for each gene) constructed by the 
#' \code{\link{construct.gene.model}} function (provided for backwards
#' compatibility). This .RData file can be generated by setting 
#' \code{save.gene.model=TRUE} when performing data analysis for the first time.
#' @param sample.list a list containing condition names and the samples under each 
#' condition. It should have the format \code{sample.list <-}
#' \code{list(ConditionA=c("Sample_A1",} \code{"Sample_A2", "Sample_A3"),} 
#' \code{ConditionB=c("Sample_B1", "Sample_B2"),} 
#' \code{ConditionC=c("Sample_C1", "Sample_C2"))}. The names of the samples in list 
#' members MUST match the column names containing the read counts in the counts 
#' file. If they do not match, the pipeline will either crash or at best, ignore 
#' several of your samples. Alternative, \code{sample.list} can be a small
#' tab-delimited file structured as follows: the first line of the external tab 
#' delimited file should contain column names (names are not important). The first 
#' column MUST contain UNIQUE sample names and the second column MUST contain the 
#' biological condition where each of the samples in the first column should belong 
#' to. In this case, the function \code{\link{make.sample.list}} is used. If the
#' \code{counts} argument is missing, the \code{sample.list} argument MUST be a 
#' targets text tab-delimited file which contains the sample names, the BAM/BED 
#' file names and the biological conditions/groups for each sample/file. The file 
#' should be text tab-delimited and structured as follows: the first line of the 
#' external tab delimited file should contain column names (names are not important).
#' The first column MUST contain UNIQUE sample names. The second column MUST contain
#' the raw BAM/BED files WITH their full path. Alternatively, the \code{path} 
#' argument should be provided (see below). The third column MUST contain the 
#' biological condition where each of the samples in the first column should belong 
#' to.
#' @param exclude.list a list of samples to exclude, in the same (list) format 
#' as \code{sample.list} above.
#' @param path an optional path where all the BED/BAM files are placed, to be 
#' prepended to the BAM/BED file names in the targets file. If not given and if 
#' the files in the second column of the targets file do not contain a path to a 
#' directory, the current directory is assumed to be the BAM/BED file container.
#' @param file.type the type of raw input files. It can be \code{"auto"} for 
#' auto-guessing, \code{"bed"} for BED files, \code{"sam"} for SAM files or 
#' \code{"bam"} for BAM files.
#' @param contrast a character vector of contrasts to be tested in the statistical 
#' testing step(s) of the metaseqr pipeline. Each element of the should STRICTLY 
#' have the format "ConditionA_vs_ConditionB_vs_...". A valid example based on the
#' \code{sample.list} above is \code{contrast <- c("ConditionA_vs_ConditionB",} 
#' \code{"ConditionA_vs_ConditionC",} \code{"ConditionA_vs_ConditionB_vs_ConditionC")}.
#' The first  element of pairwise contrasts (e.g. "ConditionA" above) MUST be the 
#' control condition or any reference that ConditionB is checked  against. metaseqr
#' uses this convention to properly calculate fold changes. If it's NULL, a contrast
#' between the first two  members of the \code{sample.list} will be auto-generated.
#' @param libsize.list an optional named list where names represent samples (MUST
#' be the same as the samples in \code{sample.list}) and members are the library
#' sizes (the sequencing depth) for each sample. For example
#' \code{libsize.list <- list(Sample_A1=32456913,} \code{Sample_A2=4346818)}.
#' @param id.col an integer denoting the column number in the file (or data frame)
#' provided with the counts argument, where the unique gene or exon accessions are.
#' Default to \code{4} which is the standard feature name column in a BED file.
#' @param gc.col an integer denoting the column number in the file (or data frame)
#' provided with the \code{counts} argument, where each gene's GC content is given.
#' If not provided, GC content normalization provided by EDASeq will not be available.
#' @param name.col an integer denoting the column number in the file (or data frame)
#' provided with the counts argument, where the HUGO gene symbols are given. If
#' not provided, it will not be available when reporting results. In addition, the
#' \code{"known"} gene filter will not be available.
#' @param bt.col an integer denoting the column number in the file (or data frame)
#' provided with the counts argument, where the gene biotypes are given. If not
#' provided, the \code{"biodetection"}, \code{"countsbio"}, \code{"saturation"},
#' \code{"filtered"} and \code{"biodist"} plots will not be available.
#' @param annotation instructs metaseqr where to find the annotation for the given
#' counts file. It can be one of i) \code{"download"} (default) for automatic
#' downloading of the annotation for the organism specified by the org parameter
#' (using biomaRt), ii) \code{"embedded"} if the annotation elements are embedded
#' in the read counts file or iv) a file specified by the user which should be as
#' similar as possible to the \code{"download"} case, in terms of column structure.
#' @param org the supported organisms by metaseqr. These can be, for human genomes
#' \code{"hg18"}, \code{"hg19"} or \code{"hg38"} for mouse genomes \code{"mm9"},
#' \code{"mm10"}, for rat genomes \code{"rn5"}, for drosophila genome \code{"dm3"}, 
#' for zebrafish genome \code{"danrer7"}, for chimpanzee genome \code{"pantro4"},
#' for pig genome \code{"susScr3"} and for Arabidopsis thaliana genome \code{"tair10"}.
#' Finally, \code{"custom"} will instruct metaseqR to completely ignore the 
#' \code{org} argument and depend solely on annotation file provided by the user.
#' @param refdb the reference annotation repository from which to retrieve annotation
#' elements to use with metaseqr. It can be one of \code{"ensembl"} (default),
#' \code{"ucsc"} or \code{"refseq"}.
#' @param count.type the type of reads inside the counts file. It can be one of 
#' \code{"gene"} or \code{"exon"}. This is a very important and mandatory parameter
#' as it defines the course of the workflow.
#' @param exon.filters a named list whose names are the names of the supported
#' exon filters and its members the filter parameters. See section "Exon filters"
#' below for details.
#' @param gene.filters a named list whose names are the names of the supported
#' gene filters and its members the filter parameters. See section "Gene filters"
#' below for details.
#' @param when.apply.filter a character string determining when to apply the exon
#' and/or gene filters, relative to normalization. It can be \code{"prenorm"} to
#' apply apply the filters and exclude genes from further processing before
#' normalization, or \code{"postnorm"} to apply the filters after normalization
#' (default). In the case of \code{when.apply.filter="prenorm"}, a first
#' normalization round is applied to a copy of the gene counts matrix in order to
#' derive the proper normalized values that will constitute the several 
#' expression-based filtering cutoffs.
#' @param normalization the normalization algorithm to be applied on the count
#' data. It can be one of \code{"edaseq"} (default) for EDASeq normalization,
#' \code{"deseq"} for the normalization algorithm (individual options specified
#' by the \code{norm.args} argument) in the DESeq package, \code{"edger"} for the
#' normalization algorithms present in the edgeR package (specified by the
#' \code{norm.args} argument), \code{"noiseq"} for the normalization algorithms
#' present in the NOISeq package (specified by the \code{norm.args} argument), 
#' \code{"nbpseq"} for the normalization algorithms present in the NBPSeq package
#' (specified by the \code{norm.args} argument) or  \code{"none"} to not normalize
#' the data (highly unrecommended). It can also be \code{"each"} where in this
#' case, the normalization applied will be specific to each statistical test used
#' (i.e. the normalization method bundled with each package and used in its
#' examples and documentation).
#' @param norm.args a named list whose names are the names of the normalization
#' algorithm parameters and its members parameter values. See section "Normalization
#' parameters" below for details. Leave \code{NULL} for the defaults of 
#' \code{normalization}. If \code{normalization="each"}, it must be a named list
#' of lists, where each sub-list contains normalization parameters specific to
#' each statistical test to be used.
#' @param statistics one or more statistical analyses to be performed by the
#' metaseqr pipeline.It can be one or more of \code{"deseq"} (default) to conduct
#' statistical test(s) implemented in the DESeq package, \code{"edger"} to conduct
#' statistical test(s) implemented in the edgeR package, \code{"limma"} to conduct
#' the RNA-Seq version of statistical test(s) implemented in the limma package, 
#' \code{"noiseq"} to conduct statistical test(s) implemented in the NOISeq package,
#' \code{"bayseq"} to conduct statistical test(s) implemented in the baySeq package
#' and \code{"nbpseq"} to conduct statistical test(s) implemented in the NBPSeq
#' package. In any case individual algorithm parameters are controlled by the
#' contents of the \code{stat.args} list.
#' @param stat.args a named list whose names are the names of the statistical
#' algorithms used in the pipeline. Each member is another named list whose names
#' are the algorithm parameters and its members are the parameter values. See
#' section "Statistics parameters" below for details. Leave \code{NULL} for the
#' defaults of \code{statistics}.
#' @param adjust.method the multiple testing p-value adjustment method. It can be
#' one of \code{\link{p.adjust.methods}} or \code{"qvalue"} from the qvalue
#' Bioconductor package. Defaults to \code{"BH"} for Benjamini-Hochberg correction.
#' @param meta.p the meta-analysis method to combine p-values from multiple
#' statistical tests \strong{(experimental! see also the second note below,
#' regarding meta-analysis)}. It can be one of \code{"simes"} (default), 
#' \code{"bonferroni"}, \code{"minp"}, \code{"maxp"}, \code{"weight"}, \code{"pandora"},
#' \code{"dperm.min"}, \code{"dperm.max"}, \code{"dperm.weight"}, \code{"fisher"},
#' \code{"fperm"}, \code{"whitlock"} or\code{"none"}. For the \code{"fisher"} and
#' \code{"fperm"} methods, see the documentation of the R package MADAM. For the
#' \code{"whitlock"} method, see the documentation of the survcomp Bioconductor
#' package. With the \code{"maxp"} option, the final p-value is the maximum p-value
#' out of those returned by each statistical test. This is equivalent to an
#' "intersection" of the results derived from each algorithm so as to have a final
#' list with the common genes returned by all statistical tests. Similarly, when
#' \code{meta.p="minp"}, is equivalent to a "union" of the results derived from
#' each algorithm so as to have a final list with all the genes returned by all
#' statistical tests. The latter can be used as a very lose statistical threshold
#' to aggregate results from all methods regardless of their False Positive Rate.
#' With the \code{"simes"} option, the method proposed by Simes (Simes, R. J., 1986)
#' is used. With the \code{"dperm.min"}, \code{"dperm.max"}, \code{"dperm.weight"}
#' options, a permutation procedure is initialed, where \code{nperm} permutations
#' are performed across the samples of the normalized counts matrix, producing
#' \code{nperm} permuted instances of the initital dataset. Then, all the chosen
#' statistical tests are re-executed for each permutation. The final p-value is
#' the number of times that the p-value of the permuted datasets is smaller than
#' the original dataset. The p-value of the original dataset is created based on
#' the choice of one of \code{dperm.min}, \code{dperm.max} or \code{dperm.weight}
#' options. In case of \code{dperm.min}, the intial p-value vector is consists of
#' the minimum p-value resulted from the applied statistical tests for each gene.
#' The maximum p-value is used with the \code{dperm.max} option. With the 
#' \code{dperm.weight} option, the \code{weight} weighting vector for each
#' statistical test is used to weight each p-value according to the power of
#' statistical tests (some might work better for a specific dataset). Be careful
#' as the permutation procedure usually requires a lot of time. However, it should
#' be the most accurate. This method will NOT work when there are no replicated
#' samples across biological conditions. In that case, use \code{meta.p="simes"}
#' instead. Finally, there are the \code{"minp"}, \code{"maxp"} and \code{"weight"}
#' options which correspond to the latter three methods but without permutations.
#' Generally, permutations would be accurate to use when the experiment includes
#' >5 samples per condition (or even better 7-10) which is rather rare in RNA-Seq
#' experiments. Finally, \code{"pandora"} is the same as \code{"weight"} and is
#' added to be in accordance with the metaseqR paper.
#' @param weight a vector of weights with the same length as the \code{statistics}
#' vector containing a weight for each statistical test. It should sum to 1. 
#' \strong{Use with caution with the} \code{dperm.weight} \strong{parameter! 
#' Theoretical background is not yet} \strong{solid and only experience shows 
#' improved results!}
#' @param nperm the number of permutations performed to derive the meta p-value
#' when \code{meta.p="fperm"} or \code{meta.p="dperm"}. It defaults to 10000.
#' @param reprod create reproducible permutations when \code{meta.p="dperm.min"},
#' \code{meta.p="dperm.max"} or \code{meta.p="dperm.weight"}. Ideally one would
#' want to create the same set of indices for a given dataset so as to create
#' reproducible p-values. If \code{reprod=TRUE}, a fixed seed is used by
#' \code{meta.perm} for all the datasets analyzed with \code{metaseqr}. If
#' \code{reprod=FALSE}, then the p-values will not be reproducible, although
#' statistical significance is not expected to change for a large number of
#' resamplings. Finally, \code{reprod} can be a numeric vector of seeds with the
#' same length as \code{nperm} so that the user can supply his/her own seeds.
#' @param pcut a p-value cutoff for exporting differentially genes, default is
#' to export all the non-filtered genes.
#' @param log.offset an offset to be added to values during logarithmic
#' transformations in order to avoid Infinity (default is \code{1}).
#' @param preset an analysis strictness preset. \code{preset} can be one of
#' \code{"all.basic"}, \code{"all.normal"}, \code{"all.full"}, \code{"medium.basic"},
#' \code{"medium.normal"}, \code{"medium.full"}, \code{"strict.basic"},
#' \code{"strict.normal"} or \code{"strict.full"}, each of which control the
#' strictness of the analysis and the amount of data to be exported. For an
#' explanation of the presets, see the section "Presets" below.
#' @param qc.plots a set of diagnostic plots to show/create. It can be one or more
#' of \code{"mds"}, \code{"biodetection"}, \code{"rnacomp"}, \code{"countsbio"},
#' \code{"saturation"}, \code{"readnoise"}, \code{"filtered"}, \code{"boxplot"},
#' \code{"gcbias"}, \code{"lengthbias"}, \code{"meandiff"}, \code{"meanvar"},
#' \code{"deheatmap"}, \code{"volcano"}, \code{"biodist"}, \code{"venn"}. The
#' \code{"mds"} stands for Mutlti-Dimensional Scaling and it creates a PCA-like
#' plot but using the MDS dimensionality reduction instead. It has been succesfully
#' used for NGS data (e.g. see the package htSeqTools) and it shows how well
#' samples from the same condition cluster together. For \code{"biodetection"},
#' \code{"countsbio"}, \code{"saturation"}, \code{"rnacomp"}, \code{"readnoise"},
#' \code{"biodist"} see the vignette of NOISeq package. The \code{"saturation"}
#' case has been rewritten in order to display more samples in a more simple 
#' way. See the help page of \code{\link{diagplot.noiseq.saturation}}. In addition,
#' the \code{"readnoise"} plots represent an older version or the RNA composition
#' plot included in older versions of NOISeq. For \code{"gcbias"},
#' \code{"lengthbias"}, \code{"meandiff"}, \code{"meanvar"} see the vignette of
#' EDASeq package. \code{"lenghtbias"} is similar to \code{"gcbias"} but using the
#' gene length instead of the GC content as covariate. The \code{"boxplot"} option
#' draws boxplots of log2 transformed gene counts. The \code{"filtered"} option
#' draws a 4-panel figure with the filtered genes per chromosome and per biotype,
#' as absolute numbers and as fractions of the genome. See also the help page of
#' \code{\link{diagplot.filtered}}. The \code{"deheatmap"} option performs
#' hierarchical clustering and draws a heatmap of differentially expressed genes.
#' In the context of diagnostic plots, it's useful to see if samples from the
#' same groups cluster together after statistical testing. The \code{"volcano"}
#' option draws a volcano plot for each contrast and if a report is requested, an
#' interactive volcano plot is presented in the HTML report. The \code{"venn"}
#' option will draw an up to 5-way Venn diagram depicting the common and specific
#' to each statistical algorithm genes and for each contrast, when meta-analysis
#' is performed. The \code{"correl"} option creates two correlation graphs: the
#' first one is a correlation heatmap (a correlation matrix which depicts all the
#' pairwise correlations between each pair of samples in the counts matrix is
#' drawn as a clustered heatmap) and the second one is a correlogram plot, which
#' summarizes the correlation matrix in the form of ellipses (for an explanation
#' please see the vignette/documentation of the R package corrplot. Set
#' \code{qc.plots=NULL} if you don't want any diagnostic plots created.
#' @param fig.format the format of the output diagnostic plots. It can be one or
#' more of \code{"png"}, \code{"jpg"}, \code{"tiff"}, \code{"bmp"}, \code{"pdf"},
#' \code{"ps"}. The native format \code{"x11"} (for direct display) is not provided
#' as an option as it may not render the proper display of some diagnostic plots
#' in some devices.
#' @param out.list a logical controlling whether to export a list with the results
#' in the current running environment.
#' @param export.where  an output directory for the project results (report, lists,
#' diagnostic plots etc.).
#' @param export.what the content of the final lists. It can be one or more of 
#' \code{"annotation"}, to bind the annoation elements for each gene, \code{"p.value"},
#' to bind the p-values of each method, \code{"adj.p.value"}, to bind the multiple
#' testing adjusted p-values, \code{"meta.p.value"}, to bind the combined p-value
#' from the meta-analysis, \code{"adj.meta.p.value"}, to bind the corrected  combined
#' p-value from the meta-analysis, \code{"fold.change"}, to bind the fold changes
#' of each requested contrast, \code{"stats"}, to bind several statistics calclulated
#' on raw and normalized counts (see the \code{export.stats} argument), \code{"counts"},
#' to bind the raw and normalized counts for each sample.
#' @param export.scale export values from one or more transformations applied to 
#' the data. It can be one or more of \code{"natural"}, \code{"log2"}, \code{"log10"},
#' \code{"vst"} (Variance Stabilizing Transormation, see the documentation of DESeq 
#' package) and \code{"rpgm"} which is ratio of mapped reads per gene model 
#' (either the gene length or the sum of exon lengths, depending on \code{count.type} 
#' argument). Note that this is not RPKM as reads are already normalized for 
#' library size using one of the supported normalization methods. Also, \code{"rpgm"} 
#' might be misleading when \code{normalization} is other than \code{"deseq"}.
#' @param export.values It can be one or more of \code{"raw"} to export raw values
#' (counts etc.) and \code{"normalized"} to export normalized counts.
#' @param export.stats calculate and export several statistics on raw and normalized
#' counts, condition-wise. It can be one or more of \code{"mean"}, \code{"median"},
#' \code{"sd"}, \code{"mad"}, \code{"cv"} for the Coefficient of Variation,
#' \code{"rcv"} for a robust version of CV where the median and the MAD are used
#' instead of the mean and the standard deviation.
#' @param export.counts.table exports also the calculated read counts table when
#' input is read from bam files and exports also the normalized count table in
#' all cases. Defaults to \code{FALSE}.
#' @param restrict.cores in case of parallel execution of several subfunctions,
#' the fraction of the available cores to use. In some cases if all available cores
#' are used (\code{restrict.cores=1} and the system does not have sufficient RAM,
#' the pipeline running machine might significantly slow down.
#' @param report a logical value controlling whether to produce a summary report
#' or not. Defaults to \code{TRUE}.
#' @param report.top a fraction of top statistically significant genes to append
#' to the HTML report. This helps in keeping the size of the report as small as
#' possible, as appending the total gene list might create a huge HTML file. Users
#' can always retrieve the whole gene lists from the report links. Defaults to
#' \code{0.1} (top 10% of statistically significant genes). Set to \code{NULL}
#' to report all the statistically significant genes.
#' @param report.template an HTML template to use for the report. Do not change
#' this unless you know what you are doing.
#' @param save.gene.model in case of exon analysis, a list with exon counts for
#' each gene will be saved to the file \code{export.where/data/gene_model.RData}.
#' This file can be used as input to metaseqR for exon count based analysis, in
#' order to avoid the time consuming step of assembling the counts for each gene
#' from its exons.
#' @param verbose print informative messages during execution? Defaults to
#' \code{TRUE}.
#' @param run.log write a log file of the \code{metaseqr} run using package log4r.
#' Defaults to \code{TRUE}. The filename will be auto-generated.
#' @param ... further arguments that may be passed to plotting functions, related
#' to \code{\link{par}}.
#' @return If \code{out.list} is \code{TRUE}, a named list whose length is the same
#' as the number of requested contrasts. Each list member is named according to
#' the corresponding contrast and contains a data frame of differentially expressed
#' genes for that contrast. The contents of the data frame are defined by the 
#' \code{export.what, export.scale, export.stats, export.values} parameters. If
#' \code{report} is \code{TRUE}, the output list contains two main elements. The 
#' first is described above (the analysis results) and the second contains the same
#' results but in HTML formatted tables.
#' @section Exon filters: The exon filters are a set of filters which are applied
#' after the gene models are assembled from the read counts of individual exons
#' and before the gene expression is summarized from the exons belonging to each
#' gene. These filters can be applied when the input read counts file contains exon
#' reads. It is not applicable when the input file already contains gene counts.
#' Such filters can be for example "accept genes where all the exons contain more
#' than x reads" or "accept genes where there is read presence in at least m/n
#' exons, n being the total exons of the gene". Such filters are NOT meant for
#' detecting differential splicing as also the whole metaseqr pipeline, thus they
#' should not be used in that context. The \code{exon.filters} argument is a named
#' list of filters, where the names are the filter names and the members are the
#' filter parameters (named lists with parameter name, parameter value). See the
#' usage of the \code{metaseqr} function for an example of how these lists are
#' structured. The supported exon filters in the current version are: i) 
#' \code{min.active.exons} which implements a filter for demanding m out of n exons
#' of a gene to have a certain read presence with parameters \code{exons.per.gene},
#' \code{min.exons} and \code{frac}. The filter is described as follows: if a gene
#' has up to \code{exons.per.gene} exons, then read presence is required in at
#' least \code{min.exons} of them, else read presence is required in a \code{frac}
#' fraction of the total exons. With the default values, the filter instructs that
#' if a gene has up to 5 exons, read presence is required in at least 2, else in
#' at least 20% of the exons, in order to be accepted. More filters will be
#' implemented in future versions and users are encouraged to propose exon filter
#' ideas to the author by mail. See \code{metaseqr} usage for the defaults. Set
#' \code{exon.filters=NULL} to not apply any exon filtering.
#' @section Gene filters: The gene filters are a set of filters applied to gene
#' expression as this is manifested through the read presence on each gene and
#' are preferably applied after normalization. These filters can be applied both
#' when the input file or data frame contains exon read counts and gene read
#' counts. Such filter can be for example "accept all genes above a certain count
#' threshold" or "accept all genes with expression above the median of the
#' normalized counts distribution" or "accept all with length above a certain
#' threshold in kb" or "exclude the 'pseudogene' biotype from further analysis".
#' The supported gene filters in the current version, which have the same structure
#' as the exon filters (named list of lists with filter names, parameter names and
#' parameter arguments)  are: i) \code{length} which implements a length filter
#' where genes are accepted for further analysis if they are above \code{length}
#' (its parameter) kb. ii) \code{avg.reads} which implements a filter where a gene
#' is accepted for further analysis if it has more average reads than the
#' \code{quantile} of the average count distribution per \code{average.per.bp} base
#' pairs. In summary, the reads of each gene are averaged per \code{average.per.bp}
#' based on each gene's length (in case of exons, input the "gene's length" is the
#' sum of the lengths of exons) and the \code{quantile} quantile of the average
#' counts distribution is calculated for each sample. Genes passing the filter
#' should have an average read count larger than the maximum of the vector of the
#' quantiles calculated above. iii) \code{expression} which implements a filter
#' based on the overall expression of a gene. The parameters of this filter are:
#' \code{median}, where genes below the median of the overall count distribution
#' are not accepted for further analysis (this filter has been used to distinguish
#' between "expressed" and "not expressed" genes in several cases, e.g. (Mokry et
#' al., NAR, 2011) with a logical as value, \code{mean} which is the same as
#' \code{median} but using the mean, \code{quantile} which is the same as the
#' previous two but using a specific quantile of the total counts distribution,
#' \code{known}, where in this case, a set of known not-expressed genes in the
#' system under investigation are used to estimate an expression cutoff. This can
#' be quite useful, as the genes are filtered based on a "true biological" cutoff
#' instead of a statistical cutoff. The value of this filter is a character vector
#' of HUGO gene symbols (MUST be contained in the annotation, thus it's better to
#' use \code{annotation="download"}) whose counts are used to build a "null" 
#' expression distribution. The 90th quantile of this distribution is then the
#' expression cutoff. This filter can be combined with any other filter. Be careful
#' with gene names as they are case sensitive and must match exactly ("Pten" is
#' different from "PTEN"!). iv) \code{biotype} where in this case, genes with a
#' certain biotype (MUST be contained in the annotation, thus it's better to use
#' \code{annotation="download"}) are excluded from the analysis. This filter is
#' a named list of logical, where names are the biotypes in each genome and values
#' are \code{TRUE} or \code{FALSE}. If the biotype should be excluded, the value
#' should be \code{TRUE} else \code{FALSE}. See the result of 
#' \code{get.defaults("biotype.filter","hg19")} for an example. Finally, in future
#' versions there will be support for user-defined filters in the form of a function.
#' @section Normalization parameters: The normalization parameters are passed again
#' as a named list where the names of the members are the normalization parameter
#' names and the values are the normalization parameter values. You should check
#' the documentation of the packages EDASeq, DESeq, edgeR, NOISeq and NBPSeq for
#' the parameter names and parameter values. There are a few exceptions in 
#' parameter names: in case of \code{normalization="edaseq"} the only parameter
#' names are \code{within.which} and \code{between.which}, controlling the within
#' lane/sample and between lanes/samples normalization algorithm. In the case of
#' \code{normalization="nbpseq"}, there is one additional parameter called 
#' \code{main.method} which can take the calues \code{"nbpseq"} or \code{"nbsmyth"}.
#' These values correspond to the two different workflows available in the NBPSeq
#' package. Please, consult the NBPSeq package documentation for further details.
#' For the rest of the algorithms, the parameter names are the same as the names
#' used in the respective packages. For examples, please use the 
#' \code{\link{get.defaults}} function.
#' @section Statistics parameters: The statistics parameters as passed to statistical
#' algorithms in metaseqr, exactly with the same way as the normalization parameters
#' above. In this case, there is one more layer in list nesting. Thus, \code{stat.args}
#' is a named list whose names are the names the algorithms used (see the 
#' \code{statistics} parameter). Each member is another named list,with parameters
#' to be used for each statistical algorithm. Again, the names of the member lists
#' are parameter names and the values of the member lists are parameter values.
#' You should check the documentations of DESeq, edgeR, NOISeq, baySeq, limma and
#' NBPSeq for these parameters. There are a few exceptions in parameter names:
#' In case of \code{statistics="edger"}, apart from the rest of the edgeR statistical
#' testing arguments, there is the argument \code{main.method} which can be either
#' \code{"classic"} or \code{"glm"}, defining whether the binomial test or GLMs
#' will be used for statistical testing. For examples, please use the'
#' \code{\link{get.defaults}} function. When \code{statistics="nbpseq"}, apart
#' from the rest arguments of the NBPSeq functions \code{estimate.disp} and
#' \code{estimate.dispersion}, there is the argument \code{main.method} which can
#' be \code{"nbpseq"} or \code{"nbsmyth"}. This argument determines the parameters
#' to be used by the \code{estimate.dispersion} function or by the
#' \code{estimate.disp} function to estimate RNA-Seq count dispersions. The
#' difference between the two is that they constitute different starting points
#' for the two workflows in the package NBPSeq. The first worklfow (with 
#' \code{main.method="nbpseq"} and the \code{estimate.dispersion} function is
#' NBPSeq package specific, while the second (with \code{main.method="nbsmyth"}
#' and the \code{estimate.disp} function is similar to the workflow of the edgeR
#' package. For additional information regarding the statistical testing in
#' NBPSeq, please consult the documentation of the NBPSeq package. 
#' \strong{Additinally, please note that there is currently a problem with the
#' NBPSeq package and the workflow that is specific to the NBPSeq package. The
#' problem has to do with function exporting as there are certain functions which
#' are not recognized from the package internally. For this reason and until it is
#' fixed, only the Smyth workflow will be available with the NBPSeq package (thus}
#' \code{stat.args$main.method="nbpseq"} \strong{will not be available)!}
#' @section Presets: The analysis presets are a set of keywords (only one can be
#' used) that predefine some of the parameters of the metaseqr pipeline. For the
#' time being they are quite simple and they control i) the strictness of
#' filtering and statistical thresholding with three basic levels ("all", "medium",
#' "strict") and ii) the data columns that are exported, again in three basic ways
#' ("basic", "normal", "full") controlling the amount of data to be exported. These
#' keywords can be combined with a dot in the middle (e.g. \code{"all.basic"} to
#' define an analysis preset. When using analysis presets, the following arguments
#' of metaseqr are overriden: \code{exon.filters}, \code{gene.filters}, \code{pcut},
#' \code{export.what}, \code{export.scale}, \code{export.values}, \code{exon.stats}.
#' If you want to explicitly control the above arguments, the \code{preset} argument
#' should be set to \code{NULL} (default). Following is a synopsis of the different
#' presets and the values of the arguments they moderate:
#' \itemize{
#'  \item \code{"all.basic"}: use all genes (do not filter) and export all genes
#'   and basic annotation and statistics elements. In this case, the above described
#'   arguments become:
#'  \itemize{
#'   \item \code{exon.filters=NULL}
#'   \item \code{gene.filters=NULL}
#'   \item \code{pcut=1}
#'   \item \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",}
#'    \code{"adj.meta.p.value","fold.change")}
#'   \item \code{export.scale=c("natural","log2")}
#'   \item \code{export.values=c("normalized")}
#'   \item \code{export.stats=c("mean")}
#'  }
#'  \item \code{"all.normal"}: use all genes (do not filter) and export all genes
#'   and normal annotation and statistics elements. In 
#'   this case, the above described arguments become:
#'  \itemize{
#'   \item \code{exon.filters=NULL}
#'   \item \code{gene.filters=NULL}
#'   \item \code{pcut=1}
#'   \item \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",}
#'    \code{"adj.meta.p.value","fold.change","stats","counts")}
#'   \item \code{export.scale=c("natural","log2")}
#'   \item \code{export.values=c("normalized")}
#'   \item \code{export.stats=c("mean","sd","cv")}
#'  }
#'  \item \code{"all.full"}: use all genes (do not filter) and export all genes and
#'   all available annotation and statistics elements. In this case, the above 
#'   described arguments become:
#'  \itemize{
#'   \item \code{exon.filters=NULL}
#'   \item \code{gene.filters=NULL}
#'   \item \code{pcut=1}
#'   \item \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",}
#'    \code{"adj.meta.p.value","fold.change","stats","counts")}
#'   \item \code{export.scale=c("natural","log2","log10","vst")}
#'   \item \code{export.values=c("raw","normalized")}
#'   \item \code{export.stats=c("mean","median","sd","mad","cv","rcv")}
#'  }
#'  \item \code{"medium.basic"}: apply a medium set of filters and and export
#'   statistically significant genes and basic annotation and statistics elements.
#'   In this case, the above described arguments become:
#'  \itemize{
#'   \item \code{exon.filters=list(min.active.exons=list(exons.per.gene=5,
#'    min.exons=2,frac=1/5))}
#'   \item \code{gene.filters=list(length=list(length=500),}
#'                  \code{avg.reads=list(average.per.bp=100,quantile=0.25),}
#'                    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,
#'                   known=NA,custom=NA),}
#'                    \code{biotype=get.defaults("biotype.filter",org[1]))}
#'   \item \code{pcut=0.05}
#'   \item \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",}
#'    \code{"adj.meta.p.value","fold.change")}
#'   \item \code{export.scale=c("natural","log2")}
#'   \item \code{export.values=c("normalized")}
#'   \item \code{export.stats=c("mean")}
#'  }
#'  \item \code{"medium.normal"}: apply a medium set of filters and and export
#'   statistically significant genes and normal annotation and statistics elements.
#'   In this case, the above described arguments become:
#'  \itemize{
#'   \item \code{exon.filters=list(min.active.exons=list(exons.per.gene=5,
#'    min.exons=2,frac=1/5))}
#'   \item \code{gene.filters=list(length=list(length=500),}
#'                  \code{avg.reads=list(average.per.bp=100,quantile=0.25),}
#'                    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,
#'                   known=NA,custom=NA),}
#'                    \code{biotype=get.defaults("biotype.filter",org[1]))}
#'   \item \code{pcut=0.05}
#'   \item \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",}
#'    \code{"adj.meta.p.value","fold.change","stats","counts")}
#'   \item \code{export.scale=c("natural","log2")}
#'   \item \code{export.values=c("normalized")}
#'   \item \code{export.stats=c("mean","sd","cv")}
#'  }
#  \item \code{"medium.full"}: apply a medium set of filters and and export 
#'  statistically significant genes and all available annotation and statistics 
#'  elements. In this case, the above described arguments become:
#'  \itemize{
#'   \item \code{exon.filters=list(min.active.exons=list(exons.per.gene=5,
#'    min.exons=2,frac=1/5))}
#'   \item \code{gene.filters=list(length=list(length=500),}
#'                  \code{avg.reads=list(average.per.bp=100,quantile=0.25),}
#'                    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,
#'                   known=NA,custom=NA),}
#'                    \code{biotype=get.defaults("biotype.filter",org[1]))}
#'   \item \code{pcut=0.05}
#'   \item \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",}
#'    \code{"adj.meta.p.value","fold.change","stats","counts")}
#'   \item \code{export.scale=c("natural","log2","log10","vst")}
#'   \item \code{export.values=c("raw","normalized")}
#'   \item \code{export.stats=c("mean","median","sd","mad","cv","rcv")}
#'  }
#'  \item \code{"strict.basic"}: apply a strict set of filters and and export
#'   statistically significant genes and basic annotation and statistics elements.
#'   In this case, the above described arguments become:
#'  \itemize{
#'   \item \code{exon.filters=list(min.active.exons=list(exons.per.gene=4,
#'    min.exons=2,frac=1/4))}
#'   \item \code{gene.filters=list(length=list(length=750),}
#'                  \code{avg.reads=list(average.per.bp=100,quantile=0.5),}
#'                    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,
#'                   known=NA,custom=NA),}
#'                    \code{biotype=get.defaults("biotype.filter",org[1]))}
#'   \item \code{pcut=0.01}
#'   \item \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",}
#'    \code{"adj.meta.p.value","fold.change")}
#'   \item \code{export.scale=c("natural","log2")}
#'   \item \code{export.values=c("normalized")}
#'   \item \code{export.stats=c("mean")}
#'  }
#'  \item \code{"strict.normal"}: apply a strict set of filters and and export
#'   statistically significant genes and normal annotation and statistics elements.
#'   In this case, the above described arguments become:
#'  \itemize{
#'   \item \code{exon.filters=list(min.active.exons=list(exons.per.gene=4,
#'    min.exons=2,frac=1/4))}
#'   \item \code{gene.filters=list(length=list(length=750),}
#'                  \code{avg.reads=list(average.per.bp=100,quantile=0.5),}
#'                    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,
#'                   known=NA,custom=NA),}
#'                    \code{biotype=get.defaults("biotype.filter",org[1]))}
#'   \item \code{pcut=0.01}
#'   \item \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",}
#'    \code{"adj.meta.p.value","fold.change","stats","counts")}
#'   \item \code{export.scale=c("natural","log2")}
#'   \item \code{export.values=c("normalized")}
#'   \item \code{export.stats=c("mean","sd","cv")}
#'  }
#  \item \code{"strict.full"}: apply a strict set of filters and and export
#'  statistically significant genes and all available annotation and statistics
#'  elements. In this case, the above described arguments become:
#'  \itemize{
#'   \item \code{exon.filters=list(min.active.exons=list(exons.per.gene=4,
#'    min.exons=2,frac=1/4))}
#'   \item \code{gene.filters=list(length=list(length=750),}
#'                  \code{avg.reads=list(average.per.bp=100,quantile=0.5),}
#'                    \code{expression=list(median=TRUE,mean=FALSE,quantile=NA,
#'                   known=NA,custom=NA),}
#'                    \code{biotype=get.defaults("biotype.filter",org[1]))}
#'   \item \code{pcut=0.01}
#'   \item \code{export.what=c("annotation","p.value","adj.p.value","meta.p.value",}
#'    \code{"adj.meta.p.value","fold.change","stats","counts")}
#'   \item \code{export.scale=c("natural","log2","log10","vst")}
#'   \item \code{export.values=c("raw","normalized")}
#'   \item \code{export.stats=c("mean","median","sd","mad","cv","rcv")}
#'  }
#' }
#' @note Please note that currently only gene and exon annotation from Ensembl
#' (http://www.ensembl.org) are supported. Thus, the unique gene or exon ids in
#' the counts files should correspond to valid Ensembl gene or exon accessions
#' for the organism of interest. If you are not sure about the source of your
#' counts file or do not know how to produce it, it's better to start from the
#' original BAM/BED files (metaseqr will use the \code{\link{read2count}} function
#' to create a counts file). Keep in mind that in the case of BED files, the
#' performance will be significantly lower and the overall running time significantly
#' higher as the R functions which are used to read BED files to proper structures
#' (GenomicRanges) and calculate the counts are quite slow. An alternative way is
#' maybe the easyRNASeq package (Delhomme et al, 2012). The \code{\link{read2count}}
#' function does not use this package but rather makes use of standard Bioconductor
#' functions to handle NGS data. If you wish to work outside R, you can work with
#' other popular read counters such as the HTSeq read counter
#' (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html). Please also
#' note that in the current version, the members of the \code{gene.filters} and
#' \code{exon.filters} lists are not checked for validity so be careful to supply 
#' with correct names otherwise the pipeline will crash or at the best case scenario,
#' will ignore the filters. Also note that when you are supplying metaseqr with
#' an exon counts table, gene annotation is always downloaded so please be sure
#' to have a working internet connection. In addition to the above, if you have
#' a multiple core system, be very careful on how you are using the 
#' \code{restrict.cores} argument and generally how many cores you are using with
#' scripts purely written in R. The analysis with exon read data can very easily
#' cause memory problems, so unless you have more than 64Gb of RAM available,
#' consider setting restrict.cores to something like 0.2 when working with exon
#' data. Finally, if you do not wish to download the same annotation again and
#' again when performing multiple analyses, it is best to use the
#' \code{\link{get.annotation}} function to download and store the resulting
#' data frames in local files and then use these files with the \code{annotation}
#' option.
#' @note Please note that the \strong{meta-analysis} feature provided by metaseqR
#' is currently experimental and does not satisfy the strict definition of 
#' "meta-analysis", which is the combination of multiple similar datasets under 
#' the same statistical methodology. Instead it is the use of mulitple statistical
#' tests applied to the same data so the results at this point are not guaranteed
#' and should be interpreted appropriately. We are working on a more solid
#' methodology for combining multiple statistical tests based on multiple testing
#' correction and Monte Carlo methods. For the Simes method, please consult also
#' "Simes, R. J. (1986). "An improved Bonferroni procedure for multiple tests of
#' significance". Biometrika 73 (3): 751–754."
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' # An example pipeline with exon counts
#' data("hg19.exon.data",package="metaseqR")
#' metaseqr(
#'  counts=hg19.exon.counts,
#'  sample.list=list(normal="normal",paracancerous="paracancerous",
#'    cancerous="cancerous"),
#'  contrast=c("normal_vs_paracancerous","normal_vs_cancerous",
#'    "normal_vs_paracancerous_vs_cancerous"),
#'  libsize.list=libsize.list.hg19,
#'  id.col=4,
#'  annotation="download",
#'  org="hg19",
#'  count.type="exon",
#'  normalization="edaseq",
#'  statistics="deseq",
#'  pcut=0.05,
#'  qc.plots=c("mds", "biodetection", "countsbio", "saturation", "rnacomp",
#'     "boxplot", "gcbias", "lengthbias", "meandiff", "readnoise","meanvar",
#'     "readnoise", "deheatmap", "volcano", "biodist", "filtered"),
#'  fig.format=c("png","pdf"),
#'  export.what=c("annotation","p.value","adj.p.value","fold.change","stats",
#'    "counts"),
#'  export.scale=c("natural","log2","log10","vst"),
#'  export.values=c("raw","normalized"),
#'  export.stats=c("mean","median","sd","mad","cv","rcv"),
#'  restrict.cores=0.8,
#'  gene.filters=list(
#'    length=list(
#'      length=500
#'    ),
#'    avg.reads=list(
#'      average.per.bp=100,
#'      quantile=0.25
#'    ),
#'    expression=list(
#'      median=TRUE,
#'      mean=FALSE
#'    ),
#'    biotype=get.defaults("biotype.filter","hg18")
#'  )
#' )
#'
#' # An example pipeline with gene counts
#' data("mm9.gene.data",package="metaseqR")
#' result <- metaseqr(
#'  counts=mm9.gene.counts,
#'  sample.list=list(e14.5=c("e14.5_1","e14.5_2"), 
#'    adult_8_weeks=c("a8w_1","a8w_2")),
#'  contrast=c("e14.5_vs_adult_8_weeks"),
#'  libsize.list=libsize.list.mm9,
#'  annotation="download",
#'  org="mm9",
#'  count.type="gene",
#'  normalization="edger",
#'  statistics=c("deseq","edger","noiseq"),
#'  meta.p="fisher",
#'  pcut=0.05,
#'  fig.format=c("png","pdf"),
#'  export.what=c("annotation","p.value","meta.p.value","adj.meta.p.value",
#'    "fold.change"),
#'  export.scale=c("natural","log2"),
#'  export.values="normalized",
#'  export.stats=c("mean","sd","cv"),
#'  export.where=getwd(),
#'  restrict.cores=0.8,
#'  gene.filters=list(
#'    length=list(
#'      length=500
#'    ),
#'    avg.reads=list(
#'         average.per.bp=100,
#'         quantile=0.25
#'    ),
#'    expression=list(
#'         median=TRUE,
#'         mean=FALSE,
#'         quantile=NA,
#'         known=NA,
#'         custom=NA
#'    ),
#'    biotype=get.defaults("biotype.filter","mm9")
#'  ),
#'  out.list=TRUE
#' )
#' head(result$data[["e14.5_vs_adult_8_weeks"]])
#' }
metaseqr <- function(
    counts,
    sample.list,
    exclude.list=NULL,
    file.type=c("auto","sam","bam","bed"),
    path=NULL,
    contrast=NULL,
    libsize.list=NULL,
    id.col=4,
    gc.col=NA,
    name.col=NA,
    bt.col=NA,
    annotation=c("download","embedded"),
    org=c("hg18","hg19","hg38","mm9","mm10","rn5","dm3","danrer7","pantro4",
        "susscr3","tair10","custom"),
    refdb=c("ensembl","ucsc","refseq"),
    count.type=c("gene","exon"),
    exon.filters=list(
        min.active.exons=list(
            exons.per.gene=5,
            min.exons=2,
            frac=1/5
        )
    ),
    gene.filters=list(
        length=list(
            length=500
        ),
        avg.reads=list(
            average.per.bp=100,
            quantile=0.25
        ),
        expression=list(
            median=TRUE,
            mean=FALSE,
            quantile=NA,
            known=NA,
            custom=NA
        ),
        biotype=get.defaults("biotype.filter",org[1])
    ),
    when.apply.filter=c("postnorm","prenorm"),
    normalization=c("edaseq","deseq","edger","noiseq","nbpseq","each","none"),
    norm.args=NULL,
    statistics=c("deseq","edger","noiseq","bayseq","limma","nbpseq"),
    stat.args=NULL,
    adjust.method=sort(c(p.adjust.methods,"qvalue")), # Brings BH first which is the default
    meta.p=if (length(statistics)>1) c("simes","bonferroni","fisher",
        "dperm.min","dperm.max","dperm.weight","fperm","whitlock","minp","maxp",
        "weight","pandora","none") else "none",
    weight=rep(1/length(statistics),length(statistics)),
    nperm=10000,
    reprod=TRUE,
    pcut=NA, # A p-value cutoff for exporting DE genes, default is to export all
    log.offset=1, # Logarithmic transformation offset to avoid +/-Inf (log2(a+offset/b+offset))
    preset=NULL, # An analysis strictness preset
    qc.plots=c(
        "mds","biodetection","countsbio","saturation","readnoise","filtered",
        "correl","pairwise", # Raw count data
        "boxplot","gcbias","lengthbias","meandiff","meanvar","rnacomp", # Pre and post normalization
        "deheatmap","volcano","biodist" # Post statistical testing
    ),
    fig.format=c("png","jpg","tiff","bmp","pdf","ps"),
    out.list=FALSE,
    export.where=NA, # An output directory for the project
    export.what=c("annotation","p.value","adj.p.value","meta.p.value",
        "adj.meta.p.value","fold.change","stats","counts","flags"),
    export.scale=c("natural","log2","log10","vst","rpgm"),
    export.values=c("raw","normalized"),
    export.stats=c("mean","median","sd","mad","cv","rcv"),
    export.counts.table=FALSE,
    restrict.cores=0.6,
    report=TRUE,
    report.top=0.1,
    report.template="default",
    save.gene.model=TRUE,
    verbose=TRUE,
    run.log=TRUE,
    ...
)

{
    # Check essential arguments
    from.raw <- from.previous <- FALSE
    if (missing(counts) && (missing(sample.list) || is.list(sample.list)))
        stop("You must provide a file with genomic region (gene, exon, etc.) ",
            "counts or an input targets file to create input from! If the ",
            "counts file is missing, sample.list cannot be missing or it must ",
            "be a targets file with at least three columns! See the ",
            "read.targets function. counts may also be a gene model list ",
            "(see the documentation)")
    if (!missing(counts) && !missing(sample.list) && is.character(counts) 
        && file.exists(counts) && length(grep(".RData$",counts))>0) 
    {
        warning("When restoring a previous analysis, sample.list argument is ",
            "not necessary! Ignoring...")
        from.previous <- TRUE
        tmp.env <- new.env()
        disp("Restoring previous analysis from ",basename(counts))
        load(counts,tmp.env)
        sample.list <- tmp.env$sample.list
        count.type <- tmp.env$count.type
    }
    if (!missing(counts) && missing(sample.list) && is.character(counts) 
        && file.exists(counts) && length(grep(".RData$",counts))>0)
    { # Time to load previous analysis if existing
        from.previous <- TRUE
        tmp.env <- new.env()
        message("Restoring previous analysis from ",basename(counts))
        load(counts,tmp.env)
        sample.list <- tmp.env$sample.list
        count.type <- tmp.env$count.type
    }
    if (missing(sample.list) && !from.previous || (!is.list(sample.list) &&
        !file.exists(sample.list)))
        stop("You must provide a list with condition names and sample names ",
            "(same as in the counts file) or an input file to create the ",
            "sample list from!")
    if (!missing(sample.list) && !is.list(sample.list) 
        && file.exists(sample.list) && !missing(counts) && !from.previous)
        sample.list <- make.sample.list(sample.list)
    if (!missing(sample.list) && !is.list(sample.list) 
        && file.exists(sample.list) && missing(counts))
    {
        counts <- NULL
        the.list <- read.targets(sample.list,path=path)
        sample.list <- the.list$samples
        file.list <- the.list$files
        if (tolower(file.type[1])=="auto")
            file.type <- the.list$type
        if (is.null(file.type))
            stop(paste("The type of the input files could not be recognized!",
                "Please specify (BAM or BED)..."))
        from.raw <- TRUE
    }

    # Initialize environmental variables
    HOME <- system.file(package="metaseqR")
    TEMPLATE <- HOME
    #if (!exists("HOME"))
    #    init.envar()
    # Globalize the project's verbosity and logger
    if (from.raw)
        PROJECT.PATH <- make.project.path(export.where)
    else
        PROJECT.PATH <- make.project.path(export.where,counts)
    assign("VERBOSE",verbose,envir=meta.env)
    if (run.log)
        logger <- create.logger(logfile=file.path(PROJECT.PATH$logs,
            "metaseqr_run.log"),level=2,logformat="%d %c %m")
    else
        logger <- NULL
    assign("LOGGER",logger,envir=meta.env)

    # Check if there are any mispelled or invalid parameters and throw a warning
    check.main.args(as.list(match.call()))
    
    # Check if sample names match in file/df and list, otherwise meaningless to proceed
    if (!from.raw && !from.previous)
    {
        if (!is.data.frame(counts) && !is.list(counts))
        {
            if (file.exists(counts))
            {
                aline <- read.delim(counts,nrows=5) # Read the 1st lines
                aline <- colnames(aline)
            }
            else
                stopwrap("The counts file you provided does not exist!")
        }
        else if (is.data.frame(counts))
            aline <- colnames(counts)
        else if (is.list(counts))
            aline <- names(counts)
        samples <- unlist(sample.list,use.names=FALSE)
        if (length(which(!is.na(match(samples,aline)))) != length(samples))
            stopwrap("The sample names provided in the counts file/list do ",
                "not match with those of the sample.list!")
    }
    
    # If exclude list given, check that it's a subset of sample.list, otherwise
    # just ignore exclude.list
    if (!is.null(exclude.list) && !is.na(exclude.list))
    {
        sl <- unlist(sample.list)
        el <- unlist(exclude.list)
        if (length(intersect(sl,el)) != length(el))
        {
            warnwrap("Some samples in exclude.list do not match those in the ",
                "initial sample.list! Ignoring...",now=TRUE)
            exclude.list <- NULL
        }
    }   

    file.type <- tolower(file.type[1])
    annotation <- tolower(annotation[1])
    org <- tolower(org[1])
    refdb <- tolower(refdb[1])
    count.type <- tolower(count.type[1])
    when.apply.filter <- tolower(when.apply.filter[1])
    normalization <- tolower(normalization[1])
    adjust.method <- adjust.method[1]
    meta.p <- tolower(meta.p[1])
    statistics <- tolower(statistics)
    fig.format <- tolower(fig.format)
    if (!is.null(qc.plots)) qc.plots <- tolower(qc.plots)
    export.what <- tolower(export.what)
    export.scale <- tolower(export.scale)
    export.values <- tolower(export.values)
    export.stats <- tolower(export.stats)
    if (!is.null(preset)) preset <- tolower(preset[1])
    
    if (from.raw)
        counts.name <- "imported sam/bam/bed files"
    else
    {
        if (!is.data.frame(counts) && !is.null(counts) && !is.list(counts))
        {
            check.file.args("counts",counts)
            if (from.previous)
                counts.name <- "previously stored project"
            else
                counts.name <- basename(counts)
        }
        else if (is.list(counts) && !is.data.frame(counts))
        {
            counts.name <- "previously stored gene model"
        }
        else
        {
            counts.name <- "imported custom data frame"
        }
    }

    if (is.list(counts) && count.type=="exon" && annotation=="embedded")
    {
        warnwrap("annotation cannot be \"embedded\" when importing a stored ",
            "gene model! Switching to \"download\"...")
        annotation <- "download"
    }

    if (meta.p %in% c("weight","pandora","dperm.weight") && 
        abs(1-sum(weight))>1e-5)
        stopwrap("The weights given for p-value combination should sum to 1!")

    check.text.args("file.type",file.type,c("auto","sam","bam","bed"),
        multiarg=FALSE)
    check.text.args("annotation",annotation,c("embedded","download"),
        multiarg=FALSE)
    check.text.args("org",org,c("hg18","hg19","hg38","mm9","mm10","rn5","dm3",
        "danrer7","pantro4","susscr3","tair10","custom"),multiarg=FALSE)
    check.text.args("refdb",refdb,c("ensembl","ucsc","refseq"),multiarg=FALSE)
    check.text.args("count.type",count.type,c("gene","exon"),multiarg=FALSE)
    check.text.args("when.apply.filter",when.apply.filter,c("postnorm",
        "prenorm"),multiarg=FALSE)
    check.text.args("normalization",normalization,c("edaseq","deseq","edger",
        "noiseq","nbpseq","each","none"),multiarg=FALSE)
    check.text.args("statistics",statistics,c("deseq","edger","noiseq","bayseq",
        "limma","nbpseq"),multiarg=TRUE)
    check.text.args("meta.p",meta.p,c("simes","bonferroni","fisher","dperm.min",
        "dperm.max","dperm.weight","fperm","whitlock","minp","maxp","weight",
        "pandora","none"),multiarg=FALSE)
    check.text.args("fig.format",fig.format,c("png","jpg","tiff","bmp","pdf",
        "ps"),multiarg=TRUE)
    check.text.args("export.what",export.what,c("annotation","p.value",
        "adj.p.value","meta.p.value","adj.meta.p.value","fold.change","stats",
        "counts","flags"),multiarg=TRUE)
    check.text.args("export.scale",export.scale,c("natural","log2","log10",
        "rpgm","vst"),multiarg=TRUE)
    check.text.args("export.values",export.values,c("raw","normalized"),
        multiarg=TRUE)
    check.text.args("export.stats",export.stats,c("mean","median","sd","mad",
        "cv","rcv"),multiarg=TRUE)
    if (!is.null(preset))
        check.text.args("preset",preset,c("all.basic","all.normal","all.full",
            "medium.basic","medium.normal","medium.full","strict.basic",
            "strict.normal","strict.full"),multiarg=FALSE)
    if (!is.null(qc.plots))
        check.text.args("qc.plots",qc.plots,c("mds","biodetection","countsbio",
            "saturation","readnoise","correl","pairwise","boxplot","gcbias",
            "lengthbias","meandiff","meanvar","rnacomp","deheatmap","volcano",
            "biodist","filtered","venn"),multiarg=TRUE)
    if (!is.na(restrict.cores)) check.num.args("restrict.cores",restrict.cores,
        "numeric",c(0,1),"botheq")
    if (!is.na(pcut)) check.num.args("pcut",pcut,"numeric",c(0,1),"botheq")
    if (!is.na(gc.col)) check.num.args("gc.col",gc.col,"numeric",0,"gt")
    if (!is.na(name.col)) check.num.args("name.col",name.col,"numeric",0,"gt")
    if (!is.na(bt.col)) check.num.args("bt.col",bt.col,"numeric",0,"gt")
    if (!is.na(log.offset)) check.num.args("log.offset",log.offset,"numeric",0,
        "gt")
    check.num.args("nperm",nperm,"numeric",10,"gt")
    if (!is.null(report.top))
        check.num.args("report.top",report.top,"numeric",c(0,1),"both")
    if (!is.null(contrast)) check.contrast.format(contrast,sample.list)
    if ("bayseq" %in% statistics) libsize.list <- check.libsize(libsize.list,
        sample.list)

    # Check main functionality packages
    check.packages(meta.p,qc.plots)
    # Check if parallel processing is available
    multic <- check.parallel(restrict.cores)
    # Check the case of embedded annotation but not given gc and gene name columns
    if (annotation=="embedded")
    {
        if (is.na(gc.col) && count.type=="gene")
            stopwrap("The column that contains the gene GC content ",
                "(\"gc.col\") argument is required when \"annotation\" is ",
                "\"embedded\"!")
        if (is.na(name.col) && !is.na(gene.filters$expression$known))
        {
            warnwrap("The column that contains the HUGO gene symbols ",
                "(\"bt.col\") is missing with embedded annotation! Gene name ",
                "expression filter will not be available...")
            gene.filters$expression$known=NA
            if ("volcano" %in% qc.plots)
                warnwrap("The column that contains the HUGO gene symbols ",
                    "(\"bt.col\") is missing with embedded annotation! ",
                    "Interactive volcano plots will not contain gene names...")
        }
        if (is.na(bt.col) && count.type=="gene")
        {
            warnwrap("The column that contains the gene biotypes (\"bt.col\") ",
                "is missing with embedded annotation! Biotype filters and ",
                "certain plots will not be available...")
            gene.filters$biotype=NULL
            to.remove <- match(c("biodetection","countsbio","saturation",
                "biodist","filtered"),qc.plots)
            no.match <- which(is.na(to.remove))
            if (length(no.match)>0)
                to.remove <- to.remove[-no.match]
            if (length(to.remove)>0)
                qc.plots <- qc.plots[-to.remove]
        }
    }
    if (org=="hg18" && (refdb %in% c("ucsc","refseq")))
    {
        warnwrap("Gene/exon biotypes cannot be retrieved when organism is ",
            "\"hg18\" and annotation database is \"ucsc\" or \"refseq\"! ",
            "Biotype filters and certain plots will not be available...")
        gene.filters$biotype=NULL
        to.remove <- match(c("biodetection","countsbio","saturation",
            "biodist","filtered"),qc.plots)
        no.match <- which(is.na(to.remove))
        if (length(no.match)>0)
            to.remove <- to.remove[-no.match]
        if (length(to.remove)>0)
            qc.plots <- qc.plots[-to.remove]
    }
    
    # Check if drawing a Venn diagram is possible
    if ("venn" %in% qc.plots && length(statistics)==1)
    {
        warnwrap("The creation of a Venn diagram is possible only when more ",
            "than one statistical algorithms are used (meta-analysis)! ",
            "Removing from figures list...")
        to.remove <- match("venn",qc.plots)
        no.match <- which(is.na(to.remove))
        if (length(no.match)>0)
            to.remove <- to.remove[-no.match]
        if (length(to.remove)>0)
            qc.plots <- qc.plots[-to.remove]
    }

    # Check additional input arguments for normalization and statistics
    alg.args <- validate.alg.args(normalization,statistics,norm.args,stat.args)
    norm.args <- alg.args$norm.args
    stat.args <- alg.args$stat.args
    
    # Override settigs if a preset is given
    if (!is.null(preset))
    {
        preset.opts <- get.preset.opts(preset,org)
        exon.filters <- preset.opts$exon.filters
        gene.filters <- preset.opts$gene.filters
        pcut <- preset.opts$pcut
        export.what <- preset.opts$export.what
        export.scale <- preset.opts$export.scale
        export.values <- preset.opts$export.values
        export.stats <- preset.opts$export.stats
    }

    if (report)
    {
        report.messages <- make.report.messages("en")
        if (!is.null(qc.plots) && !("png" %in% fig.format))
        {
            warnwrap("png format is required in order to build a report! ",
                "Adding to figure output formats...")
            fig.format <- c(fig.format,"png")
        }
    }

    # Display initialization report
    TB <- Sys.time()
    disp(strftime(Sys.time()),": Data processing started...\n")
    ############################################################################
    disp("Read counts file: ",counts.name)
    disp("Conditions: ",paste(names(sample.list),collapse=", "))
    disp("Samples to include: ",paste(unlist(sample.list),collapse=", "))
    if (!is.null(exclude.list) && !is.na(exclude.list))
        disp("Samples to exclude: ",paste(unlist(exclude.list),collapse=", "))
    else
        disp("Samples to exclude: none")
    disp("Requested contrasts: ",paste(contrast,collapse=", "))
    if (!is.null(libsize.list))
    {
        disp("Library sizes: ")
        for (n in names(libsize.list))
            disp("  ",paste(n,libsize.list[[n]],sep=": "))
    }
    disp("Annotation: ",annotation)
    disp("Organism: ",org)
    disp("Reference source: ",refdb)
    disp("Count type: ",count.type)
    if (!is.null(preset))
        disp("Analysis preset: ",preset)
    if (!is.null(exon.filters))
    {
        disp("Exon filters: ",paste(names(exon.filters),collapse=", "))
        for (ef in names(exon.filters))
        {
            disp("  ",ef,": ")
            for (efp in names(exon.filters[[ef]]))
            {
                if (length(exon.filters[[ef]][[efp]])==1 && 
                    is.function(exon.filters[[ef]][[efp]]))
                    print(exon.filters[[ef]][[efp]])
                else if (length(exon.filters[[ef]][[efp]])==1)
                    disp("    ",paste(efp,exon.filters[[ef]][[efp]],sep=": "))
                else if (length(exon.filters[[ef]][[efp]])>1)
                    disp("    ",paste(efp,paste(exon.filters[[ef]][[efp]],
                        collapse=", "),sep=": "))
            }
        }
    }
    else
        disp("Exon filters: none applied")
    if (!is.null(gene.filters))
    {
        disp("Gene filters: ",paste(names(gene.filters),collapse=", "))
        for (gf in names(gene.filters))
        {
            disp("  ",gf,": ")
            for (gfp in names(gene.filters[[gf]]))
            {
                if (length(gene.filters[[gf]][[gfp]])==1 && 
                    is.function(gene.filters[[gf]][[gfp]]))
                    print(gene.filters[[gf]][[gfp]])
                else if (length(gene.filters[[gf]][[gfp]])==1)
                    disp("    ",paste(gfp,gene.filters[[gf]][[gfp]],sep=": "))
                else if (length(gene.filters[[gf]][[gfp]])>1)
                    disp("    ",paste(gfp,paste(gene.filters[[gf]][[gfp]],
                        collapse=", "),sep=": "))
            }
        }
    }
    else
        disp("Gene filters: none applied")
    disp("Filter application: ",when.apply.filter)
    disp("Normalization algorithm: ",normalization)
    if (!is.null(norm.args))
    {
        disp("Normalization arguments: ")
        for (na in names(norm.args))
        {
            if (length(norm.args[[na]])==1 && is.function(norm.args[[na]]))
            {
                disp("  ",na,": ")
                disp(as.character(substitute(norm.args[[na]])))
            }
            else if (length(norm.args[[na]])==1)
                disp("  ",paste(na,norm.args[[na]],sep=": "))
            else if (length(norm.args[[na]])>1)
                disp("  ",paste(na,paste(norm.args[[na]],collapse=", "),
                    sep=": "))
        }
    }
    disp("Statistical algorithm: ",paste(statistics,collapse=", "))
    if (!is.null(stat.args))
    {
        disp("Statistical arguments: ")
        for (sa in names(stat.args))
        {
            if (length(stat.args[[sa]])==1 && is.function(stat.args[[sa]]))
            {
                disp("  ",sa,": ")
                disp(as.character(substitute(stat.args[[na]])))
            }
            else if (length(stat.args[[sa]])==1)
                disp("  ",paste(sa,stat.args[[sa]],sep=": "))
            else if (length(stat.args[[sa]])>1)
                disp("  ",paste(sa,paste(stat.args[[sa]],collapse=", "),
                    sep=": "))
        }
    }
    disp("Meta-analysis method: ",meta.p)
    disp("Multiple testing correction: ",adjust.method)
    if (!is.na(pcut)) disp("p-value threshold: ",pcut)
    disp("Logarithmic transformation offset: ",log.offset)
    if (!is.null(preset)) disp("Analysis preset: ",preset)
    disp("Quality control plots: ",paste(qc.plots,collapse=", "))
    disp("Figure format: ",paste(fig.format,collapse=", "))
    if (!is.na(export.where)) disp("Output directory: ",export.where)
    disp("Output data: ",paste(export.what,collapse=", "))
    disp("Output scale(s): ",paste(export.scale,collapse=", "))
    disp("Output values: ",paste(export.values,collapse=", "))
    if ("stats" %in% export.what)
        disp("Output statistics: ",paste(export.stats,collapse=", "),"\n")
    ############################################################################

    if (count.type=="exon")
    {
        # Download gene annotation anyway if not previous analysis restored
        if (!from.previous) 
        {
            disp("Downloading gene annotation for ",org,"...")
            gene.data <- get.annotation(org,"gene",refdb)
        }
        
        if (!from.previous)
        {
            if (annotation=="download")
            {
                disp("Downloading exon annotation for ",org,"...")
                exon.data <- get.annotation(org,count.type,refdb,multic)
            }
            else if (annotation=="embedded")
            {
                # The following should work if annotation elements are arranged in 
                # MeV-like data style
                # Embedded annotation can NEVER occur when receiving data from 
                # read2count, so there is no danger here
                if (!is.data.frame(counts))
                {
                    disp("Reading counts file ",counts.name,"...")
                    exon.counts <- read.delim(counts)
                }
                else
                    exon.counts <- counts
                rownames(exon.counts) <- as.character(exon.counts[,id.col])
                all.cols <- 1:ncol(exon.counts)
                sam.cols <- match(unlist(sample.list),colnames(exon.counts))
                sam.cols <- sam.cols[which(!is.na(sam.cols))]
                ann.cols <- all.cols[-sam.cols]
                exon.data <- exon.counts[,ann.cols]
                exon.counts <- exon.counts[,sam.cols]
                colnames(exon.data)[id.col] <- "exon_id"
                if (!is.na(name.col)) colnames(exon.data)[name.col] <- 
                    "gene_name"
                if (!is.na(bt.col)) colnames(exon.data)[bt.col] <- "biotype"
                exon.counts <- cbind(exon.data[rownames(exon.counts),c("start",
                    "end","exon_id","gene_id")],exon.counts)
            }
            else # Reading from external file, similar to embedded
            {
                disp("Reading external exon annotation for ",org," from ",
                    annotation,"...")
                exon.data <- read.delim(annotation)
                colnames(exon.data)[id.col] <- "exon_id"
            }
        }
        else
        {
            counts <- tmp.env$the.counts
            exon.data <- tmp.env$exon.data
            gene.data <- tmp.env$gene.data
        }

        # Else everything is provided and done
        #if (is.data.frame(counts))
        if (annotation!="embedded" & !from.previous)
        {
            if (!is.null(counts)) # Otherwise it's coming ready from read2count
            {
                if (!is.data.frame(counts) && !is.list(counts))
                {
                    disp("Reading counts file ",counts.name,"...")
                    exon.counts <- read.delim(counts)
                }
                else # Already a data frame as input
                    exon.counts <- counts
                rownames(exon.counts) <- as.character(exon.counts[,id.col])
                exon.counts <- exon.counts[,unlist(sample.list,
                    use.names=FALSE)]
            }
            else # Coming from read2count
            {
                if (from.raw) # Double check
                {
                    r2c <- read2count(the.list,exon.data,file.type,
                        multic=multic)
                    exon.counts <- r2c$counts
                    # Merged exon data!
                    exon.data <- r2c$mergedann
                    if (is.null(libsize.list))
                        libsize.list <- r2c$libsize
                    if (export.counts.table) {
                        disp("Exporting raw read counts table to ",
                            file.path(PROJECT.PATH[["lists"]],
                            "raw_counts_table.txt.gz"))
                        res.file <- file.path(PROJECT.PATH[["lists"]],
                            "raw_counts_table.txt.gz")
                        gzfh <- gzfile(res.file,"w")
                        write.table(cbind(
                            exon.data[rownames(exon.counts),],
                            exon.counts),gzfh,sep="\t",row.names=FALSE,
                            quote=FALSE)
                        close(gzfh)
                    }
                }
            }
            exon.counts <- cbind(exon.data[rownames(exon.counts),c("start","end",
                "exon_id","gene_id")],exon.counts[,unlist(sample.list,
                use.names=FALSE)])

            # Get the exon counts per gene model
            disp("Checking chromosomes in exon counts and gene annotation...")
            gene.data <- reduce.gene.data(exon.data[rownames(exon.counts),],
                gene.data)
            disp("Processing exons...")
            the.counts <- construct.gene.model(exon.counts,sample.list,gene.data,
                multic=multic)

            if (save.gene.model)
            {
                disp("Saving gene model to ",file.path(PROJECT.PATH[["data"]],
                    "gene_model.RData"))
                save(the.counts,exon.data,gene.data,sample.list,count.type,
                    file=file.path(PROJECT.PATH$data,"gene_model.RData"),
                    compress=TRUE)
            }
        }
        else # Retrieved gene model and/or previous analysis
            the.counts <- counts
            
        # Exclude any samples not wanted (when e.g. restoring a previous project
        # and having determined that some samples are of bad quality
        if (!is.null(exclude.list) && !is.na(exclude.list))
        {
            for (n in names(exclude.list)) {
                sample.list[[n]] <- setdiff(sample.list[[n]],
                    exclude.list[[n]])
                if (length(sample.list[[n]])==0) # Removed whole condition
                    sample.list[n] <- NULL
            }
            the.counts <- the.counts[unlist(sample.list)]
        }

        # Apply exon filters
        if (!is.null(exon.filters))
        {
            exon.filter.out <- filter.exons(the.counts,gene.data,sample.list,
                exon.filters)
            exon.filter.result <- exon.filter.out$result
            exon.filter.flags <- exon.filter.out$flags
        }
        else
            exon.filter.result <- exon.filter.flags <- NULL
        
        disp("Summarizing count data...")
        the.gene.counts <- the.exon.lengths <- vector("list",
            length(unlist(sample.list)))
        names(the.gene.counts) <- names(the.exon.lengths) <- names(the.counts)
        for (n in names(the.gene.counts))
        {
            the.gene.counts[[n]] <- wapply(multic,the.counts[[n]],
                function(x) return(sum(x$count)))
            the.exon.lengths[[n]] <- wapply(multic,the.counts[[n]],
                function(x) return(sum(x$length)))
            the.gene.counts[[n]] <- do.call("c",the.gene.counts[[n]])
            the.exon.lengths[[n]] <- do.call("c",the.exon.lengths[[n]])
        }
        gene.counts <- do.call("cbind",the.gene.counts)
        gene.length <- the.exon.lengths[[1]] # Based on the sum of their exon lengths
        names(gene.length) <- rownames(gene.data)
        
        # In case there are small differences between annotation data and external 
        # file, due to e.g. slightly different Ensembl versions
        gene.data <- gene.data[rownames(gene.counts),]
        total.gene.data <- gene.data # We need this for some total stats
    }
    else if (count.type=="gene")
    {
        if (!from.previous)
        {
            if (annotation=="download")
            {
                disp("Downloading gene annotation for ",org,"...")
                gene.data <- get.annotation(org,count.type,refdb)
            }
            else if (annotation=="embedded")
            {
                # The following should work if annotation elements are arranged 
                # in MeV-like data style
                if (!is.data.frame(counts))
                {
                    disp("Reading counts file ",counts.name,"...")
                    gene.counts <- read.delim(counts)
                }
                else
                    gene.counts <- counts
                rownames(gene.counts) <- as.character(gene.counts[,id.col])
                all.cols <- 1:ncol(gene.counts)
                sam.cols <- match(unlist(sample.list),colnames(gene.counts))
                sam.cols <- sam.cols[which(!is.na(sam.cols))]
                ann.cols <- all.cols[-sam.cols]
                gene.data <- gene.counts[,ann.cols]
                gene.counts <- gene.counts[,sam.cols]
                colnames(gene.data)[id.col] <- "gene_id"
                if (!is.na(gc.col))
                {
                    colnames(gene.data)[gc.col] <- "gc_content"
                    if (max(gene.data$gc_content<=1)) # Is already divided
                        gene.data$gc_content = 100*gene.data$gc_content
                }
                if (!is.na(name.col)) colnames(gene.data)[name.col] <- 
                    "gene_name"
                if (!is.na(bt.col)) colnames(gene.data)[bt.col] <- "biotype"
            }
            else # Reading from external file, similar to embedded
            {
                if (!is.data.frame(counts))
                {
                    disp("Reading counts file ",counts.name,"...")
                    gene.counts <- read.delim(counts)
                }
                else
                    gene.counts <- counts
                rownames(gene.counts) <- as.character(gene.counts[,id.col])
                disp("Reading external gene annotation for ",org," from ",
                    annotation,"...")
                gene.data <- read.delim(annotation)
                rownames(gene.data) <- as.character(gene.data$gene_id)
                gene.data <- gene.data[rownames(gene.counts),]
                if (max(gene.data$gc_content)<=1) # Is already divided
                    gene.data$gc_content = 100*gene.data$gc_content
            }
        }
        else
        {
            gene.counts <- tmp.env$gene.counts
            gene.data <- tmp.env$gene.data
        }
        
        total.gene.data <- gene.data # We need this for some total stats
        exon.filter.result <- NULL

        if (annotation!="embedded" & !from.previous) # Else everything is provided and done
        {
            if (!is.null(counts)) # Otherwise it's coming ready from read2count
            {
                if (!is.data.frame(counts)) # Else it's already here
                {
                    disp("Reading counts file ",counts.name,"...")
                    gene.counts <- read.delim(counts)
                }
                else # Already a data frame as input
                    gene.counts <- counts
                rownames(gene.counts) <- as.character(gene.counts[,id.col])
                gene.counts <- gene.counts[,unlist(sample.list,
                    use.names=FALSE)]
            }
            else # Coming from read2count
            {
                if (from.raw) # Double check
                {
                    r2c <- read2count(the.list,gene.data,file.type,
                        multic=multic)
                    gene.counts <- r2c$counts
                    if (is.null(libsize.list))
                        libsize.list <- r2c$libsize
                    if (export.counts.table) {
                        disp("Exporting raw read counts table to ",
                            file.path(PROJECT.PATH[["lists"]],
                            "raw_counts_table.txt.gz"))
                        res.file <- file.path(PROJECT.PATH[["lists"]],
                            "raw_counts_table.txt.gz")
                        gzfh <- gzfile(res.file,"w")
                        write.table(cbind(gene.data[rownames(gene.counts),],
                            gene.counts),gzfh,sep="\t",row.names=FALSE,
                            quote=FALSE)
                        close(gzfh)
                    }
                }
            }
        }

        gene.data <- gene.data[rownames(gene.counts),]
        gene.length <- gene.data$end - gene.data$start # Based on total gene lengths
        names(gene.length) <- rownames(gene.data)
        
        # Exclude any samples not wanted (when e.g. restoring a previous project
        # and having determined that some samples are of bad quality
        if (!is.null(exclude.list) && !is.na(exclude.list))
        {
            for (n in names(exclude.list)) {
                sample.list[[n]] <- setdiff(sample.list[[n]],
                    exclude.list[[n]])
                if (length(sample.list[[n]])==0) # Removed whole condition
                    sample.list[n] <- NULL
            }
            gene.counts <- gene.counts[,unlist(sample.list,use.names=FALSE)]
        }
        
        if (save.gene.model)
        {
            disp("Saving gene model to ",file.path(PROJECT.PATH[["data"]],
                "gene_model.RData"))
            save(gene.counts,gene.data,sample.list,count.type,
                file=file.path(PROJECT.PATH$data,"gene_model.RData"),
                compress=TRUE)
        }
    }

    # Transform GC-content and biotype
    gene.data$gc_content <- as.numeric(gene.data$gc_content)/100
    if (is.null(gene.data$biotype))
        gene.data$biotype <- rep(NA,nrow(gene.data))
    names(gene.length) <- rownames(gene.counts)
    attr(gene.data,"gene.length") <- gene.length

    ############################################################################
    # BEGIN FILTERING SECTION
    ############################################################################

    # GC bias is NOT alleviated if we do not remove the zeros!!!
    disp("Removing genes with zero counts in all samples...")
    the.zeros <- which(apply(gene.counts,1,filter.low,0))
    if (length(the.zeros)>0)
    {
        # Store the filtered, maybe we do some stats
        gene.counts.zero <- gene.counts[the.zeros,]
        gene.data.zero <- gene.data[the.zeros,]
        attr(gene.data.zero,"gene.length") <- gene.length[the.zeros]
        the.zero.names <- rownames(gene.data)[the.zeros]
        # Then remove
        gene.counts <- gene.counts[-the.zeros,]
        gene.data <- gene.data[-the.zeros,]
        attr(gene.data,"gene.length") <- gene.length[-the.zeros]
    }
    else
        gene.counts.zero <- gene.data.zero <- the.zero.names <- NULL

    # Store un-normalized gene counts for export purposes
    gene.counts.unnorm <- gene.counts

    # Apply filtering prior to normalization if desired
    if (when.apply.filter=="prenorm")
    {
        # However, a first round of normalization has to be performed in order to 
        # get proper expression filters
        disp("Prefiltering normalization with: ",normalization)
        switch(normalization,
            edaseq = {
                temp.genes <- normalize.edaseq(gene.counts,sample.list,
                    norm.args,gene.data,output="matrix")
            },
            deseq = {
                temp.genes <- normalize.deseq(gene.counts,sample.list,norm.args,
                    output="matrix")
            },
            edger = {
                temp.genes <- normalize.edger(gene.counts,sample.list,norm.args,
                    output="matrix")
            },
            noiseq = {
                temp.genes <- normalize.noiseq(gene.counts,sample.list,
                    norm.args,gene.data,log.offset,output="matrix")
            },
            nbpseq = {
                temp.genes <- normalize.nbpseq(gene.counts,sample.list,
                    norm.args,libsize.list,output="matrix")
            },
            none = {
                # In case some external normalization is applied (e.g. equal read 
                # counts from all samples)
                temp.genes <- gene.counts
            }
        )
        
        # Now filter
        if (!is.null(gene.filters))
        {
            gene.filter.out <- filter.genes(temp.genes,gene.data,gene.filters)
            gene.filter.result <- gene.filter.out$result
            gene.filter.cutoff <- gene.filter.out$cutoff
            gene.filter.flags <- gene.filter.out$flags
        }
        else
            gene.filter.result <- gene.filter.cutoff <-
                gene.filter.flags <- NULL

        # Unify the filters and filter
        the.dead.genes <- list(
            gene.filter.result$expression$median,
            gene.filter.result$expression$mean,
            gene.filter.result$expression$quantile,
            gene.filter.result$expression$known,
            gene.filter.result$expression$custom
        )
        the.dead <- unique(unlist(c(gene.filter.result,exon.filter.result)))
        # Some genes filtered by zero, were present in exon filters, not yet applied
        if (count.type=="exon")
            the.dead <- setdiff(the.dead,the.zero.names)
        
        # All method specific objects are row-index subsettable
        if (length(the.dead)>0)
        {
            # Store the filtered for later export or some stats
            gene.counts.dead <- gene.counts[the.dead,]
            gene.counts.unnorm <- gene.counts.unnorm[the.dead,]
            gene.data.dead <- gene.data[the.dead,]
            attr(gene.data.dead,"gene.length") <- attr(gene.data,
                "gene.length")[the.dead]
            # Now filter
            the.dead.ind <- match(the.dead,rownames(gene.counts))
            gene.counts.expr <- gene.counts[-the.dead.ind,]
            gene.data.expr <- gene.data[-the.dead.ind,]
            attr(gene.data.expr,"gene.length") <- attr(gene.data,
                "gene.length")[-the.dead.ind]
        }
        else
        {
            gene.counts.expr <- gene.counts
            gene.data.expr <- gene.data
            gene.counts.dead <- gene.data.dead <- gene.counts.unnorm <- NULL
        }

        disp("Normalizing with: ",normalization)
        switch(normalization,
            edaseq = {
                norm.genes <- normalize.edaseq(gene.counts.expr,sample.list,
                    norm.args,gene.data.expr,output="matrix")
            },
            deseq = {
                norm.genes <- normalize.deseq(gene.counts.expr,sample.list,
                    norm.args,output="native")
            },
            edger = {
                norm.genes <- normalize.edger(gene.counts.expr,sample.list,
                    norm.args,output="native")
            },
            noiseq = {
                norm.genes <- normalize.noiseq(gene.counts.expr,sample.list,
                    norm.args,gene.data.expr,log.offset,output="matrix")
            },
            nbpseq = {
                norm.genes <- normalize.nbpseq(gene.counts.expr,sample.list,
                    norm.args,libsize.list,output="native")
            },
            none = {
                norm.genes <- gene.counts.expr
            }
        )
        norm.genes.expr <- norm.genes
    }
    else if (when.apply.filter=="postnorm")
    {
        # Apply filtering prior to normalization if desired (default)
        disp("Normalizing with: ",normalization)
        switch(normalization,
            edaseq = {
                norm.genes <- normalize.edaseq(gene.counts,sample.list,
                    norm.args,gene.data,output="matrix")
            },
            deseq = {
                norm.genes <- normalize.deseq(gene.counts,sample.list,norm.args,
                    output="native")
            },
            edger = {
                norm.genes <- normalize.edger(gene.counts,sample.list,norm.args,
                    output="native")
            },
            noiseq = {
                norm.genes <- normalize.noiseq(gene.counts,sample.list,
                    norm.args,gene.data,log.offset,output="matrix")
            },
            nbpseq = {
                norm.genes <- normalize.nbpseq(gene.counts,sample.list,
                    norm.args,libsize.list,output="native")
            },
            none = {
                norm.genes <- gene.counts
            }
        )
    
        switch(class(norm.genes),
            CountDataSet = { # Has been normalized with DESeq
                temp.matrix <- round(counts(norm.genes,normalized=TRUE))
            },
            DGEList = { # Has been normalized with edgeR
                # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
                scl <- norm.genes$samples$lib.size *
                    norm.genes$samples$norm.factors
                temp.matrix <- round(t(t(norm.genes$counts)/scl)*mean(scl))
            },
            matrix = { # Has been normalized with EDASeq or NOISeq or nothing
                temp.matrix <- norm.genes
            },
            data.frame = { # Has been normalized with or nothing
                temp.matrix <- as.matrix(norm.genes)
            },
            list = { # Has been normalized with NBPSeq and main method was "nbpseq"
                temp.matrix <- as.matrix(round(sweep(norm.genes$counts,2,
                    norm.genes$norm.factors,"*")))
            },
            nbp = { # Has been normalized with NBPSeq and main method was "nbsmyth"
                 temp.matrix <- as.matrix(round(norm.genes$pseudo.counts))
            }
        )

        # Implement gene filters after normalization
        if (!is.null(gene.filters)) {
            gene.filter.out <- filter.genes(temp.matrix,gene.data,gene.filters)
            gene.filter.result <- gene.filter.out$result
            gene.filter.cutoff <- gene.filter.out$cutoff
            gene.filter.flags <- gene.filter.out$flags
        }
        else
            gene.filter.result <- gene.filter.cutoff <-
                gene.filter.flags <- NULL

        # Unify the filters and filter
        the.dead.genes <- list(
            gene.filter.result$expression$median,
            gene.filter.result$expression$mean,
            gene.filter.result$expression$quantile,
            gene.filter.result$expression$known,
            gene.filter.result$expression$custom
        )
        #gene.filter.result$expression <- Reduce("union",the.dead.genes)
        the.dead <- unique(unlist(c(gene.filter.result,exon.filter.result)))
        # Some genes filtered by zero, were present in exon filters, not yet applied
        if (count.type=="exon")
            the.dead <- setdiff(the.dead,the.zero.names)
        
        # All method specific object are row-index subsettable
        if (length(the.dead)>0)
        {
            # Store the filtered for later export or some stats
            gene.counts.dead <- temp.matrix[the.dead,]
            gene.counts.unnorm <- gene.counts.unnorm[the.dead,]
            gene.data.dead <- gene.data[the.dead,]
            attr(gene.data.dead,"gene.length") <- attr(gene.data,
                "gene.length")[the.dead]
            # Now filter
            the.dead.ind <- match(the.dead,rownames(temp.matrix))
            switch(class(norm.genes),
                CountDataSet = {
                    norm.genes.expr <- norm.genes[-the.dead.ind,]
                },
                DGEList = { # edgeR bug???
                    norm.genes.expr <- norm.genes[-the.dead.ind,]
                    norm.genes.expr$AveLogCPM <-
                        norm.genes.expr$AveLogCPM[-the.dead.ind]
                },
                matrix = { # Has been normalized with EDASeq or NOISeq
                    norm.genes.expr <- norm.genes[-the.dead.ind,]
                },
                data.frame = { # Has been normalized with EDASeq or NOISeq
                    norm.genes.expr <- as.matrix(norm.genes[-the.dead.ind,])
                },
                list = { # Has been normalized with NBPSeq, main.method="nbpseq"
                    norm.genes.expr <- norm.genes
                    norm.genes.expr$counts <-
                        as.matrix(norm.genes.expr$counts[-the.dead.ind,])
                    norm.genes.expr$rel.frequencies <- 
                        norm.genes.expr$rel.frequencies[-the.dead.ind,]
                    norm.genes.expr$tags <-
                        as.matrix(norm.genes.expr$tags[-the.dead.ind,])
                },
                nbp = {
                    norm.genes.expr <- norm.genes
                    norm.genes.expr$counts <-
                        as.matrix(norm.genes.expr$counts[-the.dead.ind,])
                    norm.genes.expr$pseudo.counts <- 
                        as.matrix(norm.genes.expr$pseudo.counts[-the.dead.ind,])
                    norm.genes.expr$pseudo.lib.sizes <- 
                        colSums(as.matrix(norm.genes.expr$pseudo.counts))*
                            rep(1,dim(norm.genes.expr$counts)[2])
                }
            )
            gene.counts.expr <- gene.counts[rownames(norm.genes.expr),]
            gene.data.expr <- gene.data[-the.dead.ind,]
            attr(gene.data.expr,"gene.length") <-
                attr(gene.data,"gene.length")[-the.dead.ind]
            
        }
        else
        {
            norm.genes.expr <- norm.genes
            gene.counts.expr <- gene.counts
            gene.data.expr <- gene.data
            gene.counts.dead <- gene.data.dead <- gene.counts.unnorm <- NULL
        }
    }
    
    # Store the final filtered, maybe we do some stats
    gene.data.filtered <- rbind(gene.data.zero,gene.data.dead)
    if (!is.null(gene.data.filtered) && nrow(gene.data.filtered)>0)
    {
        disp(nrow(gene.data.filtered)," genes filtered out")
        if (!is.null(gene.data.zero) && nrow(gene.data.zero)>0)
            attr(gene.data.filtered,"gene.length") <- c(attr(gene.data.zero,
                "gene.length"),attr(gene.data.dead,"gene.length"))
        else
            attr(gene.data.filtered,"gene.length") <-
                attr(gene.data.dead,"gene.length")
    }
    if (!is.null(gene.filters) || !is.null(exon.filters))
        disp(nrow(gene.data.expr)," genes remain after filtering")

    ############################################################################
    # END FILTERING SECTION
    ############################################################################

    # There is a small case that no genes are left after filtering...
    if(any(dim(norm.genes.expr)==0))
        stopwrap("No genes left after gene and/or exon filtering! Try again ",
            "with no filtering or less strict filter rules...")

    # Run the statistical test, norm.genes is always a method-specific object,
    # handled in the metaseqr.stat.R stat.* functions
    cp.list <- vector("list",length(contrast))
    names(cp.list) <- contrast
    contrast.list <- make.contrast.list(contrast,sample.list)
    for (n in names(cp.list))
    {
        cp.list[[n]] <- vector("list",length(statistics))
        names(cp.list[[n]]) <- statistics
    }
    for (alg in statistics)
    {
        disp("Running statistical tests with: ",alg)    
        switch(alg,
            deseq = {
                p.list <- stat.deseq(norm.genes.expr,sample.list,contrast.list,
                    stat.args[[alg]])
                if (!is.na(pcut)) {
                    for (con in names(contrast.list))
                        disp("  Contrast ",con,": found ",
                            length(which(p.list[[con]]<=pcut))," genes")
                }
            },
            edger = {
                p.list <- stat.edger(norm.genes.expr,sample.list,contrast.list,
                    stat.args[[alg]])
                if (!is.na(pcut)) {
                    for (con in names(contrast.list))
                        disp("  Contrast ",con,": found ",
                            length(which(p.list[[con]]<=pcut))," genes")
                }
            },
            noiseq = {
                p.list <- stat.noiseq(norm.genes.expr,sample.list,contrast.list,
                    stat.args[[alg]],gene.data.expr,log.offset)
                if (!is.na(pcut)) {
                    for (con in names(contrast.list))
                        disp("  Contrast ",con,": found ",
                            length(which(p.list[[con]]<=pcut))," genes")
                }
            },
            bayseq = {
                p.list <- stat.bayseq(norm.genes.expr,sample.list,contrast.list,
                    stat.args[[alg]],libsize.list)
                if (!is.na(pcut)) {
                    for (con in names(contrast.list))
                        disp("  Contrast ",con,": found ",
                            length(which(p.list[[con]]<=pcut))," genes")
                }
            },
            limma = {
                p.list <- stat.limma(norm.genes.expr,sample.list,contrast.list,
                    stat.args[[alg]])
                if (!is.na(pcut)) {
                    for (con in names(contrast.list))
                        disp("  Contrast ",con,": found ",
                            length(which(p.list[[con]]<=pcut))," genes")
                }
            },
            nbpseq = {
                p.list <- stat.nbpseq(norm.genes.expr,sample.list,contrast.list,
                    stat.args[[alg]],libsize.list)
                if (!is.na(pcut)) {
                    for (con in names(contrast.list))
                        disp("  Contrast ",con,": found ",
                            length(which(p.list[[con]]<=pcut))," genes")
                }
            }
        )
        for (n in names(p.list))
            cp.list[[n]][[alg]] <- p.list[[n]]
    }
    for (n in names(cp.list))
        cp.list[[n]] <- do.call("cbind",cp.list[[n]])

    # Create the adjusted p-value matrices (if needed)
    if ("adj.p.value" %in% export.what)
    {
        adj.cp.list <- wapply(multic,cp.list,
            function(x,a) return(apply(x,2,p.adjust,a)),adjust.method)
        for (n in names(cp.list))
        {
            noi <- grep("noiseq",colnames(cp.list[[n]]))
            if (length(noi)>0)
            {
                # DESeq has not run in this case, FDR cannot be calculated
                if (length(strsplit(n,"_vs_")[[1]])==2)
                    adj.cp.list[[n]][,noi] <- rep(NA,nrow(cp.list[[n]]))
            }
        }
    }
    else
        adj.cp.list <- NULL

    # At this point, all method-specific objects must become matrices for exporting 
    # and plotting
    switch(class(norm.genes.expr),
        CountDataSet = { # Has been processed with DESeq
            norm.genes <- round(counts(norm.genes,normalized=TRUE))
            norm.genes.expr <- round(counts(norm.genes.expr,normalized=TRUE))
        },
        DGEList = { # Has been processed with edgeR
            # Trick found at http://cgrlucb.wikispaces.com/edgeR+spring2013
            scl.r <- norm.genes$samples$lib.size*norm.genes$samples$norm.factors
            norm.genes <- round(t(t(norm.genes$counts)/scl.r)*mean(scl.r))
            scl.n <- norm.genes.expr$samples$lib.size *
                norm.genes.expr$samples$norm.factors
            norm.genes.expr <- round(t(t(norm.genes.expr$counts)/scl.n) *
                mean(scl.n))
        },
        list = {
            norm.genes <- as.matrix(round(sweep(norm.genes$counts,2,
                norm.genes$norm.factors,"*")))
            norm.genes.expr <- as.matrix(round(sweep(norm.genes.expr$counts,2,
                norm.genes$norm.factors,"*")))
        },
        nbp = {
            norm.genes <- as.matrix(round(norm.genes$pseudo.counts))
            norm.genes.expr <- as.matrix(round(norm.genes.expr$pseudo.counts))
        }
        # We don't need the matrix case
    )

    # Now that everything is a matrix, export the normalized counts if asked
    if (export.counts.table) {
        disp("Exporting and compressing normalized read counts table to ",
            file.path(PROJECT.PATH[["lists"]],"normalized_counts_table.txt"))
        expo <- cbind(
            rbind(gene.data.expr,gene.data.filtered),
            rbind(norm.genes.expr,gene.counts.zero,gene.counts.dead)
        )
        res.file <- file.path(PROJECT.PATH[["lists"]],
            "normalized_counts_table.txt.gz")
        gzfh <- gzfile(res.file,"w")
        write.table(expo,gzfh,sep="\t",row.names=FALSE,quote=FALSE)
        close(gzfh)
    }

    # Calculate meta-statistics, if more than one statistical algorithm has been used
    if (length(statistics)>1)
    {
        sum.p.list <- meta.test(
            cp.list=cp.list,
            meta.p=meta.p,
            counts=norm.genes.expr,
            sample.list=sample.list,
            statistics=statistics,
            stat.args=stat.args,
            libsize.list=libsize.list,
            nperm=nperm,
            weight=weight,
            reprod=reprod,
            multic=multic
        )
    }
    # We assign the p-values from the only statistic used to sum.p.list in order
    # to use it for stat plots
    else
        sum.p.list <- cp.list
    # Useless for one statistics but just for safety
    if ("adj.meta.p.value" %in% export.what) 
        adj.sum.p.list <- wapply(multic,sum.p.list,
            function(x,a) return(p.adjust(x,a)),adjust.method)
    else
        adj.sum.p.list <- NULL

    ############################################################################
    # BEGIN EXPORT SECTION
    ############################################################################

    # Bind all the flags
    if (count.type=="gene")
        flags <- gene.filter.flags
    else if (count.type=="exon")
    {
        if (!is.null(exon.filter.flags))
        {
            flags <- cbind(gene.filter.flags,
                as.matrix(exon.filter.flags[rownames(gene.filter.flags),]))
            nams <- c(colnames(gene.filter.flags),colnames(exon.filter.flags))
            rownames(flags) <- rownames(gene.filter.flags)
            colnames(flags) <- nams
        }
        else
            flags <- gene.filter.flags
    }
    
    disp("Building output files...")
    if (out.list) out <- make.export.list(contrast) else out <- NULL
    if (report) html <- make.export.list(contrast) else html <- NULL
    if ("rpgm" %in% export.scale)
        fa <- attr(gene.data,"gene.length")
    else
        fa <- NULL
    if ("normalized" %in% export.values) {
        fac <- fa[rownames(norm.genes.expr)]
        norm.list <- make.transformation(norm.genes.expr,export.scale,fac,
            log.offset)
    }
    else
        norm.list <- NULL
    if ("raw" %in% export.values) {
        fac <- fa[rownames(gene.counts.expr)]
        raw.list <- make.transformation(gene.counts.expr,export.scale,fac,
            log.offset)
    }
    else
        raw.list <- NULL
    if ("flags" %in% export.what)
        good.flags <- flags[rownames(norm.genes.expr),]
    else
        good.flags <- NULL

    if (!is.null(gene.counts.zero) || !is.null(gene.counts.dead))
    {
        gene.counts.filtered <- rbind(gene.counts.zero,gene.counts.dead)
        gene.counts.unnorm.filtered <- rbind(gene.counts.zero,
            gene.counts.unnorm)
        if ("normalized" %in% export.values) {
            fac <- fa[rownames(gene.counts.filtered)]
            norm.list.filtered <- make.transformation(gene.counts.filtered,
                export.scale,fac,log.offset)
        }
        else
            norm.list.filtered <- NULL
        if ("raw" %in% export.values) {
            fac <- fa[rownames(gene.counts.unnorm.filtered)]
            raw.list.filtered <- make.transformation(
                gene.counts.unnorm.filtered,export.scale,fac,log.offset)
        }
        else
            raw.list.filtered <- NULL
        if ("flags" %in% export.what && !is.null(flags))
            all.flags <- rbind(
                matrix(1,nrow(gene.counts.zero),ncol(flags)),
                    flags[rownames(gene.counts.dead),]
            )
        else
            all.flags <- NULL
    }
    
    counter <- 1
    for (cnt in contrast)
    {
        disp("  Contrast: ",cnt)
        disp("    Adding non-filtered data...")
        the.export <- build.export(
            gene.data=gene.data.expr,
            raw.gene.counts=gene.counts.expr,
            norm.gene.counts=norm.genes.expr,
            flags=good.flags,
            sample.list=sample.list,
            cnt=cnt,
            statistics=statistics,
            raw.list=raw.list,
            norm.list=norm.list,
            p.mat=cp.list[[cnt]],
            adj.p.mat=adj.cp.list[[cnt]],
            sum.p=sum.p.list[[cnt]],
            adj.sum.p=adj.sum.p.list[[cnt]],
            export.what=export.what,
            export.scale=export.scale,
            export.values=export.values,
            export.stats=export.stats,
            log.offset=log.offset,
            report=report
        )

        # Adjust the export based on what statistics have been done and a possible 
        # p-value cutoff
        export <- the.export$text.table
        if (report)
            export.html <- the.export$html.table
        if (!is.na(pcut))
        {
            if (length(statistics)>1)
            {
                switch(meta.p,
                    fisher = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    fperm = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    whitlock = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    dperm.min = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    dperm.max = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    dperm.weight = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    minp = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    maxp = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    weight = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    pandora = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    simes = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    },
                    none = {
                        cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                    }
                )
                pp <- sum.p.list[[cnt]][cut.ind]
                export <- export[cut.ind,]
                export <- export[order(pp),]
                if (report)
                {
                    export.html <- export.html[cut.ind,]
                    export.html <- export.html[order(pp),]
                }
            }
            else
            {
                cut.ind <- which(sum.p.list[[cnt]]<=pcut)
                pp <- sum.p.list[[cnt]][cut.ind,]
                export <- export[cut.ind,]
                export <- export[order(pp),]
                if (report)
                {
                    export.html <- export.html[cut.ind,]
                    export.html <- export.html[order(pp),]
                }
            }
        }
        else
        {
            pp <- sum.p.list[[cnt]]
            export <- export[order(pp),]
            if (report)
                export.html <- export.html[order(pp),]
        }

        # Final safety trigger
        na.ind <- grep("NA",rownames(export))
        if (length(na.ind)>0)
        {
            export <- export[-na.ind,]
            if (report) export.html <- export.html[-na.ind,]
        }

        res.file <- file.path(PROJECT.PATH[["lists"]],
            paste("metaseqr_sig_out_",cnt,".txt.gz",sep=""))
        disp("    Writing output...")
        gzfh <- gzfile(res.file,"w")
        write.table(export,gzfh,quote=FALSE,row.names=FALSE,sep="\t")
        close(gzfh)
        if (out.list)
            out[[cnt]] <- export
        if (report)
        {
            the.html.header <- make.html.header(the.export$headers)
            if (!is.null(report.top)) {
                topi <- ceiling(report.top*nrow(export.html))
                the.html.rows <- make.html.rows(export.html[1:topi,])
            }
            else
                the.html.rows <- make.html.rows(export.html)
            the.html.body <- make.html.body(the.html.rows)
            the.html.table <- make.html.table(the.html.body,the.html.header,
                id=paste("table_",counter,sep=""))
            html[[cnt]] <- the.html.table
            counter <- counter+1
        }

        if (!is.null(gene.counts.zero) || !is.null(gene.counts.dead))
        {
            disp("    Adding filtered data...")
            the.export.filtered <- build.export(
                gene.data=gene.data.filtered,
                raw.gene.counts=gene.counts.unnorm.filtered,
                norm.gene.counts=gene.counts.filtered,
                flags=all.flags,
                sample.list=sample.list,
                cnt=cnt,
                statistics=statistics,
                raw.list=raw.list.filtered,
                norm.list=norm.list.filtered,
                export.what=export.what,
                export.scale=export.scale,
                export.values=export.values,
                export.stats=export.stats,
                log.offset=log.offset,
                report=FALSE
            )

            # Now we should be having the.export and the.export.filtered. We do not 
            # generate html output for filtered or total results just a compressed
            # text file. We thus have to append the.export$text.table and 
            # the.export.filtered$html.table before writing the final output...
            export.all <- rbind(the.export$text.table,
                the.export.filtered$text.table)
            # ...and order them somehow... alphabetically according to row names, as
            # the annotation might not have been bundled...
            export.all <- export.all[order(rownames(export.all)),]

            res.file <- file.path(PROJECT.PATH[["lists"]],paste(
                "metaseqr_all_out_",cnt,".txt.gz",sep=""))
            disp("    Writing output...")
            gzfh <- gzfile(res.file,"w")
            write.table(export.all,gzfh,quote=FALSE,row.names=FALSE,sep="\t")
            close(gzfh)
        }
    }

    ############################################################################
    # END EXPORT SECTION
    ############################################################################

    ############################################################################
    # BEGIN PLOTTING SECTION
    ############################################################################
    
    if (!is.null(qc.plots))
    {
        disp("Creating quality control graphs...")
        plots <- list(
            raw=c("mds","biodetection","countsbio","saturation","readnoise",
                "correl","pairwise"),
            norm=c("boxplot","gcbias","lengthbias","meandiff","meanvar",
                "rnacomp"),
            stat=c("deheatmap","volcano","biodist"),
            other=c("filtered"),
            venn=c("venn")
        )
        fig.raw <- fig.unorm <- fig.norm <- fig.stat <- fig.other <- fig.venn <- 
            vector("list",length(fig.format))
        names(fig.raw) <- names(fig.unorm) <- names(fig.norm) <-
            names(fig.stat) <- names(fig.other) <- names(fig.venn) <-
            fig.format
        for (fig in fig.format)
        {
            disp("Plotting in ",fig," format...")
            fig.raw[[fig]] <- diagplot.metaseqr(gene.counts,sample.list,
                annotation=gene.data,diagplot.type=intersect(qc.plots,
                plots$raw),is.norm=FALSE,output=fig,path=PROJECT.PATH$qc)
            fig.unorm[[fig]] <- diagplot.metaseqr(gene.counts,sample.list,
                annotation=gene.data,diagplot.type=intersect(qc.plots,
                plots$norm),is.norm=FALSE,output=fig,
                path=PROJECT.PATH$normalization)
            if (when.apply.filter=="prenorm") # The annotation dimensions change...
                fig.norm[[fig]] <- diagplot.metaseqr(norm.genes,sample.list,
                    annotation=gene.data.expr,diagplot.type=intersect(qc.plots,
                    plots$norm),is.norm=TRUE,output=fig,
                    path=PROJECT.PATH$normalization) 
            else if (when.apply.filter=="postnorm")
                fig.norm[[fig]] <- diagplot.metaseqr(norm.genes,sample.list,
                    annotation=gene.data,diagplot.type=intersect(qc.plots,
                    plots$norm),is.norm=TRUE,output=fig,
                    path=PROJECT.PATH$normalization)
            fig.stat[[fig]] <- diagplot.metaseqr(norm.genes.expr,sample.list,
                annotation=gene.data.expr,contrast.list=contrast.list,
                p.list=sum.p.list,thresholds=list(p=pcut,f=1),
                diagplot.type=intersect(qc.plots,plots$stat),is.norm=TRUE,
                output=fig,path=PROJECT.PATH$statistics)
            if (!is.null(gene.data.filtered))
                fig.other[[fig]] <- diagplot.metaseqr(gene.data.filtered,
                    sample.list,annotation=total.gene.data,
                    diagplot.type=intersect(qc.plots,plots$other),
                    is.norm=FALSE,output=fig,path=PROJECT.PATH$qc)
            else fig.other[[fig]] <- NULL
            if ("venn" %in% qc.plots)
                fig.venn[[fig]] <- diagplot.metaseqr(norm.genes.expr,
                    sample.list,annotation=gene.data.expr,
                    contrast.list=contrast.list,
                    p.list=cp.list,thresholds=list(p=pcut,f=1),
                    diagplot.type=intersect(qc.plots,plots$venn),
                    output=fig,path=PROJECT.PATH$statistics)
        }
    }

    ############################################################################
    # END PLOTTING SECTION
    ############################################################################

    ############################################################################
    # BEGIN REPORTING SECTION
    ############################################################################
    
    if (report)
    {
        disp("Creating HTML report...")
        if (!is.null(qc.plots))
        {
            # First create zip archives of the figures
            disp("Compressing figures...")
            zipfiles <- file.path(PROJECT.PATH$plots,paste("metaseqr_figures_",
                fig.format,".zip",sep=""))
            names(zipfiles) <- fig.format
            for (f in fig.format)
            {
                files <- c(
                    dir(PROJECT.PATH$qc,pattern=paste(".",f,sep=""),
                        full.names=TRUE),
                    dir(PROJECT.PATH$normalization,pattern=paste(".",f,sep=""),
                        full.names=TRUE),
                    dir(PROJECT.PATH$statistics,pattern=paste(".",f,sep=""),
                        full.names=TRUE)
                )
                zip(zipfiles[f],files)
            }
            # Then create the final figure variables which brew will find...
            fig.raw <- fig.raw[["png"]]
            fig.unorm <- fig.unorm[["png"]]
            fig.norm <- fig.norm[["png"]]
            fig.stat <- fig.stat[["png"]]
            fig.other <- fig.other[["png"]]
            fig.venn <- fig.venn[["png"]]
        }
        else
            fig.raw <- fig.unorm <- fig.norm <- fig.stat <- fig.other <-
                fig.venn <- NULL

        if (tolower(report.template)=="default")
        {
            if (exists("TEMPLATE"))
            {
                report.template=list(
                    html=file.path(TEMPLATE,"metaseqr_report.html"),
                    css=file.path(TEMPLATE,"styles.css"),
                    logo=file.path(TEMPLATE,"logo.png"),
                    loader=file.path(TEMPLATE,"loader.gif")
                )
            }
            else
                report.template=list(html=NULL,css=NULL,logo=NULL,
                    loader=NULL)
        }

        if (!is.null(report.template$html))
        {
            if (file.exists(report.template$html))
            {
                template <- report.template$html
                has.template <- TRUE
            }
            else
            {
                warnwrap(paste("The template file",report.template$html,
                    "was not ","found! The HTML report will NOT be generated."))
                has.template <- FALSE
            }
        }
        else
        {
            warnwrap(paste("The report option was enabled but no template ",
                "file is provided! The HTML report will NOT be generated."))
            has.template <- FALSE
        }
        if (!is.null(report.template$css))
        {
            if (file.exists(report.template$css))
                file.copy(from=report.template$css,to=PROJECT.PATH$media)
            else
                warnwrap(paste("The stylesheet file",report.template$css,
                    "was not ","found! The HTML report will NOT be styled."))
        }
        else
            warnwrap(paste("The report stylesheet file was not provided! The ",
                "HTML report will NOT be styled."))
        if (!is.null(report.template$logo))
        {
            if (file.exists(report.template$logo))
                file.copy(from=report.template$logo,to=PROJECT.PATH$media)
            else
                warnwrap(paste("The report logo image",report.template$logo,
                    "was not found!"))
        }
        else
            warnwrap(paste("The report logo image was not provided!"))
        if (!is.null(report.template$loader))
        {
            if (file.exists(report.template$loader))
                file.copy(from=report.template$loader,to=PROJECT.PATH$media)
            else
                warnwrap(paste("The report logo image",report.template$loader,
                    "was not found!"))
        }
        else
            warnwrap(paste("The report loader image was not provided!"))
        if (has.template)
        {
            exec.time <- elap2human(TB)
            TEMP <- environment()
            brew(
                file=report.template$html,
                #output=file.path(PROJECT.PATH$main,paste(
                #    basename(PROJECT.PATH$main),"html",sep=".")),
                output=file.path(PROJECT.PATH$main,"index.html"),
                envir=TEMP
            )
        }
    }

    ############################################################################
    # END REPORTING SECTION
    ############################################################################
    
    disp("\n",strftime(Sys.time()),": Data processing finished!\n")
    exec.time <- elap2human(TB)
    disp("\n","Total processing time: ",exec.time,"\n\n")
    
    if (out.list) return(list(data=out,html=html))
} # End metaseqr

#' Assemble a gene model based on exon counts
#'
#' This function assembles gene models (single genes, not isoforms) based on the
#' input exon read counts file and a gene annotation data frame, either from an
#' external file provided by the user, or with the \code{\link{get.annotation}}
#' function. The \code{gene.data} argument should have a specific format and for
#' this reason it's better to use one of the two aforementioned ways to supply it.
#' This function is intended mostly for internal use but can be used if the
#' requirements are met.
#'
#' @param exon.counts the exon counts data frame produced by reading the exon read
#' counts file.
#' @param sample.list the list containing condition names and the samples under
#' each condition.
#' @param gene.data an annotation data frame from the same organism as
#' \code{exon.counts} (such the ones produced by \code{get.annotation}).
#' @param multic a logical value indicating the presence of multiple cores. Defaults
#' to \code{FALSE}. Do not change it if you are not sure whether package parallel
#' has been loaded or not.
#' @return A named list where names represent samples. Each list member is a also
#' a named list where names correspond to gene ids and members are named vectors.
#' Each vector is named according to the exons corresponding to each gene and
#' contains the read counts for each exon. This structure is used for exon filtering
#' and assembling final gene counts in the metaseqr pipeline.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' @examples
#' \dontrun{
#' data("hg19.exon.data",package="metaseqR")
#' gene.data <- get.annotation("hg19","gene","ensembl")
#' reduced.gene.data <- reduce.gene.data(hg19.exon.counts,gene.data)
#' multic <- check.parallel(0.4)
#' gene.model <- construct.gene.model(hg19.exon.counts,sample.list.hg19,gene.data,
#'   multic)
#'}
construct.gene.model <- function(exon.counts,sample.list,
    gene.data,multic=FALSE) {
    the.counts <- vector("list",length(unlist(sample.list)))
    names(the.counts) <- unlist(sample.list,use.names=FALSE)
    the.genes <- as.character(unique(gene.data$gene_id))
    #the.exons <- as.character(unique(exon.counts$exon_id))
    #the.exons <- as.character(unique(exon.data[,id.col]))
    for (n in names(the.counts))
    {
        disp("  Separating exons per gene for ",n,"...")
        #the.counts[[n]] <- vector("list",length(the.genes))
        the.counts[[n]] <- the.genes
        names(the.counts[[n]]) <- the.genes
        the.counts[[n]] <- wapply(multic,the.counts[[n]],function(x,d,n) {
            tmp <- d[which(d$gene_id==x),c("start","end","exon_id",n)]
            xx <- tmp[,n]
            yy <- tmp$end - tmp$start
            names(xx) <- names(yy) <- tmp$exon_id
            return(list(count=xx,length=yy))
        },exon.counts,n)
    }
    return(the.counts)
}

#' Reduce the gene annotation in case of not all chromosomes present in counts
#'
#' This function reduces the gene annotation in case of exon reads and when the
#' data to be analyzed do not contain all the standard chromosomes of the genome
#' under investigation. This can greatly reduce processing time in these cases.
#'
#' @param exon.data the exon annotation already reduced to the size of the input
#' exon counts table.
#' @param gene.data an annotation data frame from the same organism as 
#' \code{exon.counts} (such the ones produced by \code{get.annotation}).
#' @return The \code{gene.data} annotation, reduced to have the same chromosomes
#' as in \code{exon.data}, or the original \code{gene.data} if \code{exon.data}
#' does contain the standard chromosomes.
#' @author Panagiotis Moulos
#' @export
#' @examples
#' \dontrun{
#' data("hg19.exon.data",package="metaseqR")
#' gene.data <- get.annotation("hg19","gene","ensembl")
#' reduced.gene.data <- reduce.gene.data(hg19.exon.counts,gene.data)
#'}
reduce.gene.data <- function(exon.data,gene.data) {
    exon.chrs <- unique(as.character(exon.data$chromosome))
    gene.chrs <- unique(as.character(gene.data$chromosome))
    if (length(exon.chrs)!=length(gene.chrs)) {
        m <- match(gene.data$chromosome,exon.chrs)
        gene.data <- gene.data[which(!is.na(m)),]
    }
    return(gene.data)
}

# Initialize environment
#
# Initializes metaseqr environmental variables. Internal use only.
#
# @author Panagiotis Moulos
#init.envar <- function() {
#    HOME <<- system.file(package="metaseqR")
#    SCRIPT <<- file.path(HOME,"R")
#    TEMPLATE <<- HOME
#    ANNOTATION <<- file.path(HOME,"data")
#}
