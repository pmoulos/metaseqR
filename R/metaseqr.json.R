#' MDS plot JSON exporter for the metaseqR package
#'
#' Non-exportable JSON exporter for \code{\link{diagplot.mds}}.
#'
#' @param obj A list holding MDS plot data. See \code{\link{diaplot.mds}}.
#' @param jl JavaScript charting library to export. Currently only \code{"highcharts"}
#' supported.
#' @return A JSON string.
#' @author Panagiotis Moulos
mdsToJSON <- function(obj,jl=c("highcharts")) {
    jl <- tolower(jl[1])
    x <- obj$x
    y <- obj$y
    xlim <- obj$xlim
    ylim <- obj$ylim
    samples <- obj$samples
    cols <- getColorScheme()
    
    # Construct series
    counter <- 0
    series <- vector("list",length(samples))
    names(series) <- names(samples)
    for (n in names(series)) {
        counter <- counter + 1
        series[[n]] <- list()
        series[[n]]$name=n
        series[[n]]$type="scatter"
        series[[n]]$color=cols$fill[counter]
        series[[n]]$marker <- list(
            lineWidth=1,
            states=list(
                hover=list(
                    enabled=TRUE,
                    lineColor=cols$border[counter]
                ),
                select=list(
                    fillColor=cols$selected[counter],
                    lineColor=cols$border[counter],
                    lineWidth=2
                )
            )
        )
        m <- match(samples[[n]],names(x))
        if (length(m)>0) {            
            series[[n]]$data <- make.highcharts.points(x[m],y[m])
        }
    }
    
    switch(jl,
        highcharts = {
            point.format=paste("<strong>Sample name: </strong>{point.name}<br>",
                "<strong>Principal coordinate 1: </strong>{point.x}<br>",
                "<strong>Principal coordinate 2: </strong>{point.y}",sep="")
                
                json <- toJSON(
                    list(
                        chart=list(
                        type="scatter",
                        zoomType="xy"
                    ),
                    title=list(
                        text=paste("Multidimensional Scaling")
                    ),
                    xAxis=list(
                        title=list(
                            useHTML=TRUE,
                            text="1<sup>st</sup> Principal Coordinate",
                            margin=20,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        min=round(xlim[1],3),
                        max=round(xlim[2],3)
                    ),
                    yAxis=list(
                        title=list(
                            useHTML=TRUE,
                            text="2<sup>nd</sup> Principal Coordinate",
                            margin=25,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        min=round(ylim[1],3),
                        max=round(ylim[2],3)
                    ),
                    plotOptions=list(
                        scatter=list(
                            allowPointSelect=TRUE,
                            tooltip=list(
                                headerFormat=paste("<span style=",
                                    "\"font-size:1.1em;color:{series.color};",
                                    "font-weight:bold\">{series.name}<br>",
                                    sep=""),
                                pointFormat=point.format
                            )
                        )
                    ),
                    series=unname(series)
                )
            )
        }
    )
    return(json)
}

#' Biodetection counts plot JSON exporter for the metaseqR package
#'
#' Non-exportable JSON exporter for \code{\link{diaplot.noiseq}}.
#'
#' @param obj A list holding boxplot data. See \code{\link{diaplot.noiseq}}.
#' @param by Can be \code{"sample"} to create biotypes boxplots per sample or
#' \code{"biotype"} to create samples boxplots per biotype.
#' @param jl JavaScript charting library to export. Currently only \code{"highcharts"}
#' supported.
#' @return A JSON string.
#' @author Panagiotis Moulos
countsBioToJSON <- function(obj,by=c("sample","biotype"),jl=c("highcharts")) {
    by <- tolower(by[1])
    jl <- tolower(jl[1])
    samples <- obj$samples
    status <- obj$status
    altnames <- obj$altnames
    counts <- round(2^obj$user$counts - 1)
    counts[counts==0] <- 0.001
    #counts <- obj$user$counts
    
    covars <- obj$user$covars
    biotypes <- unique(as.character(covars$biotype))
    if (!is.null(altnames))
        names(altnames) <- rownames(counts)
    
    grouped <- FALSE
    if (is.null(samples)) {
        if (is.null(colnames(counts)))
            samplenames <- paste("Sample",1:ncol(counts),sep=" ")
        else
            samplenames <- colnames(counts)
        samples <- list(Samples=nams)
    }
    else if (is.list(samples)) {
        samplenames <- unlist(samples,use.names=FALSE)
        grouped <- TRUE
    }
    
    # y label formatter for logarithmic axis
    y.label.formatter <- paste('function() {if(this.value === 0.001)',
        '{return 0;} else {return Highcharts.Axis.prototype.',
        'defaultLabelFormatter.call(this);}}',sep="")
        
    tooltip.point.formatter <- paste("function() {",
        "   var min = this.low === 0.001 ? 0 : this.low;" ,
        "   var q1 = this.q1 === 0.001 ? 0 : this.q1;" ,
        "   var med = this.median === 0.001 ? 0 : this.median;",
        "   var q3 = this.q3 === 0.001 ? 0 : this.q3;",
        "   var max = this.high === 0.001 ? 0 : this.high;",
        "   var str = 'Maximum: ' + max + '<br/>' +",
        "       'Upper quartile: ' + q3 + '<br/>' +",
        "       'Median: ' + med + '<br/>' +",
        "       'Lower quartile: ' + q1 + '<br/>' +",
        "       'Minimum: ' + min + '<br/>';",
        "   return  str;",
        "}",sep="")
        
    # Legend clicker
    boxplot.onclick <- paste("function() { ",
        "   var chart =  this.chart;",
        "   var outlier_id =  chart.get(this.name);",
        "   if (!outlier_id.visible) {",
        "       outlier_id.show();",
        "   } else {",
        "       outlier_id.hide();",
        "   }",
        "}",sep="")
    
    # Outliers tooltip
    if (is.null(obj$altnames)) {
        outlier.pointformat <- paste(
            '<strong>Sample {point.category}</strong><br/>',
            'Gene ID: {point.name}<br/>',
            'Value: {point.y}<br/>',sep=""
        )
    }
    else {
        outlier.pointformat <- paste(
            '<strong>Sample {point.category}</strong><br/>',
            'Gene ID: {point.name}<br/>',
            'Gene name: {point.alt_name}<br/>',
            'Value: {point.y}<br/>',
            sep=""
        )
    }
    
    if (by=="sample") {
        cols <- getColorScheme(length(biotypes))
        box.list <- json <- vector("list",length(samplenames))
        names(box.list) <- names(json) <- samplenames
        for (n in samplenames) {
            box.list[[n]] <- vector("list",length(biotypes))
            names(box.list[[n]]) <- biotypes
            for (b in biotypes)
                box.list[[n]][[b]] <- counts[covars$biotype==b,n]
            
            B <- boxplot(box.list[[n]],plot=FALSE)$stats
            colnames(B) <- biotypes
            o.list <- lapply(names(box.list[[n]]),function(x,M,b) {
                v <- b[,x]
                o <- which(M[[x]]<v[1] | M[[x]]>v[5])
                if (length(o)>0)
                    return(M[[x]][o])
                else
                    return(NULL)
            },box.list[[n]],B)
            names(o.list) <- biotypes
    
            # Data series
            BB <- matrix(0,nrow(B),ncol(B)) # Workaround of strange problem...
            colnames(BB) <- colnames(B)
            for (jj in 1:ncol(B))
                BB[,jj] <- round(B[,jj],3)
            d <- as.data.frame(BB)
            ids <- 0:(ncol(d)-1)
            d <- rbind(ids,d)
            names(ids) <- colnames(d)
            counter <- 0
            series <- vector("list",length(biotypes))
            names(series) <- biotypes
            for (s in names(series)) {
                counter <- counter + 1
                series[[s]] <- list()
                series[[s]]$name <- s
                series[[s]]$color <- cols$fill[counter]
                #series[[s]]$turboThreshold <- 10000
                series[[s]]$data <- list(unname(as.list(d[,s])))
                r <- round(d[,s])
                series[[s]]$tooltip=list(
                    pointFormat=paste('<strong>Population: ',
                        length(box.list[[n]][[s]]),'</strong><br/>',
                        'Maximum: ',r[6],'<br/>',
                        'Upper quartile: ',r[5],'<br/>',
                        'Median: ',r[4],'<br/>',
                        'Lower quartile: ',r[3],'<br/>',
                        'Minimum: ',r[2],'<br/>',sep="")
                )
            }
            
            # Outlier series (if any)
            counter <- 0
            outliers <- vector("list",length(biotypes))
            names(outliers) <- biotypes
            for (o in names(outliers)) {
                counter <- counter + 1
                outliers[[o]] <- list()
                outliers[[o]]$id <- o
                outliers[[o]]$name <- o
                outliers[[o]]$type <- "scatter"
                outliers[[o]]$showInLegend <- FALSE
                #outliers[[o]]$turboThreshold <- 10000
                outliers[[o]]$color <- cols$fill[counter]
                outliers[[o]]$marker <- list(
                    fillColor=cols$fill[counter],
                    symbol="circle",
                    lineWidth=1,
                    lineColor=cols$border[counter]
                )
                outliers[[o]]$data <- list()
                x <- rep(d[1,o],length(o.list[[o]]))
                names(x) <- names(o.list[[o]])
                if (is.null(obj$altnames)) {
                    outliers[[o]]$data <- 
                        make.highcharts.points(x,o.list[[o]])
                }
                else {
                    outliers[[o]]$data <- 
                        make.highcharts.points(x,o.list[[o]],
                            unname(altnames[names(x)]))
                }
            }
            
            json[[n]] <- switch(jl,
                highcharts = {
                    toJSON(
                        list(
                            chart=list(
                                type="boxplot"
                            ),
                            title=list(
                                text=paste("Biotype detection for sample ",n,
                                    sep="")
                            ),
                            legend=list(
                                enabled=TRUE,
                                itemHoverStyle=list(
                                    color="#B40000"
                                )
                            ),
                            xAxis=list(
                                categories=biotypes,
                                title=list(
                                    text="Biotype",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.2em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontWeight="bold"
                                    )
                                )
                            ),
                            yAxis=list(
                                type="logarithmic",
                                showFirstLabel=FALSE,
                                min=1e-4,
                                tickInterval=1,
                                title=list(
                                    useHTML=TRUE,
                                    #text="Read count (log<sub>2</sub>)",
                                    text="Expression (read count)",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em",
                                        fontWeight="bold"
                                    ),
                                    formatter=y.label.formatter
                                )
                            ),
                            plotOptions=list(
                                boxplot=list(
                                    fillColor="#F0F0E0",
                                    lineWidth=2,
                                    medianColor="#000000",
                                    medianWidth=3,
                                    stemColor="#000000",
                                    stemDashStyle="dash",
                                    stemWidth=1,
                                    whiskerColor="#000000",
                                    whiskerLength="75%",
                                    whiskerWidth=1,
                                    grouping=FALSE,
                                    tooltip=list(
                                        headerFormat=paste(
                                            '<span style="font-size:1.1em;',
                                            'color:{series.color};',
                                            'font-weight:bold">',
                                            '\u25CF </span>',
                                            '<span style="font-size:1.1em;',
                                            'font-weight:bold">',
                                            'Biotype {series.name}</span><br/>',
                                            sep=""
                                        )
                                    ),
                                    events=list(
                                        legendItemClick=boxplot.onclick
                                    )
                                ),
                                scatter=list(
                                    allowPointSelect=TRUE,
                                    tooltip=list(
                                        headerFormat=paste(
                                            '<span style="font-weight:bold;',
                                            'color:{series.color};">',
                                            '\u25CF </span>',
                                            '<span style="font-weight:bold">',
                                            'Biotype {series.name}</span><br/>',
                                            sep=""
                                        ),
                                        pointFormat=outlier.pointformat
                                    ),
                                    states=list(
                                        hover=list(
                                            marker=list(
                                                enabled=FALSE
                                            )
                                        )
                                    )
                                )
                            ),
                            series=c(unname(series),unname(outliers))
                        )
                    )
                }
            )
        }
        return(unquote_js_fun(json))
    }
    else if (by=="biotype") {
        cols <- getColorScheme(length(samples))
        box.list <- json <- vector("list",length(biotypes))
        names(box.list) <- names(json) <- biotypes
        for (b in biotypes) {
            box.list[[b]] <- vector("list",length(samplenames))
            names(box.list[[b]]) <- samplenames
            for (n in samplenames)
                box.list[[b]][[n]] <- counts[covars$biotype==b,n]
            
            B <- boxplot(box.list[[b]],plot=FALSE)$stats
            colnames(B) <- samplenames
            o.list <- lapply(names(box.list[[b]]),function(x,M,b) {
                v <- b[,x]
                o <- which(M[[x]]<v[1] | M[[x]]>v[5])
                if (length(o)>0)
                    return(M[[x]][o])
                else
                    return(NULL)
            },box.list[[b]],B)
            names(o.list) <- samplenames
    
            # Data series
            BB <- matrix(0,nrow(B),ncol(B)) # Workaround of strange problem...
            colnames(BB) <- colnames(B)
            for (jj in 1:ncol(B))
                BB[,jj] <- round(B[,jj],3)
            d <- as.data.frame(BB)
            ids <- 0:(ncol(d)-1)
            d <- rbind(ids,d)
            names(ids) <- colnames(d)
            counter <- 0
            series <- vector("list",length(samples))
            names(series) <- names(samples)
            for (s in names(series)) {
                counter <- counter + 1
                series[[s]] <- list()
                series[[s]]$name=s
                if (grouped)
                    series[[s]]$color=cols$fill[counter]
                else
                    series[[s]]$color=cols$fill[1]
                m <- match(samples[[s]],colnames(d))
                series[[s]]$data <- unname(as.list(d[,m]))
            }
            
            # Outlier series (if any)
            counter <- 0
            outliers <- vector("list",length(samples))
            names(outliers) <- names(samples)
            for (o in names(outliers)) {
                counter <- counter + 1
                outliers[[o]] <- list()
                outliers[[o]]$id <- o
                outliers[[o]]$name <- o
                outliers[[o]]$type <- "scatter"
                outliers[[o]]$showInLegend <- FALSE
                if (grouped) {
                    outliers[[o]]$color <- cols$fill[counter]
                    outliers[[o]]$marker <- list(
                        fillColor=cols$fill[counter],
                        symbol="circle",
                        lineWidth=1,
                        lineColor=cols$border[counter]
                    )
                }
                else {
                    outliers[[o]]$color <- cols$fill[1]
                    outliers[[o]]$marker <- list(
                        fillColor=cols$fill[1],
                        symbol="circle",
                        lineWidth=1,
                        lineColor=cols$border[1]
                    )
                }
                outliers[[o]]$data <- list()
                m <- match(samples[[o]],colnames(d))
                if (length(m)>0) {
                    for (i in m) {
                        x <- rep(d[1,i],length(o.list[[i]]))
                        names(x) <- names(o.list[[i]])
                        if (is.null(obj$altnames)) {
                            outliers[[o]]$data <- 
                                make.highcharts.points(x,o.list[[i]])
                        }
                        else {
                            outliers[[o]]$data <- c(outliers[[o]]$data,
                                make.highcharts.points(x,o.list[[i]],
                                unname(altnames)))
                        }
                    }
                }
            }
            
            json[[b]] <- switch(jl,
                highcharts = {
                    toJSON(
                        list(
                            chart=list(
                                type="boxplot"
                            ),
                            title=list(
                                text=paste("Detection for biotype ",b,
                                    " (population: ",lengths(box.list[[b]])[1],
                                    ")",sep="")
                            ),
                            legend=list(
                                enabled=TRUE
                            ),
                            xAxis=list(
                                categories=samplenames,
                                title=list(
                                    text="Sample name",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.2em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontWeight="bold"
                                    )
                                )
                            ),
                            yAxis=list(
                                type="logarithmic",
                                showFirstLabel=FALSE,
                                min=1e-4,
                                tickInterval=1,
                                title=list(
                                    text="Expression (read count)",
                                    margin=25,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em"
                                    )
                                ),
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em",
                                        fontWeight="bold"
                                    ),
                                    formatter=y.label.formatter
                                )
                            ),
                            plotOptions=list(
                                boxplot=list(
                                    fillColor="#F0F0E0",
                                    lineWidth=2,
                                    medianColor="#000000",
                                    medianWidth=3,
                                    stemColor="#000000",
                                    stemDashStyle="dash",
                                    stemWidth=1,
                                    whiskerColor="#000000",
                                    whiskerLength="75%",
                                    whiskerWidth=1,
                                    grouping=FALSE,
                                    tooltip=list(
                                        headerFormat=paste(
                                            '<span style="font-size:1.1em;',
                                            'color:{series.color};',
                                            'font-weight:bold">',
                                            '\u25CF </span>',
                                            '<span style="font-size:1.1em;',
                                            'font-weight:bold">',
                                            'Condition {series.name}</span>',
                                            '<br/>',
                                            '<span style="font-weight:bold;">',
                                            'Sample {point.key}',
                                            '</span><br/>',sep=""
                                        ),
                                        pointFormatter=tooltip.point.formatter
                                    ),
                                    events=list(
                                        legendItemClick=boxplot.onclick
                                    )
                                ),
                                scatter=list(
                                    allowPointSelect=TRUE,
                                    tooltip=list(
                                        headerFormat=paste(
                                            '<span style="font-weight:bold;',
                                            'color:{series.color};">',
                                            '\u25CF </span>',
                                            '<span style="font-weight:bold">',
                                            'Condition {series.name}</span>',
                                            '<br/>',sep=""
                                        ),
                                        pointFormat=outlier.pointformat
                                    ),
                                    states=list(
                                        hover=list(
                                            marker=list(
                                                enabled=FALSE
                                            )
                                        )
                                    )
                                )
                            ),
                            series=c(unname(series),unname(outliers))
                        )
                    )
                }
            )
        }
        return(unquote_js_fun(json))
    }
}

#' Biodetection barplot JSON exporter for the metaseqR package
#'
#' Non-exportable JSON exporter for \code{\link{diaplot.noiseq}}.
#'
#' @param obj A list holding boxplot data. See \code{\link{diaplot.noiseq}}.
#' @param jl JavaScript charting library to export. Currently only \code{"highcharts"}
#' supported.
#' @return A JSON string.
#' @author Panagiotis Moulos
bioDetectionToJSON <- function(obj,jl=c("highcharts")) {
    jl <- tolower(jl[1])
    samples <- obj$samples
    status <- obj$status
    plotdata <- obj$user$plotdata
    covars <- obj$user$covars
    biotypes <- unique(as.character(covars$biotype))
    
    if (!is.null(samples)&& is.list(samples)) {
        samplenames <- unlist(samples,use.names=FALSE)
        names(plotdata$biotables) <- samplenames
    }
    # Otherwise we are using the names present in the input object
    
    cols <- getColorScheme()
    bar.list <- json <- vector("list",length(samplenames))
    names(bar.list) <- names(json) <- samplenames
    for (n in samplenames) {
        # Data series
        series <- vector("list",3)
        names(series) <- c("genome","detectionVSgenome","detectionVSsample")
        series$genome <- list()
        series$genome$name <- "% in genome"
        series$genome$color <- cols$trans[1]
        series$genome$pointPlacement <- -0.2
        series$genome$data <- round(as.numeric(plotdata$genome),3)
        #series$genome$tooltip <- ""
        series$detectionVSgenome <- list()
        series$detectionVSgenome$name <- "% detected"
        series$detectionVSgenome$color <- cols$trans[2]
        series$detectionVSgenome$pointPlacement <- 0
        series$detectionVSgenome$data <- round(as.numeric(
            plotdata$biotables[[n]][1,]),3)
        #series$detectionVSgenome$tooltip
        series$detectionVSsample <- list()
        series$detectionVSsample$name <- "% in sample"
        series$detectionVSsample$color <- cols$trans[3]
        series$detectionVSsample$pointPlacement <- 0.2
        series$detectionVSsample$data <- round(as.numeric(
            plotdata$biotables[[n]][2,]),3)
        #series$detectionVSsample$tooltip
        
        json[[n]] <- switch(jl,
            highcharts = {
                toJSON(
                    list(
                        chart=list(
                            type="column"
                        ),
                        title=list(
                            text=paste("Comparative biotype detection for ",
                                "sample ",n,sep="")
                        ),
                        legend=list(
                            enabled=TRUE,
                            itemHoverStyle=list(
                                color="#B40000"
                            )
                        ),
                        tooltip=list(
                            shared=TRUE
                        ),  
                        xAxis=list(
                            categories=biotypes,
                            title=list(
                                text="Biotype",
                                margin=25,
                                style=list(
                                    color="#000000",
                                    fontSize="1.2em"
                                )
                            ),
                            labels=list(
                                style=list(
                                    color="#000000",
                                    fontWeight="bold"
                                )
                            )
                        ),
                        yAxis=list(
                            list(
                                min=0,
                                max=80,
                                title=list(
                                    text="% of abundant features",
                                    margin=20,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em"
                                    )
                                ),                          
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em",
                                        fontWeight="bold"
                                    )
                                )
                            ),
                            list(
                                min=0,
                                max=10,
                                title=list(
                                    text="% of non-abundant features",
                                    margin=20,
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em"
                                    )
                                ),                          
                                labels=list(
                                    style=list(
                                        color="#000000",
                                        fontSize="1.1em",
                                        fontWeight="bold"
                                    )
                                ),
                                opposite=TRUE
                            )
                        ),
                        plotOptions=list(
                            column=list(
                                grouping=FALSE,
                                shadow=FALSE,
                                groupPadding=0.3,
                                pointPadding=0.25,
                                tooltip=list(
                                    headerFormat=paste(
                                        '<span style="font-size:1.1em;',
                                        'font-weight:bold">',
                                        '{point.key}</span><br/>',sep=""
                                    )
                                )
                            )
                        ),
                        series=unname(series)
                    )
                )
            }
        )
    }
    return(json)
}

#' Read noise plot JSON exporter for the metaseqR package
#'
#' Non-exportable JSON exporter for \code{\link{diagplot.noiseq}}.
#'
#' @param obj A list holding plot data. See \code{\link{diaplot.noiseq}}.
#' @param jl JavaScript charting library to export. Currently only \code{"highcharts"}
#' supported.
#' @return A JSON string.
#' @author Panagiotis Moulos
readNoiseToJSON <- function(obj,jl=c("highcharts")) {
    jl <- tolower(jl[1])
    d <- obj$user
    samples <- obj$samples
    
    # Too many points for a lot of curves of interactive data
    if (nrow(d)>1000) {
        ii <- sort(sample(1:nrow(d),998))
        ii <- c(1,ii,nrow(d))
        d <- cbind(d[ii,1],d[ii,2:ncol(d)])
    }

    if (is.null(samples)) 
        samples <- 1:(ncol(d)-1)
    if (is.numeric(samples)) 
        samplenames = colnames(dat)[samples+1]
    if (is.list(samples))
        samplenames <- unlist(samples)
    
    cols <- getColorScheme(length(samplenames))
    
    # Construct series
    counter <- 0
    series <- vector("list",length(samplenames))
    names(series) <- samplenames
    for (n in names(series)) {
        counter <- counter + 1
        series[[n]] <- list()
        series[[n]]$name=n
        series[[n]]$color=cols$fill[counter]
        series[[n]]$data <- make.highcharts.points(d[,1],d[,n])
        series[[n]]$tooltip=list(
            headerFormat=paste("<span style=",
                "\"font-size:1.1em;color:{series.color};",
                "font-weight:bold\">{series.name}<br>",
                sep=""),
            pointFormat=NULL
        )
    }
    
    switch(jl,
        highcharts = {
                json <- toJSON(list(
                    chart=list(
                        type="line",
                        zoomType="xy"
                    ),
                    title=list(
                        text=paste("RNA-Seq mapped reads noise")
                    ),
                    xAxis=list(
                        useHTML=TRUE,
                        title=list(
                            text="% detected features",
                            margin=20,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        min=0,
                        max=100
                    ),
                    yAxis=list(
                        useHTML=TRUE,
                        title=list(
                            text="% of total reads",
                            margin=25,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        tickPositions=seq(0,110,10)
                    ),
                    plotOptions=list(
                        line=list(
                            allowPointSelect=TRUE,
                            lineWidth=1,
                            marker=list(
                                enabled=FALSE
                            ),
                            tooltip=list(
                                headerFormat=paste("<span style=",
                                    "\"font-size:1.1em;color:{series.color};",
                                    "font-weight:bold\">{series.name}<br>",
                                    sep=""),
                                pointFormat=NULL
                            ),
                            turboThreshold=50000
                        )
                    ),
                    series=unname(series)
                )
            )
        }
    )
    return(json)
}

#' Boxplots JSON exporter for the metaseqR package
#'
#' Non-exportable JSON exporter for \code{\link{diaplot.boxplot}}.
#'
#' @param obj A list holding boxplot data. See \code{\link{diaplot.boxplot}}.
#' @param jl JavaScript charting library to export. Currently only \code{"highcharts"}
#' supported.
#' @return A JSON string.
#' @author Panagiotis Moulos
boxplotToJSON <- function(obj,jl=c("highcharts")) {
    jl <- tolower(jl[1])
    b <- obj$plot
    name <- obj$samples
    status <- obj$status
    altnames <- obj$altnames
    o.list <- obj$user
    
    grouped <- FALSE
    if (is.null(name)) {
        if (is.null(colnames(b$stats)))
            nams <- paste("Sample",1:ncol(b$stats),sep=" ")
        else
            nams <- colnames(b$stats)
        name <- list(Samples=nams)
    }
    else if (length(name)==1 && name=="none") {
        nams <- rep("",ncol(b$stats))
        name <- list(Samples=nams)
    }
    else if (is.list(name)) { # Is sample.list
        nams <- unlist(name,use.names=FALSE)
        grouped <- TRUE
    }
    cols <- getColorScheme()
    
    # Data series
    d <- as.data.frame(round(b$stat,3))
    ids <- 0:(ncol(d)-1)
    d <- rbind(ids,d)
    colnames(d) <- nams
    names(ids) <- colnames(d)
    counter <- 0
    series <- vector("list",length(name))
    names(series) <- names(name)
    for (n in names(series)) {
        counter <- counter + 1
        series[[n]] <- list()
        series[[n]]$name=n
        if (grouped)
            series[[n]]$color=cols$fill[counter]
        else
            series[[n]]$color=cols$fill[1]
        m <- match(name[[n]],colnames(d))
        series[[n]]$data <- unname(as.list(d[,m]))
    }
    # Outlier series (if any)
    counter <- 0
    outliers <- vector("list",length(name))
    names(outliers) <- names(name)
    for (n in names(outliers)) {
        counter <- counter + 1
        outliers[[n]] <- list()
        outliers[[n]]$id <- n
        outliers[[n]]$name <- n
        outliers[[n]]$type <- "scatter"
        outliers[[n]]$showInLegend <- FALSE
        if (grouped) {
            outliers[[n]]$color <- cols$fill[counter]
            outliers[[n]]$marker <- list(
                fillColor=cols$fill[counter],
                symbol="circle",
                lineWidth=1,
                lineColor=cols$border[counter]
            )
        }
        else {
            outliers[[n]]$color <- cols$fill[1]
            outliers[[n]]$marker <- list(
                fillColor=cols$fill[1],
                symbol="circle",
                lineWidth=1,
                lineColor=cols$border[1]
            )
        }
        outliers[[n]]$data <- list()
        m <- match(name[[n]],colnames(d))
        if (length(m)>0) {
            for (i in m)
                x <- rep(d[1,i],length(o.list[[i]]))
                names(x) <- names(o.list[[i]])
                outliers[[n]]$data <- c(outliers[[n]]$data,
                    make.highcharts.points(x,o.list[[i]],unname(altnames)))
        }
    }
        
    # Boxplot tooltip point formatter for the case of zeros
    tooltip.point.formatter <- paste("function() {",
        "   var min = this.low === 0.001 ? 0 : this.low;" ,
        "   var q1 = this.q1 === 0.001 ? 0 : this.q1;" ,
        "   var med = this.median === 0.001 ? 0 : this.median;",
        "   var q3 = this.q3 === 0.001 ? 0 : this.q3;",
        "   var max = this.high === 0.001 ? 0 : this.high;",
        "   var str = 'Maximum: ' + max + '<br/>' +",
        "       'Upper quartile: ' + q3 + '<br/>' +",
        "       'Median: ' + med + '<br/>' +",
        "       'Lower quartile: ' + q1 + '<br/>' +",
        "       'Minimum: ' + min + '<br/>';",
        "   return  str;",
        "}",sep="")
    
    # Legend clicker
    boxplot.onclick <- paste("function() {",
        "   var chart =  this.chart;",
        "   var outlier_id =  chart.get(this.name);",
        "   if (!outlier_id.visible) {",
        "       outlier_id.show();",
        "   } else {",
        "       outlier_id.hide();",
        "   }",
        "}",sep="")
    
    if (is.null(obj$altnames)) {
        outlier.pointformat=paste(
            '<strong>Sample {point.category}</strong><br/>',
            'Gene ID: {point.name}<br/>',
            'Value: {point.y}<br/>',sep=""
        )
    }
    else {
        outlier.pointformat=paste(
            '<strong>Sample {point.category}</strong><br/>',
            'Gene ID: {point.name}<br/>',
            'Gene name: {point.alt_name}<br/>',
            'Value: {point.y}<br/>',sep=""
        )
    }
    
    json <- switch(jl,
        highcharts = {
            toJSON(
                list(
                    chart=list(
                        type="boxplot"
                    ),
                    title=list(
                        text=paste("Boxplot ",status,sep="")
                    ),
                    legend=list(
                        enabled=TRUE
                    ),
                    xAxis=list(
                        categories=nams,
                        title=list(
                            text="Sample name",
                            margin=25,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontWeight="bold"
                            )
                        )
                    ),
                    yAxis=list(
                        title=list(
                            useHTML=TRUE,
                            text="Read count (log<sub>2</sub>)",
                            margin=30,
                            style=list(
                                color="#000000",
                                fontSize="1.1em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        )
                    ),
                    plotOptions=list(
                        boxplot=list(
                            fillColor="#F0F0E0",
                            lineWidth=2,
                            medianColor="#000000",
                            medianWidth=3,
                            stemColor="#000000",
                            stemDashStyle="dash",
                            stemWidth=1,
                            whiskerColor="#000000",
                            whiskerLength="75%",
                            whiskerWidth=1,
                            grouping=FALSE,
                            tooltip=list(
                                headerFormat=paste(
                                    '<span style="font-size:1.1em;',
                                    'color:{series.color};',
                                    'font-weight:bold">',
                                    '\u25CF </span>',
                                    '<span style="font-size:1.1em;',
                                    'font-weight:bold">',
                                    'Condition {series.name}</span><br/>',
                                    '<span style="font-weight:bold">',
                                    'Sample {point.key}</span><br/>',sep=""
                                ),
                                pointFormatter=tooltip.point.formatter
                            ),
                            events=list(
                                legendItemClick=boxplot.onclick
                            )
                        ),
                        scatter=list(
                            allowPointSelect=TRUE,
                            tooltip=list(
                                headerFormat=paste(
                                    '<span style="font-weight:bold;',
                                    'color:{series.color};">',
                                    '\u25CF </span>',
                                    '<span style="font-weight:bold">',
                                    'Condition {series.name}</span><br/>',
                                    sep=""
                                ),
                                pointFormat=outlier.pointformat
                            ),
                            states=list(
                                hover=list(
                                    marker=list(
                                        enabled=FALSE
                                    )
                                )
                            )
                            
                        )
                    ),
                    series=c(unname(series),unname(outliers))
                )
            )
        }
    )
    return(unquote_js_fun(json))
}

#' Volcano JSON exporter for the metaseqR package
#'
#' Non-exportable JSON exporter for \code{\link{diagplot.volcano}}.
#'
#' @param obj A list holding volcano plot data. See \code{\link{diaplot.volcano}}.
#' @param jl JavaScript charting library to export. Currently only \code{"highcharts"}
#' supported.
#' @return A JSON string.
#' @author Panagiotis Moulos
volcanoToJSON <- function(obj,jl=c("highcharts")) {
    jl <- tolower(jl[1])
    f <- obj$x
    p <- obj$y
    xlim <- obj$xlim
    ylim <- obj$ylim
    alt.names <- obj$altnames
    pcut <- obj$pcut
    fcut <- obj$fcut
    up <- obj$user$up
    down <- obj$user$down
    ff <- obj$user$unf
    pp <- obj$user$unp
    alt.names.neutral <- obj$user$ualt
    con <- obj$user$con
    
    switch(jl,
        highcharts = {
            if (is.null(alt.names))
                point.format=paste("<strong>Gene ID: </strong>{point.name}<br>",
                    "<strong>Fold change: </strong>{point.x}<br>",
                    "<strong>Significance: </strong>{point.y}",sep="")
            else
                point.format=paste("<strong>Gene name: </strong>",
                    "{point.alt_name}<br>",
                    "<strong>Gene ID: </strong>{point.name}<br>",
                    "<strong>Fold change: </strong>{point.x}<br>",
                    "<strong>Significance: </strong>{point.y}",sep="")
                json <- toJSON(
                    list(
                        chart=list(
                        type="scatter",
                        zoomType="xy"
                    ),
                    title=list(
                        text=paste("Volcano plot for",con)
                    ),
                    xAxis=list(
                        title=list(
                            text="Fold change",
                            margin=20,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        min=round(xlim[1],3),
                        max=round(xlim[2],3)
                    ),
                    yAxis=list(
                        title=list(
                            useHTML=TRUE,
                            text="Significance (-log<sub>10</sub>(p-value))",
                            margin=25,
                            style=list(
                                color="#000000",
                                fontSize="1.2em"
                            )
                        ),
                        labels=list(
                            style=list(
                                color="#000000",
                                fontSize="1.1em",
                                fontWeight="bold"
                            )
                        ),
                        startOnTick=TRUE,
                        endOnTick=TRUE,
                        showLastLabel=TRUE,
                        gridLineWidth=1,
                        min=round(ylim[1]-2,3),
                        max=round(ylim[2],3)
                    ),
                    #legend=list(
                    #    layout="vertical",
                    #    align="left",
                    #    verticalAlign="top",
                    #    floating=TRUE,
                    #    backgroundColor="#FFFFFF",
                    #    borderWidth=1
                    #),
                    plotOptions=list(
                        scatter=list(
                            allowPointSelect=TRUE,
                            marker=list(
                                radius=2,
                                states=list(
                                    hover=list(
                                        enabled=TRUE,
                                        lineColor="#333333"
                                    )
                                )
                            ),
                            states=list(
                                hover=list(
                                    marker=list(
                                        enabled=FALSE
                                    )
                                )
                            ),
                            tooltip=list(
                                headerFormat=paste("<span style=",
                                    "\"font-size:1.1em;color:{series.color};",
                                    "font-weight:bold\">{series.name}<br>",
                                    sep=""),
                                pointFormat=point.format
                            ),
                            turboThreshold=50000
                        )
                    ),
                    series=list(
                        list(
                            name="Up-regulated",
                            color="#EE0000",
                            marker=list(
                                symbol="circle"
                            ),
                            data=make.highcharts.points(f[up],-log10(p[up]),
                                unname(alt.names[up]))
                        ),
                        list(
                            name="Down-regulated",
                            marker=list(
                                symbol="circle"
                            ),
                            color="#00CD00",
                            data=make.highcharts.points(f[down],-log10(p[down]),
                                unname(alt.names[down]))
                        ),
                        list(
                            name="Unregulated",
                            marker=list(
                                symbol="circle"
                            ),
                            color="#0000EE",
                            data=make.highcharts.points(ff,-log10(pp),
                                unname(alt.names.neutral))
                        ),
                        list(
                            name="Downfold threshold",
                            color="#000000",
                            type="line",
                            dashStyle="dash",
                            marker=list(
                                enabled=FALSE
                            ),
                            tooltip=list(
                                headerFormat=paste('<strong>{series.name}',
                                    '</strong><br/>',sep=""),
                                pointFormat=paste('<strong>Threshold: ',
                                    '</strong>{point.x}<br/>',sep="")
                            ),
                            data=list(round(c(-fcut,ylim[1]-5),3),
                                round(c(-fcut,ylim[2]),3))
                        ),
                        list(
                            name="Upfold threshold",
                            color="#000000",
                            type="line",
                            dashStyle="Dash",
                            marker=list(
                                enabled=FALSE
                            ),
                            tooltip=list(
                                headerFormat=paste('<strong>{series.name}',
                                    '</strong><br/>',sep=""),
                                pointFormat=paste('<strong>Threshold: ',
                                    '</strong>{point.x}<br/>',sep="")
                            ),
                            data=list(round(c(fcut,ylim[1]-5),3),
                                round(c(fcut,ylim[2]),3))
                        ),
                        list(
                            name="Significance threshold",
                            color="#000000",
                            type="line",
                            dashStyle="DashDot",
                            marker=list(
                                enabled=FALSE
                            ),
                            tooltip=list(
                                headerFormat=paste('<strong>{series.name}',
                                    '</strong><br/>',sep=""),
                                pointFormat=paste('<strong>Threshold: ',
                                    '</strong>{point.y}<br/>',sep="")
                            ),
                            data=list(round(c(xlim[1],-log10(pcut)),3),
                                round(c(xlim[2],-log10(pcut)),3))
                        )
                    )
                )
            )
        }
    )
    return(json)
}

unquote_js_fun <- function(js) {
    if (is.list(js))
        js <- lapply(js,unquote_js_fun)
    else {
        op <- gregexpr(pattern="function",js)
        cl <- gregexpr(pattern="}\\\"",js)
        if (length(op)>0) {
            starts <- as.numeric(op[[1]])
            for (i in 1:length(starts))
                substr(js,starts[i]-1,starts[i]-1) <- " "
            ends <- as.numeric(cl[[1]])
            for (i in 1:length(starts))
                substr(js,ends[i]+1,ends[i]+1) <- " "
        }
    }
    return(js)
}

getColorScheme <- function(n=NULL) {
    if (missing(n) || is.null(n))
        return(getColors())
    else {
        cols <- getColors()
        if (n > length(cols$fill)) {
            cols$fill <- rep(cols$fill,length.out=n)
            cols$border <- rep(cols$border,length.out=n)
            cols$select <- rep(cols$select,length.out=n)
            cols$trans <- rep(cols$trans,length.out=n)
        }
        return(cols)
    }
}

getColors <- function() {
    return(list(
        fill=c("#CD0000","#00CD00","#0000EE","#FFD700","#87CEEB","#CD8500",
            "#DEB887","#FF0000","#0000FF","#00FF00","#FFA500","#A9A9A9",
            "#008B00","#313131","#FFC0CB","#A52A2A","#FF00FF","#9ACD32",
            "#8B636C","#2E8B57","#008B8B"),
        border=c("#850000","#006B00","#000085","#927C00","#156280","#5A3A00",
            "#8B7457","#935E18","#000080","#008500","#603E00","#454545",
            "#073E07","#000000","#896067","#691111","#7C007C","#3A4D14",
            "#5B1726","#0C2517","#062A2A"),
        selected=c("#FF0000","#00FF00","#0066FF","#FFD700","#FFEB77","#FFB428",
            "#FFD9A5","#FF326D","#0089FF","#B3FF00","#FFC352","#D9D9D9",
            "#00EC00","#8E8E8E","#FFDAE0","#F94444","#FF87FF","#C2FF45",
            "#EA889D","#4EE590","#00DADA"),
        trans=c("rgba(205,0,0,0.6)","rgba(0,205,0,0.6)","rgba(0,0,238,0.6)",
            "rgba(255,215,0,0.6)","rgba(135,206,235,0.6)","rgba(205,133,0,0.6)",
            "rgba(222,184,135,0.6)","rbga(255,0,0,0.5)","rgba(0,0,255,0.5)",
            "rgba(0,255,0,0.5)","rgba(255,165,0,0.6)","rgba(169,169,169,0.5)",
            "rgba(0,139,0,0.6)","rgba(49,49,49,0.6)","rgba(255,192,203,0.5)",
            "rgba(165,42,42,0.6)","rgba(255,0,255,0.6)","rgba(154,205,50,0.6)",
            "rgba(139,99,108,0.6)","rgba(46,139,87,0.6)","rgba(0,139,139,0.6)")
    ))
}
