## Authors: Gautier Stoll, Hélène Fohrer-Ting, Estelle Devêvre, Sarah LEVESQUE, Julie LE NAOUR, Juliette PAILLET, Jonathan POL
## 2019, INSERM U1138
## Version 0.10.5

##tmpIsV3p6 = (as.integer(strsplit(strsplit(version$version.string,split=" ")[[1]][3],split=".",fixed=TRUE)[[1]][1]) >= 3) & (as.integer(strsplit(strsplit(version$version.string,split=" ")[[1]][3],split=".",fixed=TRUE)[[1]][2]) >= 6) ## for testing R version

tmpIsV3p6 = TRUE

##library(flowCore)
##library(flowDensity)
##library(flowWorkspace)
##library(flowAI)
##library(FlowSOM)
##library(FlowSOMworkshop)
##library(Rtsne)
##library("pheatmap")

##library("data.table")

##library(reshape)
##library(reshape2)

##library(ggplot2)

##library(ggpubr)

#library(plyr)
#library(readr)
#library(dplyr)
#library(beeswarm)

#library(gridExtra)
#library(gplots)
#library(dunn.test)
if (tmpIsV3p6) {
    ##library(CytoML)
    }
## Internal tool: PlotStar with bigger legend
PlotStarsBigLeg <- function(fsom,
                      markers=fsom$map$colsUsed,
                      view="MST", #"grid","tSNE"
                      colorPalette=grDevices::colorRampPalette(
                        c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                          "yellow", "#FF7F00", "red", "#7F0000")),
                      starBg = "white",
                      backgroundValues = NULL,
                      backgroundColor = function(n){
                        grDevices::rainbow(n, alpha=0.3)},
                      backgroundLim = NULL,
                      backgroundBreaks = NULL,
                      backgroundSize = NULL,
                      thresholds=NULL,
                      legend=TRUE,
                      query=NULL,
                      main=""){
    # Add star chart option to iGraph
    igraph::add.vertex.shape("star", clip=igraph::igraph.shape.noclip, plot=mystarBL,
                    parameters=list(vertex.data=NULL,vertex.cP = colorPalette,
                                    vertex.scale=TRUE, vertex.bg = starBg))

    if(is.null(thresholds)){
        # Use MFIs
        data <- fsom$map$medianValues[, markers,drop=FALSE]
        scale <- TRUE
    } else {
        # scale thresholds same as data
        if(fsom$transform){
            warning("Thresholds should be given in the transformed space")
        }
        if(!is.null(fsom$scaled.center)){
          thresholds <- scale(t(thresholds),
                              center = fsom$scaled.center[markers],
                              scale = fsom$scaled.scale[markers])
        }
        # Use pctgs of cells above threshold as star plot values
        data <-
            t(sapply(seq_len(fsom$map$nNodes), function(i) {
                res = NULL
                for(m in seq_along(markers)){
                    res = c(res,
                            sum(subset(fsom$data,
                               fsom$map$mapping[,1] == i)[,
                                                  markers[m]] > thresholds[m])/
                            sum(fsom$map$mapping[,1] == i))

                }
                res
            }))
        scale <- FALSE
    }

    # Choose layout type
    switch(view,
        MST  = { layout <- fsom$MST$l
                    lty <- 1},
        grid = { layout <- as.matrix(fsom$map$grid)
                    lty <- 0},
        tSNE = { layout <- fsom$MST$l2
                    lty <- 0},
        stop("The view should be MST, grid or tSNE. tSNE will only work
                   if you specified this when building the MST.")
    )

    # Choose background colour
    if (!is.null(backgroundValues)) {
        background <- computeBackgroundColorBL(backgroundValues,backgroundColor,
                                             backgroundLim, backgroundBreaks)
        if (is.null(backgroundSize)) {
          backgroundSize <- fsom$MST$size
          backgroundSize[backgroundSize == 0] <- 3
        }
    } else {
        background <- NULL
    }

    # Save plot window settings and minimize margin
    oldpar <- graphics::par(no.readonly = TRUE)
    graphics::par(mar=c(1,1,1,1))

    # Add legend
    if(legend){
        if(!is.null(backgroundValues)){
            # Make plot with folowing layout
            # 1 3
            # 2 3
            graphics::layout(matrix(c(1,1,3,3,2,2), 3, 2, byrow = TRUE),
                    widths=c(1), heights=c(1,3,1))
        } else {
            graphics::layout(matrix(c(1,2), 1, 2, byrow = TRUE),
                   widths=c(1,2), heights=c(1))
        }

       if(is.null(query)){
            plotStarLegendBL(fsom$prettyColnames[markers],
                            colorPalette(ncol(data)))
        } else {
            plotStarQuery(fsom$prettyColnames[markers],
                            values=query == "high",
                            colorPalette(ncol(data)))
        }

        if(!is.null(backgroundValues)){
            PlotBackgroundLegendBL(backgroundValues,background)
        }
    }

    # Plot the actual graph
    igraph::plot.igraph(fsom$MST$g,
                        vertex.shape = "star",
                        vertex.label = NA,
                        vertex.size = fsom$MST$size,
                        vertex.data = data,
                        vertex.cP = colorPalette(ncol(data)),
                        vertex.scale = scale,
                        layout = layout,
                        edge.lty = lty,
                        mark.groups = background$groups,
                        mark.col = background$col[background$values],
                        mark.border = background$col[background$values],
                        mark.expand	= backgroundSize,
                        main=main
    )
    # Reset plot window
    graphics::par(oldpar)
    graphics::layout(1)
}
## Internal tool, for BigLegendPlot
computeBackgroundColorBL <- function(backgroundValues,
                                    backgroundColor,
                                    backgroundLim = NULL,
                                    backgroundBreaks = NULL){
    # Choose background colour
    backgroundList <- list()
    backgroundColors <- NULL
    if(!is.null(backgroundValues)){
        if(is.numeric(backgroundValues)){
            backgroundList <- as.list(seq_along(backgroundValues))

            if(class(backgroundColor)=="function" &
               !is.null(backgroundBreaks) &
               length(backgroundBreaks)>1)
            {
                backgroundColor <- backgroundColor(length(backgroundBreaks))
            } else if (class(backgroundColor)=="function"){
                backgroundColor <- backgroundColor(100)
                backgroundBreaks <- length(backgroundColor)
            } else if (is.null(backgroundBreaks)){
                backgroundBreaks <- length(backgroundColor)
            }

            if(length(backgroundLim) > 0){
                ids <- cut(c(backgroundLim,backgroundValues),
                           backgroundBreaks
                )[-c(seq_along(backgroundLim))]
            } else {
                ids <- cut(backgroundValues,
                           backgroundBreaks)
            }
            backgroundValues <- ids
#             backgroundColors <- backgroundColor[ids]
        } else {
            if(! is.factor(backgroundValues)){
                backgroundValues <- as.factor(backgroundValues)
            }

            backgroundList <- as.list(seq_along(backgroundValues))

            if(class(backgroundColor)=="function"){
                backgroundColor <- backgroundColor(
                    length(levels(backgroundValues)))
            }

            if(length(backgroundColor) < length(levels(backgroundValues))){
                stop("You specified less backgroundcolours than groups.")
            }

        }
    }
    backgroundColors <- backgroundColor[backgroundValues]

    list(values=backgroundValues,
         col=backgroundColor,
         groups=backgroundList)
}
## Internal tool, for BL
plotStarLegendBL <- function(labels, colors=grDevices::rainbow(length(labels)),
                            main=""){
    graphics::plot(1, type="n", xlab="", ylab="",
        xlim=c(-10, 10), ylim=c(-3, 3),asp=1,
        bty="n",xaxt="n",yaxt="n",main=main)

    graphics::stars(matrix(c(1:(2*length(labels))),nrow=2),col.segments=colors,
        locations=c(0,0),draw.segments = TRUE,add=TRUE,
        inches=FALSE)
    n <- length(labels)
    angle <- 2*pi / n
    angles <- seq(angle/2,2*pi,by=angle)

    left <- (angles > (pi/2) & angles < (3*pi/2))
    x <- c(2,-2)[left+1]
    y_tmp <- c(seq(-2,2,by= 4/(sum(!left)+1))[-c(1,sum(!left)+2)],
                seq(2,-2,by=-4/(sum(left)+1))[-c(1,sum(left)+2)])
    y <- shiftFunctionBL(y_tmp,max((cummax(y_tmp)<0)*seq_along(y_tmp)))

    for(i in seq_along(labels)){
        graphics::text(x= x[i],
            y= y[i],
            labels=labels[i],
            adj = c(as.numeric(left)[i],0.5),
            cex = 0.5)

        graphics::lines(x=c(x[i]+c(-0.2,0.2)[left[i]+1],
                c(1.5,-1.5)[left[i]+1],
                cos(angles[i])),
            y=c(y[i],
                y[i],
                sin(angles[i])),
            col=colors[i],
            lwd=2)
    }
}

##Internal tool, for Big Legend
shiftFunctionBL <- function(x,n){
    c(x[(n+1):length(x)],x[1:n])
}
##Internal tool, for Big Legend
PlotBackgroundLegendBL <- function(backgroundValues, background,
                                 main="Background"){
    graphics::plot.new()
    if(is.numeric(backgroundValues)) {
        legendContinuous(background$col,
                         as.numeric(gsub(".*,","",
                                         gsub("].*","",
                                              levels(background$values)))))
    } else {
        graphics::legend("center", legend=levels(background$values),
               fill=background$col,
               cex=0.7,
               ncol =  ceiling(length(levels(background$values)) / 10),
               bty="n",
               title=main)
    }
}
##Internal tool, for Big Legend
mystarBL <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
    }
    vertex.size    <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    data <- params("vertex", "data")
    cP <- params("vertex","cP")
    scale <- params("vertex","scale")
    bg <- params("vertex","bg")
    graphics::symbols(coords[, 1], coords[, 2], circles = vertex.size,
                      inches = FALSE, bg = bg, bty='n', add=TRUE)
    graphics::stars(data, locations = coords, labels = NULL,scale=scale,
            len = vertex.size, col.segments = cP,
            draw.segments = TRUE, mar = c(0, 0, 0, 0), add=TRUE,
            inches=FALSE)

}

## Internal tool: gating subset from FlowSOMworshop

gating_subset_toolBox<- function(flowjo_res, gate){

  if(!is.null(flowjo_res$flowSet)){
    res <- lapply(seq_len(length(flowjo_res$flowSet)),
                  function(i){
                    flowjo_res$flowSet[[i]][flowjo_res$gates[[i]][,gate],]
                  })
    names(res) <- flowCore::sampleNames(flowjo_res$flowSet)
    return(list(
      flowSet = flowCore::flowSet(res),
      gates = lapply(flowjo_res$gates, function(x)x[x[,gate], ])
    ))

  } else {
    return(list(flowFrame = flowjo_res$flowFrame[flowjo_res$gates[,gate], ],
                gates = flowjo_res$gates[flowjo_res$gates[,gate], ]))
  }
}

## Internal tool: new parse_flowjo for CytoML 1.12.0

parse_flowjo_CytoML_v12 <- function (files, wsp_file, group = "All Samples")
{
  wsp <- CytoML::open_flowjo_xml(wsp_file)
  o <- capture.output(gates <- suppressMessages(CytoML::flowjo_to_gatingset(wsp,
                                                                       group)))
  files_in_wsp <- gates@data@origSampleVector
  counts <- as.numeric(gsub(".*_([0-9]*)$", "\\1", files_in_wsp))
  result <- list()
  for (file in files) {
    print(paste0("Processing ", file))
    file_id <- grep(gsub(".*/", "", file), files_in_wsp)
    if (length(file_id) == 0) {
      stop("File not found. Files available: ", gsub("_[0-9]*$",
                                                   "\n", files_in_wsp))
    }
    gate_names <- flowWorkspace::gs_get_pop_paths(gates[[file_id]],path = "auto")
    gatingMatrix <- matrix(FALSE, nrow = counts[file_id],
      ncol = length(gate_names), dimnames = list(NULL,
      gate_names))
    for (gate in gate_names) {
      gatingMatrix[, gate] <- flowWorkspace::gh_pop_get_indices(gates[[file_id]],gate)
    }
if ((unlist(packageVersion("flowWorkspace"))[1] == 3) && (unlist(packageVersion("flowWorkspace"))[2] <= (32)))     
{ff <- flowWorkspace::getData(gates[[file_id]], "root")}
    else {ff <- flowWorkspace::gh_pop_get_data(gates[[file_id]], "root")}

    ff@exprs[, "Time"] <- ff@exprs[, "Time"] * 100
    result[[file]] <- list(flowFrame = ff, gates = gatingMatrix)
  }
  if (length(files) == 1) {
    result <- result[[1]]
  }
  else {
    result <- list(flowSet = flowCore::flowSet(lapply(result,
                function(x) x$flowFrame)), gates = lapply(result,
                     function(x) x$gates))
  }
  return(result)
}
## Internal tool: corrected parse_flowjo

    parse_flowjo_CytoML <- function (files, wsp_file, group = "All Samples", plot = FALSE)
    {
       ##wsp <- flowWorkspace::openWorkspace(wsp_file)
      wsp <- CytoML::openWorkspace(wsp_file)
        o <- capture.output(gates <- suppressMessages(CytoML::parseWorkspace(wsp,
                                                                             group)))
        files_in_wsp <- gates@data@origSampleVector
        counts <- as.numeric(gsub(".*_([0-9]*)$", "\\1", files_in_wsp))
        result <- list()
        for (file in files) {
            print(paste0("Processing ", file))
            file_id <- grep(gsub(".*/", "", file), files_in_wsp)
            if (length(file_id) == 0) {
                stop("File not found. Files available: ", gsub("_[0-9]*$",
                    "\n", files_in_wsp))
            }
            gate_names <- flowWorkspace::getNodes(gates[[file_id]],
                path = "auto")
            gatingMatrix <- matrix(FALSE, nrow = counts[file_id],
                ncol = length(gate_names), dimnames = list(NULL,
                    gate_names))
            for (gate in gate_names) {
                gatingMatrix[, gate] <- flowWorkspace::getIndiceMat(gates[[file_id]],
                    gate)
            }
            ff <- flowWorkspace::getData(gates[[file_id]], "root")
            ff@exprs[, "Time"] <- ff@exprs[, "Time"] * 100
            result[[file]] <- list(flowFrame = ff, gates = gatingMatrix)
            if (plot) {
                flowWorkspace::plot(gates[[file_id]])
            }
        }
        if (length(files) == 1) {
            result <- result[[1]]
        }
        else {
            result <- list(flowSet = flowCore::flowSet(lapply(result,
                function(x) x$flowFrame)), gates = lapply(result,
                                                          function(x) x$gates))
        }
        return(result)
    }

## Internal tool: GetClusters ?
##if(!exists("GetClusters",mode="function")) {
    GetClusters <- function(fsom) {
      if (class(fsom) == "list" & !is.null(fsom$FlowSOM)) {
        fsom <- fsom$FlowSOM
      }
      if (class(fsom) != "FlowSOM") {
        stop("fsom should be a FlowSOM object.")
      }
      return(fsom$map$mapping[,1])
    }
##}

## Internal tool: GetMetClusters ?
##if(!exists("GetMetaclusters",mode="function")) {
    GetMetaclusters <- function(fsom, meta = NULL){

      if (class(fsom) == "list" & !is.null(fsom$FlowSOM)) {
        if (is.null(meta) & !is.null(fsom$metaclustering)) {
          meta <- fsom$metaclustering
        }
        fsom <- fsom$FlowSOM
      }
      if (class(fsom) != "FlowSOM"){
        stop("fsom should be a FlowSOM object.")
      }
      if(is.null(meta)){
        stop("No metaclustering found.")
      }

      return(meta[fsom$map$mapping[,1]])
    }
##}

    ## Internal tool: Seems that PlotLabels diseapear...
##if(!exists("PlotLabels",mode="function")) {

    PlotLabels <- function(fsom,
                       labels,
                       view="MST",
                       main=NULL,
                       nodeSize=fsom$MST$size,
                       fontSize = 1,
                       backgroundValues = NULL,
                       backgroundColor = function(n){
                         grDevices::rainbow(n,alpha=0.3)},
                       backgroundLim = NULL,
                       backgroundBreaks = NULL){
        switch(view,
               MST  = { layout <- fsom$MST$l
                   lty <- 1},
               grid = { layout <- as.matrix(fsom$map$grid)
                   lty <- 0},
               tSNE = { layout <- fsom$MST$l2
                   lty <- 0},
               stop("The view should be MST, grid or tSNE. tSNE will only work
                   if you specified this when building the MST.")
               )

                                        # Choose background colour
        if(!is.null(backgroundValues)){
            background <- computeBackgroundColor(backgroundValues,backgroundColor,
                                                 backgroundLim, backgroundBreaks)
        } else {
            background <- NULL
        }

        igraph::plot.igraph(fsom$MST$graph,
                            layout=layout,
                            vertex.size=nodeSize,
                            vertex.label=labels,
                            vertex.label.cex = fontSize,
                            edge.lty=lty,
                            mark.groups=background$groups,
                            mark.col=background$col[background$values],
                            mark.border=background$col[background$values],
                            main=main)

    }
##}


## Internal tool: extract meta-clusters count ratio in percent
get_pctgsMT <- function(fSOM,metacl, meta_names = NULL){
  `%>%` <- magrittr::`%>%`
  cell_ids <- fSOM$metaData
  files <- sapply(seq_len(length(cell_ids)),
                  function(i){
                    rep(gsub(".*/", "", names(cell_ids)[i]),
                        cell_ids[[i]][2] - cell_ids[[i]][1] + 1)
                  }) %>%
    unlist()
  pctgs <- table(files, GetClusters(fSOM)) %>%
    as.matrix() %>%
    apply(1, function(x){x/sum(x)}) %>%
    t()
  pctgs_meta <- table(files, GetMetaclusters(fSOM,meta = metacl)) %>%
    as.matrix() %>%
    apply(1, function(x){x/sum(x)}) %>%
    t()
  if(!is.null(meta_names)) colnames(pctgs_meta) <- meta_names
  return(list("pctgs" = as.matrix(pctgs),
              "pctgs_meta" = as.matrix(pctgs_meta)))
}

## Internal tool: extract absolute count of meta-clusters
get_abstgsMT <- function(fSOM,metacl, meta_names = NULL){
  `%>%` <- magrittr::`%>%`
  cell_ids <- fSOM$metaData
  files <- sapply(seq_len(length(cell_ids)),
                  function(i){
                    rep(gsub(".*/", "", names(cell_ids)[i]),
                        cell_ids[[i]][2] - cell_ids[[i]][1] + 1)
                  }) %>%
      unlist()
  pctgs <- table(files, GetClusters(fSOM)) %>%
      as.matrix() %>% apply(1,as.numeric) %>% t()
  pctgs_meta <- table(files, GetMetaclusters(fSOM,meta = metacl)) %>%
      as.matrix()
  ## %>% apply(1,as.numeric) %>% t()
  if(!is.null(meta_names)) colnames(pctgs_meta) <- meta_names
  return(list("abstgs" = pctgs,
              "abstgs_meta" = pctgs_meta))
}

##Internal tool: return p-value of Tukey test, given metacluster names
TukeyTestSarah = function(fSOMTable, metaClust){TukeyHSD(aov(as.formula(paste(metaClust,"~ Treatment")),data=fSOMTable))$Treatment[,4]}

##Internal tool every kind of boxplots-heatmaps
BoxPlotMetaClustFull <- function(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab,Norm=FALSE,Marker="",Robust,ClustHeat)
{
    ## Search for the marker
    treatmentTable$Treatment=gsub(" ","_",treatmentTable$Treatment,fixed=T)
    ControlTreatment = gsub(" ","_",ControlTreatment,fixed=T)
    MarkerIndex = which(TreeMetaCl$fSOMTree$prettyColnames == Marker)
    if(length(MarkerIndex) == 1) {
        print(paste("User marker:",Marker,":",names(TreeMetaCl$fSOMTree$prettyColnames)[MarkerIndex]))
        fSOMnbrs = sapply(unique(TreeMetaCl$metaCl),function(metaClust){
            clusterList=which(TreeMetaCl$metaCl == metaClust)
            metaClustIndices=unlist(sapply(clusterList,function(cluster){which(TreeMetaCl$fSOMTree$map$mapping[,1] == cluster)}))
            sapply(TreeMetaCl$fSOMTree$metaData,function(StartEnd){
                indices = intersect((StartEnd[1]:StartEnd[2]),metaClustIndices)
                return(median(TreeMetaCl$fSOMTree$data[indices,MarkerIndex]))
            })
        })
        colnames(fSOMnbrs)=unique(TreeMetaCl$metaCl)
        row.names(fSOMnbrs)=gsub(".*/","",row.names(fSOMnbrs))
        PlotLab=Marker
    }
    else {
    ## constuct fSOMnbrs, according to Norm (false: percentage, true: normalized)

    if (Norm) {
        abstgs=get_abstgsMT(TreeMetaCl$fSOMTree,TreeMetaCl$metaCl)
        fSOMnbrs<-abstgs$abstgs_meta
        if (is.null(treatmentTable$NormalizationFactor)) {stop("No column NormalizationFactor in annotation table")}
        NormFactors = sapply(row.names(fSOMnbrs),function(fileFCS){treatmentTable$NormalizationFactor[which(treatmentTable$files == fileFCS)]})
        for (index in 1:length(NormFactors)) {fSOMnbrs[index,] = fSOMnbrs[index,]/NormFactors[index] }
        PlotLab=paste("size of ",yLab,sep="")
    }
    else {
        pctgs<-get_pctgsMT(TreeMetaCl$fSOMTree,TreeMetaCl$metaCl)
        fSOMnbrs<-pctgs$pctgs_meta
        fSOMnbrs<-fSOMnbrs*100
        PlotLab=paste("% of ",yLab,sep="")
    }
        ##colnames(fSOMnbrs)=unique(TreeMetaCl$metaCl)
        ##row.names(fSOMnbrs)=gsub(".*/","",row.names(fSOMnbrs))
        }
    treatmentsFSOM=sapply(row.names(fSOMnbrs),function(fileFCS){treatmentTable$Treatment[which(treatmentTable$files == fileFCS)]})
    Treatments=unique(treatmentTable$Treatment)
    if(length(which(Treatments == ControlTreatment)) == 0) {stop(paste("No",ControlTreatment,"in annotation table"))}
    treatmentsFSOM=factor(treatmentsFSOM,levels=c(ControlTreatment,setdiff(Treatments,ControlTreatment))) # set control treatment at first
    if(length(MarkerIndex) == 1) {
      pdf(file=gsub("/","_",paste(Title,"_BoxPlot",Marker,"Metacl.pdf",sep=""),fixed=T))
      } else {
    if (Norm) {pdf(file=paste(Title,"_BoxPlotNormMetacl.pdf",sep=""))}
    else {pdf(file=paste(Title,"_BoxPlotPercentMetacl.pdf",sep=""))}}
    metaclNumber=length(fSOMnbrs[1,])
    par(mfrow=c(6,6),las=2,mar=c(BottomMargin,3,1,.5),mgp=c(1.8,.8,0)) ## page have 6x6 boxplots
    fSOMnbrs=fSOMnbrs[,unique(unique(TreeMetaCl$metaCl))]
    cex4Title=exp(-max(sapply(colnames(fSOMnbrs),nchar))/50)
    for (metaCl in (1:metaclNumber)){ ## boxplots with no annotations
        plotDf=data.frame(PP=fSOMnbrs[,metaCl],TreatmentFSOM=treatmentsFSOM) ## dataframe for box plot
        boxplot(PP ~ TreatmentFSOM,data=plotDf,main=paste("mtcl",colnames(fSOMnbrs)[metaCl],sep="_"),xlab="",ylab=PlotLab,cex.axis=.5,cex.main=cex4Title,cex.lab=.5)
        beeswarm::beeswarm(PP ~ TreatmentFSOM,data=plotDf,main=paste("mtcl",colnames(fSOMnbrs)[metaCl],sep="_"),add=T,cex=.5,col="red")
    }
    par(mfrow=c(6,6),las=2,mar=c(BottomMargin,3,1,.5),mgp=c(1.8,.8,0))
    PvalPairwiseTable = sapply((1:metaclNumber),function(metaCl) ## construct pval table of tukey pairwise comparison test, boxplots with p-values annotation
    {
        plotDf=data.frame(PP=fSOMnbrs[,metaCl],TreatmentFSOM=treatmentsFSOM,noData=rep(1,length(treatmentsFSOM)))
        if (Robust) {
          invisible(capture.output(tmpNoData <- dunn.test::dunn.test(plotDf$noData,plotDf$TreatmentFSOM,table=F)))
          pairwisePval=tmpNoData$P
          names(pairwisePval) = gsub(" ","",tmpNoData$comparisons,fixed=T) ## create template
          if (length(unique(plotDf$TreatmentFSOM[which(is.finite(plotDf$PP))])) > 2 ){
          invisible(capture.output(tmp <- dunn.test::dunn.test(plotDf$PP,plotDf$TreatmentFSOM)
                                       ##tryCatch({dunn.test::dunn.test(plotDf$PP,plotDf$TreatmentFSOM)},
                  ##error = function(e){
                  ##  print(str(unique(treatmentsFSOM)))
                  ##  print("haha")
                  ##  data.frame(P=rep(1,length(plotDf$PP)*(length(plotDf$PP)-1)/2),comparisons = combn(unique(treatmentsFSOM),2,function(x){paste(x[1],x[2],sep=" - ")}))
                  ##  })
                  ))
            pairwiseSignPval=tmp$P
            names(pairwiseSignPval) = gsub(" ","",tmp$comparisons,fixed=T)
            for(name in names(pairwiseSignPval)){
              pairwisePval[[name]] = pairwiseSignPval[[name]]
              }}
            
        } else {
          pairwisePval=TukeyHSD(aov(noData ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM[,4]
          pairwisePval[]=1
          names(pairwisePval)=row.names(TukeyHSD(aov(noData ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM) ##create template
          if (length(unique(plotDf$TreatmentFSOM[which(is.finite(plotDf$PP))])) > 2 ){ 
          pairwiseSignPval=TukeyHSD(aov(PP ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM[,4]
          names(pairwiseSignPval)=row.names(TukeyHSD(aov(PP ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM)
          for(name in names(pairwiseSignPval)){
            pairwisePval[[name]] = pairwiseSignPval[[name]]
          }
          }
            
        }
        ListSignif=(sapply(1:length(pairwisePval),function(index){
            if (!is.finite(pairwisePval[index])){return(c())}
            if(pairwisePval[index] < 0.0001){return(c("****",strsplit(names(pairwisePval)[index],split = "-")[[1]]))}
            else if(pairwisePval[index] < 0.001){return(c("***",strsplit(names(pairwisePval)[index],split = "-")[[1]]))}
            else if(pairwisePval[index] < 0.01){return(c("**",strsplit(names(pairwisePval)[index],split = "-")[[1]]))}
            else if(pairwisePval[index] < 0.05){return(c("*",strsplit(names(pairwisePval)[index],split = "-")[[1]]))}
        }))
        ListSignif = ListSignif[which(sapply(ListSignif,length) > 0)]
        ListSignifPosIndex = lapply(ListSignif,function(hit){
            return(c(which(levels(plotDf$TreatmentFSOM) == hit[2]),which(levels(plotDf$TreatmentFSOM) == hit[3])))})
        minTr=min(plotDf$PP,na.rm=T)
        maxTr=max(plotDf$PP,na.rm=T)
        boxplot(PP ~ TreatmentFSOM,
                data=plotDf,main=paste("mtcl",colnames(fSOMnbrs)[metaCl],sep="_"),
                xlab="",
                ylab=PlotLab,
                cex.axis=.5,
                cex.main=cex4Title,
                cex.lab=.5,
                ylim=c(minTr,length(ListSignif)*abs(maxTr-minTr)*.2+maxTr)
                )
        if (length(ListSignif) > 0) {
            if (length(pairwisePval) > 1) ## more than one pair of comparison
                {
                    for (signifIndex in (1:length(ListSignif))) {
                        ##  print(maxTr+(signifIndex-.4)*abs(maxTr-minTr)*.2)
                        segments(y0=maxTr+(signifIndex-.4)*abs(maxTr-minTr)*.2,
                                 x0= ListSignifPosIndex[[signifIndex]][1],x1=ListSignifPosIndex[[signifIndex]][2])
                        text(x=(ListSignifPosIndex[[signifIndex]][1]+ListSignifPosIndex[[signifIndex]][2])/2,y=maxTr+(signifIndex-.1)*abs(maxTr-minTr)*.2,
                             labels=ListSignif[[signifIndex]][1])
                    }
                } else {
                    segments(y0=maxTr+(1-.4)*abs(maxTr-minTr)*.2,
                             x0= 1,x1=2)
                    text(x=1+1/2,y=maxTr+(1-.1)*abs(maxTr-minTr)*.2,
                             labels=ListSignif[1])
                    }
       }
        beeswarm::beeswarm(PP ~ TreatmentFSOM,data=plotDf,add=T,cex=.5,col="red")
       return(pairwisePval)
    })
    ## finish the construction of PvalTable, write csv files
    if(is.matrix(PvalPairwiseTable)) {
        PvalPairwiseTable=as.data.frame(PvalPairwiseTable)
    }
    else if (is.list(PvalPairwiseTable)) {PvalPairwiseTable = as.data.frame(do.call(cbind,PvalPairwiseTable))}
    else {PvalPairwiseTable=as.data.frame(t(PvalPairwiseTable))}
    names(PvalPairwiseTable)=paste("mtcl",colnames(fSOMnbrs)[1:metaclNumber],sep="_")
    par(mfrow=c(1,1),mar=c(3,2,3,1),cex=.5)

    if(length(MarkerIndex) == 1) {write.table(PvalPairwiseTable,gsub("/","_",paste(Title,"_PairwisePval",Marker,"Metacl.csv",sep=""),fixed=T),sep=";",col.names = NA)} else {
    if (Norm) {write.table(PvalPairwiseTable,paste(Title,"_PairwisePvalNormMetacl.csv",sep=""),sep=";",col.names = NA)}
    else {write.table(PvalPairwiseTable,paste(Title,"_PairwisePvalPercentMetacl.csv",sep=""),sep=";",col.names = NA)}}

    DF4lm = data.frame(y=c(fSOMnbrs),metaCl = c(sapply(colnames(fSOMnbrs),function(name){rep(name,length(fSOMnbrs[,1]))})),treat = rep(c(sapply(row.names(fSOMnbrs),function(name){as.character(treatmentTable$Treatment[which(treatmentTable$files == name)])})),length(fSOMnbrs[1,])))
    DF4lm$treat = factor(DF4lm$treat,levels=c(ControlTreatment,setdiff(unique(DF4lm$treat),ControlTreatment)))
    if (Robust) {
        pvalLmMatrix=as.matrix(PvalPairwiseTable)[which(sapply(row.names(PvalPairwiseTable),function(name){spName = strsplit(name,split = "-")[[1]];((spName[1] == ControlTreatment) |(spName[2] == ControlTreatment) )})),]
        row.names(pvalLmMatrix)=gsub(ControlTreatment,"",row.names(pvalLmMatrix),fixed=T)
        row.names(pvalLmMatrix)=gsub("-","",row.names(pvalLmMatrix),fixed=T)
    } else {
        pvalLmMatrix = t(do.call(rbind,lapply(colnames(fSOMnbrs),function(metaCl){
          SubDF4Lm = DF4lm[which(DF4lm$metaCl == metaCl),]
          SubDF4Lm$noData = rep(1,length(SubDF4Lm$y))
          retPval = summary(lm(noData ~ treat,data = SubDF4Lm))$coefficient[-1,4] ##create template
          if (length(unique(SubDF4Lm$treat[which(is.finite(SubDF4Lm$y))])) > 2) {
            SignRetPval = summary(lm(y ~ treat,data = DF4lm[which(DF4lm$metaCl == metaCl),]))$coefficient[-1,4]
            for(name in names(SignRetPval)){retPval[[name]] = SignRetPval[[name]]}
          }
          return(retPval)
        })))
        row.names(pvalLmMatrix) = setdiff(unique(DF4lm$treat),ControlTreatment)
    }
    colnames(pvalLmMatrix) = paste("mtcl",colnames(fSOMnbrs),sep= "_")
    pvalLmMatrix = rbind(rep(1,length(fSOMnbrs[1,])),pvalLmMatrix)
    row.names(pvalLmMatrix)[1] = ControlTreatment
    DF4lm$metaCl = factor(DF4lm$metaCl,level=unique(DF4lm$metaCl)) ## to get the right ordering after "by" function
    if (Robust) {
        meanMatrix  = t(by(DF4lm$y,list(DF4lm$metaCl,DF4lm$treat),function(x){median(x,na.rm=T)}))
       ## print(row.names(meanMatrix))
      ##  print(pvalLmMatrix)
        pvalLmMatrix=pvalLmMatrix[row.names(meanMatrix),]
    } else {
        meanMatrix  = t(by(DF4lm$y,list(DF4lm$metaCl,DF4lm$treat),function(x){mean(x,na.rm=T)})) }
    attr(meanMatrix,"class") = NULL
    attr(meanMatrix,"call") = NULL
    colnames(meanMatrix) = paste("mtcl",colnames(meanMatrix),sep= "_")
    pvalAnnotationMatrix = apply(pvalLmMatrix,c(1,2),function(x){
    if (!is.finite(x)){return("")}  
            else if (x < 0.0001){return("****")}
            else if(x < 0.001){return("***")}
            else if(x < 0.01){return("**")}
            else if(x < 0.05){return("*")} else {return("")}})
    if (length(MarkerIndex) == 1) {
        colMarginSize=20-18*exp(-max(sapply(colnames(meanMatrix[-1,]),nchar))/10)
        rowMarginSize=20-18*exp(-max(sapply(row.names(meanMatrix[-1,]),nchar))/10)
        colCex4Plot=exp(-max(sapply(colnames(meanMatrix[-1,]),nchar))/70)
        rowCex4Plot=exp(-max(sapply(row.names(meanMatrix[-1,]),nchar))/70)

        if (Robust) {heatTitle = paste("Median MFI of ",PlotLab,sep="")} else {heatTitle = paste("Mean MFI of ",PlotLab,sep="")}
        par(cex.main=exp(-nchar(heatTitle)/70))
        if (ClustHeat) {
            gplots::heatmap.2(meanMatrix,Rowv=F,Colv=T,dendrogram = "column",scale="none",col = heat.colors(100),cellnote = pvalAnnotationMatrix,
                      notecol = "black",trace = "none",cexRow = rowCex4Plot,cexCol=colCex4Plot,density.info="none",main=heatTitle,
                      notecex=.5,margins=c(colMarginSize,rowMarginSize),key.xlab = "",key.title="")
        }
        gplots::heatmap.2(meanMatrix,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = heat.colors(100),cellnote = pvalAnnotationMatrix,
                  notecol = "black",trace = "none",cexRow = rowCex4Plot,cexCol=colCex4Plot,density.info="none",main=heatTitle,
                  notecex=.5,margins=c(colMarginSize,rowMarginSize),key.xlab = "",key.title="")
    } else {
        colMarginSize=20-18*exp(-max(sapply(colnames(meanMatrix[-1,]),nchar))/10)
        rowMarginSize=20-18*exp(-max(sapply(row.names(meanMatrix[-1,]),nchar))/10)
        colCex4Plot=exp(-max(sapply(colnames(meanMatrix[-1,]),nchar))/70)
        rowCex4Plot=exp(-max(sapply(row.names(meanMatrix[-1,]),nchar))/70)
        if (Robust) { heatTitle = paste("Median ",PlotLab,sep="")}
        else {heatTitle = paste("Mean ",PlotLab,sep="")}
        par(cex.main=exp(-max(c(nchar(heatTitle),(nchar(ControlTreatment)+18)))/70))
        heatTitle=paste(heatTitle,"\n(rel. to ",ControlTreatment,", scaled)",sep="")
        meanMatrix=apply(meanMatrix,2,function(x){
          varCol = sd(x,na.rm=T)
          if (varCol == 0) {varCol = 1} ## contant column, no scale
          return((x-x[1])/varCol)
          })
        meanMatrix=meanMatrix[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")] ## get the correct ordering
        if (ClustHeat) {
          par(cex.main=.5)
            gplots::heatmap.2(meanMatrix[-1,],Rowv=F,Colv=T,dendrogram = "column",scale="none",col = gplots::bluered(100),cellnote = pvalAnnotationMatrix[-1,],
                      notecol = "black",trace = "none",cexRow = rowCex4Plot,cexCol=colCex4Plot,density.info="none",main=heatTitle,
                      ##distfun=function(x){dist(t(apply(meanMatrix,2,function(y){scale(y)})))},notecex=.5,margins=c(colMarginSize,rowMarginSize),key.xlab = "",key.title="")
                      notecex=.5,margins=c(colMarginSize,rowMarginSize),key.xlab = "",key.title="")
        }
            gplots::heatmap.2(meanMatrix[-1,],Rowv=F,Colv=F,dendrogram = "none",scale="none",col = gplots::bluered(100),cellnote = pvalAnnotationMatrix[-1,],
                      notecol = "black",trace = "none",cexRow = rowCex4Plot,cexCol=colCex4Plot,density.info="none",main=heatTitle,
                      distfun=function(x){dist(t(apply(meanMatrix,2,function(y){scale(y)})))},notecex=.5,margins=c(colMarginSize,rowMarginSize),key.xlab = "",key.title="")
    }
    par(cex.main=1)
    matrixPval4Heat=apply(as.matrix(PvalPairwiseTable)[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")],c(1,2),function(x){
      if (!is.finite(x)) {return(0)}
      else if (x < 0.0001) {return(-log10(0.0001))}
      else {return(-log10(x))}
    })
    matrixAnnot4Heat=apply(as.matrix(PvalPairwiseTable)[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")],c(1,2),function(x){
      if (!is.finite(x)){return("")}
      else if (x < 0.0001){return("****")}
      else if(x < 0.001){return("***")}
      else if(x < 0.01){return("**")}
      else if(x < 0.05){return("*")} 
      else if (x < 0.1){return(".")} else {return("")}})
    maxLogPval = max(unlist(matrixPval4Heat))
    colMarginSize=20-18*exp(-max(sapply(colnames(meanMatrix[-1,]),nchar))/10)
    rowMarginSize=20-18*exp(-max(sapply(row.names(meanMatrix[-1,]),nchar))/10)
    colCex4Plot=exp(-max(sapply(colnames(meanMatrix[-1,]),nchar))/70)
    rowCex4Plot=exp(-max(sapply(row.names(meanMatrix[-1,]),nchar))/70)
    if (Robust) {
        if (ClustHeat) {
            gplots::heatmap.2(matrixPval4Heat,Rowv=T,Colv=T,dendrogram = "both",scale="none",col = gray(1-((0:100)/100*maxLogPval/(-log10(0.0001)))),
                      trace="none",main="log10(Dunn p-values)",cexRow = rowCex4Plot,cexCol=colCex4Plot,margins=c(colMarginSize,rowMarginSize),density.info="none",
                      cellnote = matrixAnnot4Heat,notecol = "blue",notecex = .5,key.xlab = "",key.title="")
        }
            gplots::heatmap.2(matrixPval4Heat,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = gray(1-((0:100)/100*maxLogPval/(-log10(0.0001)))),
                      trace="none",cexRow = rowCex4Plot,cexCol=colCex4Plot,main="log10(Dunn p-values)",margins=c(colMarginSize,rowMarginSize),density.info="none",
                      cellnote = matrixAnnot4Heat,notecol = "blue",notecex = .5,key.xlab = "",key.title="")
    } else {
             if (ClustHeat) {
                 gplots::heatmap.2(matrixPval4Heat,Rowv=T,Colv=T,dendrogram = "both",scale="none",col = gray(1-((0:100)/100*maxLogPval/(-log10(0.0001)))),
                           trace="none",cexRow = rowCex4Plot,cexCol=colCex4Plot,main="log10(Tukey p-values)",margins=c(colMarginSize,rowMarginSize),density.info="none",
                           cellnote = matrixAnnot4Heat,notecol = "blue",notecex = .5,key.xlab = "",key.title="")
             }
             gplots::heatmap.2(matrixPval4Heat,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = gray(1-((0:100)/100*maxLogPval/(-log10(0.0001)))),
                       trace="none",cexRow = rowCex4Plot,cexCol=colCex4Plot,main="log10(Tukey p-values)",margins=c(colMarginSize,rowMarginSize),density.info="none",
                       cellnote = matrixAnnot4Heat,notecol = "blue",notecex = .5,key.xlab = "",key.title="") }
    dev.off()
    
    if (is.table(fSOMnbrs)) {
      DFSizes = as.data.frame(as.matrix.data.frame(fSOMnbrs))
      row.names(DFSizes) = row.names(fSOMnbrs)
      colnames(DFSizes) = colnames(fSOMnbrs)
    } else {DFSizes = as.data.frame(fSOMnbrs)}
    
    retData=list(DFSizes,as.data.frame(PvalPairwiseTable),as.data.frame(pvalLmMatrix))

    if(length(MarkerIndex) == 1) {write.table(pvalLmMatrix,gsub("/","_",paste(Title,"_LmPval",Marker,"Metacl.csv",sep=""),fixed=T),sep=";",col.names = NA)} else {
    if (Norm) {write.table(pvalLmMatrix,paste(Title,"_LmPvalNormMetacl.csv",sep=""),sep=";",col.names = NA)}
    else {write.table(pvalLmMatrix,paste(Title,"_LmPvalPercentMetacl.csv",sep=""),sep=";",col.names = NA)}}

    names(retData)=c("Sizes","PvalPairwise","PvalLm")
    return(retData)
}

## User tool: Plot Meta clusters labels
PlotLabelsRm <- function(fSOMObject,metaClustFactors,mainTitle,nbRm=0)
{
     fSOM4Plot=list(
        map=fSOMObject$map,
        prettyColnames=fSOMObject$prettyColnames,
       MST=fSOMObject$MST)
    if (nbRm>0)
    {
       indexKeep =  which(fSOMObject$MST$size > sort(fSOMObject$MST$size)[nbRm])
       indexRemove  = setdiff((1:length(fSOMObject$MST$size)),indexKeep)
       fSOM4Plot$MST$size=fSOMObject$MST$size[indexKeep]
       fSOM4Plot$map$medianValues=fSOMObject$map$medianValues[indexKeep,]
       fSOM4Plot$MST$graph=igraph::induced_subgraph(fSOMObject$MST$graph,indexKeep)
       fSOM4Plot$MST$l = fSOMObject$MST$l[indexKeep,]
       PlotLabels(fSOM4Plot,as.factor(metaClustFactors[indexKeep]), main=mainTitle)
    }
    else
        {PlotLabels(fSOM4Plot,as.factor(metaClustFactors), main=mainTitle)}
    }

## User tool: tree representaton of metacluster, given size and marker representation, removing a given number of smallest metacluster
PlotStarsMSTRm <- function(fSOMObject,metaClustFactors,mainTitle,nbRm=0)
{
   fSOM4Plot=list(
        map=fSOMObject$map,
        prettyColnames=fSOMObject$prettyColnames,
       MST=fSOMObject$MST)
    if (nbRm>0)
    {
       indexKeep =  which(fSOMObject$MST$size > sort(fSOMObject$MST$size)[nbRm])
       indexRemove  = setdiff((1:length(fSOMObject$MST$size)),indexKeep)
       fSOM4Plot$MST$size=fSOMObject$MST$size[indexKeep]
       fSOM4Plot$map$medianValues=fSOMObject$map$medianValues[indexKeep,]
       fSOM4Plot$MST$graph=igraph::induced_subgraph(fSOMObject$MST$graph,indexKeep)
       fSOM4Plot$MST$l = fSOMObject$MST$l[indexKeep,]
       PlotStarsBigLeg(fSOM4Plot,backgroundValues = as.factor(metaClustFactors[indexKeep]), main=mainTitle)
    }
    else
        {PlotStarsBigLeg(fSOM4Plot,backgroundValues = as.factor(metaClustFactors), main=mainTitle)}

}

## User tool: tree representaton of metacluster, given size and marker representation on a subset of samples given by an list of index, removing a given number of smallest metacluster
PlotStarsMSTCondRm <- function(fSOMObject,metaClustFactors,condIndex,mainTitle,nbRm=0)
{
    fSOM4Plot=list(
        map=fSOMObject$map,
        prettyColnames=fSOMObject$prettyColnames,
        MST=fSOMObject$MST)
    dataIndex = unlist(sapply(fSOMObject$metaData[condIndex],function(x){x[1]:x[2]}))
    clSizes = sapply(1:length(fSOMObject$map$medianValues[,1]),function(x){length(which(fSOMObject$map$mapping[dataIndex,1] == x))})
    if (nbRm>0)
    {
        indexKeep =  which(fSOMObject$MST$size > sort(fSOMObject$MST$size)[nbRm])
        indexRemove  = setdiff((1:length(fSOMObject$MST$size)),indexKeep)
        fSOM4Plot$MST$size = (sqrt(clSizes)/max(sqrt(clSizes))*15)[indexKeep]
        fSOM4Plot$map$medianValues = t(sapply(1:length(fSOMObject$map$medianValues[,1]),function(i){apply(fSOMObject$data[intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i)),,drop=F],2,function(x){median(x)})}))[indexKeep,]
        fSOM4Plot$MST$graph=igraph::induced_subgraph(fSOMObject$MST$graph,indexKeep)
        fSOM4Plot$MST$l = fSOMObject$MST$l[indexKeep,]
        PlotStarsBigLeg(fSOM4Plot,backgroundValues = as.factor(metaClustFactors[indexKeep]), main=mainTitle)
    }
    else {
        fSOM4Plot$MST$size = sqrt(clSizes)/max(sqrt(clSizes))*15
        fSOM4Plot$map$medianValues = t(sapply(1:length(fSOMObject$map$medianValues[,1]),function(i){apply(fSOMObject$data[intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i)),,drop=F],2,function(x){median(x)})}))
         PlotStarsBigLeg(fSOM4Plot,backgroundValues = as.factor(metaClustFactors), main=mainTitle)
        }
}

## User tool: marker level represented on metacluster tree, removing a given number of smallest metacluster
PlotMarkerMSTRm <- function(fSOMObject,markerName,mainTitle,nbRm=0,globalMinMax=c())
{
   fSOM4Plot=list(
        map=fSOMObject$map,
        prettyColnames=fSOMObject$prettyColnames,
       MST=fSOMObject$MST)
    if (nbRm>0)
    {
       indexKeep =  which(fSOMObject$MST$size > sort(fSOMObject$MST$size)[nbRm])
       indexRemove  = setdiff((1:length(fSOMObject$MST$size)),indexKeep)
       fSOM4Plot$MST$size=fSOMObject$MST$size[indexKeep]
       fSOM4Plot$map$medianValues=fSOMObject$map$medianValues[indexKeep,]
       fSOM4Plot$MST$graph=igraph::induced_subgraph(fSOMObject$MST$graph,indexKeep)
       fSOM4Plot$MST$l = fSOMObject$MST$l[indexKeep,]
    }

 if (length(globalMinMax) > 0)
     {
         maxGlobal = globalMinMax[2]
         minGlobal = globalMinMax[1]
         maxCond=max(fSOM4Plot$map$medianValues[,markerName])
         minCond=min(fSOM4Plot$map$medianValues[,markerName])
         colorPalette1000=grDevices::colorRampPalette(c("#00007F", "blue","#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(1000)
        minIndex=(minCond-minGlobal)/(maxGlobal-minGlobal)*1000
        maxIndex=(maxCond-minGlobal)/(maxGlobal-minGlobal)*1000
        colorIndexRound=round(minIndex+(0:8)*(maxIndex-minIndex)/8)
        colorPaletteCond=grDevices::colorRampPalette(colorPalette1000[colorIndexRound])

         FlowSOM::PlotMarker(fSOM4Plot,marker=markerName, view = "MST",main=mainTitle,colorPalette = colorPaletteCond)
     }
   else {FlowSOM::PlotMarker(fSOM4Plot,marker=markerName, view = "MST",main=mainTitle)
   }
}

## User tool: marker level represented on metacluster tree, on a subset of samples given by an list of index, removing a given number of smallest metacluster
PlotMarkerMSTCondRm <- function(fSOMObject,markerName,condIndex,mainTitle,nbRm=0,globalMinMax=c()){
    fSOM4Plot=list(
        map=fSOMObject$map,
        prettyColnames=fSOMObject$prettyColnames,
        MST=fSOMObject$MST)
    dataIndex = unlist(sapply(fSOMObject$metaData[condIndex],function(x){x[1]:x[2]}))
    clSizes = sapply(1:length(fSOMObject$map$medianValues[,1]),function(x){length(which(fSOMObject$map$mapping[dataIndex,1] == x))})
    if (nbRm>0) {
        indexKeep =  which(fSOMObject$MST$size > sort(fSOMObject$MST$size)[nbRm])
        indexRemove  = setdiff((1:length(fSOMObject$MST$size)),indexKeep)
        fSOM4Plot$MST$size = (sqrt(clSizes)/max(sqrt(clSizes))*15)[indexKeep]
        fSOM4Plot$map$medianValues = t(sapply(1:length(fSOMObject$map$medianValues[,1]),function(i){
            if(length(intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i))) == 0 )
            {print(paste("Cluster ",i," has size zero for given condition, use global median value"))
                apply(fSOMObject$data[which(fSOMObject$map$mapping[,1] == i),,drop=F],2,function(x){median(x)})}
            else {apply(fSOMObject$data[intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i)),,drop=F],2,function(x){median(x)})}
        }))[indexKeep,]
        ##print(fSOM4Plot$map$medianValues)
        fSOM4Plot$MST$graph=igraph::induced_subgraph(fSOMObject$MST$graph,indexKeep)
        fSOM4Plot$MST$l = fSOMObject$MST$l[indexKeep,]

    }
    else {
        fSOM4Plot$MST$size = sqrt(clSizes)/max(sqrt(clSizes))*15
        fSOM4Plot$map$medianValues = t(sapply(1:length(fSOMObject$map$medianValues[,1]),function(i){
            if(length(intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i))) == 0 )
            {print(paste("Cluster ",i," has size zero for given condition, use global median value"))
                apply(fSOMObject$data[which(fSOMObject$map$mapping[,1] == i),,drop=F],2,function(x){median(x)})}
            else{apply(fSOMObject$data[intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i)),,drop=F],2,function(x){median(x)})}
        }))
    }


     if (length(globalMinMax) > 0)
     {
         maxGlobal = globalMinMax[2]
         minGlobal = globalMinMax[1]
         maxCond=max(fSOM4Plot$map$medianValues[,markerName])
         minCond=min(fSOM4Plot$map$medianValues[,markerName])
         colorPalette1000=grDevices::colorRampPalette(c("#00007F", "blue","#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(1000)
        minIndex=(minCond-minGlobal)/(maxGlobal-minGlobal)*1000
        maxIndex=(maxCond-minGlobal)/(maxGlobal-minGlobal)*1000
        colorIndexRound=round(minIndex+(0:8)*(maxIndex-minIndex)/8)
        colorPaletteCond=grDevices::colorRampPalette(colorPalette1000[colorIndexRound])
        FlowSOM::PlotMarker(fSOM4Plot, marker=markerName, view = "MST", main=mainTitle,colorPalette = colorPaletteCond)
     }
     else {FlowSOM::PlotMarker(fSOM4Plot, marker=markerName, view = "MST", main=mainTitle)}
}


## User tool: Download data, given fcs files, FlowJo workspace should in in current environment, FCS directory is inside wd, given with no "/"
DownLoadCytoData <- function(dirFCS="",gatingName,fcsPattern = "Tube",compensate=FALSE ){
    flowJoWS=list.files(pattern=".wsp")
    if (length(flowJoWS) > 1)
    {
        stop("several FlowJoWorkSpace")
    }
    if (dirFCS == "") {absoluteDirFCS = getwd() }
    else {
        absoluteDirFCS=paste(getwd(),"/",dirFCS,"/",sep="")
        }
    files<-list.files(path = absoluteDirFCS , pattern=fcsPattern)
    if ((unlist(packageVersion("CytoML"))[1] == 1) & (unlist(packageVersion("CytoML"))[2] >= 12)) {
     ## print("Parse flowJo within version 1.12 of CytoML")
        data<-parse_flowjo_CytoML_v12(files,flowJoWS)}
    else {
            data<-parse_flowjo_CytoML(files,flowJoWS)
        }
    dataGated<-gating_subset_toolBox(data,gatingName)
    fSOM<-FlowSOM::ReadInput(dataGated$flowSet,compensate = compensate,transform = FALSE,scale = FALSE,scaled.center = TRUE,scaled.scale = TRUE,silent = FALSE)
    return(list(fSOMData=fSOM,flJoDataGated=dataGated))
}

## User tool: build FSOM tree with the metacluster, plot the tree
buildFSOMTree <- function(fSOMDloaded,prettyNames,clustDim,metaClNb,fSOMSeed)
{
  set.seed(fSOMSeed)
  if (length(which(names(fSOMDloaded) == "fSOMData")) > 0 )
  {fSOMData = fSOMDloaded$fSOMData} else {fSOMData = fSOMDloaded}
    ## ff<-fSOMDloaded$flJoDataGated$flowSet[[1]]
    fSOMNicePrettyColNames=gsub(" <.*", "", fSOMData$prettyColnames)
    colNamesIndices=unlist(lapply(prettyNames,function(name){which(fSOMNicePrettyColNames == name)}))
    print("Catched col indices:")
    print(colNamesIndices)
    ##channels_of_interest <-  fSOMData$prettyColnames[colNamesIndices]
    fSOM<-FlowSOM::BuildSOM(fSOMData,colsToUse = colNamesIndices,silent = FALSE,xdim=clustDim,ydim=clustDim,rlen=10,init=FALSE,distf=2)
    fSOM<-FlowSOM::BuildMST(fSOM,silent = FALSE,tSNE=FALSE)
    fSOM$prettyColnames =  fSOMNicePrettyColNames
    metacl<-FlowSOM::metaClustering_consensus(fSOM$map$codes,k=metaClNb,seed=fSOMSeed)
    PlotStarsBigLeg(fSOM,backgroundValues = as.factor(metacl))
    return(list(fSOMTree = fSOM,metaCl = metacl))
}

## User tool: plot figures, use the object created by buildFSOMTree
## the treatmentTable should be a dataframe with two column: "Treatment", "files"
plotTreeSet <- function(TreeMetacl,markers,Title,rmClNb=0,treatmentTable,globalScale=T){
    if (rmClNb>0) {
        indexKeep =  which(TreeMetacl$fSOMTree$MST$size > sort(TreeMetacl$fSOMTree$MST$size)[rmClNb])} else {indexKeep = 1:length(TreeMetacl$fSOMTree$MST$size)}
    pdf(file=paste(Title,"_TreatmentTree.pdf",sep=""))
    ## plot tree of pooled data
    PlotStarsMSTRm(TreeMetacl$fSOMTree,TreeMetacl$metaCl,paste(Title," MainTree",sep=""),rmClNb)

    Treatments=unique(treatmentTable$Treatment[which(sapply(tableTreatmentFCS$files,function(files){length(grep(files,names(TreeMetacl$fSOMTree$metaData),fixed=T))>0}))])
    print("Treatments:")
    print(Treatments)
    ## plot tree of subdata for each treatments
    for (treatName in Treatments) {
        fcsFiles=treatmentTable$files[which(treatmentTable$Treatment == treatName)]
        treatIndex = unlist(sapply(fcsFiles,function(file){grep(file,names(TreeMetacl$fSOMTree$metaData),fixed=T)}))
        print(paste("Plot tree for treatment",treatName))
        PlotStarsMSTCondRm(TreeMetacl$fSOMTree,TreeMetacl$metaCl,treatIndex,paste(Title," Treat: ",treatName,sep=""),rmClNb)
        }
    dev.off()
    pdf(file=paste(Title,"_MarkerTree.pdf",sep=""))
    markersInData = intersect(markers,TreeMetacl$fSOMTree$prettyColnames)
    print("Markers:")
    print(markersInData)
    ##plot tree for each marker
    for (marker in markersInData)
    {
        uglyName = names(which(TreeMetacl$fSOMTree$prettyColnames == marker))
        if (globalScale) ## construct the global scale for color scale of markers
        {
            listTreatmentIndices = lapply(unique(treatmentTable$Treatment),function(treatName){
                fcsFiles=treatmentTable$files[which(treatmentTable$Treatment == treatName)]
                CondIndex = unlist(sapply(fcsFiles,function(file){grep(file,names(TreeMetacl$fSOMTree$metaData),fixed=T)}))
                return(unlist(sapply(TreeMetacl$fSOMTree$metaData[CondIndex],function(x){x[1]:x[2]})))
            })
            maxMin = matrix(unlist(lapply((1:length(TreeMetacl$fSOMTree$map$medianValues[,1]))[indexKeep],function(cluster){
                medianList=sapply(listTreatmentIndices,function(indices){
                    median(TreeMetacl$fSOMTree$data[intersect(indices,which(TreeMetacl$fSOMTree$map$mapping[,1] == cluster)),uglyName])});
                c(max(medianList,na.rm=T),min(medianList,na.rm=T))
            })),nrow=2)
            maxGlobal = max(maxMin[1,])
            minGlobal = min(maxMin[2,])
            print(paste("Scale for marker",marker,":",minGlobal,"->",maxGlobal))
            globMinMax=c(minGlobal,maxGlobal)
        }
        else {globMinMax=c()}
        ##print(paste("Marker: ",marker,sep=""))
        print(paste("Plot tree for marker",marker," -- ",uglyName))
        PlotMarkerMSTRm(TreeMetacl$fSOMTree,uglyName,paste(Title," Marker: ",marker,sep=""),rmClNb,globMinMax)
    }
    for (marker in markersInData){ ## plot tree for each marher and each treatment
         uglyName = names(which(TreeMetacl$fSOMTree$prettyColnames == marker))
        if (globalScale) ## construct the global scale for color scale of markers
        {
            listTreatmentIndices = lapply(unique(treatmentTable$Treatment),function(treatName){
                fcsFiles=treatmentTable$files[which(treatmentTable$Treatment == treatName)]
                CondIndex = unlist(sapply(fcsFiles,function(file){grep(file,names(TreeMetacl$fSOMTree$metaData),fixed=T)}))
                return(unlist(sapply(TreeMetacl$fSOMTree$metaData[CondIndex],function(x){x[1]:x[2]})))
            })
            maxMin = matrix(unlist(lapply((1:length(TreeMetacl$fSOMTree$map$medianValues[,1]))[indexKeep],function(cluster){medianList=sapply(listTreatmentIndices,function(indices){median(TreeMetacl$fSOMTree$data[intersect(indices,which(TreeMetacl$fSOMTree$map$mapping[,1] == cluster)),uglyName])});c(max(medianList,na.rm=T),min(medianList,na.rm=T))})),nrow=2)
            maxGlobal = max(maxMin[1,])
            minGlobal = min(maxMin[2,])
            globMinMax=c(minGlobal,maxGlobal)
        }
        else {globMinMax=c()}
        ##print(paste("Marker: ",marker,sep=""))
        for (treatName in Treatments){
            ##print(paste("Treatment: ",treatName,sep=""))
            fcsFiles=treatmentTable$files[which(treatmentTable$Treatment == treatName)]
            treatIndex = unlist(sapply(fcsFiles,function(file){grep(file,names(TreeMetacl$fSOMTree$metaData),fixed=T)}))
             print(paste("Plot tree for marker",marker," -- ",uglyName,"and treatment",treatName))
            PlotMarkerMSTCondRm(TreeMetacl$fSOMTree,uglyName,treatIndex,paste(Title," Marker: ",marker," Treat: ",treatName,sep=""),rmClNb,globMinMax)
        }
    }
    dev.off()
}

## User tool: Box plot of metacluster, either percentage or normalized size is Norm = T.
## TreatmentTable should be a dataframe with two column: "Treatment", "files" (a third one with column "NormalizationFactor" if Norm=T).
## Robust specifies either Tukey/lm or Dunn (non adjusted p-values).
## ClustHeat=FALSE for no clustering on heatmap.
BoxPlotMetaClust = function(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab,Norm=FALSE,Robust = TRUE,ClustHeat=TRUE) {
   BoxPlotMetaClustFull(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab,Norm,Marker="",Robust,ClustHeat)
}

## User tool: Box plot of metacluster, For a given marker.
## treatmentTable should be a dataframe with two column: "Treatment", "files". Robust specifies either Tukey/lm or Dunn (non adjusted p-values).
## ClustHeat=FALSE for no clustering on heatmap.
BoxPlotMarkerMetaClust = function(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,Marker,Robust=TRUE,ClustHeat=TRUE) {
    MarkerIndex = which(TreeMetaCl$fSOMTree$prettyColnames == Marker)
    if (length(MarkerIndex) < 1) {stop("No marker ",Marker," in data")} else {
    BoxPlotMetaClustFull(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab="",Norm=FALSE,Marker,Robust,ClustHeat) }
}

MetaClusterNaming <- function(TreeMetaCl,Markers)
{
    MarkerIn = sapply(Markers,function(Marker){length(which(TreeMetaCl$fSOMTree$prettyColnames == Marker))})
    Markers = Markers[which(MarkerIn > 0)]
    print(paste("Use Marker",Markers))


    metaClustListName=lapply(unique(TreeMetaCl$metaCl),function(metaClust){
            clusterList=which(TreeMetaCl$metaCl == metaClust)
            metaClustIndices=unlist(sapply(clusterList,function(cluster){which(TreeMetaCl$fSOMTree$map$mapping[,1] == cluster)}))
            nameList = sapply(Markers,function(Marker){
                MarkerIndex=which(TreeMetaCl$fSOMTree$prettyColnames == Marker)
                metaClustMedian=median(TreeMetaCl$fSOMTree$data[metaClustIndices,MarkerIndex],na.rm=T)
                vectQuantiles=quantile(TreeMetaCl$fSOMTree$data[,MarkerIndex],na.rm=T)
                nonRobustName = c(paste(Marker,"-",sep=""),paste(Marker,"+",sep=""))[as.numeric(metaClustMedian > vectQuantiles[3])+1]
                robustName = c(paste(Marker,"-",sep=""),paste(Marker,"med",sep=""),paste(Marker,"+",sep=""))[as.numeric(metaClustMedian > vectQuantiles[2])+as.numeric(metaClustMedian > vectQuantiles[4])+1]
                robustNameNoMed= c(paste(Marker,"-",sep=""),"",paste(Marker,"+",sep=""))[as.numeric(metaClustMedian > vectQuantiles[2])+as.numeric(metaClustMedian > vectQuantiles[4])+1]
                return(c(nonRobustName,robustName,robustNameNoMed))
            })
            c(metaClust,apply(nameList,1,function(l){paste(l,collapse="")}))
    })
    metaClustDF=as.data.frame(do.call(rbind,metaClustListName),stringsAsFactors=F)
    names(metaClustDF)=c("number","nonRobustName","robustName","shortRobustName")
    return(metaClustDF)
}

## To do:  put plotlabel in plotTreeSet; automatic naming of metaclusters (median or quartile); extract sub data from a single metacluster
DataFromMetaClust <- function(FSOMData,TreeMetaCl,MetaClusters)
{
  newFSOMData=FSOMData[-2]
  Clusters = unlist(lapply(MetaClusters,function(x){which(TreeMetaCl$metaCl == x)}))
  ClusterIndices = unlist(lapply(Clusters,function(x){which(TreeMetaCl$fSOMTree$map$mapping[,1] == x)}))
  newFSOMData$fSOMData$data=FSOMData$fSOMData$data[ClusterIndices,]
  metaDataLengthKept=lapply(FSOMData$fSOMData$metaData,function(x){
    length(intersect((x[1]:x[2]),ClusterIndices))
  })
  if (length(which(metaDataLengthKept < 1))) {print(paste("Remove files:",names(which(metaDataLengthKept < 1))))}
  keepFilesIndices = which(metaDataLengthKept > 0)
  LastFilesIndex=cumsum(metaDataLengthKept[keepFilesIndices])
  FirstFilesIndex=c(1,LastFilesIndex[-length(LastFilesIndex)]+1)
  newMetaData=lapply(1:length(LastFilesIndex),function(x){unname(c(FirstFilesIndex[x],LastFilesIndex[x]))})
  names(newMetaData)=names(LastFilesIndex)
  newFSOMData$fSOMData$metaData=newMetaData
  return(newFSOMData)
}

# plot2DDensityMetaCl = function(TreeMetaCl,marker1,marker2,metaClusterNames,file)
#   {
#   pdf(file)
#   #ncolGraph = as.integer(sqrt(length(TreeMetaCl$fSOMTree$metaData)))
#   #nrowGraph = as.integer(length(TreeMetaCl$fSOMTree$metaData)/ncolGraph) + 1
#   #par(mfrow=c(nrowGraph,ncolGraph))
# for (index in 1:length(TreeMetaCl$fSOMTree$metaData))
#   {
#   
#   mData = as.data.frame(
#     TreeMetaCl$fSOMTree$data[
#       ((TreeMetaCl$fSOMTree$metaData[index][[1]][1]):(TreeMetaCl$fSOMTree$metaData[index][[1]][2])),
#       c(which(TreeMetaCl$fSOMTree$prettyColnames == marker1),which(TreeMetaCl$fSOMTree$prettyColnames == marker2))])
#   names(mData) = c("marker1","marker2")
#   print(str(mData))
#   mData$Cl = TreeMetaCl$fSOMTree$map$mapping[(TreeMetaCl$fSOMTree$metaData[index][[1]][1]):(TreeMetaCl$fSOMTree$metaData[index][[1]][2]),1]
#   mData$metaCl = TreeMetaCl$metaCl[mData$Cl]
#   print(str(mData$metaCl))
#   tmpIndex = is.element(mData$metaCl,metaClusterNames)
#   mData = mData[which(tmpIndex),]
#   mData$metaCl 
#   print(str(mData$metaCl))
#   print(ggplot2::ggplot(mData,ggplot2::aes(x=marker1,y=marker2)) + geom_density_2d(ggplot2::aes(color = metaCl)) + 
#           xlab(marker1) + ylab(marker2) +ggplot2::ggtitle(gsub(".*/","",names(TreeMetaCl$fSOMTree$metaData)[index])))
#   }
#   dev.off()
# }
