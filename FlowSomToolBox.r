## Authors: Gautier Stoll, Hélène Fohrer-Ting, Estelle Devêvre, Sarah LEVESQUE, Julie LE NAOUR, Juliette PAILLET, Jonathan POL
## 2019, INSERM U1138
## Version 0.8.1

tmpIsV3p6 = (as.integer(strsplit(strsplit(version$version.string,split=" ")[[1]][3],split=".",fixed=TRUE)[[1]][1]) >= 3) & (as.integer(strsplit(strsplit(version$version.string,split=" ")[[1]][3],split=".",fixed=TRUE)[[1]][2]) >= 6) ## for testing R version


library(flowCore)
library(flowDensity)
library(flowWorkspace)
library(flowAI)
library(FlowSOM)
library(FlowSOMworkshop)
library(Rtsne)
library("pheatmap")

library("data.table")

library(reshape)
library(reshape2)

library(ggplot2)

library(ggpubr)

library(plyr)
library(readr)
library(dplyr)
library(beeswarm)

library(gridExtra)
library(gplots)
library(dunn.test)
if (tmpIsV3p6) {
    library(CytoML)
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
    add.vertex.shape("star", clip=igraph.shape.noclip, plot=mystarBL, 
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


## Internal tool: corrected parse_flowjo

if (tmpIsV3p6) {
    parse_flowjo_CytoML <- function (files, wsp_file, group = "All Samples", plot = FALSE)
    {
        wsp <- flowWorkspace::openWorkspace(wsp_file)
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
}
## Internal tool: GetClusters ?
if(!exists("GetClusters",mode="function")) {
    GetClusters <- function(fsom) {
      if (class(fsom) == "list" & !is.null(fsom$FlowSOM)) {
        fsom <- fsom$FlowSOM
      }
      if (class(fsom) != "FlowSOM") {
        stop("fsom should be a FlowSOM object.")
      }
      return(fsom$map$mapping[,1])
    }
}
    
## Internal tool: GetMetClusters ?
if(!exists("GetMetaclusters",mode="function")) {
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
}
    
    ## Internal tool: Seems that PlotLabels diseapear...
if(!exists("PlotLabels",mode="function")) {
    
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
}
    

## Internal tool: extract meta-clusters count ratio in percent
get_pctgsMT <- function(fSOM,metacl, meta_names = NULL){
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
    treatmentsFSOM=factor(treatmentsFSOM,levels=c(ControlTreatment,setdiff(Treatments,ControlTreatment))) # set control treatment at first
    if(length(MarkerIndex) == 1) {pdf(file=paste(Title,"_BoxPlot",Marker,"Metacl.pdf",sep=""))} else {
    if (Norm) {pdf(file=paste(Title,"_BoxPlotNormMetacl.pdf",sep=""))}
    else {pdf(file=paste(Title,"_BoxPlotPercentMetacl.pdf",sep=""))}}
    metaclNumber=length(fSOMnbrs[1,])
    par(mfrow=c(6,6),las=2,mar=c(BottomMargin,3,1,.5),mgp=c(1.8,.8,0)) ## page have 6x6 boxplots
    fSOMnbrs=fSOMnbrs[,unique(unique(TreeMetaCl$metaCl))]
    cex4Title=1/(1+log(max(sapply(colnames(fSOMnbrs),nchar)))/2)
    for (metaCl in (1:metaclNumber)){ ## boxplots with no annotations
        plotDf=data.frame(PP=fSOMnbrs[,metaCl],TreatmentFSOM=treatmentsFSOM) ## dataframe for box plot
        boxplot(PP ~ TreatmentFSOM,data=plotDf,main=paste("mtcl",colnames(fSOMnbrs)[metaCl],sep="_"),xlab="",ylab=PlotLab,cex.axis=.5,cex.main=cex4Title,cex.lab=.5)
        beeswarm(PP ~ TreatmentFSOM,data=plotDf,main=paste("mtcl",colnames(fSOMnbrs)[metaCl],sep="_"),add=T,cex=.5,col="red")
    }
    par(mfrow=c(6,6),las=2,mar=c(BottomMargin,3,1,.5),mgp=c(1.8,.8,0))
    
    PvalPairwiseTable = sapply((1:metaclNumber),function(metaCl) ## construct pval table of tukey pairwise comparison test, boxplots with p-values annotation
    {
        plotDf=data.frame(PP=fSOMnbrs[,metaCl],TreatmentFSOM=treatmentsFSOM)
        if (Robust) {
            invisible(capture.output(tmp <- dunn.test(plotDf$PP,plotDf$TreatmentFSOM)))
            pairwisePval=tmp$P
            names(pairwisePval) = gsub(" ","",tmp$comparisons,fixed=T)
        } else {
            pairwisePval=TukeyHSD(aov(PP ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM[,4]
            names(pairwisePval)=row.names(TukeyHSD(aov(PP ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM)
        }
        ListSignif=(sapply(1:length(pairwisePval),function(index){
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
       beeswarm(PP ~ TreatmentFSOM,data=plotDf,add=T,cex=.5,col="red")
       return(pairwisePval)
    })
    ## finish the construction of PvalTable, write csv files 
    if(is.matrix(PvalPairwiseTable)) {
        PvalPairwiseTable=as.data.frame(PvalPairwiseTable)
    }
    else {PvalPairwiseTable=as.data.frame(t(PvalPairwiseTable))}
    names(PvalPairwiseTable)=paste("mtcl",colnames(fSOMnbrs)[1:metaclNumber],sep="_")
    par(mfrow=c(1,1),mar=c(3,2,3,1),cex=.5)
  
    if(length(MarkerIndex) == 1) {write.table(PvalPairwiseTable,paste(Title,"_PairwisePval",Marker,"Metacl.csv",sep=""),sep=";",col.names = NA)} else {
    if (Norm) {write.table(PvalPairwiseTable,paste(Title,"_PairwisePvalNormMetacl.csv",sep=""),sep=";",col.names = NA)}
    else {write.table(PvalPairwiseTable,paste(Title,"_PairwisePvalPercentMetacl.csv",sep=""),sep=";",col.names = NA)}}

    DF4lm = data.frame(y=c(fSOMnbrs),metaCl = c(sapply(colnames(fSOMnbrs),function(name){rep(name,length(fSOMnbrs[,1]))})),treat = rep(c(sapply(row.names(fSOMnbrs),function(name){as.character(treatmentTable$Treatment[which(treatmentTable$files == name)])})),length(fSOMnbrs[1,])))
    DF4lm$treat = factor(DF4lm$treat,levels=c(ControlTreatment,setdiff(unique(DF4lm$treat),ControlTreatment)))
    if (Robust) {
        pvalLmMatrix=as.matrix(PvalPairwiseTable)[which(sapply(row.names(PvalPairwiseTable),function(name){spName = strsplit(name,split = "-")[[1]];((spName[1] == ControlTreatment) |(spName[2] == ControlTreatment) )})),]
        row.names(pvalLmMatrix)=gsub(ControlTreatment,"",row.names(pvalLmMatrix),fixed=T)
        row.names(pvalLmMatrix)=gsub("-","",row.names(pvalLmMatrix),fixed=T)
    } else {
        pvalLmMatrix = t(do.call(rbind,lapply(colnames(fSOMnbrs),function(metaCl){summary(lm(y ~ treat,data = DF4lm[which(DF4lm$metaCl == metaCl),]))$coefficient[-1,4]})))
        row.names(pvalLmMatrix) = setdiff(unique(DF4lm$treat),ControlTreatment)
    }
    colnames(pvalLmMatrix) = paste("mtcl",colnames(fSOMnbrs),sep= "_")
    pvalLmMatrix = rbind(rep(1,length(fSOMnbrs[1,])),pvalLmMatrix)
    row.names(pvalLmMatrix)[1] = ControlTreatment
    DF4lm$metaCl = factor(DF4lm$metaCl,level=unique(DF4lm$metaCl)) ## to get the right ordering after "by" function
    if (Robust) {
        meanMatrix  = t(by(DF4lm$y,list(DF4lm$metaCl,DF4lm$treat),function(x){median(x,na.rm=T)}))
        pvalLmMatrix=pvalLmMatrix[row.names(meanMatrix),]
    } else {
        meanMatrix  = t(by(DF4lm$y,list(DF4lm$metaCl,DF4lm$treat),function(x){mean(x,na.rm=T)})) }
    attr(meanMatrix,"class") = NULL
    attr(meanMatrix,"call") = NULL
    colnames(meanMatrix) = paste("mtcl",colnames(meanMatrix),sep= "_")
    pvalAnnotationMatrix = apply(pvalLmMatrix,c(1,2),function(x){
    if (x < 0.0001){return("****")}
            else if(x < 0.001){return("***")}
            else if(x < 0.01){return("**")}
            else if(x < 0.05){return("*")} else {return("")}})
    if (length(MarkerIndex) == 1) {
        colMarginSize=max(c(5,.5*max(sapply(colnames(meanMatrix),nchar))))
        rowMarginSize=max(c(8,.5*max(sapply(row.names(meanMatrix),nchar))))
        if (Robust) {heatTitle = paste("Median MFI of ",PlotLab,sep="")} else {heatTitle = paste("Mean MFI of ",PlotLab,sep="")}
        if (ClustHeat) {
            heatmap.2(meanMatrix,Rowv=F,Colv=T,dendrogram = "column",scale="none",col = heat.colors(100),cellnote = pvalAnnotationMatrix,notecol = "black",trace = "none",cexRow = .8,cexCol=.8,density.info="none",main=heatTitle,notecex=.5,margins=c(colMarginSize,rowMarginSize))
        }
        heatmap.2(meanMatrix,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = heat.colors(100),cellnote = pvalAnnotationMatrix,notecol = "black",trace = "none",cexRow = .8,cexCol=.8,density.info="none",main=heatTitle,notecex=.5,margins=c(colMarginSize,rowMarginSize))
    } else {
        colMarginSize=max(c(5,.5*max(sapply(colnames(meanMatrix[-1,]),nchar))))
        rowMarginSize=max(c(8,.5*max(sapply(row.names(meanMatrix[-1,]),nchar))))
        if (Robust) { heatTitle = paste("Median ",PlotLab,sep="")}
        else {heatTitle = paste("Mean ",PlotLab,sep="")}
        heatTitle=paste(heatTitle,"\n(rel. to ",ControlTreatment,", scaled)",sep="")
        meanMatrix=apply(meanMatrix,2,function(x){(x-x[1])/sd(x,na.rm=T)})
        meanMatrix=meanMatrix[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")] ## get the correct ordering
        if (ClustHeat) {
            heatmap.2(meanMatrix[-1,],Rowv=F,Colv=T,dendrogram = "column",scale="none",col = bluered(100),cellnote = pvalAnnotationMatrix[-1,],notecol = "black",trace = "none",cexRow = .8,cexCol=.8,density.info="none",main=heatTitle,distfun=function(x){dist(t(apply(meanMatrix,2,function(y){scale(y)})))},notecex=.5,margins=c(colMarginSize,rowMarginSize))
        } 
            heatmap.2(meanMatrix[-1,],Rowv=F,Colv=F,dendrogram = "none",scale="none",col = bluered(100),cellnote = pvalAnnotationMatrix[-1,],notecol = "black",trace = "none",cexRow = .8,cexCol=.8,density.info="none",main=heatTitle,distfun=function(x){dist(t(apply(meanMatrix,2,function(y){scale(y)})))},notecex=.5,margins=c(colMarginSize,rowMarginSize))   
    }
    matrixPval4Heat=as.matrix(PvalPairwiseTable)[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")]
    colMarginSize=max(c(5,.5*max(sapply(colnames(matrixPval4Heat),nchar))))
    rowMarginSize=max(c(8,.5*max(sapply(row.names(matrixPval4Heat),nchar))))
    if (Robust) {
        if (ClustHeat) {
            heatmap.2(matrixPval4Heat,Rowv=T,Colv=T,dendrogram = "both",scale="none",col = gray((0:100)/100),trace="none",cexRow=.8,main="Dunn p-values",cexCol=.8,margins=c(colMarginSize,rowMarginSize),density.info="none")
}       
            heatmap.2(matrixPval4Heat,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = gray((0:100)/100),trace="none",cexRow=.8,main="Dunn p-values",cexCol=.8,margins=c(colMarginSize,rowMarginSize),density.info="none")
    } else {
             if (ClustHeat) {
                 heatmap.2(matrixPval4Heat,Rowv=T,Colv=T,dendrogram = "both",scale="none",col = gray((0:100)/100),trace="none",cexRow=.8,main="Tukey p-values",cexCol=.8,margins=c(colMarginSize,rowMarginSize),density.info="none")
             }
             heatmap.2(matrixPval4Heat,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = gray((0:100)/100),trace="none",cexRow=.8,main="Tukey p-values",cexCol=.8,margins=c(colMarginSize,rowMarginSize),density.info="none") }      
    dev.off()
    retData=list(fSOMnbrs,PvalPairwiseTable,pvalLmMatrix)

    if(length(MarkerIndex) == 1) {write.table(pvalLmMatrix,paste(Title,"_LmPval",Marker,"Metacl.csv",sep=""),sep=";",col.names = NA)} else {
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
       fSOM4Plot$MST$graph=induced_subgraph(fSOMObject$MST$graph,indexKeep)
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
       fSOM4Plot$MST$graph=induced_subgraph(fSOMObject$MST$graph,indexKeep)
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
        fSOM4Plot$MST$graph=induced_subgraph(fSOMObject$MST$graph,indexKeep)
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
       fSOM4Plot$MST$graph=induced_subgraph(fSOMObject$MST$graph,indexKeep)
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
   
         PlotMarker(fSOM4Plot,marker=markerName, view = "MST",main=mainTitle,colorPalette = colorPaletteCond)
     }
   else {PlotMarker(fSOM4Plot,marker=markerName, view = "MST",main=mainTitle)
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
        fSOM4Plot$MST$graph=induced_subgraph(fSOMObject$MST$graph,indexKeep)
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
        PlotMarker(fSOM4Plot, marker=markerName, view = "MST", main=mainTitle,colorPalette = colorPaletteCond)
     }
     else {PlotMarker(fSOM4Plot, marker=markerName, view = "MST", main=mainTitle)}
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
    if (tmpIsV3p6) {
        data<-parse_flowjo_CytoML(files,flowJoWS)}
    else {
            data<-parse_flowjo(files,flowJoWS)
        }
    dataGated<-gating_subset(data,gatingName)
    fSOM<-ReadInput(dataGated$flowSet,compensate = compensate,transform = FALSE,scale = FALSE,scaled.center = TRUE,scaled.scale = TRUE,silent = FALSE)
    return(list(fSOMData=fSOM,flJoDataGated=dataGated))
}
## User tool: build FSOM tree with the metacluster, plot the tree
buildFSOMTree <- function(fSOMDloaded,prettyNames,clustDim,metaClNb,fSOMSeed)
{
    ff<-fSOMDloaded$flJoDataGated$flowSet[[1]]
    fSOMNicePrettyColNames=gsub(" <.*", "", fSOMDloaded$fSOMData$prettyColnames)
    colNamesIndices=unlist(lapply(prettyNames,function(name){which(fSOMNicePrettyColNames == name)}))
    print("Catched col indices:")
    print(colNamesIndices)
    channels_of_interest <-  colnames(ff)[colNamesIndices]
    fSOM<-BuildSOM(fSOMDloaded$fSOMData,colsToUse = channels_of_interest,silent = FALSE,xdim=clustDim,ydim=clustDim,rlen=10,init=FALSE,distf=2)
    fSOM<-BuildMST(fSOM,silent = FALSE,tSNE=FALSE)
    fSOM$prettyColnames =  fSOMNicePrettyColNames
    metacl<-metaClustering_consensus(fSOM$map$codes,k=metaClNb,seed=fSOMSeed)
    PlotStarsBigLeg(fSOM,backgroundValues = as.factor(metacl))
    return(list(fSOMTree = fSOM,metaCl = metacl))
}

## User tool: plot figures, use the object created by buildFSOMTree
## the treatmentTable should be a dataframe with two column: "Treatment", "files"
plotTreeSet <- function(TreeMetacl,markers,Title,rmClNb=0,treatmentTable,globalScale=F){
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
