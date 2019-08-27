## Authors: Gautier Stoll, Hélène Fohrer-Ting, Estelle Devêvre, Sarah LEVESQUE, Julie LE NAOUR, Juliette PAILLET, Jonathan POL
## 2019, INSERM U1138
## Version 0.7

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
if (tmpIsV3p6) {
    library(CytoML)
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
      as.matrix() %>% apply(1,as.numeric) %>% t()
  if(!is.null(meta_names)) colnames(pctgs_meta) <- meta_names
  return(list("abstgs" = pctgs,
              "abstgs_meta" = pctgs_meta))
}

##Internal tool: return p-value of Tukey test, given metacluster names
TukeyTestSarah = function(fSOMTable, metaClust){TukeyHSD(aov(as.formula(paste(metaClust,"~ Treatment")),data=fSOMTable))$Treatment[,4]}

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
       PlotStars(fSOM4Plot,backgroundValues = as.factor(metaClustFactors[indexKeep]), main=mainTitle)
    }
    else
        {PlotStars(fSOM4Plot,backgroundValues = as.factor(metaClustFactors), main=mainTitle)}

}

## User tool: tree representaton of metacluster, given size and marker representation on a subset of samples given by an list of index, removing a given number of smallest metacluster
PlotStarsMSTCondRm=function(fSOMObject,metaClustFactors,condIndex,mainTitle,nbRm=0)
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
        PlotStars(fSOM4Plot,backgroundValues = as.factor(metaClustFactors[indexKeep]), main=mainTitle)
    }
    else {
        fSOM4Plot$MST$size = sqrt(clSizes)/max(sqrt(clSizes))*15
        fSOM4Plot$map$medianValues = t(sapply(1:length(fSOMObject$map$medianValues[,1]),function(i){apply(fSOMObject$data[intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i)),,drop=F],2,function(x){median(x)})}))
         PlotStars(fSOM4Plot,backgroundValues = as.factor(metaClustFactors), main=mainTitle)
        }
}

## User tool: marker level represented on metacluster tree, removing a given number of smallest metacluster
PlotMarkerMSTRm <- function(fSOMObject,markerName,mainTitle,nbRm=0)
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
       PlotMarker(fSOM4Plot,marker=markerName, view = "MST",main=mainTitle)
    }
    else
    {PlotMarker(fSOM4Plot,marker=markerName, view = "MST",main=mainTitle)}
}

## User tool: marker level represented on metacluster tree, on a subset of samples given by an list of index, removing a given number of smallest metacluster
PlotMarkerMSTCondRm <- function(fSOMObject,markerName,condIndex,mainTitle,nbRm=0,globeScale=F){
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
    PlotMarker(fSOM4Plot, marker=markerName, view = "MST", main=mainTitle)
}

## User tool: Download data, given fcs files, FlowJo workspace should in in current environment, FCS directory is inside wd, given with no "/"
DownLoadCytoData <- function(dirFCS="",gatingName,fcsPattern = "Tube"){
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
    fSOM<-ReadInput(dataGated$flowSet,compensate = FALSE,transform = FALSE,scale = FALSE,scaled.center = TRUE,scaled.scale = TRUE,silent = FALSE)
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
    PlotStars(fSOM,backgroundValues = as.factor(metacl))
    return(list(fSOMTree = fSOM,metaCl = metacl))
}

## User tool: plot figures, use the object created by buildFSOMTree
## the treatmentTable should be a dataframe with two column: "Treatment", "files"
plotTreeSet <- function(TreeMetacl,markers,Title,rmClNb=0,treatmentTable){
    pdf(file=paste(Title,"_TreatmentTree.pdf",sep=""))
    PlotStarsMSTRm(TreeMetacl$fSOMTree,TreeMetacl$metaCl,paste(Title," MainTree",sep=""),rmClNb)
    Treatments=unique(treatmentTable$Treatment)
    for (treatName in Treatments) {
        print(paste("Treatment: ",treatName,sep=""))
        fcsFiles=treatmentTable$files[which(treatmentTable$Treatment == treatName)]
        treatIndex = unlist(sapply(fcsFiles,function(file){grep(file,names(TreeMetacl$fSOMTree$metaData),fixed=T)}))
        PlotStarsMSTCondRm(TreeMetacl$fSOMTree,TreeMetacl$metaCl,treatIndex,paste(Title," Treat: ",treatName,sep=""),rmClNb)
        }
    dev.off()
    pdf(file=paste(Title,"_MarkerTree.pdf",sep=""))
    for (marker in markers)
    {
        print(paste("Marker: ",marker,sep=""))
        uglyName = names(which(TreeMetacl$fSOMTree$prettyColnames == marker))
        PlotMarkerMSTRm(TreeMetacl$fSOMTree,uglyName,paste(Title," Marker: ",marker,sep=""),rmClNb)
    }
    for (marker in markers){
        print(paste("Marker: ",marker,sep=""))
        uglyName = names(which(TreeMetacl$fSOMTree$prettyColnames == marker))
        for (treatName in Treatments){
            print(paste("Treatment: ",treatName,sep=""))
            fcsFiles=treatmentTable$files[which(treatmentTable$Treatment == treatName)]
            treatIndex = unlist(sapply(fcsFiles,function(file){grep(file,names(TreeMetacl$fSOMTree$metaData),fixed=T)}))
            PlotMarkerMSTCondRm(TreeMetacl$fSOMTree,uglyName,treatIndex,paste(Title," Marker: ",marker," Treat: ",treatName,sep=""),rmClNb)
        }
    }
    dev.off()
}

## User tool: Box plot of metacluster, either percentage or normlized size is Norm = T
## treatmentTable should be a dataframe with two column: "Treatment", "files" (a third one with column "NormalizationFactor" in Norm=T
BoxPlotMetaClust <- function(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab,Norm=FALSE,Plot=T)
{
    if (Norm) {
        abstgs=get_abstgsMT(TreeMetaCl$fSOMTree,TreeMetaCl$metaCl)
        fSOMnbrs<-abstgs$abstgs_meta
        NormFactors = sapply(row.names(fSOMnbrs),function(fileFCS){treatmentTable$NormalizationFactor[which(treatmentTable$files == fileFCS)]})
        for (index in 1:length(NormFactors)) {fSOMnbrs[index,] = fSOMnbrs[index,]/NormFactors[index] }
        PlotLab=yLab
    }
    else {
        pctgs<-get_pctgsMT(TreeMetaCl$fSOMTree,TreeMetaCl$metaCl)
        fSOMnbrs<-pctgs$pctgs_meta
        fSOMnbrs<-fSOMnbrs*100
        PlotLab=paste("% of ",yLab,sep="")
    }
    treatmentsFSOM=sapply(row.names(fSOMnbrs),function(fileFCS){treatmentTable$Treatment[which(treatmentTable$files == fileFCS)]})
    Treatments=unique(treatmentTable$Treatment)
    treatmentsFSOM=factor(treatmentsFSOM,levels=c(ControlTreatment,setdiff(Treatments,ControlTreatment))) # set control treatment in first
    if (Norm) {pdf(file=paste(Title,"_BoxPlotNormMetacl.pdf",sep=""))}
    else {pdf(file=paste(Title,"_BoxPlotPercentMetacl.pdf",sep=""))}
    metaclNumber=length(fSOMnbrs[1,])
    par(mfrow=c(6,6),las=2,mar=c(BottomMargin,3,1,.5),mgp=c(1.8,.8,0))
    for (metaCl in (1:metaclNumber)){
        plotDf=data.frame(PP=fSOMnbrs[,metaCl],TreatmentFSOM=treatmentsFSOM) ## dataframe for box plot
        boxplot(PP ~ TreatmentFSOM,data=plotDf,main=paste("mtcl",metaCl,sep=""),xlab="",ylab=PlotLab,cex.axis=.5,cex.main=.8,cex.lab=.5)
        beeswarm(PP ~ TreatmentFSOM,data=plotDf,main=paste("mtcl",metaCl,sep=""),add=T,cex=.5,col="red")
    }
    par(mfrow=c(6,6),las=2,mar=c(BottomMargin,3,1,.5),mgp=c(1.8,.8,0))
    PvalTable = sapply((1:metaclNumber),function(metaCl)
    {
        plotDf=data.frame(PP=fSOMnbrs[,metaCl],TreatmentFSOM=treatmentsFSOM)
        tukeyPval=TukeyHSD(aov(PP ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM[,4]
        names(tukeyPval)=row.names(TukeyHSD(aov(PP ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM)
        ListSignif=(sapply(1:length(tukeyPval),function(index){
            if(tukeyPval[index] < 0.0001){return(c("****",strsplit(names(tukeyPval)[index],split = "-")[[1]]))}
            else if(tukeyPval[index] < 0.001){return(c("***",strsplit(names(tukeyPval)[index],split = "-")[[1]]))}
            else if(tukeyPval[index] < 0.01){return(c("**",strsplit(names(tukeyPval)[index],split = "-")[[1]]))}
            else if(tukeyPval[index] < 0.05){return(c("*",strsplit(names(tukeyPval)[index],split = "-")[[1]]))}
        }))
        ListSignif = ListSignif[which(sapply(ListSignif,length) > 0)]
        ListSignifPosIndex = lapply(ListSignif,function(hit){
            return(c(which(levels(plotDf$TreatmentFSOM) == hit[2]),which(levels(plotDf$TreatmentFSOM) == hit[3])))})
        minTr=min(plotDf$PP)
        maxTr=max(plotDf$PP)
        boxplot(PP ~ TreatmentFSOM,
                data=plotDf,main=paste("mtcl",metaCl,sep=""),
                xlab="",
                ylab=paste("% of ",yLab,sep=""),
                cex.axis=.5,
                cex.main=.8,
                cex.lab=.5,
                ylim=c(minTr,length(ListSignif)*abs(maxTr-minTr)*.2+maxTr)
                )
        if (length(ListSignif) > 0) {
            if (length(tukeyPval) > 1) ## more than one pair of comparison
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

       return(tukeyPval)
    })
    ##print(PvalTable)
    ##tmpTable=PvalTable[,1]
    if(is.matrix(PvalTable)) {
        PvalTable=as.data.frame(PvalTable)
    }
    else {PvalTable=as.data.frame(t(PvalTable))}
    names(PvalTable)=paste("mtcl",(1:metaclNumber),sep="")
    par(mfrow=c(1,1),mar=c(3,2,3,1),cex=.5)
    plot.new()
    tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),base_size = 3*(35/metaclNumber))
    grid.table(t(PvalTable),theme = tt)
    dev.off()
    if (Norm) {write.table(PvalTable,paste(Title,"_TukeyPvalNormMetacl.csv",sep=""),sep=";",col.names = NA)}
    else {write.table(PvalTable,paste(Title,"_TukeyPvalPercentMetacl.csv",sep=""),sep=";",col.names = NA)}
    retData=list(fSOMnbrs,PvalTable)
    names(retData)=c("Sizes","Pval")
    return(retData)
}
