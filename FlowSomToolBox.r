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
library(gplots)
library(dunn.test)
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
        PlotLab=yLab
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
    for (metaCl in (1:metaclNumber)){ ## boxplots with no annotations
        plotDf=data.frame(PP=fSOMnbrs[,metaCl],TreatmentFSOM=treatmentsFSOM) ## dataframe for box plot
        boxplot(PP ~ TreatmentFSOM,data=plotDf,main=paste("mtcl",colnames(fSOMnbrs)[metaCl],sep="_"),xlab="",ylab=PlotLab,cex.axis=.5,cex.main=.8,cex.lab=.5)
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
                cex.main=.8,
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
        if (Robust) {heatTitle = paste("Median MFI of ",PlotLab,sep="")} else {heatTitle = paste("Mean MFI of ",PlotLab,sep="")}
        if (ClustHeat) {
            heatmap.2(meanMatrix,Rowv=F,Colv=T,dendrogram = "column",scale="none",col = heat.colors(100),cellnote = pvalAnnotationMatrix,notecol = "black",trace = "none",cexRow = .8,cexCol=.8,density.info="none",main=heatTitle,notecex=.5)
        }
heatmap.2(meanMatrix,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = heat.colors(100),cellnote = pvalAnnotationMatrix,notecol = "black",trace = "none",cexRow = .8,cexCol=.8,density.info="none",main=heatTitle,notecex=.5)
    } else {
        if (Robust) {if (Norm) {heatTitle = paste("Median size of ",PlotLab,sep="")} else {heatTitle = paste("Median % of ",PlotLab,sep="")}}
        else {if (Norm) {heatTitle = paste("Mean size of ",PlotLab,sep="")} else {heatTitle = paste("Mean % of ",PlotLab,sep="")}}
        heatTitle=paste(heatTitle,"(rel. to control, scaled)",sep="\n")
        meanMatrix=apply(meanMatrix,2,function(x){(x-x[1])/sd(x,na.rm=T)})
        meanMatrix=meanMatrix[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")] ## get the correct ordering
        print(meanMatrix)
        if (ClustHeat) {
            heatmap.2(meanMatrix[-1,],Rowv=F,Colv=T,dendrogram = "column",scale="none",col = bluered(100),cellnote = pvalAnnotationMatrix[-1,],notecol = "black",trace = "none",cexRow = .8,cexCol=.8,density.info="none",main=heatTitle,distfun=function(x){dist(t(apply(meanMatrix,2,function(y){scale(y)})))},notecex=.5)
        } 
            heatmap.2(meanMatrix[-1,],Rowv=F,Colv=F,dendrogram = "none",scale="none",col = bluered(100),cellnote = pvalAnnotationMatrix[-1,],notecol = "black",trace = "none",cexRow = .8,cexCol=.8,density.info="none",main=heatTitle,distfun=function(x){dist(t(apply(meanMatrix,2,function(y){scale(y)})))},notecex=.5)   
    }
    matrixPval4Heat=as.matrix(PvalPairwiseTable)[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")]
    if (Robust) {
        if (ClustHeat) {
            heatmap.2(matrixPval4Heat,Rowv=T,Colv=T,dendrogram = "both",scale="none",col = gray((0:100)/100),trace="none",cexRow=.8,main="Dunn p-values",cexCol=.8,margins=c(5,8),density.info="none")
}       
            heatmap.2(matrixPval4Heat,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = gray((0:100)/100),trace="none",cexRow=.8,main="Dunn p-values",cexCol=.8,margins=c(5,8),density.info="none")
    } else {
             if (ClustHeat) {
                 heatmap.2(matrixPval4Heat,Rowv=T,Colv=T,dendrogram = "both",scale="none",col = gray((0:100)/100),trace="none",cexRow=.8,main="Tukey p-values",cexCol=.8,margins=c(5,8),density.info="none")
             }
                 heatmap.2(matrixPval4Heat,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = gray((0:100)/100),trace="none",cexRow=.8,main="Tukey p-values",cexCol=.8,margins=c(5,8),density.info="none") }      
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
       PlotStars(fSOM4Plot,backgroundValues = as.factor(metaClustFactors[indexKeep]), main=mainTitle)
    }
    else
        {PlotStars(fSOM4Plot,backgroundValues = as.factor(metaClustFactors), main=mainTitle)}

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
        PlotStars(fSOM4Plot,backgroundValues = as.factor(metaClustFactors[indexKeep]), main=mainTitle)
    }
    else {
        fSOM4Plot$MST$size = sqrt(clSizes)/max(sqrt(clSizes))*15
        fSOM4Plot$map$medianValues = t(sapply(1:length(fSOMObject$map$medianValues[,1]),function(i){apply(fSOMObject$data[intersect(dataIndex,which(fSOMObject$map$mapping[,1] == i)),,drop=F],2,function(x){median(x)})}))
         PlotStars(fSOM4Plot,backgroundValues = as.factor(metaClustFactors), main=mainTitle)
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
    PlotStars(fSOM,backgroundValues = as.factor(metacl))
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

## To do:  put plotlabel in plotTreeSet; automatic naming of metaclusters (median or quartile); extract sub data from a single metacluster
