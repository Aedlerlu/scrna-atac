#多样本整合（包含doublet remove和常规标准化）
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(harmony)
library("DoubletFinder")
library(glmGamPoi)
#设置工作路径
setwd("E:/WOR/")
#读取数据
dirs<-list.dirs(path = "E:/prac/",recursive = FALSE,full.names = TRUE)
dirs<-subset(dirs,subset = dirs !=dirs[1])

kdod[["percentage.mt"]]=PercentageFeatureSet(kdod,pattern="^MT-")
kdod <- subset(kdod,subset = nFeature_RNA>500&percentage.mt<20&nCount_RNA>1000)
kdod <- SplitObject(kdod,split.by = "orig.ident")
normalize_and_rmdoublet<-function(x){

  nEXP<-round(0.08*length(x@meta.data$orig.ident)*length(x@meta.data$orig.ident)/10000)
  x <-doubletFinder(x,PCs = 1:20,nExp = nEXP,pK = 0.09,pN = 0.25,sct = FALSE)#没用sct就false,第一次跑不要加reuse.pann，会报错
  #去除
  x$doublet.class<- x[[paste0("DF.classifications_0.25_0.09_",nEXP)]]
  x[[paste0("DF.classifications_0.25_0.09_",nEXP)]]<-NULL
  pann<-grep(pattern = "^pANN",x=names(x@meta.data),value = TRUE)

  x[[paste0("paNN_0.25_0.09_",nEXP)]]<-x[[pann]]
  x[[pann]] <- NULL
  x <-subset(x,subset = doublet.class != "Doublet")
  return(x)
}
kdod_l<-lapply(kdod,normalize_and_rmdoublet)
kdod_l<-IntegrateLayers(kdod_l,normalization.method="SCT",method = CCAIntegration,verbose=T)
kdod_l<-ScaleData(kdod_l)
kdod_l<-RunPCA(kdod_l)
kdod_l<-FindNeighbors(kdod_l,reduction="integrated.dr")
kdod_l<-FindClusters(kdod_l,resolution= 0.6)
kdod_l<-RunUMAP(kdod_l,dims=1:30)