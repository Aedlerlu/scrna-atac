#合并多组数据
dirs=c('E:/prac/BC2/','E:/prac/BC21/')
names(dirs)=c('BC2','BC21')
kdo<-Read10X(data.dir = dirs)
scrnacom<-CreateSeuratObject(counts,min.cells=3,min.features=40,project = "sampl.1")
#依旧过滤线粒体基因
scrnacom[["percentage.mt"]]=PercentageFeatureSet(scrnacom,pattern="^MT-")
scrnacom <- subset(scrnacom,subset = nFeature_RNA>500&percentage.mt<20&nCount_RNA>1000)

#标准化
kdo<- NormalizeData(kdo)
kdo <- FindVariableFeatures(kdo)
kdo<-ScaleData(kdo)
kdo <- RunPCA(kdo)
#双细胞去除,这一步应该在标准化之后
#计算nexp
nEXP<-round(0.08*length(kdo@meta.data$orig.ident)*length(kdo@meta.data$orig.ident)/10000)
kdod <-doubletFinder(kdo,PCs = 1:20,nExp = nEXP,pK = 0.09,pN = 0.25,sct = FALSE)#没用sct就false,第一次跑不要加reuse.pann，会报错
#去除
kdod$doublet.class<- kdod[[paste0("DF.classifications_0.25_0.09_",nEXP)]]
kdod[[paste0("DF.classifications_0.25_0.09_",nEXP)]]<-NULL
pann<-grep(pattern = "^pANN",x=names(kdod@meta.data),value = TRUE)
kdod$pANN_0.25_0.09_nEXP<-kdod[[pann]]
kdod[[pann]] <- NULL
kdod <-subset(kdod,subset = doublet.class != "Doublet")
#开始整合，整合前要先将读在一块的样本分组
kdod[["RNA"]]<- split(kdod[["RNA"]],f = kdod$orig.ident)
ifnb<-IntegrateLayers(ifnb,method = CCAIntegration,orig.reduction = "pca",new.reduction="integrated.cca",verbose=T)
ifnb[["RNA"]]<-JoinLayers(ifnb[["RNA"]])
#整合后分析
ifnb<-FindNeighbors(ifnb,reduction = "integrated.cca")
ifnb<-FindClusters(ifnb,resolution = 1)
ifnb<-RunUMAP(ifnb,dims = 1:30,reduction = "integrated.cca")
#依旧画图
ployn1<-DimPlot(ifnb,reduction = "umap",group.by = c("orig.ident"))
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#harmony方法整合,完整pipline
dir<-list.dirs(path = "E:/prac/",full.names = TRUE,recursive = FALSE)
dir<-dir[-1]
dir <- c(dir[1],dir[2],dir[3],dir[4])
#分开处理每个样本
kdod_l<-list()
for (i in seq_along(dir)) {
  x<-Read10X(data.dir = dir[i])
  projectname<-paste0("sample.",i)
  x<-CreateSeuratObject(x,min.cells=3,min.features=40,project =projectname )
  x[["percentage.mt"]]=PercentageFeatureSet(x,pattern="^MT-")
  x<- subset(x,subset = nFeature_RNA>300&percentage.mt<20&nCount_RNA>1000)
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x,nFeatures=3000)
  x<-ScaleData(x)
  x<-RunPCA(x)
  nEXP<-round(0.08*length(x@meta.data$orig.ident)*length(x@meta.data$orig.ident)/10000)
  x <-doubletFinder(x,PCs = 1:30,nExp = nEXP,pK = 0.09,pN = 0.25,sct = FALSE)#没用sct就false,第一次跑不要加reuse.pann，会报错
  #去除
  x$doublet.class<- x[[paste0("DF.classifications_0.25_0.09_",nEXP)]]
  x[[paste0("DF.classifications_0.25_0.09_",nEXP)]]<-NULL
  pann<-grep(pattern = "^pANN",x=names(x@meta.data),value = TRUE)
  x[[paste0("paNN_0.25_0.09_",nEXP)]]<-x[[pann]]
  x[[pann]] <- NULL
  x <-subset(x,subset = doublet.class != "Doublet")
  kdod_l[[i]]<-x
  names(kdod_l)[i]<-paste0("c",i)
}
#features <- SelectIntegrationFeatures(object.list = kdod_l, nfeatures = 3000)
kdod_t<-merge(x=kdod_l[[1]],y=kdod_l[-1],add.cell.ids = names(kdod_l))
kdod_t[["RNA"]]<-JoinLayers(kdod_t[["RNA"]])
kdod_t<-FindVariableFeatures(kdod_t,nfeatures = 3000)
kdod_t<- ScaleData(kdod_t)
kdod_t<-RunPCA(kdod_t,npcs = 30)
#整合
kdod_t<-RunHarmony(object=kdod_t,group.by.vars="orig.ident",reduction.use="pca",dims=1:30)
kdod_t<-FindNeighbors(kdod_t,reduction="harmony")
kdod_t<-FindClusters(kdod_t,resolution= 0.6,cluster.name="harmony_cluster")
kdod_t<-RunUMAP(kdod_t,dims=1:30,reduction.name="umap.harmony",reduction = "harmony")
plot<-DimPlot(kdod_t,reduction = "umap.harmony",group.by = c("harmony_cluster","orig.ident"),label = T)
#注释，用findaallmarker
markers<-FindAllMarkers(kdod_t,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.5,test.use = "MAST")
top5<-markers%>%group_by(cluster)%>%top_n(n=5,wt = avg_log2FC)
#dotplot显示marker基因在不同cluster的表达情况,a_geneset_list为手动筛选的一些个细胞marker集的列表
DotPlot(kdod_t,features = a_geneset_list,assay = "RNA") + theme(axis.title.y = element_text(angle = 30,vjust = 0.5,hjust = 1),axis.title = element_blank())
#singler注释
celldex::HumanPrimaryCellAtlasData()
refdata<- HumanPrimaryCellAtlasData
HumanPrimaryCellAtlasData<-celldex::HumanPrimaryCellAtlasData()
refdata<- HumanPrimaryCellAtlasData
#提取layer
testdata<-LayerData(kdod_t,assay = "RNA",layer = "data")
cluster<-kdod_t@meta.data$seurat_clusters
#singler分析正式开始
cellpred<-SingleR(testdata,ref = refdata,clusters = cluster,labels = refdata$label.main,method = "cluster",assay.type.test = "logcounts",assay.type.ref = "logcounts")
celltype<-data.frame(ClusterID=rownames(cellpred),celltype=cellpred@listData$labels)
#baocun，写成csv

#将得到的注释文件返回到seurat对象
kdod_sr<-RenameIdents(kdod_t,c("0" = "Macrophage","1"="T_cells","2"="Chondrocytes","3"="Macrophage","4"="Tissue_stem_cells","5"="Endothelial_cells","6"="Monocyte","7"="Macrophage","8"="Macrophage","9"="MSC","10"="Chondrocytes","11"="Endothelial_cells"))
dp<-DimPlot(kdod_sr,label = T)
#还有方法二,写循环匹配替换
kdod_t@meta.data$celltype<-NA
  for (i in 1:nrow(celltype)) {
   ldx<- which(kdod_t@meta.data$seurat_clusters == celltype$ClusterID[i])
   kdod_t@meta.data$celltype[ldx]<-celltype$celltype[i]
  }
#三，match
kdod_t$celltype<-celltype$celltype[match(kdod_t@meta.data$seurat_clusters,celltype$ClusterID)]

#注释掉的另一种方式，用软件azimuth
#注释亚群,以某一细胞为例
chron<- subset(kdod_t,subset=celltype=="Chondrocytes")
chron<-NormalizeData(chron)
chron<-FindVariableFeatures(chron)
chron<-ScaleData(chron)
chron<-RunPCA(chron)
ts<-ElbowPlot(chron,ndims = 50)
chron<-FindNeighbors(chron,dims = 1:30)
chron<-FindClusters(chron,resolution = 1)
chron<-RunUMAP(chron,reduction.name = "umap",reduction = "pca",dims = 1:30)
chron_p<-DimPlot(chron,reduction = "umap",label = T)
clusterch<-chron@meta.data$seurat_clusters
#亚群注释建议用findallmarker，结合文献进行手动注释
markers<-FindAllMarkers(chron,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.5,test.use = "MAST")
top5<-markers%>%group_by(clusterch)%>%top_n(n=5,wt = avg_log2FC)