################ 
## VALIDATION 1: Topic Enrichment Analysis

## Load Gene Reference
LoadRef = function(){
  filename = "/Users/balubao/Documents/MATLAB/Projects/data/uniprot_wagner2018_1.tsv"
  GOref = read.table(filename, sep = "\t", header = TRUE, fill = TRUE)
  filename = "/Users/balubao/Documents/MATLAB/Projects/data/uniprot_wagner2018_2.tsv"
  GOref1 = read.table(filename, sep = "\t", header = TRUE, fill = TRUE)
  clabs = colnames(GOref)
  clabs[1] = "GOI"
  colnames(GOref) = colnames(GOref1) = clabs
  GOref = rbind(GOref,GOref1)
  return(GOref)
}
GOref = LoadRef()

### Set up functions and ref
geneRef = tolower(as.character(GOref$GOI))

TopicGO.bioprocess = function(InpTopicGenes){
  InpTopicGenes = tolower(InpTopicGenes)
  gene.select = geneRef %in% InpTopicGenes
  TopicGOnames = GOref$Gene.ontology..biological.process.[gene.select]
  TopicGOnames = sub("\\[GO","",lapply(strsplit(as.character(TopicGOnames),split = ":"),"[",1))
  TopicGOtbl = table(TopicGOnames)
  return(TopicGOtbl)}
TopicGO.cellularcomp = function(InpTopicGenes){
  InpTopicGenes = tolower(InpTopicGenes)
  gene.select = geneRef %in% InpTopicGenes
  TopicGOnames = GOref$Gene.ontology..cellular.component.[gene.select]
  TopicGOnames = sub("\\[GO","",lapply(strsplit(as.character(TopicGOnames),split = ":"),"[",1))
  TopicGOtbl = table(TopicGOnames)
  return(TopicGOtbl)}
TopicGO.molecularfunc = function(InpTopicGenes){
  InpTopicGenes = tolower(InpTopicGenes)
  gene.select = geneRef %in% InpTopicGenes
  TopicGOnames = GOref$Gene.ontology..molecular.function.[gene.select]
  TopicGOnames = sub("\\[GO","",lapply(strsplit(as.character(TopicGOnames),split = ":"),"[",1))
  TopicGOtbl = table(TopicGOnames)
  return(TopicGOtbl)}
TopicGO.GO = function(InpTopicGenes){
  InpTopicGenes = tolower(InpTopicGenes)
  gene.select = geneRef %in% InpTopicGenes
  TopicGOnames = GOref$Gene.ontology..GO.[gene.select]
  TopicGOnames = sub("\\[GO","",lapply(strsplit(as.character(TopicGOnames),split = ":"),"[",1))
  TopicGOtbl = table(TopicGOnames)
  return(TopicGOtbl)}

#GO labels
geneRef

TopicGOnames = GOref$Gene.ontology..biological.process.
TopicGOnames = sub("\\[GO","",lapply(strsplit(as.character(TopicGOnames),split = ":"),"[",1))
uniqueTopicGOnames = unique(TopicGOnames)

### Load Data
LoadTopicGeneList = function(filename){
  filename = sub(pattern = "..",replacement = "/Users/balubao/Documents/MATLAB/Projects",filename)
  # system(paste0("rm ",sub(pattern = "csv",replacement = "tsv",x = filename)," ",sub(pattern = "csv",replacement = "txt",x = filename)))
  df = read.delim(filename,
                  sep = "\t",
                  header=T, 
                  fill=T, 
                  as.is = T)
  # df = t.data.frame(df)
  # colnames(df) = df[1,]
  # df = df[-1,]
  return(df)
}

#Functions
MergeTopics = function(TopicTable){
  COL<-unique(names(TopicTable))
  ROW<-unique(unlist(lapply(TopicTable, rownames)))
  # Empty DF with all combinations
  TOTAL<-matrix(data=0, nrow=length(ROW), ncol=length(COL), dimnames=list(ROW, COL))
  # Subsetting :
  for (topic in c(1:length(TopicTable))) { 
    TOTAL[names(TopicTable[[topic]]),names(TopicTable[topic])] <- as.vector(TopicTable[[topic]])
  }
  return(TOTAL)}

filename = "/Users/balubao/Documents/MATLAB/Projects/results/filedirectory_load.txt"
filedir = read.delim(filename,header = F)
dataid = sapply(strsplit(basename(as.character(filedir$V1)),split = "\\."),"[",1)
TopicMatbp = TopicMatcc = TopicMatmf = TopicMatgo = list()
for(i in 1:length(dataid)){ #FOR ALL FILES IN DIRECTORY
  TopicGeneList = LoadTopicGeneList(as.character(filedir[i,]))
  TopicGObp=TopicGOcc=TopicGOmf=TopicGOgo=list()
  for(topic in 1:ncol(TopicGeneList)){
    TopicGObp[[paste0("Topic",topic)]] = TopicGO.bioprocess(TopicGeneList[,topic])
    TopicGOcc[[paste0("Topic",topic)]] = TopicGO.cellularcomp(TopicGeneList[,topic])
    TopicGOmf[[paste0("Topic",topic)]] = TopicGO.molecularfunc(TopicGeneList[,topic])
    TopicGOgo[[paste0("Topic",topic)]] = TopicGO.GO(TopicGeneList[,topic])
  }
  
  TopicMatbp[[dataid[i]]] = MergeTopics(TopicGObp)
  TopicMatcc[[dataid[i]]] = MergeTopics(TopicGOcc)
  TopicMatmf[[dataid[i]]] = MergeTopics(TopicGOmf)
  TopicMatgo[[dataid[i]]] = MergeTopics(TopicGOgo)
}





# Merge Topics
#libraries
library(ComplexHeatmap)
        
hgo2top.bp = hgo2top.mf = hgo2top.cc = hgo2top.go = list()  
for(i in 1:length(dataid)){
  hgo2top.bp[[dataid[i]]] = Heatmap(TopicMatbp[[dataid[i]]])
  hgo2top.cc[[dataid[i]]] = Heatmap(TopicMatcc[[dataid[i]]])
  hgo2top.mf[[dataid[i]]] = Heatmap(TopicMatmf[[dataid[i]]])
  hgo2top.go[[dataid[i]]] = Heatmap(TopicMatgo[[dataid[i]]])
}

## Visualize Topic GO
PlotTopicPlots = function(TopicGO){
  ntopics = length(TopicGO)
  if(ntopics > 16){
    ntopics = 16; par(mfrow = c(4,4))
  }else if(ntopics==10){
    par(mar=c(4,7,4,4), mfrow = c(5,2))
  }else if(ntopics==4){
    par(mfrow = c(2,2))}
  
  for(i in 1:ntopics){
    barplot(TopicGO[[i]],horiz = T, las=1, main = paste0("Topic ",i))
  }
}


i=1
PlotTopicPlots(TopicGOgo[[dataid[i]]])

TopicMatbp[[dataid[1]]][,2]

i=2
plt.data = as.data.frame(TopicMatbp[[dataid[i]]])
plt.data$GO = rownames(plt.data)
# ggplot(plt.data,aes(x=GO,y=Topic1))+
#   geom_col()+
#   coord_flip() +
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))

plt.data = melt(plt.data)
ggplot(plt.data,aes(x=GO,y=variable,fill=value))+
  geom_tile(colour = "grey50")+
  coord_flip() +
  theme_classic()+
  ggtitle(dataid[i])+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
 

# geom_density_ridges

# Heatmap(TopicMatbp[[dataid[i]]])

################ 
## VALIDATION 2: Clustering Analysis

  #libraries
  library(ggplot2)
  library(Rtsne)
  library("tidyr")
  library(reshape2)
  library(gridExtra)
  library(dbscan)
  library(umap)
library(seriation)
  
  #Functions
  LoadMatEmbeddings = function(filename){
    filename = as.character(filename)
    filename = sub(pattern = "..",replacement = "/Users/balubao/Documents/MATLAB/Projects",filename)
    mat = read.table(filename, sep = "\t", header = T, fill = TRUE)
    colnames(mat) = paste0("Topic",c(1:ncol(mat)))
    rownames(mat) = paste0("c",rownames(mat))
    return(mat)
  }
  
  ## Load Data (cell labels)
  filename = "/Users/balubao/Documents/MATLAB/Projects/results/matrices/cellnames.txt"
  CellLabels = read.table(filename,header = F, sep = "\t")
  # CellLabels = as.factor(CellLabels$V1)
  CellLabels = as.factor(unlist(lapply(strsplit(as.character(CellLabels$V1),split = " - "),"[",1)))
  
  
  
  ## Load Data (cell embeddings)
  filename = "/Users/balubao/Documents/MATLAB/Projects/results/filedirectory_mtx.txt"
  filedir = read.delim(filename,header = F)
  dataid = sapply(strsplit(basename(as.character(filedir$V1)),split = "\\."),"[",1)
  mat = list()
  # mat = mattsne = matumap = list()
  for(i in 1:length(dataid)){ #FOR ALL FILES IN DIRECTORY
    print(paste0("Loading ",dataid[i]," data "))
    mat[[dataid[i]]] = LoadMatEmbeddings(filedir[i,])
    # mattsne[[dataid[i]]] = Rtsne(mat[[i]],dims = 2)
    # matumap[[dataid[i]]] = umap(mat[[i]])
    
  }

  # jaccardIndex = intersect()/union()
  
  #Plot tSNE Embeddings
  ptsne = list()
  for(i in 1:length(dataid)){
    plt.data = as.data.frame(mattsne[[i]]$Y)
    ptsne[[i]] = ggplot(plt.data,aes(x=V1,y=V2))+
      geom_point(color=as.numeric(CellLabels), size=0.1)+
      ggtitle(dataid[i])+
      xlab("tSNE1")+
      ylab("tSNE2")+
      theme_classic()
  }

  #Plot tSNE plots
  grid.arrange(grobs= p[1:10],nrow=5,ncol=2)
  grid.arrange(grobs= p[11:20],nrow=5,ncol=2)
  grid.arrange(grobs= p[21:30],nrow=5,ncol=2)

  
  #Plot UMAP Embeddings
  pumap = list()
  for(i in 1:length(dataid)){
  plt.data = as.data.frame(matumap[[i]]$layout)
  pumap[[i]] = ggplot(plt.data,aes(x=V1,y=V2))+
    geom_point(color=as.numeric(CellLabels), size=0.1)+
    ggtitle(dataid[i])+
    xlab("UMAP1")+
    ylab("UMAP2")+
    theme_classic()
  }
  
  #Plot UMAP plots
  grid.arrange(grobs= pumap[1:10],nrow=5,ncol=2)
  grid.arrange(grobs= pumap[11:20],nrow=5,ncol=2)
  grid.arrange(grobs= pumap[21:30],nrow=5,ncol=2)


ggplot(plt.data,aes(x=V1,y=V2))+
  geom_point(color=DBSCAN_clust$cluster+1, size=0.1)+
  ggtitle(dataid[i])+
  xlab("tSNE1")+
  ylab("tSNE2")+
  theme_classic()

  DBSCAN_clust = dbscan(mattsne[[1]]$Y,eps = 0.8,minPt = 5)
  ggplot(plt.data,aes(x=V1,y=V2))+
    geom_point(color=DBSCAN_clust$cluster+1, size=0.1)+
    ggtitle(dataid[i])+
    xlab("tSNE1")+
    ylab("tSNE2")+
    theme_classic()
  
  p = list()
  for(i in 1:length(dataid)){
    plt.data = as.data.frame(mattsne[[i]]$Y)
    p[[i]] = ggplot(plt.data,aes(x=V1,y=V2))+
      geom_point(color=as.numeric(CellLabels), size=0.1)+
      ggtitle(dataid[i])+
      xlab("tSNE1")+
      ylab("tSNE2")+
      theme_classic()
  }
  
  
  
  ggplot(plt.data,aes(x=CellLabels))+
    geom_boxplot(aes(y=Topic1))+
    geom_boxplot(aes(y=Topic2))
  
  cor(mat[[1]],as.numeric(CellLabels))
  plot(CellLabels,mat$Topic1,las=2)
  
  plt.data = cbind(CellLabels,mat[[1]])
  # plt.data = reshape(plt.data,ids = "CellLabels", varying = c(2:ncol(mat)),direction = "long")
  
  # Topic enrichment in cell-types (Heatmap)
  
  # Heatmap(as.matrix(plt.data[,2:5]))
  
  ggplot(plt.data,aes(y=as.numeric(CellLabels),x=as.numeric(value),fill=variable))+
    geom_raster()+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  # Topic enrichment in cell-types (Boxplot)
  plt.data2 = melt(plt.data)
  ggplot(plt.data2,aes(x=CellLabels,y=value,color=variable))+
    geom_boxplot()+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
   
  #Get Mean Topic Value per Celltype
  neomat = list() 
  for(i in 1:length(dataid)){
    neomat_part = c()
    for(j in c(1:nlevels(CellLabels))){
      select.c = which(as.numeric(CellLabels) == j)
      neomat_row = colMeans(mat[[i]][select.c,], na.rm = T)
      neomat_part = rbind(neomat_part,neomat_row)
    }
    rownames(neomat_part) = levels(CellLabels)
    neomat[[dataid[i]]] = neomat_part
  }
  
  
  for(){
    
  }
  i=30
  multout = neomat[[dataid[i]]]%*%t(TopicMatgo[[dataid[i]]])
  
  o = seriate(max(neomat[[dataid[i]]]) - neomat[[dataid[i]]], method = "BEA_TSP")
  Heatmap(max(neomat[[dataid[i]]]) - neomat[[dataid[i]]], name = "mat", 
          row_order = get_order(o, 1), column_order = get_order(o, 2))
  
  o = seriate(max(neomat[[dataid[i]]]) - neomat[[dataid[i]]], method = "BEA_TSP")
  Heatmap(max(neomat[[dataid[i]]]) - neomat[[dataid[i]]], name = "mat", 
          row_order = get_order(o, 1), column_order = get_order(o, 2))
  
  Heatmap_seriated = function(mat,Name){
    # matrep = mat - min(mat)
    # o = seriate(matrep, method = "BEA_TSP")
    Heatmap(mat,name = Name,column_title = Name)
            # row_order = get_order(o, 1), column_order = get_order(o, 2),column_title = Name)
  }
  Heatmap_seriated2 = function(mat,Name){
    # matrep = mat - min(mat)
    # o = seriate(matrep, method = "BEA_TSP")
    Heatmap(mat,name = Name,row_title = Name)
    # row_order = get_order(o, 1), column_order = get_order(o, 2),column_title = Name)
  }

  i=10
  multout = neomat[[dataid[i]]]%*%t(TopicMatgo[[dataid[i]]])
  ht1.LDA = Heatmap_seriated(neomat[[dataid[i]]],dataid[i])
  ht2.LDA = Heatmap_seriated2(TopicMatgo[[dataid[i]]],dataid[i])
  htmult.LDA = Heatmap_seriated(multout,dataid[i])
  
  i=20
  multout = neomat[[dataid[i]]]%*%t(TopicMatgo[[dataid[i]]])
  ht1.LSA = Heatmap_seriated(neomat[[dataid[i]]],dataid[i])
  ht2.LSA = Heatmap_seriated2(TopicMatgo[[dataid[i]]],dataid[i])
  htmult.LSA = Heatmap_seriated(multout,dataid[i])
  
  i=30
  multout = neomat[[dataid[i]]]%*%t(TopicMatgo[[dataid[i]]])
  ht1.NMF = Heatmap_seriated(neomat[[dataid[i]]],dataid[i])
  ht2.NMF = Heatmap_seriated2(TopicMatgo[[dataid[i]]],dataid[i])
  htmult.NMF = Heatmap_seriated(multout,dataid[i])
  
  ht1.LDA+ht1.LSA+ht1.NMF
  ht2.LDA%v%ht2.LSA%v%ht2.NMF
  
  htmult.LDA
  
  # ht = 

   ht1 =  Heatmap(t(neomat[[1]]))
  ht2 = Heat
  
  
  NormalizeMatrix = function(mat){
    #minmax
    
  }
  
###################
  ## Heirarichal clustering of cell-types using topics
  standardize <- function(x){(x-min(x))/(max(x)-min(x))}

# By latent topics  
  i=30
  d <- dist(neomat[[dataid[i]]])   # find distance matrix 
  hc <- hclust(d)                # apply hirarchical clustering 
  plot(hc, main=dataid[i])                       # plot the dendrogram
 
# By PCA
  matPCA = data.seurat@reductions$pca@cell.embeddings[,1:20]
  neomatPCA = c()
  for(j in c(1:nlevels(CellLabels))){
    select.c = which(as.numeric(CellLabels) == j)
    neomat_row = colMeans(matPCA[select.c,], na.rm = T)
    neomatPCA = rbind(neomatPCA,neomat_row)
  }
  rownames(neomatPCA) = levels(CellLabels)
  
  d <- dist(neomatPCA)   # find distance matrix 
  hc <- hclust(d)                # apply hirarchical clustering 
  plot(hc, main="PCA20")                       # plot the dendrogram
  
  # clusterCut <- cutree(hc, 3)
  # Compute distance matrix
  res.dist <- dist(df, method = "euclidean")
  
  # Compute 2 hierarchical clusterings
  hcLDA <- hclust(res.dist, method = "average")
  hcLSA <- hclust(res.dist, method = "ward.D2")
  hcNMF <- hclust(res.dist, method = "ward.D2")
  
  # Create two dendrograms
  dend1 <- as.dendrogram (hc1)
  dend2 <- as.dendrogram (hc2)
  
  # Create a list to hold dendrograms
  dend_list <- dendlist(dend1, dend2)
  
  # Align and plot two dendrograms side by side
  dendlist(dend1, dend2) %>%
    untangle(method = "step1side") %>% # Find the best alignment layout
    tanglegram()                       # Draw the two dendrograms
  # Compute alignment quality. Lower value = good alignment quality
  dendlist(dend1, dend2) %>%
    untangle(method = "step1side") %>% # Find the best alignment layout
    entanglement()
  # Cophenetic correlation matrix
  cor.dendlist(dend_list, method = "cophenetic")
  
  # Create multiple dendrograms by chaining
  library(dendextend)
  dendmethod = "average"
  dend1 <- neomat[[dataid[10]]] %>% dist %>% hclust(dendmethod) %>% as.dendrogram
  dend2 <- neomat[[dataid[20]]] %>% dist %>% hclust(dendmethod) %>% as.dendrogram
  dend3 <- neomat[[dataid[30]]] %>% dist %>% hclust(dendmethod) %>% as.dendrogram
  dend4 <- neomatPCA %>% dist %>% hclust(dendmethod) %>% as.dendrogram
  # Compute correlation matrix
  dend_list <- dendlist("LDA20" = dend1, "LSA20" = dend2,
                        "NMF20" = dend3, "PCA20" = dend4)
  
  
  # dend_list <- dendlist("Complete" = dend1, "Single" = dend2,
                        # "Average" = dend3, "Centroid" = dend4)
  dendcors <- cor.dendlist(dend_list)
  # Print correlation matrix
  round(dendcors, 2)
  
  # Visualize the correlation matrix using corrplot package
  library(corrplot)
  corrplot(dendcors, "pie", "lower")
  title(dendmethod)
  
  
  #Ground truth
  
  
  # ggtitle(label = "Gene-Gene Pearson Correlation")
  
  #retrive embedding
  #Construct tSNE on embedding
  #Density clustering
  #Jaccard of real label?



# write.table(df,gsub(pattern = "csv",replacement = "tsv",x = filename),
#             sep = "\t",
#             row.names = F)
# system(paste0("sed 's/\"//g' ",gsub(pattern = "csv",replacement = "tsv",x = filename)," > ",gsub(pattern = "csv",replacement = "txt",x = filename)))


################ 
## State-of-the-Art 
    
    library(R.matlab)
    library(Seurat)

    filename = "/Users/balubao/Documents/MATLAB/Projects/data/DataPartition3.mat"
    InpData = readMat(filename) 
    
    #Plot PCA based tSNE and UMAP
    data.X = InpData$X
    colnames(data.X) = paste0("c",c(1:ncol(data.X)))
    gene.names = read.delim(file = "/Users/balubao/Documents/MATLAB/Projects/data/Wagner2018_GeneNames.txt",header = F)
    rownames(data.X) = gene.names$V1
    data.seurat = CreateSeuratObject(counts = data.X)
    data.seurat$CellLabels = CellLabels
    #Run Esentials to UMAP
    data.seurat = NormalizeData(data.seurat)
    data.seurat = ScaleData(data.seurat)
    data.seurat = FindVariableFeatures(data.seurat)
    data.seurat = RunPCA(data.seurat, npcs = 30)
    ElbowPlot(data.seurat, ndims=30)
    data.seurat = RunTSNE(data.seurat, dims = 1:10, n.components = 2L)
    data.seurat = RunUMAP(data.seurat, dims = 1:10, n.components = 2L)
    
    data.seurat = FindNeighbors(data.seurat, dims = 1:10)
    data.seurat = FindClusters(data.seurat, resolution = 0.5)
    plotClusterTree(data.seurat)
    
    pseurat.tsne = DimPlot(data.seurat, reduction = "tsne" ,dims = c(1:2), group.by = "CellLabels")#+theme(legend.position="bottom")
    pseurat.umap = DimPlot(data.seurat, reduction = "umap" ,dims = c(1:2), group.by = "CellLabels")#+theme(legend.position="bottom")
    
    DimPlot(data.seurat, reduction = "umap" ,dims = c(1:2), group.by = "CellLabels")
    DimPlot(data.seurat, reduction = "umap" ,dims = c(1:2), group.by = "seurat_clusters")    
    
    
    ##
    plotting.datatsne = data.frame(TSNE1 = data.seurat@reductions$tsne@cell.embeddings[,1],
                                   TSNE2 = data.seurat@reductions$tsne@cell.embeddings[,2])
    
    i=10
    plotting.data = cbind(plotting.datatsne[-1,],mat[[dataid[i]]])
    ptop = list()
    #length(mat[[dataid[i]]])){
    # topiclabs = colnames(mat[[dataid[i]]])
    for(j in c(1:4)){
      ptop[[j]]=ggplot(plotting.data)+
        geom_point(aes(x=TSNE1,y=TSNE2,colour=plotting.data[,j+2]), size=0.5)+
        theme_classic()
    }
     
    grid.arrange(grobs= ptop[c(1:4)],nrow=2,ncol=2)
    
    