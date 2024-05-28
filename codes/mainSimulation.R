# delete work space
rm(list = ls(all = TRUE))
graphics.off()

library(data.table)
library(purrr)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(MOFA2)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(cowplot)
library(caret)
library(ggpubr)
library(corrplot)
library(geigen)
library(MASS)
library(Rtsne)
library(umap)
library(pheatmap)
library(RGCCA)
library(ConsensusClusterPlus)
library(mclust)
library(aricode)
library(omicade4)
library(IntNMF)
library(r.jive)
library(tensorBSS)
library(mixOmics) # import the mixOmics library
library(PANDA)

source("C:/Projects/mixOmics/Data/SNF data/Breast/tICA.R")
source("C:/Projects/mixOmics/Data/SNF data/Breast/iCluster2.R")


library(reticulate)

load("C:/Users/MAminu/OneDrive - Inside MD Anderson/Projects/DOmics working draft/Final/Nature Communications/PANDA/data/Simulation/data_simulation.rda")

data.rna <- CreateSeuratObject(counts = as.data.frame(data_simulation[["data"]][["RNA"]]),
                               project = "simProject",
                               assay = "RNA",
                               meta.data = data_simulation[["labels"]])

data.atac <- CreateSeuratObject(counts = as.data.frame(data_simulation[["data"]][["ATAC"]]),
                                project = "simProject",
                                assay = "ATAC",
                                meta.data = data_simulation[["labels"]])


# Perform standard analysis of each modality independently RNA analysis
data.rna <- NormalizeData(data.rna)
data.rna <- FindVariableFeatures(data.rna,nfeatures = 2000)
data.rna <- ScaleData(data.rna)
data.rna <- RunPCA(data.rna)
data.rna <- RunUMAP(data.rna, dims = 1:30)

# We exclude the first dimension as this is typically correlated with sequencing depth
data.atac <- RunTFIDF(data.atac)
data.atac <- FindTopFeatures(data.atac, min.cutoff = "q0")
data.atac <- RunSVD(data.atac,n=30)
data.atac <- RunUMAP(data.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
data.atac <- NormalizeData(data.atac)
data.atac <- ScaleData(data.atac, features = rownames(data.atac))


# select highly discriminant features in RNA
rna = t(as.data.frame(data.rna@assays[["RNA"]]@scale.data))
rnapcor <- data.frame(matrix(ncol = 5, nrow = ncol(rna)))
colnames(rnapcor) <- c('features', 'pvalue', 'corcoef', 'scores', 'type of association')
rnapcor$features <- colnames(rna)

gnd <- as.numeric(as.factor(data.rna@meta.data[["group"]]))
for (i in 1:ncol(rna)) {
  cor <- cor.test(rna[,i] ,gnd)
  rnapcor[i,2] <- cor[["p.value"]]
  rnapcor[i,3] <- cor[["estimate"]][["cor"]]
}

rnapcor <- rnapcor[order(rnapcor$pvalue),]

# feature ranking
rnaFRank = rnapcor[c(1:30),]
for (i in 1:30){
  rnaFRank[i,4] <- -log(rnaFRank[i,2])
  if(rnaFRank[i,3] > 0){
    rnaFRank[i,5] <- "positive"
  }
  else{
    rnaFRank[i,5] <- "negative"
  }
}

ggdotchart(rnaFRank, x = "features", y = "scores",
           color = "type of association",
           palette = c("#FF0000", "#1B9E77"),
           sorting = "descending",
           add = "segments",
           add.params = list(color = "darkgray", size = 4),
           rotate = TRUE,
           dot.size = 6,
           #label = round(fscImaging$Scores,1),
           #font.label = list(color = "white", size = 0,
           #                  vjust = 0.4),
           ggtheme = theme_pubr()
)+
  font("x.text", size = 18, color = "black", face = "bold.italic")+
  font("y.text", size = 18, color = "black", face = "bold.italic")+
  font("xy", size = 18, color = "black", face = "bold.italic")

# select highly discriminant features in ATAC
atac <- t(as.data.frame(data.atac@assays[["ATAC"]]@scale.data))
atacpcor <- data.frame(matrix(ncol = 5, nrow = ncol(atac)))
colnames(atacpcor) <- c('features', 'pvalue', 'corcoef', 'scores', 'type of association')
atacpcor$features <- colnames(atac)

for (i in 1:ncol(atac)) {
  cor <- cor.test(atac[,i] ,gnd)
  atacpcor[i,2] <- cor[["p.value"]]
  atacpcor[i,3] <- cor[["estimate"]][["cor"]]
}

atacpcor <- atacpcor[order(atacpcor$pvalue),]

# feature ranking
atacFRank = atacpcor[c(1:30),]
for (i in 1:30){
  atacFRank[i,4] <- -log(atacFRank[i,2])
  if(atacFRank[i,3] > 0){
    atacFRank[i,5] <- "positive"
  }
  else{
    atacFRank[i,5] <- "negative"
  }
}

ggdotchart(atacFRank, x = "features", y = "scores",
           color = "type of association",
           palette = c("#FF0000", "#1B9E77"),
           sorting = "descending",
           add = "segments",
           add.params = list(color = "darkgray", size = 4),
           rotate = TRUE,
           dot.size = 6,
           #label = round(fscImaging$Scores,1),
           #font.label = list(color = "white", size = 0,
           #                  vjust = 0.4),
           ggtheme = theme_pubr()
)+
  font("x.text", size = 18, color = "black", face = "bold.italic")+
  font("y.text", size = 18, color = "black", face = "bold.italic")+
  font("xy", size = 18, color = "black", face = "bold.italic")


### compare and save performance of different methods
reducedRNA <- rna[,c(rnapcor$features[1:1500])]
reducedATAC <- atac[,c(atacpcor$features[1:3000])]

ClustPerf <- data.frame(matrix(ncol = 2, nrow = 16))
colnames(ClustPerf) <- c('ARI', 'NMI')
rownames(ClustPerf) <- c("iCluster","intNMF","JIVE","MCIA","MOFA","MEFISTO","tICA_RNA","tICA_ATAC","RGCCA_RNA","RGCCA_ATAC",
                         "DIABLO_RNA","DIABLO_ATAC","seurat_RNA","seurat_ATAC","PANDA_RNA","PANDA_ATAC")

numComponents <- 10
omics <- list(rna = as.matrix(t(reducedRNA)), atac = as.matrix(t(reducedATAC)))
omics_pos<-list()
for(j in 1:length(omics)){
  if(min(omics[[j]])<0){
    omics_pos[[j]]<-omics[[j]]+abs(min(omics[[j]]))
  }else{
    omics_pos[[j]]<-omics[[j]]
  }
  omics_pos[[j]]<-omics_pos[[j]]/max(omics_pos[[j]])
}
###icluster
factorizations_icluster<-iCluster2(lapply(omics, function(x) t(x)), k=5)
factors_icluster<-as.matrix(t(factorizations_icluster$expZ))

iClusterRed <- CreateDimReducObject(
  embeddings = as.matrix(t(factorizations_icluster$expZ)),
  loadings = matrix(),
  stdev = numeric(),
  key = "iCluster_",
  assay = "iCluster"
)

data.rna@reductions[["iCluster"]] <- iClusterRed
rownames(data.rna@reductions[["iCluster"]]@cell.embeddings) <- colnames(data.rna@assays[["RNA"]])

cols <- c("#9F8F12","#027B8E","#E58601","#4E9F50","#A0522D")
DimPlot(data.rna, reduction = "iCluster", group.by = "group", label = FALSE, pt.size = 4, cols = cols)

iclusterclust = ConsensusClusterPlus(
  t(factors_icluster[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

subgroup <- 5
iclusterlabel <- unlist(iclusterclust[[subgroup]]["consensusClass"])
iclusterclustres <- clustComp(gnd, iclusterlabel)

ClustPerf[1,1] <- iclusterclustres[["ARI"]]
ClustPerf[1,2] <- iclusterclustres[["NMI"]]

###intNMF
factorizations_intnmf<-nmf.mnnals(dat=lapply(omics_pos, function(x) t(x)), k=5)
factors_intNMF<-as.matrix(factorizations_intnmf$W)

intNMFRed <- CreateDimReducObject(
  embeddings = as.matrix(factorizations_intnmf$W),
  loadings = matrix(),
  stdev = numeric(),
  key = "intNMF_",
  assay = "intNMF"
)

data.rna@reductions[["intNMF"]] <- intNMFRed
rownames(data.rna@reductions[["intNMF"]]@cell.embeddings) <- colnames(data.rna@assays[["RNA"]])

DimPlot(data.rna, reduction = "intNMF", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend()

intNMFclust = ConsensusClusterPlus(
  t(factors_intNMF[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

intNMFlabel <- unlist(intNMFclust[[subgroup]]["consensusClass"])
intNMFclustres <- clustComp(gnd, intNMFlabel)

ClustPerf[2,1] <- intNMFclustres[["ARI"]]
ClustPerf[2,2] <- intNMFclustres[["NMI"]]

### JIVE
factorizations_jive<-jive(omics, rankJ=numComponents, rankA = rep(numComponents, length(omics)), method = "given",
                          conv = "default", maxiter = 100, showProgress=FALSE)
rankJV <- factorizations_jive$rankJ
rankIV.v <- factorizations_jive$rankA
J<-numeric(0)
ng<-0
metagenes_jive <- list()
for(j in 1:length(omics)){
  J <- rbind(J,factorizations_jive$joint[[j]])
  ng<-c(ng,dim(factorizations_jive$joint[[j]])[1])
}
svd.o <- svd(J)
jV <- svd.o$v %*% diag(svd.o$d)
for(j in 1:length(omics)){
  metagenes_jive[[j]] <- svd.o$u[(1+sum(ng[1:j])):sum(ng[1:j+1]),1:rankJV]; ###error in dimension
  rownames(metagenes_jive[[j]])<-rownames(omics[[j]])
  colnames(metagenes_jive[[j]])<-1:numComponents
}
factors_jive=jV[,1:rankJV]

jiveRed <- CreateDimReducObject(
  embeddings = as.matrix(factors_jive),
  loadings = matrix(),
  stdev = numeric(),
  key = "jive_",
  assay = "jive"
)

data.rna@reductions[["jive"]] <- jiveRed
rownames(data.rna@reductions[["jive"]]@cell.embeddings) <- colnames(data.rna@assays[["RNA"]])

DimPlot(data.rna, reduction = "jive", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend()

jiveclust = ConsensusClusterPlus(
  t(factors_jive[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

jivelabel <- unlist(jiveclust[[subgroup]]["consensusClass"])
jiveclustres <- clustComp(gnd, jivelabel)

ClustPerf[3,1] <- jiveclustres[["ARI"]]
ClustPerf[3,2] <- jiveclustres[["NMI"]]

### MCIA
factorizations_mcia<-mcia(omics_pos, cia.nf = numComponents)
factors_mcia<-as.matrix(factorizations_mcia$mcoa$SynVar)

mciaRed <- CreateDimReducObject(
  embeddings = as.matrix(factorizations_mcia$mcoa$SynVar),
  loadings = as.matrix(factorizations_mcia$mcoa$axis[1:dim(omics[[j]])[1],]),
  stdev = numeric(),
  key = "MCIA_",
  assay = "mcia"
)

data.rna@reductions[["mcia"]] <- mciaRed
rownames(data.rna@reductions[["mcia"]]@cell.embeddings) <- colnames(data.rna@assays[["RNA"]])

DimPlot(data.rna, reduction = "mcia", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend()

mciaclust = ConsensusClusterPlus(
  t(factors_mcia[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

mcialabel <- unlist(mciaclust[[subgroup]]["consensusClass"])
mciaclustres <- clustComp(gnd, mcialabel)

ClustPerf[4,1] <- mciaclustres[["ARI"]]
ClustPerf[4,2] <- mciaclustres[["NMI"]]

### MOFA
mofa <- create_mofa(list(
  "RNA" = as.matrix(t(reducedRNA)),
  "ATAC" = as.matrix(t(reducedATAC))
))

mofa_opts <- get_default_model_options(mofa)
mofa_opts$num_factors <- numComponents

# Training options: let's use default options
train_opts <- get_default_training_options(mofa)
# train_opts$seed <- 42

mofa <- prepare_mofa(
  object = mofa,
  model_options = mofa_opts,
  training_options = train_opts
)

mofa <- run_mofa(mofa)

factors_mofa <- mofa@expectations[["Z"]][["group1"]]

# Add sample metadata to the model
# metaData$sample <- paste0("sample",c(1:206))
# samples_metadata(mofa) <- metaData

mofaRed <- CreateDimReducObject(
  embeddings = mofa@expectations[["Z"]][["group1"]],
  loadings = matrix(),
  stdev = numeric(),
  key = "MOFA_",
  assay = "mofa"
)

data.rna@reductions[["mofa"]] <- mofaRed
rownames(data.rna@reductions[["mofa"]]@cell.embeddings) <- colnames(data.rna@assays[["RNA"]])

DimPlot(data.rna, reduction = "mofa", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend()

mofaclust = ConsensusClusterPlus(
  t(factors_mofa[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

mofalabel <- unlist(mofaclust[[subgroup]]["consensusClass"])
mofaclustres <- clustComp(gnd, mofalabel)

ClustPerf[5,1] <- mofaclustres[["ARI"]]
ClustPerf[5,2] <- mofaclustres[["NMI"]]

### MEFISTO
mefisto <- create_mofa(list(
  "RNA" = as.matrix(t(reducedRNA)),
  "ATAC" = as.matrix(t(reducedATAC))
))

Princomp <- prcomp(reducedRNA,rank. = 50,center = TRUE,scale. = TRUE)
umapcomp <- umap(Princomp[["x"]])
mefisto <- set_covariates(mefisto, t(umapcomp[["layout"]]))

mefisto_opts <- get_default_model_options(mefisto)
mefisto_opts$num_factors <- numComponents

mefisto <- prepare_mofa(
  object = mefisto,
  model_options = mefisto_opts,
)

mefisto <- run_mofa(mefisto)
factors_mefisto <- mefisto@expectations[["Z"]][["group1"]]

# Add sample metadata to the model
# metaData$sample <- paste0("sample",c(1:206))
# samples_metadata(mefisto) <- metaData

mefistoRed <- CreateDimReducObject(
  embeddings = factors_mefisto,
  loadings = matrix(),
  stdev = numeric(),
  key = "MEFISTO_",
  assay = "mefisto"
)

data.rna@reductions[["mefisto"]] <- mefistoRed
rownames(data.rna@reductions[["mefisto"]]@cell.embeddings) <- colnames(data.rna@assays[["RNA"]])

DimPlot(data.rna, reduction = "mefisto", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend()

mefistoclust = ConsensusClusterPlus(
  t(factors_mefisto[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

mefistolabel <- unlist(mefistoclust[[subgroup]]["consensusClass"])
mefistoclustres <- clustComp(gnd, mefistolabel)

ClustPerf[6,1] <- mefistoclustres[["ARI"]]
ClustPerf[6,2] <- mefistoclustres[["NMI"]]

### tICA
omics_tensor<-list()
for(j in 1:length(omics)){
  omics_tensor[[j]]<-cor(omics[[j]], method = "spearman")
}

S<-vector(length = dim(omics[[1]])[2]*dim(omics[[1]])[2]*length(omics))
dim(S) <- c(length(omics), dim(omics[[1]])[2], dim(omics[[1]])[2])
for(j in 1:length(omics)){
  S[j,,]<-t(omics_tensor[[j]])
}
tICA<-DoTICA(S,numComponents,method="FOBI")
# factors_tica<-as.matrix(tICA$signals)
factors_ticaRNA <- as.matrix(t(tICA[["projS"]][[1]]))
factors_ticaATAC <- as.matrix(t(tICA[["projS"]][[2]]))

tICARNA <- CreateDimReducObject(
  embeddings = factors_ticaRNA,
  loadings = matrix(),
  stdev = numeric(),
  key = "tICA_",
  assay = "tICA"
)

data.rna@reductions[["tICA"]] <- tICARNA
rownames(data.rna@reductions[["tICA"]]@cell.embeddings) <- colnames(data.rna@assays[["RNA"]])

tICAATAC <- CreateDimReducObject(
  embeddings = factors_ticaATAC,
  loadings = matrix(),
  stdev = numeric(),
  key = "tICA_",
  assay = "tICA"
)

data.atac@reductions[["tICA"]] <- tICAATAC
rownames(data.atac@reductions[["tICA"]]@cell.embeddings) <- colnames(data.atac@assays[["ATAC"]])


p1 <- DimPlot(data.rna, reduction = "tICA", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(data.atac, reduction = "tICA", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend() + ggtitle("ATAC")
p1 + p2

ticaRNAclust = ConsensusClusterPlus(
  t(factors_ticaRNA[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

ticaRNAlabel <- unlist(ticaRNAclust[[subgroup]]["consensusClass"])
ticaRNAclustres <- clustComp(gnd, ticaRNAlabel)

ClustPerf[7,1] <- ticaRNAclustres[["ARI"]]
ClustPerf[7,2] <- ticaRNAclustres[["NMI"]]

ticaATACclust = ConsensusClusterPlus(
  t(factors_ticaATAC[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

ticaATAClabel <- unlist(ticaATACclust[[subgroup]]["consensusClass"])
ticaATACclustres <- clustComp(gnd, ticaATAClabel)

ClustPerf[8,1] <- ticaATACclustres[["ARI"]]
ClustPerf[8,2] <- ticaATACclustres[["NMI"]]

### RGCCA
factorizations_RGCCA<-rgcca(lapply(omics, function(x) t(x)), ncomp = rep(numComponents, length(omics)), scheme = "centroid", scale = TRUE, init = "svd",bias = TRUE, tol = 1e-08, verbose = F)
factors_rgccaRNA <- factorizations_RGCCA[["Y"]][[1]]
factors_rgccaATAC<- factorizations_RGCCA[["Y"]][[2]]

rgccaRNA <- CreateDimReducObject(
  embeddings = factorizations_RGCCA[["Y"]][[1]],
  loadings = factorizations_RGCCA[["a"]][[1]],
  stdev = numeric(),
  key = "RGCCA_",
  assay = "rgcca"
)

data.rna@reductions[["rgcca"]] <- rgccaRNA
rownames(data.rna@reductions[["rgcca"]]@cell.embeddings) <- colnames(data.rna@assays[["RNA"]])

rgccaATAC <- CreateDimReducObject(
  embeddings = factorizations_RGCCA[["Y"]][[2]],
  loadings = factorizations_RGCCA[["a"]][[2]],
  stdev = numeric(),
  key = "RGCCA_",
  assay = "rgcca"
)

data.atac@reductions[["rgcca"]] <- rgccaATAC
rownames(data.atac@reductions[["rgcca"]]@cell.embeddings) <- colnames(data.atac@assays[["ATAC"]])

p1 <- DimPlot(data.rna, reduction = "rgcca", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(data.atac, reduction = "rgcca", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend() + ggtitle("ATAC")
p1 + p2

rgccaRNAclust = ConsensusClusterPlus(
  t(factors_rgccaRNA[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

rgccaRNAlabel <- unlist(rgccaRNAclust[[subgroup]]["consensusClass"])
rgccaRNAclustres <- clustComp(gnd, rgccaRNAlabel)

ClustPerf[9,1] <- rgccaRNAclustres[["ARI"]]
ClustPerf[9,2] <- rgccaRNAclustres[["NMI"]]

rgccaATACclust = ConsensusClusterPlus(
  t(factors_rgccaATAC[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

rgccaATAClabel <- unlist(rgccaATACclust[[subgroup]]["consensusClass"])
rgccaATACclustres <- clustComp(gnd, rgccaATAClabel)

ClustPerf[10,1] <- rgccaATACclustres[["ARI"]]
ClustPerf[10,2] <- rgccaATACclustres[["NMI"]]

### DIABLO
data <- list(rna = as.matrix(reducedRNA), atac = as.matrix(reducedATAC))

diablo <- block.plsda(data, gnd, ncomp = numComponents, design = "null") # run the DIABLO method
factors_diabloRNA <- as.matrix(diablo[["variates"]][["rna"]])
factors_diabloATAC <- as.matrix(diablo[["variates"]][["atac"]])

DIABLORNA <- CreateDimReducObject(
  embeddings = diablo[["variates"]][["rna"]],
  loadings = diablo[["loadings"]][["rna"]],
  stdev = numeric(),
  key = "DIABLO_",
  assay = "RNA"
)

data.rna@reductions[["DIABLO"]] <- DIABLORNA
rownames(data.rna@reductions[["DIABLO"]]@cell.embeddings) <- colnames(data.rna@assays[["RNA"]])
# momix.rna <- RunTSNE(momix.rna, dims = 1:10, reduction = "DIABLO")

DIABLOATAC <- CreateDimReducObject(
  embeddings = diablo[["variates"]][["atac"]],
  loadings = diablo[["loadings"]][["atac"]],
  stdev = numeric(),
  key = "DIABLO_",
  assay = "ATAC"
)

data.atac@reductions[["DIABLO"]] <- DIABLOATAC
rownames(data.atac@reductions[["DIABLO"]]@cell.embeddings) <- colnames(data.atac@assays[["ATAC"]])
# momix.atac <- RunTSNE(momix.atac, dims = 1:10, reduction = "DIABLO")

p1 <- DimPlot(data.rna, reduction = "DIABLO", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(data.atac, reduction = "DIABLO", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend() + ggtitle("ATAC")
p1 + p2

diabloRNAclust = ConsensusClusterPlus(
  t(factors_diabloRNA[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

diabloRNAlabel <- unlist(diabloRNAclust[[subgroup]]["consensusClass"])
diabloRNAclustres <- clustComp(gnd, diabloRNAlabel)

ClustPerf[11,1] <- diabloRNAclustres[["ARI"]]
ClustPerf[11,2] <- diabloRNAclustres[["NMI"]]


diabloATACclust = ConsensusClusterPlus(
  t(factors_diabloATAC[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

diabloATAClabel <- unlist(diabloATACclust[[subgroup]]["consensusClass"])
diabloATACclustres <- clustComp(gnd, diabloATAClabel)

ClustPerf[12,1] <- diabloATACclustres[["ARI"]]
ClustPerf[12,2] <- diabloATACclustres[["NMI"]]

### seurat
p1 <- DimPlot(data.rna, reduction = "umap", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(data.atac, reduction = "umap.atac", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend() + ggtitle("ATAC")
p1 + p2

seuratRNAclust = ConsensusClusterPlus(
  t(data.rna@reductions[["umap"]]@cell.embeddings),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

seuratRNAlabel <- unlist(seuratRNAclust[[subgroup]]["consensusClass"])
seuratRNAclustres <- clustComp(gnd, seuratRNAlabel)

ClustPerf[13,1] <- seuratRNAclustres[["ARI"]]
ClustPerf[13,2] <- seuratRNAclustres[["NMI"]]

seuratATACclust = ConsensusClusterPlus(
  t(data.atac@reductions[["umap.atac"]]@cell.embeddings),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

seuratATAClabel <- unlist(seuratATACclust[[subgroup]]["consensusClass"])
seuratATACclustres <- clustComp(gnd, seuratATAClabel)

ClustPerf[14,1] <- seuratATACclustres[["ARI"]]
ClustPerf[14,2] <- seuratATACclustres[["NMI"]]

### PanDA
data <- list(rna = t(as.matrix(reducedRNA)), atac = t(as.matrix(reducedATAC)))
PanDAModel <- PanDA(data,gnd,numComponents,meu=0.02)
factors_PanDARNA <- PanDAModel[["PanDAComponents"]][["rnaComponents"]]
factors_PanDAATAC <- PanDAModel[["PanDAComponents"]][["atacComponents"]]

PanDARNA <- CreateDimReducObject(
  embeddings = PanDAModel[["PanDAComponents"]][["rnaComponents"]],
  loadings = PanDAModel[["projMatrices"]][["Wrna"]],
  stdev = numeric(),
  key = "PANDA_",
  assay = "RNA"
)

data.rna@reductions[["PanDA"]] <- PanDARNA
rownames(data.rna@reductions[["PanDA"]]@cell.embeddings) <- colnames(data.rna@assays[["RNA"]])

PanDAATAC <- CreateDimReducObject(
  embeddings = PanDAModel[["PanDAComponents"]][["atacComponents"]],
  loadings = PanDAModel[["projMatrices"]][["Watac"]],
  stdev = numeric(),
  key = "PANDA_",
  assay = "ATAC"
)

data.atac@reductions[["PanDA"]] <- PanDAATAC
rownames(data.atac@reductions[["PanDA"]]@cell.embeddings) <- colnames(data.atac@assays[["ATAC"]])
# momix.atac <- RunTSNE(momix.atac, dims = 1:100, reduction = "PanDA")

p1 <- DimPlot(data.rna, reduction = "PanDA", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(data.atac, reduction = "PanDA", group.by = "group", label = FALSE, pt.size = 4, cols = cols) + NoLegend() + ggtitle("ATAC")
p1 + p2


PanDARNAclust = ConsensusClusterPlus(
  t(factors_PanDARNA[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

PanDARNAlabel <- unlist(PanDARNAclust[[subgroup]]["consensusClass"])
PanDARNAclustres <- clustComp(gnd, PanDARNAlabel)

ClustPerf[15,1] <- PanDARNAclustres[["ARI"]]
ClustPerf[15,2] <- PanDARNAclustres[["NMI"]]

PanDAATACclust = ConsensusClusterPlus(
  t(factors_PanDAATAC[,c(1:2)]),
  maxK=5,
  reps=100,
  distance="euclidean",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

PanDAATAClabel <- unlist(PanDAATACclust[[subgroup]]["consensusClass"])
PanDAATACclustres <- clustComp(gnd, PanDAATAClabel)

ClustPerf[16,1] <- PanDAATACclustres[["ARI"]]
ClustPerf[16,2] <- PanDAATACclustres[["NMI"]]


### boxplots
par(mar = c(7, 4, 2, 2) + 0.2)
# colors <- c("#008000","#008080","#FF8C00","#20B2AA","#708090","#CCCC00","#006666",
#             "#CC0066","#808000","#00BFFF","#ABA300","#924900","#006666","#006666")

colors <- c("#008000","#708090","#CCCC00","#CC0066","#808000","#265DAB","#006666","#008080","#FF8C00","#FF8C00",
                     "#00BFFF","#00BFFF","#A0522D","#A0522D","#20B2AA","#20B2AA")
                     # ARI barplot
barplot(height=ClustPerf$ARI,col=colors,names.arg=rownames(ClustPerf),las=2,font.axis=2, ylab = substitute(paste(bold("ARI"))))

# NMI barplot
barplot(height=ClustPerf$NMI,col=colors,names.arg=rownames(ClustPerf),las=2,font.axis=2, ylab = substitute(paste(bold("NMI"))))

