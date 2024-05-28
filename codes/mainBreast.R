# delete work space
rm(list = ls(all = TRUE))
graphics.off()

# setwd("C:/Projects/mixOmics/Data/SNF data/Breast")

library(corrplot)
library(kernlab)
library(IntNMF)
library(mclust)
library(aricode)
library(ConsensusClusterPlus)
library("survival")
library("survcomp")
library(survminer)
library("forestplot")
library(reticulate)
library(aricode)
library(omicade4)
library(IntNMF)
library(r.jive)
library(tensorBSS)
library(mixOmics) # import the mixOmics library
library(umap)
library(RGCCA)
library(Seurat)
library(MOFA2)
library(ggforce)
library(PANDA)

# source("PanDA2.R")
source("C:/Projects/mixOmics/Data/SNF data/Breast/tICA.R")
source("C:/Projects/mixOmics/Data/SNF data/Breast/iCluster2.R")


mrna <- read.table("C:/Users/MAminu/OneDrive - Inside MD Anderson/Projects/DOmics working draft/Final/Nature Communications/Github/Codes/PanDA/Data/TCGA Breast/gene_expression.txt", header=TRUE, sep = ",")
rownames(mrna) <- mrna[,1]
mrna <- scale(mrna[,-c(1)])

meth <- read.table("C:/Users/MAminu/OneDrive - Inside MD Anderson/Projects/DOmics working draft/Final/Nature Communications/Github/Codes/PanDA/Data/TCGA Breast/meth_expression.txt", header=TRUE, sep = ",")
rownames(meth) <- meth[,1]
meth <- scale(meth[,-c(1)])

mirna <- read.table("C:/Users/MAminu/OneDrive - Inside MD Anderson/Projects/DOmics working draft/Final/Nature Communications/Github/Codes/PanDA/Data/TCGA Breast/mirna_expression.txt", header=TRUE, sep = ",")
rownames(mirna) <- mirna[,1]
mirna <- scale(mirna[,-c(1)])

surv <- read.table("C:/Users/MAminu/OneDrive - Inside MD Anderson/Projects/DOmics working draft/Final/Nature Communications/Github/Codes/PanDA/Data/TCGA Breast/survival.txt", header=TRUE)
colnames(surv) <- c("time","event","label")

data <- list(mirna = t(as.matrix(mirna)), mrna = t(as.matrix(mrna)), methylation = t(as.matrix(meth)))

Y <- surv$label # use the subtype as the outcome variable
subtype <- factor(Y)

col <- c("#be0000","#00468BFF")
# Converting Factor to numeric
gnd <- as.numeric(as.factor(Y))
numComponents <- 20

# create a Seurat object and add the assays
breast <- CreateSeuratObject(counts = t(mrna), assay = "mRNA", meta.data = surv)
breast@assays[["miRNA"]] <- CreateSeuratObject(counts = t(mirna))
breast@assays[["methylation"]] <- CreateSeuratObject(counts = t(meth))

### compare and save performance of different methods
ClustPerf <- data.frame(matrix(ncol = 3, nrow = 18))
colnames(ClustPerf) <- c('ARI', 'NMI')
rownames(ClustPerf) <- c("iCluster","intNMF","JIVE","MCIA","MOFA","MEFISTO","tICA_mRNA","tICA_miRNA","tICA_meth",
                         "RGCCA_mRNA","RGCCA_miRNA","RGCCA_meth","DIABLO_mRNA","DIABLO_miRNA","DIABLO_meth",
                         "PANDA_mRNA","PANDA_miRNA","PANDA_meth")

CIndex <- data.frame(matrix(ncol = 1, nrow = 10))
dCIndex <- data.frame(matrix(ncol = numComponents, nrow = 10))
colnames(CIndex) <- c("C-index")
rownames(CIndex) <- c("iCluster","intNMF","jIVE","MCIA","MOFA","MEFISTO","tICA","RGCCA","DIABLO","PANDA")
pvalueCI <- data.frame(matrix(ncol = 1, nrow = 10))
rownames(dCIndex) <- c("iCluster","intNMF","jIVE","MCIA","MOFA","MEFISTO","tICA","RGCCA","DIABLO","PANDA")
rownames(pvalueCI) <- c("iCluster","intNMF","jIVE","MCIA","MOFA","MEFISTO","tICA","RGCCA","DIABLO","PANDA")
colnames(pvalueCI) <- c("pvalues")

omics <- list(mrna = as.matrix(t(mrna)), mirna = as.matrix(t(mirna)), methylation = as.matrix(t(meth)))
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
factorizations_icluster<-iCluster2(lapply(omics, function(x) t(x)), k=numComponents+1)
factors_icluster<-as.data.frame(as.matrix(t(factorizations_icluster$expZ)))
colnames(factors_icluster) <- paste(rep("iCluster",20),c(1:20), sep="_")

factors_icluster %>%
  ggplot(aes(x = `iCluster_1`,
             y = `iCluster_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("iCluster_1")+ylab("iCluster_2")+xlim(-0.4, 0.5)+ylim(-0.5, 0.6)+NoLegend()

iclusterclust = ConsensusClusterPlus(
  t(factors_icluster[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

subgroup <- 2
iclusterlabel <- unlist(iclusterclust[[subgroup]]["consensusClass"])

iclusterclustres <- clustComp(gnd, iclusterlabel)

ClustPerf[1,1] <- iclusterclustres[["ARI"]]
ClustPerf[1,2] <- iclusterclustres[["NMI"]]

iclusterForest <- cbind(factors_icluster, surv, iclusterlabel)

iclusterfit <- coxph(Surv(time/365,event)
                     ~ iclusterlabel, data = iclusterForest)

pvalueCI[1,1] <- summary(iclusterfit)$coefficients[5]

CIndex[1,1] <- iclusterfit[["concordance"]][["concordance"]]

for(k in 1:numComponents){
  if (k == 1){
    covariates <- colnames(iclusterForest)[1]
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    iclusterfit2 <- lapply( func, function(x){coxph(x, data = iclusterForest)})
    dCIndex[1,k] <- iclusterfit2[[covariates]][["concordance"]][["concordance"]]
  }else{
    covariates <- paste(colnames(iclusterForest)[1:k], collapse="+")
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    iclusterfit2 <- lapply( func, function(x){coxph(x, data = iclusterForest)})
    dCIndex[1,k] <- iclusterfit2[[covariates]][["concordance"]][["concordance"]]
  }
}

iclustergg <- coxph(Surv(time/365,event)
                    ~ iCluster_1+iCluster_2+iCluster_3+iCluster_4+iCluster_5+iCluster_6+iCluster_7+iCluster_8+iCluster_9+iCluster_10,
                    data = iclusterForest)

ggforest(iclustergg, main = "Breast cancer (iCluster components): Hazard Ratio", fontsize = 1)

cut <- surv_cutpoint(iclusterForest, variables='iCluster_1',minprop = 0.5)
iclusterForest$FactorCluster <- iclusterForest$iCluster_1 > cut$cutpoint$cutpoint
fitiCluster <- survfit(Surv(time/365, event) ~ FactorCluster, iclusterForest)
summary(coxph(Surv(time/365,event) ~ FactorCluster, data = iclusterForest))

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(size = 14, color = "black",hjust=0.5,face = "bold"),
      axis.text.x = element_text(size = 14, color = "black", face = "bold"),
      legend.text = element_text(size = 14, color = "black", face = "bold"),
      legend.title = element_text(size = 14, color = "black", face = "bold"),
      axis.text.y = element_text(size = 14, color = "black", face = "bold"),
      axis.title.x = element_text(size = 14, color = "black", face = "bold"),
      axis.title.y = element_text(size = 14, color = "black", face = "bold", angle = 90) , #angle=(90))
    )
}

ggsurvplot(fitiCluster, data = iclusterForest,title = "Breast cancer (iCluster first mRNA component)",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 16),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#be0000","#00468BFF"),

           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "",
           legend.labs = c("High risk  ",
                           "Low risk  "))

###intNMF
factorizations_intnmf<-nmf.mnnals(dat=lapply(omics_pos, function(x) t(x)), k=2)
factors_intNMF<-as.data.frame(as.matrix(factorizations_intnmf$W))
colnames(factors_intNMF) <- paste(rep("intNMF",2),c(1:2), sep="_")

factors_intNMF %>%
  ggplot(aes(x = `intNMF_1`,
             y = `intNMF_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("intNMF_1")+ylab("intNMF_2")+xlim(0, 0.15)+ylim(0, 0.2)+NoLegend()

intNMFclust = ConsensusClusterPlus(
  t(factors_intNMF[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

intNMFlabel <- unlist(intNMFclust[[subgroup]]["consensusClass"])
intNMFclustres <- clustComp(gnd, intNMFlabel)

ClustPerf[2,1] <- intNMFclustres[["ARI"]]
ClustPerf[2,2] <- intNMFclustres[["NMI"]]

intNMFForest <- cbind(factors_intNMF, surv, intNMFlabel)

intNMFfit <- coxph(Surv(time/365,event)
                   ~ intNMFlabel, data = intNMFForest)

pvalueCI[2,1] <- summary(intNMFfit)$coefficients[5]

CIndex[2,1] <- intNMFfit[["concordance"]][["concordance"]]

for(k in 1:2){
  if (k == 1){
    covariates <- colnames(intNMFForest)[1]
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    intNMFfit2 <- lapply( func, function(x){coxph(x, data = intNMFForest)})
    dCIndex[2,k] <- intNMFfit2[[covariates]][["concordance"]][["concordance"]]
  }else{
    covariates <- paste(colnames(intNMFForest)[1:k], collapse="+")
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    intNMFfit2 <- lapply( func, function(x){coxph(x, data = intNMFForest)})
    dCIndex[2,k] <- intNMFfit2[[covariates]][["concordance"]][["concordance"]]
  }
}

intNMFgg <- coxph(Surv(time/365,event)
                  ~ intNMF_1+intNMF_2,data = intNMFForest)

ggforest(intNMFgg, main = "Breast cancer (intNMF components): Hazard Ratio", fontsize = 1)

cut <- surv_cutpoint(intNMFForest, variables='intNMF_1', minprop = 0.5)
intNMFForest$FactorCluster <- intNMFForest$intNMF_1 > cut$cutpoint$cutpoint
fitintNMF <- survfit(Surv(time/365, event) ~ FactorCluster, intNMFForest)
summary(coxph(Surv(time/365,event) ~ FactorCluster, data = intNMFForest))

ggsurvplot(fitintNMF, data = intNMFForest,title = "Breast cancer (intNMF first mRNA component)",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 16),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#be0000","#00468BFF"),

           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "",
           legend.labs = c("High risk  ",
                           "Low risk  "))

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
factors_jive=as.data.frame(jV[,1:rankJV])
colnames(factors_jive) <- paste(rep("jive",20),c(1:20), sep="_")

factors_jive %>%
  ggplot(aes(x = `jive_1`,
             y = `jive_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("jive_1")+ylab("jive_2")+xlim(-0.0008, 0.0005)+ylim(-0.0005, 0.0005)+NoLegend()

jiveclust = ConsensusClusterPlus(
  t(factors_jive[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

jivelabel <- unlist(jiveclust[[subgroup]]["consensusClass"])
jiveclustres <- clustComp(gnd, jivelabel)

ClustPerf[3,1] <- jiveclustres[["ARI"]]
ClustPerf[3,2] <- jiveclustres[["NMI"]]

jiveForest <- cbind(factors_jive, surv, jivelabel)

jivefit <- coxph(Surv(time/365,event)
                 ~ jivelabel, data = jiveForest)

pvalueCI[3,1] <- summary(jivefit)$coefficients[5]

CIndex[3,1] <- jivefit[["concordance"]][["concordance"]]

for(k in 1:numComponents){
  if (k == 1){
    covariates <- colnames(jiveForest)[1]
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    jivefit2 <- lapply( func, function(x){coxph(x, data = jiveForest)})
    dCIndex[3,k] <- jivefit2[[covariates]][["concordance"]][["concordance"]]
  }else{
    covariates <- paste(colnames(jiveForest)[1:k], collapse="+")
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    jivefit2 <- lapply( func, function(x){coxph(x, data = jiveForest)})
    dCIndex[3,k] <- jivefit2[[covariates]][["concordance"]][["concordance"]]
  }
}

jivefitgg <- coxph(Surv(time/365,event)
                   ~ jive_1+jive_2+jive_3+jive_4+jive_5+jive_6+jive_7+jive_8+jive_9+jive_10, data = jiveForest)

# ggforest(jivefitgg,main = "breast (jive components): Hazard Ratio", fontsize = 0.8)

cut <- surv_cutpoint(jiveForest, variables='jive_1',minprop = 0.5)
jiveForest$FactorCluster <- jiveForest$jive_1 > cut$cutpoint$cutpoint
fitjive <- survfit(Surv(time/365, event) ~ FactorCluster, jiveForest)
summary(coxph(Surv(time/365,event) ~ FactorCluster, data = jiveForest))

ggsurvplot(fitjive, data = jiveForest,title = "Breast cancer (JIVE first mRNA component)",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 16),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#be0000","#00468BFF"),

           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "",
           legend.labs = c("High risk  ",
                           "Low risk  "))

### MCIA
factorizations_mcia<-mcia(omics_pos, cia.nf = numComponents)
factors_mcia<-as.data.frame(as.matrix(factorizations_mcia$mcoa$SynVar))
colnames(factors_mcia) <- paste(rep("MCIA",20),c(1:20), sep="_")

factors_mcia %>%
  ggplot(aes(x = `MCIA_1`,
             y = `MCIA_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("MCIA_1")+ylab("MCIA_2")+xlim(-3, 3)+ylim(-4, 3)+NoLegend()

mciaclust = ConsensusClusterPlus(
  t(factors_mcia[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

mcialabel <- unlist(mciaclust[[subgroup]]["consensusClass"])
mciaclustres <- clustComp(gnd, mcialabel)

ClustPerf[4,1] <- mciaclustres[["ARI"]]
ClustPerf[4,2] <- mciaclustres[["NMI"]]

mciaForest <- cbind(factors_mcia, surv, mcialabel)

mciafit <- coxph(Surv(time/365,event)
                 ~ mcialabel, data = mciaForest)

pvalueCI[4,1] <- summary(mciafit)$coefficients[5]

CIndex[4,1] <- mciafit[["concordance"]][["concordance"]]

for(k in 1:numComponents){
  if (k == 1){
    covariates <- colnames(mciaForest)[1]
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    mciafit2 <- lapply( func, function(x){coxph(x, data = mciaForest)})
    dCIndex[4,k] <- mciafit2[[covariates]][["concordance"]][["concordance"]]
  }else{
    covariates <- paste(colnames(mciaForest)[1:k], collapse="+")
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    mciafit2 <- lapply( func, function(x){coxph(x, data = mciaForest)})
    dCIndex[4,k] <- mciafit2[[covariates]][["concordance"]][["concordance"]]
  }
}

mciafitgg <- coxph(Surv(time/365,event)
                   ~ MCIA_1+MCIA_2+MCIA_3+MCIA_4+MCIA_5+MCIA_6+MCIA_7+MCIA_8+MCIA_9+MCIA_10, data = mciaForest)

ggforest(mciafitgg,main = "Breast cancer (MCIA mRNA components): Hazard Ratio", fontsize = 1)

cut <- surv_cutpoint(mciaForest, variables='MCIA_1',minprop = 0.5)
mciaForest$FactorCluster <- mciaForest$MCIA_1 > cut$cutpoint$cutpoint
fitmcia <- survfit(Surv(time/365, event) ~ FactorCluster, mciaForest)
summary(coxph(Surv(time/365,event) ~ FactorCluster, data = mciaForest))

ggsurvplot(fitmcia, data = mciaForest,title = "Breast cancer (MCIA first mRNA component)",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 16),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#be0000","#00468BFF"),

           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "",
           legend.labs = c("High risk  ",
                           "Low risk  "))

### MOFA
mofa <- create_mofa(list(
  "mRNA" = as.matrix(t(mrna)),
  "mirna" = as.matrix(t(mirna),
                      "methylation" = as.matrix(t(meth)))
))

# Model options: let's use only 4 factors, should be enough to distinguish the four cell lines.
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

factors_mofa <- as.data.frame(mofa@expectations[["Z"]][["group1"]])
colnames(factors_mofa) <- paste(rep("MOFA",20),c(1:20), sep="_")

factors_mofa %>%
  ggplot(aes(x = `MOFA_1`,
             y = `MOFA_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("MOFA_1")+ylab("MOFA_2")+xlim(-2, 4)+ylim(-2.5, 2.5)+NoLegend()

mofaclust = ConsensusClusterPlus(
  t(factors_mofa[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

mofalabel <- unlist(mofaclust[[subgroup]]["consensusClass"])
mofaclustres <- clustComp(gnd, mofalabel)

ClustPerf[5,1] <- mofaclustres[["ARI"]]
ClustPerf[5,2] <- mofaclustres[["NMI"]]

mofaForest <- cbind(factors_mofa, surv, mofalabel)

mofafit <- coxph(Surv(time/365,event)
                 ~ mofalabel, data = mofaForest)

pvalueCI[5,1] <- summary(mofafit)$coefficients[5]

CIndex[5,1] <- mofafit[["concordance"]][["concordance"]]

for(k in 1:numComponents){
  if (k == 1){
    covariates <- colnames(mofaForest)[1]
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    mofafit2 <- lapply( func, function(x){coxph(x, data = mofaForest)})
    dCIndex[5,k] <- mofafit2[[covariates]][["concordance"]][["concordance"]]
  }else{
    covariates <- paste(colnames(mofaForest)[1:k], collapse="+")
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    mofafit2 <- lapply( func, function(x){coxph(x, data = mofaForest)})
    dCIndex[5,k] <- mofafit2[[covariates]][["concordance"]][["concordance"]]
  }
}

mofafitgg <- coxph(Surv(time/365,event)
                   ~ MOFA_1+MOFA_2+MOFA_3+MOFA_4+MOFA_5+MOFA_6+MOFA_7+MOFA_8+MOFA_9+MOFA_10, data = mofaForest)
summary(mofafitgg)
ggforest(mofafitgg,main = "Breast cancer (MOFA mRNA components): Hazard Ratio", fontsize = 1)

cut <- surv_cutpoint(mofaForest, variables='MOFA_1',minprop = 0.5)
mofaForest$FactorCluster <- mofaForest$MOFA_1 > cut$cutpoint$cutpoint
fitmofa <- survfit(Surv(time/365, event) ~ FactorCluster, mofaForest)
summary(coxph(Surv(time/365,event) ~ FactorCluster, data = mofaForest))

ggsurvplot(fitmofa, data = mofaForest,title = "Breast cancer (MOFA first mRNA component)",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 16),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#00468BFF","#be0000"),

           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "",
           legend.labs = c("Low risk  ",
                           "High risk  "))

### MEFISTO
mefisto <- create_mofa(list(
  "mRNA" = as.matrix(t(mrna)),
  "mirna" = as.matrix(t(mirna),
                      "methylation" = as.matrix(t(meth)))
))

Princomp <- prcomp(mrna,rank. = 50,center = TRUE,scale. = TRUE)
umapcomp <- umap(Princomp[["x"]])
mefisto <- set_covariates(mefisto, t(umapcomp[["layout"]]))

mefisto_opts <- get_default_model_options(mefisto)
mefisto_opts$num_factors <- numComponents

mefisto <- prepare_mofa(
  object = mefisto,
  model_options = mefisto_opts,
)

mefisto <- run_mofa(mefisto)
factors_mefisto <- as.data.frame(mefisto@expectations[["Z"]][["group1"]])
colnames(factors_mefisto) <- paste(rep("MEFISTO",20),c(1:20), sep="_")

factors_mefisto %>%
  ggplot(aes(x = `MEFISTO_1`,
             y = `MEFISTO_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("MEFISTO_1")+ylab("MEFISTO_2")+xlim(-2.5, 4)+ylim(-2.5, 2.5)+NoLegend()


mefistoclust = ConsensusClusterPlus(
  t(factors_mefisto[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

mefistolabel <- unlist(mefistoclust[[subgroup]]["consensusClass"])
mefistoclustres <- clustComp(gnd, mefistolabel)

ClustPerf[6,1] <- mefistoclustres[["ARI"]]
ClustPerf[6,2] <- mefistoclustres[["NMI"]]

mefistoForest <- cbind(factors_mefisto, surv, mefistolabel)

mefistofit <- coxph(Surv(time/365,event)
                    ~ mefistolabel, data = mefistoForest)

pvalueCI[6,1] <- summary(mefistofit)$coefficients[5]

CIndex[6,1] <- mefistofit[["concordance"]][["concordance"]]

for(k in 1:numComponents){
  if (k == 1){
    covariates <- colnames(mefistoForest)[1]
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    mefistofit2 <- lapply( func, function(x){coxph(x, data = mefistoForest)})
    dCIndex[6,k] <- mefistofit2[[covariates]][["concordance"]][["concordance"]]
  }else{
    covariates <- paste(colnames(mefistoForest)[1:k], collapse="+")
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    mefistofit2 <- lapply( func, function(x){coxph(x, data = mefistoForest)})
    dCIndex[6,k] <- mefistofit2[[covariates]][["concordance"]][["concordance"]]
  }
}

mefistofitgg <- coxph(Surv(time/365,event)
                      ~ MEFISTO_1+MEFISTO_2+MEFISTO_3+MEFISTO_4+MEFISTO_5+MEFISTO_6+MEFISTO_7+MEFISTO_8+MEFISTO_9+MEFISTO_10,
                      data = mefistoForest)

summary(mefistofitgg)
ggforest(mefistofitgg,main = "Breast cancer (MEFISTO  mRNA components): Hazard Ratio", fontsize = 1)

cut <- surv_cutpoint(mefistoForest, variables='MEFISTO_1',minprop = 0.5)
mefistoForest$FactorCluster <- mefistoForest$MEFISTO_1 > cut$cutpoint$cutpoint
fitmefisto <- survfit(Surv(time/365, event) ~ FactorCluster, mefistoForest)
summary(coxph(Surv(time/365,event) ~ FactorCluster, data = mefistoForest))

ggsurvplot(fitmefisto, data = mefistoForest,title = "Breast cancer (MEFISTO first mRNA component)",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 16),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#00468BFF","#be0000"),

           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "",
           legend.labs = c("Low risk  ",
                           "High risk  "))

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
factors_ticamRNA <- as.data.frame(as.matrix(t(tICA[["projS"]][[1]])))
colnames(factors_ticamRNA) <- paste(rep("tICAmRNA",20),c(1:20), sep="_")

factors_ticamiRNA <- as.data.frame(as.matrix(t(tICA[["projS"]][[2]])))
colnames(factors_ticamiRNA) <- paste(rep("tICAmiRNA",20),c(1:20), sep="_")

factors_ticameth <- as.data.frame(as.matrix(t(tICA[["projS"]][[3]])))
colnames(factors_ticameth) <- paste(rep("tICAmeth",20),c(1:20), sep="_")

p1 <- factors_ticamRNA %>%
  ggplot(aes(x = `tICAmRNA_1`,
             y = `tICAmRNA_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("tICAmRNA_1")+ylab("tICAmRNA_2")+xlim(-20, 20)+ylim(-20, 20)+NoLegend()

p2 <- factors_ticamiRNA %>%
  ggplot(aes(x = `tICAmiRNA_1`,
             y = `tICAmiRNA_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("tICAmiRNA_1")+ylab("tICAmiRNA_2")+xlim(-4, 4)+ylim(-4, 4)+NoLegend()

p3 <- factors_ticameth %>%
  ggplot(aes(x = `tICAmeth_1`,
             y = `tICAmeth_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("tICAmeth_1")+ylab("tICAmeth_2")+xlim(-4, 4)+ylim(-4, 4)+NoLegend()

p1 + p2 + p3

ticamRNAclust = ConsensusClusterPlus(
  t(factors_ticamRNA[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

ticamRNAlabel <- unlist(ticamRNAclust[[subgroup]]["consensusClass"])
ticamRNAclustres <- clustComp(gnd, ticamRNAlabel)

ClustPerf[7,1] <- ticamRNAclustres[["ARI"]]
ClustPerf[7,2] <- ticamRNAclustres[["NMI"]]

ticamiRNAclust = ConsensusClusterPlus(
  t(factors_ticamiRNA[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

ticamiRNAlabel <- unlist(ticamiRNAclust[[subgroup]]["consensusClass"])
ticamiRNAclustres <- clustComp(gnd, ticamiRNAlabel)

ClustPerf[8,1] <- ticamiRNAclustres[["ARI"]]
ClustPerf[8,2] <- ticamiRNAclustres[["NMI"]]

ticamethclust = ConsensusClusterPlus(
  t(factors_ticameth[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

ticamethlabel <- unlist(ticamethclust[[subgroup]]["consensusClass"])
ticamethclustres <- clustComp(gnd, ticamethlabel)

ClustPerf[9,1] <- ticamethclustres[["ARI"]]
ClustPerf[9,2] <- ticamethclustres[["NMI"]]

ticaForest <- cbind(factors_ticamRNA, surv, ticamRNAlabel)

ticafit <- coxph(Surv(time/365,event)
                 ~ ticamRNAlabel, data = ticaForest)

pvalueCI[7,1] <- summary(ticafit)$coefficients[5]

CIndex[7,1] <- ticafit[["concordance"]][["concordance"]]

for(k in 1:numComponents){
  if (k == 1){
    covariates <- colnames(ticaForest)[1]
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    ticafit2 <- lapply( func, function(x){coxph(x, data = ticaForest)})
    dCIndex[7,k] <- ticafit2[[covariates]][["concordance"]][["concordance"]]
  }else{
    covariates <- paste(colnames(ticaForest)[1:k], collapse="+")
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    ticafit2 <- lapply( func, function(x){coxph(x, data = ticaForest)})
    dCIndex[7,k] <- ticafit2[[covariates]][["concordance"]][["concordance"]]
  }
}

ticafitgg <- coxph(Surv(time/365,event)
                   ~ tICAmRNA_1+tICAmRNA_2+tICAmRNA_3+tICAmRNA_4+tICAmRNA_5+tICAmRNA_6+tICAmRNA_7+tICAmRNA_8+tICAmRNA_9+tICAmRNA_10,
                   data = ticaForest)

summary(ticafitgg)
ggforest(ticafitgg,main = "Breast cancer (tICA mRNA components): Hazard Ratio", fontsize = 1)

cut <- surv_cutpoint(ticaForest, variables='tICAmRNA_1',minprop = 0.5)
ticaForest$FactorCluster <- ticaForest$tICAmRNA_1 > cut$cutpoint$cutpoint
fittica <- survfit(Surv(time/365, event) ~ FactorCluster, ticaForest)
summary(coxph(Surv(time/365,event) ~ FactorCluster, data = ticaForest))

ggsurvplot(fittica, data = ticaForest,title = "Breast cancer (tICA first mRNA components)",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 16),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#00468BFF","#be0000"),

           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "",
           legend.labs = c("Low risk  ",
                           "High risk  "))

### RGCCA
factorizations_RGCCA<-rgcca(lapply(omics, function(x) t(x)), ncomp = rep(numComponents, length(omics)), scheme = "centroid", scale = TRUE, init = "svd",bias = TRUE, tol = 1e-08, verbose = F)
factors_rgccamRNA <- as.data.frame(factorizations_RGCCA[["Y"]][[1]])
colnames(factors_rgccamRNA) <- paste(rep("rgccamRNA",20),c(1:20), sep="_")

factors_rgccamiRNA <- as.data.frame(factorizations_RGCCA[["Y"]][[2]])
colnames(factors_rgccamiRNA) <- paste(rep("rgccamiRNA",20),c(1:20), sep="_")

factors_rgccameth <- as.data.frame(factorizations_RGCCA[["Y"]][[3]])
colnames(factors_rgccameth) <- paste(rep("rgccameth",20),c(1:20), sep="_")

p1 <- factors_rgccamRNA %>%
  ggplot(aes(x = `rgccamRNA_1`,
             y = `rgccamRNA_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("rgccamRNA_1")+ylab("rgccamRNA_2")+xlim(-1.5, 1.5)+ylim(-1, 1)+NoLegend()

p2 <- factors_rgccamiRNA %>%
  ggplot(aes(x = `rgccamiRNA_1`,
             y = `rgccamiRNA_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("rgccamiRNA_1")+ylab("rgccamiRNA_2")+xlim(-1.5, 1.5)+ylim(-1, 1)+NoLegend()

p3 <- factors_rgccameth %>%
  ggplot(aes(x = `rgccameth_1`,
             y = `rgccameth_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("rgccameth_1")+ylab("rgccameth_2")+xlim(-2, 2)+ylim(-1, 1.5)+NoLegend()

p1 + p2 + p3

rgccamRNAclust = ConsensusClusterPlus(
  t(factors_rgccamRNA[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

rgccamRNAlabel <- unlist(rgccamRNAclust[[subgroup]]["consensusClass"])
rgccamRNAclustres <- clustComp(gnd, rgccamRNAlabel)

ClustPerf[10,1] <- rgccamRNAclustres[["ARI"]]
ClustPerf[10,2] <- rgccamRNAclustres[["NMI"]]

rgccamiRNAclust = ConsensusClusterPlus(
  t(factors_rgccamiRNA[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

rgccamiRNAlabel <- unlist(rgccamiRNAclust[[subgroup]]["consensusClass"])
rgccamiRNAclustres <- clustComp(gnd, rgccamiRNAlabel)

ClustPerf[11,1] <- rgccamiRNAclustres[["ARI"]]
ClustPerf[11,2] <- rgccamiRNAclustres[["NMI"]]

rgccamethclust = ConsensusClusterPlus(
  t(factors_rgccameth[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

rgccamethlabel <- unlist(rgccamethclust[[subgroup]]["consensusClass"])
rgccamethclustres <- clustComp(gnd, rgccamethlabel)

ClustPerf[12,1] <- rgccamethclustres[["ARI"]]
ClustPerf[12,2] <- rgccamethclustres[["NMI"]]

rgccaForest <- cbind(factors_rgccamRNA, surv, rgccamRNAlabel)

rgccafit <- coxph(Surv(time/365,event)
                  ~ rgccamRNAlabel, data = rgccaForest)

pvalueCI[8,1] <- summary(rgccafit)$coefficients[5]

CIndex[8,1] <- rgccafit[["concordance"]][["concordance"]]

for(k in 1:numComponents){
  if (k == 1){
    covariates <- colnames(rgccaForest)[1]
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    rgccafit2 <- lapply( func, function(x){coxph(x, data = rgccaForest)})
    dCIndex[8,k] <- rgccafit2[[covariates]][["concordance"]][["concordance"]]
  }else{
    covariates <- paste(colnames(rgccaForest)[1:k], collapse="+")
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    rgccafit2 <- lapply( func, function(x){coxph(x, data = rgccaForest)})
    dCIndex[8,k] <- rgccafit2[[covariates]][["concordance"]][["concordance"]]
  }
}

rgccafitgg <- coxph(Surv(time/365,event)
                    ~ rgccamRNA_1+rgccamRNA_2+rgccamRNA_3+rgccamRNA_4+rgccamRNA_5+rgccamRNA_6+rgccamRNA_7+rgccamRNA_8+rgccamRNA_9+rgccamRNA_10,
                    data = rgccaForest)

summary(rgccafitgg)
ggforest(rgccafitgg,main = "Breast cancer (RGCCA mRNA components): Hazard Ratio", fontsize = 1)

cut <- surv_cutpoint(rgccaForest, variables='rgccamRNA_1',minprop = 0.5)
rgccaForest$FactorCluster <- rgccaForest$rgccamRNA_1 > cut$cutpoint$cutpoint
fitrgcca<- survfit(Surv(time/365, event) ~ FactorCluster, rgccaForest)
summary(coxph(Surv(time/365,event) ~ FactorCluster, data = rgccaForest))

ggsurvplot(fitrgcca, data = rgccaForest,title = "Breast cancer (RGCCA first mRNA component)",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 16),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#00468BFF","#be0000"),

           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "",
           legend.labs = c("Low risk  ",
                           "High risk  "))

### DIABLO
data <- list(mrna = as.matrix(mrna), mirna = as.matrix(mirna), meth = as.matrix(meth))

diablo <- block.plsda(data, gnd, ncomp = numComponents, design = "null") # run the DIABLO method
factors_diablomRNA <- as.data.frame(as.matrix(diablo[["variates"]][["mrna"]]))
colnames(factors_diablomRNA) <- paste(rep("diablomRNA",20),c(1:20), sep="_")

factors_diablomiRNA <- as.data.frame(as.matrix(diablo[["variates"]][["mirna"]]))
colnames(factors_diablomiRNA) <- paste(rep("diablomiRNA",20),c(1:20), sep="_")

factors_diablometh <- as.data.frame(as.matrix(diablo[["variates"]][["meth"]]))
colnames(factors_diablometh) <- paste(rep("diablometh",20),c(1:20), sep="_")

p1 <- factors_diablomRNA %>%
  ggplot(aes(x = `diablomRNA_1`,
             y = `diablomRNA_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("diablomRNA_1")+ylab("diablomRNA_2")+xlim(-15, 15)+ylim(-15, 15)+NoLegend()

p2 <- factors_diablomiRNA %>%
  ggplot(aes(x = `diablomiRNA_1`,
             y = `diablomiRNA_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("diablomiRNA_1")+ylab("diablomiRNA_2")+xlim(-15, 15)+ylim(-15, 15)+NoLegend()

p3 <- factors_diablometh %>%
  ggplot(aes(x = `diablometh_1`,
             y = `diablometh_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("diablometh_1")+ylab("diablometh_2")+xlim(-30, 30)+ylim(-10, 10)+NoLegend()

p1 + p2 + p3

diablomRNAclust = ConsensusClusterPlus(
  t(factors_diablomRNA[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

diablomRNAlabel <- unlist(diablomRNAclust[[subgroup]]["consensusClass"])
diablomRNAclustres <- clustComp(gnd, diablomRNAlabel)

ClustPerf[13,1] <- diablomRNAclustres[["ARI"]]
ClustPerf[13,2] <- diablomRNAclustres[["NMI"]]


diablomiRNAclust = ConsensusClusterPlus(
  t(factors_diablomiRNA[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

diablomiRNAlabel <- unlist(diablomiRNAclust[[subgroup]]["consensusClass"])
diablomiRNAclustres <- clustComp(gnd, diablomiRNAlabel)

ClustPerf[14,1] <- diablomiRNAclustres[["ARI"]]
ClustPerf[14,2] <- diablomiRNAclustres[["NMI"]]

diablomethclust = ConsensusClusterPlus(
  t(factors_diablometh[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

diablomethlabel <- unlist(diablomethclust[[subgroup]]["consensusClass"])
diablomethclustres <- clustComp(gnd, diablomethlabel)

ClustPerf[15,1] <- diablomethclustres[["ARI"]]
ClustPerf[15,2] <- diablomethclustres[["NMI"]]

diabloForest <- cbind(factors_diablomRNA, surv, diablomRNAlabel)

diablofit <- coxph(Surv(time/365,event)
                   ~ diablomethlabel, data = diabloForest)

pvalueCI[9,1] <- summary(diablofit)$coefficients[5]

CIndex[9,1] <- diablofit[["concordance"]][["concordance"]]

for(k in 1:numComponents){
  if (k == 1){
    covariates <- colnames(diabloForest)[1]
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    diablofit2 <- lapply( func, function(x){coxph(x, data = diabloForest)})
    dCIndex[9,k] <- diablofit2[[covariates]][["concordance"]][["concordance"]]
  }else{
    covariates <- paste(colnames(diabloForest)[1:k], collapse="+")
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    diablofit2 <- lapply( func, function(x){coxph(x, data = diabloForest)})
    dCIndex[9,k] <- diablofit2[[covariates]][["concordance"]][["concordance"]]
  }
}

diablofitgg <- coxph(Surv(time/365,event)
                     ~ diablomRNA_1+diablomRNA_2+diablomRNA_3+diablomRNA_4+diablomRNA_5+diablomRNA_6+diablomRNA_7+
                       diablomRNA_8+diablomRNA_9+diablomRNA_10,data = diabloForest)

summary(diablofitgg)
ggforest(diablofitgg,main = "Breast cancer (DIABLO mRNA components): Hazard Ratio", fontsize = 1)

cut <- surv_cutpoint(diabloForest, variables='diablomRNA_1',minprop = 0.5)
diabloForest$FactorCluster <- diabloForest$diablomRNA_1 > cut$cutpoint$cutpoint
fitdiablo<- survfit(Surv(time/365, event) ~ FactorCluster, diabloForest)
summary(coxph(Surv(time/365,event) ~ FactorCluster, data = diabloForest))

ggsurvplot(fitdiablo, data = diabloForest,title = "Breast cancer (DIABLO first mRNA component)",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 16),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#00468BFF","#be0000"),

           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "",
           legend.labs = c("Low risk  ",
                           "High risk  "))

### PanDA
data <- list(mrna = t(as.matrix(mrna)), mirna = t(as.matrix(mirna)), meth = t(as.matrix(meth)))
PanDAModel <- PanDA(data,gnd,numComponents,0.02)
factors_PanDAmRNA <- as.data.frame(PanDAModel[["PanDAComponents"]][["mrnaComponents"]])
colnames(factors_PanDAmRNA) <- paste(rep("PanDAmRNA",20),c(1:20), sep="_")

factors_PanDAmiRNA <- as.data.frame(PanDAModel[["PanDAComponents"]][["mirnaComponents"]])
colnames(factors_PanDAmiRNA) <- paste(rep("PanDAmiRNA",20),c(1:20), sep="_")

factors_PanDAmeth <- as.data.frame(PanDAModel[["PanDAComponents"]][["methComponents"]])
colnames(factors_PanDAmeth) <- paste(rep("PanDAmeth",20),c(1:20), sep="_")

p1 <- factors_PanDAmRNA %>%
  ggplot(aes(x = `PanDAmRNA_1`,
             y = `PanDAmRNA_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("PandamRNA_1")+ylab("PandamRNA_2")+xlim(-1, 1)+ylim(-0.04, 0.04)+NoLegend()

p2 <- factors_PanDAmiRNA %>%
  ggplot(aes(x = `PanDAmiRNA_1`,
             y = `PanDAmiRNA_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("PandamiRNA_1")+ylab("PandamiRNA_2")+xlim(-1, 1)+ylim(-0.04, 0.04)+NoLegend()

p3 <- factors_PanDAmeth %>%
  ggplot(aes(x = `PanDAmeth_1`,
             y = `PanDAmeth_2`))+
  scale_color_manual(values=col)+
  geom_point(aes(color=subtype),size = 4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))+xlab("Pandameth_1")+ylab("Pandameth_2")+xlim(-1, 1)+ylim(-0.04, 0.04)

p1 + p2 + p3

PanDAmRNAclust = ConsensusClusterPlus(
  t(factors_PanDAmRNA[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

PanDAmRNAlabel <- unlist(PanDAmRNAclust[[subgroup]]["consensusClass"])
PanDAmRNAclustres <- clustComp(gnd, PanDAmRNAlabel)

ClustPerf[16,1] <- PanDAmRNAclustres[["ARI"]]
ClustPerf[16,2] <- PanDAmRNAclustres[["NMI"]]

PanDAmiRNAclust = ConsensusClusterPlus(
  t(factors_PanDAmiRNA[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

PanDAmiRNAlabel <- unlist(PanDAmiRNAclust[[subgroup]]["consensusClass"])
PanDAmiRNAclustres <- clustComp(gnd, PanDAmiRNAlabel)

ClustPerf[17,1] <- PanDAmiRNAclustres[["ARI"]]
ClustPerf[17,2] <- PanDAmiRNAclustres[["NMI"]]

PanDAmethclust = ConsensusClusterPlus(
  t(factors_PanDAmeth[,c(1:2)]),
  maxK=3,
  reps=100,
  distance="canberra",
  tmyPal=c("white","#135078"),
  clusterAlg="hc")

PanDAmethlabel <- unlist(PanDAmethclust[[subgroup]]["consensusClass"])
PanDAmethclustres <- clustComp(gnd, PanDAmethlabel)

ClustPerf[18,1] <- PanDAmethclustres[["ARI"]]
ClustPerf[18,2] <- PanDAmethclustres[["NMI"]]

PanDAForest <- cbind(factors_PanDAmRNA, surv, PanDAmRNAlabel)

PanDAfit <- coxph(Surv(time/365,event)
                   ~ PanDAmRNAlabel, data = PanDAForest)
PanDAfit[["concordance"]][["concordance"]]

pvalueCI[10,1] <- summary(PanDAfit)$coefficients[5]

CIndex[10,1] <- PanDAfit[["concordance"]][["concordance"]]

for(k in 1:numComponents){
  if (k == 1){
    covariates <- colnames(PanDAForest)[1]
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    PanDAfit2 <- lapply( func, function(x){coxph(x, data = PanDAForest)})
    dCIndex[10,k] <- PanDAfit2[[covariates]][["concordance"]][["concordance"]]
  }else{
    covariates <- paste(colnames(PanDAForest)[1:k], collapse="+")
    func <- sapply(covariates,
                   function(x) as.formula(paste('Surv(time/365,event)~', x)))
    PanDAfit2 <- lapply( func, function(x){coxph(x, data = PanDAForest)})
    dCIndex[10,k] <- PanDAfit2[[covariates]][["concordance"]][["concordance"]]
  }
}

PanDAfitgg <- coxph(Surv(time/365,event)
                     ~ PanDAmRNA_1+PanDAmRNA_2+PanDAmRNA_3+PanDAmRNA_4+PanDAmRNA_5+PanDAmRNA_6+PanDAmRNA_7+PanDAmRNA_8+
                       PanDAmRNA_9+PanDAmRNA_10,data = PanDAForest)

summary(PanDAfitgg)
ggforest(PanDAfitgg,main = "Breast cancer (PanDA mRNA components): Hazard Ratio", fontsize = 1)

cut <- surv_cutpoint(PanDAForest, variables='PanDAmRNA_1',minprop = 0.5)
PanDAForest$FactorCluster <- PanDAForest$PanDAmRNA_1 > cut$cutpoint$cutpoint
fitPanDA <- survfit(Surv(time/365, event) ~ FactorCluster, PanDAForest)
summary(coxph(Surv(time/365,event) ~ FactorCluster, data = PanDAForest))

ggsurvplot(fitPanDA, data = PanDAForest,title = "Breast cancer (PanDA first mrna component)",ggtheme=custom_theme(),
           conf.int = FALSE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           xlab = "Time (Years)",
           ylab = "Overall Survival (%)",
           xlim = c(0, 16),
           risk.table.fontsize =5,
           size = 2,
           linetype = "solid",
           palette = c("#be0000","#00468BFF"),

           risk.table.col = "strata",
           #legend = "bottom",
           legend.title = "",
           legend.labs = c("High risk  ",
                           "Low risk  "))

### boxplots
par(mar = c(7, 4, 2, 2) + 0.2)
# colors <- c("#008000","#008080","#FF8C00","#20B2AA","#708090","#CCCC00","#006666",
#             "#CC0066","#808000","#00BFFF","#ABA300","#924900","#006666","#006666")

colors <- c("#008000","#708090","#CCCC00","#CC0066","#808000","#265DAB","#006666","#006666","#006666",
            "#FF8C00","#FF8C00","#FF8C00","#00BFFF","#00BFFF","#00BFFF","#20B2AA","#20B2AA","#20B2AA")
colors2 <- c("#008000","#708090","#CCCC00","#CC0066","#808000","#265DAB","#006666","#FF8C00","#00BFFF","#20B2AA")
# ARI barplot
barplot(height=ClustPerf$ARI,col=colors,names.arg=rownames(ClustPerf), ylab = substitute(paste(bold("ARI"))),las=2,font.axis=2)

# NMI barplot
barplot(height=ClustPerf$NMI,col=colors,names.arg=rownames(ClustPerf), ylab = substitute(paste(bold("NMI"))),las=2,font.axis=2)

# C-Index barplot
barplot(height=CIndex$`C-index`,col=colors2,names.arg=rownames(CIndex), ylab = substitute(paste(bold("C-Index"))),las=2,font.axis=2, ylim=c(0, 0.8))

nvar <- 1:20
dCIndex <- as.data.frame(t(dCIndex))
# Create a first line
plot(nvar, dCIndex$iCluster, type = "b", frame = FALSE, pch = 1, lwd=3,
     xlim=c(0,  20), col = "#008000", xlab = substitute(paste(bold("# of components"))), ylab = substitute(paste(bold("C-Index"))), ylim=c(0.4,1),
     las=2, font=2, cex.lab=1.4, cex.axis=1.4, las=1)
# Add a second line
lines(nvar, dCIndex$intNMF, pch = 2, col = "#708090", type = "b", lty = 1, lwd=3)
# Add a second line
lines(nvar, dCIndex$jIVE, pch = 3, col = "#CCCC00", type = "b", lty = 1, lwd=3)
# Add a second line
lines(nvar, dCIndex$MCIA, pch = 4, col = "#CC0066", type = "b", lty = 1, lwd=3)
# Add a second line
lines(nvar, dCIndex$MOFA, pch = 5, col = "#808000", type = "b", lty = 1, lwd=3)
# Add a second line
lines(nvar, dCIndex$MEFISTO, pch = 6, col = "#265DAB", type = "b", lty = 1, lwd=3)
# Add a second line
lines(nvar, dCIndex$tICA, pch = 7, col = "#006666", type = "b", lty = 1, lwd=3)
# Add a second line
lines(nvar, dCIndex$RGCCA, pch = 8, col = "#FF8C00", type = "b", lty = 1, lwd=3)
# Add a second line
lines(nvar, dCIndex$DIABLO, pch = 9, col = "#00BFFF", type = "b", lty = 1, lwd=3)
# Add a second line
lines(nvar, dCIndex$PANDA, pch = 10, col = "#20B2AA", type = "b", lty = 1, lwd=3)


# Add a legend to the plot
legend("bottomright", legend=c("iCluster","intNMF","jIVE","MCIA","MOFA","MEFISTO","tICA","RGCCA","DIABLO","PANDA"),
       col=colors2, text.font = 2, pch = 1:10, lty = 1, cex=0.8)

# save.image("C:/Projects/mixOmics/Data/SNF data/Breast/Figures 3/BreastWorkspace.RData")

# cfea <- as.data.frame(cbind(rownames(factors_icluster),factors_icluster[,1],factors_intNMF[,1],factors_jive[,1],factors_mcia[,1],factors_mofa[,1],factors_mefisto[,1],
#                             factors_ticamRNA[,1],factors_rgccamRNA[,1],factors_diablomRNA[,1],factors_PanDAmRNA[,1]))
#
# write.xlsx(cfea, file = "FirstComponents_AllMoldels.xlsx", col.names = TRUE)

