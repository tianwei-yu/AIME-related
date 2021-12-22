setwd('C:/users/tianw/Onedrive/1-work/autoencoder metaanalysis/CCLE')
setwd('D:/Onedrive/1-work/autoencoder metaanalysis/CCLE')
library(AIME)
source('D:/Onedrive/1-work/autoencoder metaanalysis/fdrgamma.r')
source('C:/users/tianw/Onedrive/1-work/autoencoder metaanalysis/fdrgamma.r')


load("CCLE_miRNA_20181103.bin")
load("CCLE_RNAseq_genes_rpkm_20180929.bin")

mir<-log10(1+mir)
gene<-log10(1+gene)
gene<-as.matrix(gene)
mir<-as.matrix(mir)

mir<-mir[,which(colnames(mir) %in% colnames(gene))]
gene<-gene[,which(colnames(gene) %in% colnames(mir))]

mir<-mir[,order(colnames(mir))]
gene<-gene[,order(colnames(gene))]

sum(colnames(mir)==colnames(gene))
dim(mir)

mir<-t(mir)
gene<-t(gene)

cv.mir<-apply(mir,2,sd)/apply(mir,2,mean)
mir<-mir[,cv.mir>=0.1]

ze<-apply(gene==0,2,sum)
gene<-gene[,which(ze<=0.25*nrow(gene))]
cv.gene<-apply(gene,2,sd)/apply(gene,2,mean)
gene<-gene[,cv.gene>=0.5]

############

load("mir to gene aime layers 3 4 dropout 0.3 no confounder.bin")
m.to.g<-rec

#### important mirs

m.imp<-m.to.g[[1]]$imp
for(i in 2:length(m.to.g)) m.imp<-cbind(m.imp, m.to.g[[i]]$imp)

summary(as.dist(cor(m.imp)))

# gene specific importance

mg<-m.to.g[[1]]$pair.imp
for(i in 2:10) mg<-mg+m.to.g[[i]]$pair.imp

rownames(mg)<-colnames(gene)
colnames(mg)<-colnames(mir)

ccc<-cor(gene, mir)

mir.names<-read.table("mir names A-GEOD-16231_comments.txt",sep="\t", header=T)
for(i in 1:ncol(mir))
{
	sel<-which(mir.names[,4] == colnames(mir)[i])
	#print(c(i,sel))
	colnames(mir)[i]<-mir.names[sel,1]
}


########### my own FDR with gamma null

m.imp<-apply(m.imp,1,mean)
m.fdr<-fdrgamma(m.imp)

m.aime<-colnames(mir)[which(m.fdr<1e-3)]

load("CCA_standardized_weights.bin")
weights<-apply(abs(weights),1,sum)
#w.cut<-fdrgamma.cut(weights,use.percentile=use.percentile.1, fdr.cut=0.1)
w.cut<-quantile(weights, (700-35)/700)
m.cca<-colnames(mir)[which(weights>w.cut)]

load("MOFA_standardized_weights.bin")
weights<-apply(abs(weights[[1]]),1,sum)
#w.cut<-fdrgamma.cut(weights,use.percentile=use.percentile.1, fdr.cut=0.1)
w.cut<-quantile(weights, (700-35)/700)
m.mofa<-colnames(mir)[which(weights>w.cut)]

load("iCluster results.Rdata")
w<-results$W[1:700,]
weights<-apply(w^2,1,sum)
w.cut<-quantile(weights, (700-35)/700)
m.iCluster<-colnames(mir)[which(weights>w.cut)]

library(eulerr)
plot(venn(list(AIME=m.aime, CCA=m.cca, MOFA2=m.mofa, icluster=m.iCluster)))


### followed up by mirPath website
