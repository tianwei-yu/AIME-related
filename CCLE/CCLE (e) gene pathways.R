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
m.imp<-apply(m.to.g[[1]]$pair.imp,1,max)
for(i in 2:length(m.to.g)) m.imp<-m.imp+apply(m.to.g[[i]]$pair.imp,1,max)

m.fdr<-fdrgamma(m.imp)
m.aime<-colnames(gene)[which(m.fdr<0.05)]

###

load("CCA_normalized_results.bin")
weights.cca<-apply((ppp$loadings$Y)^2,1,sum)
w.cut<-quantile(weights.cca, (length(weights.cca)-286)/length(weights.cca))
m.cca<-colnames(gene)[weights.cca>w.cut]


######## (4) MOFA2

library(MOFA2)
load("MOFA2_model.bin")
weights.mofa <- get_weights(model, views = "all", factors = "all")
weights.mofa<-apply(abs(weights.mofa[[2]]),1,sum)
w.cut<-quantile(weights.mofa, (length(weights.mofa)-286)/length(weights.mofa))
m.mofa<-colnames(gene)[weights.mofa>w.cut]


#####

load("iCluster results.Rdata")
		weights.icluster<-results$W[701:nrow(results$W),]
		weights.icluster<-apply(weights.icluster^2,1,sum)
w.cut<-quantile(weights.icluster, (length(weights.icluster)-286)/length(weights.icluster))
m.iCluster<-colnames(gene)[weights.icluster>w.cut]

######

library(eulerr)
plot(venn(list(AIME=m.aime, CCA=m.cca, MOFA2=m.mofa, iCluster2=m.iCluster)))


###

gg2<-list(AIME=m.imp, CCA=weights.cca, MOFA2=weights.mofa, iCluster2=weights.icluster)
newnames<-names<-names(gg2[[2]])
for(i in 1:length(names)) 
{
	this.name<-unlist(mget(substr(names[i],1,15), org.Hs.egENSEMBL2EG,ifnotfound=NA))[1]
	if(is.na(this.name)) this.name<-paste("isna",i)
	newnames[i]<-this.name
}
for(i in 1:4)  names(gg2[[i]])<-newnames

pathways <- reactomePathways(names(gg2[[1]]))

pa<-new("list")
for(i in 1:4)
{
	f<-fgseaSimple(pathways, gg2[[i]],nperm=5000)
	pa[[i]]<-f
}

save(pa, file="four methods all genes fgsea.bin")

pa2<-pa[[1]][,c(1,7,3)]
for(i in 2:4) pa2<-cbind(pa2, pa[[i]][,3])
colnames(pa2)[3:6]<-names(gg2)
pa2[is.na(pa2)]<-1
l<-apply(pa2[,3:6],1,min)
pa2<-pa2[l<=0.05,]
pa2<-pa2[which(pa2[,2]>=100 & pa2[,2]<=500),]
pa2[,1]<-as.vector(pa2[,1])
rownames.pa2<-unlist(as.vector(pa2[,1]))
pa2<-as.matrix(pa2[,-1])
rownames(pa2)<-rownames.pa2

panet<-c(0,0)
for(i in 1:nrow(pa2))
{
	for(j in 2:5)
	{
		if(pa2[i,j]<=0.05) panet<-rbind(panet, c(rownames(pa2)[i], colnames(pa2)[j]))
	}
}
panet<-panet[-1,]
write.table(panet, "pathway net 4 methods.txt", sep="\t", quote=F, col.names=F, row.names=F)

##### to be followed up using cytoscape
