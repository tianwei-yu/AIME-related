setwd("D:\\OneDrive\\1-work\\autoencoder metaanalysis\\CCLE")
setwd("C:\\Users\\tianw\\OneDrive\\1-work\\autoencoder metaanalysis\\CCLE")

library(AIME)
source("dist_group.R")

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

cells<-rownames(gene)
cell.type<-cells
for(i in 1:length(cells))
{
	this<-strsplit(cells[i], "_")[[1]][-1]
	this<-paste(this, collapse="_")
	cell.type[i]<-this
}

ttt<-table(cell.type)
sel<-which(cell.type %in% names(ttt)[ttt<10])
cell.type[sel]<-NA

anno<-read.table("Cell_lines_annotations_20181226.txt",header=T,sep="\t")
anno<-anno[which(anno[,1] %in% cells),]
anno<-anno[order(anno[,1]),]

##### limit to cells with annotations
s.anno<-which(cells %in% anno[,1])
cells<-cells[s.anno]
cell.type<-cell.type[s.anno]
sum(cells == anno[,1])
mir<-mir[s.anno,]
gene<-gene[s.anno,]

mir2<-mir
gene2<-gene

for(i in 1:ncol(mir2)) mir2[,i]<-(mir2[,i]-mean(mir2[,i]))/sd(mir2[,i])
for(i in 1:ncol(gene2)) gene2[,i]<-(gene2[,i]-mean(gene2[,i]))/sd(gene2[,i])

####

all.colors<-c("black","green","blue","cyan", "yellow","orange","red","grey50", "wheat","purple","darkblue","dodgerblue4","darkred","darkorange","darkcyan","magenta","firebrick","khaki4","cornsilk3","darkgoldenrod4")

load("mofa_embeded.Rdata")
load("canon_embeded.Rdata")
load("mir to gene aime layers 3 4 dropout 0.3 no confounder no pair.bin")
load("iCluster_embeded.Rdata")
load("jSVD_embeded.Rdata")
load("SNF_embeded.Rdata")

for(j in 1:10)
{
	b<-rec[[j]]
	aime.embeded<-b$embeded
	r2<-dist.grp(as.matrix(dist(aime.embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])
	if(j == 1) r<-r2
	else r<-cbind(r, r2[,2])
}

this<-dist.grp(as.matrix(dist(canon.embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(sca.embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(mofa.embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(icluster.embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(SNF.embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])
r<-cbind(r,this[,2]) 

colnames(r)<-c("k",rep("AIME (10 repeats)",10),"CCA","jSVD","MOFA2","iCluster","SNF")

cols<-c(rep("chocolate1",10),"cyan1","blue1","darkorchid1","red","green","black","black")
pchs<-c(rep(19,10),1,2,3,4,6)
cexes<-c(rep(1,10),1.5,1,1,1,1)

par(mfrow=c(1,2))
plot(r[,1], r[,2], xlab="k", ylab="proportion of neighbors in the same class", ylim=range(r[,-1]), type="b", col=cols[1], pch=pchs[1], cex=1)
for(i in 3:ncol(r)) lines(r[,1], r[,i], pch=pchs[i-1], col=cols[i-1], type="b", cex=cexes[i-1])

plot(1:8,1:8, type="n", axes=F, xlab="", ylab="")
for(i in 11:ncol(r))
{
	points(1,i-9, pch=pchs[i-1],col=cols[i-1], cex=cexes[i-1])
	text(1,i-9,colnames(r)[i],pos=4)
}

#####################

pdf("six method pairs.pdf", width=8, height=8)

colnames(aime.embeded)<-colnames(canon.embeded)<-colnames(sca.embeded)<-colnames(mofa.embeded)<-colnames(icluster.embeded)<-colnames(SNF.embeded)<-paste("Factor",1:4)
pairs(aime.embeded[!is.na(cell.type),], cex=0.5, col=all.colors[as.factor(cell.type[!is.na(cell.type)])], pch=19,main="AIME")
pairs(canon.embeded[!is.na(cell.type),], cex=0.5, col=all.colors[as.factor(cell.type[!is.na(cell.type)])], pch=19,main="CCA")
pairs(sca.embeded[!is.na(cell.type),], cex=0.5, col=all.colors[as.factor(cell.type[!is.na(cell.type)])], pch=19,main="jSVD")
pairs(mofa.embeded[!is.na(cell.type),], cex=0.5, col=all.colors[as.factor(cell.type[!is.na(cell.type)])], pch=19,main="MOFA2")
pairs(icluster.embeded[!is.na(cell.type),], cex=0.5, col=all.colors[as.factor(cell.type[!is.na(cell.type)])], pch=19,main="iCluster2")
pairs(SNF.embeded[!is.na(cell.type),], cex=0.5, col=all.colors[as.factor(cell.type[!is.na(cell.type)])], pch=19, main="SNF")

dev.off()