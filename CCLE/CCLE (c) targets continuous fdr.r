setwd('C:/users/tianw/Onedrive/1-work/autoencoder metaanalysis/CCLE')
setwd('D:/Onedrive/1-work/autoencoder metaanalysis/CCLE')
library(AIME)
library(MLmetrics)
library(pROC)
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

########### my own FDR with gamma null

m.imp<-apply(m.imp,1,mean)
m.fdr<-fdrgamma(m.imp)

rownames(mg)<-substr(rownames(mg),1,15)

mir.names<-read.table("mir names A-GEOD-16231_comments.txt",sep="\t", header=T)
for(i in 1:ncol(mg))
{
	sel<-which(mir.names[,4] == colnames(mg)[i])
	#print(c(i,sel))
	colnames(mg)[i]<-mir.names[sel,1]
}


##########################

mg.2<-mg

library(multiMiR)
gg<-get_multimir(mirna = colnames(mg.2))
g.known<-gg@data
g.known<-g.known[which(g.known[,6] %in% rownames(mg.2)),]
gg2<-get_multimir(mirna = colnames(mg.2),table="predicted")
g.known2<-gg2@data
g.known2<-g.known2[which(g.known2[,6] %in% rownames(mg.2)),]
g.known<-rbind(g.known[,c(2:6)], g.known2[,c(2:6)])
g.known<-unique(g.known)

mg.2<-mg.2*0
for(i in 1:nrow(g.known))
{
	mg.2[which(rownames(mg.2) == g.known[i,5]),which(colnames(mg.2) == g.known[i,1])]<-1
}
mg.target.all<-mg.2

mg.fdr<-mg
for(i in 1:ncol(mg)) 
{
	zzz<-fdrgamma(mg[,i])
	sel<-which(mg[,i]>=quantile(mg[,i], 1-10/nrow(mg)))
	zzz[sel]<- 0
	mg.fdr[,i]<-zzz
}
save(mg.fdr, file="mg_fdr_one_by_one.Rdata")

load("mg_fdr_one_by_one.Rdata")
load("mg_target_all.Rdata")

m.fdr.cuts<-c(1e-7, 1e-5, 0.001, 0.01, 0.1, 0.2)
fdr.cols<-c("blue","green","cyan","purple","orange","red","grey")
for(j in 1:length(m.fdr.cuts))
{
	s<-which(m.fdr<m.fdr.cuts[j])  
	print(c(m.fdr.cuts[j], length(s)))
	s2<-which(apply(mg.target.all,1,sum)>0)
	x<-as.vector(mg.fdr[,s])
	y<-as.vector(mg.target.all[,s])
	xrand<-sample(x,length(x),replace=F)

	fdr.cuts<-seq(log(1e-5), log(1), length.out=100)
	fdr.cuts<-exp(fdr.cuts)
	rr<-fdr.cuts
	for(i in 1:length(fdr.cuts))
	{
		this.x<-1*(x<=fdr.cuts[i])
		rr[i]<-(sum(this.x*y)/sum(this.x))/(sum(y)/length(y))
	}
	
	if(j == 1)
	{
		plot(fdr.cuts, rr, type="l",col=fdr.cols[j],ylim=c(1, 4),log="x",xlab="miRNA-gene pair fdr cutoff", ylab="fold-change of %validated miRNA-gene pairs over random pairs")
	}else{
		lines(fdr.cuts, rr, col=fdr.cols[j])
	}
}

top.pos=4.1
shift.pos=0.15
x.pos<-c(1e-4, 2e-4)

for(i in 1:length(m.fdr.cuts))
{
	top.pos<-top.pos-shift.pos
	lines(x.pos, c(top.pos,top.pos), col=fdr.cols[i])
	text(x.pos[2], top.pos, pos=4, paste("miRNA fdr cutoff: ", m.fdr.cuts[i],", selected miRNAs: ", sum(m.fdr<m.fdr.cuts[i])))
}


##### generate graph

s<-which(m.fdr<1e-7)  
x<-mg.fdr[,s]
y<-mg.target.all[,s]
this.x<-1*(x<=1e-5)

mg.3<-this.x
m.g.graph<-c(0,0,0)
for(i in 1:ncol(mg.3))
{
	for(j in 1:nrow(mg.3))
	{
		if(mg.3[j,i]==1 & y[j,i]==1) m.g.graph<-rbind(m.g.graph, c(1, colnames(mg.3)[i], rownames(mg.3)[j]))
		if(mg.3[j,i]==1 & y[j,i]==0) m.g.graph<-rbind(m.g.graph, c(0, colnames(mg.3)[i], rownames(mg.3)[j]))
	}
}
m.g.graph<-m.g.graph[-1,]
m.g.table<-m.g.graph

### output with gene name

library(org.Hs.eg.db)
g<-m.g.table[,3]
for(i in 1:length(g))
{
 g[i]<-strsplit(g[i],"\\.")[[1]][1]
}
m.g.table[,3]<-g
ents<-mget(g,org.Hs.egENSEMBL2EG, ifnotfound=NA)
for(i in 1:length(ents)) if(length(ents[[i]])>1) ents[[i]]<-ents[[i]][1]
ents<-unlist(ents)
genes<-mget(ents[!is.na(ents)], org.Hs.egSYMBOL, ifnotfound=NA)
for(i in 1:length(genes)) if(length(genes[[i]])>1) genes[[i]]<-genes[[i]][1]
genes<-unlist(genes)
ents[!is.na(ents)]<-genes
m.g.table<-cbind(m.g.table, unlist(ents))
m.g.table<-m.g.table[!is.na(m.g.table[,4]),]

write.table(m.g.table, "miRNA fdr 1e-7_gene_fdr 1e-5.txt",sep="\t", col.names=F, row.names=T, quote=F) # to be followed up by cytoscape




