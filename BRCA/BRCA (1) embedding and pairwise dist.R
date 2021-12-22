library(AIME)
library(MLmetrics)
library(pROC)
source('fdrgamma.r')

library(keras)
library(ggplot2)
library(tidyverse)
library(e1071)

if(Sys.info()["nodename"] == "DESKTOP-0L4M44G") 
{
	prefix<-"D:/OneDrive"
	temp.path<-"C:\\Users\\temp\\Desktop"
	use_python("C:\\ProgramData\\Anaconda3\\envs\\tf-keras")
}


load("BRCA.miRseq_mature_RPM_log2.bin")
load("BRCA.rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.bin")

mir2<-mir
gene2<-gene

for(i in 1:ncol(mir2)) mir2[,i]<-(mir2[,i]-mean(mir2[,i]))/sd(mir2[,i])
for(i in 1:ncol(gene2)) gene2[,i]<-(gene2[,i]-mean(gene2[,i]))/sd(gene2[,i])

############
load("sample_info_mapped.bin")
load("clinical data mapped.bin")

T1<-rep(NA, nrow(clin))
T1[clin$Tumor..T1.Coded=="T1"]<-1
T1[clin$Tumor..T1.Coded=="T_Other"]<-0

ER<-rep(NA, nrow(clin))
ER[clin$ER.Status == "Positive"]<-1
ER[clin$ER.Status == "Negative"]<-0

Age<-clin$Diagnosis.Age

conf<-cbind(clin$Diagnosis.Age, T1, ER)

sel<-which(apply(is.na(conf),1,sum) == 0)

all.colors<-c("red","green","blue","cyan", "yellow","orange","grey50", "red","white","purple","darkblue","dodgerblue4","darkred","darkorange","darkcyan","magenta","firebrick","khaki4")


#g<-aime.select(data.in=gene[sel,], data.out=mir[sel,], confounder=conf[sel,], all.in.layers=2:5, all.out.layers=2:5, all.dropouts=c(0.2, 0.3,0.4, 0.5), repeats=3, temp.path=temp.path, col=all.colors[as.numeric(as.factor(ER))], cor.cut=0.5, kurtosis.cut=0.5, skew.cut=0.5)

#cbind(g,rank(-apply(g[,4:5],1,sum)))

#save(g, file="selection results.bin")

rec<-new("list")
for(n in 1:10)
{
   b<-aime(data.in=mir, data.out=gene, in.layers=4, out.layers=5,  max.epochs=100, max.dropout=0.2, importance.permutations=2, ncomp=4, pairwise.importance = TRUE)
   rec[[n]]<-b
   save(rec, file="mir to gene aime layers 4 5 dropout 0.2 no confounder.bin")
}

rec<-new("list")
for(n in 1:10)
{
   b<-aime(data.in=mir, data.out=gene, confounder=conf, in.layers=5, out.layers=4,  max.epochs=100, max.dropout=0.2, importance.permutations=2, ncomp=4, pairwise.importance = TRUE)
   rec[[n]]<-b
   save(rec, file="mir to gene aime layers 5 4 dropout 0.2 with confounder.bin")
}


######################## other methods
 

library(MOFA2)

data <- list(mir=t(mir2), gene=t(gene2))

lapply(data,dim)
N = ncol(data[[1]])
MOFAobject <- create_mofa(data)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors<-4
train_opts <- get_default_training_options(MOFAobject)
data_opts <- get_default_data_options(MOFAobject)
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
outfile = file.path(getwd(),"model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

model<-MOFAobject.trained

plot_factors(model, 
  factors = 1:4,
  color_by = anno[,5]
)
weights <- get_weights(MOFAobject.trained, views = "all", factors = "all")
x<-attributes(model)$expectations
mofa.embeded<-mir %*% x$W[[1]]


save(MOFAobject.trained, mofa.embeded, file="MOFA2 results.Rdata")

##############

library(mixOmics)
ppp<-pls(mir2, gene2, ncomp = 4, mode="canonical")
zzz<-mir2 %*% ppp$loadings$X
canon.embeded<-zzz
save(ppp, canon.embeded, file="CCA results.Rdata")

################### iCluster

library(iCluster)
data <- list(mir=(mir), gene=(gene))
results<-iCluster2(data, k=5)
icluster.embeded<-t(results$expZ)
save(results, icluster.embeded, file="iCluster results.Rdata")

################### jSVD
library(multiway) #for SCA == jSVD
data <- list(mir=t(mir2), gene=t(gene2))
s<-sca(data, 4)
sca.embeded<-s$B
save(s, sca.embeded, file="jSVD results.Rdata")

################### SNF

library(SNFtool)

## First, set all the parameters:
K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 20; 	# Number of Iterations, usually (10~20)

Data1<-mir2
Data2<-gene2

Dist1 = (dist2(as.matrix(Data1),as.matrix(Data1)))^(1/2)
Dist2 = (dist2(as.matrix(Data2),as.matrix(Data2)))^(1/2)

## next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)

## next, we fuse all the graphs
## then the overall matrix can be computed by similarity network fusion(SNF):
W = SNF(list(W1,W2), K, T)
e<-eigen(W)
SNF.embeded<-e$vec[,1:4]

save(W, SNF.embeded, file="SNF results.Rdata")

#####################################

########## plot the new methods
source('dist_group.r')

load("SNF results.Rdata")
load("iCluster results.Rdata")
load("jSVD results.Rdata")
load("MOFA2 results.Rdata")
load("CCA results.Rdata")
load("mir to gene aime layers 4 5 dropout 0.2 no confounder.bin")

ER2<-ER

for(j in 1:10)
{
	b<-rec[[j]]
	aime.embeded<-b$embeded
	r2<-dist.grp(as.matrix(dist(aime.embeded[!is.na(ER2),])), ER2[!is.na(ER2)])
	if(j == 1) r<-r2
	else r<-cbind(r, r2[,2])
}

this<-dist.grp(as.matrix(dist(canon.embeded[!is.na(ER2),])), ER2[!is.na(ER2)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(sca.embeded[!is.na(ER2),])), ER2[!is.na(ER2)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(mofa.embeded[!is.na(ER2),])), ER2[!is.na(ER2)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(icluster.embeded[!is.na(ER2),])), ER2[!is.na(ER2)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(SNF.embeded[!is.na(ER2),])), ER2[!is.na(ER2)])
r<-cbind(r,this[,2]) 

colnames(r)<-c("k",rep("AIME (10 repeats)",10),"CCA","jSVD","MOFA2","iCluster","SNF")

cols<-c(rep("chocolate1",10),"cyan1","blue1","darkorchid1","red","green","black","black")
pchs<-c(rep(19,10),1,2,3,4,6)
cexes<-c(rep(1,10),2,2,2,2,2)

par(mfrow=c(1,3))
plot(r[,1], r[,2], xlab="k", ylab="proportion of neighbors in the same class", ylim=range(r[,-1]), type="b", col=cols[1], pch=pchs[1], cex=1)
for(i in 3:ncol(r)) lines(r[,1], r[,i], pch=pchs[i-1], col=cols[i-1], type="b", cex=cexes[i-1])


T1<-clin$PAM50.Subtype
T1[T1==""]<-NA

for(j in 1:10)
{
	b<-rec[[j]]
	aime.embeded<-b$embeded
	r2<-dist.grp(as.matrix(dist(aime.embeded[!is.na(T1),])), T1[!is.na(T1)])
	if(j == 1) r<-r2
	else r<-cbind(r, r2[,2])
}

this<-dist.grp(as.matrix(dist(canon.embeded[!is.na(T1),])), T1[!is.na(T1)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(sca.embeded[!is.na(T1),])), T1[!is.na(T1)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(mofa.embeded[!is.na(T1),])), T1[!is.na(T1)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(icluster.embeded[!is.na(T1),])), T1[!is.na(T1)])
r<-cbind(r,this[,2]) 
this<-dist.grp(as.matrix(dist(SNF.embeded[!is.na(T1),])), T1[!is.na(T1)])
r<-cbind(r,this[,2]) 

colnames(r)<-c("k",rep("AIME (10 repeats)",10),"CCA","jSVD","MOFA2","iCluster","SNF")

cols<-c(rep("chocolate1",10),"cyan1","blue1","darkorchid1","red","green","black","black")
pchs<-c(rep(19,10),1,2,3,4,6)
cexes<-c(rep(1,10),2,2,2,2,2)

plot(r[,1], r[,2], xlab="k", ylab="proportion of neighbors in the same class", ylim=range(r[,-1]), type="b", col=cols[1], pch=pchs[1], cex=1)
for(i in 3:ncol(r)) lines(r[,1], r[,i], pch=pchs[i-1], col=cols[i-1], type="b", cex=cexes[i-1])


plot(1:8,1:8, type="n", axes=F, xlab="", ylab="")
for(i in 11:ncol(r))
{
	points(1,i-9, pch=pchs[i-1],col=cols[i-1], cex=cexes[i-1])
	text(1.25,i-9,colnames(r)[i],pos=4)
}


#####################

pdf("six method pairs.pdf", width=8, height=8)

all.colors<-c("orange","cyan","orange","blue", "yellow","orange","grey50", "red","white","purple","darkblue","dodgerblue4","darkred","darkorange","darkcyan","magenta","firebrick","khaki4")

colnames(aime.embeded)<-colnames(canon.embeded)<-colnames(sca.embeded)<-colnames(mofa.embeded)<-colnames(icluster.embeded)<-colnames(SNF.embeded)<-paste("Factor",1:4)
pairs(aime.embeded[!is.na(ER2),], cex=0.5, col=all.colors[as.factor(ER2[!is.na(ER2)])], pch=19,main="AIME")
pairs(canon.embeded[!is.na(ER2),], cex=0.5, col=all.colors[as.factor(ER2[!is.na(ER2)])], pch=19,main="CCA")
pairs(sca.embeded[!is.na(ER2),], cex=0.5, col=all.colors[as.factor(ER2[!is.na(ER2)])], pch=19,main="jSVD")
pairs(mofa.embeded[!is.na(ER2),], cex=0.5, col=all.colors[as.factor(ER2[!is.na(ER2)])], pch=19,main="MOFA2")
pairs(icluster.embeded[!is.na(ER2),], cex=0.5, col=all.colors[as.factor(ER2[!is.na(ER2)])], pch=19,main="iCluster2")
pairs(SNF.embeded[!is.na(ER2),], cex=0.5, col=all.colors[as.factor(ER2[!is.na(ER2)])], pch=19, main="SNF")

all.colors<-c("red","green","orange","blue", "yellow","orange","grey50", "red","white","purple","darkblue","dodgerblue4","darkred","darkorange","darkcyan","magenta","firebrick","khaki4")

pairs(aime.embeded[!is.na(T1),], cex=0.5, col=all.colors[as.factor(T1[!is.na(T1)])], pch=19,main="AIME")
pairs(canon.embeded[!is.na(T1),], cex=0.5, col=all.colors[as.factor(T1[!is.na(T1)])], pch=19,main="CCA")
pairs(sca.embeded[!is.na(T1),], cex=0.5, col=all.colors[as.factor(T1[!is.na(T1)])], pch=19,main="jSVD")
pairs(mofa.embeded[!is.na(T1),], cex=0.5, col=all.colors[as.factor(T1[!is.na(T1)])], pch=19,main="MOFA2")
pairs(icluster.embeded[!is.na(T1),], cex=0.5, col=all.colors[as.factor(T1[!is.na(T1)])], pch=19,main="iCluster2")
pairs(SNF.embeded[!is.na(T1),], cex=0.5, col=all.colors[as.factor(T1[!is.na(T1)])], pch=19, main="SNF")

dev.off()


