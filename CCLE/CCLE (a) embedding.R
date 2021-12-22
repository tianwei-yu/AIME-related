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

############ AIME 

#g<-aime.select(data.in=mir, data.out=gene, all.in.layers=2:5, all.out.layers=2:5, all.dropouts=c(0.2, 0.3,0.4, 0.5), repeats=3, cor.cut=0.5, kurtosis.cut=0.5, skew.cut=0.5)
#cbind(g,rank(-apply(g[,4:5],1,sum)))
#save(g, file="mir gene selection results.bin")

rec<-new("list")
for(n in 1:10)
{
   b<-aime(data.in=mir, data.out=gene, in.layers=3, out.layers=4,  max.epochs=100, max.dropout=0.3, importance.permutations=2, ncomp=4, pairwise.importance = TRUE)
   rec[[n]]<-b
   save(rec, file="mir to gene aime layers 3 4 dropout 0.3 no confounder.bin")
}

b<-rec[[1]]

aime.embeded<-b$embeded

pairs(b$embeded, cex=0.5,pch=19, col=all.colors[as.numeric(as.factor(cell.type))])
dist.grp(as.matrix(dist(b$embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])

#########


#############
col.codes<-data.frame(all.colors[as.numeric(as.factor(cell.type))], as.factor(cell.type))
col.codes<-unique(col.codes)
col.codes<-col.codes[!is.na(col.codes[,2]),]

plot(1:20,1:20,type="n")
for(i in 1:nrow(col.codes))
{
	points(1,i,pch=19, cex=2, col=col.codes[i,1])
	text(1.25, i, col.codes[i,2], pos=4)
}

#### limit to lung and blood (two largest types),taking sub-types with 10 or more samples

colors<-rep(NA,length(cell.type))

colors[anno[,10] == "acute_lymphoblastic_B_cell_leukaemia" & anno[,5]=="haematopoietic_and_lymphoid_tissue"]<-"firebrick1" 
colors[anno[,10] == "acute_lymphoblastic_T_cell_leukaemia" & anno[,5]=="haematopoietic_and_lymphoid_tissue"]<-"deeppink4"
colors[anno[,10] == "acute_myeloid_leukaemia" & anno[,5]=="haematopoietic_and_lymphoid_tissue"]<-"orange"
colors[anno[,10] == "blast_phase_chronic_myeloid_leukaemia" & anno[,5]=="haematopoietic_and_lymphoid_tissue"]<-"yellow"
colors[anno[,10] == "Burkitt_lymphoma" & anno[,5]=="haematopoietic_and_lymphoid_tissue"]<-"burlywood3"
colors[anno[,10] == "diffuse_large_B_cell_lymphoma" & anno[,5]=="haematopoietic_and_lymphoid_tissue"]<-"darkgoldenrod"
colors[anno[,10] == "plasma_cell_myeloma" & anno[,5]=="haematopoietic_and_lymphoid_tissue"]<-"darkred"
colors[anno[,10] == "adenocarcinoma" & anno[,5]=="lung"]<-"blue"
colors[anno[,10] == "large_cell_carcinoma" & anno[,5]=="lung"]<-"darkslateblue"
colors[anno[,10] == "non_small_cell_carcinoma" & anno[,5]=="lung"]<-"midnightblue"
colors[anno[,10] == "small_cell_carcinoma" & anno[,5]=="lung"]<-"purple"
colors[anno[,10] == "squamous_cell_carcinoma" & anno[,5]=="lung"]<-"cyan"

pairs(b$embeded[!is.na(colors),], cex=0.5, col=colors[!is.na(colors)], pch=19)

library(rgl)
plot3d(b$embeded[!is.na(colors),1], b$embeded[!is.na(colors),2], b$embeded[!is.na(colors),4], col=colors[!is.na(colors)],size=6,axes=FALSE)

pairs((mir %*% ppp$loadings$X)[!is.na(colors),], cex=0.5, col=colors[!is.na(colors)], pch=19)
#pairs((gene %*% ppp$loadings$Y)[!is.na(colors),], cex=0.75, col=colors[!is.na(colors)], pch=19)

col.codes<-matrix(c(
"acute_lymphoblastic_B_cell_leukaemia","firebrick1",
"acute_lymphoblastic_T_cell_leukaemia","deeppink4",
"acute_myeloid_leukaemia","orange",
"blast_phase_chronic_myeloid_leukaemia","yellow",
"Burkitt_lymphoma","burlywood3",
"diffuse_large_B_cell_lymphoma","darkgoldenrod",
"plasma_cell_myeloma","darkred",
"adenocarcinoma","blue",
"large_cell_carcinoma","darkslateblue",
"non_small_cell_carcinoma","midnightblue",
"small_cell_carcinoma","purple",
"squamous_cell_carcinoma","cyan"
),nrow=2)

col.codes<-cbind(col.codes[2,], col.codes[1,])

plot(1:20,1:20,type="n")
for(i in 1:nrow(col.codes))
{
	points(1,i,pch=19, cex=2, col=col.codes[i,1])
	text(1.25, i, col.codes[i,2], pos=4)
}


########

library(MOFA2)

data <- list(mir=t(mir2), gene=t(gene2))

lapply(data,dim)
N = ncol(data[[1]])
#groups = c(rep("A",N/2), rep("B",N/2))
#groups=anno[,5]
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

pairs(mofa.embeded[!is.na(cell.type),], cex=0.5, col=colors[!is.na(cell.type)], pch=19)

dist.grp(as.matrix(dist(mofa.embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])


##############

library(mixOmics)
ppp<-pls(mir2, gene2, ncomp = 4, mode="canonical")
zzz<-mir2 %*% ppp$loadings$X
pairs(zzz, cex=.5, pch=1, col=all.colors[as.numeric(as.factor(cell.type))])
dist.grp(as.matrix(dist((mir2 %*% ppp$loadings$X)[!is.na(cell.type),])), cell.type[!is.na(cell.type)])

################### iCluster

library(iCluster)
data <- list(mir=(mir2), gene=(gene2))
results<-iCluster2(data, k=5)
#results<-iCluster(data, k=5, lambda=c(.05,.05))
icluster.embeded<-t(results$expZ)
pairs(icluster.embeded[!is.na(cell.type),], cex=0.5, col=all.colors[as.factor(cell.type[!is.na(cell.type)])], pch=19)

dist.grp(as.matrix(dist(icluster.embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])

################### jSVD

library(multiway) #for SCA == jSVD
data <- list(mir=t(mir2), gene=t(gene2))
s<-sca(data, 4)
sca.embeded<-s$B
dist.grp(as.matrix(dist(sca.embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])
pairs(sca.embeded[!is.na(cell.type),], cex=0.5, col=all.colors[as.factor(cell.type[!is.na(cell.type)])], pch=19)

################### SNF

library(SNFtool)

K = 20;		# number of neighbors, usually (10~30)
alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
T = 20; 	# Number of Iterations, usually (10~20)

Data1<-mir2
Data2<-gene2

Dist1 = (dist2(as.matrix(Data1),as.matrix(Data1)))^(1/2)
Dist2 = (dist2(as.matrix(Data2),as.matrix(Data2)))^(1/2)
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)

W = SNF(list(W1,W2), K, T)
e<-eigen(W)
SNF.embeded<-e$vec[,1:4]
dist.grp(as.matrix(dist(SNF.embeded[!is.na(cell.type),])), cell.type[!is.na(cell.type)])
pairs(SNF.embeded[!is.na(cell.type),], cex=0.5, col=all.colors[as.factor(cell.type[!is.na(cell.type)])], pch=19)
