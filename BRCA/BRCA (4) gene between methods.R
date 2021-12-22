source('fdrgamma.r')

load("BRCA.miRseq_mature_RPM_log2.bin")
load("BRCA.rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.bin")

for(i in 1:ncol(mir))
{
	this<-strsplit(colnames(mir)[i], "\\|")
	colnames(mir)[i]<-this[[1]][1]
}


for(i in 1:ncol(gene))
{
	this<-strsplit(colnames(gene)[i], "\\|")
	colnames(gene)[i]<-this[[1]][2]

}

############
load("sample_info_mapped.bin")
load("clinical data mapped.bin")

###########################


all.colors<-c("black","green","orange","blue","yellow","cyan",  "red","grey50", "purple","darkblue","khaki4","darkred","darkorange","darkcyan","magenta","firebrick","dodgerblue4")


################
clin$ER.Status[!(clin$ER.Status %in% c("Positive", "Negative"))]<-NA
clin$PAM50.Subtype[clin$PAM50.Subtype==""]<-NA

for(i in 1:ncol(mir))
{
	this<-strsplit(colnames(mir)[i], "\\|")
	colnames(mir)[i]<-this[[1]][1]
}


####### (1) with confounders


load("mir to gene aime layers 5 4 dropout 0.2 with confounder.bin")
m.to.g<-rec
m.imp<-apply(m.to.g[[1]]$pair.imp,1,sum)
for(i in 2:length(m.to.g)) m.imp<-m.imp+apply(m.to.g[[i]]$pair.imp,1,sum)

m.fdr<-fdrgamma(m.imp)
sum(m.fdr<0.1)

AIME_adjusted<-colnames(gene)[m.fdr<0.1]

####### (2) no confounders

load("mir to gene aime layers 4 5 dropout 0.2 no confounder.bin")
m.to.g<-rec
m.imp2<-apply(m.to.g[[1]]$pair.imp,1,sum)
for(i in 2:length(m.to.g)) m.imp2<-m.imp2+apply(m.to.g[[i]]$pair.imp,1,sum)

m.cut<-quantile(m.imp2, (length(m.imp2)-270)/length(m.imp2))
sum(m.imp2>m.cut)

AIME_unadjusted<-colnames(gene)[m.imp2>m.cut]

###

load("CCA_normalized_results.bin")
weights.cca<-apply((ppp$loadings$Y)^2,1,sum)
w.cut<-quantile(weights.cca, (length(weights.cca)-270)/length(weights.cca))
CCA<-colnames(gene)[weights.cca>w.cut]

######## (4) MOFA2

load("MOFA2_model.bin")
weights.mofa <- get_weights(model, views = "all", factors = "all")
weights.mofa<-apply(abs(weights.mofa[[2]]),1,sum)
w.cut<-quantile(weights.mofa, (length(weights.mofa)-270)/length(weights.mofa))
MOFA2<-colnames(gene)[weights.mofa>w.cut]


#####

load("iCluster results.Rdata")
		weights.icluster<-results$W[243:nrow(results$W),]
		weights.icluster<-apply(weights.icluster^2,1,sum)
w.cut<-quantile(weights.icluster, (length(weights.icluster)-270)/length(weights.icluster))
iCluster<-colnames(gene)[weights.icluster>w.cut]

######


######## compare


library(eulerr)
plot(venn(list(AIME_adjusted=AIME_adjusted,MOFA2=MOFA2, AIME_unadjusted=AIME_unadjusted, iCluster2=iCluster)))

###


gg2<-list(AIME_adjusted=m.imp,MOFA2=weights.mofa, AIME_unadjusted=m.imp2, iCluster=weights.icluster)

newnames<-names<-names(gg2[[2]])
for(i in 1:length(names)) 
{
	this.name<-strsplit(names[i],'\\|')[[1]][2]
	if(is.na(this.name)) this.name<-paste("isna",i)
	newnames[i]<-this.name
}
for(i in 1:4)  names(gg2[[i]])<-newnames

library(fgsea)
library(reactome.db)
pathways <- reactomePathways(names(gg2[[1]]))

pa<-new("list")
for(i in 1:4)
{
	f<-fgseaSimple(pathways, gg2[[i]],nperm=5000)
	pa[[i]]<-f
}
names(pa)<-names(gg2)

save(pa, file="four methods all genes fgsea.bin")

pa2<-pa[[1]][,c(1,7,3)]
for(i in 2:4) pa2<-cbind(pa2, pa[[i]][,3])
colnames(pa2)[3:6]<-names(gg2)
pa2[is.na(pa2)]<-1
l<-apply(pa2[,3:6],1,min)
pa2<-pa2[l<=0.1,]
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
		if(pa2[i,j]<=0.1) panet<-rbind(panet, c(rownames(pa2)[i], colnames(pa2)[j]))
	}
}
panet<-panet[-1,]
write.table(panet, "pathway net 4 methods FDR 0.1.txt", sep="\t", quote=F, col.names=F, row.names=F)

### to be followed by cytoscape