source('fdrgamma.r')

load("BRCA.miRseq_mature_RPM_log2.bin")
load("BRCA.rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.bin")

############
load("sample_info_mapped.bin")
load("clinical data mapped.bin")

###########################

T1<-rep(NA, nrow(clin))
T1[clin$Tumor..T1.Coded=="T1"]<-1
T1[clin$Tumor..T1.Coded=="T_Other"]<-0

ER<-rep(NA, nrow(clin))
ER[clin$ER.Status == "Positive"]<-1
ER[clin$ER.Status == "Negative"]<-0

Age<-clin$Diagnosis.Age

conf<-cbind(clin$Diagnosis.Age, T1, ER)


all.colors<-c("black","green","orange","blue","yellow","cyan",  "red","grey50", "purple","darkblue","khaki4","darkred","darkorange","darkcyan","magenta","firebrick","dodgerblue4")


################
clin$ER.Status[!(clin$ER.Status %in% c("Positive", "Negative"))]<-NA
clin$PAM50.Subtype[clin$PAM50.Subtype==""]<-NA

for(i in 1:ncol(mir))
{
	this<-strsplit(colnames(mir)[i], "\\|")
	colnames(mir)[i]<-this[[1]][1]
}


########

##### (1) with confounders

load("mir to gene aime layers 5 4 dropout 0.2 with confounder.bin")
m.imp<-rec[[1]]$imp
for(i in 2:length(rec)) m.imp<-m.imp+rec[[i]]$imp

m.fdr<-fdrgamma(m.imp)
sum(m.fdr<0.1)

AIME_adjusted<-colnames(mir)[m.fdr<0.1]

####### (2) no confounders

load("mir to gene aime layers 4 5 dropout 0.2 no confounder.bin")
m.imp<-rec[[1]]$imp
for(i in 2:length(rec)) m.imp<-m.imp+rec[[i]]$imp

#m.cut=fdrgamma.cut(m.imp,use.percentile=NA, fdr.cut=0.1)
m.cut<-quantile(m.imp, (length(m.imp)-16)/length(m.imp))
sum(m.imp>m.cut)

AIME_unadjusted<-colnames(mir)[m.imp>m.cut]


######## (3) CCA

load("CCA_normalized_results.bin")
weights<-apply((ppp$loadings$X)^2,1,sum)
w.cut<-quantile(weights, (length(weights)-16)/length(weights))
CCA<-colnames(mir)[weights>w.cut]


######## (4) MOFA2

load("MOFA2_model.bin")
weights <- get_weights(model, views = "all", factors = "all")
weights<-apply(abs(weights[[1]]),1,sum)
w.cut<-quantile(weights, (length(weights)-16)/length(weights))
MOFA2<-colnames(mir)[weights>w.cut]


######## (5) iCluster2

load("iCluster results.Rdata")
		weights.icluster<-results$W[1:242,]
		weights.icluster<-apply(weights.icluster^2,1,sum)
w.cut<-quantile(weights.icluster, (length(weights.icluster)-16)/length(weights.icluster))
iCluster<-colnames(mir)[weights.icluster>w.cut]


######## compare

library(eulerr)
plot(venn(list(AIME_adjusted=AIME_adjusted,MOFA2=MOFA2, AIME_unadjusted=AIME_unadjusted, iCluster=iCluster)))


###

write.table(AIME_adjusted, "mir_sel_AIME_adjusted.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(AIME_unadjusted, "mir_sel_AIME_unadjusted.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(CCA, "mir_sel_CCA.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(MOFA2, "mir_sel_MOFA2.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(iCluster, "mir_sel_iCluster.txt", sep="\t", quote=F, col.names=F, row.names=F)

############################## the following is after submiting to mirPath and downloading results

a<-read.csv("mirpath mir_sel_AIME_adjusted.csv")

BHFDR<-p.adjust(a[,2], method="BH")
a[,2]<-BHFDR
a<-a[,c(1,4,2)]
met<-rep("AIME_adjusted",nrow(a))
a<-cbind(met,a)

b<-read.csv("mirpath mir_sel_AIME_unadjusted.csv")
BHFDR<-p.adjust(b[,2], method="BH")
b[,2]<-BHFDR
b<-b[,c(1,4,2)]
met<-rep("AIME_unadjusted",nrow(b))
b<-cbind(met,b)
a<-rbind(a,b)

b<-read.csv("mirpath mir_sel_CCA.csv")
BHFDR<-p.adjust(b[,2], method="BH")
b[,2]<-BHFDR
b<-b[,c(1,4,2)]
met<-rep("CCA",nrow(b))
b<-cbind(met,b)
a<-rbind(a,b)

b<-read.csv("mirpath mir_sel_iCluster.csv")
BHFDR<-p.adjust(b[,2], method="BH")
b[,2]<-BHFDR
b<-b[,c(1,4,2)]
met<-rep("iCluster2",nrow(b))
b<-cbind(met,b)
a<-rbind(a,b)

b<-read.csv("mirpath mir_sel_MOFA2.csv")
BHFDR<-p.adjust(b[,2], method="BH")
b[,2]<-BHFDR
b<-b[,c(1,4,2)]
met<-rep("MOFA2",nrow(b))
b<-cbind(met,b)
a<-rbind(a,b)

a2<-a[a[,4]<=0.025,]

write.table(a2, "five methods mir pathway FDR 0.025.txt",sep="\t",col.names=T, row.names=F, quote=F)
