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

T1<-rep(NA, nrow(clin))
T1[clin$Tumor..T1.Coded=="T1"]<-1
T1[clin$Tumor..T1.Coded=="T_Other"]<-0

ER<-rep(NA, nrow(clin))
ER[clin$ER.Status == "Positive"]<-1
ER[clin$ER.Status == "Negative"]<-0

Age<-clin$Diagnosis.Age

conf<-cbind(clin$Diagnosis.Age, T1, ER)

sel<-which(apply(is.na(conf),1,sum) == 0)



all.colors<-c("black","green","orange","blue","yellow","cyan",  "red","grey50", "purple","darkblue","khaki4","darkred","darkorange","darkcyan","magenta","firebrick","dodgerblue4")


################

load("mir to gene aime layers 5 4 dropout 0.2 with confounder.bin")
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


# overall imiportance

m.imp<-apply(m.imp,1,mean)
m.fdr<-fdrgamma(m.imp)

#mg.fdr<-fdrgamma(mg)
#mg.fdr<-matrix(mg.fdr, nrow=nrow(mg))
#colnames(mg.fdr)<-colnames(mg)
#rownames(mg.fdr)<-rownames(mg)
#mg.fdr.together<-mg.fdr
#save(mg.fdr.together, file="mg_fdr_all_together.Rdata")

#mg.fdr<-mg
#for(i in 1:ncol(mg)) 
#{
#	zzz<-fdrgamma(mg[,i])
#	sel<-which(mg[,i]>=quantile(mg[,i], 1-10/nrow(mg)))
#	zzz[sel]<- 0
#	mg.fdr[,i]<-zzz
#}
#save(mg.fdr, file="mg_fdr_one_by_one.Rdata")

#mg.2<-mg

#library(multiMiR)
#gg<-get_multimir(mirna = colnames(mg.2))
#gg2<-get_multimir(mirna = colnames(mg.2),table="predicted")
#g.known<-gg@data
#g.known2<-gg2@data
#g.known<-rbind(g.known[,1:6], g.known2[,1:6])
#g.known<-unique(g.known)

#mg.2<-mg.2*0
#for(i in 1:nrow(g.known))
#{
#	mg.2[which(rownames(mg.2) == g.known[i,5]),which(colnames(mg.2) == g.known[i,3])]<-1
#}
#mg.target.all<-mg.2
#save(mg.target.all, file="mg_target_all.Rdata")

load("mg_fdr_one_by_one.Rdata")
load("mg_fdr_all_together.Rdata")
load("mg_target_all.Rdata")

m.fdr.cuts<-c(0.01, 0.05, 0.1, 0.2)
fdr.cols<-c("blue","green","cyan","purple","orange","red")
for(j in 1:length(m.fdr.cuts))
{
	s<-which(m.fdr<m.fdr.cuts[j])  ### 0.01 yields 30, 0.001 yields 8
	print(c(m.fdr.cuts[j], length(s)))
	s2<-which(apply(mg.target.all,1,sum)>0)
	x<-as.vector(mg.fdr[,s])
	y<-as.vector(mg.target.all[,s])
	xrand<-sample(x,length(x),replace=F)

	fdr.cuts<-seq(log(1e-3), log(1), length.out=100)
	fdr.cuts<-exp(fdr.cuts)
	rr<-fdr.cuts
	for(i in 1:length(fdr.cuts))
	{
		this.x<-1*(x<=fdr.cuts[i])
		print(c(fdr.cuts[i], sum(this.x), sum(this.x*y), (sum(this.x*y)/sum(this.x))/(sum(y)/length(y))))
		rr[i]<-(sum(this.x*y)/sum(this.x))/(sum(y)/length(y))
	}
	
	if(j == 1)
	{
		plot(fdr.cuts, rr, type="l",col=fdr.cols[j],ylim=c(1, 2),xlab="miRNA-gene pair fdr cutoff", ylab="fold-change of %validated miRNA-gene pairs over random pairs")
	}else{
		lines(fdr.cuts, rr, col=fdr.cols[j])
	}
}

top.pos=2
shift.pos=0.12
x.pos<-c(0.4, 0.45)

for(i in 1:length(m.fdr.cuts))
{
	top.pos<-top.pos-shift.pos
	lines(x.pos, c(top.pos,top.pos), col=fdr.cols[i])
	text(x.pos[2], top.pos, pos=4, paste("miRNA fdr cutoff: ", m.fdr.cuts[i],", selected miRNAs: ", sum(m.fdr<m.fdr.cuts[i])))
}


##### generate graph

s<-which(m.fdr<1e-2)  
x<-mg.fdr[,s]
y<-mg.target.all[,s]
this.x<-1*(x<=1e-3)

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

library(org.Hs.eg.db)
genes<-mget(m.g.table[,3], org.Hs.egSYMBOL, ifnotfound=NA)
for(i in 1:length(genes)) if(length(genes[[i]])>1) genes[[i]]<-genes[[i]][1]
m.g.table<-cbind(m.g.table, unlist(genes))
m.g.table<-m.g.table[!is.na(m.g.table[,4]),]
write.table(m.g.table, "miRNA_gene_fdr 0.01 0.001.txt",sep="\t", col.names=T, row.names=F, quote=F)

mir.valid.label<-rbind(cbind(unique(m.g.table[,2]), rep("mir", length(unique(m.g.table[,2])))), m.g.table[m.g.table[,1]=="1",c(4,1)])
write.table(mir.valid.label, "mir_valid_relation_labels_for_cytoscape.txt",sep="\t", col.names=T, row.names=F, quote=F)
