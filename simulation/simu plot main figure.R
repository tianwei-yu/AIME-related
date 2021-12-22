
library(MLmetrics)
library(pROC)

nfac=3
n.simu=10
nx=c(1000) # dimension of X
ny=c(1000) # dimension of Y
kx=c(40,20,10) # number of effective variables in X
ky=c(40,20,10) # number of effective variables in Y per each effective variable in X
links=c(2,1) # types of link function, 0 is linear only, 1 is half linear half mix, 2 is all mixed
sn=c(3.16) # signal to noise
N=c(200,1000, 5000,10000) # sample size
rho=c(0,0.3)

combos<-expand.grid(nx, ny, kx, ky, links, sn, N, rho)
colnames(combos)<-c("nx", "ny", 'kx', 'ky', 'links', 'sn', 'N',"rho")
to.delete<-which(combos[,3]*combos[,4]>combos[,2])
if(length(to.delete)>0) combos<-combos[-to.delete,]

sel<-which((combos[,3]==10 & combos[,4] == 10) | (combos[,3]==20 & combos[,4] == 20) |(combos[,3]==20 & combos[,4] == 40) |(combos[,3]==40 & combos[,4] == 20))
combos<-combos[sel,]

combos2<-cbind(combos, matrix(0, ncol=10, nrow=nrow(combos)))
colnames(combos2)[9:18]<-c("AIME AUC","Canonical AUC", "PLS AUC","iCluster AUC", "jSVD AUC", "AIME PRAUC","Canonical PRAUC", "PLS PRAUC","iCluster PRAUC", "jSVD PRAUC")

pdf("hahaha 3.pdf",width=12, height=12)
par(mfrow=c(4,4))

d<-new("list")


for(nnn in 1:nrow(combos))
{
	nx<-combos[nnn,1]
	ny<-combos[nnn,2]
	kx<-combos[nnn,3]
	ky<-combos[nnn,4]
	links<-combos[nnn,5]
	sn<-combos[nnn,6]
	N<-combos[nnn,7]
	rho<-combos[nnn,8]

	filename<-paste("nx", nx, 'ny', ny, 'kx', kx, 'ky', ky, 'links', links, 'sn', sn, 'N', N, "rho",rho, ".bin")

	load(filename)
	rec0<-rec
	
	filename2<-paste("MOFA dim 3 nx", nx, 'ny', ny, 'kx', kx, 'ky', ky, 'links', links, 'sn', sn, 'N', N, "rho",rho, ".bin")
	load(filename2)
	for(i in 1:length(rec0)) rec0[[i]]<-cbind(rec0[[i]], rec[[i]][,2])
	
	filename3<-paste("iCluster dim 3 nx", nx, 'ny', ny, 'kx', kx, 'ky', ky, 'links', links, 'sn', sn, 'N', N, "rho",rho, ".bin")
	load(filename3)
	for(i in 1:length(rec0)) rec0[[i]]<-cbind(rec0[[i]], rec[[i]][,2])
	
	filename4<-paste("sca nx", nx, 'ny', ny, 'kx', kx, 'ky', ky, 'links', links, 'sn', sn, 'N', N, "rho",rho, ".bin")
	load(filename4)
	for(i in 1:length(rec0)) rec0[[i]]<-cbind(rec0[[i]], rec[[i]][,2])

	#filename5<-paste("AIME dim 3 nx", nx, 'ny', ny, 'kx', kx, 'ky', ky, 'links', links, 'sn', sn, 'N', N, "rho",rho, ".bin")
	#load(filename5)
	#for(i in 1:length(rec0)) rec0[[i]][,1]<-rec[[i]]

	filename6<-paste("PLS CCA SS nx", nx, 'ny', ny, 'kx', kx, 'ky', ky, 'links', links, 'sn', sn, 'N', N, "rho",rho, ".bin")
	load(filename6)
	for(i in 1:length(rec0)) rec0[[i]][,2:3]<-rec[[i]][,3:4]

	
	rec<-rec0

	r.auc<-matrix(NA, ncol=6, nrow=length(rec))
	colnames(r.auc)<-colnames(rec[[1]])
	r.pr<-r.auc

	y<-rep(0, nrow(rec[[1]]))
	y[1:kx]<-1

	for(n in 1:length(rec))
	{
		a<-rec[[n]]
		if(!is.na(rec[[n]][1]))
		{
			for(i in c(1,2,3,4,5,6))
			{
				r.auc[n,i]<-auc(roc(y~a[,i], direction="<"))
				if(r.auc[n,i]<0.5) r.auc[n,i]<-0.5
				r.pr[n,i]<-PRAUC(a[,i],y)
			}
		}
	}
	main=paste(nx,ny,kx,ky,links,sn,N,rho)

	boxplot(r.auc, main=main)
	boxplot(r.pr,main=main)

	combos2[nnn,9:13]<-apply(r.auc,2,mean,na.rm=T)
	combos2[nnn,14:18]<-apply(r.pr,2,mean,na.rm=T)

	d[[nnn]]<-list(r.auc=r.auc, r.pr=r.pr)
}

dev.off()



n.simu=10
nx=c(1000) # dimension of X
ny=c(1000) # dimension of Y
kx=c(10, 20, 40) # number of effective variables in X
ky=c(10, 20, 40) # number of effective variables in Y per each effective variable in X
links=c(1,2) # types of link function, 0 is linear only, 1 is half linear half mix, 2 is all mixed
sn=c(3.16) # signal to noise
N=c(200, 1000, 5000,10000) # sample size
rho=c(0, 0.3)

kxky<-cbind(c(20,20,40), c(20,40,20))


#pdf(paste("hohoho",nx,".pdf"),width=3*4, height=1+(length(ky)*length(kx))*4)
#par(mfrow=c(1+(length(ky)*length(kx)), 3))
pdf(paste("hohoho 2",nx,".pdf"),width=3*2.5, height=4*2.5)
par(mfrow=c(4,3))
par(mar=c(4.5, 4, 3.5, 1.5) + 0.1)
for(k in length(links):1)
{
	for(m in 1:length(rho))
	{
		thisname<-"mixed"
		if(links[k]==2) thisname="nonlinear"
		for(i in 1:nrow(kxky))
		{
				kx<-kxky[i,1]
				ky<-kxky[i,2]


				a<-b<-new("list")
				p<-NULL
				labs<-NULL

				for(n in N)
				{
					s<-which(combos2[,1]==nx & combos2[,2]==ny & combos2[,3]==kx & combos2[,4]==ky & combos2[,5]==links[k] & combos2[,7]==n & combos2[,8]==rho[m])
					
					a[[length(a)+1]]<-NA
					a[[length(a)+1]]<-NA
					a[[length(a)+1]]<-d[[s]][[1]][,1]
					a[[length(a)+1]]<-d[[s]][[1]][,2]
					a[[length(a)+1]]<-d[[s]][[1]][,3]
					a[[length(a)+1]]<-d[[s]][[1]][,4]
					a[[length(a)+1]]<-d[[s]][[1]][,5]
					a[[length(a)+1]]<-d[[s]][[1]][,6]

					
					b[[length(b)+1]]<-NA
					b[[length(b)+1]]<-NA
					b[[length(b)+1]]<-d[[s]][[2]][,1]
					b[[length(b)+1]]<-d[[s]][[2]][,2]
					b[[length(b)+1]]<-d[[s]][[2]][,3]
					b[[length(b)+1]]<-d[[s]][[2]][,4]
					b[[length(b)+1]]<-d[[s]][[2]][,5]
					b[[length(b)+1]]<-d[[s]][[2]][,6]
										
					labs<-c(labs, paste("N",n,"rho",rho[m]))
					
					p<-c(p, NA, NA)
					message("---------------------------------------------------------------")
					message(c(kx,"  ",ky, "  ",links[k], "  ",n,"  ",rho[m]))
					for(j in 1:6)
					{
						message(wilcox.test(d[[s]][[2]][,j], d[[s]][[2]][,2])$p.value)
						p<-c(p, wilcox.test(d[[s]][[2]][,j], d[[s]][[2]][,2])$p.value)
					}
				}
				a<-a[-1:-2]
				b<-b[-1:-2]
				p<-p[-1:-2]
				
				border.cols<-c("chocolate4","cyan4","blue4","darkorchid4","darkred","darkgreen","black","black")
				fill.cols<-c("chocolate1","cyan1","blue1","darkorchid1","red","green","black","black")
				labs<-c(NA, NA, NA, 200, NA, NA, NA, NA, NA, NA, NA,1000,NA, NA, NA, NA, NA,NA, NA,5000,NA, NA, NA, NA, NA,NA, NA,10000,NA, NA)
				
				boxplot(b, col=fill.cols, border=border.cols, horizontal=F, ylim=c(0, 1.05),main=paste("rho=", rho[m], ", k=", kx,", m=",ky,", ",thisname,sep=""),names=labs,xlab="sample size",ylab="PR-AUC")
				
				mx<-lapply(b,max)
				lab<-rep("",length(b))
				lab[p<=0.01]<-"*"
				lab[p<=0.001]<-"**"
				#lab[p<=0.0001]<-"***"
				
				text(1:length(p),mx,lab, pos=3,col=border.cols, offset=0, cex=0.9)
				
				abline(v=c(7.5, 15.5, 23.5), col="grey")
		}

}}
dev.off()


