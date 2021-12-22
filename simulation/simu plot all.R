library(MLmetrics)
library(pROC)

nfac=3
n.simu=10
nx=c(1000) # dimension of X
ny=c(1000) # dimension of Y
kx=c(40,20,10) # number of effective variables in X
ky=c(40,20,10) # number of effective variables in Y per each effective variable in X
links=c(2,1,0) # types of link function, 0 is linear only, 1 is half linear half mix, 2 is all mixed
sn=c(3.16) # signal to noise
N=c(200,1000, 5000,10000) # sample size
rho=c(0,0.3,0.6)

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
links=c(0,1,2) # types of link function, 0 is linear only, 1 is half linear half mix, 2 is all mixed
sn=c(3.16) # signal to noise
N=c(200, 1000, 5000,10000) # sample size
rho=c(0, 0.3, 0.6)

kxky<-cbind(c(10,20,20,40), c(10,20,40,20))


#pdf(paste("hohoho",nx,".pdf"),width=3*4, height=1+(length(ky)*length(kx))*4)
#par(mfrow=c(1+(length(ky)*length(kx)), 3))
pdf(paste("hohoho 2",nx,".pdf"),width=3*4, height=6*4)
par(mfrow=c(5,3))

for(i in 1:nrow(kxky))
{
	kx<-kxky[i,1]
	ky<-kxky[i,2]
		for(k in length(links):1)
		{
			

			a<-b<-new("list")
			labs<-NULL
			for(n in N)
			{
				for(m in 1:length(rho))
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
				}
			}

			a<-a[-1]
			b<-b[-1]
			#names(a)<-rep(c("AIME","CCA","PLS",""),100)[1:length(a)]
			#names(b)<-rep(c("AIME","CCA","PLS",""),100)[1:length(a)]

			thisname<-"all linear"
			if(links[k]==1) thisname<-"mixed"
			if(links[k]==2) thisname<-"all nonlinear"

			#border.cols<-c("chocolate4","cyan4","blue4","black")
			#fill.cols<-c("chocolate1","cyan1","blue1","black")
			#boxplot(a, col=fill.cols, border=border.cols, horizontal=TRUE, ylim=c(0.35, 1),main=paste("ROC-AUC, k=", kx[i],", m=",ky[j],", ",thisname,sep=""))
			#abline(h=seq(4,length(a), by=4),col="grey")
			#text(rep(0.3, length(labs)), seq(2, length(a),by=4), labs, pos=4)


			border.cols<-c("black","chocolate4","cyan4","blue4","darkorchid4","darkred","darkgreen","black")
			fill.cols<-c("black","chocolate1","cyan1","blue1","darkorchid1","red","green","black")
			#minpos<- -0.4* max(unlist(b),na.rm=T)
			#boxplot(b, col=fill.cols, border=border.cols, horizontal=TRUE, ylim=c(minpos, max(unlist(b),na.rm=T)),main=paste("PR-AUC, k=", kx,", m=",ky,", ",thisname,sep=""),names=NA,xlab="PR-AUC",ylab="Settings")
			minpos<- -0.4
			boxplot(b, col=fill.cols, border=border.cols, horizontal=TRUE, ylim=c(minpos, 1),main=paste("PR-AUC, k=", kx,", m=",ky,", ",thisname,sep=""),names=NA,xlab="PR-AUC",ylab="Settings")

			abline(h=seq(8,length(b), by=8),col="grey75",lty=3)
			abline(h=seq(24, length(b),by=24),col="black")
			text(rep(minpos, length(labs)), seq(2.5, length(a),by=8), labs, pos=4)

		}

}


plot(1:5,1:5, type="n", xlab="",ylab="",main="",axes=FALSE)
points(2.8,5,pch=15,cex=3,col="chocolate1")
text(3,5,"AIME",pos=4)
points(4.3,5,pch=15,cex=3,col="cyan1")
text(4.5,5,"CCA",pos=4)
plot(1:5,1:5, type="n", xlab="",ylab="",main="",axes=FALSE)
points(1.3,5,pch=15,cex=3,col="blue1")
text(1.5,5,"PLS",pos=4)
points(2.8,5,pch=15,cex=3,col="darkorchid1")
text(3,5,"MOFA2",pos=4)
plot(1:5,1:5, type="n", xlab="",ylab="",main="",axes=FALSE)
points(1.3,5,pch=15,cex=3,col="red")
text(1.5,5,"iCluster",pos=4)
points(2.8,5,pch=15,cex=3,col="green")
text(3,5,"jSVD",pos=4)
dev.off()









######################


n.simu=10
nx=c(1000) # dimension of X
ny=c(1000) # dimension of Y
kx=c(10, 20) # number of effective variables in X
ky=c(10, 20) # number of effective variables in Y per each effective variable in X
links=c(0,1,2) # types of link function, 0 is linear only, 1 is half linear half mix, 2 is all mixed
sn=c(5) # signal to noise
N=c(200,1000) # sample size
rho=c(0.6, 0.3, 0)


pdf("hohoho.pdf",width=length(ky)*5, height=(1+length(kx))*2.5)
par(mfrow=c(length(kx)+1, 2*length(ky)))

for(i in 1:length(kx))
{
	for(j in 1:length(ky))
	{
		a<-combos2[combos2[,3]==kx[i] & combos2[,4]==ky[j],]
		for(ind in c(9,12))
		{
			ylim=c(0,1)
			ylab="PR-AUC"
			if(ind == 9)
			{
				ylim=c(0.5,1)
				ylab="ROC-AUC"
			}
			plot(a[,7], a[,ind], log="x", ylim=ylim, ylab=ylab, xlab="Sample size", type="n",main=paste(ylab, ", K=", kx[i], ", m=", ky[j], sep=""))
			cols=c("red","blue","green")
			pchs<-c(19, 2, 3)
			rho.cex<-c(3,2,1)

			for(m in 1:length(rho))
			{
				a2<-a[a[,8]==rho[m],]
				for(addi in c(2,1,0))
				{
					points(a2[,7], a2[,ind+addi], pch=pchs[addi+1], col=cols[addi+1], cex=rho.cex[m])
				}

				for(addi in c(2,1,0))
				{
					for(k in length(links):1)
					{
						b<-a2[a2[5]==links[k],]
						lines(b[,7], b[,ind+addi], col=cols[addi+1],lty=k)
					}
				}
			}
		}
	}
}


plot(1,1,type="n",axes=F,xlab="",ylab="")

plot(1:6,1:6,type="n",axes=F,xlab="",ylab="")
lines(c(1,2),c(2.5,2.5), lty=1)
lines(c(1,2),c(3.5,3.5), lty=2)
lines(c(1,2),c(4.5,4.5), lty=3)
text(2,2.5,"100%", pos=4, cex=1.25)
text(2,3.5," 50%", pos=4, cex=1.25)
text(2,4.5,"  0%", pos=4, cex=1.25)
text(1,5.5,"Percent non-linear",pos=4, cex=1.5)

plot(1:6,1:6,type="n",axes=F,xlab="",ylab="")
#points(1.1,1.5,cex=2,pch=1,col="red")
points(1.1,2.5,cex=2,pch=19,col="red")
points(1.1,3.5,cex=2,pch=2,col="blue")
points(1.1,4.5,cex=2,pch=3,col="green")
#text(1.5,1.5,"AIME iterative", pos=4, cex=1.25)
text(1.5,2.5,"AIME", pos=4, cex=1.25)
text(1.5,3.5,"Canonical", pos=4, cex=1.25)
text(1.5,4.5,"PLS", pos=4, cex=1.25)
text(1.5,5.5,"Method",pos=4, cex=1.5)


dev.off()

