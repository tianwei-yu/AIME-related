
dist.grp.majority<-function(distmat, grp, k=1:20)
{
	
	grpmat<-matrix(0,nrow=length(grp),ncol=length(grp))
	for(i in 1:length(grp)) for(j in 1:length(grp)) if(grp[i]==grp[j]) grpmat[i,j]<-1

	rec<-cbind(k,k)
	
	for(j in 1:length(k))
	{
		r<-NULL
		diag(distmat)<-Inf
		for(i in 1:nrow(distmat))
		{
			sel<-which(distmat[i,]<=quantile(distmat[i,], k[j]/ncol(distmat)))
			this<-grpmat[i,sel]
			this<-1*(sum(this) >= length(this)*0.5)
			r<-c(r,this)
		}
		rec[j,2]<-sum(r==1)/length(r)
	}
	rec
}

dist.grp<-function(distmat, grp, k=1:20)
{
	
	grpmat<-matrix(0,nrow=length(grp),ncol=length(grp))
	for(i in 1:length(grp)) for(j in 1:length(grp)) if(grp[i]==grp[j]) grpmat[i,j]<-1

	rec<-cbind(k,k)
	
	for(j in 1:length(k))
	{
		r<-NULL
		diag(distmat)<-Inf
		for(i in 1:nrow(distmat))
		{
			sel<-which(distmat[i,]<=quantile(distmat[i,], k[j]/ncol(distmat)))
			r<-c(r,grpmat[i,sel])
		}
		rec[j,2]<-sum(r==1)/length(r)
	}
	rec
}
