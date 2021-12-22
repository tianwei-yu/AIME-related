
fdrgamma<-function(x, use.percentile=NA, min.use.percentile=0.8)
{
	library(MASS)
	
	x0<-x
	if(length(x)>50000) x<-sample(x, 50000, replace=F)
	
	all.use.percentile<-use.percentile
	if(is.na(use.percentile)) all.use.percentile<-seq(min.use.percentile,0.99,by=0.01)
	
	d.x<-density(x,from=0)
	x.mode<-d.x$x[which(d.x$y==max(d.x$y))[1]]
	rec<-rep(0, length(all.use.percentile))
	
	for(i in 1:length(all.use.percentile))
	{
		use.percentile<-all.use.percentile[i]
		params<-fitdistr(x[x<quantile(x,use.percentile)],"gamma")
	
		d.x.lim<-density(x,from=0,to=quantile(x,use.percentile))
		d.null.lim<-dgamma(d.x.lim$x, shape=params$estimate[1], rate=params$estimate[2])*use.percentile
		sel<-which(d.x.lim$x <= x.mode)
		rec[i]<-sum((d.x.lim$y[sel]-d.null.lim[sel])^2)
	
	}
	
	sel<-which(rec==min(rec))[1]
	use.percentile=all.use.percentile[sel]
	cat('null %', use.percentile, "; ")

	d.x<-density(x,from=0, n=min(512, length(x)/2))
	params<-fitdistr(x[x<quantile(x,use.percentile)],"gamma")
	d.null<-dgamma(d.x$x, shape=params$estimate[1], rate=params$estimate[2])
	
	null.mode.x<-d.x$x[which(d.null==max(d.null))[1]]

	#plot(d.x, type="l")
	hist(x,nclass=min(round(length(x)/10),100),freq=FALSE)
	lines(d.x$x, d.null*use.percentile,col="red")
	lines(d.x, col="blue")

	lfdr<-d.null/d.x$y*use.percentile
	lfdr[lfdr>1]<-1
	lfdr[is.na(lfdr)]<-1
	lfdr[d.x$x<=null.mode.x]<-1
	sel<-which(lfdr==max(lfdr))[1]
	lfdr[1:sel]<-max(lfdr)
	lfdr<-cummin(lfdr)
	
	all.lfdr<-approx(d.x$x, lfdr, x0, rule=2)
	
	return(all.lfdr[[2]])
}
