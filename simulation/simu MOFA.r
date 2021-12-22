
library(MOFA2)
library(ggplot2)
library(mvtnorm)

library(BiocParallel)


nfac=3
n.simu=10
nx=c(1000) # dimension of X
ny=c(1000) # dimension of Y
kx=c(40,20,10) # number of effective variables in X
ky=c(40,20,10) # number of effective variables in Y per each effective variable in X
links=c(2,1,0) # types of link function, 0 is linear only, 1 is half linear half mix, 2 is all mixed
sn=c(3.16) # signal to noise
N=c(200,1000, 5000, 10000) # sample size
rho=c(0.6, 0.3, 0)

combos<-expand.grid(nx, ny, kx, ky, links, sn, N, rho)
colnames(combos)<-c("nx", "ny", 'kx', 'ky', 'links', 'sn', 'N',"rho")
to.delete<-which(combos[,3]*combos[,4]>combos[,2])
if(length(to.delete)>0) combos<-combos[-to.delete,]

sel<-which((combos[,3]==10 & combos[,4] == 10) | (combos[,3]==20 & combos[,4] == 20) |(combos[,3]==20 & combos[,4] == 40) |(combos[,3]==40 & combos[,4] == 20))
combos<-combos[sel,]

combos2<-new("list")
for(i in 1:nrow(combos)) combos2[[i]]<-combos[i,]

simu.fun<-function(vec,n.simu=10,nfac=3)
{

library(MOFA2)
library(mvtnorm)
	vec<-as.numeric(vec)
	
	nx<-vec[1]
	ny<-vec[2]
	kx<-vec[3]
	ky<-vec[4]
	links<-vec[5]
	sn<-vec[6]
	N<-vec[7]
	rho<-vec[8]
	
   filename<-paste("MOFA dim 3 nx", nx, 'ny', ny, 'kx', kx, 'ky', ky, 'links', links, 'sn', sn, 'N', N,"rho",rho, ".bin")
  is.there<-dir(pattern=filename)

  if(length(is.there) == 0)
  {
  rec<-new("list")
  for(mmm in 1:n.simu)
  {
  	sig<-matrix(rho, ncol=nx, nrow=nx)
  	diag(sig)<-1
	X<-rmvnorm(N, mean=rep(0,nx), sigma=sig)
	X<-pnorm(X)-0.5

 	sig<-matrix(rho, ncol=ny, nrow=ny)
  	diag(sig)<-1
	Y<-rmvnorm(N, mean=rep(0,ny), sigma=sig)
	Y<-pnorm(Y)-0.5

	z<-matrix(0, ncol=nfac, nrow=N)
	
	X0<-X

	
	for(k in 1:nfac)
	{
		X<-X0
		
		betas<-runif(kx, min=1, max=2)
		betas<-betas*sample(c(-1,1), kx, replace=T)
		for(i in 1:kx) X[,i]<-X[,i]*betas[i]

		z[,k]<-apply(X[,1:kx],1,sum)
		#z[,k]<-z[,k]-min(z[,k])
		#z[,k]<-z[,k]/max(z[,k])-0.5
	}

	X<-X0
	
	for(i in 1:(kx*ky))
	{
		
		betas<-runif(nfac, min=1, max=2)
		betas<-betas*sample(c(-1,1), nfac, replace=T)
		this.z<-z %*% betas
		
		this.z<-(this.z-mean(this.z))/sd(this.z)
		this.z<-this.z/3
		
		#this.z<-this.z-min(this.z)
		#this.z<-this.z/max(this.z)-0.5
		
		this.fun.type<-sample(2:5, 1)
		
		if(links == 0)
		{
			this.fun.type<-1
		}else if(links==1){
			if(runif(1)>0.5) this.fun.type<-1
		}

		if(this.fun.type == 2) Y[,i]<- abs(this.z)
        if(this.fun.type == 3) Y[,i]<- sin((5*(this.z+0.5)*pi))
		if(this.fun.type == 4) Y[,i]<- 1*(this.z>quantile(this.z, 0.25) & this.z<quantile(this.z, 0.75))
        #if(this.fun.type == 4) Y[,i]<- round((this.z-min(this.z))/(max(this.z)-min(this.z))*2) %% 2
        if(this.fun.type == 5) Y[,i]<- (2*this.z)^2
		if(this.fun.type == 1) Y[,i]<- this.z
	}
	

	
	#X<-X2
	
	for(i in 1:ncol(X)) X[,i]<- (X[,i]-mean(X[,i]))/sd(X[,i]) + rnorm(nrow(X), mean=0, sd=1/sn)
	for(i in 1:ncol(Y)) Y[,i]<- (Y[,i]-mean(Y[,i]))/sd(Y[,i]) + rnorm(nrow(Y), mean=0, sd=1/sn)

	############# all X are used as predictors

	data <- list(mir=t(X), gene=t(Y))

	MOFAobject <- create_mofa(data)

	model_opts <- get_default_model_options(MOFAobject)
	model_opts$num_factors<-3
	train_opts <- get_default_training_options(MOFAobject)
	data_opts <- get_default_data_options(MOFAobject)
	MOFAobject <- prepare_mofa(
	object = MOFAobject,
	data_options = data_opts,
	model_options = model_opts,
	training_options = train_opts
	)
	outfile = file.path(getwd(), paste("model",nx,ny,kx, ky,links, sn, N,rho,".hdf5"))
	MOFAobject.trained <- run_mofa(MOFAobject, outfile)
	weights <- get_weights(MOFAobject.trained, views = "all", factors = "all")[[1]]
	weights2<-weights3<-weights
	for(q in 1:ncol(weights2)) weights2[,q]<-rank(abs(weights2[,q]))
	for(q in 1:ncol(weights3)) weights3[,q]<-((weights3[,q])-mean(weights3[,q]))/sd(weights3[,q])
	mofa.imp<-apply(abs(weights),1,max)
	mofa.imp2<-apply(abs(weights2),1,max)
	mofa.imp3<-apply(abs(weights3),1,max)
	
	rec[[mmm]]<-cbind(mofa.imp,mofa.imp2,mofa.imp3)
    save(rec, file=filename)
  }


  }
}


z<-bplapply(combos2, simu.fun, n.simu=10)
