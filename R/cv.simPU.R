#'@import doParallel
#'@import Matrix
#'@export
#'
cv.simPU<-function(Xlu,zlu,step_size,
                method = c("PAVA","Kernel","true"),eps = 1e-3,epoch=10,warmstart=T,
                penalty = c("l0","l1"),lambdaseq = NULL, l2projection=TRUE,
                initial= NULL,pr1=NULL,betas=NULL,verbose=F, nfolds = 10, nCores=1){
  
  nl = sum(zlu); nu = length(zlu)-nl; n = nl+nu ; p = ncol(Xlu)
  method = match.arg(method,c("PAVA","Kernel","true"))
  penalty = match.arg(penalty, c("l0","l1"))
  k = nfolds; 
  
  if(is.null(lambdaseq)){
    stop("currently default lambda seq not implemented")
  }
  if(verbose){cat('Fitting full data\n')}
  if(is.null(initial)){initial = rep(0,p)}else{beta=initial}
  betaMat<-matrix(NA,nrow=p,ncol=length(lambdaseq))
  for(j in 1:length(lambdaseq)){
    res.full<-simPU(Xlu = Xlu,zlu = zlu,step_size = step_size,method = method,
                    epoch = epoch,warmstart = warmstart,penalty = penalty,
                    lambda = lambdaseq[j],l2projection = l2projection,
                    eps = eps,initial = initial,pr1 = pr1,betas = betas,verbose = verbose)
    betaMat[,j]<-res.full$beta
    initial = res.full$beta # warm start
  }
  colnames(betaMat) = paste("l",1:length(lambdaseq),sep="")
  if(is.null(colnames(Xlu))){
    rownames(betaMat) = paste("X",1:p,sep="")
  }else{
    rownames(betaMat) = colnames(Xlu)
  }
  
  
  # shuffle Xlu
  pl <- sample(1:nl); pu <- sample(1:nu)
  X_l <- Xlu[1:nl,];X_l <- X_l[pl,]
  X_u <- Xlu[(nl+1):(nl+nu),];X_u <- X_u[pu,]
  
  # k=10
  # Create k equally size folds
  lfolds <- cut(seq(1,nl),breaks=k,labels=FALSE)
  ufolds <- cut(seq(1,nl),breaks=k,labels=FALSE)
  
  if(nCores==0){nCores = max(1, detectCores() - 1)
  }else{
    nCores = min(nCores,detectCores())
  }
  
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())    
  if(verbose){cat('Cross-Validation with',getDoParWorkers(), 'workers\n')}
  
  res<-foreach(i=1:k,
          .packages = "simPU",
          .combine = cbind,
          .multicombine = TRUE)  %dopar%  {
            
            vlidx <- which(lfolds==i,arr.ind=TRUE)
            vuidx <- which(ufolds==i,arr.ind=TRUE)
            train_X=rbind(X_l[-vlidx,],rbind(X_u[-vuidx,]))
            train_z = c(rep(1,nl-length(vlidx)),rep(0,nu-length(vuidx)))
            
            # betaMat<-matrix(NA,nrow=p,ncol=length(lambdaseq))
            sqloss <-c()
            for(j in 1:length(lambdaseq)){
              res.i<-simPU(Xlu = train_X,zlu = train_z,step_size = step_size,method = method,
                           epoch = epoch,warmstart = warmstart,penalty = penalty,
                           lambda = lambdaseq[j],l2projection = l2projection,
                           eps = eps,initial = initial,pr1 = pr1,betas = betas,verbose = verbose)
              # betaMat[,j]<-res.i$beta
              sqloss[j]<-res.i$sqloss
              initial = res.i$beta
            }
            return(sqloss)
  
          }

  cvm = apply(res,1,mean); cvsd = apply(res,1,sd)
  indmin <- min(which(cvm==min(cvm)))
  lambda.min <- lambdaseq[indmin]
  
  ind <-  intersect(which(cvm>=cvm[indmin]+cvsd[indmin]),(1:indmin))
  if(length(ind)==0){ind1se <-  indmin
  } else {
    ind1se <-  max(ind)
  }
  lambda.1se <- lambdaseq[ind1se]
  
  
  stopCluster(cl)
  return(list(lambdaseq=lambdaseq,beta=betaMat,cvm=cvm,cvsd=cvsd,lambda.min=lambda.min,lambda.1se=lambda.1se))  
}
