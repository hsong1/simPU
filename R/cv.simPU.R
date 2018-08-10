#'@import doParallel
#'@import Matrix
#'@export
#'
cv.simPU<-function(Xlu,zlu,step_size,decaying_factor=0.8,
                method = c("Kernel","PAVA","true"),eps = 1e-3,epoch=10,warmstart=T,
                penalty = c("l0","l1"),lambda = NULL, l2projection=FALSE,
                backtracking = FALSE, backtracking_tol=1e-3, 
                initial= NULL,pr1=NULL,betas=NULL,verbose=F){
  
  nl = sum(zlu); nu = length(zlu)-nl; n = nl+nu ; p = ncol(Xlu)
  
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
          .combine = list,
          .multicombine = TRUE)  %dopar%  {
            vlidx <- which(lfolds==i,arr.ind=TRUE)
            vuidx <- which(ufolds==i,arr.ind=TRUE)
            train_X=rbind(X_l[-vlidx,],rbind(X_u[-vuidx,]))
            train_z = c(rep(1,nl-length(vlidx)),rep(0,nu-length(vuidx)))
            
            res.i<-simPU(Xlu = train_X,zlu = train_z,step_size = step_size,
                    decaying_factor = decaying_factor,method = method,
                    epoch = epoch,warmstart = warmstart,penalty = penalty,
                    lambda = lambda,l2projection = l2projection,
                    backtracking = backtracking,backtracking_tol = backtracking_tol,
                    eps = eps,initial = initial,pr1 = pr1,betas = betas,verbose = verbose)

            return(res.i$beta)
  
          }

  stopCluster(cl)
    
    
  }
  
  
  
  
  

