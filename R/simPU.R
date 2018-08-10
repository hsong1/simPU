#'@import Matrix
#'@export
#'
simPU<-function(Xlu,zlu,step_size,
                method = c("PAVA","Kernel","true"),eps = 1e-4,epoch=10,warmstart=T,
                penalty = c("l0","l1","l0proj"),lambda = 0,cardinality =ncol(Xlu),l2projection=TRUE,
                initial= NULL,pr1=NULL,betas=NULL,verbose=F){
  
  method = match.arg(method,c("PAVA","Kernel","true"))
  penalty = match.arg(penalty, c("l0","l1","l0proj"))
  eta = step_size
  if(penalty=="l0proj"&&cardinality<ncol(Xlu)){lambda=1}
  
  l2err = c()
  nl = sum(zlu); nu = length(zlu)-nl; n = nl+nu ; p = ncol(Xlu)
  
  #If no initialization available, do a random initialization
  if(is.null(initial)){beta= rep(0,p)}else{beta=initial}
  
  # Start iteration
  ff = c()
  # ff_upper = c(); ff_upper[1]=0
  # use lazy evaluation
  eval_fn_g <- function(beta,block=NULL,bandwidth=NULL){
    x= Xlu%*%beta
    # estimate mean using (x,zlu)
    gg = approx_mean(x = x,z = zlu,method = method,
                     warmstart = warmstart,block = block,
                     bandwidth = bandwidth,
                     pr1 = pr1)
    gt = gg$fit; 
    if(method=="PAVA"&&warmstart){block = gg$block}
    if(method=="Kernel"){bandwidth = gg$bandwidth}
    # descending direction = -vt
    mean_diff = gt-zlu # expected - observed
    vt = (t(Xlu)%*%mean_diff)/n # upper bound of gradient
    # estimated loss
    ft = crossprod(mean_diff)/(2*n); ff[1] = ft; #f(beta)
    return(list(ft=as.numeric(ft),vt=vt,block=block,bandwidth=bandwidth))
  }
  
  # when t=1
  fn_gt = eval_fn_g(beta); ft = fn_gt$ft; vt = fn_gt$vt; ff[1] = ft;
  
  if(method=="PAVA"&&warmstart){
    gt.block = fn_gt$block
  }else{
    gt.block = NULL
  }
  if(method=="Kernel"){
    bandwidth = fn_gt$bandwidth
  }else{
    bandwidth = NULL
  }
  
  if(!is.null(betas)){
    l2err[1]= sum((proj(beta,1,"l2")-betas)^2)
  }
  
  converged = FALSE; t=2; 
  betaMat = beta
  while(t<=(n*epoch) && (!converged)){
    
    if(t/n==floor(t/n)){
      if(verbose){cat("epoch = ",t/n,"\n")}
      if(method=="Kernel"){
        bandwidth = NULL # update Kernel bandwidth
      }
    }
    # update beta
    if(lambda>0){
      # if penalty="l0proj", project a gradient step onto l0 <= lambda ball
      # otherwise do prox_{\eta\lambda*penalty}(gradient step)
      if(penalty =="l0proj"){
        alpha = cardinality
      }else{
        alpha = eta*lambda
      }
      # take a proximal step
      beta1 = prox(beta - eta*vt,alpha = alpha,func = penalty)
      if(l2projection){beta1 = proj(beta1,1,"l2")}
    }else{
      beta1 = beta - eta*vt
      if(l2projection){beta1 = proj(beta1,1,"l2")}
    }
    
    # difference between beta_(t+1)-beta_t
    beta_diff = max(abs(beta1-beta)); 
    converged = beta_diff<eps
    
    if(!converged){
      fn_gt = eval_fn_g(beta1,block = gt.block,bandwidth = bandwidth); ft = fn_gt$ft; 
      vt = fn_gt$vt; ff[t] = ft;
      if(method=="PAVA"&&warmstart){gt.block = fn_gt$block}
      if(method=="Kernel"&&is.null(bandwidth)){
        if(verbose){cat("update bandwidth\n")}
        bandwidth = fn_gt$bandwidth}
      
      beta = beta1; betaMat = cbind(betaMat,beta)
      
      # if true beta is provided
      if(!is.null(betas)){
        l2err[t]= sum((proj(beta,1,"l2")-betas)^2)
      }
      t=t+1
    }
  } # end of while loop
  
  if(verbose&&converged){
    cat("converged at iter =",t,"\n")
  }
  
  # Finally scale to unit vector
  if(!l2projection){
    beta = proj(beta,1,"l2")
  }
  # Return
  
  if(!is.null(betas)){
    return(structure(list(beta=beta,sqloss=ff,l2err=l2err,betaMat=betaMat,convergence=converged,call=match.call()),class="simPU"))
  }else{
    return(structure(list(beta=beta,sqloss=ff,betaMat=betaMat,convergence=converged,call=match.call()),class="simPU"))
  }
}
