#'@import np
#'@import isotone
#'@import Matrix
#'@export
#'
approx_mean<-function(x,z,method=c("Kernel","PAVA","true"),
                      warmstart=T,block=NULL,bandwidth=NULL,pr1=NULL){
  
  #Initialization
  if(is.null(block)){
    pava_start=NULL;
  }else{
    pava_start=1;
    gt.block = block
  }
  
  method = match.arg(method,c("Kernel","PAVA","true"))
  h<-function(eta,py1){
    log(nl/(py1*nu))+eta-log(1+exp(eta))
  }
  
  if(method=="true"){
    if(is.null(pr1)){stop("valid pr1 needs to be provided")}
    gt = 1/(1+exp(-h(x,pr1)))
    
  }else if(method=="PAVA"){
    if(warmstart){
      if(is.null(pava_start)){
        # save break
        pavafit = pava(y = z,x = x)
        gt = pavafit$fit
        gt.block = pavafit$block
        pava_start = 1
        
      }else{
        pavafit = pava(y = z,x = x,block = gt.block)
        gt = pavafit$fit;
        
        # gt.ref = gpava(x,z)$x
        # print(max(abs(gt-gt.ref)))
      }
    }else{
      gt = gpava(x,z)$x
    }

  }else{
    if(is.null(bandwidth)){
      # first use CV to choose one
      bw = npregbw(xdat=x,ydat = z,regtype="lc") #Nadaraya-Watson
      bandwidth = bw$bw
    }else{
      bw <- npregbw(xdat=x,ydat = z,regtype="lc",bandwidth.compute = FALSE,bws = bandwidth)
    }
    g <- npreg(bws = bw, gradients = TRUE)
    gt <- fitted(g)
  }
  
  if(method=="PAVA"&&warmstart){
    return(list(fit=gt,block=gt.block))
  }else if(method=="Kernel"){
    return(list(fit=gt,bandwidth=bandwidth))
  }else{
    return(list(fit=gt))
  }
  
}
