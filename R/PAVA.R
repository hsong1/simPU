#'@import Matrix
#'@useDynLib simPU
#'@export
pava <-function(y,x=NULL,block=NULL){
  
  n = length(y)
  if(!is.null(x)){
    pm = order(x); ipm = invPerm(pm)
    y = y[pm]}
  
  if(!is.null(block)){
    nblocks = block$nblocks
    sidx = block$sidx
    size = block$size
    # Error check
    if(nblocks<0||nblocks>n){stop("nblocks must be in 1:n")}
    if(!all(sidx%in%(1:n))){stop("block starting index must be in 1:n")}
    if(!all(sidx==c(sidx[1],sidx[1]+cumsum(size)[-nblocks]))){stop("sidx and block size are not compatible")}
    
    sidx = sidx-1 # to make it comparable with c++ indexing
    fit = pava_main(y,nblocks,sidx,size,TRUE)
  }else{
    fit = pava_main(y,1,1,1,FALSE)
  }
  fitted = fit$fit
  if(!is.null(x)){fitted = fitted[ipm]}
  block=list("nblocks"=fit$nblocks,"sidx"=fit$sidx+1,"size"=fit$size)
  return(list(fit=fitted,block=block))
}
