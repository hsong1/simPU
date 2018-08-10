#'@export
soft_threshold<-function(x,th){
  if(th < 0){stop("threshold must be >=0")}
  pmax(abs(x)-th,0)*sign(x)
}

#'@export
hard_threshold<-function(x,th){
  if(th < 0){stop("threshold must be >=0")}
  ind = which(abs(x)<th)
  x[ind] = 0
  x
}


#' proximal function
#' objective = f(w)+alpha*g(w), g=l0 or l1
#'@export
prox<-function(x,alpha,norm=c("l0","l1")){
  norm = match.arg(norm,c("l0","l1"))
  if(norm=="l0"){
    px= hard_threshold(x,alpha)
  }else{
    px=soft_threshold(x,sqrt(2*alpha))
  }
  return(px)
}

#'@export
proj<-function(x,alpha,norm=c("l2")){
  if(all(x==0)){
    return(x)
  }else{
    return(x/sqrt(as.numeric(crossprod(x))))
  }
}
