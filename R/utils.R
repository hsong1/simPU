#'@export
soft_threshold<-function(x,th){
  if(th < 0){stop("threshold must be >=0")}
  pmax(abs(x)-th,0)*sign(x)
}

#'@export
hard_threshold<-function(x,cardinality=NULL,threshold=NULL){
  if(!xor(is.null(cardinality),is.null(threshold))){
    stop("only one of cardinality or threshold should be not NULL")
  }
  if(!is.null(cardinality)){
    if(cardinality < 0){stop("cardinality must be >=0")}
    k = floor(cardinality)
    if(k==0){
      x = rep(0,length(x))
    }else{
      ind = order(abs(x),decreasing = T)[1:k]
      x[-ind] = 0
    }
    return(x)
  }else if(!is.null(threshold)){
    th = threshold
    if(th < 0){stop("threshold must be >=0")}
    ind = which(abs(x)<th)
    x[ind] = 0
    return(x)
  }

}


#' proximal function
#' 
#' objective = f(w)+g(w);
#' g(w) = alpha*norm(w)_0(func = "l0"), alpha*norm(w)_1 (func="l1") or I_(norm(w)_0<=alpha)(w) (func="l0proj")
#' 
#'@export
prox<-function(x,alpha,func=c("l0","l1","l0proj")){
  func = match.arg(func,c("l0","l1","l0proj"))
  if(func=="l0"){
    px = hard_threshold(x,threshold = sqrt(2*alpha))
  }else if(func=="l1"){
    px = soft_threshold(x,alpha)
  }else if(func=="l0proj"){
    px = hard_threshold(x,cardinality = alpha)
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
