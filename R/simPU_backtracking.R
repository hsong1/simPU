#' #'@import Matrix
#' #'
#' simPU_backtracking<-function(Xlu,zlu,step_size,
#'                 method = c("Kernel","PAVA","true"),eps = 1e-3,epoch=10,warmstart=T,
#'                 penalty = c("l0","l1"),lambda = 0, l2projection=FALSE,
#'                 decaying_factor=0.8, step_size_min = 1e-3,
#'                 initial= NULL,pr1=NULL,betas=NULL,verbose=F){
#'   
#'   method = match.arg(method,c("Kernel","PAVA","true"))
#'   penalty = match.arg(penalty, c("l0","l1"))
#'   eta = step_size 
#'   
#'   l2err = c()
#'   nl = sum(zlu); nu = length(zlu)-nl; n = nl+nu ; p = ncol(Xlu)
#'   
#'   #If no initialization available, do a random initialization
#'   if(is.null(initial)){beta = rnorm(p)}else{beta=initial}
#'   
#'   # Start iteration
#'   ff = c()
#'   ff_upper = c(); ff_upper[1]=0
#'   # use lazy evaluation
#'   eval_fn_g <- function(beta,block=NULL){
#'     x= Xlu%*%beta
#'     # estimate mean using (x,zlu)
#'     gg = approx_mean(x = x,z = zlu,method = method,
#'                      warmstart = warmstart,block = block,pr1 = pr1)
#'     gt = gg$fit; 
#'     if(method=="PAVA"&&warmstart){block = gg$block}
#'     
#'     # descending direction = -vt
#'     mean_diff = gt-zlu # expected - observed
#'     vt = (t(Xlu)%*%mean_diff)/n # upper bound of gradient
#'     # estimated loss
#'     ft = crossprod(mean_diff)/(2*n); ff[1] = ft; #f(beta)
#'     return(list(ft=as.numeric(ft),vt=vt,block=block))
#'   }
#'   
#'   # when t=1
#'   fn_gt = eval_fn_g(beta); ft = fn_gt$ft; vt = fn_gt$vt; ff[1] = ft;
#'   if(method=="PAVA"&&warmstart){
#'     gt.block = fn_gt$block
#'   }else{
#'     gt.block =NULL
#'   }
#'   
#'   if(!is.null(betas)){
#'     l2err[1]= sum((proj(beta,1,"l2")-betas)^2)
#'   }
#'   
#'   converged = FALSE; t=2; backtr_tol_factor = 1
#'   if(is.null(backtracking_tol)){
#'     set_backtr_tol = TRUE; 
#'     backtracking_tol = 100 # will be reset after 1 epoch
#'   }else{
#'     set_backtr_tol = NULL
#'   }
#'   
#'   while(t<=(n*epoch) && (!converged)){
#'     success_count = 0;
#'     if(verbose){if(t/n==floor(t/n)){cat("epoch = ",t/n,"\n")}}
#'     # update beta
#'     if(lambda>0){
#'       beta1 = prox(beta - eta*vt,alpha = eta*lambda,norm = penalty)
#'       if(l2projection){beta1 = proj(beta,1,"l2")}
#'     }else{
#'       beta1 = beta - eta*vt
#'       if(l2projection){beta1 = proj(beta,1,"l2")}
#'     }
#'     
#'     beta_diff = max(abs(beta1-beta)); 
#'     converged = beta_diff<eps
#'     
#'     if(!converged){
#'       fn_gt = eval_fn_g(beta1,block = gt.block); ft = fn_gt$ft; 
#'       vt = fn_gt$vt; ff[t] = ft;
#'       if(method=="PAVA"&&warmstart){gt.block = fn_gt$block}
#'       beta = beta1
#'       t=t+1
#'     }
#'     if(!backtracking){
#'       
#'       
#'     }else{
#'       descending= FALSE
#'       if(!is.null(set_backtr_tol)){
#'         if(t/n ==floor(t/n)){
#'           backtracking_tol = backtr_tol_factor*sd(ff_upper);
#'           cat("backtracking_tol",backtracking_tol,"at:",t/n,"\n")
#'         }}
#'       while(!descending){
#'         # update beta
#'         if(lambda>0){
#'           beta1 = prox(beta - eta*vt,alpha = eta*lambda,norm = penalty)
#'           if(l2projection){beta1 = proj(beta,1,"l2")}
#'         }else{
#'           beta1 = beta - eta*vt
#'           if(l2projection){beta1 = proj(beta,1,"l2")}
#'         }
#'         
#'         fn_gt = eval_fn_g(beta1,block = gt.block); ft = fn_gt$ft; 
#'         f_upper = ff[t-1]+crossprod(vt,beta1-beta)+(1/(2*eta))*crossprod(beta1-beta)
#'         ff_upper[t] = f_upper
#'         
#'         if(ft<=f_upper+backtracking_tol){
#'           descending = TRUE; success_count = success_count+1
#'           beta_diff = max(abs(beta1-beta)); 
#'           converged = beta_diff<eps;
#'           if(!converged){
#'             vt = fn_gt$vt; ff[t] = ft;
#'             if(method=="PAVA"&&warmstart){gt.block = fn_gt$block}
#'             beta = beta1
#'             t=t+1
#'             if(success_count>round(n/10)){
#'               success_count = 0; eta = eta/decaying_factor
#'               if(verbose){cat("step size(*):",round(eta,4),"->",round(eta/decaying_factor,4),"\n")}
#'             }
#'           }
#'         }else{
#'           #if eta gets too small, consider allowing larger gap
#'           if(decaying_factor*eta<1e-3){
#'             backtr_tol_factor = 2*backtr_tol_factor;
#'             if(verbose){cat("retry with a larger margin; scale factor=",backtr_tol_factor,"\n")}
#'             backtracking_tol = backtr_tol_factor*sd(ff_upper);
#'           }else{
#'             if(verbose){cat("step size:",round(eta/decaying_factor,4),"->",round(eta,4),"\n")}
#'             eta = decaying_factor*eta  
#'           } 
#'         }
#'       }
#'     }
#'     # if true beta is provided
#'     if(!is.null(betas)){
#'       l2err[t]= sum((proj(beta,1,"l2")-betas)^2)
#'     }
#'   }
#'   
#'   if(verbose&&converged){
#'     cat("converged at iter =",t,"\n")
#'   }
#'   
#'   #Finally scale to unit vector
#'   if(!l2projection){
#'     beta = beta/sqrt(as.numeric(crossprod(beta)))
#'   }
#'   
#'   if(!is.null(betas)){
#'     return(list(beta=beta,l2err=l2err,diff=diff,fval=ff,convergence=converged))
#'   }else{
#'     return(beta,fval=ff,convergence=converged)
#'   }
#' }
