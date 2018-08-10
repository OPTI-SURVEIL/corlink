linkd <- function(d, initial_m=NULL, initial_u=NULL, p_init=0.5,fixed_col=NULL,alg="m"){
  colnames(d)[ncol(d)] <- "counts"
  L <- ncol(d) - 1
  if(is.null(initial_m)) initial_m <- rep(0.8,L)
  if(is.null(initial_u)) initial_u <- rep(0.2,L)
  if(is.null(fixed_col)) fixed_col <- c()
  themat <- create_01mat(L)
  out <- imputemissing(d,themat)
  out <- cbind(themat,counts=out$counts)
  out <- as.data.frame(out)
  colnames(out)[1:L] = colnames(d)[1:L]
  allterms <- apply(combn(colnames(d)[1:L],2),2,paste,collapse=':')
  
  if(alg == "i" | alg == "a"){
    results_independence <- EM_match_independence_v3(out,m=initial_m,u=initial_u,p_init=p_init,tol=10^-5, fixedcol=fixed_col)
    probs_independence <- reassign_probs(d[,1:L], out, results_independence$probs)
  }
  if(alg == "b" | alg == "a"){
    results_loglinear <- EM_match_modelsearch(out,m=initial_m,u=initial_u,p_init=p_init,tol=10^-5, fixedcol=fixed_col,allterms=allterms)
    probs_loglinear <- reassign_probs(d[,1:L], out,results_loglinear$probs)
  }
  if(alg == "m" | alg == "a"){
    results_loglinear_iu <- EM_match_modelsearch_iu(out,m=initial_m,u=initial_u,p_init=p_init,tol=10^-5, fixedcol=fixed_col,allterms=allterms)
    probs_loglinear_iu <- reassign_probs(d[,1:L], out,results_loglinear_iu$probs)
  }
  if(alg == "i") return(list(fitted_probs=cbind(d,"fitted_prob_match"=probs_independence),fitted_models=results_independence,
                             imputed_freqs = out))
  if(alg == "b") return(list(fitted_probs=cbind(d,"fitted_prob_match"=probs_loglinear),fitted_models=results_loglinear,
                             imputed_freqs = out))
  if(alg == "m") return(list(fitted_probs=cbind(d,"fitted_prob_match"=probs_loglinear_iu),fitted_models=results_loglinear_iu,
                             imputed_freqs = out))
  if(alg == "a") return(list(fitted_probs=cbind(d,"independence"=probs_independence,"log_linear"=probs_loglinear,"log_linear_iu"=probs_loglinear_iu),model_loglinear_iu=results_loglinear_iu,model_loglinear=results_loglinear,model_independence=results_independence,imputed_freqs = out))
  
}