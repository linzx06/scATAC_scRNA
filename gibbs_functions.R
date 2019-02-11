getClusterGibbs <- function(data_acc, data_exp, overlap_seq_acc, overlap_seq_exp, 
                            mu0=0.5, mu1=0.5, nu0=2, nCluster,
                            cluster_acc_Ini=NULL, cluster_exp_Ini=NULL, 
                            u_acc_tilda_Ini=NULL, v_exp_tilda_Ini=NULL, 
                            omega_acc_Ini=NULL, omega_exp_Ini=NULL, nu1_Ini=1,
                            q_acc_Ini=NULL, q_exp_Ini=NULL, pi_exp_Ini=NULL,
                            niter, quiet=F, eval_logPost=F){
  ## estimate the mixture density
  if (!quiet){
    cat("\n", "Estimating mixture distribution", "\n")
  }
  f1_acc <- (data_acc!=0) + 0
  f0_acc <- 1 - f1_acc
  data_exp_non0 <- data_exp[which(data_exp!=0)]
  if (length(data_exp_non0)>50000){
    data_exp_non0 <- sample(data_exp_non0, 50000)
  }
  data_exp_non0 <- log(data_exp_non0 + 1)
  ## fit gamma-gamma mixture for gene expression data
  fit <- gammamixEM(x=data_exp_non0, k = 2, epsilon=2, verb=F)
  data_exp <- log(data_exp+1)
  if (!quiet){
    cat("\n", "The mixture components for scRNA-Seq", "\n")
    print(fit$gamma.pars)
  }
  dens_scRNA_gg <- cal_density_gg_thres(gg_fit=fit, data=data_exp, thres_low=-1, thres_high=20)
  f1_exp <- dens_scRNA_gg$f1
  f0_exp <- dens_scRNA_gg$f0
  
  ## initialization
  # z_acc
  if (is.null(cluster_acc_Ini)){
    cluster_acc_Ini <- sample(nCluster, nrow(f1_acc), replace=T)
    z_acc <- expandCluster(cluster_acc_Ini)  
  } else {
    z_acc <- expandCluster(cluster_acc_Ini)  
  }
  
  # z_exp
  if (is.null(cluster_exp_Ini)){
    cluster_exp_Ini <- sample(nCluster, nrow(f1_exp), replace=T)
    z_exp <- expandCluster(cluster_exp_Ini)  
  } else {
    z_exp <- expandCluster(cluster_exp_Ini)  
  }
  
  # u_acc_tilda
  if (is.null(u_acc_tilda_Ini)){
    u_acc_tilda <- get_v_Ini(f1_acc, f0_acc, 0.5)
  } else {
    u_acc_tilda <- u_acc_tilda_Ini
  }
  
  # v_exp_tilda
  if (is.null(v_exp_tilda_Ini)){
    v_exp_tilda <- get_v_Ini(f1_exp, f0_exp, 0.5)
  } else {
    v_exp_tilda <- v_exp_tilda_Ini
  }
  
  # omega_acc
  if (is.null(omega_acc_Ini)){
    omega_acc <- get_omega_Ini(u_acc_tilda, cluster_acc_Ini)
    omega_acc <- omega_acc/max(omega_acc)
    omega_acc[which(omega_acc<=0.05)] <- 0.05
    omega_acc[which(omega_acc>=0.95)] <- 0.95  
  } else {
    omega_acc <- omega_acc_Ini
  }
  
  # omega_exp
  if (is.null(omega_exp_Ini)){
    omega_exp <- get_omega_Ini(v_exp_tilda, cluster_exp_Ini)
    omega_exp <- omega_exp/max(omega_exp)
    omega_exp[which(omega_exp<=0.05)] <- 0.05
    omega_exp[which(omega_exp>=0.95)] <- 0.95  
  } else {
    omega_exp <- omega_exp_Ini
  }
  
  # q_acc
  if (!is.null(q_acc_Ini)){
    q_acc <- q_acc_Ini
  } else {
    q_acc <- rep(0.5, nrow(f1_acc))
  }
  
  # q_exp
  if (!is.null(q_exp_Ini)){
    q_exp <- q_exp_Ini
  } else {
    q_exp <- rep(0.5, nrow(f1_exp))
  }

  # pi_exp
  if (!is.null(pi_exp_Ini)){
    pi_exp <- pi_exp_Ini
  } else {
    pi_exp <- matrix(0.3, nrow=2, ncol=nrow(f1_exp))
    pi_exp[1,] <- 0.8
  }
  
  # nu1
  nu1 <- nu1_Ini
    
  # compuate eta_acc
  eta_acc <- calEta_acc(omega_acc, q_acc)
  # compuate EtaLambda_exp
  EtaLambda_exp <- calEtaLambda_exp(omega_exp, q_exp, pi_exp, overlap_seq_exp)
  
  ## traces
  z_acc_Trace <- array(dim=c(niter, dim(z_acc)))
  z_exp_Trace <- array(dim=c(niter, dim(z_exp)))
  q_acc_Trace <- matrix(nrow=niter, ncol=length(q_acc))
  q_exp_Trace <- matrix(nrow=niter, ncol=length(q_exp))                          
  pi_exp_Trace <- array(dim=c(niter, dim(pi_exp)))                      
  omega_acc_Trace <- array(dim=c(niter, dim(omega_acc)))
  omega_exp_Trace <- array(dim=c(niter, dim(omega_exp)))
  nu1_Trace <- rep(NA, niter)

  # the auxiliary variable
  h <- 1:nCluster
  h_Trace <- c()
  logpost_Trace <- c()
  cat("\n", "MCMC:", "\n")
  for (iter in 1:niter){
    if (!quiet & iter%%round(niter/10)==0){
      cat("\r", round(iter/niter*100), "%", "completed", "\n")
    }
    ## update u_acc_tilda
    u_acc_tilda <- update_v(f1_acc, f0_acc, z_acc, eta_acc)
    
    ## update v_exp_tilda
    v_exp_tilda <- update_v(f1_exp, f0_exp, z_exp, EtaLambda_exp)
    
    ## update z_acc
    z_acc <- update_x(u_acc_tilda, eta_acc)  
  
    ## update z_exp
    z_exp <- update_x(v_exp_tilda, EtaLambda_exp)  
    
    ## update omega_acc
    omega_acc <- update_omega_acc(omega_acc, z_acc, u_acc_tilda, q_acc, omega_exp, 
                                  nu1, mu0, nu0, overlap_seq_acc, overlap_seq_exp)
    
    ## update omega_exp
    omega_exp <- update_omega_exp(omega_exp, z_exp, v_exp_tilda, q_exp, pi_exp, 
                                  omega_acc, nu1, nu0, mu1, overlap_seq_acc, overlap_seq_exp)

    ## update q_acc
    q_acc <- update_q_acc(q_acc, u_acc_tilda, z_acc, omega_acc)
    
    ## update q_exp
    q_exp <- update_q_exp(q_exp, pi_exp, v_exp_tilda, z_exp, omega_exp, overlap_seq_exp)
    
    ## update pi_exp
    pi_exp <- update_pi_exp(q_exp, pi_exp, v_exp_tilda[,overlap_seq_exp], z_exp, omega_exp[,overlap_seq_exp])
    
    ## update nu1
    nu1 <- update_nu1(nu1, as.numeric(omega_acc[,overlap_seq_acc]), 
                      as.numeric(omega_exp[,overlap_seq_exp]), min=0, max=50, iter)
    nu1_Trace[iter] <- nu1
    
    ## update h
    h <- update_h(h, omega_acc[,overlap_seq_acc], omega_exp[,overlap_seq_exp], nu1) 
    h_Trace <- rbind(h_Trace, h)
    omega_exp <- omega_exp[h,]
    z_exp <- z_exp[,h]
    
    ## compuate lambda
    eta_acc <- calEta_acc(omega_acc, q_acc)
    EtaLambda_exp <- calEtaLambda_exp(omega_exp, q_exp, pi_exp, overlap_seq_exp)
    
    ## calculate log posterior
    if (eval_logPost){
      logpost_Trace <- c( logpost_Trace, cal_logpost(log(f1_acc), log(f0_acc), log(f1_exp), log(f0_exp), 
                                                     omega_acc, omega_exp, pi_exp,  
                                                     eta_acc, EtaLambda_exp,
                                                     u_acc_tilda, z_acc, v_exp_tilda, z_exp,
                                                     mu0, mu1, nu0, nu1, 
                                                     atac_bin=F, rna_bin=F, overlap_seq_acc, overlap_seq_exp) )  
    }
    
    ## track
    z_acc_Trace[iter,,] <- z_acc
    z_exp_Trace[iter,,] <- z_exp
    q_acc_Trace[iter,] <- q_acc
    q_exp_Trace[iter,] <- q_exp 
    pi_exp_Trace[iter,,] <- pi_exp
    omega_acc_Trace[iter,,] <- omega_acc
    omega_exp_Trace[iter,,] <- omega_exp
  }
  ## correct label switching
  lab_perm <- ecr.iterative.1(z=apply(z_acc_Trace, c(1,2), which.max), K=nCluster, threshold=1e-20)
  omega_acc_Trace <- permute.mcmc(omega_acc_Trace, lab_perm$permutations)[[1]]
  omega_exp_Trace <- permute.mcmc(omega_exp_Trace, lab_perm$permutations)[[1]]
  z_acc_Trace <- permute.mcmc(aperm(z_acc_Trace, c(1, 3, 2)), lab_perm$permutations)[[1]]
  z_exp_Trace <- permute.mcmc(aperm(z_exp_Trace, c(1, 3, 2)), lab_perm$permutations)[[1]]
  z_acc_Trace <- aperm(z_acc_Trace, c(1, 3, 2))
  z_exp_Trace <- aperm(z_exp_Trace, c(1, 3, 2))
  
  ## cluster probability and cluster assignment
  z_acc_prob <- apply(z_acc_Trace[(niter/2+1):niter,,], c(2, 3), mean)
  z_exp_prob <- apply(z_exp_Trace[(niter/2+1):niter,,], c(2, 3), mean)
  cluster_acc <- apply(z_acc_prob, 1, which.max)
  cluster_exp <-  apply(z_exp_prob, 1, which.max)
  
  return(list(z_acc_prob=z_acc_prob, z_exp_prob=z_exp_prob, cluster_acc=cluster_acc, cluster_exp=cluster_exp, gamma_mixture_exp=fit,
              z_acc_Trace=z_acc_Trace, q_acc_Trace=q_acc_Trace, omega_acc_Trace=omega_acc_Trace, 
              z_exp_Trace=z_exp_Trace, q_exp_Trace=q_exp_Trace, omega_exp_Trace=omega_exp_Trace, pi_exp_Trace=pi_exp_Trace,
              h_Trace=h_Trace, logpost_Trace=logpost_Trace, nu1_Trace=nu1_Trace))
}

calEta_acc <- function(omega_acc, q_acc){
  eta_acc <- outer(omega_acc, q_acc) 
  return(eta_acc)
}

calEtaLambda_exp <- function(omega_exp, q_exp, pi_exp, overlap_seq_exp){
  EtaLambda_exp <- array( dim=c(dim(omega_exp), length(q_exp)) )
  EtaLambda_exp[,overlap_seq_exp,] <- outer(omega_exp[,overlap_seq_exp], q_exp)
  EtaLambda_exp[,-overlap_seq_exp,] <- outer(omega_exp[,-overlap_seq_exp], q_exp*pi_exp[1,]) + 
    outer(1 - omega_exp[,-overlap_seq_exp], q_exp*pi_exp[2,])
  return(EtaLambda_exp)
}

expandCluster <- function(cluster){
  nCluster <- max(cluster)
  X <- matrix(0, nrow=length(cluster), ncol=nCluster)
  for (i in 1:nCluster){
    X[which(cluster==i), i] <- 1  
  }
  return(X)
}

get_v_Ini <- function(f1, f0, prob){
  v <- ((prob*f1) > ((1-prob)*f0)) + 0
  return(v)
}

get_omega_Ini <- function(uDnase, clusterDnase){
  nComp <- max(clusterDnase)
  n <- length(clusterDnase)
  pi_pK <- matrix(nrow=nComp, ncol=ncol(uDnase))
  for (i in 1:nComp){
    if (length(which(clusterDnase==i))==1){
      pi_pK[i,] <- uDnase[which(clusterDnase==i),]
    } else {
      pi_pK[i,] <- colSums(uDnase[which(clusterDnase==i),])/length(which(clusterDnase==i))  
    }
  }
  return(pi_pK)
}

update_v <- function(f1, f0, xRna, lambda){
  #logp1 <- logf1 + apply(log(lambda), 3, function(lambdaJ){colSums(lambdaJ*t(xRna))})
  #logp0 <- logf0 + apply(log(1-lambda), 3, function(lambdaJ){colSums(lambdaJ*t(xRna))})
  t0m1 <- t(colSums(sweep(log(1-lambda), c(1, 3), t(xRna), "*"), 1)) - t(colSums(sweep(log(lambda), c(1, 3), t(xRna), "*"), 1))
  prob1 <- f1/( f1 + f0*exp(t0m1) )
  vRna <- (prob1>runif(length(prob1))) + 0
  return(vRna)   
}

update_x <- function(vRna, lambda){
  vRnat <- t(vRna)
  logp <- apply(log(lambda), 1, function(lambdaK){colSums(lambdaK*vRnat)}) + apply(log(1-lambda), 1, function(lambdaK){colSums(lambdaK*(1-vRnat))})
  logp <- logp - apply(logp, 1, max)
  probs <- exp(logp)
  xRna <- t(apply(probs, 1, rmultinom, n=1, size=1))
  return(xRna)   
}

update_omega_acc <- function(pi_pK_atac, x_atac, v_atac, q_acc, pi_pK_rna, 
                             theta, pi_p, phi, overlap_seq_atac, overlap_seq_rna){
  #
  pi_pK_atac_New <- matrix(runif(length(pi_pK_atac)), nrow=nrow(pi_pK_atac), ncol=ncol(pi_pK_atac))
  #
  p_o2n <- matrix(1, nrow=nrow(pi_pK_atac), ncol=ncol(pi_pK_atac))
  p_n2o <- matrix(1, nrow=nrow(pi_pK_atac), ncol=ncol(pi_pK_atac))
  # 
  x_atacE <- array(rep(x_atac, ncol(pi_pK_atac)), dim=c(dim(x_atac), ncol(pi_pK_atac)))
  v_atacE <- array(rep(v_atac, nrow(pi_pK_atac)), dim=c(dim(v_atac), nrow(pi_pK_atac)))
  x_atacE <- aperm(x_atacE, c(2,3,1))
  v_atacE <- aperm(v_atacE, c(3,2,1))
  #
  lambda_atac <- calEta_acc(pi_pK_atac, q_acc)
  logPold <- rowSums(log(lambda_atac)*v_atacE*x_atacE + log(1-lambda_atac)*(1-v_atacE)*x_atacE, dims=2) + 
    (phi*pi_p-1)*log(pi_pK_atac) + (-phi*pi_p+phi-1)*log(1-pi_pK_atac) 
  logPold[,overlap_seq_atac] <- logPold[,overlap_seq_atac] + (theta*pi_pK_atac[,overlap_seq_atac]-1)*log(pi_pK_rna[,overlap_seq_rna]) + 
    (-theta*pi_pK_atac[,overlap_seq_atac]+theta-1)*log(1-pi_pK_rna[,overlap_seq_rna]) -
    lbeta(theta*pi_pK_atac[,overlap_seq_atac], -theta*pi_pK_atac[,overlap_seq_atac]+theta)
  #
  lambda_atac <- calEta_acc(pi_pK_atac_New, q_acc)
  logPnew <- rowSums(log(lambda_atac)*v_atacE*x_atacE + log(1-lambda_atac)*(1-v_atacE)*x_atacE, dims=2) + 
    (phi*pi_p-1)*log(pi_pK_atac_New) + (-phi*pi_p+phi-1)*log(1-pi_pK_atac_New) 
  logPnew[,overlap_seq_atac] <- logPnew[,overlap_seq_atac] + (theta*pi_pK_atac_New[,overlap_seq_atac]-1)*log(pi_pK_rna[,overlap_seq_rna]) + 
    (-theta*pi_pK_atac_New[,overlap_seq_atac]+theta-1)*log(1-pi_pK_rna[,overlap_seq_rna]) -
    lbeta(theta*pi_pK_atac_New[,overlap_seq_atac], -theta*pi_pK_atac_New[,overlap_seq_atac]+theta)
  #
  alpha <- matrix(pmin(1, exp(logPnew-logPold)*p_n2o/p_o2n), nrow=nrow(pi_pK_atac), ncol=ncol(pi_pK_atac) ) 
  updated_lab <- which(alpha>=runif(length(alpha)))
  pi_pK_atac[updated_lab] <- pi_pK_atac_New[updated_lab]
  return(pi_pK_atac)  
}

update_omega_exp <- function(pi_pK_rna, x_rna, v_rna, q_exp, pi_exp, pi_pK_atac, theta, phi, pi_p, overlap_seq_atac, overlap_seq_rna){
  # wrong: logP should not include the prior pi_p and phi
  #
  pi_pK_rna_New <- matrix(runif(length(pi_pK_rna)), nrow=nrow(pi_pK_rna), ncol=ncol(pi_pK_rna))
  #
  p_o2n <- matrix(1, nrow=nrow(pi_pK_rna), ncol=ncol(pi_pK_rna))
  p_n2o <- matrix(1, nrow=nrow(pi_pK_rna), ncol=ncol(pi_pK_rna))
  # 
  x_rnaE <- array(rep(x_rna, ncol(pi_pK_rna)), dim=c(dim(x_rna), ncol(pi_pK_rna)))
  v_rnaE <- array(rep(v_rna, nrow(pi_pK_rna)), dim=c(dim(v_rna), nrow(pi_pK_rna)))
  x_rnaE <- aperm(x_rnaE, c(2,3,1))
  v_rnaE <- aperm(v_rnaE, c(3,2,1))
  #
  lambda_rna <- calEtaLambda_exp(pi_pK_rna, q_exp, pi_exp, overlap_seq_rna)
  logPold <- rowSums(log(lambda_rna)*v_rnaE*x_rnaE + log(1-lambda_rna)*(1-v_rnaE)*x_rnaE, dims=2) 
  logPold[,overlap_seq_rna] <- logPold[,overlap_seq_rna] + (theta*pi_pK_atac[,overlap_seq_atac]-1)*log(pi_pK_rna[,overlap_seq_rna]) + 
    (-theta*pi_pK_atac[,overlap_seq_atac]+theta-1)*log(1-pi_pK_rna[,overlap_seq_rna]) 
  logPold[,-overlap_seq_rna] <- logPold[,-overlap_seq_rna] + 
    (phi*pi_p-1)*log(pi_pK_rna[,-overlap_seq_rna]) + (-phi*pi_p+phi-1)*log(1-pi_pK_rna[,-overlap_seq_rna])
  #
  lambda_rna <- calEtaLambda_exp(pi_pK_rna_New, q_exp, pi_exp, overlap_seq_rna)
  logPnew <- rowSums(log(lambda_rna)*v_rnaE*x_rnaE + log(1-lambda_rna)*(1-v_rnaE)*x_rnaE, dims=2) 
  logPnew[,overlap_seq_rna] <- logPnew[,overlap_seq_rna] + (theta*pi_pK_atac[,overlap_seq_atac]-1)*log(pi_pK_rna_New[,overlap_seq_rna]) + 
    (-theta*pi_pK_atac[,overlap_seq_atac]+theta-1)*log(1-pi_pK_rna_New[,overlap_seq_rna]) 
  logPnew[,-overlap_seq_rna] <- logPnew[,-overlap_seq_rna] + 
    (phi*pi_p-1)*log(pi_pK_rna_New[,-overlap_seq_rna]) + (-phi*pi_p+phi-1)*log(1-pi_pK_rna_New[,-overlap_seq_rna])
  #
  alpha <- matrix(pmin(1, exp(logPnew-logPold)*p_n2o/p_o2n), nrow=nrow(pi_pK_rna), ncol=ncol(pi_pK_rna) ) 
  updated_lab <- which(alpha>=runif(length(alpha)))
  pi_pK_rna[updated_lab] <- pi_pK_rna_New[updated_lab]
  return(pi_pK_rna)    
}

update_q_acc <- function(q_acc, v_acc, x_acc, omega_acc){
  # sample from uniform
  q_accNew <- runif(length(q_acc))
  #p_o2n <- tmp$p_o2n # these are 1
  #p_n2o <- tmp$p_n2o # these are 1
  ##
  logPnew <- get_q_acc_LogProb(q_accNew, omega_acc, x_acc, v_acc)
  logPold <- get_q_acc_LogProb(q_acc, omega_acc, x_acc, v_acc)
  ##
  alpha <- pmin( 1, exp(logPnew-logPold) )
  updated_lab <- which(alpha>=runif(length(alpha)))
  q_acc[updated_lab] <- q_accNew[updated_lab]
  return(q_acc)  
}

get_q_acc_LogProb <- function(q_acc, omega_acc, x_acc, v_acc){
  eta <- calEta_acc(omega_acc, q_acc) 
  tmp1 <- colSums(sweep(log(eta), c(1, 3), t(x_acc), "*"), dims=1)
  tmp2 <- colSums(sweep(log(1-eta), c(1, 3), t(x_acc), "*"), dims=1)
  LogProb <- rowSums(t(tmp1)*v_acc + t(tmp2)*(1-v_acc))
  return(LogProb)
}

update_q_exp <- function(q_exp, pi_exp, v_exp, x_exp, omega_exp, overlap_seq_exp){
  # sample from uniform
  q_expNew <- runif(length(q_exp))
  #p_o2n <- tmp$p_o2n # these are 1
  #p_n2o <- tmp$p_n2o # these are 1
  ##
  logPnew <- get_q_exp_LogProb(q_expNew, pi_exp, omega_exp, x_exp, v_exp, overlap_seq_exp)
  logPold <- get_q_exp_LogProb(q_exp, pi_exp, omega_exp, x_exp, v_exp, overlap_seq_exp)
  ##
  alpha <- pmin( 1, exp(logPnew-logPold) )
  updated_lab <- which(alpha>=runif(length(alpha)))
  q_exp[updated_lab] <- q_expNew[updated_lab]
  return(q_exp)  
}

get_q_exp_LogProb <- function(q_exp, pi_exp, omega_exp, x_exp, v_exp, overlap_seq_exp){
  etaLambda <- calEtaLambda_exp(omega_exp, q_exp, pi_exp, overlap_seq_exp)
  tmp1 <- colSums(sweep(log(etaLambda), c(1, 3), t(x_exp), "*"), dims=1)
  tmp2 <- colSums(sweep(log(1-etaLambda), c(1, 3), t(x_exp), "*"), dims=1)
  LogProb <- rowSums(t(tmp1)*v_exp + t(tmp2)*(1-v_exp))
  return(LogProb)  
}

update_pi_exp <- function(q_exp, pi_exp, v_exp_sub, x_exp, omega_exp_sub){
  ## update pi_1
  tmp <- getnewPigA1_unif(pi_exp)
  p_o2n <- tmp$p_o2n
  p_n2o <- tmp$p_n2o
  pi_expNew <- tmp$pi_gANew
  logPnew <- get_pi_exp_LogProb(t(t(pi_expNew)*q_exp), omega_exp_sub, x_exp, v_exp_sub)
  logPold <- get_pi_exp_LogProb(t(t(pi_exp)*q_exp), omega_exp_sub, x_exp, v_exp_sub)
  alpha <- pmin(1, exp(logPnew-logPold)*p_n2o/p_o2n)
  updated_lab <- which(alpha>=runif(length(alpha)))
  pi_exp[1, updated_lab] <- pi_expNew[1, updated_lab]
  
  ## update pi_0
  tmp <- getnewPigA0_unif(pi_exp)
  p_o2n <- tmp$p_o2n
  p_n2o <- tmp$p_n2o
  pi_expNew <- tmp$pi_gANew
  logPnew <- get_pi_exp_LogProb(t(t(pi_expNew)*q_exp), omega_exp_sub, x_exp, v_exp_sub)
  logPold <- get_pi_exp_LogProb(t(t(pi_exp)*q_exp), omega_exp_sub, x_exp, v_exp_sub)
  alpha <- pmin(1, exp(logPnew-logPold)*p_n2o/p_o2n)
  updated_lab <- which(alpha>=runif(length(alpha)))
  pi_exp[2, updated_lab] <- pi_expNew[2, updated_lab]
  
  return(pi_exp) 
}

update_h <- function(h, pi_pK_atac, pi_pK_rna, theta){
  #
  label2swap <- sample(length(h), 2)
  h_New <- h
  h_New[label2swap] <- rev(h_New[label2swap])
  # only ratio matters
  p_o2n <- 1
  p_n2o <- 1
  ##
  pi_pK_rna_mean <- pi_pK_atac
  pi_pK_rna_New <- pi_pK_rna[h_New,]
  ##
  logPold <- (theta*pi_pK_rna_mean-1)*log(pi_pK_rna) + (-theta*pi_pK_rna_mean+theta-1)*log(1-pi_pK_rna) 
  logPold <- sum(logPold)
  ##
  logPnew <- (theta*pi_pK_rna_mean-1)*log(pi_pK_rna_New) + (-theta*pi_pK_rna_mean+theta-1)*log(1-pi_pK_rna_New) 
  logPnew <- sum(logPnew)
  #
  alpha_prob <- min(1, exp(logPnew-logPold)*p_n2o/p_o2n)
  if (alpha_prob>=runif(1)){
    return(h_New)  
  } else {
    return(h)    
  }
} 

getnewPigA1_unif <- function(pi_gA){
  pi_gANew <- pi_gA
  pi_gANew[1,] <- pi_gA[2,] + (1-pi_gA[2,])*runif(length(pi_gA[1,]))
  p_o2n <- 1/(1-pi_gA[2,])
  p_n2o <- p_o2n
  return(list(pi_gANew=pi_gANew, p_o2n=p_o2n, p_n2o=p_n2o))  
}

getnewPigA0_unif <- function(pi_gA){
  pi_gANew <- pi_gA
  pi_gANew[2,] <- pi_gA[1,]*runif(length(pi_gA[2,]))
  p_o2n <- 1/pi_gA[1,]
  p_n2o <- p_o2n
  return(list(pi_gANew=pi_gANew, p_o2n=p_o2n, p_n2o=p_n2o))  
}

get_pi_exp_LogProb <- function(pi_gA, pi_pK, xRna, vRna){
  lambda <- outer(pi_pK, pi_gA[1,]) + outer(1-pi_pK, pi_gA[2,])
  tmp1 <- colSums(sweep(log(lambda), c(1, 3), t(xRna), "*"), dims=1)
  tmp2 <- colSums(sweep(log(1-lambda), c(1, 3), t(xRna), "*"), dims=1)
  LogProb <- rowSums(t(tmp1)*vRna + t(tmp2)*(1-vRna))
  return(LogProb)
}

cal_logpost <- function(logf1_atac, logf0_atac, logf1_rna, logf0_rna, 
                               pi_pK_atac, pi_pK_rna, pi_exp, 
                               lambda_atac, lambda_rna, 
                               v_atac, x_atac, v_rna, x_rna, 
                               mu0, mu1, nu0, nu1, 
                               atac_bin=F, rna_bin=F, overlap_seq_atac, overlap_seq_rna){
  ##
  x_atacE <- array(rep(x_atac, ncol(pi_pK_atac)), dim=c(dim(x_atac), ncol(pi_pK_atac)))
  v_atacE <- array(rep(v_atac, nrow(pi_pK_atac)), dim=c(dim(v_atac), nrow(pi_pK_atac)))
  x_atacE <- aperm(x_atacE, c(2,3,1))
  v_atacE <- aperm(v_atacE, c(3,2,1))
  ##
  x_rnaE <- array(rep(x_rna, ncol(pi_pK_rna)), dim=c(dim(x_rna), ncol(pi_pK_rna)))
  v_rnaE <- array(rep(v_rna, nrow(pi_pK_rna)), dim=c(dim(v_rna), nrow(pi_pK_rna)))
  x_rnaE <- aperm(x_rnaE, c(2,3,1))
  v_rnaE <- aperm(v_rnaE, c(3,2,1))
  ##
  term1 <- 0
  if (!atac_bin){
    tmp1 <- v_atac*logf1_atac 
    tmp1 <- tmp1[!is.na(tmp1)]
    tmp1 <- sum(tmp1)
    ##
    tmp2 <- (1-v_atac)*logf0_atac 
    tmp2 <- tmp2[!is.na(tmp2)]
    tmp2 <- sum(tmp2)
    ##
    term1 <- term1 + tmp1 + tmp2
  }
  if (!rna_bin){
    tmp1 <- v_rna*logf1_rna
    tmp1 <- tmp1[!is.na(tmp1)]
    tmp1 <- sum(tmp1)
    ##
    tmp2 <- (1-v_rna)*logf0_rna 
    tmp2 <- tmp2[!is.na(tmp2)]
    tmp2 <- sum(tmp2)
    ##
    term1 <- term1 + tmp1 + tmp2
  }
  ##
  term2 <- sum( rowSums(log(lambda_atac)*v_atacE*x_atacE + log(1-lambda_atac)*(1-v_atacE)*x_atacE, dims=2) )
  term3 <- sum( rowSums(log(lambda_rna)*v_rnaE*x_rnaE + log(1-lambda_rna)*(1-v_rnaE)*x_rnaE, dims=2) )
  ## cross data type term
  term4 <- sum( (nu1*pi_pK_atac[,overlap_seq_atac]-1)*log(pi_pK_rna[,overlap_seq_rna]) + (-nu1*pi_pK_atac[,overlap_seq_atac]+nu1-1)*log(1-pi_pK_rna[,overlap_seq_rna]) -
                  lbeta(nu1*pi_pK_atac[,overlap_seq_atac], -nu1*pi_pK_atac[,overlap_seq_atac]+nu1) )
  ## the prior, pi_pK atac
  term5 <- (nu0*mu0-1)*sum(log(pi_pK_atac)) + (-nu0*mu0+nu0-1)*sum(log(1-pi_pK_atac))
  ## the prior, pi_pK rna
  term5_2 <- (nu0*mu1-1)*sum(log(pi_pK_rna[,-overlap_seq_rna])) + (-nu0*mu1+nu0-1)*sum(log(1-pi_pK_rna[,-overlap_seq_rna]))
  ##
  term6 <- -sum( log(1-pi_exp[2,]) )
  ##
  logPost <- term1 + term2 + term3 + term4 + term5 + term5_2 + term6
  return(logPost)  
}

simData <- function(n1, n2, p, Xprobs, muG1, sigmaG1, muG2, sigmaG2,
                    mu0, mu1, nu0, nu1, q_acc, q_exp, pi_exp,            
                    theta, pi_p, phi, pi_gA_n1, pi_gA_n2, 
                    overlap_prop=0.5, diff_prop, diff_prop_atac, diff_prop_rna, cutoff=10^-6, 
                    high1=0.9, low1=0.1, high2=0.9, low2=0.1){
  nComp <- length(Xprobs)
  ## simulate pi_pk_atac
  pi_pk_atac <- matrix(rep(rbeta(n=p, shape1=mu0*nu0, shape2=-mu0*nu0+nu0), each=nComp), nrow=nComp, ncol=p)
  high <- high1
  low <- low1
  num <- round(p*diff_prop*overlap_prop/nComp)
  pi_pk_atac[,1:(num*nComp)] <- low
  for (i in 1:nComp){
    pi_pk_atac[i, (num*i-num+1):(num*i)] <- high
  }
  #
  num_atac <- round(p*(1-overlap_prop)*diff_prop_atac/nComp)
  if (num_atac!=0){
    pi_pk_atac[,round(overlap_prop*p):(round(overlap_prop*p)+num_atac*nComp-1)] <- low
    for (i in 1:nComp){
      pi_pk_atac[i, (overlap_prop*p-1 + num_atac*i-num_atac+1):(overlap_prop*p-1 + num_atac*i)] <- high
    }  
  } 
  ## simulate pi_pk_rna
  pi_pk_rna <- pi_pk_atac
  high <- high2
  low <- low2
  pi_pk_rna[,1:(num*nComp)] <- rbeta(n=nComp*num*nComp, shape1=as.numeric(pi_pk_atac[,1:(num*nComp)])*nu1, 
                                     shape2=as.numeric(-pi_pk_atac[,1:(num*nComp)])*nu1+nu1)
  #
  pi_pk_rna[,(num*nComp+1):round(p*overlap_prop)] <- rep(rbeta(n=round(p*overlap_prop)-num*nComp, 
                                                               shape1=pi_pk_atac[1,(num*nComp+1):round(p*overlap_prop)]*nu1, 
                                                               shape2=-pi_pk_atac[1,(num*nComp+1):round(p*overlap_prop)]*nu1+nu1), 
                                                         each=nComp)
  num_rna <- round(p*(1-overlap_prop)*diff_prop_rna/nComp)
  if (num_rna!=0){
    pi_pk_rna[,round(overlap_prop*p):(round(overlap_prop*p)+num_rna*nComp-1)] <- low
    for (i in 1:nComp){
      pi_pk_rna[i, (round(overlap_prop*p) + num_rna*i-num_rna):(round(overlap_prop*p)-1 + num_rna*i)] <- high
    } 
  }
  num1 <- round(overlap_prop*p)+num_rna*nComp
  if (num1<p){
    pi_pk_rna[,num1:p] <- rep(rbeta(n=p-num1+1, shape1=mu1*nu0, shape2=-mu1*nu0+nu0), each=nComp)  
  }
  
  ## implement the cutoff
  pi_pk_atac[which(pi_pk_atac<=cutoff)] <- cutoff
  pi_pk_rna[which(pi_pk_rna<=cutoff)] <- cutoff
  pi_pk_atac[which(pi_pk_atac>=(1-cutoff))] <- 1-cutoff
  pi_pk_rna[which(pi_pk_rna>=(1-cutoff))] <- 1-cutoff
  
  ## simulate x_atac
  x_atac <- t(rmultinom(n=n1, size=1, prob=Xprobs))
  cluster_atac <- rep(0, n1)
  for (i in 1:nComp){
    cluster_atac[which(x_atac[,i]==1)] <- i  
  }
  
  ## simulate x_rna
  x_rna <- t(rmultinom(n=n2, size=1, prob=Xprobs))
  cluster_rna <- rep(0, n2)
  for (i in 1:nComp){
    cluster_rna[which(x_rna[,i]==1)] <- i  
  }
  
  ## simulate u_atac
  pi_pk_atac_x <- x_atac%*%pi_pk_atac
  u_atac <- (pi_pk_atac_x>runif(length(pi_pk_atac_x))) + 0
  
  ## simulate u_rna
  pi_pk_rna_x <- x_rna%*%pi_pk_rna
  u_rna <- (pi_pk_rna_x>runif(length(pi_pk_rna_x))) + 0
  
  ## simulate v_atac
  pi_gA_1 <- matrix(rep(q_acc, p), nrow=length(q_acc))
  v_atac <- u_atac
  v_atac[which(u_atac==1)] <- (runif(length(which(u_atac==1))) < pi_gA_1[which(u_atac==1)]) + 0
  v_atac[which(u_atac==0)] <-  0
  
  ## simulate v_rna
  pi_gA_1 <- matrix(rep(pi_exp[1,]*q_exp, p), nrow=ncol(pi_exp))
  pi_gA_0 <- matrix(rep(pi_exp[2,]*q_exp, p), nrow=ncol(pi_exp))
  num_nonoverlap <- round(p*(1-overlap_prop))
  pi_gA_1[,num_nonoverlap:p] <- matrix(rep(q_exp, p-num_nonoverlap+1), nrow=ncol(pi_exp))
  pi_gA_0[,num_nonoverlap:p] <- 0
  v_rna <- u_rna
  v_rna[which(u_rna==1)] <- (runif(length(which(u_rna==1))) < pi_gA_1[which(u_rna==1)]) + 0
  v_rna[which(u_rna==0)] <- (runif(length(which(u_rna==0))) < pi_gA_0[which(u_rna==0)]) + 0
  
  ## simulate atac data 
  Data <- matrix(nrow=n1, ncol=p)
  Data[which(v_atac==1)] <- rnorm(length(which(v_atac==1)), mean=muG1[1], sd=sigmaG1[1])
  Data[which(v_atac==0)] <- rnorm(length(which(v_atac==0)), mean=muG1[2], sd=sigmaG1[2])
  f1 <- dnorm(Data, mean=muG1[1], sd=sigmaG1[1])
  f0 <- dnorm(Data, mean=muG1[2], sd=sigmaG1[2])
  Data_atac <- Data
  
  ## simulate rna data 
  Data <- matrix(nrow=n2, ncol=p)
  Data[which(v_rna==1)] <- rnorm(length(which(v_rna==1)), mean=muG2[1], sd=sigmaG2[1])
  Data[which(v_rna==0)] <- rnorm(length(which(v_rna==0)), mean=muG2[2], sd=sigmaG2[2])
  g1 <- dnorm(Data, mean=muG2[1], sd=sigmaG2[1])
  g0 <- dnorm(Data, mean=muG2[2], sd=sigmaG2[2])
  Data_rna <- Data
  
  ##
  return(list(f1=f1, f0=f0, 
              g1=g1, g0=g0,
              cluster_atac=cluster_atac, cluster_rna=cluster_rna,
              x_atac=x_atac, x_rna=x_rna, 
              u_atac=u_atac, u_rna=u_rna, 
              v_atac=v_atac, v_rna=v_rna,
              Data_atac=Data_atac, Data_rna=Data_rna,
              pi_pk_atac=pi_pk_atac, pi_pk_rna=pi_pk_rna,
              q_acc=q_acc, q_exp=q_exp, pi_exp=pi_exp))
}

update_nu1 <- function(theta, pi_pK_atac, pi_pK_rna, min, max, iter){
  if (iter%%3==0){
    # big move
    thetaNew <- runif(n=1, min = min, max = max) 
    p_o2n <- 1
    p_n2o <- 1
  } else {
    # small move
    stepsize <- theta/5
    thetaNew <- runif(n=1, min = max(min, theta-stepsize), max = min(max, theta+stepsize) )
    p_o2n <- 1/( min(max, theta+stepsize) - max(min, theta-stepsize) )
    p_n2o <- 1/( min(max, thetaNew+stepsize) - max(min, thetaNew-stepsize) )
  }
  logPold <- (theta*pi_pK_atac-1)*log(pi_pK_rna) + (-theta*pi_pK_atac+theta-1)*log(1-pi_pK_rna) -
    lbeta(theta*pi_pK_atac, -theta*pi_pK_atac+theta)
  logPold <- sum(logPold)
  #
  logPnew <- (thetaNew*pi_pK_atac-1)*log(pi_pK_rna) + (-thetaNew*pi_pK_atac+thetaNew-1)*log(1-pi_pK_rna) -
    lbeta(thetaNew*pi_pK_atac, -thetaNew*pi_pK_atac+thetaNew)
  logPnew <- sum(logPnew)
  #
  alpha <- min(1, exp(logPnew-logPold)*p_n2o/p_o2n)
  if (alpha>=runif(1)){
    return(thetaNew)  
  } else {
    return(theta)    
  }
}

update_omega_acc_nCluster1 <- function(pi_pK_atac, v_atac, q_acc, pi_pK_rna, 
                                       theta, pi_p, phi, overlap_seq_atac, overlap_seq_rna){
  #
  pi_pK_atac_New <- runif(length(pi_pK_atac))
  #
  p_o2n <- rep(1, length(pi_pK_atac))
  p_n2o <- p_o2n
  #
  p1 <- outer(q_acc, pi_pK_atac) 
  logP_vpart_old <- colSums(log(v_atac*p1 + (1-v_atac)*(1-p1)))
  #
  p1_new <- outer(q_acc, pi_pK_atac_New) 
  logP_vpart_new <- colSums(log(v_atac*p1_new + (1-v_atac)*(1-p1_new)))
  #
  logPold <- (phi*pi_p-1)*log(pi_pK_atac) + (-phi*pi_p+phi-1)*log(1-pi_pK_atac) 
  logPold[overlap_seq_atac] <- logPold[overlap_seq_atac] + (theta*pi_pK_atac[overlap_seq_atac]-1)*log(pi_pK_rna[overlap_seq_rna]) + 
    (-theta*pi_pK_atac[overlap_seq_atac]+theta-1)*log(1-pi_pK_rna[overlap_seq_rna]) -
    lbeta(theta*pi_pK_atac[overlap_seq_atac], -theta*pi_pK_atac[overlap_seq_atac]+theta)
  #
  logPnew <- (phi*pi_p-1)*log(pi_pK_atac_New) + (-phi*pi_p+phi-1)*log(1-pi_pK_atac_New) 
  logPnew[overlap_seq_atac] <- logPnew[overlap_seq_atac] + (theta*pi_pK_atac_New[overlap_seq_atac]-1)*log(pi_pK_rna[overlap_seq_rna]) + 
    (-theta*pi_pK_atac_New[overlap_seq_atac]+theta-1)*log(1-pi_pK_rna[overlap_seq_rna]) -
    lbeta(theta*pi_pK_atac_New[overlap_seq_atac], -theta*pi_pK_atac_New[overlap_seq_atac]+theta)
  #
  alpha <- pmin(1, exp(logPnew-logPold + logP_vpart_new - logP_vpart_old)*p_n2o/p_o2n)
  updated_lab <- which(alpha>=runif(length(alpha)))
  pi_pK_atac[updated_lab] <- pi_pK_atac_New[updated_lab]
  return(pi_pK_atac)  
}

update_omega_exp_nCluster1 <- function(pi_pK_rna, v_rna, q_exp, pi_exp, pi_pK_atac, theta, 
                                       phi, pi_p, overlap_seq_atac, overlap_seq_rna){
  #
  pi_pK_rna_New <- runif(length(pi_pK_rna))
  #
  p_o2n <- rep(1, length(pi_pK_rna))
  p_n2o <- p_o2n
  #
  pi_gA_rna <- pi_exp
  pi_gA_rna[1,] <- pi_gA_rna[1,]*q_exp
  pi_gA_rna[2,] <- pi_gA_rna[2,]*q_exp
  #
  p1 <- outer(pi_gA_rna[1,], pi_pK_rna) + outer(pi_gA_rna[2,], 1-pi_pK_rna)
  p1[,overlap_seq_rna] <- outer(q_exp, pi_pK_rna[overlap_seq_rna])
  logP_vpart_old <- colSums(log(v_rna*p1 + (1-v_rna)*(1-p1)))
  #
  p1_new <- outer(pi_gA_rna[1,], pi_pK_rna_New) + outer(pi_gA_rna[2,], 1-pi_pK_rna_New)
  p1_new[,overlap_seq_rna] <- outer(q_exp, pi_pK_rna_New[overlap_seq_rna])
  logP_vpart_new <- colSums(log(v_rna*p1_new + (1-v_rna)*(1-p1_new)))
  #
  logPold <- rep(1, length(pi_pK_rna))
  logPold[overlap_seq_rna] <- logPold[overlap_seq_rna] + (theta*pi_pK_atac[overlap_seq_atac]-1)*log(pi_pK_rna[overlap_seq_rna]) + 
    (-theta*pi_pK_atac[overlap_seq_atac]+theta-1)*log(1-pi_pK_rna[overlap_seq_rna]) 
  logPold[-overlap_seq_rna] <- logPold[-overlap_seq_rna] + (phi*pi_p-1)*log(pi_pK_rna[-overlap_seq_rna]) + (-phi*pi_p+phi-1)*log(1-pi_pK_rna[-overlap_seq_rna]) 
  #
  logPnew <- rep(1, length(pi_pK_rna))
  logPnew[overlap_seq_rna] <- logPnew[overlap_seq_rna] + (theta*pi_pK_atac[overlap_seq_atac]-1)*log(pi_pK_rna_New[overlap_seq_rna]) + 
    (-theta*pi_pK_atac[overlap_seq_atac]+theta-1)*log(1-pi_pK_rna_New[overlap_seq_rna]) 
  logPnew[-overlap_seq_rna] <- logPnew[-overlap_seq_rna] + (phi*pi_p-1)*log(pi_pK_rna_New[-overlap_seq_rna]) + (-phi*pi_p+phi-1)*log(1-pi_pK_rna_New[-overlap_seq_rna]) 
  #
  alpha <- pmin(1, exp(logPnew-logPold+logP_vpart_new-logP_vpart_old)*p_n2o/p_o2n) 
  updated_lab <- which(alpha>=runif(length(alpha)))
  pi_pK_rna[updated_lab] <- pi_pK_rna_New[updated_lab]
  return(pi_pK_rna)    
}

update_q_acc_nCluster1 <- function(q_acc, v, pi_pK){
  q_acc_New <- runif(length(q_acc))
  ##
  p1 <- outer(q_acc, pi_pK) 
  p1_new <- outer(q_acc_New, pi_pK)
  logP_old <- rowSums(log(v*p1 + (1-v)*(1-p1)))
  logP_new <- rowSums(log(v*p1_new + (1-v)*(1-p1_new)))
  ##
  alpha <- pmin(1, exp(logP_new-logP_old))
  updated_lab <- which(alpha>=runif(length(alpha)))
  q_acc[updated_lab] <- q_acc_New[updated_lab]
  return(q_acc)  
}

update_q_exp_nCluster1 <- function(q_exp, pi_exp, v, pi_pK, overlap_seq_exp){
  q_exp_New <- runif(length(q_exp))
  ##
  p1 <- outer(pi_exp[1,]*q_exp, pi_pK) + outer(pi_exp[2,]*q_exp, 1-pi_pK)
  p1[,overlap_seq_exp] <- outer(q_exp, pi_pK[overlap_seq_exp])
  ##
  p1_new <- outer(pi_exp[1,]*q_exp_New, pi_pK) + outer(pi_exp[2,]*q_exp_New, 1-pi_pK)
  p1_new[,overlap_seq_exp] <- outer(q_exp_New, pi_pK[overlap_seq_exp])
  ##
  logP_old <- rowSums(log(v*p1 + (1-v)*(1-p1)))
  logP_new <- rowSums(log(v*p1_new + (1-v)*(1-p1_new)))
  ##
  alpha <- pmin(1, exp(logP_new-logP_old))
  updated_lab <- which(alpha>=runif(length(alpha)))
  q_exp[updated_lab] <- q_exp_New[updated_lab]
  return(q_exp)  
}

update_pi_exp_nCluster1 <- function(q_exp, pi_exp, v, pi_pK){
  ## update pi_1
  tmp <- getnewPigA1_unif(pi_exp)
  p_o2n <- tmp$p_o2n
  p_n2o <- tmp$p_n2o
  pi_expNew <- tmp$pi_gANew
  p1 <- outer(pi_exp[1,]*q_exp, pi_pK) + outer(pi_exp[2,]*q_exp, 1-pi_pK)
  p1_new <- outer(pi_expNew[1,]*q_exp, pi_pK) + outer(pi_expNew[2,]*q_exp, 1-pi_pK)
  logP_old <- rowSums(log(v*p1 + (1-v)*(1-p1)))
  logP_new <- rowSums(log(v*p1_new + (1-v)*(1-p1_new)))
  alpha <- pmin(1, exp(logP_new-logP_old)*p_n2o/p_o2n)
  updated_lab <- which(alpha>=runif(length(alpha)))
  pi_exp[1, updated_lab] <- pi_expNew[1, updated_lab]
  
  ## update pi_0
  tmp <- getnewPigA0_unif(pi_exp)
  p_o2n <- tmp$p_o2n
  p_n2o <- tmp$p_n2o
  pi_expNew <- tmp$pi_gANew
  p1 <- outer(pi_exp[1,]*q_exp, pi_pK) + outer(pi_exp[2,]*q_exp, 1-pi_pK)
  p1_new <- outer(pi_expNew[1,]*q_exp, pi_pK) + outer(pi_expNew[2,]*q_exp, 1-pi_pK)
  logP_old <- rowSums(log(v*p1 + (1-v)*(1-p1)))
  logP_new <- rowSums(log(v*p1_new + (1-v)*(1-p1_new)))
  alpha <- pmin(1, exp(logP_new-logP_old)*p_n2o/p_o2n)
  updated_lab <- which(alpha>=runif(length(alpha)))
  pi_exp[2, updated_lab] <- pi_expNew[2, updated_lab]
  return(pi_exp) 
}

cal_density_gg_thres <- function(gg_fit, data, thres_low=-1, thres_high=20){
  lambda <- gg_fit$lambda
  shapes <- gg_fit$gamma.pars[1,]
  scales <- gg_fit$gamma.pars[2,]
  or <- order(shapes*scales)
  lambda <- lambda[or]
  shapes <- shapes[or]
  scales <- scales[or]
  prob <- lambda[1]
  ##
  perc0 <- sum(data==0)/length(data)
  p1 <- prob*(1-perc0)
  p0 <- perc0
  ##
  f1 <- dgamma(data, shape=shapes[2], scale=scales[2])
  f1[which(data==0)] <- 0
  f0 <- dgamma(data, shape=shapes[1], scale=scales[1])*p1/(p0+p1)
  f0[which(data==0)] <- 10^4*p0/(p0+p1) 
  ##
  f1[which(data>thres_high)] <- 10^4
  f1[which(data<thres_low)] <- 0
  f0[which(data>thres_high)] <- 0
  f0[which(data<thres_low)] <- 10^4
  ##
  prob0 <- p1 + p0
  return(list(f1=f1, f0=f0, prob0=prob0))
}