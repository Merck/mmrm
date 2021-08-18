## Packages ----
library(dplyr)
library(reshape2)
library(nlme)
library(mvtnorm)
library(e1071)
library(knitr)
library(car)
library(lme4)
library(Matrix)
library(magic)
library(MASS)
library(MCMCpack)
library(xtable)

## Generate data ----
sim_one_group <- function(n, k, P, sigma, beta, pi = NULL, pi_par = NULL, trt = 1, missing_type, distribution){
  p <- P
  n_pattern <- length(pi)
  covar_x <- cbind(rnorm(n), rbinom(n, 1, 0.3))
  colnames(covar_x) <- paste0("x",1:k)
  xmat <- cbind(rep(1,n),covar_x)
  y_mat <- matrix(0, n, p)
  if(distribution == "MVN")
    eps <- rmvnorm(n, mean = rep(0,p), sigma = sigma)
  else if(distribution == "MVT")
    # eps <- rmvt(n, delta = rep(0,p), sigma = sigma, df = 3)
    eps <- rmvt(n, delta = rep(0,p), sigma = sigma/5*3, df = 5)
  else if(distribution == "MVGamma"){
    u <- rmvnorm(n, mean = rep(0,p), sigma = cov2cor(sigma))
    var_eps <- diag(sigma)
    eps <- matrix(0, n, p)
    for(j in 1:p){
      shape_para <- 2
      scale_para <- sqrt(var_eps[j]/shape_para)
      eps[,j] <- qgamma(pnorm(u[,j]), shape = shape_para, scale = scale_para) - shape_para*scale_para
    }
  }
  y_mean <- apply(beta, 1, function(x) as.vector(xmat%*%x))
  y_mat <- y_mean + eps
  mean_y <- colMeans(y_mat)
  cov_eps <- cov(eps)
  cov_y <- cov(y_mat)
  db <- data.frame(y_mat)
  colnames(db) <- paste0("y",1:p)
  db <- cbind(db, covar_x)
  db$num <- c(1:n)
  db$id <- paste0(trt, "-", 1:n)
  
  if(missing_type == "MCAR"){
    db$pattern <- sample(1:p, size = n, replace = TRUE, prob = pi)
  }
  
  if(missing_type %in% c("MAR","MNAR")){
    if(pi_par[3] != 0 & missing_type == "MAR"){
      stop("This is not a MAR setting")
    }
    logit_inv <- function(x){
      exp(x) / (1 + exp(x))
    }
    
    .pattern <- rep(1, n)
    for(i_missing in p:2){
      .score <- as.matrix(data.frame(1, db[,c(i_missing - 1,i_missing)])) %*% pi_par
      .pi <- logit_inv(as.numeric(.score))
      .pattern <- ifelse( rbinom(n = n, size = 1, prob = .pi) == 1, i_missing, .pattern)
    }
    db$pattern <- .pattern
    
  }
  
  db_comb <- db
  db_comb$trt <- trt
  
  for(i in 1:nrow(db_comb)){
    for(j in 2:p){
      if(db_comb$pattern[i] <= j & db_comb$pattern[i] > 1){
        db_comb[i, j] <- NA
      }
    }
  }
  
  db_long <- reshape2::melt(db_comb, id.vars = c("id","pattern", "num", paste0("x", 1:k), "trt"),
                            variable.name = c("time") , value.name = "aval")
  db_long <- db_long %>% group_by(id) %>%
    mutate(
      time = as.numeric(time)) %>%
    ungroup()
  
  return(list(db_comb = db_comb, db_long = db_long,
              mean_y = mean_y, cov_eps = cov_eps, cov_y = cov_y))
}

sim_one_extreme <- function(n, k, P, sigma, beta, pi = NULL, pi_par = NULL, trt = 1, missing_type, distribution, outlier = "none"){
  p <- P
  n_pattern <- length(pi)
  covar_x <- cbind(rnorm(n), rbinom(n, 1, 0.3))
  colnames(covar_x) <- paste0("x",1:k)
  xmat <- cbind(rep(1,n),covar_x)
  y_mat <- matrix(0, n, p)
  if(distribution == "MVN")
    eps <- rmvnorm(n, mean = rep(0,p), sigma = sigma)
  else if(distribution == "MVT")
    # eps <- rmvt(n, delta = rep(0,p), sigma = sigma, df = 3)
    eps <- rmvt(n, delta = rep(0,p), sigma = sigma/5*3, df = 5)
  else if(distribution == "MVGamma"){
    u <- rmvnorm(n, mean = rep(0,p), sigma = cov2cor(sigma))
    var_eps <- diag(sigma)
    eps <- matrix(0, n, p)
    for(j in 1:p){
      shape_para <- 2
      scale_para <- sqrt(var_eps[j]/shape_para)
      eps[,j] <- qgamma(pnorm(u[,j]), shape = shape_para, scale = scale_para) - shape_para*scale_para
    }
  }
  y_mean <- apply(beta, 1, function(x) as.vector(xmat%*%x))
  y_mat <- y_mean + eps
  mean_y <- colMeans(y_mat)
  cov_eps <- cov(eps)
  cov_y <- cov(y_mat)
  db <- data.frame(y_mat)
  colnames(db) <- paste0("y",1:p)
  db <- cbind(db, covar_x)
  db$num <- c(1:n)
  db$id <- paste0(trt, "-", 1:n)
  
  if(missing_type == "MCAR"){
    db$pattern <- sample(1:p, size = n, replace = TRUE, prob = pi)
  }
  
  if(missing_type %in% c("MAR","MNAR")){
    if(pi_par[3] != 0 & missing_type == "MAR"){
      stop("This is not a MAR setting")
    }
    logit_inv <- function(x){
      exp(x) / (1 + exp(x))
    }
    
    .pattern <- rep(1, n)
    for(i_missing in p:2){
      .score <- as.matrix(data.frame(1, db[,c(i_missing - 1,i_missing)])) %*% pi_par
      .pi <- logit_inv(as.numeric(.score))
      .pattern <- ifelse( rbinom(n = n, size = 1, prob = .pi) == 1, i_missing, .pattern)
    }
    db$pattern <- .pattern
    
  }
  
  db_comb <- db
  db_comb$trt <- trt
  
  for(i in 1:nrow(db_comb)){
    for(j in 2:p){
      if(db_comb$pattern[i] <= j & db_comb$pattern[i] > 1){
        db_comb[i, j] <- NA
      }
    }
  }
  if(outlier != "none"){
    ind_max <- which(db_comb[,p] == max(db_comb[,p], na.rm = TRUE))
    if(outlier == "Type 1"){
      db_comb[ind_max, p] <- 3*db_comb[ind_max, p]
    }
    if(outlier == "Type 2"){
      db_comb[ind_max, 2:p] <- 3*db_comb[ind_max, 2:p]
    }
  }
  
  db_long <- reshape2::melt(db_comb, id.vars = c("id","pattern", "num", paste0("x", 1:k), "trt"),
                            variable.name = c("time") , value.name = "aval")
  db_long <- db_long %>% group_by(id) %>%
    mutate(
      time = as.numeric(time),
      trt = trt) %>%
    ungroup()
  
  return(list(db_comb = db_comb, db_long = db_long))
}

## Likelihood-based ----
j2r_llh <- function(db_comb, db_long, M){
  n1 <- length(which(db_comb$trt == 1))
  n2 <- length(which(db_comb$trt == 2))
  p <- length(unique(db_long$time))
  obs_pi <- 1 - c(sum(is.na(db_comb[,p] & db_comb$trt == 1))/n1,
                  sum(is.na(db_comb[,p] & db_comb$trt == 2))/n2)
  names_dep_var <- colnames(db_comb[, c(1:p)])
  names_covar <- colnames(db_comb[, (p+1):(p+k)])
  
  db_comb_ctl <- db_comb[which(db_comb$trt == 1),]
  db_comb_trt <- db_comb[which(db_comb$trt == 2),]
  adj_mean <- colMeans(db_comb[,(p+1):(p+k)])
  
  names_reg_var <- c(names_covar, names_dep_var[1])
  
  dep_var <- list()
  for(j in 1:p){
    if(j == 1){
      dep_var[[j]] <- names_covar
    }
    else{
      dep_var[[j]] <- c(dep_var[[j-1]], names_dep_var[j-1])
    }
  }
  name_all_var <- colnames(db_comb)
  
  
  imp_fn <- function(){
    fit_time_ctl <- list()
    fit_time_trt <- list()
    beta_time_ctl <- list()
    beta_time_trt <- list()
    beta_cov_ctl <- list()
    beta_cov_trt <- list()
    sigma_time_ctl <- c(0)
    sigma_time_trt <- c(0)
    db_imp_ctl <- db_comb_ctl
    db_imp_trt <- db_comb_trt
    for(j in 1:p){
      # control
      fit_time_ctl[[j]] <- fit_model(reformulate(dep_var[[j]],
                                                 name_all_var[j]),
                                     data = db_imp_ctl)
      beta_time_ctl[[j]] <- fit_time_ctl[[j]]$coefficients
      beta_cov_ctl[[j]] <- vcov(fit_time_ctl[[j]])
      
      # treatment
      fit_time_trt[[j]] <- fit_model(reformulate(dep_var[[j]],
                                                 name_all_var[j]),
                                     data = db_imp_trt)
      beta_time_trt[[j]] <- fit_time_trt[[j]]$coefficients
      beta_cov_trt[[j]] <- vcov(fit_time_trt[[j]])
      
      if(!identical(fit_model, Rfit::rfit)){
        sigma_time_ctl[j] <- summary(fit_time_ctl[[j]])$sigma^2
        sigma_time_trt[j] <- summary(fit_time_trt[[j]])$sigma^2
        if(j > 1){
          df_fit_ctl <- summary(fit_time_ctl[[j]])$df[1]
          df_fit_trt <- summary(fit_time_trt[[j]])$df[1]
          imp_mean_ctl <- predict(fit_time_ctl[[j]], newdata = db_imp_ctl)[which(db_imp_ctl$pattern <= j & db_imp_ctl$pattern > 1)]
          imp_mean_trt <- predict(fit_time_trt[[j]], newdata = db_imp_trt)[which(db_imp_trt$pattern <= j & db_imp_trt$pattern > 1)]
        }
      }
      
      else{
        sigma_time_ctl[j] <- fit_time_ctl[[j]]$tauhat^2
        sigma_time_trt[j] <- fit_time_trt[[j]]$tauhat^2
        df_fit_ctl <- sum(!is.na(db_comb_ctl[,j])) - (k+j)
        df_fit_trt <- sum(!is.na(db_comb_trt[,j])) - (k+j)
        if(j == 1){
          design_ctl <- cbind(rep(1, n1), data.matrix(db_imp_ctl[,(p+1):(p+k)]))
          design_trt <- cbind(rep(1, n2), data.matrix(db_imp_trt[,(p+1):(p+k)]))
        }
        else{
          design_ctl <- cbind(rep(1, n1), data.matrix(db_imp_ctl[,(p+1):(p+k)]), data.matrix(db_imp_ctl[,1:(j-1)]))
          design_trt <- cbind(rep(1, n2), data.matrix(db_imp_trt[,(p+1):(p+k)]), data.matrix(db_imp_trt[,1:(j-1)]))
        }
        pred_ctl <- as.vector(design_ctl%*%beta_time_ctl[[j]])
        imp_mean_ctl <- pred_ctl[which(db_imp_ctl$pattern <= j & db_imp_ctl$pattern > 1)]
        pred_trt <- as.vector(design_trt%*%beta_time_trt[[j]])
        imp_mean_trt <- pred_trt[which(db_imp_trt$pattern <= j & db_imp_trt$pattern > 1)]
      }
      
      if(j > 1){
        sigma_ini <- sigma_time_ctl[j]
        sigma_imp_ctl <- (df_fit_ctl-2)*sigma_ini/rchisq(1, df = df_fit_ctl)
        imp_sd_ctl <- rep(sqrt(sigma_imp_ctl), length(which(db_imp_ctl$pattern <= j & db_imp_ctl$pattern > 1)))
        db_imp_ctl[which(db_imp_ctl$pattern <= j & db_imp_ctl$pattern > 1),j] <- rnorm(length(which(db_imp_ctl$pattern <= j & db_imp_ctl$pattern > 1)), mean = imp_mean_ctl, sd = imp_sd_ctl)
        
        sigma_ini <- sigma_time_trt[j]
        sigma_imp_trt <- (df_fit_trt-2)*sigma_ini/rchisq(1, df = df_fit_trt)
        imp_sd_trt <- rep(sqrt(sigma_imp_trt), length(which(db_imp_trt$pattern <= j & db_imp_trt$pattern > 1)))
        db_imp_trt[which(db_imp_trt$pattern <= j & db_imp_trt$pattern > 1),j] <- rnorm(length(which(db_imp_trt$pattern <= j & db_imp_trt$pattern > 1)), mean = imp_mean_trt, sd = imp_sd_trt)
      }
    }
    imp_value <- c(db_imp_ctl[,p], db_imp_trt[,p])
    return(imp_value)
  }
  imp_value <- matrix(0, n1+n2, M)
  for(m in 1:M){
    imp_value[,m] <- imp_fn()
  }
  base_value <- db_comb[,1]
  chg_imp <- apply(imp_value, 2, function(x) x - base_value)
  
  # (ii) robust regression for each group to get var(\hat \beta)
  db_reg_ctl <- db_comb_ctl[,(p+1):(p+k)]
  db_reg_trt <- db_comb_trt[,(p+1):(p+k)]
  adj_mean <- colMeans(db_comb[,(p+1):(p+k)])
  coef1_mat <- matrix(c(1, adj_mean), 1, k+1)
  coef2_mat <- matrix(c(1, adj_mean), 1, k+1)
  var_xbar <- apply(db_comb[,(p+1):(p+k)], 2, var)/(n1+n2)
  
  rr_fn <- function(y){
    z1 <- y[1:n1]
    z2 <- y[(n1+1):(n1+n2)]
    db_reg_ctl$wt <- z1
    db_reg_trt$wt <- z2
    
    fit_ctl <- fit_model(wt ~ x1 + x2, data = db_reg_ctl)
    beta_ctl <- coef(fit_ctl)
    fit_trt <- fit_model(wt ~ x1 + x2, data = db_reg_trt)
    beta_trt <- coef(fit_trt)
    cov_beta_ctl <- vcov(fit_ctl)
    cov_beta_trt <- vcov(fit_trt)
    return(list(beta_ctl = beta_ctl,
                beta_trt = beta_trt,
                cov_beta_ctl = cov_beta_ctl,
                cov_beta_trt = cov_beta_trt))
  }
  mi_reg_mat <- apply(chg_imp, 2, rr_fn)
  beta_m_ctl <- sapply(mi_reg_mat, function(x) x$beta_ctl)
  beta_m_trt <- sapply(mi_reg_mat, function(x) x$beta_trt)
  beta_ctl <- rowMeans(beta_m_ctl)
  beta_trt <- rowMeans(beta_m_trt)
  cov_m_ctl <- sapply(mi_reg_mat, function(x) x$cov_beta_ctl)
  cov_m_trt <- sapply(mi_reg_mat, function(x) x$cov_beta_trt)
  
  # Rubin's estimate
  wm_ctl <- matrix(rowMeans(cov_m_ctl), 3, 3)
  bm_ctl <- cov(t(beta_m_ctl))
  var_rubin_ctl <- wm_ctl + (1+1/M)*bm_ctl
  wm_trt <- matrix(rowMeans(cov_m_trt), 3, 3)
  bm_trt <- cov(t(beta_m_trt))
  var_rubin_trt <- wm_trt + (1+1/M)*bm_trt
  
  # point estimate
  mut_ctl <- sum(coef1_mat*beta_ctl)
  mut_trt <- sum(coef2_mat*beta_trt)
  mu1 <- mut_ctl
  mu2 <- mut_trt*obs_pi[2] + mut_ctl*(1 - obs_pi[2])
  llh_est <- c(mu1, mu2, mu2 - mu1)
  
  # Variance estimate
  cov4_ctl <- var_rubin_ctl
  cov4_trt <- var_rubin_trt
  var_xbar <- cov(db_comb[,(p+1):(p+k)])/(n1+n2)
  coef_mat <- matrix(c(1, adj_mean), nrow = 1)
  coef_beta_mat_ctl <- matrix(beta_ctl[2:(k+1)], nrow = 1)
  coef_beta_mat_trt <- matrix(beta_trt[2:(k+1)]*obs_pi[2] + beta_ctl[2:(k+1)]*(1 - obs_pi[2]), nrow = 1)
  coef_beta_mat_diff <- matrix((beta_trt[2:(k+1)] - beta_ctl[2:(k+1)])*obs_pi[2], nrow = 1)
  var_theta_ctl <- as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat)) +
    as.numeric(coef_beta_mat_ctl%*%var_xbar%*%t(coef_beta_mat_ctl))
  var_theta_trt <- as.numeric(coef_mat%*%cov4_trt%*%t(coef_mat))*obs_pi[2]^2 +
    as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat))*(1 - obs_pi[2])^2 +
    (mut_trt - mut_ctl)^2*obs_pi[2]*(1 - obs_pi[2])/n2 +
    as.numeric(coef_beta_mat_trt%*%var_xbar%*%t(coef_beta_mat_trt))
  var_theta_diff <- as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat))*obs_pi[2]^2 +
    as.numeric(coef_mat%*%cov4_trt%*%t(coef_mat))*obs_pi[2]^2 +
    (mut_trt - mut_ctl)^2*obs_pi[2]*(1 - obs_pi[2])/n2 +
    as.numeric(coef_beta_mat_diff%*%var_xbar%*%t(coef_beta_mat_diff))
  var_llh <- c(var_theta_ctl, var_theta_trt, var_theta_diff)
  
  return(list(llh_est = llh_est,
              var_llh = var_llh))
}

## nonparametric bootstrap ----
nonpara_llh_fn <- function(db_comb, db_long, B, M){
  db_comb_ctl <- db_comb[which(db_comb$trt == 1),]
  db_comb_trt <- db_comb[which(db_comb$trt == 2),]
  llh_boot <- matrix(0, 3, B)
  
  for(b in 1:B){
    set.seed(b)
    id_boot <- sample(unique(db_comb$id), replace = TRUE)
    data_boot <- data.frame(id_new = 1:length(id_boot), id = id_boot)
    data_boot <- merge(data_boot, db_comb, all.x = TRUE)
    data_boot <- cbind(data_boot[,-c(1:2)], data_boot[,1:2])
    colnames(data_boot) <- c(colnames(data_boot)[1:(ncol(data_boot)-2)],"id_old", "id")
    data_long <- data.frame(id_new = 1:length(id_boot), id = id_boot)
    data_long <- merge(data_long, db_long, all.x = TRUE)
    colnames(data_long) <- c("id_old", "id", colnames(data_long)[3:ncol(data_long)])
    llh_res <- j2r_llh(data_boot, data_long, M)
    llh_boot[,b] <- llh_res$llh_est
  }
  var_llh_boot <- apply(llh_boot, 1, var)
  
  return(list(var_llh_boot = var_llh_boot))
}

## main function ----
main<-function(seed, case){
  #the main fucntion to run for one simulation
  set.seed(seed)
  if(missing_type == "MCAR"){
    tmp1 <- sim_one_group(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi = Pi_ctl, trt = 1, missing_type = missing_type, distribution = distribution)
    # tmp1 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi = Pi_ctl, trt = 1, missing_type = missing_type, distribution = distribution)
    tmp2 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_trt, pi = Pi_trt, trt = 2, missing_type = missing_type, distribution = distribution)
  }
  if(missing_type %in% c("MAR", "MNAR")){
    if(case == "Case0"){
      tmp1 <- sim_one_group(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi_par = phi_ctl, trt = 1, missing_type = missing_type, distribution = distribution)
      tmp2 <- sim_one_group(n = N, k = k, P = p, sigma = Sigma, beta = beta_trt, pi_par = phi_trt, trt = 2, missing_type = missing_type, distribution = distribution)
    }
    if(case == "Case1"){
      tmp1 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi_par = phi_ctl, trt = 1, missing_type = missing_type, distribution = distribution, outlier = "Type 1")
      tmp2 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_trt, pi_par = phi_trt, trt = 2, missing_type = missing_type, distribution = distribution, outlier = "Type 1")
    }
    if(case == "Case2"){
      tmp1 <- sim_one_group(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi_par = phi_ctl, trt = 1, missing_type = missing_type, distribution = distribution)
      tmp2 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_trt, pi_par = phi_trt, trt = 2, missing_type = missing_type, distribution = distribution, outlier = "Type 1")
    }
    if(case == "Case3"){
      tmp1 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi_par = phi_ctl, trt = 1, missing_type = missing_type, distribution = distribution, outlier = "Type 1")
      tmp2 <- sim_one_group(n = N, k = k, P = p, sigma = Sigma, beta = beta_trt, pi_par = phi_trt, trt = 2, missing_type = missing_type, distribution = distribution)
    }
    if(case == "Case4"){
      tmp1 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi_par = phi_ctl, trt = 1, missing_type = missing_type, distribution = distribution, outlier = "Type 2")
      tmp2 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_trt, pi_par = phi_trt, trt = 2, missing_type = missing_type, distribution = distribution, outlier = "Type 2")
    }
    if(case == "Case5"){
      tmp1 <- sim_one_group(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi_par = phi_ctl, trt = 1, missing_type = missing_type, distribution = distribution)
      tmp2 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_trt, pi_par = phi_trt, trt = 2, missing_type = missing_type, distribution = distribution, outlier = "Type 2")
    }
    if(case == "Case6"){
      tmp1 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi_par = phi_ctl, trt = 1, missing_type = missing_type, distribution = distribution, outlier = "Type 2")
      tmp2 <- sim_one_group(n = N, k = k, P = p, sigma = Sigma, beta = beta_trt, pi_par = phi_trt, trt = 2, missing_type = missing_type, distribution = distribution)
    }
  }
  db_comb <- rbind(tmp1[[1]], tmp2[[1]]) # to fit lm
  db_long <- rbind(tmp1[[2]], tmp2[[2]]) # to fit mmrm
  
  ## Frank's likelihood-based method
  llh_res <- j2r_llh(db_comb, db_long, M)
  llh_boot <- nonpara_llh_fn(db_comb, db_long, B, M)
  
  return(list(llh_est = llh_res$llh_est,
              llh_var = llh_res$var_llh,
              llh_boot_var = llh_boot$var_llh_boot
  ))
}

## Data ----
N <- 100
k <- 2 # dimension of covariates (omit intercept)
p <- 5 # number of visits
mu_beta_ctl <- c(0, 1, 2, 3, 4)
# mu_beta_trt <- c(0, 1.3, 2.8, 4, 5.5)
mu_beta_trt <- c(0, 1.3, 2.3, 3.5, 4.8)
set.seed(123)
beta_ctl <- rbind(rnorm(k+1, mu_beta_ctl[1], 1), rnorm(k+1,mu_beta_ctl[2], 1),
                  rnorm(k+1, mu_beta_ctl[3], 1), rnorm(k+1,mu_beta_ctl[4], 1),
                  rnorm(k+1, mu_beta_ctl[5], 1))
# beta_trt <- rbind(rnorm(k+1, mu_beta_trt[1], 1), rnorm(k+1,mu_beta_trt[2], 1),
#                   rnorm(k+1, mu_beta_trt[3], 1), rnorm(k+1,mu_beta_trt[4], 1),
#                   rnorm(k+1, mu_beta_trt[5], 1))
# beta_trt[1,] <- beta_ctl[1,]
## Under H0
beta_trt <- beta_ctl

sd <- c(2.0, 1.8, 2.0, 2.1, 2.2)
corr <- matrix(
  c(1, 0.6, 0.3, 0.2, 0.1,
    0.6, 1, 0.7, 0.5, 0.2,
    0.3, 0.7, 1, 0.6, 0.4,
    0.2, 0.5, 0.6, 1, 0.5,
    0.1, 0.2, 0.4, 0.5, 1), 5, 5)
Sigma <- diag(sd) %*% corr %*% diag(sd)

## Missing types
# (1) MCAR
# missing_type = "MCAR"
# (i) case 1
# Pi_ctl <- c(80, rep(5,4)) / 100
# Pi_trt <- c(80, rep(5,4)) / 100
# (ii) case 2
# Pi_ctl <- c(80, rep(5,4)) / 100
# Pi_trt <- c(70, 9, rep(7,3)) / 100
# (iii) case 3
# Pi_ctl <- c(70, 9, rep(7,3)) / 100
# Pi_trt <- c(80, rep(5,4)) / 100
# (iv) case 4
# Pi_ctl <- c(60, rep(10,4)) / 100
# Pi_trt <- c(60, rep(10,4)) / 100
# (2) MAR
missing_type = "MAR"
## (i) ctl > trt, approx 0.8, 0.7
# phi_ctl <- c(-3.2, 0.2, 0)
# phi_trt <- c(-3.5, 0.2, 0)
## (ii) ctl = trt, approx 0.8
# phi_ctl <- c(-3.5, 0.2, 0)
# phi_trt <- c(-3.6, 0.2, 0)
# ## (iii) trt > ctl, approxiamte 0.7, 0.8
# phi_ctl <- c(-2.8, 0.2, 0)
# phi_trt <- c(-4.0, 0.2, 0)
## Under H0
phi_ctl <- c(-3.5, 0.2, 0)
phi_trt <- c(-3.5, 0.2, 0)

## Distribution
distribution = "MVN"
# distribution = "MVT"
# distribution = "MVGamma"

## True value ----
if(distribution == "MVN"){
  # true_mean <- c(5.526485, 6.604056, 1.077571) # MVN
  true_mean <- c(5.526485, 5.526485, 0)
}
if(distribution == "MVT"){
  # true_mean <- c(5.526554, 6.607928, 1.081374) # MVT
  true_mean <- c(5.526554, 5.526554, 0)
}
if(distribution == "MVGamma"){
  # true_mean <- c(5.526218, 6.603291, 1.077073) # MVGamma
  true_mean <- c(5.526218, 5.526218, 0)
}

## Simulation results ----
M <- 10
B <- 100
# fit_model <- lm
# fit_model <- MASS::rlm
fit_model <- Rfit::rfit
# Get task id for each simulation job (not run in Rstudio Serve)
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
res <- main(seed = task_id, case = "Case6")

# Save all the objects into i.Rdata
save(res, file = paste0(task_id, ".Rdata"))
