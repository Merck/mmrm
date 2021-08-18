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

## j2r imputation (proper MI) ----
j2r_imp <- function(db_comb, db_long, M){
  n1 <- length(which(db_comb$trt == 1))
  n2 <- length(which(db_comb$trt == 2))
  p <- length(unique(db_long$time))
  obs_pi <- 1 - c(sum(is.na(db_comb[,p] & db_comb$trt == 1))/n1,
                  sum(is.na(db_comb[,p] & db_comb$trt == 2))/n2)
  names_dep_var <- colnames(db_comb[, c(1:p)])
  names_covar <- colnames(db_comb[, c((p+1):(p+k))])
  db_comb_ctl <- db_comb[which(db_comb$trt == 1),]
  db_comb_trt <- db_comb[which(db_comb$trt == 2),]
  db_long_ctl <- db_long[which(db_long$trt == 1),]
  db_long_trt <- db_long[which(db_long$trt == 2),]
  db_avaiable <- na.omit(db_long)
  db_avaiable_ctl <- na.omit(db_long_ctl)
  db_avaiable_trt <- na.omit(db_long_trt)
  child_ctl <- factor(db_avaiable_ctl$id)
  child_trt <- factor(db_avaiable_trt$id)
  fit_lmer_ctl <- lmer(aval ~ x1 + x2 + factor(time) +
                         x1:factor(time) + x2:factor(time) +
                         (0 + factor(time)|id),
                       control = lmerControl(check.nobs.vs.nRE = "ignore",
                                             check.conv.grad = "ignore",
                                             check.conv.hess = "ignore",
                                             check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),
                       data = db_avaiable_ctl)
  fit_lmer_trt <- lmer(aval ~ x1 + x2 + factor(time) +
                         x1:factor(time) + x2:factor(time) +
                         (0 + factor(time)|id),
                       control = lmerControl(check.nobs.vs.nRE = "ignore",
                                             check.conv.grad = "ignore",
                                             check.conv.hess = "ignore",
                                             check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),
                       data = db_avaiable_trt)
  # may take more time using lsmeans()
  # lsmeans(fit_lmer_ctl, list(pairwise ~ factor(time)), df = nrow(db_long_ctl) - 1, adjust = NULL)
  beta_lmer_ctl <- fixef(fit_lmer_ctl)
  cov_beta_ctl <- vcov(fit_lmer_ctl, full = TRUE, ranpar = "var")
  cov_beta_ctl <- data.matrix(cov_beta_ctl)
  # draw beta
  beta_mi_ctl <- rmvnorm(M, mean = beta_lmer_ctl, sigma = cov_beta_ctl)
  beta_transform_ctl <- apply(beta_mi_ctl, 1, function(x) c(x[1:3],
                                                            x[1:3] + c(x[4], x[8], x[12]),
                                                            x[1:3] + c(x[5], x[9], x[13]),
                                                            x[1:3] + c(x[6], x[10], x[14]),
                                                            x[1:3] + c(x[7], x[11], x[15])))
  beta_lmer_trt <- fixef(fit_lmer_trt)
  cov_beta_trt <- vcov(fit_lmer_trt, full = TRUE, ranpar = "var")
  cov_beta_trt <- data.matrix(cov_beta_trt)
  # draw beta
  beta_mi_trt <- rmvnorm(M, mean = beta_lmer_trt, sigma = cov_beta_trt)
  beta_transform_trt <- apply(beta_mi_trt, 1, function(x) c(x[1:3],
                                                            x[1:3] + c(x[4], x[8], x[12]),
                                                            x[1:3] + c(x[5], x[9], x[13]),
                                                            x[1:3] + c(x[6], x[10], x[14]),
                                                            x[1:3] + c(x[7], x[11], x[15])))
  # random effect covariance matrix D
  vc.a <- VarCorr(fit_lmer_ctl)
  vc.da <- as.data.frame(vc.a,order="lower.tri")
  var_D_ctl <- matrix(c(vc.da[1:5,4],
                        vc.da[2,4], vc.da[6:9,4],
                        vc.da[3,4], vc.da[7,4], vc.da[10:12,4],
                        vc.da[4,4], vc.da[8,4], vc.da[11,4],vc.da[13:14,4],
                        vc.da[5,4], vc.da[9,4], vc.da[12,4],vc.da[14,4], vc.da[15,4]
  ),p,p,byrow=TRUE)
  var_R <- summary(fit_lmer_ctl)$sigma^2*diag(p)
  Z_i <- diag(p)
  V_ctl <- var_R + t(Z_i)%*%var_D_ctl%*%Z_i
  
  db_imp <- db_comb
  imp_fn <- function(beta){
    beta_imp_ctl <- matrix(beta[1:15], 1+k, p)
    beta_imp_trt <- matrix(beta[16:30], 1+k, p)
    V_mi_ctl <- rwish(p+1, V_ctl)/(p+1)
    
    control_term1 <- as.vector(beta_imp_ctl[1,])
    control_term2 <- beta_imp_ctl[-1,]
    trt_term1 <- as.vector(beta_imp_trt[1,])
    trt_term2 <- beta_imp_trt[-1,]
    imp_value_ctl <- db_comb_ctl[,p]
    imp_value_trt <- db_comb_trt[,p]
    sigma_pattern <- c(0)
    mean_ctl <- rep(0, n1)
    mean_trt <- rep(0, n2)
    for(j in 2:p){
      sigma_mo <- V_mi_ctl[j:p, 1:(j-1)]%*%solve(V_mi_ctl[1:(j-1), 1:(j-1)])
      sigma_pattern[j-1] <- (V_mi_ctl[j:p, j:p] - sigma_mo%*%V_mi_ctl[1:(j-1),j:p])[p-j+1, p-j+1]
      if(length(which(db_comb_ctl$pattern == j)) > 0){
        pattern_mat_ctl <- db_comb_ctl[which(db_comb_ctl$pattern == j), ]
        mu_miss_ctl <- control_term1[j:p] + apply(pattern_mat_ctl[,(p+1):(p+k)], 1, function(x) apply(matrix(control_term2[,j:p], k, (p-j)+1), 2, function(y) sum(x*y)))
        mu_obs_ctl <- control_term1[1:(j-1)] + apply(pattern_mat_ctl[,(p+1):(p+k)], 1, function(x) apply(matrix(control_term2[,1:(j-1)], k, j-1), 2, function(y) sum(x*y)))
        mean_pattern_ctl <- mu_miss_ctl + apply(matrix(data.matrix(pattern_mat_ctl[,1:(j-1)] - t(mu_obs_ctl)), nrow(pattern_mat_ctl), j-1), 1, function(x) sigma_mo%*%x)
        mean_pattern_ctl <- matrix(mean_pattern_ctl, p-j+1, nrow(pattern_mat_ctl))
        mean_ctl[which(db_comb_ctl$pattern == j)] <- mean_pattern_ctl[p-j+1,]
        
        imp_value_ctl[which(db_comb_ctl$pattern == j)] <- sapply(mean_pattern_ctl[p-j+1,], function(x) rnorm(1, x, sqrt(sigma_pattern[j-1])))
      }
      if(length(which(db_comb_trt$pattern == j)) > 0){
        pattern_mat_trt <- db_comb_trt[which(db_comb_trt$pattern == j), ]
        mu_miss_trt <- control_term1[j:p] + apply(pattern_mat_trt[,(p+1):(p+k)], 1, function(x) apply(matrix(control_term2[,j:p], k, (p-j)+1), 2, function(y) sum(x*y)))
        mu_obs_trt <- trt_term1[1:(j-1)] + apply(pattern_mat_trt[,(p+1):(p+k)], 1, function(x) apply(matrix(trt_term2[,1:(j-1)], k, j-1), 2, function(y) sum(x*y)))
        mean_pattern_trt <- mu_miss_trt + apply(matrix(data.matrix(pattern_mat_trt[,1:(j-1)] - t(mu_obs_trt)), nrow(pattern_mat_trt), j-1), 1, function(x) sigma_mo%*%x)
        mean_pattern_trt <- matrix(mean_pattern_trt, p-j+1, nrow(pattern_mat_trt))
        mean_trt[which(db_comb_trt$pattern == j)] <- mean_pattern_trt[p-j+1,]
        imp_value_trt[which(db_comb_trt$pattern == j)] <- sapply(mean_pattern_trt[p-j+1,], function(x) rnorm(1, x, sqrt(sigma_pattern[j-1])))
      }
    }
    
    imp_value <- c(imp_value_ctl, imp_value_trt)
    
    return(imp_value)
  }
  imp_value <- apply(rbind(beta_transform_ctl, beta_transform_trt), 2, imp_fn)
  
  # change from baseline
  base_value <- c(db_comb_ctl[,1], db_comb_trt[,1])
  chg_imp <- t(apply(imp_value, 2, function(x) x - base_value))
  chg_imp_ctl <- chg_imp[,1:n1]
  chg_imp_trt <- chg_imp[,(n1+1):(n1+n2)]
  
  return(chg_imp = chg_imp)
}

rubin_est <- function(db_comb, db_long, chg_imp, M, fit_model){
  n1 <- length(which(db_comb$trt == 1))
  n2 <- length(which(db_comb$trt == 2))
  p <- length(unique(db_long$time))
  obs_pi <- 1 - c(sum(is.na(db_comb[,p] & db_comb$trt == 1))/n1,
                  sum(is.na(db_comb[,p] & db_comb$trt == 2))/n2)
  names_dep_var <- colnames(db_comb[, 1:p])
  names_covar <- colnames(db_comb[, (p+1):(p+k)])
  db_comb_ctl <- db_comb[which(db_comb$trt == 1),]
  db_comb_trt <- db_comb[which(db_comb$trt == 2),]
  chg_imp_ctl <- chg_imp[,1:n1]
  chg_imp_trt <- chg_imp[,(n1+1):(n1+n2)]
  
  ## point estimates
  # (1) mean estimator
  # (i) simple average
  theta1 <- mean(chg_imp_ctl)
  theta2 <- mean(chg_imp_trt)
  theta_diff <-  theta2 - theta1
  mean_est <- c(theta1, theta2, theta_diff)
  
  # Rubin's estimate
  wm_ctl <- mean(apply(chg_imp_ctl, 1, var)/n1)
  bm_ctl <- var(apply(chg_imp_ctl, 1, mean))
  var_rubin_ctl <- wm_ctl + (1+1/M)*bm_ctl
  wm_trt <- mean(apply(chg_imp_trt, 1, var)/n2)
  bm_trt <- var(apply(chg_imp_trt, 1, mean))
  var_rubin_trt <- wm_trt + (1+1/M)*bm_trt
  wm_diff <- wm_ctl + wm_trt
  bm_diff <- var(apply(chg_imp_trt, 1, mean) - apply(chg_imp_ctl, 1, mean))
  var_rubin_diff <- wm_diff + (1+1/M)*bm_diff
  var_mean_rubin <- c(var_rubin_ctl, var_rubin_trt, var_rubin_diff)
  
  # (ii) robust regression for each group
  db_reg_ctl <- db_comb_ctl[,c(1,(p+1):(p+k))]
  db_reg_trt <- db_comb_trt[,c(1,(p+1):(p+k))]
  adj_mean <- colMeans(db_comb[,(p+1):(p+k)])
  coef1_mat <- matrix(c(1, adj_mean), 1, k+1)
  coef2_mat <- matrix(c(1, adj_mean), 1, k+1)
  var_xbar <- apply(db_comb[,(p+1):(p+k)], 2, var)/(n1+n2)
  
  rr_fn <- function(y){
    z1 <- y[1:n1]
    z2 <- y[(n1+1):(n1+n2)]
    db_reg_ctl$wt <- z1
    db_reg_trt$wt <- z2
    
    # using Huber weight
    # fit_ctl <- rlm(wt ~ x1 + x2, data = db_reg_ctl)
    fit_ctl <- fit_model(wt ~ x1 + x2, data = db_reg_ctl)
    # using bisquare weight
    # fit_ctl <- rlm(wt ~ x1 + x2, data = db_reg_ctl,
    #                psi = psi.bisquare)
    beta_ctl <- coef(fit_ctl)
    # fit_trt <- rlm(wt ~ x1 + x2, data = db_reg_trt)
    fit_trt <- fit_model(wt ~ x1 + x2, data = db_reg_trt)
    beta_trt <- coef(fit_trt)
    theta1 <- beta_ctl[1] + sum(beta_ctl[2:(k+1)]*adj_mean)
    theta2 <- beta_trt[1] + sum(beta_trt[2:(k+1)]*adj_mean)
    theta_diff <- theta2 - theta1
    mean_est <- c(theta1, theta2, theta_diff)
    cov_beta_ctl <- vcov(fit_ctl)
    cov_beta_trt <- vcov(fit_trt)
    var_theta1 <- as.numeric(coef1_mat%*%cov_beta_ctl%*%t(coef1_mat)) +
      sum(beta_ctl[2:(k+1)]^2*var_xbar)
    var_theta2 <- as.numeric(coef2_mat%*%cov_beta_trt%*%t(coef2_mat)) +
      sum(beta_trt[2:(k+1)]^2*var_xbar)
    var_diff <- as.numeric(coef1_mat%*%cov_beta_ctl%*%t(coef1_mat)) +
      as.numeric(coef2_mat%*%cov_beta_trt%*%t(coef2_mat)) +
      sum((beta_trt[2:(k+1)] - beta_ctl[2:(k+1)])^2*var_xbar)
    var_est <- c(var_theta1, var_theta2, var_diff)
    return(c(mean_est, var_est))
  }
  mi_reg_mat <- apply(chg_imp, 1, rr_fn)
  reg_est <- rowMeans(mi_reg_mat)[1:3]
  
  # Rubin's estimate
  wm_ctl <- rowMeans(mi_reg_mat)[4]
  bm_ctl <- var(mi_reg_mat[1,])
  var_rubin_ctl <- wm_ctl + (1+1/M)*bm_ctl
  wm_trt <- rowMeans(mi_reg_mat)[5]
  bm_trt <- var(mi_reg_mat[2,])
  var_rubin_trt <- wm_trt + (1+1/M)*bm_trt
  wm_diff <- rowMeans(mi_reg_mat)[6]
  bm_diff <- var(mi_reg_mat[3,])
  var_rubin_diff <- wm_diff + (1+1/M)*bm_diff
  var_reg_rubin <- c(var_rubin_ctl, var_rubin_trt, var_rubin_diff)
  
  # (iii) robust regression for the whole data
  db_reg <- db_comb[,c("x1", "x2", "trt")]
  coef1_mat <- matrix(c(1, adj_mean, 0), 1, k+2)
  coef2_mat <- matrix(c(1, adj_mean, 1), 1, k+2)
  coef3_mat <- coef2_mat - coef1_mat
  var_xbar <- apply(db_comb[,(p+1):(p+k)], 2, var)/(n1+n2)
  
  rr2_fn <- function(y){
    z <- y
    db_reg$wt <- z
    
    # (1)
    # using Huber weight
    # fit_rr <- rlm(wt ~ x1 + x2 + factor(trt), data = db_reg)
    fit_rr <- fit_model(wt ~ x1 + x2 + factor(trt), data = db_reg)
    # using bisquare weight
    # fit_ctl <- rlm(wt ~ x1 + x2, data = db_reg_ctl,
    #                psi = psi.bisquare)
    beta_rr <- coef(fit_rr)
    theta1 <- beta_rr[1] + sum(beta_rr[2:(k+1)]*adj_mean)
    theta2 <- beta_rr[1] + sum(beta_rr[2:(k+1)]*adj_mean) + beta_rr[k+2]
    theta_diff <- theta2 - theta1
    mean_est <- c(theta1, theta2, theta_diff)
    cov_beta_rr <- vcov(fit_rr)
    var_theta1 <- as.numeric(coef1_mat%*%cov_beta_rr%*%t(coef1_mat)) +
      sum(beta_rr[2:(k+1)]^2*var_xbar)
    var_theta2 <- as.numeric(coef2_mat%*%cov_beta_rr%*%t(coef2_mat)) +
      sum(beta_rr[2:(k+1)]^2*var_xbar)
    var_diff <- as.numeric(coef3_mat%*%cov_beta_rr%*%t(coef3_mat))
    var_est <- c(var_theta1, var_theta2, var_diff)
    return(c(mean_est, var_est))
  }
  mi_reg2_mat <- apply(chg_imp, 1, rr2_fn)
  reg2_est <- rowMeans(mi_reg2_mat)[1:3]
  
  # Rubin's estimate
  wm_ctl <- rowMeans(mi_reg2_mat)[4]
  bm_ctl <- var(mi_reg2_mat[1,])
  var_rubin_ctl <- wm_ctl + (1+1/M)*bm_ctl
  wm_trt <- rowMeans(mi_reg2_mat)[5]
  bm_trt <- var(mi_reg2_mat[2,])
  var_rubin_trt <- wm_trt + (1+1/M)*bm_trt
  wm_diff <- rowMeans(mi_reg2_mat)[6]
  bm_diff <- var(mi_reg2_mat[3,])
  var_rubin_diff <- wm_diff + (1+1/M)*bm_diff
  var_reg2_rubin <- c(var_rubin_ctl, var_rubin_trt, var_rubin_diff)
  
  return(list(mean_est = mean_est,
              var_mean_rubin = var_mean_rubin,
              reg_est = reg_est,
              var_reg_rubin = var_reg_rubin,
              reg2_est = reg2_est,
              var_reg2_rubin = var_reg2_rubin))
}

rubin_mean <- function(db_comb, db_long, chg_imp, fit_model){
  n1 <- length(which(db_comb$trt == 1))
  n2 <- length(which(db_comb$trt == 2))
  p <- length(unique(db_long$time))
  obs_pi <- 1 - c(sum(is.na(db_comb[,p] & db_comb$trt == 1))/n1,
                  sum(is.na(db_comb[,p] & db_comb$trt == 2))/n2)
  names_dep_var <- colnames(db_comb[, c(1:p)])
  names_covar <- colnames(db_comb[, c((p+1):(p+k))])
  db_comb_ctl <- db_comb[which(db_comb$trt == 1),]
  db_comb_trt <- db_comb[which(db_comb$trt == 2),]
  chg_imp_ctl <- chg_imp[,1:n1]
  chg_imp_trt <- chg_imp[,(n1+1):(n1+n2)]
  
  ## point estimates
  # (1) mean estimator
  # (i) simple average
  theta1 <- mean(chg_imp_ctl)
  theta2 <- mean(chg_imp_trt)
  theta_diff <-  theta2 - theta1
  mean_est <- c(theta1, theta2, theta_diff)
  
  # (ii) robust regression for each group
  db_reg_ctl <- db_comb_ctl[,(p+1):(p+k)]
  db_reg_trt <- db_comb_trt[,(p+1):(p+k)]
  adj_mean <- colMeans(db_comb[,(p+1):(p+k)])
  
  rr_fn <- function(y){
    z1 <- y[1:n1]
    z2 <- y[(n1+1):(n1+n2)]
    db_reg_ctl$wt <- z1
    db_reg_trt$wt <- z2
    
    # using Huber weight
    fit_ctl <- fit_model(wt ~ x1 + x2, data = db_reg_ctl)
    # using bisquare weight
    # fit_ctl <- rlm(wt ~ x1 + x2, data = db_reg_ctl,
    #                psi = psi.bisquare)
    beta_ctl <- coef(fit_ctl)
    fit_trt <- fit_model(wt ~ x1 + x2, data = db_reg_trt)
    beta_trt <- coef(fit_trt)
    theta1 <- beta_ctl[1] + sum(beta_ctl[2:(k+1)]*adj_mean)
    theta2 <- beta_trt[1] + sum(beta_trt[2:(k+1)]*adj_mean)
    theta_diff <- theta2 - theta1
    mean_est <- c(theta1, theta2, theta_diff)
    return(mean_est)
  }
  mi_reg_mat <- apply(chg_imp, 1, rr_fn)
  reg_est <- rowMeans(mi_reg_mat)
  
  # (iii) robust regression for the whole data
  db_reg <- db_comb[,c("x1", "x2", "trt")]
  coef1_mat <- matrix(c(1, adj_mean, 0), 1, k+2)
  coef2_mat <- matrix(c(1, adj_mean, 1), 1, k+2)
  coef3_mat <- coef2_mat - coef1_mat
  var_xbar <- apply(db_comb[,(p+1):(p+k)], 2, var)/(n1+n2)
  
  rr2_fn <- function(y){
    z <- y
    db_reg$wt <- z
    
    # (1)
    # using Huber weight
    # fit_rr <- rlm(wt ~ x1 + x2 + factor(trt), data = db_reg)
    fit_rr <- fit_model(wt ~ x1 + x2 + factor(trt), data = db_reg)
    beta_rr <- coef(fit_rr)
    theta1 <- beta_rr[1] + sum(beta_rr[2:(k+1)]*adj_mean)
    theta2 <- beta_rr[1] + sum(beta_rr[2:(k+1)]*adj_mean) + beta_rr[k+2]
    theta_diff <- theta2 - theta1
    mean_est <- c(theta1, theta2, theta_diff)
    return(mean_est)
  }
  mi_reg2_mat <- apply(chg_imp, 1, rr2_fn)
  reg2_est <- rowMeans(mi_reg2_mat)
  
  return(list(mean_est = mean_est,
              reg_est = reg_est,
              reg2_est = reg2_est))
}

nonpara_mi_fn <- function(db_comb, db_long, M, B, fit_model){
  db_comb_ctl <- db_comb[which(db_comb$trt == 1),]
  db_comb_trt <- db_comb[which(db_comb$trt == 2),]
  mean_boot <- matrix(0, 3, B)
  reg_boot <- matrix(0, 3, B)
  reg2_boot <- matrix(0, 3, B)
  
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
    chg_imp <- j2r_imp(data_boot, data_long, M)
    mi_res <- rubin_mean(data_boot, data_long, chg_imp, fit_model)
    mean_boot[,b] <- mi_res$mean_est
    reg_boot[,b] <- mi_res$reg_est
    reg2_boot[,b] <- mi_res$reg2_est
  }
  var_mean_boot <- apply(mean_boot, 1, var)
  var_reg_boot <- apply(reg_boot, 1, var)
  var_reg2_boot <- apply(reg2_boot, 1, var)
  
  return(list(var_mean_boot = var_mean_boot,
              var_reg_boot = var_reg_boot,
              var_reg2_boot = var_reg2_boot))
}

## Simulate function ----
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
  
  ## multiple imputation
  # chg_imp <- rr_seq_imp(db_comb, db_long, M)
  chg_imp <- j2r_imp(db_comb, db_long, M)
  mi_res <- rubin_est(db_comb, db_long, chg_imp, M, fit_model)
  minp_res <- nonpara_mi_fn(db_comb, db_long, M, B, fit_model)
  mean_mi <- mi_res$mean_est
  reg_mi <- mi_res$reg_est
  reg2_mi <- mi_res$reg2_est
  var_rubin_mean <- mi_res$var_mean_rubin
  var_rubin_reg <- mi_res$var_reg_rubin
  var_rubin_reg2 <- mi_res$var_reg2_rubin
  
  var_boot_mean <- minp_res$var_mean_boot
  var_boot_reg <- minp_res$var_reg_boot
  var_boot_reg2 <- minp_res$var_reg2_boot
  
  return(list(mean_mi = mean_mi,
              reg_mi = reg_mi,
              reg2_mi = reg2_mi,
              var_rubin_mean = var_rubin_mean,
              var_rubin_reg = var_rubin_reg,
              var_rubin_reg2 = var_rubin_reg2,
              var_boot_mean = var_boot_mean,
              var_boot_reg = var_boot_reg,
              var_boot_reg2 = var_boot_reg2))
}

## Data ----
N <- 100
k <- 2 # dimension of covariates (omit intercept)
p <- 5 # number of visits
# mu_beta_ctl <- c(0, 1, 2, 3, 4)
# mu_beta_trt <- c(0, 1.3, 2.8, 4, 5.5)
mu_beta_ctl <- c(0, 1, 2, 3, 4)
mu_beta_trt <- c(0, 1.3, 2.3, 3.5, 4.8)
set.seed(123)
beta_ctl <- rbind(rnorm(k+1, mu_beta_ctl[1], 1), rnorm(k+1,mu_beta_ctl[2], 1),
                  rnorm(k+1, mu_beta_ctl[3], 1), rnorm(k+1,mu_beta_ctl[4], 1),
                  rnorm(k+1, mu_beta_ctl[5], 1))
beta_trt <- rbind(rnorm(k+1, mu_beta_trt[1], 1), rnorm(k+1,mu_beta_trt[2], 1),
                  rnorm(k+1, mu_beta_trt[3], 1), rnorm(k+1,mu_beta_trt[4], 1),
                  rnorm(k+1, mu_beta_trt[5], 1))
beta_trt[1,] <- beta_ctl[1,]
## Under H0
# beta_trt <- beta_ctl

sd  <- c(2.0, 1.8, 2.0, 2.1, 2.2)
corr   <- matrix(
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
phi_ctl <- c(-3.5, 0.2, 0)
phi_trt <- c(-3.6, 0.2, 0)
# ## (iii) trt > ctl, approxiamte 0.7, 0.8
# phi_ctl <- c(-2.8, 0.2, 0)
# phi_trt <- c(-4.0, 0.2, 0)
## Under H0, ctl = trt, approx 0.8
# phi_ctl <- c(-3.5, 0.2, 0)
# phi_trt <- c(-3.5, 0.2, 0)

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
  true_mean <- c(5.526554, 6.607928, 1.081374) # MVT
  # true_mean <- c(5.526554, 5.526554, 0)
}
if(distribution == "MVGamma"){
  true_mean <- c(5.526218, 6.603291, 1.077073) # MVGamma
  # true_mean <- c(5.526218, 5.526218, 0)
}

## Analysis model ----
fit_model <- lm
# fit_model <- MASS::rlm
# fit_model <- Rfit::rfit

## Simulation results ----
M <- 10
B <- 100
# Get task id for each simulation job (not run in Rstudio Serve)
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
res <- main(seed = task_id, case = "Case0")

# Save all the objects into i.Rdata
save(res, file = paste0(task_id, ".Rdata"))
