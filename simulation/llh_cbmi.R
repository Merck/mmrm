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
    eps <- rmvt(n, delta = rep(0,p), sigma = sigma/3, df = 5)
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
j2r_llh <- function(db_comb, db_long){
  n1 <- length(which(db_comb$trt == 1))
  n2 <- length(which(db_comb$trt == 2))
  p <- length(unique(db_long$time))
  obs_pi <- 1 - c(sum(is.na(db_comb[,p] & db_comb$trt == 1))/n1,
                  sum(is.na(db_comb[,p] & db_comb$trt == 2))/n2)
  names_dep_var <- colnames(db_comb[, c(1:p)])
  names_covar <- colnames(db_comb[, c((p+1):(p+k))])
  db_comb_ctl <- db_comb[which(db_comb$trt == 1),]
  db_comb_trt <- db_comb[which(db_comb$trt == 2),]
  ## Compute relative change directly
  db_long <- db_long %>% group_by(id) %>%
    mutate(
      base = aval[time == 1],
      chg = aval - base) %>%
    ungroup()
  db_0 <- subset(db_long, time != 1)
  db_long_ctl <- db_0[which(db_0$trt == 1),]
  db_long_trt <- db_0[which(db_0$trt == 2),]
  db_avaiable <- na.omit(db_0)
  db_avaiable_ctl <- na.omit(db_long_ctl)
  db_avaiable_trt <- na.omit(db_long_trt)
  child_ctl <- factor(db_avaiable_ctl$id)
  child_trt <- factor(db_avaiable_trt$id)
  adj_mean <- colMeans(db_comb[,(p+1):(p+k)])
  
  # Step 1: for all observed data in each group, MMRM (under MAR)
  fit_lmer_ctl <- lmer(chg ~ x1 + x2 + factor(time) +
                         x1:factor(time) + x2:factor(time) +
                         (0 + factor(time)|id),
                       control = lmerControl(check.nobs.vs.nRE = "ignore",
                                             check.conv.grad = "ignore",
                                             check.conv.hess = "ignore",
                                             check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),
                       data = db_avaiable_ctl)
  fit_lmer_trt <- lmer(chg ~ x1 + x2 + factor(time) +
                         x1:factor(time) + x2:factor(time) +
                         (0 + factor(time)|id),
                       control = lmerControl(check.nobs.vs.nRE = "ignore",
                                             check.conv.grad = "ignore",
                                             check.conv.hess = "ignore",
                                             check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)),
                       data = db_avaiable_trt)
  beta_lmer_ctl <- fixef(fit_lmer_ctl)
  cov_beta_ctl <- vcov(fit_lmer_ctl, full = TRUE, ranpar = "var")
  cov_beta_ctl <- data.matrix(cov_beta_ctl)
  beta_transform_ctl <- c(beta_lmer_ctl[1:3],
                          beta_lmer_ctl[1:3] + c(beta_lmer_ctl[4], beta_lmer_ctl[7], beta_lmer_ctl[10]),
                          beta_lmer_ctl[1:3] + c(beta_lmer_ctl[5], beta_lmer_ctl[8], beta_lmer_ctl[11]),
                          beta_lmer_ctl[1:3] + c(beta_lmer_ctl[6], beta_lmer_ctl[9], beta_lmer_ctl[12]))
  beta_lmer_trt <- fixef(fit_lmer_trt)
  cov_beta_trt <- vcov(fit_lmer_trt, full = TRUE, ranpar = "var")
  cov_beta_trt <- data.matrix(cov_beta_trt)
  beta_transform_trt <- c(beta_lmer_trt[1:3],
                          beta_lmer_trt[1:3] + c(beta_lmer_trt[4], beta_lmer_trt[7], beta_lmer_trt[10]),
                          beta_lmer_trt[1:3] + c(beta_lmer_trt[5], beta_lmer_trt[8], beta_lmer_trt[11]),
                          beta_lmer_trt[1:3] + c(beta_lmer_trt[6], beta_lmer_trt[9], beta_lmer_trt[12]))
  
  ## Step 2: For completer in each group, MMRM -> ANCOVA at last time point
  db_complete_ctl <- db_comb_ctl[which(is.na(db_comb_ctl[,p]) == FALSE),]
  db_complete_trt <- db_comb_trt[which(is.na(db_comb_trt[,p]) == FALSE),]
  
  chg_ctl <- db_complete_ctl$y5 - db_complete_ctl$y1
  chg_trt <- db_complete_trt$y5 - db_complete_trt$y1
  fit_complete_ctl <- lm(chg_ctl ~ x1 + x2, data = db_complete_ctl)
  fit_complete_trt <- lm(chg_trt ~ x1 + x2, data = db_complete_trt)
  beta_complete_ctl <- coef(fit_complete_ctl)
  beta_complete_trt <- coef(fit_complete_trt)
  cov_complete_trt <- vcov(fit_complete_trt)
  
  # point estimate
  mu_complete_ctl <- beta_complete_ctl[1] + sum(beta_complete_ctl[2:3]*adj_mean)
  mu_mar_ctl <- beta_transform_ctl[10] + sum(beta_transform_ctl[11:12]*adj_mean)
  mu_ctl <- mu_mar_ctl
  mu_complete_trt <- beta_complete_trt[1] + sum(beta_complete_trt[2:3]*adj_mean)
  mu_trt <- obs_pi[2]*mu_complete_trt + (1 - obs_pi[2])*mu_mar_ctl
  llh_est <- c(mu_ctl, mu_trt, mu_trt - mu_ctl)
  
  # Variance estimate
  mat_vec <- diag(length(beta_lmer_ctl))
  mat_vec[1, (k+2):(k+p-1)] <- 1
  mat_vec[2, (k+p):(k+2*p-3)] <- 1
  mat_vec[3, (k+2*p-2):length(beta_lmer_ctl)] <- 1
  coef4_ctl <- rbind(mat_vec[,6], mat_vec[,9], mat_vec[,12])
  cov4_ctl <- coef4_ctl%*%cov_beta_ctl%*%t(coef4_ctl)
  cov4_trt <- cov_complete_trt
  var_xbar <- cov(db_comb[,(p+1):(p+k)])/(n1+n2)
  coef_mat <- matrix(c(1, adj_mean), nrow = 1)
  
  coef_beta_mat_ctl <- matrix(beta_transform_ctl[11:12], nrow = 1)
  coef_beta_mat_trt <- matrix(beta_complete_trt[2:3]*obs_pi[2] + beta_transform_ctl[11:12]*(1 - obs_pi[2]), nrow = 1)
  coef_beta_mat_diff <- matrix((beta_complete_trt[2:3] - beta_transform_ctl[11:12])*sqrt(obs_pi[2]^2 + obs_pi[2]*(1 - obs_pi[2])/n2), nrow = 1)
  var_theta_ctl <- as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat)) +
    as.numeric(coef_beta_mat_ctl%*%var_xbar%*%t(coef_beta_mat_ctl))
  var_theta_trt <- as.numeric(coef_mat%*%cov4_trt%*%t(coef_mat))*(obs_pi[2]^2 + obs_pi[2]*(1 - obs_pi[2])/n2) +
    as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat))*((1 - obs_pi[2])^2 + obs_pi[2]*(1 - obs_pi[2])/n2) +
    (mu_complete_trt - mu_mar_ctl)^2*obs_pi[2]*(1 - obs_pi[2])/n2 +
    as.numeric(coef_beta_mat_trt%*%var_xbar%*%t(coef_beta_mat_trt))
  var_theta_diff <- as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat))*(obs_pi[2]^2 + obs_pi[2]*(1 - obs_pi[2])/n2) +
    as.numeric(coef_mat%*%cov4_trt%*%t(coef_mat))*(obs_pi[2]^2 + obs_pi[2]*(1 - obs_pi[2])/n2) +
    (mu_complete_trt - mu_mar_ctl)^2*obs_pi[2]*(1 - obs_pi[2])/n2 +
    as.numeric(coef_beta_mat_diff%*%var_xbar%*%t(coef_beta_mat_diff))
  var_llh <- c(var_theta_ctl, var_theta_trt, var_theta_diff)
  
  return(list(llh_est = llh_est,
              var_llh = var_llh
              # obs_pi = obs_pi,
              # mut_ctl = mut_ctl,
              # mut_trt = mut_trt,
              # beta_transform_ctl = beta_transform_ctl,
              # beta_transform_trt = beta_transform_trt
  ))
}

## nonparametric bootstrap ----
nonpara_llh_fn <- function(db_comb, db_long, B){
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
    llh_res <- j2r_llh(data_boot, data_long)
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
      tmp1 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi_par = phi_ctl, trt = 1, missing_type = missing_type, distribution = distribution)
      tmp2 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_trt, pi_par = phi_trt, trt = 2, missing_type = missing_type, distribution = distribution)
    }
    if(case == "Case2"){
      tmp1 <- sim_one_group(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi_par = phi_ctl, trt = 1, missing_type = missing_type, distribution = distribution)
      tmp2 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_trt, pi_par = phi_trt, trt = 2, missing_type = missing_type, distribution = distribution)
    }
    if(case == "Case3"){
      tmp1 <- sim_one_extreme(n = N, k = k, P = p, sigma = Sigma, beta = beta_ctl, pi_par = phi_ctl, trt = 1, missing_type = missing_type, distribution = distribution)
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
  llh_res <- j2r_llh(db_comb, db_long)
  llh_boot <- nonpara_llh_fn(db_comb, db_long, B)
  
  return(list(llh_est = llh_res$llh_est,
              llh_var = llh_res$var_llh,
              llh_boot_var = llh_boot$var_llh_boot
  ))
}

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
# Under H0
beta_trt <- beta_ctl

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
# phi_ctl <- c(-3.5, 0.2, 0)
# phi_trt <- c(-3.6, 0.2, 0)
# ## (iii) trt > ctl, approxiamte 0.7, 0.8
# phi_ctl <- c(-2.8, 0.2, 0)
# phi_trt <- c(-4.0, 0.2, 0)
# Under H0, ctl = trt, approx 0.8
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
  true_mean <- c(5.526554, 6.607928, 1.081374) # MVT
  # true_mean <- c(5.526554, 5.526554, 0)
}
if(distribution == "MVGamma"){
  # true_mean <- c(5.526218, 6.603291, 1.077073) # MVGamma
  true_mean <- c(5.526218, 5.526218, 0)
}

## Simulation results ----
B <- 100
# Get task id for each simulation job (not run in Rstudio Serve)
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
res <- main(seed = task_id, case = "Case6")

# Save all the objects into i.Rdata
save(res, file = paste0(task_id, ".Rdata"))
