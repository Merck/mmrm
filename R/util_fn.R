# Copyright (c) 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
#
# This file is part of the mmrm program.
#
# mmrm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Simulate longitudinal data with monotone missingness for one group
#'
#' Generate longitudinal data with monotone missingness
#'
#' @param n sample size
#' @param k dimension of the covariates
#' @param P the number of time point
#' @param sigma Covariance matrix for the longitudinal response
#' @param beta regression coefficient for the baseline covariates
#' @param pi Observed probability for each missing pattern under MCAR
#' @param pi_par Parameters for the missing mechanism model under MAR
#' @param trt name of the treatment
#' @param missing_type Type of missing. Available: MCAR, MAR, MNAR
#' @param distribution Underlying multivariate distribution for the longitudinal responses
#' @param outlier The types of outlier may exists. Available: none, Type 1, Type 2
#'
#' @return a list including a wide-form data frame and a long-form data frame
#' \itemize{
#'   \item{db_comb}{Wide form of the original longitudinal data}
#'   \item{db_long}{Long form of the original longitudinal data}
#'   \item{db_comb0}{Wide form of the change from baseline longitudinal data}
#'   \item{db_long0}{Long form of the change from baseline longitudinal data}
#' }
#' @export
sim_one_group <- function(n, k, P, sigma, beta, pi = NULL, pi_par = NULL, trt = 1, missing_type, distribution, outlier = "none"){
  p <- P
  n_pattern <- length(pi)
  covar_x <- cbind(stats::rnorm(n), stats::rbinom(n, 1, 0.3))
  colnames(covar_x) <- paste0("x",1:k)
  xmat <- cbind(rep(1,n),covar_x)
  y_mat <- matrix(0, n, p)
  if(distribution == "MVN")
    eps <- mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = sigma)
  else if(distribution == "MVT")
    # eps <- rmvt(n, delta = rep(0,p), sigma = sigma, df = 3)
    eps <- mvtnorm::rmvt(n, delta = rep(0,p), sigma = sigma/5*3, df = 5)
  else if(distribution == "MVGamma"){
    u <- mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = stats::cov2cor(sigma))
    var_eps <- diag(sigma)
    eps <- matrix(0, n, p)
    for(j in 1:p){
      shape_para <- 2
      scale_para <- sqrt(var_eps[j]/shape_para)
      eps[,j] <- stats::qgamma(stats::pnorm(u[,j]), shape = shape_para, scale = scale_para) - shape_para*scale_para
    }
  }
  y_mean <- apply(beta, 1, function(x) as.vector(xmat%*%x))
  y_mat <- y_mean + eps
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
      .pattern <- ifelse(stats::rbinom(n = n, size = 1, prob = .pi) == 1, i_missing, .pattern)
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

  db_long <- stats::reshape(db_comb, varying = paste0("y", 1:p), v.names = "aval",
                            idvar = "id", times = c(1:p), timevar = "time",
                            direction = "long")
  rownames(db_long) <- 1:(n*p)

  # compute relative change
  rela_chg <- db_comb[, 2:p] - db_comb[,1]
  colnames(rela_chg) <- paste0("chg", 1:(p-1))
  db_comb0 <- cbind(db_comb, rela_chg)
  db_comb0 <- db_comb0[,-(2:p)]
  db_long0 <- stats::reshape(db_comb0, varying = paste0("chg", 1:(p-1)),
                             v.names = "chg",
                             idvar = "id", times = c(2:p), timevar = "time",
                             direction = "long")
  rownames(db_long0) <- 1:(n*(p-1))


  return(list(db_comb = db_comb, db_long = db_long,
              # change from baseline data
              db_comb0 = db_comb0, db_long0 = db_long0))
}

#' Transform the regression coefficients for lmer() object
#'
#' @param x regression coefficient returns from lmer()
#' @param k dimension of the covariates
#' @param p the number of time point
#'
#' @return transformed estimated coefficients at each time point
#' @export
transform_coef <- function(x, p, k){
  vec_bind <- c(x[1:(k+1)])
  for(j in 0:(p-2)){
    for(i in 0:k){
      vec_bind <- c(vec_bind, x[i+1] + x[k+2 + j +i*(p-1)])
    }
  }
  vec_bind
}

#' Transform the population covariance matrix for lmer() object
#'
#' @param x a data frame of the lower triangle elements obtained from lme4::VarCorr()
#' @param p the number of time point
#'
#' @return Estimated population covariance matrix
#' @export
transform_cov <- function(x, p){
  mat <- matrix(0, p, p)
  mat[c(lower.tri(mat, diag = TRUE))] <- x[-nrow(x),4]
  mat <- t(mat)
  mat[c(lower.tri(mat, diag = TRUE))] <- x[-nrow(x),4]
  mat
}


#' fit MMRM
#'
#' Fit MMRM for the long form data
#'
#' @param mmrm_formula formula for MMRM, including fixed and random effect
#' @param data_long long form of the data
#' @param p number of time points
#' @param k dimension of baseline covariates
#'
#' @return A list of model estimation results.
#' \itemize{
#'   \item{beta_lmer}{Estimated coefficients for the fixed effect for the raw
#'   formula.}
#'   \item{beta}{Estimated coefficients for the fixed effect for each time point.}
#'   \item{cov_beta}{Covariance matrix of the estimated cofficients for the
#'   fixed effect for the raw formula.}
#'   \item{cov4}{Estimated covariance matrix of the estimated cofficients at the
#'   last time point.}
#'   \item{V}{Estimated covariance matrix of the longitudinal response.}
#' }
#'
#' @export
mmrm_fit <- function(mmrm_formula, data_long, p, k){
  fit_lmer <- lme4::lmer(mmrm_formula,
                        control = lme4::lmerControl(check.nobs.vs.nRE = "ignore",
                                                    check.conv.grad = "ignore",
                                                    check.conv.hess = "ignore",
                                                    check.conv.singular = lme4::.makeCC(action = "ignore",  tol = 1e-4)),
                        data = data_long)
  beta_lmer <- lme4::fixef(fit_lmer)
  cov_beta <- stats::vcov(fit_lmer, full = TRUE, ranpar = "var")
  cov_beta <- data.matrix(cov_beta)
  beta_transform <- transform_coef(beta_lmer, p, k)

  coef4 <- cov_beta_last(beta_lmer, p, k)
  cov4 <- coef4%*%cov_beta%*%t(coef4)

  # random effect covariance matrix D for control group
  vc.a <- lme4::VarCorr(fit_lmer)
  vc.da <- as.data.frame(vc.a,order="lower.tri")
  var_D <- transform_cov(vc.da, p)
  var_R <- summary(fit_lmer)$sigma^2*diag(p)
  Z_i <- diag(p)
  V <- var_R + t(Z_i)%*%var_D%*%Z_i

  return(list(beta_lmer = beta_lmer,
              beta = beta_transform,
              cov_beta = cov_beta,
              cov4 = cov4,
              V = V))
}

#' MI-based Sequential regression
#'
#' Conduct sequential imputation under MAR to get parameter estimates for MMRM.
#'
#' @param seq_formula a list of formulas for sequential regression
#' @param data_wide wide form of the longitudinal data
#' @param trt_name name of treatment in the data frame
#' @param M imputation size
#' @param fit_model type of the analysis model. Available: lm, Rfit::rfit
#' (rank regression), MASS::rlm (robust regression)
#'
#' @return A list of estimation results from sequential regression:
#' \itemize{
#'   \item{obs_pi}{marginal probability of the observed data for each treatment.}
#'   \item{beta_ctl_imp}{M sets of the estimated MMRM coefficients used in MI
#'   for control group.}
#'   \item{beta_trt_imp}{M sets of the estimated MMRM coefficients used in MI
#'   for treatment group.}
#'   \item{sigma_ctl}{Estimated population covariance matrix of the longitudinal
#'   response.}
#'   \item{beta_rubin_ctl}{Estimated coefficients at the lsat time point for
#'   control group.}
#'   \item{beta_rubin_trt}{Estimated coefficients at the lsat time point for
#'   treatment group.}
#'   \item{var_rubin_ctl}{Estimated covariance matrix of the coefficients at
#'   the lsat time point for control group.}
#'   \item{var_rubin_trt}{Estimated covariance matrix of the coefficients at
#'   the lsat time point for treatment group.}
#' }
#' @export
seq_mi_fit <- function(seq_formula, data_wide, trt_name, M, fit_model){
  p <- length(seq_formula)
  covar_name <- all.vars(seq_formula[[1]])[-1]
  k <- length(covar_name)
  y_name <- sapply(seq_formula, function(x) all.vars(x)[1])
  y_name <- sort(y_name)
  base_name <- y_name[1]
  covar_name <- all.vars(seq_formula[[1]])[-1]
  last_name <- y_name[p]
  # trt_level <- unique(data_wide[,trt_name])
  data_wide_ctl <- data_wide[which(data_wide[,trt_name] == 1),]
  data_wide_trt <- data_wide[which(data_wide[,trt_name] == 2),]
  n1 <- nrow(data_wide_ctl)
  n2 <- nrow(data_wide_trt)
  obs_pi <- 1 - c(sum(is.na(data_wide_ctl[,last_name]))/n1,
                  sum(is.na(data_wide_trt[,last_name]))/n2)

  imp_fn <- function(){
    fit_time_ctl <- list()
    fit_time_trt <- list()
    beta_time_ctl <- list()
    beta_time_trt <- list()
    beta_cov_ctl <- list()
    beta_cov_trt <- list()
    sigma_time_ctl <- c(0)
    sigma_time_trt <- c(0)
    db_imp_ctl <- data_wide_ctl
    db_imp_trt <- data_wide_trt
    for(j in 1:p){
      # control
      fit_time_ctl[[j]] <- fit_model(seq_formula[[j]],
                                     data = db_imp_ctl)
      beta_time_ctl[[j]] <- fit_time_ctl[[j]]$coefficients
      beta_cov_ctl[[j]] <- stats::vcov(fit_time_ctl[[j]])

      # treatment
      fit_time_trt[[j]] <- fit_model(seq_formula[[j]],
                                     data = db_imp_trt)
      beta_time_trt[[j]] <- fit_time_trt[[j]]$coefficients
      beta_cov_trt[[j]] <- stats::vcov(fit_time_trt[[j]])

      if(!identical(fit_model, Rfit::rfit)){
        sigma_time_ctl[j] <- summary(fit_time_ctl[[j]])$sigma^2
        sigma_time_trt[j] <- summary(fit_time_trt[[j]])$sigma^2
        df_fit_ctl <- summary(fit_time_ctl[[j]])$df[1]
        df_fit_trt <- summary(fit_time_trt[[j]])$df[1]
        imp_mean_ctl <- stats::predict(fit_time_ctl[[j]], newdata = db_imp_ctl)[which(db_imp_ctl$pattern <= j + 1 & db_imp_ctl$pattern > 1)]
        imp_mean_trt <- stats::predict(fit_time_trt[[j]], newdata = db_imp_trt)[which(db_imp_trt$pattern <= j + 1 & db_imp_trt$pattern > 1)]
      }

      else{
        sigma_time_ctl[j] <- fit_time_ctl[[j]]$tauhat^2
        sigma_time_trt[j] <- fit_time_trt[[j]]$tauhat^2
        df_fit_ctl <- sum(!is.na(data_wide_ctl[,y_name[j]])) - (k+j)
        df_fit_trt <- sum(!is.na(data_wide_trt[,y_name[j]])) - (k+j)
        if(j == 1){
          design_ctl <- cbind(rep(1, n1), data.matrix(db_imp_ctl[,covar_name]))
          design_trt <- cbind(rep(1, n2), data.matrix(db_imp_trt[,covar_name]))
        }
        else{
          design_ctl <- cbind(rep(1, n1), data.matrix(db_imp_ctl[,covar_name]), data.matrix(db_imp_ctl[,y_name[1:(j-1)]]))
          design_trt <- cbind(rep(1, n2), data.matrix(db_imp_trt[,covar_name]), data.matrix(db_imp_trt[,y_name[1:(j-1)]]))
        }
        pred_ctl <- as.vector(design_ctl%*%beta_time_ctl[[j]])
        imp_mean_ctl <- pred_ctl[which(db_imp_ctl$pattern <= j+1 & db_imp_ctl$pattern > 1)]
        pred_trt <- as.vector(design_trt%*%beta_time_trt[[j]])
        imp_mean_trt <- pred_trt[which(db_imp_trt$pattern <= j+1 & db_imp_trt$pattern > 1)]
      }
      sigma_ini <- sigma_time_ctl[j]
      sigma_imp_ctl <- (df_fit_ctl-2)*sigma_ini/stats::rchisq(1, df = df_fit_ctl)
      imp_sd_ctl <- rep(sqrt(sigma_imp_ctl), length(which(db_imp_ctl$pattern <= j+1 & db_imp_ctl$pattern > 1)))
      db_imp_ctl[which(db_imp_ctl$pattern <= j+1 & db_imp_ctl$pattern > 1),y_name[j]] <- stats::rnorm(length(which(db_imp_ctl$pattern <= j+1 & db_imp_ctl$pattern > 1)), mean = imp_mean_ctl, sd = imp_sd_ctl)

      sigma_ini <- sigma_time_trt[j]
      sigma_imp_trt <- (df_fit_trt-2)*sigma_ini/stats::rchisq(1, df = df_fit_trt)
      imp_sd_trt <- rep(sqrt(sigma_imp_trt), length(which(db_imp_trt$pattern <= j+1 & db_imp_trt$pattern > 1)))
      db_imp_trt[which(db_imp_trt$pattern <= j+1 & db_imp_trt$pattern > 1),y_name[j]] <- stats::rnorm(length(which(db_imp_trt$pattern <= j+1 & db_imp_trt$pattern > 1)), mean = imp_mean_trt, sd = imp_sd_trt)

    }

    beta_imp_ctl <- list()
    beta_imp_trt <- list()
    sigma_imp_ctl <- c(0)
    control_term1 <- c(0)
    trt_term1 <- c(0)
    control_term2 <- matrix(0, k, p)
    trt_term2 <- matrix(0, k, p)
    for(j in 1:p){
      beta_imp_ctl[[j]] <- as.vector(mvtnorm::rmvnorm(1, mean = beta_time_ctl[[j]], sigma = beta_cov_ctl[[j]]))
      beta_imp_trt[[j]] <- as.vector(mvtnorm::rmvnorm(1, mean = beta_time_trt[[j]], sigma = beta_cov_trt[[j]]))

      if(j == 1){
        control_term1[j] <- beta_imp_ctl[[j]][1]
        control_term2[,j] <- beta_imp_ctl[[j]][2:(k+1)]
        trt_term1[j] <- beta_imp_trt[[j]][1]
        trt_term2[,j] <- beta_imp_trt[[j]][2:(k+1)]
      }
      else{
        sum_term1 <- 0
        sum_term2 <- rep(0,k)
        sum_term3 <- 0
        sum_term4 <- rep(0,k)
        for(i in (k+2):(k+j)){
          sum_term1 <- beta_imp_ctl[[j]][i]*control_term1[i-k-1] + sum_term1
          sum_term2 <- beta_imp_ctl[[j]][i]*control_term2[,i-k-1] + sum_term2

          sum_term3 <- beta_imp_trt[[j]][i]*trt_term1[i-k-1] + sum_term3
          sum_term4 <- beta_imp_trt[[j]][i]*trt_term2[,i-k-1] + sum_term4
        }
        control_term1[j] <- beta_imp_ctl[[j]][1] + sum_term1
        control_term2[,j] <- beta_imp_ctl[[j]][2:(k+1)] + sum_term2
        trt_term1[j] <- beta_imp_trt[[j]][1] + sum_term3
        trt_term2[,j] <- beta_imp_trt[[j]][2:(k+1)] + sum_term4
      }
    }
    sigma_coef_ctl <- diag(p)
    sigma_coef_trt <- diag(p)
    for(j in 1:p){
      if(j < p){
        for(ind in (j+1):p){
          sum_term <- 0
          sum_term2 <- 0
          for(i in j:(ind-1)){
            sum_term <- beta_time_ctl[[ind]][i+k+1]*sigma_coef_ctl[j,i] + sum_term
            sum_term2 <- beta_time_trt[[ind]][i+k+1]*sigma_coef_trt[j,i] + sum_term2
          }
          sigma_coef_ctl[j,ind] <- sum_term
          sigma_coef_trt[j,ind] <- sum_term2
        }
      }
    }
    component_mat_ctl <- apply(sigma_coef_ctl, 2, function(x) x^2*sigma_time_ctl)
    sigma_total_ctl <- apply(component_mat_ctl, 2, sum)
    component_mat_trt <- apply(sigma_coef_trt, 2, function(x) x^2*sigma_time_trt)
    sigma_total_trt <- apply(component_mat_trt, 2, sum)

    sigma_ctl <- diag(sigma_total_ctl)
    sigma_trt <- diag(sigma_total_trt)
    for(j1 in 1:(p-1)){
      for(j2 in (j1+1):p){
        for(ind in 1:(j2-1)){
          sigma_ctl[j1, j2] <- beta_time_ctl[[j2]][ind+1+k]*sigma_ctl[j1, ind] + sigma_ctl[j1, j2]
          sigma_trt[j1, j2] <- beta_time_trt[[j2]][ind+1+k]*sigma_trt[j1, ind] + sigma_trt[j1, j2]
        }
        sigma_ctl[j2, j1] <- sigma_ctl[j1, j2]
        sigma_trt[j2, j1] <- sigma_trt[j1, j2]
      }
    }
    # sigma_mi_ctl <- MCMCpack::rwish(p+1, sigma_ctl)/(p+1)

    beta_mi_ctl <- as.vector(rbind(control_term1, control_term2))
    beta_mi_trt <- as.vector(rbind(trt_term1, trt_term2))

    # get cov(hat beta)
    db_reg_ctl <- data_wide_ctl[,covar_name]
    db_reg_trt <- data_wide_trt[,covar_name]
    db_reg_ctl$wt <- db_imp_ctl[,last_name]
    db_reg_trt$wt <- db_imp_trt[,last_name]

    fit_ctl <- fit_model(stats::reformulate(covar_name, "wt"),
                         data = db_reg_ctl)
    beta_last_ctl <- stats::coef(fit_ctl)
    fit_trt <- fit_model(stats::reformulate(covar_name, "wt"),
                         data = db_reg_trt)
    beta_last_trt <- stats::coef(fit_trt)
    cov_beta_ctl <- stats::vcov(fit_ctl)
    cov_beta_trt <- stats::vcov(fit_trt)

    return(list(sigma_ctl = sigma_ctl,
                sigma_trt = sigma_trt,
                beta_mi_ctl = beta_mi_ctl,
                beta_mi_trt = beta_mi_trt,
                beta_last_ctl = beta_last_ctl,
                beta_last_trt = beta_last_trt,
                cov_last_ctl = cov_beta_ctl,
                cov_last_trt = cov_beta_trt))
  }
  # sigma_ctl_imp <- array(0, dim = c(p, p, M))
  beta_ctl_imp <- matrix(0, (k+1)*p, M)
  beta_trt_imp <- matrix(0, (k+1)*p, M)
  beta_last_ctl <- matrix(0, k+1, M)
  beta_last_trt <- matrix(0, k+1, M)
  cov_last_ctl <- array(0, dim = c(k+1, k+1, M))
  cov_last_trt <- array(0, dim = c(k+1, k+1, M))
  for(m in 1:M){
    imp_res <- imp_fn()
    beta_ctl_imp[,m] <- imp_res$beta_mi_ctl
    beta_trt_imp[,m] <- imp_res$beta_mi_trt
    beta_last_ctl[,m] <- imp_res$beta_last_ctl
    beta_last_trt[,m] <- imp_res$beta_last_trt
    cov_last_ctl[,,m] <- imp_res$cov_last_ctl
    cov_last_trt[,,m] <- imp_res$cov_last_trt
  }
  sigma_ctl <- imp_res$sigma_ctl
  sigma_trt <- imp_res$sigma_trt
  beta_rubin_ctl <- rowMeans(beta_last_ctl)
  beta_rubin_trt <- rowMeans(beta_last_trt)

  wm_ctl <- apply(cov_last_ctl, c(1,2), mean)
  bm_ctl <- stats::cov(t(beta_last_ctl))
  var_rubin_ctl <- wm_ctl + (1+1/M)*bm_ctl
  wm_trt <- apply(cov_last_trt, c(1,2), mean)
  bm_trt <- stats::cov(t(beta_last_trt))
  var_rubin_trt <- wm_trt + (1+1/M)*bm_trt

  return(list(obs_pi = obs_pi,
              beta_ctl_imp = beta_ctl_imp,
              beta_trt_imp = beta_trt_imp,
              sigma_ctl = sigma_ctl,
              sigma_trt = sigma_trt,
              # last time point
              beta_rubin_ctl = beta_rubin_ctl,
              beta_rubin_trt = beta_rubin_trt,
              var_rubin_ctl = var_rubin_ctl,
              var_rubin_trt = var_rubin_trt))
}

#' J2R imputation
#'
#' Conduct J2R imputation, using the estimated coefficients from MMRM
#'
#' @param data_wide wide form of the data
#' @param beta_imp_ctl M sets of the estimated MMRM coefficients used in MI
#'   for control group
#' @param beta_imp_trt M sets of the estimated MMRM coefficients used in MI
#'   for treatment group
#' @param sigma_ctl Estimated population covariance matrix of the longitudinal
#'   response for the control group
#' @param sigma_trt Estimated population covariance matrix of the longitudinal
#'   response for the treatment group
#' @param covar_name name of the baseline covariates
#' @param y_name name of the longitudinal response
#' @param trt_name name of the treatment
#'
#' @return imputation set at last time point
#' @export
j2r_imp <- function(data_wide, beta_imp_ctl, beta_imp_trt, sigma_ctl, sigma_trt,
                    covar_name, y_name, trt_name){
  p <- nrow(sigma_ctl)
  k <- nrow(beta_imp_ctl)/p - 1
  data_wide_ctl <- data_wide[which(data_wide[,trt_name] == 1),]
  data_wide_trt <- data_wide[which(data_wide[,trt_name] == 2),]
  n1 <- nrow(data_wide_ctl)
  n2 <- nrow(data_wide_trt)
  intermis_ctl <- which(data_wide_ctl$pattern == 0)
  intermis_trt <- which(data_wide_trt$pattern == 0)
  if((length(intermis_ctl) > 0) | (length(intermis_trt) > 0) ){
    print("Intermittent missing pattern is detected in the data,
              impute intermittent missing data under MAR")
  }
  imp_fn <- function(beta){
    beta_imp_ctl <- matrix(beta[1:(length(beta)/2)], 1+k, p)
    beta_imp_trt <- matrix(beta[(length(beta)/2 + 1):length(beta)], 1+k, p)
    V_mi_ctl <- MCMCpack::rwish(p+1, sigma_ctl)/(p+1)
    V_mi_trt <- MCMCpack::rwish(p+1, sigma_trt)/(p+1)

    control_term1 <- as.vector(beta_imp_ctl[1,])
    control_term2 <- beta_imp_ctl[-1,]
    trt_term1 <- as.vector(beta_imp_trt[1,])
    trt_term2 <- beta_imp_trt[-1,]
    imp_value_ctl <- data_wide_ctl[,y_name[p]]
    imp_value_trt <- data_wide_trt[,y_name[p]]
    sigma_pattern <- c(0)
    mean_ctl <- rep(0, n1)
    mean_trt <- rep(0, n2)
    # intermittent missing

    if(length(intermis_ctl) > 0){
      index_ctl <- 1:p
      for(i in 1:length(intermis_ctl)){
        missingInd <- which(is.na(data_wide_ctl[intermis_ctl[i], y_name]))
        obsInd <- index_ctl[-missingInd]
        sigma_mo <- V_mi_ctl[missingInd, obsInd]%*%
          solve(V_mi_ctl[obsInd, obsInd])
        V_mat <- V_mi_ctl[missingInd, missingInd] - sigma_mo%*%
          V_mi_ctl[obsInd,missingInd]
        sigma_pattern <- diag(V_mat)[length(diag(V_mat))]
        design_mat <- data_wide_ctl[intermis_ctl[i], covar_name]
        mu_miss_ctl <- control_term1[missingInd] + apply(matrix(control_term2[,missingInd], k, length(missingInd)), 2, function(y) sum(design_mat*y))
        mu_obs_ctl <- control_term1[obsInd] + apply(matrix(control_term2[,obsInd], k, length(obsInd)), 2, function(y) sum(design_mat*y))

        mean_pattern_ctl <- mu_miss_ctl + sigma_mo%*%(matrix(data.matrix(data_wide_ctl[intermis_ctl[i],y_name[obsInd]] - mu_obs_ctl), length(obsInd), 1))
        mean_ctl <- mean_pattern_ctl[length(mean_pattern_ctl)]
        imp_value_ctl[intermis_ctl[i]] <- stats::rnorm(1, mean_ctl, sqrt(sigma_pattern))
      }
    }
    if(length(intermis_trt) > 0){
      index_trt <- 1:p
      for(i in 1:length(intermis_trt)){
        missingInd <- which(is.na(data_wide_trt[intermis_trt[i], y_name]))
        obsInd <- index_trt[-missingInd]
        sigma_mo <- V_mi_trt[missingInd, obsInd]%*%
          solve(V_mi_trt[obsInd, obsInd])
        V_mat <- V_mi_trt[missingInd, missingInd] - sigma_mo%*%
          V_mi_trt[obsInd,missingInd]
        sigma_pattern <- diag(V_mat)[length(diag(V_mat))]
        design_mat <- data_wide_trt[intermis_trt[i], covar_name]
        mu_miss_trt <- trt_term1[missingInd] + apply(matrix(trt_term2[,missingInd], k, length(missingInd)), 2, function(y) sum(design_mat*y))
        mu_obs_trt <- trt_term1[obsInd] + apply(matrix(trt_term2[,obsInd], k, length(obsInd)), 2, function(y) sum(design_mat*y))

        mean_pattern_trt <- mu_miss_trt + sigma_mo%*%(matrix(data.matrix(data_wide_trt[intermis_trt[i],y_name[obsInd]] - mu_obs_trt), length(obsInd), 1))
        mean_trt <- mean_pattern_trt[length(mean_pattern_trt)]
        imp_value_trt[intermis_trt[i]] <- stats::rnorm(1, mean_trt, sqrt(sigma_pattern))
      }
    }

    # monotone missing
    for(j in 2:(p+1)){
      if(length(which(data_wide_ctl$pattern == j)) > 0){
        pattern_mat_ctl <- data_wide_ctl[which(data_wide_ctl$pattern == j), ]
        mu_miss_ctl <- control_term1[(j-1):p] + apply(pattern_mat_ctl[,covar_name], 1, function(x) apply(matrix(control_term2[,(j-1):p], nrow = k), 2, function(y) sum(x*y)))

        if(j == 2){
          sigma_pattern[j-1] <- V_mi_ctl[p,p]
          mean_pattern_ctl <- mu_miss_ctl
        }
        if(j > 2){
          sigma_mo <- V_mi_ctl[(j-1):p, 1:(j-2)]%*%solve(V_mi_ctl[1:(j-2), 1:(j-2)])
          sigma_pattern[j-1] <- (V_mi_ctl[(j-1):p, (j-1):p] - sigma_mo%*%V_mi_ctl[1:(j-2),(j-1):p])[p-j+2, p-j+2]

          mu_obs_ctl <- control_term1[1:(j-2)] + apply(pattern_mat_ctl[,covar_name], 1, function(x) apply(matrix(control_term2[,1:(j-2)], nrow = k), 2, function(y) sum(x*y)))
          mean_pattern_ctl <- mu_miss_ctl + apply(matrix(data.matrix(pattern_mat_ctl[,y_name[1:(j-2)]] - t(mu_obs_ctl)), nrow = nrow(pattern_mat_ctl)), 1, function(x) sigma_mo%*%x)
        }
        mean_pattern_ctl <- matrix(mean_pattern_ctl, ncol = nrow(pattern_mat_ctl))
        mean_ctl[which(data_wide_ctl$pattern == j)] <- mean_pattern_ctl[nrow(mean_pattern_ctl),]

        imp_value_ctl[which(data_wide_ctl$pattern == j)] <- sapply(mean_ctl[which(data_wide_ctl$pattern == j)], function(x) stats::rnorm(1, x, sqrt(sigma_pattern[j-1])))
      }
      if(length(which(data_wide_trt$pattern == j)) > 0){
        pattern_mat_trt <- data_wide_trt[which(data_wide_trt$pattern == j), ]
        mu_miss_trt <- control_term1[(j-1):p] + apply(pattern_mat_trt[,covar_name], 1, function(x) apply(matrix(control_term2[,(j-1):p], nrow = k), 2, function(y) sum(x*y)))

        if(j == 2){
          mean_pattern_trt <- mu_miss_trt
        }
        if(j > 2){
          sigma_mo <- V_mi_ctl[(j-1):p, 1:(j-2)]%*%solve(V_mi_ctl[1:(j-2), 1:(j-2)])
          sigma_pattern[j-1] <- (V_mi_ctl[(j-1):p, (j-1):p] - sigma_mo%*%V_mi_ctl[1:(j-2),(j-1):p])[p-j+2, p-j+2]
          mu_obs_trt <- trt_term1[1:(j-2)] + apply(pattern_mat_trt[,covar_name], 1, function(x) apply(matrix(trt_term2[,1:(j-2)], nrow = k), 2, function(y) sum(x*y)))
          mean_pattern_trt <- mu_miss_trt + apply(matrix(data.matrix(pattern_mat_trt[,y_name[1:(j-2)]] - t(mu_obs_trt)), nrow = nrow(pattern_mat_trt)), 1, function(x) sigma_mo%*%x)
        }
        mean_pattern_trt <- matrix(mean_pattern_trt, ncol = nrow(pattern_mat_trt))
        mean_trt[which(data_wide_trt$pattern == j)] <- mean_pattern_trt[nrow(mean_pattern_trt),]
        imp_value_trt[which(data_wide_trt$pattern == j)] <- sapply(mean_trt[which(data_wide_trt$pattern == j)], function(x) stats::rnorm(1, x, sqrt(sigma_pattern[j-1])))
      }
    }
    imp_value <- c(imp_value_ctl, imp_value_trt)

    return(imp_value)
  }
  imp_value <- apply(rbind(beta_imp_ctl, beta_imp_trt), 2, imp_fn)

  # # change from baseline
  # base_value <- c(data_wide_ctl[,1], data_wide_trt[,1])
  # chg_imp <- t(apply(imp_value, 2, function(x) x - base_value))

  return(imp_value)
}

#' Summary of long-form data frame given sequential regression formula
#'
#' Given a long-form data and a formula from MMRM, summarize the information
#' from the data
#'
#' @param mmrm_formula formula for fixed effect in MMRM
#' @param data_long long form of the data
#' @param time_name name of the time variable
#' @param id_name name of the id variable
#' @param trt_name name of the treatment variable
#'
#' @return A list of summary of the long form data:
#' \itemize{
#'   \item{p} number of time points
#'   \item{covar_name}{name of the covariates}
#'   \item{y_name}{name of the response}
#'   \item{k}{dimension of the covariates}
#'   \item{adj_mean}{overall covariate mean}
#'   \item{var_xbar}{covariance matrix of the covariates}
#'   \item{obs_pi}{marginal probability of the observed data for each treatment}
#'   \item{db_avaiable_ctl}{long-form available data in the control group}
#'   \item{db_avaiable_trt}{long-form available data in the treatment group}
#'   \item{data_wide}{wide-form data}
#'   \item{data_wide_ctl}{wide-form data in the control group}
#'   \item{data_wide_trt}{wide-form data in the treatment group}
#'   \item{n1}{number of subjects in the control group}
#'   \item{n2}{number of subjects in the treatment group}
#'   \item{data_complete_ctl}{complete data in the control group}
#'   \item{data_complete_trt}{complete data in the treatment group}
#' }
#' @export
summaryLong <- function(mmrm_formula, data_long, time_name, id_name, trt_name){
  p <- length(unique(data.matrix(data_long[,time_name])))
  var_name <- all.vars(mmrm_formula)
  y_name <- var_name[1]
  ind_time <- which(var_name == time_name)
  covar_name <- all.vars(mmrm_formula)[-c(1, ind_time)]
  k <- length(covar_name)

  data_long_ctl <- data_long[which(data_long[,trt_name] == 1),]
  data_long_trt <- data_long[which(data_long[,trt_name] == 2),]
  db_avaiable <- stats::na.omit(data_long)
  db_avaiable_ctl <- stats::na.omit(data_long_ctl)
  db_avaiable_trt <- stats::na.omit(data_long_trt)

  data_wide <- stats::reshape(data_long, v.names= y_name, idvar=id_name,
                        timevar=time_name, direction="wide")
  data_wide_ctl <- data_wide[which(data_wide[,trt_name] == 1),]
  data_wide_trt <- data_wide[which(data_wide[,trt_name] == 2),]
  n1 <- nrow(data_wide_ctl)
  n2 <- nrow(data_wide_trt)
  data_complete_ctl <- stats::na.omit(data_wide_ctl)
  data_complete_trt <- stats::na.omit(data_wide_trt)

  y_name <- paste0(y_name, ".", as.vector(unique(data.matrix(data_long[,time_name]))))
  y_name <- sort(y_name)
  adj_mean <- colMeans(data_wide[,covar_name])
  var_xbar <- stats::cov(data_wide[,covar_name])/nrow(data_wide)
  obs_pi <- 1 - c(sum(is.na(data_wide[,y_name[p]] & data_wide[,trt_name] == 1))/n1,
                  sum(is.na(data_wide[,y_name[p]] & data_wide[,trt_name] == 2))/n2)

  return(list(p = p,
              covar_name = covar_name,
              y_name = y_name,
              k = k,
              adj_mean = adj_mean,
              var_xbar = var_xbar,
              obs_pi = obs_pi,
              db_avaiable_ctl = db_avaiable_ctl,
              db_avaiable_trt = db_avaiable_trt,
              data_wide = data_wide,
              data_wide_ctl = data_wide_ctl,
              data_wide_trt = data_wide_trt,
              n1 = n1, n2 = n2,
              data_complete_ctl = data_complete_ctl,
              data_complete_trt = data_complete_trt))
}

#' Create formula that fits the MMRM model (fixed + random effect)
#'
#' @param mmrm_formula formula for fixed effect in MMRM
#' @param time_name name of the time variable
#' @param id_name name of the id variable
#'
#' @export
mmrmformula <- function(mmrm_formula, time_name, id_name){
  random.part <- paste0("+ (0 + factor(", time_name, ")|", id_name, ")")
  string_temp <- format(mmrm_formula)
  fix.part <- NULL
  for(i in 1:length(string_temp)){
    fix.part <- paste(fix.part,string_temp[i])
  }
  new_formula <- stats::formula(paste(fix.part, random.part, collapse = " + "))

  new_formula
}


#' Point estimation using naive method, under J2R
#'
#' Obtain estimation of the ATE using naive method by Rubin's rule
#'
#' @param mmrm_formula a symbolic description of the model to be fitted
#' @param data_long long-form of the data frame
#' @param time_name name of the time variable
#' @param id_name name of the id variable
#' @param trt_name name of the treatment variable
#' @param M imputation size
#' @param fit_model type of the analysis model. Available: lm, Rfit::rfit
#' (rank regression), MASS::rlm (robust regression)
#' @param estimator Estimator type in Rubin's method
#'
#' @return Estimation results using naive method.
#' \itemize{
#'   \item{ate_est}{MI point estimation of the ATE}
#'   \item{var_rubin}{MI variance estimation by Rubin's rule}
#' }
#' @export
naive_est <- function(mmrm_formula, data_long, time_name, id_name, trt_name,
                      M = 10, fit_model, estimator = "reg with interaction"){
  info_data <- summaryLong(mmrm_formula, data_long, time_name, id_name, trt_name)
  p <- info_data$p
  k <- info_data$k
  y_name <- sort(info_data$y_name)
  covar_name <- info_data$covar_name
  db_avaiable_ctl <- info_data$db_avaiable_ctl
  db_avaiable_trt <- info_data$db_avaiable_trt
  data_wide <- info_data$data_wide

  new_mmrm_formula <- mmrmformula(mmrm_formula, time_name, id_name)

  fit_mmrm_ctl <- mmrm_fit(new_mmrm_formula, db_avaiable_ctl, p, k)
  fit_mmrm_trt <- mmrm_fit(new_mmrm_formula, db_avaiable_trt, p, k)
  beta_lmer_ctl <- fit_mmrm_ctl$beta_lmer
  cov_beta_ctl <- fit_mmrm_ctl$cov_beta
  # draw beta
  beta_mi_ctl <- mvtnorm::rmvnorm(M, mean = beta_lmer_ctl, sigma = cov_beta_ctl)
  beta_transform_ctl <- apply(beta_mi_ctl, 1, function(x) transform_coef(x, p, k))
  beta_lmer_trt <- fit_mmrm_trt$beta_lmer
  cov_beta_trt <- fit_mmrm_trt$cov_beta
  # draw beta
  beta_mi_trt <- mvtnorm::rmvnorm(M, mean = beta_lmer_trt, sigma = cov_beta_trt)
  beta_transform_trt <- apply(beta_mi_trt, 1, function(x) transform_coef(x, p, k))
  V_ctl <- fit_mmrm_ctl$V
  V_trt <- fit_mmrm_trt$V

  imp_value <- j2r_imp(data_wide,
                       beta_transform_ctl, beta_transform_trt,
                       V_ctl, V_trt,
                       covar_name, y_name, trt_name)

  mi_res <- rubin_est(t(imp_value), data_wide = data_wide,
                      covar_name = covar_name, trt_name = trt_name,
                      M = M, fit_model = fit_model, estimator = estimator)

  return(list(ate_est = mi_res$ate_est,
              var_rubin = mi_res$var_rubin))
}

#' Inference using naive method, under J2R
#'
#' Obtain estimation of the ATE using naive method by Rubin's rule or bootstrap
#'
#' @param mmrm_formula a symbolic description of the model to be fitted
#' @param data_long long-form of the data frame
#' @param time_name name of the time variable
#' @param id_name name of the id variable
#' @param trt_name name of the treatment variable
#' @param M imputation size
#' @param fit_model type of the analysis model. Available: lm, Rfit::rfit
#' (rank regression), MASS::rlm (robust regression)
#' @param estimator Estimator type in Rubin's method
#' @param bootstrap Whether to use bootstrap variance
#' @param B the number of bootstrap replicate
#'
#' @return point estimation, variance estimation (Rubin + bootstrap)
#' \itemize{
#'   \item{ate_est}{MI point estimation of the ATE}
#'   \item{var_rubin}{MI variance estimation by Rubin's rule}
#'   \item{var_boot}{MI variance estimation by bootstrap}
#' }
#' @export
naiveMI <- function(mmrm_formula, data_long, time_name, id_name, trt_name,
                    M = 10, fit_model, estimator = "reg with interaction",
                    bootstrap = FALSE, B = 100){

  naive_res <- naive_est(mmrm_formula, data_long, time_name, id_name, trt_name,
                         M = M, fit_model = fit_model, estimator = estimator)
  if(bootstrap){
    var_boot <- nonpara_fn(formula = mmrm_formula, data_long = data_long,
                           time_name = time_name, id_name = id_name,
                           trt_name = trt_name, M = M, B = B,
                           fit_model = fit_model,
                           estimator = estimator,
                           method = "naive MI")
  }
  return(list(ate_est = naive_res$ate_est,
              var_rubin = naive_res$var_rubin,
              var_boot = var_boot))
}


#' Get Rubin's estimate
#'
#' Get estimation results by Rubin's rule
#'
#' @param chg_imp a matrix of imputed response at the last time point
#' @param data_wide wide-form of the data frame
#' @param covar_name names of the baseline covariates
#' @param trt_name name of the treatment variable
#' @param M imputation size
#' @param fit_model type of the analysis model. Available: lm, Rfit::rfit
#' (rank regressoin), MASS::rlm (robust regression)
#' @param estimator type of the estimator. Available: sample mean,
#' reg with interaction, reg without interaction. Default: reg with interaction.
#'
#' @return Rubin's point and variance estimate
#' \itemize{
#'   \item{ate_est}{MI point estimation of the ATE}
#'   \item{var_rubin}{MI variance estimation by Rubin's rule}
#' }
#' @export
rubin_est <- function(chg_imp, data_wide, covar_name,
                      trt_name, M, fit_model, estimator){
  k <- length(covar_name)
  data_wide_ctl <- data_wide[which(data_wide[,trt_name] == 1),]
  data_wide_trt <- data_wide[which(data_wide[,trt_name] == 2),]
  n1 <- nrow(data_wide_ctl)
  n2 <- nrow(data_wide_trt)
  chg_imp_ctl <- chg_imp[,1:n1]
  chg_imp_trt <- chg_imp[,(n1+1):(n1+n2)]
  adj_mean <- colMeans(data_wide[,covar_name])
  var_xbar <- stats::cov(data_wide[,covar_name])/(n1+n2)

  ## point estimates
  # (1) mean estimator
  # (i) simple average
  if(estimator == "sample mean"){
    theta1 <- mean(chg_imp_ctl)
    theta2 <- mean(chg_imp_trt)
    theta_diff <-  theta2 - theta1
    ate_est <- c(theta1, theta2, theta_diff)

    # Rubin's estimate
    wm_ctl <- mean(apply(chg_imp_ctl, 1, stats::var)/n1)
    bm_ctl <- stats::var(apply(chg_imp_ctl, 1, mean))
    var_rubin_ctl <- wm_ctl + (1+1/M)*bm_ctl
    wm_trt <- mean(apply(chg_imp_trt, 1, stats::var)/n2)
    bm_trt <- stats::var(apply(chg_imp_trt, 1, mean))
    var_rubin_trt <- wm_trt + (1+1/M)*bm_trt
    wm_diff <- wm_ctl + wm_trt
    bm_diff <- stats::var(apply(chg_imp_trt, 1, mean) -
                            apply(chg_imp_ctl, 1, mean))
    var_rubin_diff <- wm_diff + (1+1/M)*bm_diff
    var_rubin <- c(var_rubin_ctl, var_rubin_trt, var_rubin_diff)
  }

  if(estimator == "reg with interaction"){
    # (ii) robust regression for each group
    db_reg_ctl <- data_wide_ctl[,covar_name]
    db_reg_trt <- data_wide_trt[,covar_name]
    coef1_mat <- matrix(c(1, adj_mean), 1, k+1)
    coef2_mat <- matrix(c(1, adj_mean), 1, k+1)

    rr_fn <- function(y){
      z1 <- y[1:n1]
      z2 <- y[(n1+1):(n1+n2)]
      db_reg_ctl$wt <- z1
      db_reg_trt$wt <- z2

      fit_ctl <- fit_model(stats::reformulate(covar_name, "wt"),
                           data = db_reg_ctl)
      beta_ctl <- stats::coef(fit_ctl)
      fit_trt <- fit_model(stats::reformulate(covar_name, "wt"),
                           data = db_reg_trt)
      beta_trt <- stats::coef(fit_trt)
      theta1 <- beta_ctl[1] + sum(beta_ctl[2:(k+1)]*adj_mean)
      theta2 <- beta_trt[1] + sum(beta_trt[2:(k+1)]*adj_mean)
      theta_diff <- theta2 - theta1
      mean_est <- c(theta1, theta2, theta_diff)
      cov_beta_ctl <- stats::vcov(fit_ctl)
      cov_beta_trt <- stats::vcov(fit_trt)
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
    ate_est <- rowMeans(mi_reg_mat)[1:3]

    # Rubin's estimate
    wm_ctl <- rowMeans(mi_reg_mat)[4]
    bm_ctl <- stats::var(mi_reg_mat[1,])
    var_rubin_ctl <- wm_ctl + (1+1/M)*bm_ctl
    wm_trt <- rowMeans(mi_reg_mat)[5]
    bm_trt <- stats::var(mi_reg_mat[2,])
    var_rubin_trt <- wm_trt + (1+1/M)*bm_trt
    wm_diff <- rowMeans(mi_reg_mat)[6]
    bm_diff <- stats::var(mi_reg_mat[3,])
    var_rubin_diff <- wm_diff + (1+1/M)*bm_diff
    var_rubin <- c(var_rubin_ctl, var_rubin_trt, var_rubin_diff)
  }


  # (iii) robust regression for the whole data
  if(estimator == "reg without interaction"){
    db_reg <- data_wide[,c(covar_name, trt_name)]
    coef1_mat <- matrix(c(1, adj_mean, 0), 1, k+2)
    coef2_mat <- matrix(c(1, adj_mean, 1), 1, k+2)
    coef3_mat <- coef2_mat - coef1_mat
    rr2_fn <- function(y){
      z <- y
      db_reg$wt <- z

      fit_rr <- fit_model(stats::reformulate(c(covar_name, trt_name), "wt"), data = db_reg)
      beta_rr <- stats::coef(fit_rr)
      theta1 <- beta_rr[1] + sum(beta_rr[2:(k+1)]*adj_mean)
      theta2 <- beta_rr[1] + sum(beta_rr[2:(k+1)]*adj_mean) + beta_rr[k+2]
      theta_diff <- theta2 - theta1
      mean_est <- c(theta1, theta2, theta_diff)
      cov_beta_rr <- stats::vcov(fit_rr)
      var_theta1 <- as.numeric(coef1_mat%*%cov_beta_rr%*%t(coef1_mat)) +
        sum(beta_rr[2:(k+1)]^2*var_xbar)
      var_theta2 <- as.numeric(coef2_mat%*%cov_beta_rr%*%t(coef2_mat)) +
        sum(beta_rr[2:(k+1)]^2*var_xbar)
      var_diff <- as.numeric(coef3_mat%*%cov_beta_rr%*%t(coef3_mat))
      var_est <- c(var_theta1, var_theta2, var_diff)
      return(c(mean_est, var_est))
    }
    mi_reg2_mat <- apply(chg_imp, 1, rr2_fn)
    ate_est <- rowMeans(mi_reg2_mat)[1:3]

    # Rubin's estimate
    wm_ctl <- rowMeans(mi_reg2_mat)[4]
    bm_ctl <- stats::var(mi_reg2_mat[1,])
    var_rubin_ctl <- wm_ctl + (1+1/M)*bm_ctl
    wm_trt <- rowMeans(mi_reg2_mat)[5]
    bm_trt <- stats::var(mi_reg2_mat[2,])
    var_rubin_trt <- wm_trt + (1+1/M)*bm_trt
    wm_diff <- rowMeans(mi_reg2_mat)[6]
    bm_diff <- stats::var(mi_reg2_mat[3,])
    var_rubin_diff <- wm_diff + (1+1/M)*bm_diff
    var_rubin <- c(var_rubin_ctl, var_rubin_trt, var_rubin_diff)
  }

  return(list(ate_est = ate_est,
              var_rubin = var_rubin))
}

#' Get bootstrap variance estimate
#'
#' For a given dataset, compute the bootstrap variance estimate by a
#' prespecified method
#'
#' @param formula a symbolic description of the model to be fitted
#' @param data_wide wide-form of the data frame
#' @param data_long long-form of the data frame
#' @param time_name name of the time variable (needed for long-form data)
#' @param id_name name of the id variable
#' @param trt_name name of the treatment variable
#' @param M imputation size
#' @param B number of bootstrap replicates
#' @param fit_model type of the analysis model. Available: lm, Rfit::rfit
#' (rank regressoin), MASS::rlm (robust regression)
#' @param estimator type of the estimator. Available: sample mean,
#' reg with interaction, reg without interaction
#' @param method several methods to use. Available: naive MI, sequential MI,
#' naive likelihood, sequential likelihood
#' @param llh_method name of the likelihood-based method. Available: PMM, CBMI
#'
#' @return variance estimation by nonparametric bootstrap
#' @export
nonpara_fn <- function(formula, data_wide = NULL, data_long = NULL,
                       time_name = NULL, id_name, trt_name,
                       M, B, fit_model,
                       estimator,
                       method, llh_method = NULL){
  # if(!method %in% c("naive MI", "sequential MI",
  #                   "naive likelihood", "sequential likelihood"))
  #   stop("Invalid method. Available methods: naive MI, sequential MI,
  #        naive likelihood, sequential likelihood")
  # if(is.null(data_wide) & is.null(data_long))
  #   stop("No available data")
  # if(is.null(id_name))
  #   stop("Need to input id name")
  ate_boot <- matrix(0, 3, B)
  for(b in 1:B){
    set.seed(b)
    if(method == "naive MI"){
      # if(is.null(data_long))
      #   stop("Missing valid data. Need to input long form of the data")
      # if(is.null(time_name) | is.null(id_name))
      #   stop("Missing time and id variable name. For long form data,
      #        need to input time_name and id_name.")
      id_boot <- sample(unique(data_long[,id_name]), replace = TRUE)
      data_boot_long <- data.frame(id_new = 1:length(id_boot), id = id_boot)
      data_boot_long <- merge(data_boot_long, data_long, all.x = TRUE)
      colnames(data_boot_long) <- c("id_old", id_name,
                                    colnames(data_boot_long)[3:ncol(data_boot_long)])
      naive_res <- naive_est(formula, data_boot_long, time_name, id_name,
                             trt_name, M, fit_model = fit_model, estimator)
      ate_boot[,b] <- naive_res$ate_est
    }
    else if(method == "sequential MI"){
      id_boot <- sample(unique(data_wide[,id_name]), replace = TRUE)
      data_boot <- data.frame(id_new = 1:length(id_boot), id = id_boot)
      data_boot <- merge(data_boot, data_wide, all.x = TRUE)
      data_boot <- cbind(data_boot[,-c(1:2)], data_boot[,1:2])
      colnames(data_boot) <- c(colnames(data_boot)[1:(ncol(data_boot)-2)],
                               "id_old", "id")
      seq_mi_res <- rr_seq_imp(formula, data_boot, trt_name, M, fit_model,
                               estimator = estimator)
      ate_boot[,b] <- seq_mi_res$ate_est
    }
    else if(method == "naive likelihood"){
      # if(is.null(llh_method))
      #   stop("Missing specific likelihood method to use. Need to input
      #        likelihood-based method in `llh_method`. Available choice:
      #        PMM or CBMI")
      # if(is.null(data_long))
      #   stop("Missing valid data. Need to input long form of the data")
      # if(is.null(time_name) | is.null(id_name))
      #   stop("Missing time and id variable name. For long form data,
      #        need to input time_name and id_name.")
      id_boot <- sample(unique(data_long[,id_name]), replace = TRUE)
      data_boot_long <- data.frame(id_new = 1:length(id_boot), id = id_boot)
      data_boot_long <- merge(data_boot_long, data_long, all.x = TRUE)
      colnames(data_boot_long) <- c("id_old", id_name,
                                    colnames(data_boot_long)[3:ncol(data_boot_long)])
      llh_res <- j2r_llh(formula, data_boot_long, time_name, id_name,
                         trt_name, method = llh_method)
      ate_boot[,b] <- llh_res$llh_est
    }
    else if(method == "sequential likelihood"){
      # if(is.null(llh_method))
      #   stop("Missing specific likelihood method to use. Need to input
      #        likelihood-based method in `llh_method`. Available choice:
      #        PMM or CBMI")
      id_boot <- sample(unique(data_wide[,id_name]), replace = TRUE)
      data_boot <- data.frame(id_new = 1:length(id_boot), id = id_boot)
      data_boot <- merge(data_boot, data_wide, all.x = TRUE)
      data_boot <- cbind(data_boot[,-c(1:2)], data_boot[,1:2])
      colnames(data_boot) <- c(colnames(data_boot)[1:(ncol(data_boot)-2)],
                               "id_old", "id")
      llh_res <- j2r_seq_llh(formula, data_boot, trt_name, M, fit_model,
                             method = llh_method)
      ate_boot[,b] <- llh_res$llh_est
    }
  }
  var_boot <- apply(ate_boot, 1, stats::var)

  return(var_boot)
}

#' Sequential imputation
#'
#' For a given dataset with monotone missingness, conduct sequential imputation
#' based on linear / robust / rank regression
#'
#' @param seq_formula a list of formulas for sequential regression
#' @param data_wide wide-form of the data frame
#' @param trt_name name of the treatment variable
#' @param M imputation size
#' @param fit_model type of the analysis model. Available: lm, Rfit::rfit
#' (rank regression), MASS::rlm (robust regression)
#' @param estimator type of the estimator. Available: sample mean,
#' reg with interaction, reg without interaction
#'
#' @return Rubin's point and variance estimate
#' \itemize{
#'   \item{ate_est}{MI point estimation of the ATE}
#'   \item{var_rubin}{MI variance estimation by Rubin's rule}
#' }
#' @export
rr_seq_imp <- function(seq_formula, data_wide, trt_name, M, fit_model,
                       estimator = "reg with interaction"){
  base_name <- all.vars(seq_formula[[1]])[1]
  covar_name <- all.vars(seq_formula[[1]])[-1]
  y_name <- sapply(seq_formula, function(x) all.vars(x)[1])
  seq_fit <- seq_mi_fit(seq_formula, data_wide, trt_name, M, fit_model)
  beta_ctl <- seq_fit$beta_ctl_imp
  beta_trt <- seq_fit$beta_trt_imp
  sigma_ctl <- seq_fit$sigma_ctl
  sigma_trt <- seq_fit$sigma_trt
  imp_value <- j2r_imp(data_wide, beta_ctl, beta_trt, sigma_ctl, sigma_trt,
                       covar_name, y_name, trt_name)
  mi_res <- rubin_est(t(imp_value), data_wide = data_wide,
                      covar_name = covar_name, trt_name = trt_name,
                      M = M, fit_model = fit_model, estimator = estimator)

  return(list(ate_est = mi_res$ate_est,
              var_rubin = mi_res$var_rubin))
}

#' Inference using MI-based sequential based method, under J2R
#'
#' Given a wide-form data, compute the estimate of the ATE under J2R using
#' MI-based sequential method.
#'
#' @param seq_formula a list of formulas for sequential regression
#' @param data_wide wide-form of the data frame
#' @param trt_name name of the treatment variable
#' @param M imputation size
#' @param fit_model type of the analysis model. Available: lm, Rfit::rfit
#' (rank regression), MASS::rlm (robust regression)
#' @param id_name name of the id variable
#' @param estimator type of the estimator. Available: sample mean,
#' reg with interaction, reg without interaction. Default: reg with interaction
#' @param bootstrap Whether to use bootstrap variance
#' @param B the number of bootstrap replicate
#'
#' @return A list of MI estimation of the ATE
#' \itemize{
#'   \item{ate_est}{MI point estimation of the ATE}
#'   \item{var_rubin}{MI variance estimation by Rubin's rule}
#'   \item{var_boot}{MI variance estimation by nonparametric bootstrap}
#' }
#' @export
seqMI <- function(seq_formula, data_wide, trt_name, M, fit_model,
                  id_name, estimator = "reg with interaction",
                  bootstrap = FALSE, B = 100){

  seqMI_res <- rr_seq_imp(seq_formula, data_wide, trt_name, M, fit_model)
  if(bootstrap){
    var_boot <- nonpara_fn(formula = seq_formula, data_wide = data_wide,
                           id_name = id_name, trt_name = trt_name,
                           M = M, B = B,
                           fit_model = fit_model,
                           estimator = estimator,
                           method = "sequential MI")
  }
  return(list(ate_est = seqMI_res$ate_est,
              var_rubin = seqMI_res$var_rubin,
              var_boot = var_boot))
}


#' Coefficient matrix used to get variance of the transformed coefficients
#'
#' @param beta coefficients from lmer()
#' @param k dimension of the covariates
#' @param p the number of time point
#'
#' @return Coefficient matrix
#' @export
cov_beta_last <- function(beta, p, k){
  coef_mat <- matrix(0, k+1, length(beta))
  for(i in 1:(k+1)){
    coef_mat[i,i] <- 1
    coef_mat[i, i + (k+1) + (p-2)*i] <- 1
  }
  coef_mat
}


#' fit ANCOVA for complete data
#'
#' Fit ANCOVA for the complete data
#' @param data_complete_trt complete wide form data of the treatment group
#' @param last_name name of the response at the last time point
#' @param covar_name name of the covariates
#' @param fit_model type of the analysis model. Available: lm, Rfit::rfit
#' (rank regression), MASS::rlm (robust regression)
#'
#' @return A list of estimates:
#' \itemize{
#'   \item{beta_complete_trt}{estimated regression coefficients}
#'   \item{cov_complete_trt}{covariance of the estimated regression coefficients}
#' }
#'
#' @export
ancova_fit <- function(data_complete_trt, last_name, covar_name, fit_model){
  chg_trt <- data_complete_trt[,last_name]
  fit_complete_trt <- fit_model(stats::reformulate(covar_name, "chg_trt"), data = data_complete_trt)
  beta_complete_trt <- stats::coef(fit_complete_trt)
  cov_complete_trt <- stats::vcov(fit_complete_trt)

  return(list(beta_complete_trt = beta_complete_trt,
              cov_complete_trt = cov_complete_trt))
}


#' Likelihood-based method with estimation from delta method
#'
#' For a given dataset with monotone missingness, find the likelihood-based
#' inference of the ATE
#'
#' @param mmrm_formula formula for fixed effect MMRM under MAR
#' @param data_long long-form of the data frame
#' @param time_name name of the time variable
#' @param id_name name of the id variable
#' @param trt_name name of the treatment variable
#' @param method specific method for likelihood-based approach. Available: PMM (default),
#' CBMI
#'
#' @return Point estimation and variance estimation obtained by delta method:
#' \itemize{
#'   \item{llh_est}{point estimate}
#'   \item{var_llh}{variance estimate approximated by delta method}
#' }
#' @export
j2r_llh <- function(mmrm_formula, data_long, time_name, id_name, trt_name,
                    method = "PMM"){

  info_data <- summaryLong(mmrm_formula, data_long, time_name, id_name, trt_name)
  p <- info_data$p
  k <- info_data$k
  y_name <- info_data$y_name
  base_name <- y_name[1]
  last_name <- y_name[p]
  covar_name <- info_data$covar_name
  db_avaiable_ctl <- info_data$db_avaiable_ctl
  db_avaiable_trt <- info_data$db_avaiable_trt
  data_wide <- info_data$data_wide
  adj_mean <- info_data$adj_mean
  var_xbar <- info_data$var_xbar
  obs_pi <- info_data$obs_pi
  n2 <- info_data$n2
  db_complete_ctl <- info_data$data_complete_ctl
  db_complete_trt <- info_data$data_complete_trt

  new_mmrm_formula <- mmrmformula(mmrm_formula, time_name, id_name)

  fit_mmrm_ctl <- mmrm_fit(new_mmrm_formula, db_avaiable_ctl, p, k)
  fit_mmrm_trt <- mmrm_fit(new_mmrm_formula, db_avaiable_trt, p, k)

  cov_beta_ctl <- fit_mmrm_ctl$cov_beta
  beta_transform_ctl <- fit_mmrm_ctl$beta
  cov_beta_trt <- fit_mmrm_trt$cov_beta
  beta_transform_trt <- fit_mmrm_trt$beta

  cov4_ctl <- fit_mmrm_ctl$cov4
  cov4_trt <- fit_mmrm_trt$cov4

  # complete data
  fit_complete <- ancova_fit(data_complete_trt = db_complete_trt,
                             last_name = last_name,
                             covar_name = covar_name, fit_model = stats::lm)
  beta_complete_trt <- fit_complete$beta_complete_trt
  cov4_complete_trt <- fit_complete$cov_complete_trt
  approx_res <- approx_inference(beta_last_ctl = beta_transform_ctl[((k+1)*(p-1) + 1):length(beta_transform_ctl)],
                                 beta_last_trt = beta_transform_trt[((k+1)*(p-1) + 1):length(beta_transform_trt)],
                                 cov4_ctl = cov4_ctl, cov4_trt = cov4_trt,
                                 beta_complete_trt = beta_complete_trt,
                                 cov4_complete_trt = cov4_complete_trt,
                                 obs_pi = obs_pi, var_xbar = var_xbar,
                                 adj_mean = adj_mean, n2 = n2, method = method)
  llh_est <- approx_res$llh_est
  var_llh <- approx_res$var_llh

  return(list(llh_est = llh_est,
              var_llh = var_llh))
}

#' Likelihood-based method
#'
#' For a given dataset, find the likelihood-based inference of the ATE by
#' delta method or nonparametric bootstrap
#'
#' @param mmrm_formula formula for fixed effect MMRM under MAR
#' @param data_long long-form of the data frame
#' @param time_name name of the time variable
#' @param id_name name of the id variable
#' @param trt_name name of the treatment variable
#' @param method specific method for likelihood-based approach. Available: PMM (default),
#' CBMI
#' @param bootstrap whether to use bootstrap variance
#' @param B the number of bootstrap replicate. Default: B = 100
#'
#' @return Point estimation and variance estimation obtained by delta method:
#' \itemize{
#'   \item{ate_est}{point estimate}
#'   \item{var_llh}{variance estimate approximated by delta method}
#'   \item{var_boot}{variance estimate approximated by delta method}
#' }
#' @export
j2rLlh <- function(mmrm_formula, data_long, time_name, id_name, trt_name,
                   method = "PMM", bootstrap = FALSE, B = 100){

  llh_res <- j2r_llh(mmrm_formula, data_long,
                     time_name = time_name, id_name = id_name,
                     trt_name = trt_name, method = method)
  if(bootstrap){
    var_boot <- nonpara_fn(formula = mmrm_formula, data_long = data_long,
                           time_name = time_name, id_name = id_name,
                           trt_name = trt_name, B = B,
                           method = "naive likelihood", llh_method = method)
  }
  return(list(ate_est = llh_res$llh_est,
              var_llh = llh_res$var_llh,
              var_boot = var_boot))
}


#' Summary of wide-form data frame given sequential regression formula
#'
#' @param seq_formula a list of formulas for sequential regression
#' @param data_wide wide-form of the data
#' @param trt_name name of the treatment variable
#'
#' @return A list of summary of the long form data:
#' \itemize{
#'   \item{p} number of time points
#'   \item{covar_name}{name of the covariates}
#'   \item{k}{dimension of the covariates}
#'   \item{last_name}{name of the response at last time point}
#'   \item{adj_mean}{overall covariate mean}
#'   \item{var_xbar}{covariance matrix of the covariates}
#'   \item{obs_pi}{marginal probability of the observed data for each treatment}
#'   \item{n1}{number of subjects in the control group}
#'   \item{n2}{number of subjects in the treatment group}
#'   \item{data_complete_ctl}{complete data in the control group}
#'   \item{data_complete_trt}{complete data in the treatment group}
#' }
#' @export
summaryWide <- function(seq_formula, data_wide, trt_name){
  p <- length(seq_formula)
  base_name <- all.vars(seq_formula[[1]])[1]
  covar_name <- all.vars(seq_formula[[1]])[-1]
  k <- length(covar_name)
  last_name <- all.vars(seq_formula[[p]])[1]
  data_wide_ctl <- data_wide[which(data_wide[,trt_name] == 1),]
  data_wide_trt <- data_wide[which(data_wide[,trt_name] == 2),]
  n1 <- nrow(data_wide_ctl)
  n2 <- nrow(data_wide_trt)
  obs_pi <- 1 - c(sum(is.na(data_wide[,last_name] & data_wide[,trt_name] == 1))/n1,
                  sum(is.na(data_wide[,last_name] & data_wide[,trt_name] == 2))/n2)

  data_complete_ctl <- data_wide_ctl[which(!is.na(data_wide_ctl[,last_name])),]
  data_complete_trt <- data_wide_trt[which(!is.na(data_wide_trt[,last_name])),]

  adj_mean <- colMeans(data_wide[,covar_name])
  var_xbar <- stats::cov(data_wide[,covar_name])/nrow(data_wide)

  return(list(p = p,
              base_name = base_name,
              covar_name = covar_name,
              k = k,
              last_name = last_name,
              obs_pi = obs_pi,
              adj_mean = adj_mean,
              var_xbar = var_xbar,
              n1 = n1, n2 = n2,
              data_complete_ctl = data_complete_ctl,
              data_complete_trt = data_complete_trt))
}

#' Variance estimation of likelihood-based inference approximated by delta
#' method
#'
#' Compute variance estimate of the likelihood-based inference using delta
#' method
#'
#' @param beta_last_ctl regression coefficient at last time point in control
#' @param beta_last_trt regression coefficient at last time point in treatment
#' @param cov4_ctl covariance of the coefficient at last time point in control
#' @param cov4_trt covariance of the coefficient at last time point in treatment
#' @param beta_complete_trt estimated regression coefficients for complete data
#' in treatment
#' @param cov4_complete_trt covariance of the estimated regression coefficients
#' for complete data in treatment
#' @param obs_pi marginal observed probability in each group
#' @param var_xbar covariance of the covariates
#' @param adj_mean overall mean of the covariates
#' @param n2 number of subjects in treatment
#' @param method name of the likelihood-based method. Available: PMM, CBMI
#'
#' @return Estimation results:
#' \itemize{
#'   \item{llh_est}{point estimate}
#'   \item{var_llh}{variance estimate approximated by delta method}
#' }
#' @export
approx_inference <- function(beta_last_ctl, beta_last_trt, cov4_ctl, cov4_trt,
                             beta_complete_trt = NULL, cov4_complete_trt = NULL,
                             obs_pi, var_xbar, adj_mean, n2, method){
  mut_ctl <- beta_last_ctl[1] + sum(beta_last_ctl[-1]*adj_mean)
  mut_trt <- beta_last_trt[1] + sum(beta_last_trt[-1]*adj_mean)
  coef_mat <- matrix(c(1, adj_mean), nrow = 1)
  coef_beta_mat_ctl <- matrix(beta_last_ctl[-1], nrow = 1)
  var_theta_ctl <- as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat)) +
    as.numeric(coef_beta_mat_ctl%*%var_xbar%*%t(coef_beta_mat_ctl))

  if(method == "PMM"){
    # point estimate
    mu_ctl <- mut_ctl
    mu_trt <- mut_trt*obs_pi[2] + mut_ctl*(1 - obs_pi[2])

    # variance estimate
    coef_beta_mat_trt <- matrix(beta_last_ctl[-1]*obs_pi[2] +
                                  beta_last_ctl[-1]*(1 - obs_pi[2]), nrow = 1)
    coef_beta_mat_diff <- matrix((beta_last_trt[-1] -
                                    beta_last_ctl[-1])*obs_pi[2], nrow = 1)
    var_theta_trt <- as.numeric(coef_mat%*%cov4_trt%*%t(coef_mat))*obs_pi[2]^2 +
      as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat))*(1 - obs_pi[2])^2 +
      (mut_trt - mut_ctl)^2*obs_pi[2]*(1 - obs_pi[2])/n2 +
      as.numeric(coef_beta_mat_trt%*%var_xbar%*%t(coef_beta_mat_trt))
    var_theta_diff <- as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat))*obs_pi[2]^2 +
      as.numeric(coef_mat%*%cov4_trt%*%t(coef_mat))*obs_pi[2]^2 +
      (mut_trt - mut_ctl)^2*obs_pi[2]*(1 - obs_pi[2])/n2 +
      as.numeric(coef_beta_mat_diff%*%var_xbar%*%t(coef_beta_mat_diff))

  }

  if(method == "CBMI"){
    if(is.null(beta_complete_trt) | is.null(cov4_complete_trt))
      stop("Need to input the estimates for completers in treatment")
    mu_complete_trt <- beta_complete_trt[1] + sum(beta_complete_trt[-1]*adj_mean)
    cov4_trt <- cov4_complete_trt
    # point estimate
    mu_ctl <- mut_ctl
    mu_trt <- obs_pi[2]*mu_complete_trt + (1 - obs_pi[2])*mut_ctl

    # Variance estimate
    coef_beta_mat_trt <- matrix(beta_complete_trt[-1]*obs_pi[2] +
                                  beta_last_ctl[-1]*(1 - obs_pi[2]), nrow = 1)
    coef_beta_mat_diff <- matrix((beta_complete_trt[-1] - beta_last_ctl[-1])*
                                   sqrt(obs_pi[2]^2 + obs_pi[2]*(1 - obs_pi[2])/n2),
                                 nrow = 1)
    var_theta_trt <- as.numeric(coef_mat%*%cov4_trt%*%t(coef_mat))*(obs_pi[2]^2 + obs_pi[2]*(1 - obs_pi[2])/n2) +
      as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat))*((1 - obs_pi[2])^2 + obs_pi[2]*(1 - obs_pi[2])/n2) +
      (mu_complete_trt - mut_ctl)^2*obs_pi[2]*(1 - obs_pi[2])/n2 +
      as.numeric(coef_beta_mat_trt%*%var_xbar%*%t(coef_beta_mat_trt))
    var_theta_diff <- as.numeric(coef_mat%*%cov4_ctl%*%t(coef_mat))*(obs_pi[2]^2 + obs_pi[2]*(1 - obs_pi[2])/n2) +
      as.numeric(coef_mat%*%cov4_trt%*%t(coef_mat))*(obs_pi[2]^2 + obs_pi[2]*(1 - obs_pi[2])/n2) +
      (mu_complete_trt - mut_ctl)^2*obs_pi[2]*(1 - obs_pi[2])/n2 +
      as.numeric(coef_beta_mat_diff%*%var_xbar%*%t(coef_beta_mat_diff))
  }

  llh_est <- c(mu_ctl, mu_trt, mu_trt - mu_ctl)
  var_llh <- c(var_theta_ctl, var_theta_trt, var_theta_diff)

  return(list(llh_est = llh_est,
              var_llh = var_llh))
}

#' Sequential likelihood-based inference from delta method
#'
#' Compute variance estimate of the likelihood-based inference using delta
#' method
#'
#' @param seq_formula a list of formulas for sequential regression
#' @param data_wide wide-form of the data
#' @param trt_name name of the treatment variable
#' @param M imputation size
#' @param fit_model type of the analysis model. Available: lm, Rfit::rfit
#' (rank regression), MASS::rlm (robust regression)
#' @param method name of the likelihood-based method. Available: PMM, CBMI
#'
#' @return Point estimation and variance estimation obtained by delta method:
#' \itemize{
#'   \item{llh_est}{point estimate}
#'   \item{var_llh}{variance estimate approximated by delta method}
#' }
#' @export
j2r_seq_llh <- function(seq_formula, data_wide, trt_name, M = 10,
                        fit_model, method = "PMM"){

  info_data <- summaryWide(seq_formula, data_wide, trt_name)
  obs_pi <- info_data$obs_pi
  covar_name <- info_data$covar_name
  k <- info_data$k
  adj_mean <- info_data$adj_mean
  var_xbar <- info_data$var_xbar
  n2 <- info_data$n2
  data_complete_ctl <- info_data$data_complete_ctl
  data_complete_trt <- info_data$data_complete_trt
  # available data
  seq_fit <- seq_mi_fit(seq_formula, data_wide, trt_name, M, fit_model)
  beta_ctl <- seq_fit$beta_rubin_ctl
  beta_trt <- seq_fit$beta_rubin_trt
  cov4_ctl <- seq_fit$var_rubin_ctl
  cov4_trt <- seq_fit$var_rubin_trt
  # complete data
  fit_complete <- ancova_fit(data_complete_trt,
                             last_name = info_data$last_name,
                             covar_name = covar_name, fit_model = fit_model)
  beta_complete_trt <- fit_complete$beta_complete_trt
  cov4_complete_trt <- fit_complete$cov_complete_trt
  approx_res <- approx_inference(beta_last_ctl = beta_ctl,
                                 beta_last_trt = beta_trt,
                                 cov4_ctl = cov4_ctl, cov4_trt = cov4_trt,
                                 beta_complete_trt = beta_complete_trt,
                                 cov4_complete_trt = cov4_complete_trt,
                                 obs_pi = obs_pi, var_xbar = var_xbar,
                                 adj_mean = adj_mean, n2 = n2, method = method)
  llh_est <- approx_res$llh_est
  var_llh <- approx_res$var_llh

  return(list(llh_est = llh_est,
              var_llh = var_llh))
}

#' Inference using MI-based sequential based method, under J2R
#'
#' For a given dataset with monotone missingness, impute the missing data from
#' the conditional distribution given observed under J2R by proper MI, output
#' the imputed dataset
#'
#' @param seq_formula a list of formulas for sequential regression
#' @param data_wide wide-form of the data frame
#' @param trt_name name of the treatment variable
#' @param M imputation size
#' @param fit_model type of the analysis model. Available: lm, Rfit::rfit
#' (rank regression), MASS::rlm (robust regression)
#' @param id_name name of the id variable
#' @param method name of the likelihood-based method. Available: PMM, CBMI
#' @param bootstrap Whether to use bootstrap variance
#' @param B the number of bootstrap replicate
#'
#' @return Point estimation and variance estimation obtained by delta method:
#' \itemize{
#'   \item{llh_est}{point estimate}
#'   \item{var_llh}{variance estimate approximated by delta method}
#'   \item{var_boot}{variance estimate approximated by delta method}
#' }
#' @export
seqLlh <- function(seq_formula, data_wide, trt_name, M, fit_model,
                  id_name, method = "PMM",
                  bootstrap = FALSE, B = 100){

  seqLlh_res <- j2r_seq_llh(seq_formula, data_wide,
                            trt_name, M, fit_model, method)
  if(bootstrap){
    var_boot <- nonpara_fn(formula = seq_formula, data_wide = data_wide,
                           id_name = id_name, trt_name = trt_name,
                           M = M, B = B,
                           fit_model = fit_model,
                           method = "sequential likelihood",
                           llh_method = method)
  }
  return(list(ate_est = seqLlh_res$llh_est,
              var_llh = seqLlh_res$var_llh,
              var_boot = var_boot))
}

#' Assign missing patterns
#'
#' For wide form of the data, assign the missing pattern indicator
#' @param data_wide wide form data
#' @param y_name names of the longitudinal response
#'
#' @return a wide-form data frame with the assigned missing patterns
#' @export
add_pattern <- function(data_wide, y_name){
  p <- length(y_name)
  # find id with intermittent missing
  inter_ind <- c(0)
  count <- 1
  for(i in 1:nrow(data_wide)){
    for(j in 1:(p-1)){
      if(is.na(data_wide[i,y_name[j]]) & !is.na(data_wide[i,y_name[j+1]])){
        inter_ind[count] <- i
        count <- count + 1
        break
      }
    }
  }
  mis_pattern <- apply(data_wide, 1, function(x) sum(is.na(x[y_name])))
  pattern <- rep(1, nrow(data_wide))
  pattern[inter_ind] <- 0
  for(i in 1:length(mis_pattern[-inter_ind])){
    if(mis_pattern[-inter_ind][i] == 0)
      pattern[-inter_ind][i] <- 1
    for(j in 1:p){
      if(mis_pattern[-inter_ind][i] == j)
        pattern[-inter_ind][i] <- p + 1 - j
    }
  }
  data_wide$pattern <- pattern

  data_wide
}



