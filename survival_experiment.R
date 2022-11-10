source("survival_fcs.R")
library(parallel)
library(sjstats)
library(car)
library(survtmle)

surv_causal <- function(fit) {
  surv_adj <- fit[grepl(adjust_gp, names(fit$strata), fixed = T)]
  surv_tar <- fit[! grepl(adjust_gp, names(fit$strata), fixed = T)]
  ate <- summary(surv_tar, times = time_seq)$surv - summary(surv_adj, times = time_seq)$surv
  return(ate)
}

tmle <- function(formula, data) {
  mf <- model.frame(formula(formula), data)
  is_treat <- as.numeric(mf[, 2])
  covM <- model.matrix(formula(formula), data)
  covM <- data.frame(covM[, - c(1, 2), drop = F])
  form1 <- paste(colnames(covM), collapse = "+")
  form2 <- paste0("trt+", paste(colnames(covM), collapse = "+"))
  fit1 <- lapply(as.list(time_seq), 
                 survtmle, 
                 ftime = data$time, 
                 ftype = data$status, 
                 trt = is_treat, 
                 adjustVars = covM, 
                 glm.trt = form1,
                 glm.ftime = form2,
                 glm.ctime = form2,
                 method = "mean")
  fit2 <- do.call(Map, c(f = list, fit1))
  fit <- do.call(Map, c(f = c, fit2$est))
  ate <- fit[[1]] - fit[[2]]
  return(ate)
}

cov_test <- function(w, data) {
  mf <- model.frame(formula(formula), data)
  X9X10 <- mf$X9 * mf$X10
  X4X5 <- mf$X4 * mf$X5
  mf <- data.frame(is_treat = mf$is_treat, X5 = mf$X5, X10 = mf$X10, X4X5, X9X10)
  data <- cbind(mf, data.frame(w))
  lhs <- names(mf)[2 : length(mf)]
  forms <- as.list(paste0(lhs, "~", names(mf)[1], "+w"))
  formulas <- lapply(forms, as.formula)
  list_test <- lapply(formulas, weighted_mannwhitney, data)
  list_p <- do.call(Map, c(f = c, list_test))
  return(list_p[[8]])
}

leve_test <- function(data) {
  lhs <- names(data)[2 : length(data)]
  forms <- as.list(paste0(lhs, "~", names(data)[1]))
  formulas <- lapply(forms, as.formula)
  list_test <- lapply(formulas, leveneTest, data)
  list_p <- na.omit(do.call(rbind, list_test))
  return(list_p[, 3])
}

calc_rmse <- function(x) {
  rmse <- sqrt(mean((x - hr_treat) ^ 2, na.rm = T))
  return(rmse)
}

calc_bias <- function(x) {
  bias <- mean(x - hr_treat, na.rm = T)
  return(bias)
}

adjust_gp <- "no"
wu <- 10
gamma_seq = c(0.0001, 0.001, 0.01, 0.1, 1)
hr_treat <- 0
sigma_seq <- c(0.3, 0.7)
time_seq <- c(1, 2, 4, 8, 16)
path <- "./survival_simdata/"
file_num <- 1000
feature_num <- 10

cov_names <- paste0("+ X", 1 : feature_num)
cov_form <- paste(cov_names, collapse = " ")
formula <- paste("Surv(time, status) ~ is_treat", cov_form)
cores <- detectCores()
for (gamma in gamma_seq) {
  for (sigma in sigma_seq) {
    filedir = paste0(path, sigma, "/")
    data_list <- mclapply(as.list(paste0(filedir, 1 : file_num, ".csv")),
                          read.csv,
                          stringsAsFactors = T,
                          mc.cores = cores)
    
    surv_mksd <- mclapply(data_list,
                          causal_survfit,
                          formula = formula,
                          method = "mksd",
                          adjust_gp = adjust_gp,
                          gamma = gamma,
                          wu = wu,
                          mc.cores = cores)
    list_mksd <- do.call(Map, c(f = list, surv_mksd))
    result_mksd <- mclapply(list_mksd[[1]],
                            surv_causal,
                            mc.cores = cores)
    df_mksd <- do.call(rbind, result_mksd)
    rmse_mksd <- apply(df_mksd, 2, calc_rmse)
    bias_mksd <- apply(df_mksd, 2, calc_bias)
    mksd_cb <- mcmapply(cov_test,
                        list_mksd[[3]],
                        data_list,
                        mc.cores = cores)
    meth <- rep("mksd", file_num)
    mksd_cb <- cbind(meth, t(mksd_cb))
    leve_mksd <- cbind(meth = factor(meth), data.frame(df_mksd))
    
    # surv_csw <- mclapply(data_list,
    #                      causal_survfit,
    #                      formula = formula,
    #                      method = "csw",
    #                      adjust_gp = adjust_gp,
    #                      gamma = gamma,
    #                      wu = wu,
    #                      mc.cores = cores)
    # list_csw <- do.call(Map, c(f = list, surv_csw))
    # result_csw <- mclapply(list_csw[[1]],
    #                        surv_causal,
    #                        mc.cores = cores)
    # df_csw <- do.call(rbind, result_csw)
    # rmse_csw <- apply(df_csw, 2, calc_rmse)
    # bias_csw <- apply(df_csw, 2, calc_bias)
    # csw_cb <- mcmapply(cov_test,
    #                    list_csw[[3]],
    #                    data_list,
    #                    mc.cores = cores)
    # meth <- rep("csw", file_num)
    # csw_cb <- cbind(meth, t(csw_cb))
    
    surv_ncsw <- mclapply(data_list,
                          causal_survfit,
                          formula = formula,
                          method = "ncsw",
                          adjust_gp = adjust_gp,
                          gamma = gamma,
                          wu = wu,
                          mc.cores = cores)
    list_ncsw <- do.call(Map, c(f = list, surv_ncsw))
    result_ncsw <- mclapply(list_ncsw[[1]],
                            surv_causal,
                            mc.cores = cores)
    df_ncsw <- do.call(rbind, result_ncsw)
    rmse_ncsw <- apply(df_ncsw, 2, calc_rmse)
    bias_ncsw <- apply(df_ncsw, 2, calc_bias)
    ncsw_cb <- mcmapply(cov_test,
                        list_ncsw[[3]],
                        data_list,
                        mc.cores = cores)
    meth <- rep("ncsw", file_num)
    ncsw_cb <- cbind(meth, t(ncsw_cb))
    leve_ncsw <- cbind(meth = factor(meth), data.frame(df_ncsw))
    
    surv_ipw <- mclapply(data_list,
                         causal_survfit,
                         formula = formula,
                         method = "ipw",
                         mc.cores = cores)
    list_ipw <- do.call(Map, c(f = list, surv_ipw))
    result_ipw <- mclapply(list_ipw[[1]],
                           surv_causal,
                           mc.cores = cores)
    df_ipw <- do.call(rbind, result_ipw)
    rmse_ipw <- apply(df_ipw, 2, calc_rmse)
    bias_ipw <- apply(df_ipw, 2, calc_bias)
    ipw_cb <- mcmapply(cov_test,
                       list_ipw[[3]],
                       data_list,
                       mc.cores = cores)
    meth <- rep("ipw", file_num)
    ipw_cb <- cbind(meth, t(ipw_cb))
    
    surv_ow <- mclapply(data_list,
                        causal_survfit,
                        formula = formula,
                        method = "ow",
                        mc.cores = cores)
    list_ow <- do.call(Map, c(f = list, surv_ow))
    result_ow <- mclapply(list_ow[[1]],
                          surv_causal,
                          mc.cores = cores)
    df_ow <- do.call(rbind, result_ow)
    rmse_ow <- apply(df_ow, 2, calc_rmse)
    bias_ow <- apply(df_ow, 2, calc_bias)
    ow_cb <- mcmapply(cov_test,
                      list_ow[[3]],
                      data_list,
                      mc.cores = cores)
    meth <- rep("ow", file_num)
    ow_cb <- cbind(meth, t(ow_cb))
    
    surv_CBPS <- mclapply(data_list,
                          causal_survfit,
                          formula = formula,
                          method = "CBPS",
                          mc.cores = cores)
    list_CBPS <- do.call(Map, c(f = list, surv_CBPS))
    result_CBPS <- mclapply(list_CBPS[[1]],
                            surv_causal,
                            mc.cores = cores)
    df_CBPS <- do.call(rbind, result_CBPS)
    rmse_CBPS <- apply(df_CBPS, 2, calc_rmse)
    bias_CBPS <- apply(df_CBPS, 2, calc_bias)
    CBPS_cb <- mcmapply(cov_test,
                        list_CBPS[[3]],
                        data_list,
                        mc.cores = cores)
    meth <- rep("CBPS", file_num)
    CBPS_cb <- cbind(meth, t(CBPS_cb))
    
    result_tmle <- mclapply(data_list,
                            tmle,
                            formula = formula,
                            mc.cores = cores)
    df_tmle <- do.call(rbind, result_tmle)
    rmse_tmle <- apply(df_tmle, 2, calc_rmse)
    bias_tmle <- apply(df_tmle, 2, calc_bias)
    
    surv_crude <- mclapply(data_list,
                           survfit,
                           formula = as.formula("Surv(time, status) ~ is_treat"),
                           mc.cores = cores)
    result_crude <- mclapply(surv_crude,
                             surv_causal,
                             mc.cores = cores)
    df_crude <- do.call(rbind, result_crude)
    rmse_crude <- apply(df_crude, 2, calc_rmse)
    bias_crude <- apply(df_crude, 2, calc_bias)
    crude_cb <- mcmapply(cov_test,
                         rep(list(rep(1, nrow(data_list[[1]]))), file_num),
                         data_list,
                         mc.cores = cores)
    meth <- rep("crude", file_num)
    crude_cb <- cbind(meth, t(crude_cb))
    
    leve_data <- rbind(leve_mksd, leve_ncsw)
    leve_p <- leve_test(leve_data)
    rmse <- rbind(rmse_mksd, rmse_ncsw, rmse_ipw, rmse_ow, rmse_CBPS, rmse_crude, rmse_tmle, leve_p)
    bias <- rbind(bias_mksd, bias_ncsw, bias_ipw, bias_ow, bias_CBPS, bias_crude, bias_tmle)
    p_value <- rbind(mksd_cb, ncsw_cb, ipw_cb, ow_cb, CBPS_cb, crude_cb)
    write.csv(rmse, file = paste0(filedir, gamma, "rmse.csv"))
    write.csv(bias, file = paste0(filedir, gamma, "bias.csv"))
    write.csv(p_value, file = paste0(filedir, gamma, "p_value.csv"), row.names = F)
  }
}