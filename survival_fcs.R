library(survival)
library(kernlab)
library(quadprog)
library(Matrix)
library(MASS)
library(CBPS)

causal_survfit <- function(formula, data, method = "mksd", adjust_gp = NULL, gamma = NULL, wu = 10, eps = NULL, adj_cen = T) {
  mf <- model.frame(formula(formula), data)
  data <- data[rownames(mf), ]
  if (method == "mksd")
    w <- mksd_weight(formula, data, adjust_gp, gamma, wu)
  else if (method == "csw")
    w <- csw_weight(formula, data, adjust_gp, gamma, wu, eps)
  else if (method == "ncsw")
    w <- ncsw_weight(formula, data, adjust_gp, gamma, wu, eps)
  else if (method == "psatt")
    w <- psatt_weight(formula, data, adjust_gp)
  else if (method == "ipw")
    w <- ipw_weight(formula, data)
  else if (method == "ow")
    w <- ow_weight(formula, data)
  else if (method == "CBPS")
    w <- CBPS_weight(formula, data)
  if (adj_cen)
    w <- w * ipcw(formula, data)
  new_form <- formula(paste0(names(mf)[1], "~", names(mf)[2]))
  surv_object <- survfit(new_form, data, w)
  surv_object$call$formula <- new_form
  return(list(surv_object, mf, w))
}

ipcw <- function(formula, data) {
  mf <- model.frame(formula(formula), data)
  is_treat <- mf[, 2]
  z <- as.numeric(is_treat)
  covM <- model.matrix(formula(formula), data)
  covM <- covM[, - 2, drop = F]
  covM1 <- covM[z == 1, ]
  covM2 <- covM[z == 2, ]
  x1 <- data[z == 1, ]
  x2 <- data[z == 2, ]
  censor <- tail(as.numeric(mf[, 1]), nrow(mf))
  timevec <- head(as.numeric(mf[, 1]), nrow(mf))
  time1 <- timevec[z == 1]
  time2 <- timevec[z == 2]
  form <- formula(paste0(names(mf)[1], "~", paste(names(mf)[- c(1, 2)], collapse = "+")))
  model1 <- survreg(form, data = x1, dist = 'weibull', score = T, control = survreg.control(maxiter = 350))
  theta1 = - model1$coefficients / model1$scale
  gamma1 = 1 / model1$scale
  kc1 = exp(- exp(c(covM1 %*% theta1)) * time1 ^ gamma1)
  model2 <- survreg(form, data = x2, dist = 'weibull', score = T, control = survreg.control(maxiter = 350))
  theta2 = - model2$coefficients / model2$scale
  gamma2 = 1 / model2$scale
  kc2 = exp(- exp(c(covM2 %*% theta2)) * time2 ^ gamma2)
  kc <- rep(1, nrow(data))
  kc[z == 1] <- kc1
  kc[z == 2] <- kc2
  kc <- cbind(1 / (1 - kc), 1 / kc)
  cs <- rep(1, nrow(data))
  cs[censor == unique(censor)[1]] <- kc[censor == unique(censor)[1], 1]
  cs[censor == unique(censor)[2]] <- kc[censor == unique(censor)[2], 2]
  return(cs)
}

mksd_weight <- function(formula, data, adjust_gp, gamma = NULL, wu = 10) {
  mf <- model.frame(formula(formula), data)
  is_treat <- mf[, 2]
  covM <- model.matrix(formula(formula), data)
  covM <- covM[, - c(1, 2), drop = F]
  x1 <- covM[is_treat != adjust_gp, ]
  x2 <- covM[is_treat == adjust_gp, ]
  w <- rep(1, nrow(data))
  if(is.null(gamma)) gamma <- kernelwidth(covM)
  w_adjust <- mksd(x1, x2, gamma, wu)
  w[is_treat == adjust_gp] <- w_adjust
  return(w)
}

csw_weight <- function(formula, data, adjust_gp, gamma = NULL, B = 10, eps = NULL) {
  mf <- model.frame(formula(formula), data)
  is_treat <- mf[, 2]
  covM <- model.matrix(formula(formula), data)
  covM <- covM[, - c(1, 2), drop = F]
  x1 <- covM[is_treat != adjust_gp, ]
  x2 <- covM[is_treat == adjust_gp, ]
  w <- rep(1, nrow(data))
  if(is.null(gamma)) gamma <- kernelwidth(covM)
  if(is.null(eps)) eps <- eps_compute(nrow(x2))
  w_adjust <- kmm(x1, x2, gamma, B, eps)
  w[is_treat == adjust_gp] <- w_adjust
  return(w)
}

ncsw_weight <- function(formula, data, adjust_gp, gamma = NULL, B = 10, eps = NULL) {
  mf <- model.frame(formula(formula), data)
  is_treat <- mf[, 2]
  covM <- model.matrix(formula(formula), data)
  covM <- scale(covM[, - c(1, 2), drop = F])
  covM <- covM[, colSums(is.nan(covM)) == 0]
  x1 <- covM[is_treat != adjust_gp, ]
  x2 <- covM[is_treat == adjust_gp, ]
  w <- rep(1, nrow(data))
  # wfunc <- densratio(x1, x2, verbose = F)
  # w_adjust <- wfunc$compute_density_ratio(x2)
  if(is.null(gamma)) gamma <- kernelwidth(covM)
  if(is.null(eps)) eps <- eps_compute(nrow(x2))
  w_adjust <- kmm(x1, x2, gamma, B, eps)
  w[is_treat == adjust_gp] <- w_adjust
  return(w)
}

psatt_weight <- function(formula, data, adjust_gp) {
  mf <- model.frame(formula(formula), data)
  is_treat <- mf[, 2]
  form <- formula(paste0(names(mf)[2], "~", paste(names(mf)[- c(1, 2)], collapse = "+")))
  fit <- glm(formula = form, data = data, family = binomial(link = "logit"))
  e <- fit$fitted.values
  w <- rep(1, nrow(data))
  w_adjust <- if(levels(is_treat)[1] == adjust_gp) e / (1 - e) else (1 - e) / e
  w[is_treat == adjust_gp] <- w_adjust[is_treat == adjust_gp]
  return(w)
}

ipw_weight <- function(formula, data) {
  mf <- model.frame(formula(formula), data)
  is_treat <- mf[, 2]
  form <- formula(paste0(names(mf)[2], "~", paste(names(mf)[- c(1, 2)], collapse = "+")))
  fit <- glm(formula = form, data = data, family = binomial(link = "logit"))
  e <- fit$fitted.values
  e <- cbind(1 / (1 - e), 1 / e)
  w <- rep(1, nrow(data))
  z <- as.numeric(is_treat)
  w[z == 1] <- e[z == 1, 1]
  w[z == 2] <- e[z == 2, 2]
  return(w)
}

ow_weight <- function(formula, data) {
  mf <- model.frame(formula(formula), data)
  is_treat <- mf[, 2]
  form <- formula(paste0(names(mf)[2], "~", paste(names(mf)[- c(1, 2)], collapse = "+")))
  fit <- glm(formula = form, data = data, family = binomial(link = "logit"))
  e <- fit$fitted.values
  e <- cbind(e, 1 - e)
  w <- rep(1, nrow(data))
  z <- as.numeric(is_treat)
  w[z == 1] <- e[z == 1, 1]
  w[z == 2] <- e[z == 2, 2]
  return(w)
}

CBPS_weight <- function(formula, data) {
  mf <- model.frame(formula(formula), data)
  covM <- model.matrix(formula(formula), data)
  covM <- covM[, - c(1, 2), drop = F]
  form <- paste0(names(mf)[2], "~", paste(colnames(covM), collapse = "+"))
  mf_t <- cbind(mf[2], data.frame(covM))
  CB_fit <- CBPS(form, mf_t, ATT = 0)
  w <- CB_fit$weights
  return(w)
}

mksd <- function(x1, x2, gamma, wu) {
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  sigma <- cov(rbind(x1, x2))
  md <- mddot(gamma = gamma, sigma = sigma)
  K <- kernelMatrix(md, x2)
  Dmat <- as.matrix(nearPD(K)$mat)
  kappa <- kernelMatrix(md, x2, x1)
  dvec <- n2 / n1 * rowSums(kappa)
  Amat <- cbind(rep(1, n2), diag(n2), - diag(n2))
  bvec <- c(n2, rep(0, n2), rep(- wu, n2))
  meq <- 1
  w <- solve.QP(Dmat, dvec, Amat, bvec, meq)$solution
  w <- ifelse(w < 0, 0, w)
  return(w)
}

mddot <- function(gamma, sigma) {
  rval <- function(x, y) {
    md <- try(mahalanobis(x, y, sigma), silent = T)
    if ('try-error' %in% class(md)) {
      inv_sigma <- ginv(sigma)
      dot <- exp(- gamma * mahalanobis(x, y, inv_sigma, T))
    }
    else {
      dot <- exp(- gamma * md)
    }
    return(dot)
  }
  return(new("kernel", .Data = rval, kpar = list(gamma = gamma, sigma = sigma)))
}

kersd <- function(formula, data, w, gamma) {
  covM <- model.matrix(formula(formula), data)
  z <- covM[, 2]
  covM <- covM[, - c(1, 2), drop = F]
  x1 <- covM[z == 0, ]
  x2 <- covM[z == 1, ]
  w1 <- w[z == 0]
  w2 <- w[z == 1]
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  sigma <- cov(covM)
  md <- mddot(gamma = gamma, sigma = sigma)
  kc <- kernelMatrix(md, x2)
  ktc <- kernelMatrix(md, x1, x2)
  kt <- kernelMatrix(md, x1)
  ksd <- 1 / n2 ^ 2 * rowSums(crossprod(w2, kc) * w2) - 2 / n1 / n2 * rowSums(crossprod(w1, ktc) * w2) + 1 / n1 ^ 2 * rowSums(crossprod(w1, kt) * w1)
  return(sqrt(ksd))
}

kmm <- function(x1, x2, gamma, B, eps) {
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  rbf <- rbfdot(sigma = gamma)
  K <- kernelMatrix(rbf, x2)
  Dmat <- as.matrix(nearPD(K)$mat)
  kappa <- kernelMatrix(rbf, x2, x1)
  dvec <- n2 / n1 * rowSums(kappa)
  Amat <- cbind(rep(- 1, n2), rep(1, n2), diag(n2), - diag(n2))
  bvec <- c(- n2 * (eps + 1), - n2 * (eps - 1), rep(0, n2), rep(- B, n2))
  beta <- solve.QP(Dmat, dvec, Amat, bvec)$solution
  beta <- ifelse(beta < 0, 0, beta)
  return(beta)
}

kernelwidth <- function(data) {
  dis <- as.numeric(dist(data))
  width <- median(dis)
  return(width)
}

eps_compute <- function(m) {
  eps <- (sqrt(m) - 1) / sqrt(m)
  return(eps)
}