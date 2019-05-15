#' Stepwise model search
#'
#' @description \code{model_search} performs a stepwise log-linear model search.
#' @param x a data frame with the incomplete lists and, if present, covariates.
#' For the specification of the variables in \code{x}, see Details.
#' @param lists a vector with the column numbers of \code{x} with the lists.
#' @param k_max a scalar restricting the \emph{k}-order interaction terms to \code{k-max}.
#' @param year character string with the name of the variable in \code{x} representing the
#' years of the observations, if present.
#' @param degree_year a scalar denoting the polynomial degree of the variable \code{year},
#' if present.
#' @return Prints a table with the degrees of freedom, deviance, AIC, AICc, BIC and population
#' size estimates of the fitted models. The AIC, AICc and BIC are scaled in such a way that
#' there minimum values equal 0. As input argument for the function \code{\link{boot_mse}} a
#' list with the following objects is returned:
#' \item{models}{the table with the fit statistics of the fitted models.}
#' \item{formulas}{the formulas of the fitted models.}
#' \item{coefs}{the parameter estimates of the fitted models.}
#' \item{fits}{a list with the fitted frequencies of the models in the search path.}
#' \item{minima}{a vector with the respective unscaled minima of the AIC, AICc and BIC.}
#' \item{d}{the data frame \code{x} in the form of a contingency table.}
#' \item{d_year}{the data frame \code{d} with \code{year} (if present) coded as polynomial.}
#' \item{obs}{a vector distinguishing between observations and structural zeros in \code{d}.}
#' \item{lists}{a vector with the column numbers of \code{d} with the lists.}
#' \item{year}{a scalar indicating the column number of the variable \code{year} in \code{d}.}
#' @details \code{x} should include at least two lists, and the lists should be coded 1 if the
#' observation is in that list, and 0 if it is not.
#'
#' \code{x} can be either in an \eqn{n x p} format, where \eqn{n} is the number of observed
#' individuals and \eqn{p} the number of variables (lists plus covariates), or in an
#' \eqn{m x p + 1} format, where \eqn{m} is the number of all possible combination of the variable
#' levels, and an additional variable called \code{Freq} with the frequencies (including the
#' structural zeros).
#'
#' @importFrom stats coef extractAIC formula glm poisson poly predict step update
#' @export
#' @examples
#' # Model search with quadratic effect of year
#'
#' search <- model_search(x = westeros, lists = 1:4, year = "Y", degree_year = 2)
#'
#' # result is equivalent for x = as.data.frame(table(westeros)))



model_search <- function(x, lists = NULL, k_max = 1, year = NULL, degree_year = NULL){


  if(is.null(lists))stop("unspecified argument `lists`")
  if(length(lists) < 2)stop("`lists` should have at least two elements")

  if ("Freq" %in% colnames(x)){
    nr       <- which(colnames(x)=="Freq")
    x[, -nr] <- lapply(x[, -nr], as.factor)
    lev      <- sapply(x[, -nr], levels)
    if(nrow(x) != prod(sapply(lev, length)))stop("`x` has incorrect number of rows")
  } else {
    x <- lapply(x, as.factor)
    x <- as.data.frame(table(x))
  }

  lev <- sapply(x[, -ncol(x)], levels)

  if(!(max(as.numeric(unlist(lev[lists]))) == 1 &
       min(as.numeric(unlist(lev[lists]))) == 0))stop("the lists should all be coded 0 and 1")

  d <- x

  if (!is.null(year) & !is.null(degree_year)){

    if(length(levels(x[, paste(year)])) <= degree_year){
      stop("`degree_year` equals/exceeds the number of years")
    }

    x[, paste(year)] <- poly(as.numeric(x[, paste(year)]), degree=degree_year, simple = T)
  }

  n   <- sum(x$Freq)
  obs <- ifelse(rowSums(x[, lists] == "1") == 0, 0, 1)
  m0  <- glm(Freq ~ ., poisson, data=subset(x, obs==1))
  mx  <- glm(paste("Freq ~ . ^", k_max + 1), poisson, data=subset(x, obs==1))

  m_bic <- step(m0, scope = list(lower = formula(m0), upper = formula(mx)), k=log(n), trace=0)
  m_aic <- step(m_bic, scope = list(lower = formula(m0), upper = formula(mx)), trace=0)
  m_nxt <- step(m_aic, scope = list(lower = formula(m_aic), upper = formula(mx)), k=1, trace=0)

  B         <- data.frame(m_bic$anova,
                          AICc = m_bic$anova$AIC,
                          BIC  = m_bic$anova$AIC,
                          Nhat = rep(0,nrow(m_bic$anova)))
  B$Nhat[1] <- sum(predict(m0, newdata = x, type="response"))
  B$AIC[1]  <- extractAIC(m0)[2]
  B$AICc[1] <- extractAIC(m0)[2] + 2*(length(coef(m0))^2 + length(coef(m0)))/(n - length(coef(m0)) -1)
  f         <- formula(m0)
  fx        <- list(f)
  cx        <- list(coef(m0))
  nx        <- list(predict(m0, newdata = x, type="response"))

  if(nrow(B) > 1){
    for (j in 2:nrow(B)){
      f         <- update(f, paste("~.", B$Step[j]))
      m         <- glm(f, poisson, data=subset(x, obs==1))
      nxx       <- predict(m, newdata=x, type="response")
      B$AIC[j]  <- extractAIC(m, k=2)[2]
      B$AICc[j] <- extractAIC(m, k=2)[2] + 2*(length(m$coef)^2 + length(m$coef))/(n - length(m$coef) - 1)
      B$Nhat[j] <- sum(nxx)
      fx[[length(fx) + 1]] <- f
      cx[[length(cx) + 1]] <- coef(m)
      nx[[length(nx) + 1]] <- nxx
    }
  }
  if (nrow(m_aic$anova) > 1){
    A         <- data.frame(m_aic$anova, AICc = m_aic$anova$AIC, BIC = m_aic$anova$AIC, Nhat=rep(0,nrow(m_aic$anova)))[-1, ]
    for (j in 1:nrow(A)){
      f         <- update(f, paste("~.", A$Step[j]))
      m         <- glm(f, poisson, data=subset(x, obs==1))
      nxx       <- predict(m, newdata=x, type="response")
      A$BIC[j]  <- extractAIC(m, k=log(n))[2]
      A$AICc[j] <- extractAIC(m, k=2)[2] + 2*(length(m$coef)^2 + length(m$coef))/(n - length(m$coef) - 1)
      A$Nhat[j] <- sum(nxx)
      fx[[length(fx) + 1]] <- f
      cx[[length(cx) + 1]] <- coef(m)
      nx[[length(nx) + 1]] <- nxx
    }
    B <- rbind(B, A)
  }
  if (nrow(m_nxt$anova) > 1){
    A         <- data.frame(m_nxt$anova, AICc = m_nxt$anova$AIC, BIC = m_nxt$anova$AIC, Nhat=rep(0,nrow(m_nxt$anova)))[-1, ]
    for (j in 1:nrow(A)){
      f         <- update(f, paste("~.", A$Step[j]))
      m         <- glm(f, poisson, data=subset(x, obs==1))
      nxx       <- predict(m, newdata=x, type="response")
      A$AIC[j]  <- extractAIC(m, k=2)[2]
      A$AICc[j] <- extractAIC(m, k=2)[2] + 2*(length(m$coef)^2 + length(m$coef))/(n - length(m$coef) - 1)
      A$BIC[j]  <- extractAIC(m, k=log(n))[2]
      A$Nhat[j] <- sum(nxx)
      fx[[length(fx) + 1]] <- f
      cx[[length(cx) + 1]] <- coef(m)
      nx[[length(nx) + 1]] <- nxx
    }
    B <- rbind(B, A)

  }

  B[, -1]     <- round(B[, -1], 1)
  min_aic     <- min(B$AIC)
  min_aicc    <- min(B$AICc)
  min_bic     <- min(B$BIC)
  B$AIC       <- B$AIC - min_aic
  B$BIC       <- B$BIC - min_bic
  B$AICc      <- B$AICc - min_aicc
  rownames(B) <- 1:nrow(B)

  options(scipen=999)
  print(B)

  invisible(list(models     = B,
                 formulas   = fx,
                 coefs      = cx,
                 fits       = nx,
                 minima     = c(AIC = min_aic, AICc = min_aicc, BIC = min_bic),
                 d          = d,
                 d_year     = x,
                 obs        = obs,
                 lists      = lists,
                 year       = year))
}



#' Estimation of 95\% bootstrap confidence intervals.
#'
#' @description \code{boot_mse} performs a parametric bootstrap for a selected model.
#' @param object an object created with the function \code{\link{model_search}}.
#' @param modelnr the number of the model in \code{object} that is to be bootstrapped.
#' If left unspecified, the  model with the lowest BIC will be selected.
#' @param rep the number of bootstrap replications.
#' @param seed a scalar to set the random seed.
#' @return A lists with tables and, if the variable \code{year} is specified, with plots.
#' \item{tables}{table(s) with the number of observations, the population size estimate with
#' the 95\% bootstrap confidence interval. If covariates are present, the tables are presented
#' for the levels of the covariate and, with more than one covariae, for the combinations of
#' the levels of pairs of covariates. If the variable \code{year} is specified, these
#' tables are presented by year.}
#' \item{plots}{plots of the population size estimates and their 95\% confidence intervals,
#' split out by the levels of the covariates and with the variable \code{year} on the x-axis.
#' Only available if the variable \code{year} is specified.}
#' @import doParallel
#' @import ggplot2
#' @importFrom foreach %dopar%
#' @importFrom iterators icount
#' @importFrom stats quantile rmultinom xtabs
#' @export


boot_mse <- function(object, modelnr = NULL, rep = 1000, seed = 1){

  if(is.null(modelnr)){
    modelnr <- which.min(object$models$BIC)
  }

  cl <- parallel::makeCluster(parallel::detectCores()/2)

  doParallel::registerDoParallel(cl)

  RNGkind("L'Ecuyer-CMRG")

  parallel::clusterSetRNGStream(cl = cl, iseed=seed)
  boot_fitted <- foreach::foreach(icount(rep), .combine=cbind) %dopar% {
    bootx(b0  = object$coefs[[modelnr]],
          f   = object$formulas[[modelnr]],
          N   = object$models$Nhat[[modelnr]],
          f0  = object$fits[[modelnr]],
          d   = object$d_year,
          obs = object$obs)
  }

  parallel::stopCluster(cl)

  invisible(
    generate_tables(d     = object$d,
                    lists = object$lists,
                    m     = object$fits[[modelnr]],
                    year  = object$year,
                    boot_fitted = boot_fitted)
  )

}



