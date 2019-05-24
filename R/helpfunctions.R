
bootx <- function(b0, f, N, f0, d, obs){

  d$Freq <- rmultinom(n=1, size=N, prob=f0/N)
  mx     <- glm(f, poisson, data=subset(d, obs==1), start=b0)
  freqs  <- predict(mx, newdata = d, type="response")

  return(freqs)

}

ci_95 <- function(x){
  v <- quantile(x, probs = c(0.025, 0.975))
  names(v) <- c("min95", "max95")
  v
}


gen_plot <- function(temp, ttl , gr = NULL){

  temp <- as.data.frame(temp)

  limx <- c(min(temp[, 1]), max(temp[, 1]))
  labx <- unique(temp[, 1])
  if(is.null(gr)){
    p <- ggplot2::ggplot(temp, aes(x = temp$Year, y = temp$Nhat))+
      geom_line() +
      geom_ribbon(aes(x = temp$Year, ymin = temp$min95, ymax = temp$max95), linetype = 0, alpha = .2) +
      ggtitle(ttl) +
      ylim(0, NA) +
      scale_x_continuous(name = "year", breaks = labx, labels = labx) +
      theme_light() +
      theme(legend.title = element_blank())
  }

  if(!is.null(gr)){
    p <- ggplot2::ggplot(temp, aes(x = temp$Year, y = temp$Nhat, group = gr, col = gr))+
      geom_line() +
      geom_ribbon(aes(x = temp$Year, ymin = temp$min95, ymax = temp$max95, fill = gr), linetype = 0, alpha = .2) +
      ggtitle(ttl) +
      ylim(0, NA) +
      scale_x_continuous(name = "year", breaks = labx, labels = labx) +
      theme_light() +
      theme(legend.title=element_blank())
  }

  p
}



generate_tables <- function(d, lists, m, year, boot_fitted){

  est       <- list()
  est[[1]]  <- round(cbind(nobs = sum(d$Freq),
                            Nhat = sum(m),
                            t(ci_95(colSums(boot_fitted)))))
  names(est)[[1]]    <- "x"
  rownames(est[[1]]) <- "all"

  if(!is.null(year)){

    plots <- list()

    yearnr  <- which(colnames(d) == year)
    temp    <- matrix(0, length(levels(d[, yearnr])), 5,
                      dimnames = list(levels(d[, yearnr]),
                                      c("Year", "nobs", "Nhat", "min95", "max95")))

    for(j in 1:length(levels(d[, yearnr]))){

      nrs        <- which(d[, yearnr] == levels(d[, yearnr])[j])
      temp[j, ]  <- round(cbind(Year = as.numeric(levels(d[, yearnr])[j]),
                                nobs = sum(d$Freq[nrs]),
                                Nhat = sum(m[nrs]),
                                t(ci_95(colSums(boot_fitted[nrs, ])))))
    }

    est[[length(est) + 1]]    <- temp
    plots[[1]]                <- gen_plot(temp = temp, ttl = names(est)[length(est)])
    names(est)[[length(est)]] <- names(plots)[[1]] <- colnames(d)[yearnr]

    }




  if(ncol(d) > length(lists) + 1 & is.null(year)){    #no year and only one covariate

    covs <- d[, -lists]
    lev  <- lapply(covs[, -ncol(covs), drop = F], levels)

    for(i in 1:length(lev)){
      temp <- matrix(0, length(lev[[i]]), 4,
                     dimnames = list(lev[[i]],
                                     c("nobs", "Nhat", "min95", "max95")))
      for(j in 1:length(lev[[i]])){
        nrs <- which(covs[, i] == lev[[i]][j])
        temp[j, ]  <- round(cbind(nobs = sum(covs$Freq[nrs]),
                                  Nhat = sum(m[nrs]),
                                  t(ci_95(colSums(boot_fitted[nrs, ])))))
      }
      est[[length(est) + 1]]     <- temp
      names(est)[[length(est)]]  <- names(lev)[i]
    }


  }

  if(ncol(d) > length(lists) + 2 & is.null(year)){     #no year and at least two covariates


      for(i in 1:(length(lev) - 1)){
        for(j in (i+1):length(lev)){

          z         <- xtabs(boot_fitted ~ ., covs[, c(i, j)])
          temp      <- NULL
          tempnames <- NULL

          for(k in 1:length(lev[[i]])){
            for(l in 1:length(lev[[j]])){

              nrs <- which(covs[, i] == lev[[i]][k] & covs[, j] == lev[[j]][l])

              vtemp <- round(cbind(nobs = sum(covs$Freq[nrs]),
                                   Nhat = sum(m[nrs]),
                                   t(ci_95(colSums(boot_fitted[nrs, ])))))
              temp <- rbind(temp, vtemp)
              tempnames <- c(tempnames, paste(lev[[i]][k], lev[[j]][l], sep = ":"))

            }
          }
          rownames(temp)            <- tempnames
          est[[length(est) + 1]]    <- temp
          names(est)[[length(est)]] <-  paste(names(lev)[i], names(lev)[j], sep = "x")
        }
      }
    }

  if(ncol(d) > length(lists) + 2 & !is.null(year)){       #year and at leat one covariate

    covs    <- d[, -lists]
    lev     <- lapply(covs[, -ncol(covs), drop = F], levels)

    yearnr  <- which(colnames(covs) == year)
    covnr   <- (1:ncol(covs))[-c(yearnr, ncol(covs))]

    temp   <- matrix(0, length(lev[[yearnr]]), 5,
                   dimnames = list(lev[[yearnr]],
                                   c("Year", "nobs", "Nhat", "min95", "max95")))


    for(i in covnr){

      temp      <- NULL
      tempnames <- NULL

      for(j in 1:length(lev[[i]])){
        for(k in 1:length(lev[[yearnr]])){
          nrs        <- which(covs[, i] == lev[[i]][j]  & covs[, yearnr] == lev[[yearnr]][k])
          vtemp      <- round(cbind(Year = as.numeric(lev[[yearnr]][k]),
                                    nobs = sum(covs$Freq[nrs]),
                                    Nhat = sum(m[nrs]),
                                    t(ci_95(colSums(boot_fitted[nrs, ])))))
          temp       <- rbind(temp, vtemp)
          tempnames  <- c(tempnames, paste(lev[[i]][j]))
        }
      }
      rownames(temp)                <- tempnames
      est[[length(est) + 1]]        <- temp
      names(est)[[length(est)]]     <- names(lev)[i]
      plots[[length(plots) + 1]]    <- gen_plot(temp = temp, ttl = names(est)[length(est)], gr = as.factor(rownames(temp)))
      names(plots)[[length(plots)]] <- names(lev)[i]
    }
  }

  if(ncol(d) > length(lists) + 3 & !is.null(year)){       #year and at leat two covariates

    yearnr  <- which(colnames(covs) == year)
    covnr   <- (1:ncol(covs))[-c(yearnr, ncol(covs))]


    for(i in covnr[-length(covnr)]){
      for(j in (i + 1):length(covnr)){

        temp      <- NULL
        tempnames <- NULL

        for(k in 1:length(lev[[i]])){
          for(l in 1:length(lev[[j]])){
            for(y in 1:length(lev[[yearnr]])){

              nrs <- which(covs[, i] == lev[[i]][k] & covs[, j] == lev[[j]][l] & covs[, yearnr]  == lev[[yearnr]][y])

              vtemp <- round(cbind(Year = as.numeric(lev[[yearnr]][y]),
                                   nobs = sum(covs$Freq[nrs]),
                                   Nhat = sum(m[nrs]),
                                   t(ci_95(colSums(boot_fitted[nrs, ])))))
              temp      <- rbind(temp, vtemp)
              tempnames <- c(tempnames, paste(lev[[i]][k], lev[[j]][l], sep = ":"))
            }
          }
        }
        rownames(temp)                <- tempnames
        est[[length(est) + 1]]        <- temp
        names(est)[[length(est)]]     <- paste(names(lev)[i], names(lev)[j], sep = "x")
        plots[[length(plots) + 1]]    <- gen_plot(temp = temp, ttl = names(est)[length(est)], gr = as.factor(rownames(temp)))
        names(plots)[[length(plots)]] <- paste(names(lev)[i], names(lev)[j], sep = "x")
      }
    }
  }

  if(is.null(year)){

    return(est)

  }else{

    return(list(tables = est,
                plots = plots))
  }

}
