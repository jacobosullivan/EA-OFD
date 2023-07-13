fitTimeScales <- function(site.occ_array, catchment, x0) {
  ## This function takes as input a species-site-time array with the dimnames(site.occ_array)[3] a vector of time-point names (e.g. years)
  ## It will return range estimates of d, popDetect, alpha, al, metapopDetect, pRet and pRet.p - the model estimates
  ## Not all datasets are amenable to analysis using these fitting functions for example if the decay in compositional similarity is cannot be fit by an exponential model.
  ## Here we skip such catchments; updates are possible to model these differently
  ## pRet and pRet.p are not used in the current analysis

  ## Build a data set of pairs of years:
  time_spans_full <-
    data.frame(y1=rep(as.numeric(dimnames(site.occ_array)$year),times=dim(site.occ_array)[3]),
               y2=rep(as.numeric(dimnames(site.occ_array)$year),each=dim(site.occ_array)[3]),
               span=NA, retained=NA, weight=NA, retained2 = NA, weight2=NA)

  ## Compute time spans, similarity metrics of weights for data:
  for(y in seq(nrow(time_spans_full))) {
    y1 <- as.character(time_spans_full$y1[y])
    y2 <- as.character(time_spans_full$y2[y])
    keepAllSites <- F
    sites <-  ( colSums(site.occ_array[,,y1]) & colSums(site.occ_array[,,y2]) ) | keepAllSites
    time_spans_full$retained[y] <- sum(site.occ_array[,sites,y1] & site.occ_array[,sites,y2])/sum(site.occ_array[,sites,y1])
    time_spans_full$retained2[y] <- sum(rowSums(site.occ_array[,,y1]) & rowSums(site.occ_array[,,y2]))/sum(rowSums(site.occ_array[,,y1])>0)
    time_spans_full$weight[y] <- sqrt(sum(site.occ_array[,sites,y1])*sum(site.occ_array[,sites,y2]))
    if(time_spans_full$weight[y]==0) time_spans_full$retained[y] <- 0
    time_spans_full$weight2[y] <- sqrt(sum(rowSums(site.occ_array[,,y1])>0)*sum(rowSums(site.occ_array[,,y2])>0))
    time_spans_full$span[y] <- abs(diff(as.numeric(c(y1,y2))))  ## !! try both, signed values and abs() of this
  }

  time_spans_full <- time_spans_full[time_spans_full$span!=0,]

  fit_table <- data.frame(d=rep(NA,dim(site.occ_array)[3]),popDetect=NA,alpha=NA,al=NA,metapopDetect=NA, pRet=NA, pRet.p=NA)

  alpha <- colSums(site.occ_array > 0)# gives number of sampled species for each site and year
  if(!keepAllSites) alpha[alpha==0] <- NA # site-years without species are not sampled, discard
  alphas <- colMeans(alpha,na.rm = T) # gives mean number of species per site
  gammas <- colSums(apply(site.occ_array,MARGIN = c(1,3),sum) > 0)
  mean_a_pRet <- mean(alphas/gammas,na.rm = T)

  for(y in seq(1+dim(site.occ_array)[3])) { # leave each year out once, and finally do all years
    dropped_year <- as.character(dimnames(site.occ_array)$year[y])
    years <- dimnames(site.occ_array)$year[-y]
    year_names <- as.character(years)

    if(is.na(dropped_year)){ # do all years
      time_spans <- time_spans_full
    } else { # drop one year
      time_spans <- time_spans_full[with(time_spans_full,y1 != dropped_year & y2 != dropped_year),]
    }

    if (cor.test(time_spans$retained, time_spans$span)$p.value > 0.05) {
      print(paste("y =", y, "Fail cor.test"))
    }

    plotit <- F

    ## Fit turnover rate d of local communities:
    if(plotit){
      pdf(file = paste0("/tmp/plots/",gsub(" ","_",catchment),".pdf"),width = 5,height = 8)
      par(mfrow = c(2,1))
      sel <- time_spans$weight > 0
      with(time_spans[sel,],plot(span,retained,pch=16,ylim=c(0,1),col=gray(1-weight/max(weight))))
      with(time_spans[sel,],points(span,retained),col=gray(0.9,0.1),lwd=0.2)
      title(catchment)
    }
    mod <- lm(retained ~ span,data=time_spans[time_spans$span>0,],weights = weight)
    if(plotit) lines(abline(mod))
    fit <-
      tryCatch(
        nls(retained ~ (a-mean_a_pRet)*exp(-d*span) + mean_a_pRet, weights = weight, data = time_spans[time_spans$span>0,],
               start = c(a=as.numeric(mod$coefficients[1]),
                         d=-as.numeric(mod$coefficients[2])),
               control=nls.control(maxiter=5000, minFactor = 1e-200)),
           error = function(x) NA)

    if (is.na(fit) || (fit$m$getPars()["d"] < 0)) {
      print(paste("y =", y, "Fail fit()"))
      if(plotit) dev.off()
      next
    }
    if(plotit){
      lines(predict(fit,newdata=data.frame(span=seq(1,20))),col="red")
      legend(x="topright", col=c("black","red"),legend = c("linear","exponential"),lty=1,bt="n")
    }
    fit_table[y,"popDetect"] <- as.list(fit$m$getPars()["a"])
    fit_table[y,"d"] <- with(as.list(fit$m$getPars()),d*(1-mean_a_pRet/a))

    pRet <- mean_a_pRet
    pRet.p <- 1
    fit_table[y,c("pRet","pRet.p")] <- list(pRet=pRet, pRet.p=pRet.p)

    if(1){
      if(plotit){
        with(time_spans,plot(span,retained2,pch=16,ylim=c(0,1),col=gray(1-weight/max(weight))))
        with(time_spans,points(span,retained2),col=gray(0.9,0.1))
      }

      d_to_use <- fit_table[y,"d"]
      with_q <- F
      q <- 1
      if(with_q){
        PReObs <- function(t,m,q=1){
          -(log(1 - (m*(1 + m)*q^2)/(exp(t/(1 + m))*(1 + m*q)^2))/
              log(1 + m*q))
        }
        PReObsX <- function(t,m,q=1){
          PReObs(t,m,q)+exp(-t)*(1-PReObs(0,m,q))
        }
        fit2 <- tryCatch(
          with(list(d=d_to_use),
               nls(retained2 ~ a*PReObsX(t=span*d,m=pnorm(al)/(1-pnorm(al))/q,q=q),
                   weights = weight,  ## try both: weight and weight2
                   data = time_spans[with(time_spans,span>0),],
                   start = c(a=0.5,al=qnorm(0.5)),
                   control=nls.control(maxiter=500, minFactor = 1e-200))),
               error = function(x) NA)
      }else{# not with_q
        fit2 <- tryCatch(
          with(list(d=d_to_use),
               nls(retained2 ~ a*(log(1-pnorm(al)*exp((pnorm(al)-1)*d*span))/log(1-pnorm(al))),
                   weights = weight,  ## try both: weight and weight2
                   data = time_spans[with(time_spans,span>0),],
                   start = c(a=0.5,al=qnorm(0.5)),
                   control=nls.control(maxiter=500, minFactor = 1e-200))),
          error = function(x) NA)
      }

      if (is.na(fit2)) {
        print(paste("y =", y, "Fail fit2()"))
        if(plotit) dev.off()
        next
      }
      if(with_q){
        print(fit2$m$getPars())
        fit_table[y,c("alpha")] <- list(alpha=pnorm(fit2$m$getPars()["al"]))
        fit_table[y,c("al","metapopDetect")] <- as.list(fit2$m$getPars()[c("al","a")])
      }else{
        fit_table[y,c("alpha")] <- list(alpha=pnorm(fit2$m$getPars()["al"]))
        fit_table[y,c("al","metapopDetect")] <- as.list(fit2$m$getPars()[c("al","a")])
      }
      if(plotit){
        with(as.list(fit$m$getPars()["d"]),
             lines(predict(fit2,newdata=data.frame(span=seq(1,20))),col="red") )
        dev.off()
      }
    }

  }

  if (nrow(fit_table) < x0) {
    res <- "Filtered following fitting"
  } else {
    ## Evaluate using jackknife:
    fit_table0 <- fit_table[nrow(fit_table),]  # these are the directly computed means
    fit_table1 <- fit_table[-nrow(fit_table),] # these are the leave-one-out means
    n <- nrow(fit_table1)
    covar <- (n-1)/n * var(fit_table1,na.rm = T)
    biasCorrectedMean <- n*fit_table0 - (n-1)*colMeans(fit_table1)
    ranges <- outer(sqrt(diag(covar)),c(0,qnorm(c(0.025,0.975))))+as.double(biasCorrectedMean) # should be better, but it's not
    ranges <- outer(sqrt(diag(covar)),c(0,qnorm(c(0.025,0.975))))+as.double(fit_table0)
    colnames(ranges) <- c("mean","lower","upper")
    ranges[,"mean"] <- as.double(fit_table0)
    rownames(ranges) <- colnames(fit_table)
    al2m <- function(al){pnorm(al)/(1-pnorm(al))}
    ranges <- rbind(ranges,al2m(ranges["al",]))
    rownames(ranges)[nrow(ranges)] <- "m"
    if (!class(fit2)=="logical") {
      max_diff=unname(predict(fit2, newdata=data.frame(span=max(time_spans$span))))
      min_diff=unname(predict(fit2, newdata=data.frame(span=0)))
    } else {
      max_diff=NA
      min_diff=NA
    }
    
    res <- list(catchment=catchment, 
                ranges=ranges, 
                covar=covar, 
                core=data.frame(pRet=fit_table$pRet, pRet.p=fit_table$pRet.p),
                max_diff=max_diff,
                min_diff=min_diff)
  }

  return(res)
}
