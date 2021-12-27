fitTimeScales <- function(site.occ_array, x0) {
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
    time_spans_full$retained[y] <- sum(site.occ_array[,,y1] & site.occ_array[,,y2])/sum(site.occ_array[,,y1])
    time_spans_full$retained2[y] <- sum(rowSums(site.occ_array[,,y1]) & rowSums(site.occ_array[,,y2]))/sum(rowSums(site.occ_array[,,y1])>0)
    time_spans_full$weight[y] <- sqrt(sum(site.occ_array[,,y1])*sum(site.occ_array[,,y2]))
    time_spans_full$weight2[y] <- sqrt(sum(rowSums(site.occ_array[,,y1])>0)*sum(rowSums(site.occ_array[,,y2])>0))
    time_spans_full$span[y] <- abs(diff(as.numeric(c(y1,y2))))  ## !! try both, signed values and abs() of this
  }
  
  time_spans_full <- time_spans_full[time_spans_full$span!=0,]
  
  fit_table <- data.frame(d=rep(NA,dim(site.occ_array)[3]),popDetect=NA,alpha=NA,al=NA,metapopDetect=NA, pRet=NA, pRet.p=NA)
  
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
      next
    }
    
    ## Fit turnover rate d of local communities:
    fit <- tryCatch(
      nls(retained ~ a*exp(-d*span), weights = weight, data = time_spans[time_spans$span>0,], start = c(a=0.5,d=1/5)),
      error = function(x) NA)
    
    if (is.na(fit) || (fit$m$getPars()["d"] < 0)) {
      next
    }
    
    ## Fit community level turnover to determine alpha
    ## Least square fit
    ### The constraint that 0 < alpha < 1 is implemented by defining al = qnorm(alpha), alpha=pnorm(al)
    fit2 <- tryCatch(
      with(as.list(fit$m$getPars()["d"]),
           nls(retained2 ~ a*(1-log((1-pnorm(al))*exp((1-pnorm(al))*d*span)/(exp((1-pnorm(al))*d*span)-pnorm(al)))/log(1-pnorm(al))), 
               weights = weight,  ## try both: weight and weight2
               data = time_spans[with(time_spans,span>0),], 
               start = c(a=0.3,al=qnorm(0.9)) )),
      error = function(x) NA)
    
    if (is.na(fit2)) {
      next
    }
    
    if (0) {
      ## Estimate non-zero convergence of exponential decay: proportion of 'core' species - Not currently used in downstream analysis
      ## Least square fit
      ### The constraint that 0 < alpha < 1 is implemented by defining al = qnorm(alpha), alpha=pnorm(al)
      fit3 <- tryCatch(
        nls(retained ~ a*((1-p)*exp(-d*span)+p), weights = weight, data = time_spans[time_spans$span>0,], start = c(a=0.5,d=1/5,p=0.1)),
        error = function(x) NA)
      
      if (is.na(fit3)) {
        pRet <- pRet.p <- NA
      } else {
        pRet <- fit3$m$getPars()[3]
        pRet.p <- summary(fit3)$coefficients["p",4]
      }
    } else {
      pRet <- pRet.p <- 0
    }
    
    if (is.na(fit) || (fit$m$getPars()["d"] < 0)) {
      next
    }
    
    fit_table[y,] <- c(as.list(fit$m$getPars()[c("d","a")]),
                       list(alpha=pnorm(fit2$m$getPars()["al"])),
                       fit2$m$getPars()[c("al","a")],
                       pRet,
                       pRet.p)
  }
  
  # Early return for debugging:
  res <- fit_table
  return (res)
  
  fit_table <- fit_table[complete.cases(fit_table),]
  
  if (nrow(fit_table) < x0) {
    res <- NULL
  } else {
    ## Evaluate using jackknife:
    fit_table0 <- fit_table[nrow(fit_table),]  # these are the directly computed means
    fit_table1 <- fit_table[-nrow(fit_table),] # these are the leave-one-out means
    n <- nrow(fit_table1)
    covar <- (n-1)/n * var(fit_table1)
    biasCorrectedMean <- n*fit_table0 - (n-1)*colMeans(fit_table1)
    ranges <- outer(sqrt(diag(covar)),c(0,qnorm(c(0.025,0.975))))+as.double(biasCorrectedMean) # should be better, but it's not
    ranges <- outer(sqrt(diag(covar)),c(0,qnorm(c(0.025,0.975))))+as.double(fit_table0)
    colnames(ranges) <- c("mean","lower","upper")
    rownames(ranges) <- colnames(fit_table)
    al2m <- function(al){pnorm(al)/(1-pnorm(al))}
    ranges <- rbind(ranges,al2m(ranges["al",]))
    rownames(ranges)[nrow(ranges)] <- "m"
    res <- list(catchment=sort(unique(dat[,catch_index]))[i], ranges=ranges, covar=covar, core=data.frame(pRet=fit_table$pRet, pRet.p=fit_table$pRet.p))
  }
  
  # return(res)
}
