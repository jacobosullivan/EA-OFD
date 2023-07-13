require(data.table)

regularizedIncompleteGamma <- function(a,z){ # the regularized incomplete gamma function Q(a,z)
  pgamma(1,rate = z,shape = a, lower.tail = FALSE)
}

plainNeutralTheory <- TRUE

nicheNumberDist <- function(lambda,mixing=0.0001,nMax=15){
  S <- rep(NA,nMax)
  S[1] <- 1/(lambda+mixing)
  expFac <- exp((lambda-1)*lambda/mixing)
  mixFac <- lambda/mixing +1
  mixProp <- mixing/(lambda+mixing)
  mixX <- (lambda-1)/mixProp
  if(plainNeutralTheory){
    if(lambda!=1)stop("lambda must be 1")
    for(n in 2:nMax){
      S[n] <- (mixing/(1+mixing))^n/mixing/n
    }
  }else{
    for(n in 2:nMax){
      S[n] <- 1/(n*lambda)*
        ( 1-(mixing/(lambda+mixing))^n -
            regularizedIncompleteGamma(n,lambda-1) +
            expFac * mixFac * mixProp^n * regularizedIncompleteGamma(n,mixX)
        )
    }
  }
  return(S)
}

consistentLambda <- function(mixing=0.0001,nMax=15){
  if(plainNeutralTheory){
    return(1)
  }else{
    l <-
      uniroot(function(l)sum(nicheNumberDist(l,mixing,nMax))*(l-1)-1,
              lower=1,upper=2,extendInt="upX")
    return(l$root)
  }
}

fitMixingSS <- function(data){
  ## simple sum of squared differences
  data <- data/sum(data*(1:length(data)))
  goalFunction <- function(m){ 
    l <- consistentLambda(m,length(data))
    nND <- nicheNumberDist(l,m,length(data))
    nND <- nND/sum(nND*(1:length(data)))
    sum(abs(data - nND)^2)
  }
  o <- optim(3,goalFunction,method="Brent",lower=0.006,upper = 100)
  return(o)
}

fitMixingML <- function(data){  
  ## maximum likelihood fit
  if(any(!(as.integer(as.matrix(data))==as.vector(as.matrix(data))))) stop("Input to fitMixingML must be array of integer counts.")
  patchSum <- sum(data*(1:length(data)))
  goalFunction <- function(m){
    l <- consistentLambda(m,length(data))
    nND <- nicheNumberDist(l,m,length(data))
    nND <- nND*(patchSum/sum(nND*(1:length(data))))
    -sum(-nND + data*log(nND)-lgamma(data+1))
  }
  o <- optim(3,goalFunction,method="Brent",lower=0.006,upper = 100)
  l <- consistentLambda(o$par,length(data))
  nND <- nicheNumberDist(l,o$par,length(data)) * patchSum
  # Compute G2 statistic with grouping of counts according to Wood (2002): 
  # https://doi.org/10.1081/STA-120015014
  c <- 1.6008
  group_size <- 0
  group_mu <- 0
  group_n <- 0
  G2Sum <- 0
  degrees_of_freedom <- -2
  for(i in 1:length(data)){
    mu <- nND[i] 
    minimumGroupSize <- (c - mu + sqrt(c^2+mu^2))/(2*mu) 
    group_size <- group_size + 1
    group_mu <- group_mu + mu
    group_n <- group_n + data[i]
    finalG2Sum <- NULL
    if(group_size > minimumGroupSize){
      # G2 Sum in case the next group does not finalize:
      if(i < length(data)){
        final_group_n <- group_n + sum(data[i:length(data)])
        final_group_mu <- group_mu + sum(nND[i:length(data)])
      }else{
        final_group_n <- group_n 
        final_group_mu <- group_mu
      }
      finalG2Sum <- G2Sum + 
        2 * ( ifelse( final_group_n>0, final_group_n*log(final_group_n/final_group_mu), 0) 
              - final_group_n + final_group_mu )
      # finalise current group
      G2Sum <- G2Sum +
        2 * ( ifelse( group_n>0, group_n*log(group_n/group_mu), 0) - group_n + group_mu )
      degrees_of_freedom <-
        degrees_of_freedom + 1
      group_size <- 0
      group_mu <- 0
      group_n <- 0
    }
  }
  if (!is.null(finalG2Sum)) {
    goodness_of_fit_stats <-
      list(G2=finalG2Sum,dof=degrees_of_freedom,
           p_fit = pchisq(finalG2Sum,degrees_of_freedom, lower.tail = FALSE))      
  } else {
    goodness_of_fit_stats <-
      list(G2=NA,dof=NA,p_fit = NA)      
  }

  return(c(o,goodness_of_fit_stats))
}

f <- function(site.occ) {
  # wrapper for the fitting procedure, expects either a binary species by site matrix
  # or a vector of length N with elements equal to the (average) occupancy of each site class
  
  if (length(dim(site.occ))==2) {
    specs <- data.frame(alpha=NA,alpha2=NA,gamma=NA)
    specs["alpha"] <- mean(colSums(site.occ==1));
    specs["alpha2"] <- mean(colSums(site.occ==1)^2);
    specs["gamma"] <- nrow(site.occ);  
    tt <- table(rowSums(site.occ==1)); ## rowSums gives a measure of species range size
    if ("0" %in% names(tt)) {
      tt <- tt[-1]  
    }
    ss <- rep(0, ncol(site.occ))
    ss[as.numeric(names(tt))] <- tt ## add zeros to patch occupancy distribution
    N <- ncol(site.occ)
  } else {
    N <- length(site.occ)
    # ss <- site.occ[1:max(which(site.occ>0))] # remove zeros - this may lead to under-estimation of nNiches...
    ss <- site.occ
  }
  
  # print(ss)
  
  names(ss) <- 1:length(ss)
  obs <- ss
  
  if (N > 1) {
    NMAX_ss <- length(ss)
    nNiches_ss <- sum((1:length(ss))*ss)
    mixFit_ss <- fitMixingSS(ss/nNiches_ss)
    mix_ss <- mixFit_ss$par
    pred_ss <- nicheNumberDist(consistentLambda(mix_ss,nMax = NMAX_ss),mix_ss,nMax = NMAX_ss)
    nNiches <- nNiches_ss
    
    res <- list(SOD=data.frame(site.occ = 1:N,
                               taxa.obs = c(obs, rep(0, N-length(obs))),
                               taxa.pred_ss = c(pred_ss, rep(0, N-length(pred_ss)))),
                estimate=mix_ss,
                nNiches=nNiches)
  } else {
    res <- list(SOD=data.frame(site.occ = 1,
                        taxa.obs = obs,
                        taxa.pred_ss = NA),
                 estimate=NA,
                 nNiches=NA)
  }
  
  return(res)
}
