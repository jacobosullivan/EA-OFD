## This script will numerically simulate the LSPOM model to demonstrate the convergence on 
## a right-skewed time-invariant steady state OFD

simulateLSPOM <- function(N=50, n=10, m=10, invasions=1e4, storeS=T, storeO=F, neutral=T) {
  
  # List of N vectors of length n:
  occupancy_list <- apply(t(replicate(N,as.array(1:n))),1,as.list)

  # N by n array recordSing species' integer IDs
  occupancy_table <- t(replicate(N,as.array(1:n)))

  if (storeO) {
    storeS=T
  }
  
  # Storage object for intermediate richness
  if (storeS) {
    recordS <- data.frame(i=1:invasions,S=NA)  
  }
  
  # Storage object for intermediate OFDs
  if (storeO) {
    OFD_mat <- matrix(0, nrow=invasions, ncol=N)
  }
  
  # Initialise counter for integer species IDs
  last_i <- length(unique(as.vector(occupancy_table)))
  
  for(i in 1:invasions){
    S <- length(unique(as.vector(occupancy_table)))
    if (storeS) {
      recordS$S[i] <- S
    }
    new_i <- last_i+1
    
    ## Regional invasion + internal colonisation by new species
    # Sample random colonised patch
    invaded_patch <- sample.int(nrow(occupancy_table),1)
    
    # Sample additional colonisations
    if (neutral) {
      invade_patches <- rep(FALSE, nrow(occupancy_table))
      invade_patches[invaded_patch] <- TRUE  
    } else {
      invade_other_probability <- n/S # probability of additional colonisation
      invade_patches <- as.logical(rbinom(nrow(occupancy_table),1,invade_other_probability))
      invade_patches[invaded_patch] <- TRUE
    }
    
    # Update occupancy table
    for(patch in which(invade_patches)){
      occupancy_table[patch,sample.int(n,1)] <- new_i 
    }
    
    ## Internal colonisation by resident species
    # Sample species whose range expands
    mixers <- occupancy_table[rbinom(length(occupancy_table), 1, m/length(occupancy_table))==1]
    
    # Update occupancy table
    for (j in mixers) {
      target_patch <- sample.int(nrow(occupancy_table),1)
      if(j %in% occupancy_table[target_patch,]) next;
      occupancy_table[target_patch,sample.int(n,1)] <- j
    }
    
    if (storeO) {
      OFD_i <- table(table(as.numeric(as.factor(occupancy_table))))
      OFD_mat[i,as.numeric(names(OFD_i))] <- OFD_i
    }
    
    last_i <- new_i
  }
  
  OFD_table <- table(table(as.numeric(as.factor(occupancy_table))))
  OFD <- data.frame(site.occ=1:max(as.numeric(names(OFD_table))),
                    freq=0)
  OFD$freq[as.numeric(names(OFD_table))] <- OFD_table
  
  if ((storeS) && (storeO)) {
    return(list(S=recordS, OFD_mat=OFD_mat))
  } else if (storeS) {
    return(list(S=recordS, OFD=OFD))
  } else {
    return(OFD)
  }
}

N <- 50 # set number of sites
n <- 10 # set richness limit for all sites
invasions <- 1e4 # set number of invasions

par(mfrow=c(3,2))

m <- 1 # set mixing rate
res <- simulateLSPOM(N=N, n=n, m=m, invasions=invasions, storeS=T)
plot(res$S, type='l', xlab="Invasion", ylab="Richness")
plot(res$OFD, xlab="Site occupancy", ylab="Number of species", type='l')

m <- 5 # set mixing rate
res <- simulateLSPOM(N=N, n=n, m=m, invasions=invasions, storeS=T)
plot(res$S, type='l', xlab="Invasion", ylab="Richness")
plot(res$OFD, xlab="Site occupancy", ylab="Number of species", type='l')

m <- 10 # set mixing rate
res <- simulateLSPOM(N=N, n=n, m=m, invasions=invasions, storeS=T)
plot(res$S, type='l', xlab="Invasion", ylab="Richness")
plot(res$OFD, xlab="Site occupancy", ylab="Number of species", type='l')

if (0) {
  # Run in parallel, high CPU machine recommended
  require(foreach)
  require(doParallel)
  require(tidyverse)
  
  cores=detectCores()  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  N <- 50 # set number of sites
  n <- 10 # set richness limit for all sites
  invasions <- 1e4 # set number of invasions
  no_reps <- 1000
  
  res <- c()
  m <- c(2.5,5,10)
  
  for (mm in m) {
    res_m <- foreach (r=1:no_reps, .combine=rbind) %dopar% { 
      res_r <- simulateLSPOM(N=N, n=n, m=mm, invasions=invasions, storeS=F)
      res_r$m <- mm
      res_r
    }
    res <- rbind(res, res_m)
  }
  stopCluster(cl)
  
  write.csv(res, "/data/home/btx718/Paleolimnology/datasets/EA/summary_data/LSPOM_res_neutral.csv", row.names=F)
  
  res.mn <- res %>%
    group_by(m, site.occ) %>%
    summarise(freq.mn=mean(freq), freq.sd=sd(freq))
  
  write.csv(res.mn, "/data/home/btx718/Paleolimnology/datasets/EA/summary_data/LSPOM_res_neutral_mean.csv", row.names=F)
}

if (0) {
  # Test stationarity in mean, variance and skewness of OFD for simulated metacommunity
  # Generates LARGE matrix object (dimensions: invasion X N)
  
  momentExp <- function(OFD_i) {
    OFD_tab <- data.frame(x=which(OFD_i!=0), Freq=OFD_i[which(OFD_i!=0)])
    OFD <- apply(OFD_tab, MARGIN=1, FUN=function(X) rep(X[1], X[2]))
    if (class(OFD)[1]=="list") {
      OFD <- unlist(OFD)
    } else if ((class(OFD)[1]=="matrix")) {
      OFD <-as.vector(OFD)
    }
    mu1 <- mean(OFD)
    mu2 <- var(OFD)
    mu3 <- skewness(OFD)
    return(c(mu1=mu1,mu2=mu2,mu3=mu3))
  }
  
  require(e1071)
  require(tseries)
  
  N <- 50 # set number of sites
  n <- 10 # set richness limit for all sites
  invasions <- 1e4 # set number of invasions
  
  m <- 1 # set mixing rate
  res <- simulateLSPOM(N=N, n=n, m=m, invasions=invasions, storeO=T)
  
  par(mfrow=c(1,1))
  plot(res$S, type='l', xlab="Invasion", ylab="Richness")
  
  dat_momentExp <- as.data.frame(t(apply(res$OFD_mat, MAR=1, FUN=momentExp)))
  dat_momentExp$i <- 1:nrow(dat_momentExp)
  
  par(mfrow=c(3,1))
  plot(mu1 ~ i, dat_momentExp, type='l')
  plot(mu2 ~ i, dat_momentExp, type='l')
  plot(mu3 ~ i, dat_momentExp, type='l')
  
  thresh <- 5000
  
  plot(mu1 ~ i, dat_momentExp[(nrow(dat_momentExp)-thresh):nrow(dat_momentExp),], type='l')
  plot(mu2 ~ i, dat_momentExp[(nrow(dat_momentExp)-thresh):nrow(dat_momentExp),], type='l')
  plot(mu3 ~ i, dat_momentExp[(nrow(dat_momentExp)-thresh):nrow(dat_momentExp),], type='l')
  
  adf.test(dat_momentExp$mu1[(nrow(dat_momentExp)-thresh):nrow(dat_momentExp)])
  adf.test(dat_momentExp$mu2[(nrow(dat_momentExp)-thresh):nrow(dat_momentExp)])
  adf.test(dat_momentExp$mu3[(nrow(dat_momentExp)-thresh):nrow(dat_momentExp)])
  
}

