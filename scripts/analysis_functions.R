## Flexible analysis functions for EA dataset

require(tidyverse) # relational database management
require(Matrix) # sparse matrix ss tables
require(foreach) # filling large ss tables
require(doParallel) # filling large ss tables
require(e1071) # compute skewness
require(imputeTS) # impute missing data in time series
require(tseries) # adf and kpss tests
require(trend) # mann-kendall test and sen's slope

source("scripts/fitMCPD.R")

## Catchment-year geometry
area_dist <- function(dat, 
                      cores=NULL,
                      catch_max=NULL,
                      catch="management") {
  
  ## This function will, for each management catchment, record the total area and mean distance between samples
  if (is.null(cores)) {
    cores=detectCores()
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% { 
    
    require(tidyverse) # relational database management
    require(geometry) # compute convex hull
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      select(site,year,Northing, Easting)
    
    res <- c()
    for (y in sort(unique(dat_mc$year))) {
      
      dat_mc_y <- dat_mc %>%
        filter(year == y) %>%
        select(Easting, Northing)
      
      dat_mc_y <- dat_mc_y[!duplicated(dat_mc_y),]
      
      if (nrow(dat_mc_y)>3) {
        ch <- convhulln(subset(dat_mc_y, select=c(Easting, Northing)), output.options=TRUE)
        a <- ch$vol/1e6 # area of convex hull
        dsub <- as.matrix(dist(subset(dat_mc_y, select=c(Easting, Northing))))
        d <- mean(dsub[upper.tri(dsub)])/1e3
        res <- rbind(res, data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                                     year=y,
                                     area=a,
                                     dist.mn=d))
      } else if (nrow(dat_mc_y) == 2) {
        dsub <- dist(subset(dat_mc_y, select=c(Easting, Northing)))
        d <- mean(dsub[upper.tri(dsub)])/1e3
        res <- rbind(res, data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                                     year=y,
                                     area=NA,
                                     dist.mn=d))
      } else {
        res <- rbind(res, data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                                     year=y,
                                     area=NA,
                                     dist.mn=NA))
      }
    }
    res
  }
  stopCluster(cl)
  
  return(dat_comb)
}

## Compute skewness OFD management catchment
skewnessSOD <- function(dat, 
                        cores=NULL,
                        catch_max=NULL,
                        catch="management") {
  
  ## This function will compute temporal mean and standard deviation of the OFD
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {
    
    require(tidyverse) # relational database management
    require(Matrix) # sparse matrix ss tables
    require(e1071) # compute skewness 
    require(vegan) # test for bimodality
    source("scripts/fitMCPD.R")
    
    tokeshiTest <- function(x, h=0.1) {
      require(VeryLargeIntegers)
      S <- length(x)
      x <- x/max(x)
      x <- as.numeric(cut(x, breaks=seq(0,1,by=h)))
      x <- table(x)
      
      P <- c()
      for (i in c(1,length(x))) {
        p <- 0
        for (j in x[i]:(S-1)) {
          logSf <- sum(log(1:S))
          logjf <- sum(log(1:j))
          logSmjf <- sum(log(1:(S-j)))
          bc <- exp(logSf-(logjf + logSmjf))
          p <- p + (bc * (h^j) * ((1-h)^(S-j)))
        }
        P[length(P)+1] <- sum(p)
      }
      return(c(P_l=P[1], P_r=P[2]))
    }
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    
    res <- c()
    for (y in sort(unique(dat_mc$year))) {
      dat_mc_y <- subset(dat_mc, year == y)
      dat_mc_y$site <- as.numeric(factor(dat_mc_y$site))
      dat_mc_y$taxon <- as.numeric(factor(dat_mc_y$taxon))
      no_site <- max(dat_mc_y$site)
      no_taxa <- max(dat_mc_y$taxon)
      site.occ <- Matrix(0, nrow=no_taxa, ncol=no_site, sparse = T)
      for (j in 1:no_taxa) {
        site.occ[j,subset(dat_mc_y, taxon==j)$site] <- subset(dat_mc_y, taxon==j)$count
      }
      
      SOD <- rowSums(1*site.occ>0)
      SOD.prop <- rowSums(1*site.occ>0)/no_site
      mn <- mean(SOD)
      vr <- var(SOD)
      sk <- skewness(SOD)
      mn.prop <- mean(SOD.prop)
      vr.prop <- var(SOD.prop)
      sk.prop <- skewness(SOD.prop)
      mod <- tryCatch(MOStest(as.numeric(names(table(SOD))), table(SOD)),
                      error = function(c) NA)
      if (is.na(mod)) {
        bimod <- NA
      } else {
        bimod <- 1*mod$isHump
      }
      tokeshiP <- tokeshiTest(SOD)
      
      res_mcpd <- f(1*site.occ>0)
      
      res <- rbind(res, data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                                   year=y,
                                   no_site=no_site,
                                   no_taxa=no_taxa,
                                   mn=mn,
                                   vr=vr,
                                   sk=sk,
                                   mn.prop=mn.prop,
                                   vr.prop=vr.prop,
                                   sk.prop=sk.prop,
                                   bimod=bimod,
                                   mixing=res_mcpd$estimate,
                                   P_l=tokeshiP[1],
                                   P_r=tokeshiP[2]))
    }
    res
  }
  stopCluster(cl)
  
  return(dat_comb)
}

temporalSOD <- function(dat, 
                        filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10), 
                        cores=NULL,
                        catch_max=NULL,
                        catch="management") {
  
  ## This function will compute temporal mean and standard deviation of the OFD
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {
    
    require(tidyverse) # relational database management
    require(Matrix) # sparse matrix ss tables
    source("scripts/fitMCPD.R")
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    dat_mc$site <- as.numeric(factor(dat_mc$site))
    dat_mc$taxon <- as.numeric(factor(dat_mc$taxon))
    no_site <- max(dat_mc$site)
    no_taxa <- max(dat_mc$taxon)
    
    dat_mc <- dat_mc %>%
      filter(year > filter[1])
    
    SOD <- list()
    for (y in sort(unique(dat_mc$year))) {
      dat_mc_y <- subset(dat_mc, year == y)
      dat_mc_y$site <- as.numeric(factor(dat_mc_y$site))
      dat_mc_y$taxon <- as.numeric(factor(dat_mc_y$taxon))
      no_site <- max(dat_mc_y$site)
      no_taxa <- max(dat_mc_y$taxon)
      if ((no_site < filter[3])||(no_taxa < filter[2])) {
        next
      }
      site.occ <- Matrix(0, nrow=no_taxa, ncol=no_site, sparse = T)
      for (j in 1:no_taxa) {
        site.occ[j,subset(dat_mc_y, taxon==j)$site] <- subset(dat_mc_y, taxon==j)$count
      }
      SOD[[length(SOD)+1]] <- rowSums(1*site.occ>0)
      SOD[[length(SOD)]] <- SOD[[length(SOD)]][SOD[[length(SOD)]]>0]
    }
    
    if (length(SOD) < filter[4]) {
      res = NULL
    } else {
      SOD_mat <- matrix(0, nrow=length(SOD), ncol=max(sapply(SOD, max)))
      for (y in 1:length(SOD)) {
        SOD_y <- table(SOD[[y]])
        SOD_mat[y, as.numeric(names(SOD_y))] <- SOD_y
      }
      
      SOD_mn <- apply(SOD_mat, MAR=2, FUN=mean)
      if (1) {
        # uncorrected sample standard deviation
        SOD_sd <- apply(SOD_mat, MAR=2, FUN=sd)
      } else {
        # corrected sample standard deviation
        SOD_sd <- apply(SOD_mat, MAR=2, FUN=function(x) sqrt(mean((x-mean(x))^2)))  
      }
      SOD_se <- apply(SOD_mat, MAR=2, FUN=function(x) sqrt(var(x)/length(x)))
      SOD_ci5 <- apply(SOD_mat, MAR=2, FUN=function(x) quantile(x, 0.05))
      SOD_ci95 <- apply(SOD_mat, MAR=2, FUN=function(x) quantile(x, 0.95))
      res_mcpd <- f(SOD_mn)
      
      res <- data.frame(catchment=rep(sort(unique(dat[,catch_index]))[i], length(SOD_mn)),
                        site=1:length(SOD_mn),
                        SOD_mn=SOD_mn,
                        SOD_sd=SOD_sd,
                        SOD_se=SOD_se,
                        SOD_ci5=SOD_ci5,
                        SOD_ci95=SOD_ci95,
                        taxa.pred_ss=res_mcpd$nNiches[1] * res_mcpd$SOD$taxa.pred_ss,
                        mixing=rep(res_mcpd$estimate, length(SOD_mn)))
    }
    res
  }
  stopCluster(cl)
  
  return(dat_comb)
}

skew_accumulation <- function(dat,
                              filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10), 
                              no_rand=1000,
                              cores=NULL,
                              catch_max=NULL,
                              catch="management") {
  
  ## This function will compute the first (and second) derivative skewness OFD under random accumulation of the catchment year
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {
    
    require(tidyverse) # relational database management
    require(Matrix) # sparse matrix ss tables
    require(e1071) # compute skewness 
    require(sfsmisc) # first (and second) derivative
    require(stats) # smooth spline
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    
    dat_mc <- dat_mc %>%
      filter(year > filter[1])
    
    res <- c()
    for (y in sort(unique(dat_mc$year))) {
      dat_mc_y <- subset(dat_mc, year == y)
      dat_mc_y$site <- as.numeric(factor(dat_mc_y$site))
      dat_mc_y$taxon <- as.numeric(factor(dat_mc_y$taxon))
      no_site <- max(dat_mc_y$site)
      no_taxa <- max(dat_mc_y$taxon)
      if (no_site < filter[3] || no_taxa < filter[2]) {
        next
      }
      site.occ <- Matrix(0, nrow=no_taxa, ncol=no_site, sparse = T)
      for (j in 1:no_taxa) {
        site.occ[j,subset(dat_mc_y, taxon==j)$site] <- subset(dat_mc_y, taxon==j)$count
      }
      
      res_rand <- c()
      for (r in 1:no_rand) {
        no_spp <- sample(1:nrow(site.occ),1)
        site.occ_rand <- site.occ[sample(1:nrow(site.occ), no_spp, replace=F),,drop=F]
        sod <- rowSums(1*site.occ_rand>0)
        if (length(unique(sod))==1) {
          sk <- 0
        } else {
          sk <- skewness(sod)
        }
        res_rand <- rbind(res_rand, data.frame(no_spp=no_spp, sk=sk))
      }
      
      res_rand <- res_rand[order(res_rand$no_spp),]
      res_rand <- res_rand[complete.cases(res_rand),]
      
      if (length(unique(res_rand$sk)) < 4) { ## at least 4 unique values required for estimating derivatives
        res <- rbind(res, data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                                     year=as.numeric(y),
                                     no_site=no_site,
                                     no_taxa=no_taxa,
                                     sk_s_d1=1,
                                     sk_s_d2=1))
        next
      }
      
      f1 <- D1ss(res_rand$no_spp, res_rand$sk, spar.off=0.0)
      f2 <- D2ss(res_rand$no_spp, res_rand$sk, spar.off=0.0)
      
      res <- rbind(res, data.frame(catchment=sort(unique(dat[,catch_index]))[i], 
                                   year=as.numeric(y), 
                                   no_site=no_site,
                                   no_taxa=no_taxa,
                                   sk_s_d1=mean(f1[(length(f1)-floor(no_rand/10)):length(f1)]),
                                   sk_s_d2=mean(f2$y[(length(f1)-floor(no_rand/10)):length(f1)])))
    }
    res
  }
  stopCluster(cl)
  
  return(dat_comb)
}

temporalMoments <- function(dat_mom, 
                            dat_rand,
                            filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 8, thresh = -2)) {
  
  # This function will analyse temporal fluctations in the moments of the OFD
  # dat_mom: moments of the OFD or each catchment year
  # dat_rand: outcome of skewness accumulation experiment
  
  # filter by curvature of skew accumulation experiment
  filtered_catchment_years <- do.call(paste,subset(dat_rand, log10(abs(dat_rand$sk_s_d1)) > filter[5])[,c(grep("catchment", names(dat_rand)),grep("year", names(dat_rand)))])
  dat_mom <- dat_mom[!do.call(paste,dat_mom[,c(grep("catchment", names(dat_mom)),grep("year", names(dat_mom)))]) %in% filtered_catchment_years,]
  
  # analysis of OFD (abs and prop) moment time series
  dat_ret <- c()
  for (cc in sort(unique(dat_mom$catchment))) {
    dat_sub <- subset(dat_mom, catchment==cc & year>=filter[1] & no_taxa>=filter[2] & no_site>=filter[3])
    if (nrow(dat_sub) < filter[4]) {
      next
    }
    dat_sub <- dat_sub[order(dat_sub$year), ]
    
    # mean
    mn_na_impute <- rep(NA, 1 + max(dat_sub$year) - min(dat_sub$year))
    mn_na_impute[1+dat_sub$year-min(dat_sub$year)] <- dat_sub$mn
    tt <- min(dat_sub$year):max(dat_sub$year)
    mod.lin.mn <- lm(mn_na_impute ~ tt)
    mn_na_impute <- na_interpolation(mn_na_impute)
    mod.mk.mn <- mk.test(mn_na_impute)
    mod.sens.mn <- sens.slope(mn_na_impute)
    int.sens.mn <- median(mn_na_impute)-(mod.sens.mn$estimates)*median(tt)
    mod.adf.mn <- adf.test(mn_na_impute)
    mod.kpss.mn <- kpss.test(mn_na_impute)
    
    # variance
    vr_na_impute <- rep(NA, 1 + max(dat_sub$year) - min(dat_sub$year))
    vr_na_impute[1+dat_sub$year-min(dat_sub$year)] <- dat_sub$vr
    tt <- min(dat_sub$year):max(dat_sub$year)
    mod.lin.vr <- lm(vr_na_impute ~ tt)
    vr_na_impute <- na_interpolation(vr_na_impute)
    mod.mk.vr <- mk.test(vr_na_impute)
    mod.sens.vr <- sens.slope(vr_na_impute)
    int.sens.vr <- median(vr_na_impute)-(mod.sens.vr$estimates)*median(tt)
    mod.adf.vr <- adf.test(vr_na_impute)
    mod.kpss.vr <- kpss.test(vr_na_impute)
    
    # skewness
    sk_na_impute <- rep(NA, 1 + max(dat_sub$year) - min(dat_sub$year))
    sk_na_impute[1+dat_sub$year-min(dat_sub$year)] <- dat_sub$sk
    tt <- min(dat_sub$year):max(dat_sub$year)
    mod.lin.sk <- lm(sk_na_impute ~ tt)
    sk_na_impute <- na_interpolation(sk_na_impute)
    mod.mk.sk <- mk.test(sk_na_impute)
    mod.sens.sk <- sens.slope(sk_na_impute)
    int.sens.sk <- median(sk_na_impute)-(mod.sens.sk$estimates)*median(tt)
    mod.adf.sk <- adf.test(sk_na_impute)
    mod.kpss.sk <- kpss.test(sk_na_impute)

    # mean (proportional)
    mn.prop_na_impute <- rep(NA, 1 + max(dat_sub$year) - min(dat_sub$year))
    mn.prop_na_impute[1+dat_sub$year-min(dat_sub$year)] <- dat_sub$mn.prop
    tt <- min(dat_sub$year):max(dat_sub$year)
    mod.lin.mn.prop <- lm(mn.prop_na_impute ~ tt)
    mn.prop_na_impute <- na_interpolation(mn.prop_na_impute)
    mod.mk.mn.prop <- mk.test(mn.prop_na_impute)
    mod.sens.mn.prop <- sens.slope(mn.prop_na_impute)
    int.sens.mn.prop <- median(mn.prop_na_impute)-(mod.sens.mn.prop$estimates)*median(tt)
    mod.adf.mn.prop <- adf.test(mn.prop_na_impute)
    mod.kpss.mn.prop <- kpss.test(mn.prop_na_impute)

    # variance (proportional)
    vr.prop_na_impute <- rep(NA, 1 + max(dat_sub$year) - min(dat_sub$year))
    vr.prop_na_impute[1+dat_sub$year-min(dat_sub$year)] <- dat_sub$vr.prop
    tt <- min(dat_sub$year):max(dat_sub$year)
    mod.lin.vr.prop <- lm(vr.prop_na_impute ~ tt)
    vr.prop_na_impute <- na_interpolation(vr.prop_na_impute)
    mod.mk.vr.prop <- mk.test(vr.prop_na_impute)
    mod.sens.vr.prop <- sens.slope(vr.prop_na_impute)
    int.sens.vr.prop <- median(vr.prop_na_impute)-(mod.sens.vr.prop$estimates)*median(tt)
    mod.adf.vr.prop <- adf.test(vr.prop_na_impute)
    mod.kpss.vr.prop <- kpss.test(vr.prop_na_impute)

    # skewness (proportional)
    sk.prop_na_impute <- rep(NA, 1 + max(dat_sub$year) - min(dat_sub$year))
    sk.prop_na_impute[1+dat_sub$year-min(dat_sub$year)] <- dat_sub$sk
    tt <- min(dat_sub$year):max(dat_sub$year)
    mod.lin.sk.prop <- lm(sk.prop_na_impute ~ tt)
    sk.prop_na_impute <- na_interpolation(sk.prop_na_impute)
    mod.mk.sk.prop <- mk.test(sk.prop_na_impute)
    mod.sens.sk.prop <- sens.slope(sk.prop_na_impute)
    int.sens.sk.prop <- median(sk.prop_na_impute)-(mod.sens.sk.prop$estimates)*median(tt)
    mod.adf.sk.prop <- adf.test(sk.prop_na_impute)
    mod.kpss.sk.prop <- kpss.test(sk.prop_na_impute)
    
    dat_ret <- rbind(dat_ret, data.frame(catchment=cc,
                                         no_years=nrow(dat_sub),
                                         S.mn=mean(dat_sub$no_taxa),
                                         N.mn=mean(dat_sub$no_site),
                                         # mean
                                         mn.t.lin=coefficients(mod.lin.mn)[2],
                                         mn.t.lin.r2=summary(mod.lin.mn)$adj.r.squared,
                                         p.lin.mn=summary(mod.lin.mn)$coefficients[2,4],
                                         p.mkt.mn=mod.mk.mn$p.value,
                                         mn.t.sens=mod.sens.mn$estimates,
                                         mn.t.sens.int=int.sens.mn,
                                         p.adf.mn=mod.adf.mn$p.value,
                                         mn.adf=mod.adf.mn$statistic,
                                         p.kpss.mn=mod.kpss.mn$p.value,
                                         # variance
                                         vr.t.lin=coefficients(mod.lin.vr)[2],
                                         vr.t.lin.r2=summary(mod.lin.vr)$adj.r.squared,
                                         p.lin.vr=summary(mod.lin.vr)$coefficients[2,4],
                                         p.mkt.vr=mod.mk.vr$p.value,
                                         vr.t.sens=mod.sens.vr$estimates,
                                         vr.t.sens.int=int.sens.vr,
                                         p.adf.vr=mod.adf.vr$p.value,
                                         vr.adf=mod.adf.vr$statistic,
                                         p.kpss.vr=mod.kpss.vr$p.value,
                                         # skewness
                                         sk.t.lin=coefficients(mod.lin.sk)[2],
                                         sk.t.lin.r2=summary(mod.lin.sk)$adj.r.squared,
                                         p.lin.sk=summary(mod.lin.sk)$coefficients[2,4],
                                         p.mkt.sk=mod.mk.sk$p.value,
                                         sk.t.sens=mod.sens.sk$estimates,
                                         sk.t.sens.int=int.sens.sk,
                                         p.adf.sk=mod.adf.sk$p.value,
                                         sk.adf=mod.adf.sk$statistic,
                                         p.kpss.sk=mod.kpss.sk$p.value,
                                         # mean (proportional)
                                         mn.prop.t.lin=coefficients(mod.lin.mn.prop)[2],
                                         mn.prop.t.lin.r2=summary(mod.lin.mn.prop)$adj.r.squared,
                                         p.lin.mn.prop=summary(mod.lin.mn.prop)$coefficients[2,4],
                                         p.mkt.mn.prop=mod.mk.mn.prop$p.value,
                                         mn.prop.t.sens=mod.sens.mn.prop$estimates,
                                         mn.prop.t.sens.int=int.sens.mn.prop,
                                         p.adf.mn.prop=mod.adf.mn.prop$p.value,
                                         mn.prop.adf=mod.adf.mn.prop$statistic,
                                         p.kpss.mn.prop=mod.kpss.mn.prop$p.value,
                                         # variance (proportional)
                                         vr.prop.t.lin=coefficients(mod.lin.vr.prop)[2],
                                         vr.prop.t.lin.r2=summary(mod.lin.vr.prop)$adj.r.squared,
                                         p.lin.vr.prop=summary(mod.lin.vr.prop)$coefficients[2,4],
                                         p.mkt.vr.prop=mod.mk.vr.prop$p.value,
                                         vr.prop.t.sens=mod.sens.vr.prop$estimates,
                                         vr.prop.t.sens.int=int.sens.vr.prop,
                                         p.adf.vr.prop=mod.adf.vr.prop$p.value,
                                         vr.prop.adf=mod.adf.vr.prop$statistic,
                                         p.kpss.vr.prop=mod.kpss.vr.prop$p.value,
                                         # skewness (proportional)
                                         sk.prop.t.lin=coefficients(mod.lin.sk.prop)[2],
                                         sk.prop.t.lin.r2=summary(mod.lin.sk.prop)$adj.r.squared,
                                         p.lin.sk.prop=summary(mod.lin.sk.prop)$coefficients[2,4],
                                         p.mkt.sk.prop=mod.mk.sk.prop$p.value,
                                         sk.prop.t.sens=mod.sens.sk.prop$estimates,
                                         sk.prop.t.sens.int=int.sens.sk.prop,
                                         p.adf.sk.prop=mod.adf.sk.prop$p.value,
                                         sk.prop.adf=mod.adf.sk.prop$statistic,
                                         p.kpss.sk.prop=mod.kpss.sk.prop$p.value))
  }
  return(dat_ret)
}

turnover <- function(dat,
                     filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10), 
                     cores=NULL,
                     catch_max=NULL,
                     catch="management") {
  
  ## This function will compute the decay in similarity with temporal distance averaged over all start years
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {

    require(tidyverse) # relational database management
    require(Matrix) # sparse matrix ss tables
    require(vegan) # compute beta diversity
    require(trend)
    require(imputeTS)
    require(tseries)
  
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    
    dat_mc <- dat_mc %>%
      filter(year > filter[1])
    
    dat_mc$site <- as.numeric(factor(dat_mc$site))
    dat_mc$taxon <- as.numeric(factor(dat_mc$taxon))
    no_taxa <- max(dat_mc$taxon)
    no_year <- length(unique(dat_mc$year))
    no_site <- max(dat_mc$site)
    
    if ((no_year < filter[4]) || (no_taxa < filter[2]) || (no_site < filter[3])) {
      res = NULL
    } else {
      
      site.occ <- Matrix(0, nrow=no_taxa, ncol=no_year, sparse = T)
      
      for (y in 1:no_year) {
        site.occ[subset(dat_mc, year==sort(unique(dat_mc$year))[y])$taxon,y] <- 1
      }
      J <- as.matrix(vegdist(t(site.occ), method="jaccard"))
      J.t <- data.frame(J = 1-as.vector(J),
                        y = as.vector(as.matrix(dist(as.numeric(sort(unique(dat_mc$year)))))),
                        y0 = rep(1:length(unique(dat_mc$year)), each=nrow(J)))
      
      J.t.mn <- J.t %>% 
        group_by(y) %>%
        summarise(J=mean(J)) %>%
        filter(y > min(y))
      
      # Jaccard similarity
      J_na_impute <- rep(NA, 1 + max(J.t.mn) - min(J.t.mn$y))
      J_na_impute[1+J.t.mn$y-min(J.t.mn$y)] <- J.t.mn$J
      tt <- min(J.t.mn$y):max(J.t.mn$y)
      mod.lin.J <- lm(J_na_impute ~ tt)
      J_na_impute <- na_interpolation(J_na_impute)
      mod.mk.J <- mk.test(J_na_impute[1:floor(length(J_na_impute)/2)])
      mod.sens.J <- sens.slope(J_na_impute[1:floor(length(J_na_impute)/2)])
      mod.adf.J <- adf.test(J_na_impute[1:floor(length(J_na_impute)/2)])
      mod.kpss.J <- kpss.test(J_na_impute[1:floor(length(J_na_impute)/2)])
      
      res <- data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                        no_year = no_year,
                        no_taxa = no_taxa,
                        J.t.lin=coefficients(mod.lin.J)[2],
                        J.t.lin.r2=summary(mod.lin.J)$adj.r.squared,
                        p.lin.J=summary(mod.lin.J)$coefficients[2,4],
                        p.mkt.J=mod.mk.J$p.value,
                        J.t.sens=mod.sens.J$estimates,
                        p.adf.J=mod.adf.J$p.value,
                        J.adf=mod.adf.J$statistic,
                        p.kpss.J=mod.kpss.J$p.value)
    }
    res
  }
  stopCluster(cl)
  return(dat_comb)
}

turnoverLocal <- function(dat,
                          filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10), 
                          cores=NULL,
                          catch_max=NULL,
                          catch="management") {
  
  ## This function will compute the decay in similarity with temporal distance averaged over all start years
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {
    require(tidyverse) # relational database management
    require(vegan) # compute beta diversity
    require(Matrix) # sparse matrix ss tables
    require(trend)
    require(imputeTS)
    require(tseries)
    require(betapart)
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    
    dat_mc <- dat_mc %>%
      filter(year > filter[1])
    
    dat_mc$site <- as.numeric(factor(dat_mc$site))
    dat_mc$taxon <- as.numeric(factor(dat_mc$taxon))
    no_taxa <- max(dat_mc$taxon)
    no_year <- length(unique(dat_mc$year))
    no_site <- max(dat_mc$site)
    
    if ((no_year < filter[4]) || (no_taxa < filter[2]) || (no_site < filter[3])) {
      res = NULL
    } else {
      
      res <- c()
      for (x in 1:no_site) {
        dat_mc_x <- subset(dat_mc, site==x)
        dat_mc_x$taxon <- as.numeric(as.factor(dat_mc_x$taxon))
        # dat_mc_x$year <- as.numeric(as.factor(dat_mc_x$year))
        no_taxa <- max(dat_mc_x$taxon)
        no_year <- length(unique(dat_mc_x$year))
        if (no_year < filter[4]) {
          next
        }
          
        site.occ <- Matrix(0, nrow=no_taxa, ncol=no_year, sparse = T)
        for (y in 1:no_year) {
          site.occ[subset(dat_mc_x, year==sort(unique(dat_mc_x$year))[y])$taxon,y] <- 1
        }
        
        J <- as.matrix(vegdist(t(site.occ), method="jaccard"))
        J.t <- data.frame(J = 1-as.vector(J),
                          y = as.vector(as.matrix(dist(as.numeric(sort(unique(dat_mc_x$year)))))),
                          y0 = rep(1:length(unique(dat_mc_x$year)), each=nrow(J)))
        
        J.t.mn <- J.t %>% 
          group_by(y) %>%
          summarise(J=mean(J)) %>%
          filter(y > min(y))
        J.t.mn$x <- x
        J.t.mn$no_years <- 1+nrow(J.t.mn)
        
        B.t <- betapart.core(t(site.occ))
        B.t.res <-beta.pair(B.t, index.family="jac")
        
        B.t_turn <- mean(as.matrix(B.t.res$beta.jtu))
        B.t_nest <- mean(as.matrix(B.t.res$beta.jne))
        
        res <- rbind(res, cbind(J.t.mn, turn.mn=mean(as.matrix(B.t.res$beta.jtu)), 
                                nest.mn=mean(as.matrix(B.t.res$beta.jne)), tot.mn=mean(as.matrix(B.t.res$beta.jac))))
      }
      
      res$catchment <- sort(unique(dat[,catch_index]))[i]
    }
    if (is.data.frame(res)) {
      res  
    }
  }
  stopCluster(cl)
  return(dat_comb)
}

richnessLocal <- function(dat,
                          filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10, i0=100), 
                          cores=NULL,
                          catch_max=NULL,
                          catch="management") {
  
  ## This function will compute the temporal fluctuations in local richness for available sites
  
  require(tidyverse) # relational database management
  require(vegan) # compute beta diversity
  require(Matrix) # sparse matrix ss tables
  require(trend)
  require(imputeTS)
  require(tseries)

  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {
    require(tidyverse) # relational database management
    require(vegan) # compute beta diversity
    require(Matrix) # sparse matrix ss tables
    require(trend)
    require(imputeTS)
    require(tseries)
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    
    dat_mc <- dat_mc %>%
      filter(year > filter[1])
    
    dat_mc$site <- as.numeric(factor(dat_mc$site))
    dat_mc$taxon <- as.numeric(factor(dat_mc$taxon))
    no_taxa <- max(dat_mc$taxon)
    no_year <- length(unique(dat_mc$year))
    no_site <- max(dat_mc$site)
        
    if ((no_year < filter[4]) || (no_taxa < filter[2]) || (no_site < filter[3])) {
      res = NULL
    } else {
      
      res <- c()
      for (x in 1:no_site) {
        # print(x)
        dat_mc_x <- subset(dat_mc, site==x)
        dat_mc_x$taxon <- as.numeric(as.factor(dat_mc_x$taxon))
        no_taxa <- max(dat_mc_x$taxon)
        
        richness <- dat_mc_x %>%
          group_by(year) %>%
          summarise(S=length(unique(taxon)),
                    count.tot=sum(count))
        
        dat_mc_x <- dat_mc_x %>%
          filter(year %in% subset(richness, count.tot >= filter[5])$year)
        
        no_year <- length(unique(dat_mc_x$year))
        
        if (no_year < filter[4]) {
          next
        }
        
        if (nrow(dat_mc_x)==0) {
          next
        }
        
        richness <- richness %>% 
          filter(count.tot >= filter[5])
        
        dat_mc_x <- dat_mc_x %>% 
          left_join(richness[,-2], by="year")
        
        dat_mc_x$count.rare <- 0
        
        rep_rare <- 10
        
        S_rare <- matrix(0, nrow(richness), rep_rare)
        
        for (r in 1:rep_rare) {
          for (obs in 1:nrow(dat_mc_x)) {
            dat_mc_x$count.rare[obs] <- rbinom(1, size=dat_mc_x$count[obs], prob=min(richness$count.tot)/dat_mc_x$count.tot[obs])
          }
          S_rare[,r] <- (dat_mc_x %>% 
            group_by(year) %>%
            summarise(S=length(which(count.rare>0))))$S
        }
        
        richness$S_rare.mn <- rowMeans(S_rare)
        richness$S_rare.sd <- apply(S_rare, 1, sd)
        
        res <- rbind(res, cbind(site=x, richness))
      }
      res$catchment <- sort(unique(dat[,catch_index]))[i]
      if (0) {
        ggplot(res, aes(x=year, y=S_rare.mn, col=factor(site))) +
          geom_point() +
          theme(legend.position="none") +
          geom_smooth(method='lm', se=F)
      }
    }
    
    if (is.data.frame(res)) {
      res  
    }
  }
  stopCluster(cl)
  return(dat_comb)
}

turnoverSOD <- function(dat,
                        filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10), 
                        cores=NULL,
                        catch_max=NULL,
                        catch="management") {
  
  ## This function will compute the decay in similarity of OFD with temporal distance averaged over all start years
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {
    
    require(tidyverse) # relational database management
    require(Matrix) # sparse matrix ss tables
    require(e1071) # compute skewness 
    require(vegan) # compute beta diversity
    require(trend)
    require(imputeTS)
    require(tseries)
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    
    dat_mc <- dat_mc %>%
      filter(year > filter[1])
    
    dat_mc$site <- as.numeric(factor(dat_mc$site))
    dat_mc$taxon <- as.numeric(factor(dat_mc$taxon))
    no_taxa <- max(dat_mc$taxon)
    no_year <- length(unique(dat_mc$year))
    no_site <- max(dat_mc$site)
    
    if ((no_year < filter[4]) || (no_taxa < filter[2]) || (no_site < filter[3])) {
      res = NULL
    } else {
      
      site.occ <- Matrix(0, nrow=no_taxa, ncol=no_site, sparse = T)
      SOD <- Matrix(0, nrow=no_taxa, ncol=no_year, sparse = T)
      
      for (y in 1:no_year) {
        dat_y <- subset(dat_mc, year==sort(unique(dat_mc$year))[y])
        site.occ <- Matrix(0, nrow=no_taxa, ncol=no_site, sparse = T)
        for (j in unique(dat_y$taxon)) {
          site.occ[j,subset(dat_y, taxon==j)$site] <- 1  
        }
        SOD[,y] <- rowSums(site.occ)/length(unique(dat_y$site))
      }
      B <- as.matrix(vegdist(t(SOD), method="bray"))
      B.t <- data.frame(B = 1-as.vector(B),
                        y = as.vector(as.matrix(dist(as.numeric(sort(unique(dat_mc$year)))))),
                        y0 = rep(1:length(unique(dat_mc$year)), each=nrow(B)))
      
      B.t.mn <- B.t %>% 
        group_by(y) %>%
        summarise(B=mean(B)) %>%
        filter(y > min(y))
      
      # Bray-Curtis similarity
      B_na_impute <- rep(NA, 1 + max(B.t.mn) - min(B.t.mn$y))
      B_na_impute[1+B.t.mn$y-min(B.t.mn$y)] <- B.t.mn$B
      tt <- min(B.t.mn$y):max(B.t.mn$y)
      mod.lin.B <- lm(B_na_impute ~ tt)
      B_na_impute <- na_interpolation(B_na_impute)
      mod.mk.B <- mk.test(B_na_impute[1:floor(length(B_na_impute)/2)])
      mod.sens.B <- sens.slope(B_na_impute[1:floor(length(B_na_impute)/2)])
      mod.adf.B <- adf.test(B_na_impute[1:floor(length(B_na_impute)/2)])
      mod.kpss.B <- kpss.test(B_na_impute[1:floor(length(B_na_impute)/2)])
      
      res <- data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                        no_year = no_year,
                        no_taxa = no_taxa,
                        B.t.lin=coefficients(mod.lin.B)[2],
                        B.t.lin.r2=summary(mod.lin.B)$adj.r.squared,
                        p.lin.B=summary(mod.lin.B)$coefficients[2,4],
                        p.mkt.B=mod.mk.B$p.value,
                        B.t.sens=mod.sens.B$estimates,
                        p.adf.B=mod.adf.B$p.value,
                        B.adf=mod.adf.B$statistic,
                        p.kpss.B=mod.kpss.B$p.value)
    }
    res
  }
  stopCluster(cl)  
  return(dat_comb)
}

turnoverSOD_raw <- function(dat,
                            filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10), 
                            cores=NULL,
                            catch_max=NULL,
                            catch="management",
                            method="bray") {
  
  ## This function will compute the decay in similarity of OFD with temporal distance averaged over all start years
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {
    
    require(tidyverse) # relational database management
    require(Matrix) # sparse matrix ss tables
    require(vegan) # compute beta diversity
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    
    dat_mc <- dat_mc %>%
      filter(year > filter[1])
    
    dat_mc$site <- as.numeric(factor(dat_mc$site))
    dat_mc$taxon <- as.numeric(factor(dat_mc$taxon))
    no_taxa <- max(dat_mc$taxon)
    no_year <- length(unique(dat_mc$year))
    no_site <- max(dat_mc$site)
    
    if ((no_year < filter[4]) || (no_taxa < filter[2]) || (no_site < filter[3])) {
      res = NULL
    } else {
      
      site.occ <- Matrix(0, nrow=no_taxa, ncol=no_site, sparse = T)
      SOD <- Matrix(0, nrow=no_taxa, ncol=no_year, sparse = T)
      
      for (y in 1:no_year) {
        dat_y <- subset(dat_mc, year==sort(unique(dat_mc$year))[y])
        site.occ <- Matrix(0, nrow=no_taxa, ncol=no_site, sparse = T)
        for (j in unique(dat_y$taxon)) {
          site.occ[j,subset(dat_y, taxon==j)$site] <- 1  
        }
        SOD[,y] <- rowSums(site.occ)/length(unique(dat_y$site))
      }
      if (method=="bray") {
        B <- as.matrix(vegdist(t(SOD), method="bray"))
        B.t <- data.frame(B = 1-as.vector(B),
                          y = as.vector(as.matrix(dist(as.numeric(sort(unique(dat_mc$year)))))),
                          y0 = rep(1:length(unique(dat_mc$year)), each=nrow(B)))
        
        B.t.mn <- B.t %>% 
          group_by(y) %>%
          summarise(B=mean(B))
        
        res <- data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                          dt=B.t.mn$y,
                          BC.mn=B.t.mn$B)
      } else if (method=="jaccard") {
        J <- as.matrix(vegdist(t(1*SOD>0), method="jaccard"))
        J.t <- data.frame(J = 1-as.vector(J),
                          y = as.vector(as.matrix(dist(as.numeric(sort(unique(dat_mc$year)))))),
                          y0 = rep(1:length(unique(dat_mc$year)), each=nrow(J)))
        
        J.t.mn <- J.t %>% 
          group_by(y) %>%
          summarise(J=mean(J))
        
        res <- data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                          dt=J.t.mn$y,
                          J.mn=J.t.mn$J)
      }
    }
    res
  }
  stopCluster(cl) 
  return(dat_comb)
}

timescales <- function(dat,
                       dat_rand,
                       filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 8, thresh=-2), 
                       cores=NULL,
                       catch_max=NULL,
                       catch="management") {
  
  ## This function will compute the decay in similarity of OFD with temporal distance averaged over all start years
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  # filter by curvature of skew accumulation experiment
  filtered_catchment_years <- do.call(paste,subset(dat_rand, log10(abs(dat_rand$sk_s_d1)) > filter[5])[,c(grep("catchment", names(dat_rand)),grep("year", names(dat_rand)))])
  dat <- dat[!do.call(paste,dat[,c(catch_index, grep("year", names(dat)))]) %in% filtered_catchment_years,]
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms) %dopar% {
    
    require(tidyverse) # relational database management
    source("scripts/fitTimeScales.R")
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    
    dat_mc <- dat_mc %>%
      filter(year > filter[1])
    
    dat_mc_y <- dat_mc %>%
      group_by(year) %>%
      summarise(no_taxa = length(unique(taxon)), no_site = length(unique(site))) %>%
      filter(no_taxa >= filter[2] & no_site >= filter[3])
    
    if ((nrow(dat_mc_y) < filter[4])) {
      res = "Filtered"
    } else {
    
      dat_mc$site <- as.numeric(factor(dat_mc$site))
      dat_mc$taxon <- as.numeric(factor(dat_mc$taxon))
      no_taxa <- max(dat_mc$taxon)
      no_year <- length(unique(dat_mc$year))
      no_site <- max(dat_mc$site)
        
      site.occ_array <- array(0, dim=c(no_taxa, no_site, no_year))

      for (y in 1:no_year) {
        dat_y <- subset(dat_mc, year==sort(unique(dat_mc$year))[y])
        for (j in unique(dat_y$taxon)) {
          site.occ_array[j,subset(dat_y, taxon==j)$site,y] <- 1
        }
      }

      dimnames(site.occ_array) <- list(i=NULL, x=NULL, year=as.character(sort(unique(dat_mc$year))))

      res <- fitTimeScales(site.occ_array, sort(unique(dat[,catch_index]))[i], filter[4])
    }
    res
  }
  stopCluster(cl)
  return(dat_comb)
}

l_n_occupancy <- function(dat,
                          cores=NULL,
                          catch_max=NULL,
                          catch="management") {
  
  ## This function will compute the catchment average and national occupancy for each species
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {
    
    require(tidyverse) # relational database management
    require(Matrix) # sparse matrix ss tables
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,TAXON_RANK,site,year) %>%
      summarise(count=sum(count))
    
    res <- c()
    for (y in sort(unique(dat_mc$year))) {
      dat_mc_y <- subset(dat_mc, year == y)
      dat_mc_y$site <- as.numeric(factor(dat_mc_y$site))
      dat_mc_y$taxon2 <- as.numeric(factor(dat_mc_y$taxon))
      no_site <- max(dat_mc_y$site)
      no_taxa <- max(dat_mc_y$taxon2)
      site.occ <- Matrix(0, nrow=no_taxa, ncol=no_site, sparse = T)
      row.names(site.occ) <- as.character(unique(paste0(dat_mc_y$taxon, "_", dat_mc_y$TAXON_RANK)))
      for (j in 1:no_taxa) {
        site.occ[j,subset(dat_mc_y, taxon2==j)$site] <- subset(dat_mc_y, taxon2==j)$count
      }
      
      SOD <- rowSums(1*site.occ>0)
      SOD.prop <- rowSums(1*site.occ>0)/no_site
      
      res <- rbind(res, data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                                   year=y,
                                   no_site=no_site,
                                   no_taxa=no_taxa,
                                   taxon=rownames(site.occ),
                                   SO_i=SOD,
                                   SO_i.prop=SOD.prop))
    }
    res
  }
  stopCluster(cl)
  
  catch.avg <- dat_comb %>%
    group_by(year, taxon) %>%
    summarise(CA_SO_i=mean(SO_i),
              CA_SO_i.prop=mean(SO_i.prop))
  
  national <- dat_comb %>%
    group_by(year, taxon) %>%
    summarise(N_SO_i=sum(SO_i),
              N_SO_i.prop=sum(SO_i)/sum(no_site))
  
  dat_merge <- catch.avg %>%
    left_join(national, by=c("taxon", "year"))
  
  return(dat_merge)
}

invasionRate <- function(dat,
                         filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10), 
                         cores=NULL,
                         catch_max=NULL,
                         catch="management") {
  
  ## This function will estimate metacommunity invasion rate for parameterising LSPOM
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {
    
    require(tidyverse) # relational database management
    require(Matrix) # sparse matrix ss tables
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    
    dat_mc <- dat_mc %>%
      filter(year > filter[1])
    
    dat_mc$site <- as.numeric(factor(dat_mc$site))
    dat_mc$taxon <- as.numeric(factor(dat_mc$taxon))
    no_taxa <- max(dat_mc$taxon)
    no_year <- length(unique(dat_mc$year))
    no_site <- max(dat_mc$site)
    
    if ((no_year < filter[4]) || (no_taxa < filter[2]) || (no_site < filter[3])) {
      res = NULL
    } else {
      
      site.occ <- Matrix(0, nrow=no_taxa, ncol=no_year, sparse = T)
      
      for (y in 1:no_year) {
        site.occ[subset(dat_mc, year==sort(unique(dat_mc$year))[y])$taxon,y] <- 1
      }
      
      site.occ2 <- Matrix(0, nrow=no_taxa, ncol=no_year, sparse = T)
      
      for (j in 1:nrow(site.occ2)) {
        site.occ2[j, which(site.occ[j,]>0)[1]:ncol(site.occ2)] <- 1
      }
      site.occ2 <- site.occ2[order(apply(as.matrix(site.occ2), MAR=1, function(x) which(x>0)[1])),]
      image(as.matrix(site.occ2))
      plot(colSums(site.occ2) ~ sort(unique(dat_mc$year)))
      mod <- lm(colSums(site.occ2) ~ sort(unique(dat_mc$year)))
      abline(mod)
      
      res <- data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                        no_year = no_year,
                        no_taxa = no_taxa,
                        first_detect_py = coefficients(summary(mod))[2,1])
    }
    res
  }
  stopCluster(cl)
  return(dat_comb)
}

chaoEstimatorMC <- function(dat,
                            dat_rand,
                            filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10), 
                            cores=NULL,
                            catch_max=NULL,
                            catch="management") {
  
  # This function will compute the proportional observation rate [S_obs / Chao estimator (\hat{S})] for each year in a given catchment
  
  chaoObs <- function(presAbsVec) {
    # This function will compute the proportional observation rate using the Chao estimator
    D <- length(which(presAbsVec>0))
    f1 <- length(which(presAbsVec==1))
    f2 <- length(which(presAbsVec==2))
    
    if (f1==0) {
      return(NA)
    } else if (f2==0) {
      S <- D + ((f1*(f1-1))/(2*(f2+1)))
      return(D/S)
    } else {
      S <- D + ((f1^2)/(2*f2))
      return(D/S)
    }
  }
    
  ## Compute temporal mean and standard deviation of the OFD
  if (catch == "management") {
    catch_index <- which(colnames(dat) == "Management.catchment")  
  } else if (catch == "operational") {
    catch_index <- which(colnames(dat) == "Operational.catchment")  
  } else if (catch == "rbd") {
    catch_index <- which(colnames(dat) == "River.basin.district")
  }
  
  if (!is.null(catch_max)) {
    no_metacomms <- catch_max  
  } else {
    no_metacomms <- length(unique(dat[,catch_index]))
  }
  
  if (is.null(cores)) {
    cores=detectCores()  
  }
  
  if (!is.null(catch_max)) {
    cores=min(c(cores,catch_max))
  }
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # filter by curvature of skew accumulation experiment
  filtered_catchment_years <- do.call(paste,subset(dat_rand, log10(abs(dat_rand$sk_s_d1)) > filter[5])[,c(grep("catchment", names(dat_rand)),grep("year", names(dat_rand)))])
  dat <- dat[!do.call(paste,dat[,c(catch_index, grep("year", names(dat)))]) %in% filtered_catchment_years,]
  
  dat_comb <- foreach (i=1:no_metacomms, .combine=rbind) %dopar% {
    
    require(tidyverse) # relational database management
    require(Matrix) # sparse matrix ss tables
    
    dat_mc <- subset(dat, dat[,catch_index]==sort(unique(dat[,catch_index]))[i]) %>%
      group_by(taxon,site,year) %>%
      summarise(count=sum(count))
    dat_mc$site <- as.numeric(factor(dat_mc$site))
    dat_mc$taxon <- as.numeric(factor(dat_mc$taxon))
    no_site <- max(dat_mc$site)
    no_taxa <- max(dat_mc$taxon)
    
    dat_mc <- dat_mc %>%
      filter(year > filter[1])
    
    chaoObsMC <- c()
    for (y in sort(unique(dat_mc$year))) {
      dat_mc_y <- subset(dat_mc, year == y)
      dat_mc_y$site <- as.numeric(factor(dat_mc_y$site))
      dat_mc_y$taxon <- as.numeric(factor(dat_mc_y$taxon))
      no_site <- max(dat_mc_y$site)
      no_taxa <- max(dat_mc_y$taxon)
      if ((no_site < filter[3])||(no_taxa < filter[2])) {
        next
      }
      site.occ <- Matrix(0, nrow=no_taxa, ncol=no_site, sparse = T)
      for (j in 1:no_taxa) {
        site.occ[j,subset(dat_mc_y, taxon==j)$site] <- subset(dat_mc_y, taxon==j)$count
      }
      presAbsMc <- rowSums(site.occ)
      chaoObsMC[length(chaoObsMC)+1] <- chaoObs(presAbsMc)
    }
    
    chaoObsMC <- chaoObsMC[!is.na(chaoObsMC)]
    
    if (length(chaoObsMC) == 0) {
      res = NULL
    } else {
      
      res <- data.frame(catchment=sort(unique(dat[,catch_index]))[i],
                        chaoObs.med=median(chaoObsMC))
    }
    res
  }
  stopCluster(cl)
  
  return(dat_comb)
}
