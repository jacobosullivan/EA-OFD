## This script will compute the skew of the OFD at the management catchment scale and run various spatio-temporal analyses on result
## Bulk data here: https://environment.data.gov.uk/ecology-fish/downloads/
## Analysis of invert data

require(tidyverse)
setwd("/data/home/btx718/Paleolimnology/datasets/EA/") # set working directory to EA-OFD repository
source("scripts/analysis_functions.R")

################################################################################
################################### Diatoms ####################################
################################################################################

if (0) {
  
  ## Assemble dataset (first unzip files)
  metr_dit <- as_tibble(read.csv("Diatoms/DIAT_OPEN_DATA_METRICS.csv"))
  site_dit <- as_tibble(read.csv("Diatoms/DIAT_OPEN_DATA_SITE.csv"))
  taxa_dit <- as_tibble(read.csv("Diatoms/DIAT_OPEN_DATA_TAXA.csv"))
  taxaRef <- as_tibble(read.csv("TaxonInfo/OPEN_DATA_TAXON_INFO.csv"))
  
  setwd("Status/")
  
  # Load classifications for all River Basin Districts, England
  class_rbd_files <- grep("classifications_RBD*", list.files(recursive = TRUE), value=T)
  wfd_id <- c()
  for (i in 1:length(class_rbd_files)) {
    wfd_id <- rbind(wfd_id, read.csv(class_rbd_files[i]))
  }
  wfd_id <- as_tibble(wfd_id)
  wfd_id <- wfd_id %>% 
    select(River.basin.district, 
           Management.catchment, 
           Operational.catchment, 
           Waterbody.ID, 
           water.body,
           Easting,
           Northing)
  wfd_id <- wfd_id[!duplicated(wfd_id$Waterbody.ID),]
  
  ## Skewness POD macrophytes (species+genus)
  dat_dit <- taxa_dit %>% # select relevant response data
    select(ANALYSIS_ID, TAXON_LIST_ITEM_KEY, UNITS_FOUND_OR_SEQUENCE_READS)
  dat_dit <- dat_dit %>% # add corresponding sample ID
    left_join(subset(metr_dit, select=c(ANALYSIS_ID, SAMPLE_ID, SAMPLE_DATE), by="ANALYSIS_ID"))
  dat_dit <- dat_dit %>% # add corresponding site ID
    left_join(subset(metr_dit, select=c(SAMPLE_ID, SITE_ID), by="SAMPLE_ID"))
  dat_dit <- dat_dit %>% # add corresponding catchment
    left_join(subset(site_dit, select=c(SITE_ID, CATCHMENT, WATER_BODY, FULL_EASTING, FULL_NORTHING, WFD_WATERBODY_ID), by="SITE_ID"))
  dat_dit <- dat_dit %>% # add rank of taxonomic record
    left_join(subset(taxaRef, select=c(TAXON_LIST_ITEM_KEY, TAXON_RANK), by="TAXON_LIST_ITEM_KEY"))

  # Remove water bodies with no Management catchment allocation - only needed if management catchment information required for downstream analysis
  dat_dit <- dat_dit[dat_dit$WFD_WATERBODY_ID %in% wfd_id$Waterbody.ID, ]
  
  dat_dit <- dat_dit %>%
    filter(TAXON_RANK %in% c("Species", "Genus"))
  dat_dit$site <- as.numeric(factor(dat_dit$SITE_ID))
  dat_dit$taxon <- as.numeric(factor(dat_dit$TAXON_LIST_ITEM_KEY))
  dat_dit$UNITS_FOUND_OR_SEQUENCE_READS[is.na(dat_dit$UNITS_FOUND_OR_SEQUENCE_READS)] <- 0
  dat_dit$year <- stringr::str_extract(dat_dit$SAMPLE_DATE, "\\d{4}")
  colnames(dat_dit)[grep("WFD_WATERBODY_ID", colnames(dat_dit))] <- "Waterbody.ID"
  dat_dit <- dat_dit %>%
    left_join(subset(wfd_id, select=c(River.basin.district, 
                                      Management.catchment, 
                                      Operational.catchment, 
                                      Waterbody.ID,
                                      Easting,
                                      Northing)),
              by = "Waterbody.ID")
  dat_dit$count <- dat_dit[,which(colnames(dat_dit) %in% c("TOTAL_ABUNDANCE", "PERCENTAGE_COVER_BAND", "UNITS_FOUND_OR_SEQUENCE_READS"))]
  write.csv(dat_dit, "Diatoms/DIAT_DATA_CLEANED.csv", row.names=F)
} else {
  dat <- read.csv("Diatoms/DIAT_DATA_CLEANED.csv")
}

## Compute area and mean distance of catchment-year surveys
if (0) {
  dat_area_dist <- area_dist(dat, 
                             cores=NULL,
                             catch_max=NULL,
                             catch="operational")
  str(dat_area_dist)
  write.csv(dat_area_dist, "summary_data/diatom_dat_area_dist_op.csv", row.names=F)
} else {
  dat_area_dist <- read.csv("summary_data/diatom_dat_area_dist_op.csv")
}

## Compute temporal mean and standard deviation of the SOD
if (0) {
  dat_temporalSOD <- temporalSOD(dat,
                                 filter=c(y0 = 1990, t0 = 0, n0 = 0, x0 = 10), 
                                 cores=NULL,
                                 catch_max=NULL,
                                 catch="operational")
  str(dat_temporalSOD)
  write.csv(dat_temporalSOD, "summary_data/diatom_temporalSOD_op.csv", row.names=F)
} else {
  dat_temporalSOD <- read.csv("summary_data/diatom_temporalSOD_op.csv")
}

## Compute skewness SOD management catchment
if (0) {
  dat_skewnessSOD <- skewnessSOD(dat, 
                                 cores=NULL,
                                 catch_max=NULL,
                                 catch="rbd")
  str(dat_skewnessSOD)
  write.csv(dat_skewnessSOD, "summary_data/diatom_dat_skewnessSOD_RB.csv", row.names=F)
} else {
  dat_skewnessSOD <- read.csv("summary_data/diatom_dat_skewnessSOD.csv")
}

if (0) {
  dat_skew_accumulation <- skew_accumulation(dat,
                                             filter=c(y0 = 1990, t0 = 0, n0 = 0, x0 = 10), 
                                             no_rand=1000,
                                             cores=NULL,
                                             catch_max=NULL,
                                             catch="rbd")
  str(dat_skew_accumulation)
  write.csv(dat_skew_accumulation, "summary_data/diatom_dat_skew_accumulation_RB.csv", row.names=F)
} else {
  dat_skew_accumulation <- read.csv("summary_data/diatom_dat_skew_accumulation.csv")
}

if (0) {
  dat_temporalMoments <- temporalMoments(dat_mom = dat_skewnessSOD,
                                         dat_rand = dat_skew_accumulation,
                                         filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 8, thresh = -2))
  str(dat_temporalMoments)
  write.csv(dat_temporalMoments, "summary_data/diatom_dat_temporalMoments_RB.csv", row.names=F)
} else {
  dat_temporalMoments <- read.csv("summary_data/diatom_dat_temporalMoments_OC.csv")
}

if (0) {
  dat_turnover <- turnover(dat,
                           filter=c(y0 = 1990, t0 = 0, n0 = 0, x0 = 10), 
                           cores=NULL,
                           catch_max=NULL,
                           catch="operational")
  
  str(dat_turnover)
  write.csv(dat_turnover, "summary_data/diatom_dat_turnover.csv", row.names=F)
  
  dat_turnoverSOD_raw <- turnoverSOD_raw(dat,
                                         filter=c(y0 = 1990, t0 = 0, n0 = 0, x0 = 10), 
                                         cores=NULL,
                                         catch_max=NULL,
                                         catch="operational")
  str(dat_turnoverSOD_raw)
  write.csv(dat_turnoverSOD_raw, "summary_data/diatom_dat_turnoverSOD_raw_jaccard.csv", row.names=F)
} else {
  dat_turnover <- read.csv("summary_data/diatom_dat_turnover.csv")
}

if (0) {
  dat_turnoverLocal <- turnoverLocal(dat,
                                     filter=c(y0 = 1990, t0 = 10, n0 = 0, x0 = 5), 
                                     cores=NULL,
                                     catch_max=NULL,
                                     catch="management")
  str(dat_turnoverLocal)
  write.csv(dat_turnoverLocal, "summary_data/diatom_dat_turnoverLocal_betapart.csv", row.names=F)
} else {
  dat_turnoverLocal <- read.csv("summary_data/diatom_dat_turnoverLocal.csv")
}

if (0) {
  dat_turnoverSOD <- turnoverSOD(dat,
                                 filter=c(y0 = 1990, t0 = 0, n0 = 0, x0 = 10), 
                                 cores=NULL,
                                 catch_max=NULL,
                                 catch="operational")
  str(dat_turnoverSOD)
  write.csv(dat_turnoverSOD, "summary_data/diatom_dat_turnoverSOD.csv", row.names=F)
} else {
  dat_turnoverSOD <- read.csv("summary_data/diatom_dat_turnoverSOD.csv")
}

if (0) {
  if (0) {
    # dat <- subset(dat, Management.catchment %in% dat_temporalMoments$catchment)
    dat_fitTimescales <- timescales(dat,
                                    dat_rand = dat_skew_accumulation,
                                    filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 8, thresh = -2),
                                    cores=NULL,
                                    catch_max=NULL,
                                    catch="management")
    str(dat_fitTimescales[[1]])
    write_rds(dat_fitTimescales, "summary_data/diatoms_dat_fitTimescales_alpha_omega_c.RDS")  
  } else {
    dat_fitTimescales <- read_rds("summary_data/diatoms_dat_fitTimescales_alpha_omega_c.RDS")
  }
  
  dat_fitTimescales_ranges <- data.frame(catchment=rep(NA, length(dat_fitTimescales)),
                                         d.mn=NA, d.lw=NA, d.up=NA,
                                         popDetect.mn=NA, popDetect.lw=NA, popDetect.up=NA,
                                         alpha.mn=NA, alpha.lw=NA, alpha.up=NA,
                                         al.mn=NA, al.lw=NA, al.up=NA,
                                         metapopDetect.mn=NA, metapopDetect.lw=NA, metapopDetect.up=NA,
                                         pRet.mn=NA, pRet.lw=NA, pRet.up=NA, 
                                         pRet.p.mn=NA, pRet.p.lw=NA, pRet.p.up=NA,
                                         m.mn=NA, m.lw=NA, m.up=NA,
                                         cor.al = NA)
  
  # dat_pRet <- c()
  for (cc in 1:length(dat_fitTimescales)) {
    if (typeof(dat_fitTimescales[[cc]])=="list") {
      # dat_fitTimescales_ranges$catchment[cc] <- as.character(dat_fitTimescales[[cc]]$catchment)
      dat_fitTimescales_ranges[cc, 2:(ncol(dat_fitTimescales_ranges)-1)] <- as.vector(t(dat_fitTimescales[[cc]]$ranges))
      dat_fitTimescales_ranges[cc, ncol(dat_fitTimescales_ranges)] <- dat_fitTimescales[[cc]]$covar[4,4]
    }
    # dat_pRet <- rbind(dat_pRet, data.frame(catchment=as.character(dat_fitTimescales[[cc]]$catchment),
    #                                        pRet=dat_fitTimescales[[cc]]$core$pRet,
    #                                        pRet.p=dat_fitTimescales[[cc]]$core$pRet.p))
  }
  dat_fitTimescales_ranges <- dat_fitTimescales_ranges[complete.cases(dat_fitTimescales_ranges),]
  write.csv(dat_fitTimescales_ranges, "summary_data/diatom_dat_fitTimescales_ranges_alpha_omega_c.csv", row.names=F)
  # write.csv(dat_pRet, "summary_data/diatom_dat_pRet.csv", row.names=F)
} else {
  dat_fitTimescales <- read.csv("summary_data/diatom_dat_fitTimescales_alpha_omega_c.csv")
}

if (0) {
  dat_l_n_occupancy <- l_n_occupancy(dat,
                                     cores=NULL,
                                     catch_max=NULL,
                                     catch="operational")
  str(dat_l_n_occupancy)
  write.csv(dat_l_n_occupancy, "summary_data/diatom_dat_l_n_occupancy.csv", row.names=F)
} else {
  dat_l_n_occupancy <- read.csv("summary_data/diatom_dat_l_n_occupancy.csv")
}

## Time series species richness
if (0) {
  dat_richnessLocal <- richnessLocal(dat,
                                     filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 5, i=20))
  str(dat_richnessLocal)
  write.csv(dat_richnessLocal, "summary_data/diatom_dat_richnessLocal_rare.csv", row.names=F)
  
  ggplot(subset(dat_richnessLocal, catchment %in% sample(unique(dat_richnessLocal$catchment), 10)), 
         aes(x=count.tot, y=S, col=factor(paste(catchment, site)))) +
    geom_point() +
    geom_smooth(method='lm', se=F, lwd=0.5) +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') +
    theme(legend.position = "none")
  
  ggplot(subset(dat_richnessLocal, catchment %in% sample(unique(dat_richnessLocal$catchment), 10)), 
         aes(x=year, y=count.tot, col=factor(paste(catchment,site)))) +
    geom_point() +
    geom_smooth(method='lm', se=F) +
    scale_y_continuous(trans='log10') +
    theme(legend.position = "none")
  
  if (0) {
    dat_richnessLocal <- read.csv("summary_data/diatom_dat_richnessLocal.csv")
  }
  
  dat_richnessLocal$rep <- paste0(dat_richnessLocal$catchment, dat_richnessLocal$site)
  dat_richnessLocal$rep <- as.numeric(factor(dat_richnessLocal$rep))
  
  png("figures/local_richness_diatom.png", 5, 5, units = "in", res=300)
  set.seed(3)
  ggplot(subset(dat_richnessLocal, rep %in% sample(unique(dat_richnessLocal$rep), 100)), 
         aes(x=year, y=S, col=factor(paste(catchment, site)))) +
    geom_point() +
    geom_smooth(method="lm", se=F) +
    # geom_smooth(se=F) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = "none") +
    labs(x="Year", y="Richness", title="Diatoms")
  dev.off()

  sl <- c()
  for (x in 1:length(unique(dat_richnessLocal$rep))) {
    dat_sub <- subset(dat_richnessLocal, rep==x)
    dat_sub <- subset(dat_sub, year < min(dat_sub$year)+10)
    # dat_sub <- subset(dat_richnessLocal, rep==x & year < 2010)
    if (nrow(dat_sub) < 5) {
      next
    }
    mod <- lm(S~year, dat_sub)
    sl <- rbind(sl, data.frame(p=summary(mod)$coefficients[2,4],
                               sl=summary(mod)$coefficients[2,1]))
  }
  
  length(which(sl$p<0.05))/nrow(sl)
  mean(sl$sl)
  par(mfrow=c(1,2))
  hist(sl$p)
  hist(sl$sl)
  
  ## 22% of sites show significant linear trend in richness over time when years > 2010 included
  ## mean linear slope 3.74 suggesting an increase in 4 species every year or so when years > 2010 included
  ## 32% of sites show significant linear trend in richness over time when years > 2010 excluded
  ## mean linear slope 1.15 suggesting an increase in one species every year or so when years > 2010 excluded
  ## 18% of sites show significant linear trend in richness over time when first 10 years only used
  ## mean linear slope 3.42 suggesting an increase in one species every year or so when first 10 years only used
  
  mkS <- c()
  for (r in sort(unique(dat_richnessLocal$rep))) {
    dat_sub <- subset(dat_richnessLocal, rep==r)
    dat_sub <- subset(dat_sub, year < min(dat_sub$year)+10)
    # dat_sub <- subset(dat_richnessLocal, rep==r & year < 2010)
    
    if(nrow(dat_sub) < 5) {
      next
    }
    
    S_na_impute <- rep(NA, 1 + max(dat_sub$year) - min(dat_sub$year))
    S_na_impute[1+dat_sub$year-min(dat_sub$year)] <- dat_sub$S
    tt <- min(dat_sub$year):max(dat_sub$year)
    
    S_na_impute <- na_interpolation(S_na_impute)
    mod.mk <- mk.test(S_na_impute)
    mod.sens <- sens.slope(S_na_impute)
    int.sens <- median(S_na_impute)-(mod.sens$estimates)*median(tt-min(tt))
    mkS <- rbind(mkS, 
                 data.frame(rep=r,
                            mk=mod.mk$p.value,
                            sens=int.sens))
  }
  
  length(which(mkS$mk < 0.05))/nrow(mkS)
  mean(mkS$sens)
  hist(mkS$mk)
  hist(mkS$sens)
  
  ## 47% of sites significant trend according to MK (check that p value is sufficient here) when years > 2010 included
  ## 31% of sites significant trend according to MK (check that p value is sufficient here) when years > 2010 excluded
  ## 29% of sites significant trend according to MK (check that p value is sufficient here) when first 10 years are included
  
} else {
  
}

if (0) {
  dat_chaoEstimatorMC <- chaoEstimatorMC(dat,
                                         dat_rand=dat_skew_accumulation,
                                         filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10),
                                         catch="management")
  str(dat_chaoEstimatorMC)
  write.csv(dat_chaoEstimatorMC, "summary_data/diatom_dat_chaoEstimatorMC.csv", row.names=F)
} else {
  dat_chaoEstimatorMC <- read.csv("summary_data/diatom_dat_chaoEstimatorMC.csv")
}