## This script will compute the skew of the OFD at the management catchment scale and run various spatio-temporal analyses on result
## Bulk data here: https://environment.data.gov.uk/ecology-fish/downloads/
## Analysis of invert data

require(tidyverse)
setwd("/data/home/btx718/Paleolimnology/datasets/EA/") # set working directory to EA-OFD repository
source("scripts/analysis_functions.R")

################################################################################
################################ Macroinverts ##################################
################################################################################

if (0) {
  
  ## Assemble dataset (first unzip files)
  metr_inv <- as_tibble(read.csv("Macroinverts/INV_OPEN_DATA_METRICS.csv"))
  site_inv <- as_tibble(read.csv("Macroinverts/INV_OPEN_DATA_SITE.csv"))
  taxa_inv <- as_tibble(read.csv("Macroinverts/INV_OPEN_DATA_TAXA.csv"))
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
  
  ## Skewness POD macroinverts (species+genus)
  dat_inv <- taxa_inv %>% # select relevant response data
    select(ANALYSIS_ID, TAXON_LIST_ITEM_KEY, TOTAL_ABUNDANCE)
  dat_inv <- dat_inv %>% # add corresponding sample ID
    left_join(subset(metr_inv, select=c(ANALYSIS_ID, SAMPLE_ID, SAMPLE_DATE), by="ANALYSIS_ID"))
  dat_inv <- dat_inv %>% # add corresponding site ID
    left_join(subset(metr_inv, select=c(SAMPLE_ID, SITE_ID), by="SAMPLE_ID"))
  dat_inv <- dat_inv %>% # add corresponding catchment
    left_join(subset(site_inv, select=c(SITE_ID, CATCHMENT, WATER_BODY, FULL_EASTING, FULL_NORTHING, WFD_WATERBODY_ID), by="SITE_ID"))
  dat_inv <- dat_inv %>% # add rank of taxonomic record
    left_join(subset(taxaRef, select=c(TAXON_LIST_ITEM_KEY, TAXON_RANK), by="TAXON_LIST_ITEM_KEY"))
  
  # Remove water bodies with no Management catchment allocation - only needed if management catchment information required for downstream analysis
  dat_inv <- dat_inv[dat_inv$WFD_WATERBODY_ID %in% wfd_id$Waterbody.ID, ]
  
  dat_inv <- dat_inv %>%
    filter(TAXON_RANK %in% c("Species", "Genus"))
  dat_inv$site <- as.numeric(factor(dat_inv$SITE_ID))
  dat_inv$taxon <- as.numeric(factor(dat_inv$TAXON_LIST_ITEM_KEY))
  dat_inv$TOTAL_ABUNDANCE[is.na(dat_inv$TOTAL_ABUNDANCE)] <- 1
  dat_inv$year <- stringr::str_extract(dat_inv$SAMPLE_DATE, "\\d{4}")
  colnames(dat_inv)[grep("WFD_WATERBODY_ID", colnames(dat_inv))] <- "Waterbody.ID"
  dat_inv <- dat_inv %>%
    left_join(subset(wfd_id, select=c(River.basin.district, 
                                      Management.catchment, 
                                      Operational.catchment, 
                                      Waterbody.ID,
                                      Easting,
                                      Northing)),
              by = "Waterbody.ID")
  dat_inv$count <- dat_inv[,which(colnames(dat_inv) %in% c("TOTAL_ABUNDANCE", "PERCENTAGE_COVER_BAND", "UNITS_FOUND_OR_SEQUENCE_READS"))]
  write.csv(dat_inv, "Macroinverts/INV_DATA_CLEANED.csv", row.names=F)
} else {
  dat <- read.csv("Macroinverts/INV_DATA_CLEANED.csv")
}

## Compute area and mean distance of catchment-year surveys
if (0) {
  dat_area_dist <- area_dist(dat, 
                             cores=NULL,
                             catch_max=NULL,
                             catch="management")
  str(dat_area_dist)
  write.csv(dat_area_dist, "summary_data/macroinvert_dat_area_dist_op.csv", row.names=F)
} else {
  dat_area_dist <- read.csv("summary_data/macroinvert_dat_area_dist.csv")
}

## Compute temporal mean and standard deviation of the SOD
if (0) {
  dat_temporalSOD <- temporalSOD(dat,
                                 filter=c(y0 = 1990, t0 = 0, n0 = 0, x0 = 10), 
                                 cores=NULL,
                                 catch_max=NULL,
                                 catch="management")
  str(dat_temporalSOD)
  write.csv(dat_temporalSOD, "summary_data/macroinvert_dat_temporalSOD.csv", row.names=F)
} else {
  dat_temporalSOD <- read.csv("summary_data/macroinvert_dat_temporalSOD.csv")
}

## Compute skewness SOD management catchment
if (0) {
  dat_skewnessSOD <- skewnessSOD(dat, 
                                 cores=NULL,
                                 catch_max=NULL,
                                 catch="rbd")
  str(dat_skewnessSOD)
  write.csv(dat_skewnessSOD, "summary_data/macroinvert_dat_skewnessSOD.csv", row.names=F)
} else {
  dat_skewnessSOD <- read.csv("summary_data/macroinvert_dat_skewnessSOD.csv")
}

if (0) {
  dat_skew_accumulation <- skew_accumulation(dat,
                                             filter=c(y0 = 1990, t0 = 0, n0 = 0, x0 = 10), 
                                             no_rand=1000,
                                             cores=NULL,
                                             catch_max=NULL,
                                             catch="rbd")
  str(dat_skew_accumulation)
  write.csv(dat_skew_accumulation, "summary_data/macroinvert_dat_skew_accumulation.csv", row.names=F)
} else {
  dat_skew_accumulation <- read.csv("summary_data/macroinvert_dat_skew_accumulation.csv")
}

if (0) {
  dat_temporalMoments <- temporalMoments(dat_mom = dat_skewnessSOD,
                                         dat_rand = dat_skew_accumulation,
                                         filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 8, thresh = -2))
  str(dat_temporalMoments)
  write.csv(dat_temporalMoments, "summary_data/macroinvert_dat_temporalMoments.csv", row.names=F)
} else {
  dat_temporalMoments <- read.csv("summary_data/macroinvert_dat_temporalMoments.csv")
}

if (0) {
  dat_turnover <- turnover(dat,
                           filter=c(y0 = 1990, t0 = 0, n0 = 0, x0 = 10), 
                           cores=NULL,
                           catch_max=NULL,
                           catch="management")
  str(dat_turnover)
  write.csv(dat_turnover, "summary_data/macroinvert_dat_turnover.csv", row.names=F)
} else {
  dat_turnover <- read.csv("summary_data/macroinvert_dat_turnover.csv")
}

if (0) {
  dat_turnoverLocal <- turnoverLocal(dat,
                                     filter=c(y0 = 1990, t0 = 10, n0 = 0, x0 = 5), 
                                     cores=NULL,
                                     catch_max=NULL,
                                     catch="management")
  str(dat_turnoverLocal)
  write.csv(dat_turnoverLocal, "summary_data/macroinvert_dat_turnoverLocal.csv", row.names=F)
} else {
  dat_turnoverLocal <- read.csv("summary_data/macroinvert_dat_turnoverLocal.csv")
}

if (0) {
  dat_turnoverSOD <- turnoverSOD(dat,
                                 filter=c(y0 = 1990, t0 = 0, n0 = 0, x0 = 10), 
                                 cores=NULL,
                                 catch_max=NULL,
                                 catch="management")
  str(dat_turnoverSOD)
  write.csv(dat_turnoverSOD, "summary_data/macroinvert_dat_turnoverSOD.csv", row.names=F)
  
  dat_turnoverSOD_raw <- turnoverSOD_raw(dat,
                                         filter=c(y0 = 1990, t0 = 0, n0 = 0, x0 = 10), 
                                         cores=NULL,
                                         catch_max=NULL,
                                         catch="management",
                                         method="jaccard")
  str(dat_turnoverSOD_raw)
  write.csv(dat_turnoverSOD_raw, "summary_data/macroinvert_dat_turnoverSOD.csv", row.names=F)
} else {
  dat_turnoverSOD <- read.csv("summary_data/macroinvert_dat_turnoverSOD.csv")
}

if (0) {
  if (0) {
    dat_fitTimescales <- timescales(dat,
                                    dat_rand = dat_skew_accumulation,
                                    filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 8, thresh = -2),
                                    cores=20,
                                    catch_max=NULL,
                                    catch="management")
    str(dat_fitTimescales[[1]])
    write_rds(dat_fitTimescales, "summary_data/macroinvert_dat_fitTimescales_alpha_omega_c2.RDS")
  } else {
    dat_fitTimescales <- read_rds("summary_data/macroinvert_dat_fitTimescales_alpha_omega.RDS")
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
                                         cor.al = NA,
                                         min_diff=NA,
                                         max_diff=NA)
  
  for (cc in 1:length(dat_fitTimescales)) {
    if (typeof(dat_fitTimescales[[cc]])=="list") {
      dat_fitTimescales_ranges$catchment[cc] <- as.character(dat_fitTimescales[[cc]]$catchment)
      dat_fitTimescales_ranges[cc, 2:(ncol(dat_fitTimescales_ranges)-3)] <- as.vector(t(dat_fitTimescales[[cc]]$ranges))
      dat_fitTimescales_ranges[cc, ncol(dat_fitTimescales_ranges)-2] <- dat_fitTimescales[[cc]]$covar[4,4]
      dat_fitTimescales_ranges[cc, (ncol(dat_fitTimescales_ranges)-1)] <- dat_fitTimescales[[cc]]$min_diff
      dat_fitTimescales_ranges[cc, (ncol(dat_fitTimescales_ranges))] <- dat_fitTimescales[[cc]]$max_diff
    }
  }
  dat_fitTimescales_ranges <- dat_fitTimescales_ranges[complete.cases(dat_fitTimescales_ranges),]
  write.csv(dat_fitTimescales_ranges, "summary_data/macroinvert_dat_fitTimescales_ranges", row.names=F)
} else {
  dat_fitTimescales_ranges <- read.csv("summary_data/macroinvert_dat_fitTimescales_ranges.csv")
}

if (0) {
  dat_l_n_occupancy <- l_n_occupancy(dat,
                                     cores=NULL,
                                     catch_max=NULL,
                                     catch="management")
  str(dat_l_n_occupancy)
  write.csv(dat_l_n_occupancy, "summary_data/macroinvert_dat_l_n_occupancy.csv", row.names=F)
} else {
  dat_l_n_occupancy <- read.csv("summary_data/macroinvert_dat_l_n_occupancy.csv")
}

## Time series species richness
if (0) {
  dat_richnessLocal <- richnessLocal(dat,
                                     filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 5, i0=20))
  str(dat_richnessLocal)
  write.csv(dat_richnessLocal, "summary_data/macroinvert_dat_richnessLocal.csv", row.names=F)
  
} else {
  dat_richnessLocal <- read.csv("summary_data/macroinvert_dat_richnessLocal.csv")
}

if (0) {
  dat_chaoEstimatorMC <- chaoEstimatorMC(dat,
                                         dat_rand=dat_skew_accumulation,
                                         filter=c(y0 = 1990, t0 = 10, n0 = 5, x0 = 10),
                                         catch="management")
  str(dat_chaoEstimatorMC)
  write.csv(dat_chaoEstimatorMC, "summary_data/macroinvert_dat_chaoEstimatorMC.csv", row.names=F)
} else {
  dat_chaoEstimatorMC <- read.csv("summary_data/macroinvert_dat_chaoEstimatorMC.csv")
}