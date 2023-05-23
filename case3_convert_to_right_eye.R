###Start with a clean space
rm(list = ls())

###Load packages
library(stringr)
library(dplyr)
library(womblR)
library(tidyverse)
###Load Rotterdam data and combine
vf <- read.csv("glaucoma_data/VFPoints.csv")
vfs <- read.csv("glaucoma_data/VisualFields.csv")
pats <- read.csv("glaucoma_data/Patients.csv")
vf <- left_join(vf, vfs, by = "FIELD_ID")
vf <- left_join(vf, pats, by = "STUDY_ID")
vf$XY <- paste0(vf$X, "_", vf$Y)

###Split data into right and left eyes to process before combining

###Process right eyes
vf_od <- vf[vf$SITE == "OD", ]
od_xy <- unique(cbind(vf_od$X, vf_od$Y))
od_xy <- od_xy[order(od_xy[, 1], decreasing = FALSE), ]
od_xy <- od_xy[order(od_xy[, 2], decreasing = TRUE), ]
conv_od <- data.frame(Location = 1:54, XY = paste0(od_xy[, 1], "_", od_xy[, 2]))
vf_od <- left_join(vf_od, conv_od, by = "XY")

###Process left eyes
vf_os <- vf[vf$SITE == "OS", ]
os_xy <- unique(cbind(vf_os$X, vf_os$Y))
os_xy <- os_xy[order(os_xy[, 1], decreasing = TRUE), ]
os_xy <- os_xy[order(os_xy[, 2], decreasing = TRUE), ]
conv_os <- data.frame(Location = 1:54, XY = paste0(os_xy[, 1], "_", os_xy[, 2]))
vf_os <- left_join(vf_os, conv_os, by = "XY")

###Combine left and right eyes
vf <- rbind(vf_od, vf_os)
vf$X <- vf$Y <- vf$XY <- NULL
vf <- vf[order(vf$FIELD_ID, vf$Location), ]

###Output data (all data have been converted to right eyes)
dat <- with(vf, data.frame(STUDY_ID, FIELD_ID, SITE, AGE, SEX, IOP, MD, LOCATION = Location, THRESHOLD, TOTAL_DEVIATION))

###Add visual field region
region <- numeric(length = 54)
region[c(1:5, 9, 10, 17)] <- "in" # inferior nasal
region[c(6:8, 11:16, 19:22)] <- "it" # inferior temporal
region[c(23:25, 32:34)] <- "t" # temporal
region[c(29:31, 38:42, 47:48)] <- "st" # superior temporal
region[c(28, 37, 43, 45:46, 49:54)] <- "sn" # superior nasal
region[c(18, 27, 36, 44)] <- "n" # nasal
region[region == "0"] <- NA
dat <- left_join(dat, data.frame(LOCATION = 1:54, REGION = region, DEGREE = womblR::GarwayHeath), by = "LOCATION")

###Order properly and compute a eye specific id
dat <- dat[order(dat$STUDY_ID, dat$SITE, dat$AGE, dat$LOCATION), ]
dat$ID <- cumsum(!duplicated(dat[, c("STUDY_ID", "SITE")]))

###Compute time
dat$TIME <- unlist(with(dat, tapply(AGE, ID, function(x) x - x[1]))) / 365

###Create data with good names
dat <- with(dat, data.frame(id = ID, 
                    patient = STUDY_ID,
                    eye = SITE,
                    age = AGE,
                    sex = SEX,
                    iop = IOP,
                    md = MD,
                    location = LOCATION,
                    dls = THRESHOLD,
                    td = TOTAL_DEVIATION,
                    time = TIME,
                    region = REGION,
                    degree = DEGREE))

###Now that the locations are integers 1:54, you can use the adjacency matrices available in womblR, for example
W <- womblR::HFAII_Rook

dat<- dat %>%
  mutate(years = age/365.25)
