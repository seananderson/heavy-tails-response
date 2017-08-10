###################
#Code to produce results from Youngflesh and Lynch 'Black-swan events - population crashes or temporary emigration?'
#Reply to Anderson et al. 2017 PNAS 'Black-swan events in animal populations'
#
#Primary author: Casey Youngflesh
#Script contains data and some code taken from https://github.com/seananderson/heavy-tails
#
###################



# Clear environment -------------------------------------------------------


# rm(list = ls())


# Load packages -----------------------------------------------------------


require(dplyr)



# Create data -------------------------------------------------------------

setwd('Data')

#load in cleaned GPDD data from Anderson et al. 2017
gpdd <- readRDS("gpdd-clean.rds")

#list of time series IDs
ts_id <- unique(gpdd$main_id)

#calculate max r [where r = log (N_t/N_t-1)] and doubling time for each time series
gpdd_att <- data.frame()
for (i in 1:length(ts_id))
{
  #i <- 1
  temp <- filter(gpdd, main_id == ts_id[i])
  temp_pop <- temp$population_untransformed
  
  temp_r <- rep(NA, nrow(temp))
  for (j in 1:(nrow(temp)-1))
  {
    #j <- 1
    temp_r[j] <- log(temp_pop[j+1]/temp_pop[j])
  }
  
  temp_r_max <- max(temp_r, na.rm = TRUE)
  temp_double_time <- log(2)/temp_r_max
  
  t2 <- data.frame(main_id = ts_id[i],
                   r_max = temp_r_max,
                   double_time = temp_double_time)
  
  gpdd_att <- rbind(gpdd_att, t2)
}





#merge master data with model output p10 (probability nu < 10 - as noted by Anderson et al. 2017 Fig. 2)
#used to ID 'black swan species'
gomp_hat_base <- readRDS("gomp-base-hat.rds")
gpdd_att_p10 <- inner_join(gpdd_att, 
                           gomp_hat_base[,c('main_id', 'p10')], 
                           by = 'main_id')






#get taxonomic data from the gpdd
lookup <- gpdd[,c("main_id", "common_name", "taxonomic_class",
                  "taxon_name", "taxonomic_order", "taxonomic_family")]
lookup <- lookup[!duplicated(lookup), ]

#merge master data with taxonomic data
gpdd_att_p10_taxon <- inner_join(gpdd_att_p10, 
                                 lookup, 
                                 by = "main_id")






#load in life history data from Brook et al. 2006 dataset
brook <- read.csv("brook-etal.csv", stringsAsFactors = FALSE)

brook <- transform(brook, taxon_name = paste(Genus, Species))

#From Brook et al. 2006 supplmentary information
#MinAge = age at first mating in months
#Lifesp = maximum age attained in wild in months
#Fert = number of eggs laid or young born per female per year

#merge master data with life history data
#568 of 609 total time series (93%) used in Anderson et al. 2017 have life history info available 
gpdd_att_p10_taxon_lh <- inner_join(gpdd_att_p10_taxon, 
                                    brook[,c("taxon_name", "MinAge", "Lifesp", "Fert")], 
                                    by = "taxon_name")






# Calculate Rho -----------------------------------------------------------

#add Rho to master data - physiological maximum per-capita growth rate (as defined by Cole 1954)
f_data <- data.frame()
for (i in 1:nrow(gpdd_att_p10_taxon_lh))
{
  
  #filter by each time series
  temp <- gpdd_att_p10_taxon_lh[i,]
  
  a <- temp$MinAge/12 #first age at breeding in months
  b <- temp$Fert/2 #number of females produced per female per year (assume 50:50 sex ratio)
  m <- temp$Lifesp/12 #lifespan of animal in months - max age obtained in wild
  
  
  #from Cole 1954  
  f <- function(r)
  {
    exp(-r) + b*exp(-r*a) - b*exp(-r*(m+1)) - 1
  }
  
  #find root
  Rho <- uniroot(f, lower = 0.00000001, upper = 500, tol = 0.0001)$root
  
  t2 <- data.frame(temp, Rho)
  
  f_data <- rbind(f_data, t2)
}



#f_data is completed data.frame



# Analyze data ------------------------------------------------------------


#which time series have max r value that is greater than Rho
unphysical <- f_data[f_data$r_max > f_data$Rho,]

#proportion of time series that are unphysical
prop_unphysical <- nrow(unphysical)/nrow(f_data) #16% of times series unphysical



#time series that are 'black swans'
bs_data <- filter(f_data, p10 > 0.5)

#black swan time series that are unphysical
bs_unphysical <- bs_data[bs_data$r_max > bs_data$Rho,]

#proportion of black swan time series that are unphysical
bs_prop_unphysical <- nrow(bs_unphysical)/nrow(bs_data) #41% of black swan time series unphysical






# Lifespan ----------------------------------------------------------------


#compare mean life span to max life span
#Lynch and Fagan 2009 Appendix A

mean_vals <- c(2.2, 1.36, 2.02, 4.59, 0.59, 0.45, 0.21, 1.25, 0.26, 3.46, 1.17, 2.76,
               1.14, 0.81, 0.98, 1.21, 1.30, 0.84, 1.20, 0.25, 0.685, 1.07, 1.07, 0.36,
               1.13, 1.57, 1.68, 1.77, 2.28, 1.47, 3.76, 2.50, 4.09, 29.35, 14.53, 4.64, 3.43, 2.53,
               3.47, 5.25, 2.16, 3.36, 1.79, 2.08, 4.89, 5.81, 5.04, 12.1, 3.85, 2.31, 5.19,
               7.95, 6.39, 3.05, 16.25, 3.79, 2.45, 22.44)

max_vals <- c(8.0, 9.5, 11.0, 20.1, 3.9, 1.8, 1.58, 7.0, 1.58, 18.0, 6.0, 8.0, 8.0, 7.0,
              4.0, 5.5, 10.0, 4.75, 7.0, 1.1, 1.75, 6.0, 8.0, 3.3, 12.0, 9.0, 8.0, 6.0,
              11.0, 12.0, 21.0, 10.0, 11.0, 59.0, 56.0, 24.0, 13.0, 13.0, 13.0, 14.0, 11.1,
              19.0, 17.0, 7.0, 21.0, 16.0, 21.0, 36.0, 13.0, 10.0, 16.0, 21.0, 32.8, 15.0, 
              24.0, 19.0, 24.0, 44.0)

mean(max_vals/mean_vals) 

#max lifespan is on average 5 times greater than mean lifespan for mammals in the wild
#one reason for gross overestimate of rho (physiological max per-capita growth rate)
#other reason is that females don't give birth on first day of each year





# Red Grouse -------------------------------------------------------

#Red grouse - 10128 - 16-fold increase in abundance - Potts et al. 1984 source
RG <- filter(gpdd, main_id == 10128)
# plot(RG$series_step, RG$population_untransformed, type ='b')
# cbind(RG$sample_year, RG$population_untransformed)

save(f_data, unphysical, bs_data, bs_unphysical, file = "../../generated-data/youngflesh.rda")
