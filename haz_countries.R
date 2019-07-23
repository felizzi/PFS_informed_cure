### script to extract background mortality for a given country


country_map <- new.env(hash=T, parent=emptyenv())
c_mort_list <- read.csv("libraries/CountryList.csv")

code3 <- c('AUS','BEL','CAN','CHL','CZE','EST','FIN','FRA','DEU','HUN','ITA','JPN','NLD','NOR','POL','PRT','SVK','SVN','ESP','SWE','TWN','USA','GBR')

for (i in 1:length(c_mort_list$Code)){
  country_use <- as.character(c_mort_list$Code3[i])
  country_map[[country_use]] <- as.character(c_mort_list$Country[i])
}


### values for 
UTILITY_INTERCEPT <- 0.9508566 
UTILITY_AGE <- -0.0002587
UTILITY_AGE_2 <- -0.0000332  
UTILITY_MALE <- 0.0212126   
#male_utility <-  UTILITY_INTERCEPT +  UTILITY_AGE*Age + UTILITY_AGE_2*Age^2 + UTILITY_MALE


## add additional keys to the map NED SUI GER
country_map[["SUI"]] <- "Switzerland"
country_map[["NED"]] <- "Netherlands"
country_map[["GER"]] <- "Germany"
source("sasmacro/functions_long_term_survival.r")
hazard_time <- function(table_, evttme, sex, age, year, country_trial, country_output){
  # load in the file the tables 
  ## loop over the countries in the trial 
  Np <- length(table_[,sex]) ## number of rows in the datatable
  
  rate_vec <- matrix(,nrow = Np, ncol = 1);
  surv_vec <- matrix(,nrow = Np, ncol = 1);
  
  xty <-7*seq(0,3824, by = 1)/365.25 ## timespan -- 30 years 
  surv_mat <- matrix(,nrow = Np, ncol = length(xty));
  utility_vec <- matrix(,nrow = Np, ncol = length(xty));
  
  
  Nxty <- length(xty)
  #surv_mean_table_all <- data.frame(week = numeric(Nxty), residual_surv = numeric(Nxty))
  surv_mean_table <- data.frame(week = numeric(Nxty), residual_surv = numeric(Nxty))
  
  for (i in 1:Np){ ### loop over the subjects in the trial 
    ## manual adaptations of the country-codes
    if (!is.null(country_output)){
      country_in_use <- country_output
    }else{
      country_in_use <- as.character(table_[i,country_trial])  
    }
    if (!(country_in_use %in% code3)){ country_in_use <- "USA"}
    year_in_use   <- table_[i, year]
    gender_in_use <- table_[i,sex]

    ## load the tables 
    if (gender_in_use == "MALE" | gender_in_use == "M"){
      table_in_use <- read.csv(paste0("libraries/",country_map[[country_in_use]],"/Male/Mortality.csv", sep = ""))
      max_year <- max(table_in_use$Year) ## selection of the max year in case the year in the trial is missing
      if (year_in_use > max_year){ year_in_use <- max_year}
      table_in_use <- subset(table_in_use, Year == year_in_use)
      rate_vec[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/12, table_in_use);
      scvy_all <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use);
      
      #calculate here the utility for the person, starting ate table_[i,age] and continuing up to xty 
      utility_vec[i,] <- UTILITY_AGE*(as.numeric(table_[i,age]) + xty) + UTILITY_AGE_2*(as.numeric(table_[i,age]) + xty)^2 + UTILITY_MALE + UTILITY_INTERCEPT 
    }
    if (gender_in_use == "FEMALE" | gender_in_use == "F"){
      table_in_use <- read.csv(paste0("libraries/",country_map[[country_in_use]],"/Female/Mortality.csv", sep = ""))
      max_year <- max(table_in_use$Year) ## selection of the max year in case the year in the trial is missing
      if (year_in_use > max_year){ year_in_use <- max_year}
      table_in_use <- subset(table_in_use, Year == year_in_use)
      rate_vec[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/12, table_in_use);
      
      scvy_all <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use);
      utility_vec[i,] <-  UTILITY_AGE*(as.numeric(table_[i,age]) + xty) + UTILITY_AGE_2*(as.numeric(table_[i,age]) + xty)^2 + UTILITY_INTERCEPT 
    }
   
    surv_mat[i,] <- scvy_all/scvy_all[1];
    #  print("MALE")
    #rate_vec[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/365.25, table_in_use_M);
    #    rate_vec_USA[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/12, USA_table_male_2013)
    #   surv_vec_USA[i] <- Surv_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/12, USA_table_male_2013)/Surv_c_years(as.numeric(table_[i,age]), USA_table_male_2013)
    #scvy_all <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use_M);
    #    scvy_usa <- Surv_cv_years( as.numeric(table_[i,age]) + xty, USA_table_male_2013)
    #scvy_uk <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use_UK_M);
    
    # }else{
    #    print("FEMALE")
    #rate_vec[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/365.25, table_in_use_F)
    #    rate_vec_USA[i] <- rate_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/12, USA_table_female_2013)
    #   surv_vec_USA[i] <- Surv_c_years(as.numeric(table_[i,age]) + as.numeric(table_[i,evttme])/12, USA_table_female_2013)/Surv_c_years(as.numeric(table_[i,age]), USA_table_female_2013)
    #scvy_all <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use_F)
    #    scvy_usa <- Surv_cv_years( as.numeric(table_[i,age]) + xty, USA_table_female_2013)
    #scvy_uk <- Surv_cv_years( as.numeric(table_[i,age]) + xty, table_in_use_UK_F);
  }
  
  #surv_mat_all[i,] <- scvy_all/scvy_all[1];
  # surv_mat_usa[i,] <- scvy_usa/scvy_usa[1];
  #surv_mat_uk[i,] <- scvy_uk/scvy_uk[1];
  #}
  
  
  
  ###########################################################
  #mean survival curve 
  #surv_mean_all <- colSums(surv_mat_all)/Np
  #surv_mean_usa <- colSums(surv_mat_usa)/Np
  #surv_mean_uk <- colSums(surv_mat_uk)/Np
  #print("dimension of SURV_MAT UK")
  #print(dim(surv_mat_uk))
  
  ret_obj <- list()
  
  had_mat <- surv_mat*utility_vec
  ret_obj$rate_vec <- rate_vec
  ret_obj$surv_mat <- surv_mat
  
  ###########################################################
  #mean survival curve 
  surv_mean <- colSums(surv_mat)/length(table_[,sex])
  ret_obj$surv_mean <- surv_mean
  ret_obj$xty <- xty
  ret_obj$utility_mat <- utility_vec
  ret_obj$had_mat <- had_mat
  
  ret_obj$had_mean <- colSums(had_mat)/length(table_[,sex])
  #ret_obj$surv_mean_all <- surv_mean_all
  #ret_obj$surv_mean_uk <- surv_mean_uk
  #ret_obj$surv_mean_usa <- surv_mean_usa
  #ret_obj$time <- xty
  return(ret_obj)
}