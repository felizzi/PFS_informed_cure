# file for informed use 
# select distinct endpoints 
require(MASS)
require(reshape2)
require(tidyr)

end_point <- "OS"
end_point_info <- "PFSINV"

data_cut <- 'i29365h'
#data_cut <- "ash_2019"
data_file <- 'ate.sas7bdat'
#data_cut_m <- "2018 October"
data_cut_m <- "2019 March"

censor_transpl_fl <- FALSE

d_frame <- data.frame("t" = numeric(), "S_control" = numeric(), "S_intervention" = numeric(), "dist" = character(), stringsAsFactors = F)

### if datacut is March 2019, then select data from entimICE
SAICE::initialize_connection(entimice_env = "PROD") # WebUI prod version
event_enti <- SAICE::read_entimice("root/clinical_studies/RO5541077/CDPT7898/GO29365/data_analysis/ema_questions_90dsu/qa/outdata_vad/ateupdt.sas7bdat")
event_ash <- SAICE::read_entimice("root/clinical_studies/RO5541077/CDPT7898/GO29365/data_analysis/ash_2019/qa/outdata_vad/ate.sas7bdat")

data_loc <- paste('/opt/BIOSTAT/qa/cdpt7898', data_cut, 'libraries',  data_file, sep = '/')

event <- read_bce(data_loc) %>% dplyr::filter(PARAMCD == end_point) %>% 
  dplyr::select(AVAL, CNSR, ARMCD, AGE, SEX, COUNTRY, USUBJID ) %>% dplyr::mutate(CNSR = (1+CNSR)%%2) %>% filter(ARMCD %in% c("RAN DL BR", "RAN DL BR+POV"))

if (end_point == "PFSIRC"){
  event <- event_enti %>% dplyr::filter(PARAMCD == "PFSIRC3") %>% 
    dplyr::select(AVAL, CNSR, ARMCD, AGE, SEX, COUNTRY, USUBJID ) %>% dplyr::mutate(CNSR = (1+CNSR)%%2) %>% filter(ARMCD %in% c("RAN DL BR", "RAN DL BR+POV"))
}
if (data_cut == "ash_2019"){
  event <- event_ash %>% dplyr::filter(PARAMCD == end_point) %>% 
    dplyr::select(AVAL, CNSR, ARMCD, AGE, SEX, COUNTRY, USUBJID ) %>% dplyr::mutate(CNSR = (1+CNSR)%%2) %>% filter(ARMCD %in% c("RAN DL BR", "RAN DL BR+POV"))
}
if (end_point == "OS"  && censor_transpl_fl){
  ## in this case, modify the survival
  subj_list <- c("GO29365-273002-4046", "GO29365-273002-4058", "GO29365-272994-4038", "GO29365-272994-4044")
  event <- event %>% dplyr::mutate(CNSR = ifelse(USUBJID %in% subj_list, 0, CNSR)) 
  event <- event %>% dplyr::mutate(AVAL = ifelse(USUBJID == subj_list[1],  198/30.4375, AVAL)) 
  event <- event %>% dplyr::mutate(AVAL = ifelse(USUBJID == subj_list[2],  342/30.4375, AVAL))
  event <- event %>% dplyr::mutate(AVAL = ifelse(USUBJID == subj_list[3],  288/30.4375, AVAL))
  event <- event %>% dplyr::mutate(AVAL = ifelse(USUBJID == subj_list[4],  210/30.4375, AVAL))
}


#source("Pola_cure_endpinv.R")
parms_fit_frame_1 <- data.frame("Function" = character(0), "Param_Name" = character(0), "Value_mix" = numeric(0), "Cov_mix_p1" = numeric(0), "Cov_mix_p2" = numeric(0), "Cov_mix_p3" = numeric(0), "Cov_mix_PI_1" = numeric(0), "Cov_mix_PI_2" = numeric(0), stringsAsFactors=FALSE)
parms_fit_frame_2 <- data.frame("Function" = character(0), "Param_Name" = character(0), "Value_mix" = numeric(0), "Cov_mix_p1" = numeric(0), "Cov_mix_p2" = numeric(0), "Cov_mix_p3" = numeric(0), "Cov_mix_PI_1" = numeric(0), "Cov_mix_PI_2" = numeric(0), stringsAsFactors=FALSE)

event$YEAR <- 2016
event$COUNTRY <- "USA"
haz_res_os <- hazard_time(table_ = as.data.frame(event), evttme = "AVAL", sex = "SEX" , age = "AGE", year = "YEAR", country_trial =  "COUNTRY", country_output = NULL)



Nfr <- length(models)
event <- as.data.frame(event)
Nfr_cov <- 16*4 + 25 + 9 
out_frame_intervention <- data.frame("STUDY NAME" = character(Nfr), "POPULATION" = character(Nfr),	"TREATMENT ARM" = character(Nfr), 
                     "Endpoint" = character(Nfr), "Endpoint_def" = character(Nfr), "Model" = character(Nfr), "Data_cut" = character(Nfr),
                     "Distribution" = character(Nfr),	
                     "rate" = numeric(Nfr), 	"shape" = numeric(Nfr),	"scale" = numeric(Nfr),	"meanlog" = numeric(Nfr), "sdlog" = numeric(Nfr),
              
                     "mu" = numeric(Nfr), "sigma" = numeric(Nfr), "Q" = numeric(Nfr), "PI_1" = numeric(Nfr), "PI_2" = numeric(Nfr),        "aic" = numeric(Nfr),	"bic" = numeric(Nfr), 
                     "Info" = character(Nfr),
                     "Source" = character(Nfr), stringsAsFactors = F)

out_frame_control <- out_frame_intervention


cov_frame_intervention <- data.frame("STUDY NAME" = character(Nfr_cov), "POPULATION" = character(Nfr_cov),	"TREATMENT ARM" = character(Nfr_cov), 
                     "Endpoint" = character(Nfr_cov), "Endpoint_def" = character(Nfr_cov), "Data_cut" = character(Nfr_cov), "Hazard" = character(Nfr_cov),
                     "Distribution" = character(Nfr_cov),	"row_name" = numeric(Nfr_cov),
                     "column_name" = numeric(Nfr_cov), 	"row_num" = numeric(Nfr_cov),	"col_num" = numeric(Nfr_cov),	
                     "data_point" = numeric(Nfr_cov), stringsAsFactors = F)

cov_frame_control <- cov_frame_intervention

#oder_distr <- EXPONENTIAL, WEIBULL, LNORMAL, Gen GAMMA, LLOGISTIC, GOMPERTZ

#models <- c("exp",       "weibullPH", "llogis",    "lnorm",     "gengamma",  "gompertz" )
#standars order preferred by Thuresson 
models <- c("exp",       "weibullPH", "lnorm",     "gengamma", "llogis",  "gompertz" )

event$rate_mod <- haz_res_os$rate_vec


event <- event %>% dplyr::mutate(ARMCDN = ifelse(ARMCD == "RAN DL BR",1,2))
event <- as.data.frame(event)
##### likelihood function 
ll.mix_f <- function(table_, parms, pi_, time_col, cen_col, rate_col, obj){  
  time_v <- as.numeric(table_[,time_col])/12;
  cen_v <- as.numeric(table_[,cen_col]);
  rate_v <- as.numeric(table_[,rate_col]);
  #pi_ <- parms[length(parms)];
  
  #shape <- parms[1];  scale <- parms[2]; pi_ <- parms[3];  
  f_t <- do.call(obj$dfns$d, args = c(as.list(parms[1:length(parms)]),
                                      list (x = time_v, log = FALSE))
  )
  S_t <- do.call(obj$dfns$p, args = c(as.list(parms[1:length(parms)]),
                                      list(q = time_v, lower.tail = FALSE))
  )
  h_ <- rate_v + ((1-pi_)*f_t)/(pi_ + (1-pi_)*S_t);  
  s_ <- pi_ + (1 - pi_)*S_t;  
  ret_val <- sum( cen_v*log(h_) + log(s_));  
  return(-ret_val)
}

# load the relevant file with the info on the right endpoint 
parms_fit_frame_endp <- read.csv(paste0("reports/Parms_mixed_",end_point_info,"_intervention_and_control_", data_cut,".csv"))

out_input_endp <- parms_fit_frame_endp %>% dplyr::filter(Param_Name %in% c("PI_1","PI_2")) %>% dplyr::select("Function","Param_Name", "Value_mix")
out_input_endp <-out_input_endp %>% group_by("Function") %>% spread(Param_Name, Value_mix)
out_input_endp <- as.data.frame(out_input_endp)

### loop over the different distributions 
init_pos <- 1
j0 <- 1
for (i in 1:(length(models) )){
  #for (i in 2:2){
  
  ## apply the method for each arm separately 
  mat_in_use <- parms_fit_frame_endp %>% dplyr::filter(Function == models[i])
  mean_vals <- mat_in_use$Value_mix
  var_vals <- mat_in_use[,4:(4+length(mat_in_use$Param_Name) - 1)]
  parms_corr_1 <- matrix(,nrow = 323, ncol = length(mean_vals) )
  parms_corr_2 <- matrix(,nrow = 323, ncol = length(mean_vals) )
  #determine svd of var_vals. Set the stuff to 0 if the matrix 
  sing_dec <- svd(var_vals)
  
  if (min(sing_dec$d) < 1000*.Machine$double.eps){
    print("entered if small")
    sing_dec$d[which(sing_dec$d<.Machine$double.eps)] <- 1100*.Machine$double.eps
    var_vals <- sing_dec$u%*%diag(sing_dec$d)%*%t(sing_dec$v)
  }
  if (max(sing_dec$d) > 10){ 
    print("entered if large")
    sing_dec$d[which(sing_dec$d>10)] <- 1 #this needs to be clinically justified 
    var_vals <- sing_dec$u%*%diag(sing_dec$d)%*%t(sing_dec$v)
  }
  v_post <- mvrnorm(N_sample, mean_vals, var_vals, tol = 1e-2)
  
  
  #determine here the values that lead to non-plausible values on the cure PI for SOC
  #trt_v <- c(1,0)
  #print(c(b1,b2))
  #logit_arg <- v_post[,3:4]%*%rbind(1,trt_v) 
  #logit_arg <- mean_vals[3:4]%*%rbind(1,trt_v) 
  #pi_v <- exp(logit_arg)/(1 + exp(logit_arg))
  parms_corr_1[,(length(mat_in_use$Param_Name)-1):length(mat_in_use$Param_Name)] <- v_post[1:323,(length(mat_in_use$Param_Name)-1):length(mat_in_use$Param_Name)]
  parms_corr_2[,(length(mat_in_use$Param_Name)-1):length(mat_in_use$Param_Name)] <- v_post[1:323,(length(mat_in_use$Param_Name)-1):length(mat_in_use$Param_Name)]
  
  
  #modify here the cure for the control
  
  v_post[,length(mat_in_use$Param_Name)] <- v_post[,length(mat_in_use$Param_Name)-1] + v_post[,length(mat_in_use$Param_Name)]
  
  v_cure_1 <- exp(v_post[,length(mat_in_use$Param_Name)-1])/(1+exp(v_post[,length(mat_in_use$Param_Name)-1]))
  v_cure_2 <- exp(v_post[,length(mat_in_use$Param_Name)])/(1+exp(v_post[,length(mat_in_use$Param_Name)]))
  
  #if (i >= 4){
  #  v_cure_1 <-rep(0,123)
  #}
  
  #b_use <- as.numeric(out_input_endp[i,3:4])
  b_use <- mean_vals[3:4]
  logit_arg <- b_use%*%rbind(c(1,1),c(0,1)) 
  pi_v <- exp(logit_arg)/(1 + exp(logit_arg))
  print(paste0("  model  ", models[i],    "    ",pi_v  ))
  
  
  for (nt in 1:323){
    
    ### within each model, run a series of simulations varying the uncertainty around cures 
    
    fsr_fits_1 <- flexsurvreg(Surv(as.numeric(AVAL)/12, CNSR)~1, data = subset(event,ARMCDN == 1), dist = models[i])
    fsr_fits_2 <- flexsurvreg(Surv(as.numeric(AVAL)/12, CNSR)~1, data = subset(event,ARMCDN == 2), dist = models[i])
    fsr_use_1 <- fsr_fits_1$res[,'est'];
    fsr_use_2 <- fsr_fits_2$res[,'est'];
    print(models[i])
    
    ## calculate there the values of the input function for each arm separately 
    b_use <- as.numeric(out_input_endp %>% dplyr::filter(Function == models[i]))[3:4]
    logit_arg <- b_use%*%rbind(c(1,1),c(0,1)) 
    pi_v <- exp(logit_arg)/(1 + exp(logit_arg))
    
    #print(paste0(" Cure values ",pi_v))
    opt_obj_1 <- try(optim(par=c(fsr_use_1[1:(length(fsr_use_1)) ]), #dummy staring parameter for the intercept and coefficient of the cure fraction  
                           ll.mix_f,
                           pi_ = v_cure_1[nt],
                           table_ = subset(event, ARMCDN == 1),
                           time_col = "AVAL",
                           cen_col = "CNSR",
                           rate_col = "rate_mod",
                           method = "Nelder-Mead", 
                           obj = fsr_fits_1,
                           hessian = TRUE
                           #lower=c(rep(-Inf,length(fsr_use)),-Inf,-Inf), upper = c(rep(Inf,length(fsr_use)),Inf,Inf)
    ) , silent = TRUE)
    #print(paste0("cure  ",v_cure_1[nt]))
    #print(opt_obj_1)
    
    opt_obj_2 <- try(optim(par=c(.1*fsr_use_2[1:(length(fsr_use_2)) ]), #dummy staring parameter for the intercept and coefficient of the cure fraction  
                           ll.mix_f,
                           pi_ = v_cure_2[nt],
                           table_ = subset(event, ARMCDN == 2),
                           time_col = "AVAL",
                           cen_col = "CNSR",
                           rate_col = "rate_mod",
                           method = "Nelder-Mead", 
                           obj = fsr_fits_2,
                           hessian = TRUE), 
                     #lower=c(rep(-Inf,length(fsr_use)),-Inf,-Inf), upper = c(rep(Inf,length(fsr_use)),Inf,Inf)) , 
                     silent = FALSE)
    #print(opt_obj_2)
    if ((class(opt_obj_1)[1] != "try-error") ) {
      #print (paste0(models[i], "  ",nt))
      parms_corr_1[nt,1:(length(mat_in_use$Param_Name)-2)] <- opt_obj_1$par
      
    }else{
      print(" cure 1 not converged")
      parms_corr_1[nt,] <- NA
    }
    
    if ((class(opt_obj_2)[1] != "try-error") ) {
      #print (paste0(models[i], "  ",nt))
      parms_corr_2[nt,1:(length(mat_in_use$Param_Name)-2)] <- opt_obj_2$par
    }else{
      print(" cure 2 not converged")
      parms_corr_2[nt,] <- NA
    }
  }
  
  
  #removal of the missing values 
  parms_corr_1 <- na.omit(parms_corr_1)
  parms_corr_2 <- na.omit(parms_corr_2)
  
  parm_names <- rownames(fsr_fits_1$res)
  print(parm_names)
  row_util <- length(parm_names) + 2
  
  parms_means_1 <- colMeans(parms_corr_1)
  parms_means_2 <- colMeans(parms_corr_2)
  
  
  miu <- models[i]
  parms_fit_frame_1[init_pos:(init_pos+row_util -1),"Function"] <- miu
  parms_fit_frame_1[init_pos:(init_pos+row_util -3),"Param_Name"] <- parm_names
  parms_fit_frame_1[(init_pos+row_util-2),"Param_Name"] <- "PI_1"
  parms_fit_frame_1[(init_pos+row_util -1),"Param_Name"] <- "PI_2"
  parms_fit_frame_1[init_pos:(init_pos+row_util-3),"Value_mix"] <-  parms_means_1[1:(length(mean_vals)-2)]
  parms_fit_frame_1[(init_pos+row_util-2),"Value_mix"] <- b_use[1]
  parms_fit_frame_1[(init_pos+row_util-1),"Value_mix"] <- b_use[2]
  parms_fit_frame_1[init_pos:(init_pos+row_util-1),4:(4+row_util - 1)] <-  cov(parms_corr_1);
  
  parms_fit_frame_2[init_pos:(init_pos+row_util -1),"Function"] <- miu
  parms_fit_frame_2[init_pos:(init_pos+row_util -3),"Param_Name"] <- parm_names
  parms_fit_frame_2[(init_pos+row_util-2),"Param_Name"] <- "PI_1"
  parms_fit_frame_2[(init_pos+row_util -1),"Param_Name"] <- "PI_2"
  parms_fit_frame_2[init_pos:(init_pos+row_util-3),"Value_mix"] <-  parms_means_2[1:(length(mean_vals)-2)]
  parms_fit_frame_2[(init_pos+row_util-2),"Value_mix"] <- b_use[1]
  parms_fit_frame_2[(init_pos+row_util-1),"Value_mix"] <- b_use[2]
  parms_fit_frame_2[init_pos:(init_pos+row_util-1),4:(4+row_util - 1)] <-  cov(parms_corr_2);
  
  out_frame_control[i,c(parm_names, "PI_1","PI_2")] <- parms_means_1[1:(length(mean_vals))]
  out_frame_control$Distribution[i] <- models[i]
  out_frame_control$STUDY.NAME[i] <- "GO29365"
  out_frame_control$POPULATION[i] <- "ITT"
  out_frame_control$TREATMENT.ARM[i] <- "Control"
  out_frame_control$Endpoint[i] <- end_point
  out_frame_control$Model[i] <- "Separate"
  out_frame_control$Data_cut <- data_cut_m
  out_frame_control$aic[i] <- 2*(length(parm_names) + 2) + 2*opt_obj_1$value
  out_frame_control$bic[i] <- sum(sf_pola$n)*log(length(parm_names) + 2) + 2*opt_obj_1$value
   
  
  out_frame_intervention[i,c(parm_names, "PI_1","PI_2")] <- parms_means_2[1:(length(mean_vals))]
  out_frame_intervention$Distribution[i] <- models[i]
  out_frame_intervention$STUDY.NAME[i] <- "GO29365"
  out_frame_intervention$POPULATION[i] <- "ITT"
  out_frame_intervention$TREATMENT.ARM[i] <- "Intervention"
  out_frame_intervention$Endpoint[i] <- end_point
  out_frame_intervention$Model[i] <- "Separate"
  out_frame_intervention$Data_cut <- data_cut_m
  out_frame_intervention$aic[i] <- 2*(length(parm_names) + 2) + 2*opt_obj_2$value
  out_frame_intervention$bic[i] <- sum(sf_pola$n)*log(length(parm_names) + 2) + 2*opt_obj_2$value
  

  l_mat <- length(mean_vals)
  for ( irow in 1:l_mat){
    for ( icol in 1:l_mat){
      cov_frame_intervention$STUDY.NAME[j0] <- "GO293650"
      cov_frame_intervention$POPULATION[j0] <- "ITT"
      cov_frame_intervention$TREATMENT.ARM[j0] <- "Intervention"
      cov_frame_intervention$Distribution[j0] <- models[i]
      cov_frame_intervention$Endpoint[j0] <- end_point
      cov_frame_intervention$Data_cut[j0] <- data_cut_m
      cov_frame_intervention$row_name[j0] <- (c(parm_names, "cure param 1","cure param 2"))[irow]
      cov_frame_intervention$column_name[j0] <- (c(parm_names, "cure param 1","cure param 2"))[icol]
      cov_frame_intervention$row_num[j0] <- irow
      cov_frame_intervention$col_num[j0] <- icol
      cov_frame_intervention$data_point[j0] <- (cov(parms_corr_2))[irow,icol]
      
      cov_frame_control$STUDY.NAME[j0] <- "GO293650"
      cov_frame_control$POPULATION[j0] <- "ITT"
      cov_frame_control$TREATMENT.ARM[j0] <- "Control"
      cov_frame_control$Distribution[j0] <- models[i]
      cov_frame_control$Endpoint[j0] <- end_point
      cov_frame_control$Data_cut[j0] <- data_cut_m
      cov_frame_control$row_name[j0] <- (c(parm_names, "cure param 1","cure param 2"))[irow]
      cov_frame_control$column_name[j0] <- (c(parm_names, "cure param 1","cure param 2"))[icol]
      cov_frame_control$row_num[j0] <- irow
      cov_frame_control$col_num[j0] <- icol
      cov_frame_control$data_point[j0] <- (cov(parms_corr_1))[irow,icol]
      
      
      j0 <- j0 + 1
    }
  }
  

  

  init_pos <- init_pos + row_util 
  ### create now plots for 1 and 2
  
  
  surv_mean <- haz_res_romu$surv_mean
  xty <- haz_res_romu$xty
  S_1_crude <- do.call(fsr_fits_1$dfns$p, args = c(as.list(parms_means_1[1:(length(fsr_use_1))]),
                                                   list(q = xty, lower.tail = FALSE)))
  
  S_2_crude <-  do.call(fsr_fits_1$dfns$p, args = c(as.list(parms_means_2[1:(length(fsr_use_2))]),
                                                    list(q = xty, lower.tail = FALSE)))
  
  
  S_1 <- surv_mean*((1-pi_v[1])*do.call(fsr_fits_1$dfns$p, args = c(as.list(parms_means_1[1:(length(fsr_use_1))]),
                                                                    list(q = xty, lower.tail = FALSE))) + pi_v[1])
  
  S_2 <- surv_mean*((1-pi_v[2])*do.call(fsr_fits_2$dfns$p, args = c(as.list(parms_means_2[1:(length(fsr_use_2))]),
                                                                    list(q = xty, lower.tail = FALSE))) + pi_v[2])
  
  d_frame <- rbind(d_frame, cbind(xty,S_1, S_2, models[i]))
  
  #################################
  ## a PLOOAAAT 
  ################################
  
  sf_pola1 <- survfit(Surv(as.numeric(AVAL)/12, CNSR)~1, data = subset(event, ARMCD == "RAN DL BR"))
  sf_pola2 <- survfit(Surv(as.numeric(AVAL)/12, CNSR)~1, data = subset(event, ARMCD == "RAN DL BR+POV"))
  if (censor_transpl_fl){
    pdf(paste0("reports/Plot_informed_",end_point,"_from_",end_point_info,"_",models[i],"_ad_hoc_cens_",data_cut, ".pdf"))
  }else{ 
    pdf(paste0("reports/Plot_informed_",end_point,"_from_",end_point_info,"_",models[i],"_",data_cut, ".pdf"))
  }
  plot(sf_pola2, lwd = 4, col = "blue", conf.int = F, xaxis = "years", xlim = c(0,7), xlab = "Time (years)", ylab = "Survival", cex.axis  =1.2, cex.lab = 1.3, main = paste0("Distribution = ", models[i]))
  lines(xty, S_1, col = 'red', lwd = 5, lty = 3)
  lines(xty, S_2, col = 'blue', lwd = 5, lty = 3)
    lines(sf_pola1, col = "red", lwd= 4, conf.int = F)
  legend("topright", c("POLA TRIAL -- BR + POLA", "POLA + BR cure","POLA TRIAL -- BR ", "BR cure","Romulus"), col = c("blue","blue","red","red", "purple"), lwd = c(5,5,5,5,5), lty = c(1,3,1,3,1))
  text(6,pi_v[1], paste0(100*signif(pi_v[1],3),"%"), cex = 1)
  text(6,pi_v[2], paste0(100*signif(pi_v[2],3),"%"), cex = 1)
  dev.off()
  
}


## output the informed uncertainties
if (censor_transpl_fl){
  write.csv(parms_fit_frame_1, file = paste0("reports/Parms_mixed_informed_",end_point,"_from_",end_point_info, "_uncertainty_control_ad_hoc_cens_", data_cut,".csv"), sep = "")
  write.csv(parms_fit_frame_2, file = paste0("reports/Parms_mixed_informed_",end_point,"_from_",end_point_info, "_uncertainty_intervention_ad_hoc_cens_", data_cut,".csv"), sep = "")

  write.csv(rbind(out_frame_control, out_frame_intervention), file = paste0("reports/Parms_mixed_informed_HE_friendly_", end_point,"_from_", end_point_info,"_uncertainty_ad_hoc_cens_",data_cut,".csv"))
  write.csv(rbind(cov_frame_control, cov_frame_intervention), file = paste0("reports/Variance_mixed_informed_HE_friendly_", end_point,"_from_", end_point_info,"_uncertainty_ad_hoc_cens_",data_cut,".csv"))
  write.csv(d_frame, file = paste0("reports/Surv_times_",end_point,"_from_", end_point_info,"_intervention_and_control_ad_hoc_cens_",data_cut,".csv"), sep = "")
}else{
  write.csv(parms_fit_frame_1, file = paste0("reports/Parms_mixed_informed_",end_point,"_from_",end_point_info, "_uncertainty_control_", data_cut,".csv"), sep = "")
  write.csv(parms_fit_frame_2, file = paste0("reports/Parms_mixed_informed_",end_point,"_from_",end_point_info, "_uncertainty_intervention_", data_cut,".csv"), sep = "")
  
  write.csv(rbind(out_frame_control, out_frame_intervention), file = paste0("reports/Parms_mixed_informed_HE_friendly_", end_point,"_from_", end_point_info,"_uncertainty_",data_cut,".csv"))
  write.csv(rbind(cov_frame_control, cov_frame_intervention), file = paste0("reports/Variance_mixed_informed_HE_friendly_", end_point,"_from_", end_point_info,"_uncertainty_",data_cut,".csv"))
  write.csv(d_frame, file = paste0("reports/Surv_times_",end_point,"_from_", end_point_info,"_intervention_and_control_",data_cut,".csv"), sep = "")
  
}