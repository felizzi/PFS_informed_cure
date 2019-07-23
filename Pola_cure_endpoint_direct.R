# cures 
require(dplyr)
require(flexsurv)
require(stem2r)
require(optimx)
require(MASS)
require(optimr)

source("haz_countries.R")

data_cut <- 'i29365h'
#data_cut <- "ash_2019"
data_file <- 'ate.sas7bdat'
#data_cut_m <- "2018 October"
data_cut_m <- "2019 March"

end_point <- "PFSINV"

censor_transpl_fl <- FALSE

SAICE::initialize_connection(entimice_env = "PROD") # WebUI prod version
event_enti <- SAICE::read_entimice("root/clinical_studies/RO5541077/CDPT7898/GO29365/data_analysis/ema_questions_90dsu/qa/outdata_vad/ateupdt.sas7bdat")
event_ash <- SAICE::read_entimice("root/clinical_studies/RO5541077/CDPT7898/GO29365/data_analysis/ash_2019/qa/outdata_vad/ate.sas7bdat")
data_loc <- paste('/opt/BIOSTAT/qa/cdpt7898', data_cut, 'libraries',  data_file, sep = '/')

event <- read_bce(data_loc) %>% dplyr::filter(PARAMCD == end_point) %>% 
  dplyr::select(AVAL, CNSR, ARMCD, AGE, SEX, COUNTRY, USUBJID ) %>% dplyr::mutate(CNSR = (1+CNSR)%%2) %>% filter(ARMCD %in% c("RAN DL BR", "RAN DL BR+POV"))

if (end_point == "PFSIRC" & data_cut != "ash_2019"){
  print(" entered IRC PFS sel")
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


event$YEAR <- 2016
haz_res_romu <- hazard_time(table_ = as.data.frame(event), evttme = "AVAL", sex = "SEX" , 
                            age = "AGE", year = "YEAR", country_trial =  "COUNTRY", country_output = NULL)

models <- c("exp",       "weibullPH", "lnorm",     "gengamma", "llogis",  "gompertz" )



event$rate_mod <- haz_res_romu$rate_vec

##### likelihood function 

parms_fit_frame_endp <- data.frame("Function" = character(0), "Param_Name" = character(0), "Value_mix" = numeric(0), "Cov_mix_p1" = numeric(0), "Cov_mix_p2" = numeric(0), "Cov_mix_p3" = numeric(0), "Cov_mix_PI_1" = numeric(0), "Cov_mix_PI_2" = numeric(0), stringsAsFactors=FALSE)
est_AIC_BIC <- data.frame("Function" = character(6), "AIC" = numeric(6), "BIC" = numeric(6), stringsAsFactors = F)

Nfs <- length(models)

out_frame_comb <- data.frame("STUDY NAME" = character(Nfr), "POPULATION" = character(Nfr),	"TREATMENT ARM" = character(Nfr), 
                                     "Endpoint" = character(Nfr), "Endpoint_def" = character(Nfr), "Model" = character(Nfr), "Data_cut" = character(Nfr),
                                     "Distribution" = character(Nfr),	
                                     "rate" = numeric(Nfr), 	"shape" = numeric(Nfr),	"scale" = numeric(Nfr),	"meanlog" = numeric(Nfr), "sdlog" = numeric(Nfr),
                                     
                                     "mu" = numeric(Nfr), "sigma" = numeric(Nfr), "Q" = numeric(Nfr), "PI_1" = numeric(Nfr), "PI_2" = numeric(Nfr),        "aic" = numeric(Nfr),	"bic" = numeric(Nfr), 
                                     "Info" = character(Nfr),
                                     "Source" = character(Nfr), stringsAsFactors = F)

Nfr_cov <- 16*4 + 25 + 9 


cov_frame_comb <- data.frame("STUDY NAME" = character(Nfr_cov), "POPULATION" = character(Nfr_cov),	"TREATMENT ARM" = character(Nfr_cov), 
                                     "Endpoint" = character(Nfr_cov), "Endpoint_def" = character(Nfr_cov), "Data_cut" = character(Nfr_cov), "Hazard" = character(Nfr_cov),
                                     "Distribution" = character(Nfr_cov),	"row_name" = numeric(Nfr_cov),
                                     "column_name" = numeric(Nfr_cov), 	"row_num" = numeric(Nfr_cov),	"col_num" = numeric(Nfr_cov),	
                                     "data_point" = numeric(Nfr_cov), stringsAsFactors = F)



ll.mix.b <- function(table, parms, time_col, cen_col, rate_col, trt_col, obj){  
  time_v <- as.numeric(table[,time_col])/12;
  cen_v <- as.numeric(table[,cen_col]);
  rate_v <- as.numeric(table[,rate_col]);
  trt_v <- as.numeric(table[,trt_col]) - 1; # the treatment value is coded 0, 1
  lp <- length(parms)
  b <- parms[(lp-1):lp];
  
  #construct the vector PI; 
  #shape <- parms[1];  scale <- parms[2]; pi_ <- parms[3];  
  f_t <- do.call(obj$dfns$d, args = c(as.list(parms[1:(length(parms)-2)]),
                                      list (x = time_v, log = FALSE))
  )
  S_t <- do.call(obj$dfns$p, args = c(as.list(parms[1:(length(parms)-2)]),
                                      list(q = time_v, lower.tail = FALSE))
  )
  logit_arg <- b%*%rbind(1,trt_v) 
  pi_v <- exp(logit_arg)/(1 + exp(logit_arg))
  h_ <- rate_v + ((1-pi_v)*f_t)/(pi_v + (1-pi_v)*S_t);  
  s_ <- pi_v + (1 - pi_v)*S_t;  
  ret_val <- sum( cen_v*log(h_) + log(s_));  
  return(-ret_val)
  
}              



event <- event %>% dplyr::mutate(ARMCDN = ifelse(ARMCD == "RAN DL BR",1,2))
event<- as.data.frame(event)



rnd <- 2*runif(2) -1
### loop over many things 

Nfr <- length(models)



out_endp <- data.frame(
  "Distribution" = character(Nfr),	"intercept" = numeric(Nfr),
  "scale" = numeric(Nfr), 	"shape" = numeric(Nfr),	"theta" = numeric(Nfr),
  "b1" = numeric(Nfr), "b2" = numeric(Nfr)
  , stringsAsFactors = F)


init_pos <- 1
j0 <- 1

### loop over the different distributions 
for (i in 1:length(models)){
#for (i in 5:5){
  fsr_fits_ <- flexsurvreg(Surv(as.numeric(AVAL), CNSR)~ARMCD, data = event, dist = models[i])
  fsr_use <- fsr_fits_$res[,'est'];
  print(models[i])
  attempts <- 1
  opt_obj_ <- try(optimx(par=c(fsr_use[1:(length(fsr_use) -1)], rnd), #dummy staring parameter for the intercept and coefficient of the cure fraction  
                         ll.mix.b,                      
                         table = event,
                         time_col = "AVAL",
                         cen_col = "CNSR",
                         rate_col = "rate_mod",
                         trt_col = "ARMCDN",
                         method = "Nelder-Mead", 
                         obj = fsr_fits_,
                         hessian = TRUE), 
                  #lower=c(rep(-Inf,length(fsr_use)),-Inf,-Inf), upper = c(rep(Inf,length(fsr_use)),Inf,Inf)) , 
                  silent = TRUE)
  
  
  opt_obj_n <- try(optim(par=as.numeric(opt_obj_)[1:(length(fsr_use) +1 )], #dummy staring parameter for the intercept and coefficient of the cure fraction  
                          ll.mix.b,                      
                          table = event,
                          time_col = "AVAL",
                          cen_col = "CNSR",
                          rate_col = "rate_mod",
                          trt_col = "ARMCDN",
                          method = "SANN", 
                          obj = fsr_fits_,
                          hessian = TRUE), 
                   #lower=c(rep(-Inf,length(fsr_use)),-Inf,-Inf), upper = c(rep(Inf,length(fsr_use)),Inf,Inf)) , 
                   silent = TRUE)
  
  

  
  invHuse <- ginv(opt_obj_n$hessian)
  #invHuse_2 <- ginv(opt_obj_2$hessian)
  parm_names <- names(fsr_fits_$res[,'est'])
  print(parm_names)
  parm_names <- parm_names[1:(length(parm_names) - 1)]
  row_util <- length(parm_names) + 2 ## number of parameter + 2 cure coefficients 
  miu <- models[i]
  parms_fit_frame_endp[init_pos:(init_pos+row_util -1),"Function"] <- miu
  parms_fit_frame_endp[init_pos:(init_pos+row_util -3),"Param_Name"] <- parm_names
  parms_fit_frame_endp[(init_pos+row_util-2),"Param_Name"] <- "PI_1"
  parms_fit_frame_endp[(init_pos+row_util -1),"Param_Name"] <- "PI_2"
  parms_fit_frame_endp[init_pos:(init_pos+row_util-3),"Value_mix"] <- opt_obj_n$par[1:(length(parm_names))]
  parms_fit_frame_endp[(init_pos+row_util-2),"Value_mix"] <- opt_obj_n$par[(length(parm_names)+1):(length(parm_names)+1)]
  parms_fit_frame_endp[(init_pos+row_util-1),"Value_mix"] <- opt_obj_n$par[(length(parm_names)+2):(length(parm_names)+2)]
  parms_fit_frame_endp[init_pos:(init_pos+row_util-1),4:(4+row_util - 1)] <- invHuse
  init_pos <- init_pos + row_util
  
  # modify invHuse depending on the signular values. Reduce the very large singular values, as we do for the informed approach 
  
  sing_dec <- svd(invHuse)
  
  if (min(sing_dec$d) < 1000*.Machine$double.eps){
    print("entered if small")
    sing_dec$d[which(sing_dec$d<.Machine$double.eps)] <- 1100*.Machine$double.eps
    invHuse <- sing_dec$u%*%diag(sing_dec$d)%*%t(sing_dec$v)
  }
  if (max(sing_dec$d) > 10){ 
    print("entered if large")
    sing_dec$d[which(sing_dec$d>10)] <- 1 #this needs to be clinically justified 
    invHuse <- sing_dec$u%*%diag(sing_dec$d)%*%t(sing_dec$v)
  }

  out_frame_comb[i,c(parm_names, "PI_1","PI_2")] <-opt_obj_n$par[1:(length(parm_names) + 2 )]
  out_frame_comb$Distribution[i] <- models[i]
  out_frame_comb$STUDY.NAME[i] <- "GO29365"
  out_frame_comb$POPULATION[i] <- "ITT"
  out_frame_comb$TREATMENT.ARM[i] <- "Combined"
  out_frame_comb$Endpoint[i] <- end_point
  out_frame_comb$Model[i] <- "Combined"
  out_frame_comb$Data_cut[i] <- data_cut_m
  out_frame_comb$aic[i] <- 2*(length(parm_names) + 2) + 2*opt_obj_n$value
  out_frame_comb$bic[i] <- sum(sf_pola$n)*log(length(parm_names) + 2) + 2*opt_obj_n$value
  
  
  l_mat <- length(parm_names) + 2 
  for ( irow in 1:l_mat){
    for ( icol in 1:l_mat){
      cov_frame_comb$STUDY.NAME[j0] <- "GO293650"
      cov_frame_comb$POPULATION[j0] <- "ITT"
      cov_frame_comb$TREATMENT.ARM[j0] <- "Combined"
      cov_frame_comb$Distribution[j0] <- models[i]
      cov_frame_comb$Endpoint[j0] <- end_point
      cov_frame_comb$Data_cut[j0] <- data_cut_m
      cov_frame_comb$row_name[j0] <- (c(parm_names, "cure param 1","cure param 2"))[irow]
      cov_frame_comb$column_name[j0] <- (c(parm_names, "cure param 1","cure param 2"))[icol]
      cov_frame_comb$row_num[j0] <- irow
      cov_frame_comb$col_num[j0] <- icol
      cov_frame_comb$data_point[j0] <- invHuse[irow,icol]
      
      
      
      j0 <- j0 + 1
    }
  }
  
  
  if (class(opt_obj_)[1] != "try-error"){
    out_endp[i,"Distribution"] <- models[i]
    trt_v <- unique(event$ARMCDN -1)
    b <- as.numeric(opt_obj_[(length(fsr_use) ):(length(fsr_use) + 1)] )
    logit_arg <- b%*%rbind(1,trt_v) 
    pi_v <- exp(logit_arg)/(1 + exp(logit_arg))
    
    surv_mean <- haz_res_romu$surv_mean
    
    xty <- haz_res_romu$xty
    
    
    S_1 <- surv_mean*((1-pi_v[1])*do.call(fsr_fits_$dfns$p, args = c(as.list(opt_obj_[1:(length(fsr_use) -1)]),
                                                                     list(q = xty, lower.tail = FALSE))) + pi_v[1])
    
    S_2 <- surv_mean*((1-pi_v[2])*do.call(fsr_fits_$dfns$p, args = c(as.list(opt_obj_[1:(length(fsr_use) -1)]),
                                                                     list(q = xty, lower.tail = FALSE))) + pi_v[2])
    print("trt+v")
    
    trt_v <- unique(event$ARMCDN -1)
    
    b <- as.numeric(opt_obj_[(length(fsr_use) ):(length(fsr_use) + 1)] )
    logit_arg <- b%*%rbind(1,trt_v) 
    pi_v <- exp(logit_arg)/(1 + exp(logit_arg))
    
    surv_mean <- haz_res_romu$surv_mean
    
    xty <- haz_res_romu$xty
    
    
    S_1 <- surv_mean*((1-pi_v[1])*do.call(fsr_fits_$dfns$p, args = c(as.list(opt_obj_[1:(length(fsr_use) -1)]),
                                                                     list(q = xty, lower.tail = FALSE))) + pi_v[1])
    
    S_2 <- surv_mean*((1-pi_v[2])*do.call(fsr_fits_$dfns$p, args = c(as.list(opt_obj_[1:(length(fsr_use) -1)]),
                                                                     list(q = xty, lower.tail = FALSE))) + pi_v[2])
    
    
    sf_pola <- survfit(Surv(as.numeric(AVAL)/12, CNSR)~ARMCD, data = event)
    sf_pola1 <- survfit(Surv(as.numeric(AVAL)/12, CNSR)~1, data = subset(event, ARMCD == "RAN DL BR"))
    sf_pola2 <- survfit(Surv(as.numeric(AVAL)/12, CNSR)~1, data = subset(event, ARMCD == "RAN DL BR+POV"))
    sf_pola_2_m <- survfit(Surv(as.numeric(AVAL), CNSR)~1, data = subset(event, ARMCD == "RAN DL BR+POV"))
    if (censor_transpl_fl){
      pdf(paste0("reports/Plot_direct_",end_point,"_",models[i],"_ad_hoc_cens_",data_cut,".pdf"))
    }else{
      pdf(paste0("reports/Plot_direct_",end_point,"_",models[i],"_",data_cut,".pdf"))
    }
    plot(sf_pola2, lwd = 4, col = "blue", conf.int = F,  xlim = c(0,7),  xlab = "Time (years)", ylab = "Survival", cex.axis  =1.2, cex.lab = 1.3, main = paste0("Distribution = ", models[i]))
    lines(xty, S_1, col = 'red', lwd = 5, lty = 3)
    lines(xty, S_2, col = 'blue', lwd = 5, lty = 3)
    lines(sf_pola1, col = "red", lwd= 4, conf.int = F)
    legend("topright", c("POLA TRIAL -- BR + POLA", "POLA + BR cure","POLA TRIAL -- BR ", "BR cure"), col = c("blue","blue","red","red", "purple"), lwd = c(5,5,5,5,5), lty = c(1,3,1,3,1))
    text(6,pi_v[1], paste0(100*signif(pi_v[1],3),"%"), cex = 1)
    text(6,pi_v[2], paste0(100*signif(pi_v[2],3),"%"), cex = 1)
    dev.off()
    
    ###############################################################################################################
    ## contruct here the hazards fsr_fd
    ###############################################################################################################
    Sb_1 <- do.call(fsr_fits_$dfns$p, args = c(as.list(opt_obj_[1:(length(fsr_use) -1)]),
                                               list(q = xty, lower.tail = FALSE))) 
    
    fb_1 <- do.call(fsr_fits_$dfns$d, args = c(as.list(opt_obj_[1:(length(fsr_use) -1)]),
                                               list(x = xty, log = FALSE))) 
    
    h_1 <- ((1 - pi_v[1]) * fb_1)/(pi_v[1] + (1 - pi_v[1])*Sb_1 )
    
    
    Sb_2 <- do.call(fsr_fits_$dfns$p, args = c(as.list(opt_obj_[1:(length(fsr_use) -1)]),
                                               list(q = xty, lower.tail = FALSE))) 
    
    fb_2 <- do.call(fsr_fits_$dfns$d, args = c(as.list(opt_obj_[1:(length(fsr_use) -1)]),
                                               list(x = xty, log = FALSE))) 
    
    h_2 <- ((1 - pi_v[2]) * fb_2)/(pi_v[2] + (1 - pi_v[2])*Sb_2 )
    
    #######################################################################
#    x_der <- haz_res_romu$xty[1:(length(haz_res_romu$xty) - 1)]
#    h_star <- -diff(log(haz_res_romu$surv_mean)) / diff(haz_res_romu$xty[1:(length(haz_res_romu$xty) - 1)])
    
    #################################
    ## a PLOOAAAT 
    ################################
    
    #sf_pola <- survfit(Surv(as.numeric(AVAL)/12, CNSR)~ARMCD, data = event)
    #sf_pola1 <- survfit(Surv(as.numeric(AVAL)/12, CNSR)~1, data = subset(event, ARMCD == "RAN DL BR"))
    #sf_pola2 <- survfit(Surv(as.numeric(AVAL)/12, CNSR)~1, data = subset(event, ARMCD == "RAN DL BR+POV"))
#    sf_pola1_fsr <- flexsurvreg(Surv(as.numeric(AVAL)/12, CNSR)~1, data = subset(event, ARMCD == "RAN DL BR"), dist = "gengamma")
#    sf_pola2_fsr <- flexsurvreg(Surv(as.numeric(AVAL)/12, CNSR)~1, data = subset(event, ARMCD == "RAN DL BR+POV"), dist = "gengamma")
    
#    plot(sf_pola1_fsr, type = "hazard", ylim = c(0,2))
    
#    pdf(paste0("reports/Hazard_evaluation_plot_",endpoint,"_",models[i],"_",data_cut, ".pdf"))
#    par(mfrow=c(2,1))
#    plot(sf_pola1_fsr, type = "hazard", est = F, lwd.ci = 0, lwd = 4, col = "blue", conf.int = F, xaxis = "years", xlim = c(0,7), ylim = c(0,3.5), xlab = "Time (years)", ylab = "Survival", cex.axis  =1.2, cex.lab = 1.3)
#    lines(xty, h_1, col = 'red', lwd = 5, lty = 1)
#    lines(x_der, h_star, col = "black")
    
#    plot(sf_pola2_fsr, type = "hazard", est = F, lwd.ci = 0, lwd = 4, col = "red", conf.int = F, xaxis = "years", xlim = c(0,7), ylim = c(0,3.5), xlab = "Time (years)", ylab = "Survival", cex.axis  =1.2, cex.lab = 1.3)
#    lines(xty, h_2, col = 'red', lwd = 5, lty = 1)
#    lines(x_der, h_star, col = "black")
    
 
    
#    dev.off()
    
    
    
    
    est_AIC_BIC[i,"Function"] <- models[i]
    est_AIC_BIC[i,"AIC"] <- 2*(length(parm_names) + 2) + 2*opt_obj_n$value
    est_AIC_BIC[i,"BIC"] <- sum(sf_pola$n)*log(length(parm_names) + 2) + 2*opt_obj_n$value
    
  }
  
}

# export parms_fit_frame_endp 

if (censor_transpl_fl){
  write.table(parms_fit_frame_endp, file = paste0("reports/Parms_mixed_",end_point,"_intervention_and_control_ad_hoc_cens_", data_cut,".csv", sep = ""), append = F, row.names = F, sep = ',')
  write.table(data_loc, file = paste0("reports/Parms_mixed_",end_point,"_intervention_and_control_ad_hoc_cens_", data_cut,".csv",sep = ""), append = T, row.names = F)
  #write.csv(est_AIC_BIC, file = paste0("reports/AIB_BIC_pfsinv_",data_cut,".csv",sep = ""))
  write.csv(out_frame_comb, file = paste0("reports/Parms_mixed_direct_HE_friendly_", end_point,"_uncertainty_ad_hoc_cens_",data_cut,".csv"))
  write.csv(cov_frame_comb, file = paste0("reports/Variance_mixed_direct_HE_friendly_", end_point,"_uncertainty_ad_hoc_cens_",data_cut,".csv"))
}else{
  write.table(parms_fit_frame_endp, file = paste0("reports/Parms_mixed_",end_point,"_intervention_and_control_", data_cut,".csv", sep = ""), append = F, row.names = F, sep = ',')
  write.table(data_loc, file = paste0("reports/Parms_mixed_",end_point,"_intervention_and_control_", data_cut,".csv",sep = ""), append = T, row.names = F)
  #write.csv(est_AIC_BIC, file = paste0("reports/AIB_BIC_pfsinv_",data_cut,".csv",sep = ""))
  write.csv(out_frame_comb, file = paste0("reports/Parms_mixed_direct_HE_friendly_", end_point,"_uncertainty_",data_cut,".csv"))
  write.csv(cov_frame_comb, file = paste0("reports/Variance_mixed_direct_HE_friendly_", end_point,"_uncertainty_",data_cut,".csv"))
}
