setwd("/gpfs/gibbs/project/zhao/qh63/R/summary_statistics")

# Initialize vectors to store means and standard errors for multiple MR methods and PRS results
mr_sisVIVE_mean <- mr_sisVIVE_se <- c()
mr_2SLS_mean <- mr_2SLS_se <- c()
mr_ivw_mean <- mr_ivw_se <- c()
mr_median_mean <- mr_median_se <- c()
mr_mode_mean <- mr_mode_se <- c()
mr_mix_mean <- mr_mix_se <- c()
mr_egger_mean <- mr_egger_se <- c()
mr_raps_mean <- mr_raps_se <- c()
PRS_PT_mean <- PRS_PT_se <- c()
PRS_Oracle_mean <- PRS_Oracle_se <- c()
PRS_ldpred_mean <- PRS_ldpred_se <- c()



# Loop through iterations to load data and calculate means and standard errors for each method
# Although a loop could simplify the code, I list methods individually here to facilitate debugging, allowing quick adjustments if an issue arises with any specific method.
for(count in 1:5){
  
  mr_sisVIVE_mean <- cbind(mr_sisVIVE_mean, lapply(Sys.glob(paste0("mr_sisVIVE_",count,".txt")),read.table)[[1]]$beta[1:100])
  mr_sisVIVE_se <- cbind(mr_sisVIVE_se, lapply(Sys.glob(paste0("mr_sisVIVE_",count,".txt")),read.table)[[1]]$se[1:100])
  mr_2SLS_mean <- cbind(mr_2SLS_mean, lapply(Sys.glob(paste0("mr_2SLS_",count,".txt")),read.table)[[1]]$beta)
  mr_2SLS_se <- cbind(mr_2SLS_se, lapply(Sys.glob(paste0("mr_2SLS_",count,".txt")),read.table)[[1]]$se)
  mr_ivw_mean <- cbind(mr_ivw_mean, lapply(Sys.glob(paste0("mr_ivw_",count,".txt")),read.table)[[1]]$beta)
  mr_ivw_se <- cbind(mr_ivw_se, lapply(Sys.glob(paste0("mr_ivw_",count,".txt")),read.table)[[1]]$se)
  mr_median_mean <- cbind(mr_median_mean, lapply(Sys.glob(paste0("mr_median_",count,".txt")),read.table)[[1]]$beta)
  mr_median_se <- cbind(mr_median_se, lapply(Sys.glob(paste0("mr_median_",count,".txt")),read.table)[[1]]$se)
  mr_mode_mean <- cbind(mr_mode_mean, lapply(Sys.glob(paste0("mr_mode_",count,".txt")),read.table)[[1]]$beta)
  mr_mode_se <- cbind(mr_mode_se, lapply(Sys.glob(paste0("mr_mode_",count,".txt")),read.table)[[1]]$se)
  mr_mix_mean <- cbind(mr_mix_mean, lapply(Sys.glob(paste0("mr_mix_",count,".txt")),read.table)[[1]]$beta)
  mr_mix_se <- cbind(mr_mix_se, lapply(Sys.glob(paste0("mr_mix_",count,".txt")),read.table)[[1]]$se)
  mr_egger_mean <- cbind(mr_egger_mean, lapply(Sys.glob(paste0("mr_egger_",count,".txt")),read.table)[[1]]$beta)
  mr_egger_se <- cbind(mr_egger_se, lapply(Sys.glob(paste0("mr_egger_",count,".txt")),read.table)[[1]]$se)
  mr_raps_mean <- cbind(mr_raps_mean, lapply(Sys.glob(paste0("mr_raps_",count,".txt")),read.table)[[1]]$beta)
  mr_raps_se <- cbind(mr_raps_se, lapply(Sys.glob(paste0("mr_raps_",count,".txt")),read.table)[[1]]$se)
  PRS_PT_mean <- cbind(PRS_PT_mean, lapply(Sys.glob(paste0("PRS_PT_",count,".txt")),read.table)[[1]]$beta)
  PRS_PT_se <- cbind(PRS_PT_se, lapply(Sys.glob(paste0("PRS_PT_",count,".txt")),read.table)[[1]]$se)
  PRS_Oracle_mean <- cbind(PRS_Oracle_mean, lapply(Sys.glob(paste0("PRS_Oracle_",count,".txt")),read.table)[[1]]$beta)
  PRS_Oracle_se <- cbind(PRS_Oracle_se, lapply(Sys.glob(paste0("PRS_Oracle_",count,".txt")),read.table)[[1]]$se)
  PRS_ldpred_mean <- cbind(PRS_ldpred_mean, lapply(Sys.glob(paste0("PRS_ldpred_",count,".txt")),read.table)[[1]]$beta)
  PRS_ldpred_se <- cbind(PRS_ldpred_se, lapply(Sys.glob(paste0("PRS_ldpred_",count,".txt")),read.table)[[1]]$se)
  
  
}

# Function to calculate mean and standard error excluding outliers
cal_without_outlier <- function(mr_mean_column){
  Q1 <- quantile(mr_mean_column, 0.25)
  Q3 <- quantile(mr_mean_column, 0.75)
  IQR_value <- IQR(mr_mean_column)
  
  # Calculate interquartile range and define bounds for outliers
  lower_bound <- Q1 - 2 * IQR_value
  upper_bound <- Q3 + 2 * IQR_value
  
  # Filter out outliers and compute summary statistics
  filtered_data <- mr_mean_column[mr_mean_column >= lower_bound & mr_mean_column <= upper_bound]
  mean_value <- mean(filtered_data)
  standard_error <- sd(filtered_data) 
  return(c(mean_value,standard_error))
}


library(ggplot2)

# Prepare a dataframe for plotting beta estimates without outliers
circle_num <- 5
colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=T)
df <- data.frame(pi_x = rep(c(0.0005,0.0016,0.005,0.016,0.05),11),beta = c(apply(mr_sisVIVE_mean, 2, cal_without_outlier)[1,],apply(mr_2SLS_mean, 2, cal_without_outlier)[1,],apply(mr_ivw_mean, 2, cal_without_outlier)[1,],
                                                                           apply(mr_median_mean, 2, cal_without_outlier)[1,],apply(mr_mode_mean, 2, cal_without_outlier)[1,],apply(mr_mix_mean, 2, cal_without_outlier)[1,],
                                                                           apply(mr_egger_mean, 2, cal_without_outlier)[1,],apply(mr_raps_mean, 2, cal_without_outlier)[1,],apply(PRS_PT_mean, 2, cal_without_outlier)[1,],
                                                                           apply(PRS_Oracle_mean, 2, cal_without_outlier)[1,],apply(PRS_Oracle_mean, 2, cal_without_outlier)[1,]),
                 group = c(rep("mr_sisVIVE",circle_num),rep("mr_2SLS",circle_num),rep("mr_ivw",circle_num),rep("mr_median",circle_num),rep("mr_mode",circle_num),
                           rep("mr_mix",circle_num),rep("mr_egger",circle_num),rep("mr_raps",circle_num),rep("PRS_PT",circle_num),rep("PRS_Oracle",circle_num),rep("PRS_LDPred",circle_num)),
                 sd = c(apply(mr_sisVIVE_mean, 2, cal_without_outlier)[2,],apply(mr_2SLS_mean, 2, cal_without_outlier)[2,],apply(mr_ivw_mean, 2, cal_without_outlier)[2,],
                        apply(mr_median_mean, 2, cal_without_outlier)[2,],apply(mr_mode_mean, 2, cal_without_outlier)[2,],apply(mr_mix_mean, 2, cal_without_outlier)[2,],
                        apply(mr_egger_mean, 2, cal_without_outlier)[2,],apply(mr_raps_mean, 2, cal_without_outlier)[2,],apply(PRS_PT_mean, 2, cal_without_outlier)[2,],
                        apply(PRS_Oracle_mean, 2, cal_without_outlier)[2,],apply(PRS_Oracle_mean, 2, cal_without_outlier)[2,]))



# Generate a plot of beta estimates versus log(pi_x) using ggplot
ggplot(df, aes(x=log(pi_x), y=beta, group=group, color=group)) + 
  geom_errorbar(aes(ymin=beta-1.96*sd, ymax=beta+1.96*sd), width=1,position=position_dodge(0.2)) +
  geom_line(position=position_dodge(0.2)) + geom_point(position=position_dodge(0.2))+
  scale_x_continuous(name="log(pi_x)", breaks=seq(log(0.0005),log(0.05),length.out = 5),labels = c("log(0.0005)","log(0.0016)","log(0.005)","log(0.016)","log(0.05)"))+
  labs(title = "Plot of beta_hat versus log(pi_x)", x = "log(pi_x)", y = "beta_hat") + theme(plot.title = element_text(hjust = 0.5))

############################

#Function for Calculating Type-1 error
t1e_without_outlier <- function(experi_mean_100,experi_sd_100){
  
  Q1 <- quantile(experi_mean_100, 0.25)
  Q3 <- quantile(experi_mean_100, 0.75)
  IQR_value <- IQR(experi_mean_100)
  
  lower_bound <- Q1 - 2 * IQR_value
  upper_bound <- Q3 + 2 * IQR_value
  
  filtered_mean <- experi_mean_100[experi_mean_100 >= lower_bound & experi_mean_100 <= upper_bound]
  filtered_se <- experi_sd_100[experi_mean_100 >= lower_bound & experi_mean_100 <= upper_bound]
  
  pval <- 2*pnorm(abs(filtered_mean/filtered_se),lower.tail = FALSE)
  temp_t1e <- sum(pval<=0.05)/length(pval)
  return(temp_t1e)
}



cal_type_1_error <- function(mean_matrix, se_matrix){
  t1e <- c()
  for(i in 1:ncol(mean_matrix)){
    t1e[i] <- t1e_without_outlier(mean_matrix[,i],se_matrix[,i])
  }
  return(t1e)
}


# Prepare a dataframe for plotting Type-1 error
df_2 <- data.frame(pi_x = rep(c(0.0005,0.0016,0.005,0.016,0.05),11),
                   type_1_error = c(cal_type_1_error(mr_sisVIVE_mean,mr_sisVIVE_se),cal_type_1_error(mr_2SLS_mean,mr_2SLS_se),cal_type_1_error(mr_ivw_mean,mr_ivw_se),
                                    cal_type_1_error(mr_median_mean,mr_median_se),cal_type_1_error(mr_mode_mean,mr_mode_se),cal_type_1_error(mr_mix_mean,mr_mix_se),
                                    cal_type_1_error(mr_egger_mean,mr_egger_se),cal_type_1_error(mr_raps_mean,mr_raps_se),cal_type_1_error(PRS_PT_mean,PRS_PT_se),
                                    cal_type_1_error(PRS_Oracle_mean,PRS_Oracle_se),cal_type_1_error(PRS_ldpred_mean,PRS_ldpred_se)),
                   group = c(rep("mr_sisVIVE",circle_num),rep("mr_2SLS",circle_num),rep("mr_ivw",circle_num),rep("mr_median",circle_num),rep("mr_mode",circle_num),
                           rep("mr_mix",circle_num),rep("mr_egger",circle_num),rep("mr_raps",circle_num),rep("PRS_PT",circle_num),rep("PRS_Oracle",circle_num),rep("PRS_LDPred",circle_num)))


# Generate a plot of Type-1 error versus log(pi_x) using ggplot
ggplot(df_2, aes(x=log(pi_x), y=type_1_error, group=group, color=group)) + 
  geom_line(position=position_dodge(0.002)) + geom_point(position=position_dodge(0.002))+
  scale_x_continuous(name="log(pi_x)", breaks=seq(log(0.0005),log(0.05),length.out = 5),labels = c("log(0.0005)","log(0.0016)","log(0.005)","log(0.016)","log(0.05)"))+
  labs(title = "Plot of type_1_error", x = "log(pi_x)", y = "type_1_error") + theme(plot.title = element_text(hjust = 0.5))

