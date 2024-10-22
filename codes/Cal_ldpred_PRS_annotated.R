#This script uses the output from ldpred_1.sh and calculates the estimated Mendelian randomization effect when polygenic risk scores (PRS) are computed using the LDPred method.

# Initialize vectors to store LDpred scores and other statistical values
library(TwoSampleMR)
ldpred_score_1 <- c()
ldpred_score_2 <- c()
test_p <- c()
each_count_b_exp <- c()  # Beta estimates for the exposure (X)
each_count_se_exp <- c() # Standard errors for the exposure estimates
each_count_b_out <- c()  # Beta estimates for the outcome (Y)
each_count_se_out <- c() # Standard errors for the outcome estimates

count <- 1 
cycle_num <- 100 

file_name <- paste0("/gpfs/gibbs/project/zhao/qh63/R/PLink_data/output_",count,"/output")
setwd(file_name)

# Main loop for processing data files and computing LDpred scores
for(cycle_num in 1:cycle_num){
  
  datafiles_1 <- lapply(Sys.glob(paste0("output_score_1_",cycle_num,"_*.txt")),read.csv)
  datafiles_2 <- lapply(Sys.glob(paste0("output_score_2_",cycle_num,"_*.txt")),read.csv)
  X_i_1 <- read.table(paste0("/gpfs/gibbs/project/zhao/qh63/R/summary_statistics/X_i/X_i_1_",count,"_",cycle_num,".txt"))
  Y_i_2 <- read.table(paste0("/gpfs/gibbs/project/zhao/qh63/R/summary_statistics/Y_i/Y_i_2_",count,"_",cycle_num,".txt"))
  
  # Calculate correlations between exposure data (X_i_1) and LDpred scores (second column of datafiles)
  cor_2 <- c()
  for(i in 1:length(datafiles_1)){
    cor_2 <- append(cor_2,cor(X_i_1,as.data.frame(datafiles_1[i])[,2])^2)
  }
  
  # Identify the index of the dataset with the minimum correlation
  each_cycle_ldpred_score_index <- which(cor_2 == min(cor_2))[1]
  
  # Extract the corresponding LDpred scores for the selected dataset
  each_cycle_ldpred_score_1 <- as.data.frame(datafiles_1[each_cycle_ldpred_score_index])[2]
  ldpred_score_1 <- append(ldpred_score_1,each_cycle_ldpred_score_1)
  
  ##############################################
  
  # Perform linear regression of outcome (Y) on LDpred score and store coefficients and standard errors
  ldpred_lm_X_result <- lm(unlist(X_i_1) ~ unlist(each_cycle_ldpred_score_1))
  each_count_b_exp <- append(each_count_b_exp,summary(ldpred_lm_X_result)$coefficients[,1][2])
  each_count_se_exp <- append(each_count_se_exp,summary(ldpred_lm_X_result)$coefficients[,2][2])
  
  each_cycle_ldpred_score_2 <- as.data.frame(datafiles_2[each_cycle_ldpred_score_index])[2]
  ldpred_score_2 <- append(ldpred_score_2,each_cycle_ldpred_score_2)
  
  ldpred_lm_Y_result <- lm(unlist(Y_i_2) ~ unlist(each_cycle_ldpred_score_2))
  each_count_b_out <- append(each_count_b_out,summary(ldpred_lm_Y_result)$coefficients[,1][2])
  each_count_se_out <- append(each_count_se_out,summary(ldpred_lm_Y_result)$coefficients[,2][2])
  
}

# Convert results to unlisted format for easier handling
each_count_b_exp <- unname(each_count_b_exp)
each_count_se_exp <- unname(each_count_se_exp)
each_count_b_out <- unname(each_count_b_out)
each_count_se_out <- unname(each_count_se_out)

PRS_ldpred_beta <- c()
PRS_ldpred_beta_se <- c()

# Calculate Wald ratios for each cycle and store beta coefficients and standard errors
for(i in 1:cycle_num){
  temp_PRS_ldpred_result <- mr_wald_ratio(each_count_b_exp[i], each_count_b_out[i], each_count_se_exp[i], each_count_se_out[i])
  PRS_ldpred_beta <- append(PRS_ldpred_beta, temp_PRS_ldpred_result$b)
  PRS_ldpred_beta_se <- append(PRS_ldpred_beta_se,temp_PRS_ldpred_result$se)
  
}

# Combine the results into a data frame and write the results to a file
PRS_ldpred_result <- cbind(PRS_ldpred_beta,PRS_ldpred_beta_se)
colnames(PRS_ldpred_result) <- c("beta","se")
write.table(PRS_ldpred_result, file = paste0("/gpfs/gibbs/project/zhao/qh63/R/summary_statistics/PRS_ldpred_",count,".txt"), sep = " ",row.names = T, col.names = T ,quote = TRUE)



#########################################

