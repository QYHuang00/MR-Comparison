#! /bin/Rscript

# This script simulates data from scratch and calculates the estimated Mendelian randomization(MR) effect.
# The goal is to compare the performance of different methods in estimating the effect under various conditions.
# The MR methods and PRS Methods used include: mr_sisVIVE, mr_2SLS, mr_ivw, mr_median, mr_mode, mr_mix, mr_egger, mr_raps, PRS_PT, PRS_Oracle and PRS_LDPred.
# For the LDPred method, this script only handles preparing the data into PLINK format, while the actual LDPred calculations will be done in a subsequent script.

# Set working directory
setwd("/gpfs/gibbs/project/zhao/qh63/R")

# Load necessary libraries for Mendelian randomization analysis and data generation
library(CorBin)
library(MASS)
library(MendelianRandomization)
library(sisVIVE)
library(RobustIV)
library(ivreg)
library(genio)
library(tibble)


# Initialize vectors to store means and standard errors for multiple MR methods and PRS results
var_names <- c("VIVE", "mr_2SLS", "mr_ivw", "mr_ivw2", "mr_median", 
               "mr_mode", "mr_mode2", "mr_mix", "mr_egger", "mr_raps", 
               "PRS_PT", "PRS_Oracle")

for (name in var_names) {
  assign(paste0(name, "_beta"), c())
  assign(paste0(name, "_beta_se"), c())
}

count <- 1  ## Since there are many scripts running simultaneously, this 'count' variable helps users keep track of and manage the output files more effectively.
output_name <- paste0("output_",count)   
pval <- c()

# The loop will iterate to simulate data and run different MR methods.
k <- 0
while(k<100){
  k<- k+1
  print(k)
  
  # Simulate genotype data (G_ij), generate standardized genotype matrix, and other variables for MR analysis
  individual_num <- 5000*2 
  block_size <- 100
  block_num <- 100 
  SNP_Num <- block_size*block_num 
  h_x2 <- 0.5  # Heritability of X
  h_y2 <- 0    # Heritability of Y
  beta <- 0    # True causal effect
  
  pi_x <- 0.0005
  pi_10 <- pi_x
  pi_00 <- 1 - pi_10
  rho_alphagamma <- 0 
  
  ##########################################################################
  
  # Generate allele frequencies and genotype matrix (G_ij)
  # The genotype matrix will be standardized to ensure proper scaling

  A_ij <- 0
  allele_f <- c()
  for(b in 1:block_num){
    f_b <- runif(1,0.05,0.95)
    rho <- 0
    A_ij_block <- cBern(2*individual_num, rep(f_b,block_size), rho, type="DCP")
    A_ij <- cbind(A_ij,A_ij_block)
    allele_f <- append(allele_f,f_b)
  }
  A_ij <- A_ij[,-1]
  
  A_ij_add <- matrix(0,individual_num,dim(A_ij)[2])
  A_ij_add <- A_ij[2*(1:individual_num),]+A_ij[2*(1:individual_num)-1,]
  G_ij <- A_ij_add
  
  for(i in 1:ncol(G_ij)){
    if(var(G_ij[,i])==0){
      k <- k-1
      next
    }
  }

  G_ij_standardized <- matrix(0, nrow = dim(G_ij)[1], ncol = dim(G_ij)[2])
  for(b in 1:block_num){
    left <- (b-1)*block_size+1
    right <- (b-1)*block_size+block_size
    
    G_ij_standardized[,left:right] <- (G_ij[,left:right]-2*allele_f[b])/sqrt((2*allele_f[b]*(1-allele_f[b])))
  }
  
  G_ij_standardized_1 <- G_ij_standardized[1:(nrow(G_ij)/2),]
  G_ij_standardized_2 <- G_ij_standardized[((nrow(G_ij)/2)+1):nrow(G_ij),]
  
  # GenerateU_i,eps_xi,eps_yi
  U_i_1 <- as.matrix(rnorm(individual_num/2, mean = 0, sd = sqrt((1-h_x2)/2)))
  eps_xi_1 <- as.matrix(rnorm(individual_num/2, mean = 0, sd = sqrt((1-h_x2)/2)))
  eps_yi_1 <- as.matrix(rnorm(individual_num/2, mean = 0, sd = sqrt(0.5+0.5*h_x2-h_y2-beta^2)))
  
  U_i_2 <- as.matrix(rnorm(individual_num/2, mean = 0, sd = sqrt((1-h_x2)/2)))
  eps_xi_2 <- as.matrix(rnorm(individual_num/2, mean = 0, sd = sqrt((1-h_x2)/2)))
  eps_yi_2 <- as.matrix(rnorm(individual_num/2, mean = 0, sd = sqrt(0.5+0.5*h_x2-h_y2-beta^2)))
  
  #Generate(sigma_x)^2,(sigma_y)^2
  sigma_x2 <- h_x2/(SNP_Num*pi_10) #(sigma_x)^2
  
  #Generate alpha_j,gamma_j
  alpha <- matrix(0,SNP_Num,1)
  gamma <- matrix(0,SNP_Num,1)
  compIdx <- sample.int(2, SNP_Num, replace=T, prob=c(pi_00, 1-pi_00))
  alpha[compIdx==2,1] <- rnorm(sum(compIdx==2), mean=0, sd=sqrt(sigma_x2))
  
  #Generate X_i,Y_i
  X_i_1 <- t(crossprod(alpha,t(G_ij_standardized_1))) + U_i_1 + eps_xi_1
  Y_i_1 <- t(crossprod(gamma,t(G_ij_standardized_1))) + beta*X_i_1 + U_i_1 + eps_yi_1
  print(sd(X_i_1))
  print(sd(Y_i_1))
  
  X_i_2 <- t(crossprod(alpha,t(G_ij_standardized_2))) + U_i_2 + eps_xi_2
  Y_i_2 <- t(crossprod(gamma,t(G_ij_standardized_2))) + beta*X_i_2 + U_i_2 + eps_yi_2
  print(sd(X_i_2))
  print(sd(Y_i_2))
  
  
  
  ###########################################################################
  
  #Select independent instrument variables
  
  ###### step 1: lm regression
  
  G_ij_p <- function(my_model) {
    p_values <- summary(my_model)$coefficients[,4]
    return(p_values[2])         
  }
  
  # Initialize matrices to store results for regression on dataset 1
  lm_reg_p_1 <- matrix(-1,SNP_Num,1)
  p_bound <- 0.05*10/SNP_Num 
  lm_reg_beta_1 <- matrix(-1,SNP_Num,1)
  lm_reg_se_1 <- matrix(-1,SNP_Num,1)
  
  # Initialize matrices to store results for regression on dataset 2
  lm_reg_p_2 <- matrix( -1,SNP_Num,1)
  p_bound <- 0.05*10/SNP_Num 
  lm_reg_beta_2 <- matrix(-1,SNP_Num,1)
  lm_reg_se_2 <- matrix(-1,SNP_Num,1)
  
  # Loop over all SNPs and perform linear regression on both datasets
  for(i in 1:SNP_Num){
    temp_lm_1 <- lm(X_i_1 ~ G_ij_standardized_1[,i])
    lm_reg_p_1[i,1] <- summary(temp_lm_1)$coefficients[,4][2]
    lm_reg_se_1[i,1] <- summary(temp_lm_1)$coefficients[,2][2]
    lm_reg_beta_1[i,1] <- summary(temp_lm_1)$coefficients[,1][2]
    
    temp_lm_2 <- lm(X_i_2 ~ G_ij_standardized_2[,i])
    lm_reg_p_2[i,1] <- summary(temp_lm_2)$coefficients[,4][2]
    lm_reg_se_2[i,1] <- summary(temp_lm_2)$coefficients[,2][2]
    lm_reg_beta_2[i,1] <- summary(temp_lm_2)$coefficients[,1][2]
    
  }
  
  correlated_X_index <- as.matrix(which(lm_reg_p_1[,1] <= p_bound))
  
  # Check if any SNPs are identified as significant
  if(length(correlated_X_index)<1){
    k <- k-1
    next
  }
  
  
  ###### step 2: clumping
  
  # Initialize variables for clumping procedure
  S <- 0
  S_X_p <- 0
  lead_SNP <- 0
  ind_IV <- 0 
  
  # Extract the SNPs identified in the previous step (correlated_X_index) and their p-values
  S <- G_ij_standardized_1[,correlated_X_index]
  S_X_p <- as.matrix(lm_reg_p_1[correlated_X_index,1])
  
  # Identify the lead SNP (the SNP with the smallest p-value)
  lead_SNP <- as.matrix(G_ij_standardized_1[,which(lm_reg_p_1[,1] == min(S_X_p[,1]))])[,1]
  ind_IV <- as.matrix(lead_SNP)
  ind_IV_index <- c(which(lm_reg_p_1[,1] == min(S_X_p[,1])))[1]
  
  # Flag to determine if the clumping process should be skipped due to errors
  skip_to_next_clumping <- FALSE
  
  # Begin clumping process with a threshold set to 0
  threshold = 0 
  while(length(S) > 0){
    # If the number of SNPs exceeds half the individual number, perform clumping
    if(length(S) > (individual_num/2) + 1){
      threshold <- 1
      delete_index <- c()
      
      # Identify SNPs highly correlated with the lead SNP (r^2 >= 0.1) to remove them
      for(i in 1:dim(S_X_p)[1]){
        if((cor(lead_SNP,S)[i])^2 >= 0.1){
          delete_index <- append(delete_index,i)
        }
      }
      
      # Debugging output
      print("test")
      print(delete_index)
      print(dim(S))
      
      # If delete_index is empty, output error information and skip clumping
      if(is.null(delete_index)){
          print("error")
          print(dim(lead_SNP))
          print(head(lead_SNP))
          print((cor(lead_SNP,S))^2)
          print(dim(S_X_p))
          print(length(S))
          print((individual_num/2) + 1)
          skip_to_next_clumping <- TRUE
          break
      }else{
        # Remove correlated SNPs from the dataset and update the SNPs and p-values
          S <- S[,-delete_index]
          S_X_p <- as.matrix(S_X_p[-delete_index,1])
      }
      
      # If there are remaining SNPs, update lead SNP and add to independent IV set
      if(length(S)!=0){
        lead_SNP <- as.matrix(G_ij_standardized_1[,which(lm_reg_p_1[,1] == min(S_X_p[,1]))])[,1]
        ind_IV <- cbind(ind_IV,lead_SNP)  
        ind_IV_index <- append(ind_IV_index,which(lm_reg_p_1[,1] == min(S_X_p[,1]))[1])

      }
    }else if(length(S) == (individual_num/2) & threshold == 0){
      # If number of SNPs equals half of individuals and threshold not set, add lead SNP to IV
      lead_SNP <- as.matrix(G_ij_standardized_1[,which(lm_reg_p_1[,1] == min(S_X_p[,1]))])[,1]
      ind_IV <- cbind(ind_IV,lead_SNP)  
      ind_IV_index <- append(ind_IV_index,which(lm_reg_p_1[,1] == min(S_X_p[,1]))[1])
      S <- c() # Clear SNPs to end the loop
    }else{
      S <- c() # Clear SNPs if condition is not met
    }
  }
  
  # If the clumping process encounters issues, skip this iteration and continue to the next simulation
  if(skip_to_next_clumping){
      k <- k-1
      next
  }
  
  # Extract corresponding SNPs for dataset 2 based on the indices from dataset 1 clumping
  ind_IV_2 <- G_ij_standardized_2[,ind_IV_index]
  
  # Error handling: ensure column names are properly set and catch any errors
   tryCatch(
        {
        colnames(ind_IV) <- ind_IV_index      
        colnames(ind_IV_2) <- ind_IV_index
        },
        error=function(e) {
            message('Ind IV Error Occurred')
            print(dim(ind_IV))
            print(head(ind_IV))
            print(ind_IV_index)
            break
            #skip_to_next <<- TRUE
        }
    )
  

  ################################################################################
  
  # The following section begins the implementation of various MR methods and PRS methods.
  # It calculates the estimated effect and error for each MR method.
  
  # Initialize matrices to store exposure (X) and outcome (Y) data for Mendelian Randomization (MR) analysis
  exposure <- matrix(-1,dim(ind_IV)[2],4)  #X
  colnames(exposure) <- c("SNP","beta.exposure","se.exposure","pval.exposure")
  outcome <- matrix(-1,dim(ind_IV)[2],4) #Y
  colnames(outcome) <- c("SNP","beta.exposure","se.exposure","pval.exposure")
  
  # Populate exposure and outcome matrices with data for each SNP
  for(i in 1:dim(ind_IV)[2]){
    exposure[i,"SNP"] <- ind_IV_index[i] 
    exposure[i,"beta.exposure"] <- lm_reg_beta_1[ind_IV_index[i], 1] 
    exposure[i,"se.exposure"] <- lm_reg_se_1[ind_IV_index[i], 1] 
    exposure[i,"pval.exposure"] <- lm_reg_p_1[ind_IV_index[i], 1] 
    
    outcome[i,"SNP"] <- ind_IV_index[i] 
    temp_lm2 <- summary(lm(Y_i_2 ~ind_IV_2[,i]))$coefficients
    outcome[i,"beta.exposure"] <- temp_lm2[,1][2] 
    outcome[i,"se.exposure"] <- temp_lm2[,2][2] 
    outcome[i,"pval.exposure"] <- temp_lm2[,4][2] 
  }
  
  # Prepare MR input data using the exposure and outcome matrices
  MR_dat <- mr_input(bx = exposure[,"beta.exposure"],bxse = exposure[,"se.exposure"],
                     by = outcome[,"beta.exposure"],byse = outcome[,"se.exposure"],
                     snps = exposure[,"SNP"])
  
  # Initialize flag to control flow in case of exceptions
  skip_to_next <- FALSE
  
  # Try-catch block for robust MR analysis using TSHT 
  tryCatch(
        {
        TSHT_result <- TSHT(Y = Y_i_1, D = X_i_1, Z = ind_IV)
        mr_sisVIVE_beta <- append(mr_sisVIVE_beta, TSHT_result$betaHat[1])
        mr_sisVIVE_beta_se <- append(mr_sisVIVE_beta_se, TSHT_result$beta.sdHat[1])
  
        },
        error=function(e) {
            message('An Error Occurred')
            print(e)
            print(dim(ind_IV))
            skip_to_next <<- TRUE
        },
        warning=function(w) {
            message('A Warning Occurred')
            print(w)
            skip_to_next <<- TRUE
        }
    )


  mr_2SLS_result <- ivreg(Y_i_1 ~ X_i_1 | ind_IV , method="OLS")
  mr_2SLS_beta[k]<-summary(mr_2SLS_result)$coefficients[,1][2]
  mr_2SLS_beta_se[k]<-summary(mr_2SLS_result)$coefficients[,2][2]
  
  mr_ivw_result <- MendelianRandomization:::mr_ivw(MR_dat, model='default', weights='delta')
  mr_ivw_beta[k] <- mr_ivw_result@Estimate
  mr_ivw_beta_se[k] <- mr_ivw_result@StdError
  
  # Try-catch block for robust MR analysis using MR Median  
  tryCatch({
    mr_median_result <- MendelianRandomization:::mr_median(MR_dat)
    mr_median_beta[k]<- mr_median_result@Estimate
    mr_median_beta_se[k] <- mr_median_result@StdError
  },
  error = function(e) {skip_to_next <<- TRUE}
  )
  
  # Try-catch block for robust MR analysis using MR Mode
  tryCatch({
    mr_mode_result <- MendelianRandomization:::mr_mbe(MR_dat, weighting='weighted', stderror='delta')  #mr_mode
    mr_mode_beta[k]<- mr_mode_result@Estimate
    mr_mode_beta_se[k] <- mr_mode_result@StdError
  },
  error = function(e) {skip_to_next <<- TRUE}
  )

  
  if(skip_to_next){
    k <- k-1
    next
  }
  
  
  library(MRMix)
  mr_Mix_result <- MRMix(exposure[,"beta.exposure"], outcome[,"beta.exposure"], 
                         exposure[,"se.exposure"], outcome[,"se.exposure"])
  mr_mix_beta <- append(mr_mix_beta, mr_Mix_result$theta)
  mr_mix_beta_se <- append(mr_mix_beta_se, mr_Mix_result$SE_theta)
  
  mr_egger_result <- mr_egger(MR_dat)
  mr_egger_beta <- append(mr_egger_beta, mr_egger_result@Estimate)
  mr_egger_beta_se <- append(mr_egger_beta_se, mr_egger_result@StdError.Est)
  
  library(mr.raps)
  mr_raps_result <- mr.raps(exposure[,"beta.exposure"], outcome[,"beta.exposure"], 
                            exposure[,"se.exposure"], outcome[,"se.exposure"])
  mr_raps_beta <- append(mr_raps_beta, mr_raps_result$beta.hat)
  mr_raps_beta_se <- append(mr_raps_beta_se, mr_raps_result$beta.se)
  
  
  ##############################################
  
  # Calculate Mendelian Randomization using PRS as independent IVs
  
  # PRS_P+T
  library(TwoSampleMR)
  PRS_PT_score_1 <- t(as.matrix(exposure[,"beta.exposure"])) %*% t(G_ij_standardized_1[,ind_IV_index])
  PRS_PT_score_2 <- t(as.matrix(exposure[,"beta.exposure"])) %*% t(G_ij_standardized_2[,ind_IV_index])
  PT_lm_X_result <- lm(X_i_1 ~ t(PRS_PT_score_1))
  PT_lm_Y_result <- lm(Y_i_2 ~ t(PRS_PT_score_2))
  PRS_PT_result <- mr_wald_ratio(summary(PT_lm_X_result)$coefficients[,1][2],
                                 summary(PT_lm_Y_result)$coefficients[,1][2],
                                 summary(PT_lm_X_result)$coefficients[,2][2],
                                 summary(PT_lm_Y_result)$coefficients[,2][2])
  PRS_PT_beta <- append(PRS_PT_beta, PRS_PT_result$b)
  PRS_PT_beta_se <- append(PRS_PT_beta_se, PRS_PT_result$se)
  
  # PRS_Oracle
  PRS_Oracle_score_1 <- t(as.matrix(alpha[ind_IV_index])) %*% t(G_ij_standardized_1[,ind_IV_index])
  PRS_Oracle_score_2 <- t(as.matrix(alpha[ind_IV_index])) %*% t(G_ij_standardized_2[,ind_IV_index])
  Oracle_lm_X_result <- lm(X_i_1 ~ t(PRS_Oracle_score_1))
  Oracle_lm_Y_result <- lm(Y_i_2 ~ t(PRS_Oracle_score_2))
  PRS_Oracle_result <- mr_wald_ratio(summary(Oracle_lm_X_result)$coefficients[,1][2],
                                     summary(Oracle_lm_Y_result)$coefficients[,1][2],
                                     summary(Oracle_lm_X_result)$coefficients[,2][2],
                                     summary(Oracle_lm_Y_result)$coefficients[,2][2])
  PRS_Oracle_beta <- append(PRS_Oracle_beta, PRS_Oracle_result$b)
  PRS_Oracle_beta_se <- append(PRS_Oracle_beta_se, PRS_Oracle_result$se)
  
  # Preparing the data in PLINK format to facilitate subsequent LDPred calculations
  write_plink(paste0("./PLink_data/",output_name,"/test_1_",k), t(G_ij[1:(nrow(G_ij)/2),]))

  tib <- tibble(
    chr = 1,
    id = 1:SNP_Num,
    posg = 0,
    pos = 1:SNP_Num,
    alt = 'G',
    ref = 'A'
  )
  
  write_bim(paste0("./PLink_data/",output_name,"/test_1_",k,".bim"), tib)
  
  summary_data_1 <- data.frame(CHR = rep.int(1,SNP_Num),
                               POS = 1:SNP_Num,
                               SNP_ID = 1:SNP_Num,
                               REF = rep.int("A",SNP_Num),
                               ALT = rep.int("G",SNP_Num),
                               REF_FRQ = rep(allele_f,each = block_size),
                               PVAL = lm_reg_p_1,
                               BETA = lm_reg_beta_1,
                               SE = lm_reg_se_1,
                               N = rep(individual_num/2,SNP_Num))
  
  write.table(summary_data_1, file = paste0("./PLink_data/",output_name,"/ssf_1_",k,".txt"), sep = " ",row.names = FALSE, quote = FALSE)
  
  write_plink(paste0("./PLink_data/",output_name,"/test_2_",k), t(G_ij[((nrow(G_ij)/2)+1):nrow(G_ij),]))
  write_bim(paste0("./PLink_data/",output_name,"/test_2_",k,".bim"), tib)
  summary_data_2 <- data.frame(CHR = rep.int(1,SNP_Num),
                               POS = 1:SNP_Num,
                               SNP_ID = 1:SNP_Num,
                               REF = rep.int("A",SNP_Num),
                               ALT = rep.int("G",SNP_Num),
                               REF_FRQ = rep(allele_f,each = block_size),
                               PVAL = lm_reg_p_2,
                               BETA = lm_reg_beta_2,
                               SE = lm_reg_se_2,
                               N = rep(individual_num/2,SNP_Num))
  
  write.table(summary_data_2, file = paste0("./PLink_data/",output_name,"/ssf_2_",k,".txt"), sep = " ",row.names = FALSE, quote = FALSE)
  
  # Save the values of X_i_1, Y_i_1, X_i_2, and Y_i_2 for subsequent data processing
  write.table(X_i_1,file = paste0("./summary_statistics/X_i/X_i_1_",count,"_",k,".txt"))
  write.table(Y_i_1,file = paste0("./summary_statistics/Y_i/Y_i_1_",count,"_",k,".txt"))
  
  write.table(X_i_2,file = paste0("./summary_statistics/X_i/X_i_2_",count,"_",k,".txt"))
  write.table(Y_i_2,file = paste0("./summary_statistics/Y_i/Y_i_2_",count,"_",k,".txt"))
  
}

# Write out the estimated effect and standard error for each MR method

mr_sisVIVE <- data.frame(beta = unlist(mr_sisVIVE_beta), se = unlist(mr_sisVIVE_beta_se))
write.table(mr_sisVIVE, file = paste0("./summary_statistics/mr_sisVIVE_",count,".txt"))

mr_2SLS <- data.frame(beta = unlist(mr_2SLS_beta), se = unlist(mr_2SLS_beta_se))
write.table(mr_2SLS, file = paste0("./summary_statistics/mr_2SLS_",count,".txt"))

mr_ivw_output <- data.frame(beta = unlist(mr_ivw_beta), se = unlist(mr_ivw_beta_se))
write.table(mr_ivw_output, file = paste0("./summary_statistics/mr_ivw_",count,".txt"))

mr_median_output <- data.frame(beta = unlist(mr_median_beta), se = unlist(mr_median_beta_se))
write.table(mr_median_output, file = paste0("./summary_statistics/mr_median_",count,".txt"))


mr_mode_output <- data.frame(beta = unlist(mr_mode_beta), se = unlist(mr_mode_beta_se))
write.table(mr_mode_output, file = paste0("./summary_statistics/mr_mode_",count,".txt"))

mr_mix <- data.frame(beta = unlist(mr_mix_beta), se = unlist(mr_mix_beta_se))
write.table(mr_mix, file = paste0("./summary_statistics/mr_mix_",count,".txt"))

mr_egger_output <- data.frame(beta = unlist(mr_egger_beta), se = unlist(mr_egger_beta_se))
write.table(mr_egger_output, file = paste0("./summary_statistics/mr_egger_",count,".txt"))

mr_raps_output <- data.frame(beta = unlist(mr_raps_beta), se = unlist(mr_raps_beta_se))
write.table(mr_raps_output, file = paste0("./summary_statistics/mr_raps_",count,".txt"))

PRS_PT <- data.frame(beta = unlist(PRS_PT_beta), se = unlist(PRS_PT_beta_se))
write.table(PRS_PT, file = paste0("./summary_statistics/PRS_PT_",count,".txt"))

PRS_Oracle <- data.frame(beta = unlist(PRS_Oracle_beta), se = unlist(PRS_Oracle_beta_se))
write.table(PRS_Oracle, file = paste0("./summary_statistics/PRS_Oracle_",count,".txt"))

# Save the session for future use
save.image(file = paste0("./RData/output",count,".RData"))
