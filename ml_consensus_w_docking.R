##########################################
# Fisher 80 for any target
##########################################
library(dplyr)
library(neuralnet)
library(randomForest)
library(caret)
library(e1071)

############################
# Parameters
############################
num_compounds <- 36273

###########################
# Read Sig2Lead File
###########################

dat_s2l <- read.csv(file="./s2l_files/my_Candidates.csv")
dat_s2l <- dat_s2l[order(dat_s2l$User.added.Candidate),]

################################
# Read LINCS Candidates File
###############################

lincs <- read.csv(file="./s2l_files/lincs_Candidates.csv")
max_conc_tail <- max(lincs$Concordance)

##########################
# Read Docking File
##########################

dat_docking <- read.csv("./s2l_files/docking_file.csv")
dat_docking <- dat_docking[,c("compound", "autodock")]
dat_docking <- dat_docking[order(dat_docking$compound),]
dat_docking[dat_docking==0] <- .00001

###########################
# Merge Docking with S2L
###########################

dat <- data.frame(dat_s2l, dat_docking$autodock)
colnames(dat)[7] <- "autodock"

#############################################
# Compute 11 and 13 Features for Test Target
#############################################

s2l <- dat %>% dplyr::arrange(desc(Concordance)) %>% dplyr::arrange(desc(Similarity)) %>% head(nrow(dat))

ind <- which(is.na(s2l), arr.ind = T)
s2l[ind] <- 0  


s2l <- s2l %>% mutate("neg_log10_autodock" = -log10((s2l$autodock)))# %>% 
s2l <- s2l %>% #mutate("neg_log10_autodock" = -log10(s2l$autodock)) %>% 
  mutate("neg_log10_autodock" = -log10(s2l$autodock)) %>% 
  mutate("median_normalized_concordance" = s2l$Concordance/median(s2l$Concordance)) %>% 
  mutate("max_over_one_minus_conc" = (max_conc_tail / (1 - s2l$Concordance))) %>% 
  mutate("max_minus_conc_over_max" = ((max_conc_tail - s2l$Concordance) / (max_conc_tail))) %>% 
  #mutate("median_normalized_tanimoto" = s2l$Tanimoto/median(s2l$Tanimoto[i])) %>%
  mutate("median_normalized_tanimoto" = s2l$Similarity/median(s2l$Similarity)) %>%
  #mutate("median_normalized_autodock" = s2l$neg_log10_autodock/ (-log10(median(dat$autodock))) ) 
  mutate("median_normalized_autodock" = s2l$neg_log10_autodock/ (-log10(median(s2l$autodock))) ) 
#

num_compounds <- nrow(dat)
s2l$Rank_S2l <- 1:num_compounds
s2l <- s2l[order(s2l$neg_log10_autodock, decreasing = TRUE),]
s2l$Rank_Autodock <- 1:num_compounds 
s2l <-  s2l[order(s2l$Concordance, decreasing = TRUE),]
s2l$Rank_Concordance <- 1:num_compounds 

s2l$Fraction_Rank_S2L <- s2l$Rank_S2l / num_compounds
s2l$Fraction_Rank_Autodock <- s2l$Rank_Autodock / num_compounds
s2l$Fraction_Rank_Concordance <- s2l$Rank_Concordance / num_compounds

s2l$Fisher_2 <- -log10(s2l$Fraction_Rank_S2L) + -log10(s2l$Fraction_Rank_Autodock)
s2l$Fisher_3 <-  -log10(s2l$Fraction_Rank_S2L) + -log10(s2l$Fraction_Rank_Autodock) + -log10(s2l$Fraction_Rank_Concordance)

colnames(s2l)[6] <- "Tanimoto"
s2l <- s2l[order(s2l$Rank_S2l, decreasing=FALSE),]

#################################################################################################
########################################
# Test Data
########################################
test_dat <- as.matrix(s2l[,c("Tanimoto",
                   "Concordance",
                   "neg_log10_autodock",
                   "median_normalized_concordance",
                   "median_normalized_tanimoto",
                   "median_normalized_autodock",
                   "Fraction_Rank_S2L",
                   "Fraction_Rank_Autodock",
                   "Fraction_Rank_Concordance",
                   "Fisher_3",
                   "Fisher_2",
                   "max_over_one_minus_conc",
                   "max_minus_conc_over_max")])


test_dat[which(!is.finite(test_dat))] <- 6



#######################################
# Load Models
########################################

model_rf_11 <- list()
model_rf_13_b <- list()
model_nn_11 <- list()
model_nn_13 <- list()

for (i in 1:20){
  model_rf_11[[i]] <- readRDS(file=paste("./models/individual_models/rf_11/rf_11_",i,".rds",sep=""))
  model_rf_13_b[[i]] <- readRDS(file=paste("./models/individual_models/rf_13b/rf_13b_", i,".rds",sep=""))
  model_nn_11[[i]] <- readRDS(file=paste("./models/individual_models/nn_11/nn_11_",i,".rds",sep=""))
  model_nn_13[[i]] <- readRDS(file=paste("./models/individual_models/nn_13/nn_13_",i,".rds",sep=""))
}

##########################################
# Compute Predicitions 4 models:20 targets 
##########################################
predict_rf_11 <- list()
predict_rf_13 <- list()
predict_nn_11 <- list()
predict_nn_13 <- list()

rf_13_o <- list()
nn_11_o <- list()
nn_13_o <- list()
rf_11_o_f <- list()
rf_13_o_f <- list()
nn_11_o_f <- list()
nn_13_o_f <- list()
fisher_4_ml <- list()
rf_11_o <- list()

for (i in 1:20){
  #i <- 1
  ####### RF11
  predict_rf_11[[i]] <- predict(model_rf_11[[i]], newdata = test_dat[,c(1:11)], type="prob")
  rf_11_o[[i]] <- data.frame(s2l[,1], predict_rf_11[[i]])
 
  ########RF13
  predict_rf_13[[i]] <- predict(model_rf_13_b[[i]], newdata = test_dat[,c(1:13)], type="prob")
  rf_13_o[[i]] <- data.frame(s2l[,1], predict_rf_13[[i]])
  
  ########NN11
  predict_nn_11[[i]] <- compute(model_nn_11[[i]], test_dat[,c(1:11)])
  nn_11_o[[i]] <- data.frame(s2l[,1], predict_nn_11[[i]]$net.result)
  
  ########NN13
  predict_nn_13[[i]] <- compute(model_nn_13[[i]], test_dat[,c(1:13)])
  nn_13_o[[i]] <- data.frame(s2l[,1], predict_nn_13[[i]]$net.result)
  #######################################################################################################################
  
  
  ##################### Order Data and Assign Rank ###################################################################################################
   
  rf_11_o[[i]] <- rf_11_o[[i]][order(rf_11_o[[i]]$TRUE., decreasing=TRUE),]
  rf_13_o[[i]] <- rf_13_o[[i]][order(rf_13_o[[i]]$TRUE., decreasing=TRUE),]
  nn_11_o[[i]] <- nn_11_o[[i]][order(nn_11_o[[i]]$X2, decreasing=TRUE),]
  nn_13_o[[i]] <- nn_13_o[[i]][order(nn_13_o[[i]]$X2, decreasing=TRUE),]
  ##################### Add Rank Column #############################################################################################################################
  
  rf_11_o[[i]]$Rank <- 1:nrow(rf_11_o[[i]])
  rf_13_o[[i]]$Rank <- 1:nrow(rf_13_o[[i]])
  nn_11_o[[i]]$Rank <- 1:nrow(nn_11_o[[i]])
  nn_13_o[[i]]$Rank <- 1:nrow(nn_13_o[[i]])
  
  ################## Compute Fisher Probabilities ####################################################
  
  rf_11_o[[i]]$Prob_rf_11 <- -log10(rf_11_o[[i]]$Rank / nrow(rf_11_o[[i]]))
  rf_13_o[[i]]$Prob_rf_13 <- -log10(rf_13_o[[i]]$Rank / nrow(rf_13_o[[i]]))
  nn_11_o[[i]]$Prob_nn_11 <- -log10(nn_11_o[[i]]$Rank / nrow(nn_11_o[[i]]))
  nn_13_o[[i]]$Prob_nn_13 <- -log10(nn_13_o[[i]]$Rank / nrow(nn_13_o[[i]]))
  

  #################### Reorder to Orginal Order ##########################################
  
  
  rf_11_o_f[[i]] <- rf_11_o[[i]][order(rf_11_o[[i]][,1], decreasing = FALSE),]
  rf_13_o_f[[i]] <- rf_13_o[[i]][order(rf_13_o[[i]][,1], decreasing = FALSE),]
  nn_11_o_f[[i]] <- nn_11_o[[i]][order(nn_11_o[[i]][,1], decreasing = FALSE),]
  nn_13_o_f[[i]] <- nn_13_o[[i]][order(nn_13_o[[i]][,1], decreasing = FALSE),]
  

  
  
  ############### Summation of Fisher Probabiites #####################################################
  colnames(rf_11_o_f[[1]])[1] <- "Compound"
  fisher_4_ml[[i]] <- data.frame(rf_11_o_f[[i]]$Prob_rf_11,  rf_13_o_f[[i]]$Prob_rf_13, nn_11_o_f[[i]]$Prob_nn_11,nn_13_o_f[[i]]$Prob_nn_13) 
 
  print(paste("Finished iteration", i))
}

fisher_80_ml <- do.call("cbind", fisher_4_ml)
fisher_80_ml$Fisher_Sum <- rowSums(fisher_80_ml)
fisher_80_ml$Compound <- rf_11_o_f[[1]]$Compound
fisher_80_ml_o <- fisher_80_ml[order(fisher_80_ml$Fisher_Sum, decreasing=TRUE),]
fisher_80_ml_o <- fisher_80_ml[,c(82,1:81)]

write.csv(fisher_80_ml_o, file="./results/ml_consensus_output.csv")
