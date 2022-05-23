require(dplyr)
require(Morpho)
require(microbenchmark)
require(usedist)
require(stringr)
require(ggplot2)
require(svDialogs)
require(otuSummary)
require(fs)
require(Rmisc)

#reads MODEL data .txt file exported from ID unknown script
model_data <- read.table(choose.files(caption = "Choose data file (.txt) with MODEL", multi = FALSE),sep="\t",header = T)
names(model_data)[1] <- "ind"
names(model_data)[2] <- "gr"
model_data <- data.frame(row.names = model_data$ind,model_data)

model_data_alpha <- model_data

cross_val_no <- dlg_input("Please specify the number of cross-validation iterations:", Sys.info()["cross_val_no"])$res
cross_val_no <- as.numeric(cross_val_no)

## 60% of the model size for known cross_val dataset
smp_size <- floor(0.6 * nrow(model_data_alpha))

#create empty tables for final scores
accuracy_val <- vector("numeric", length(cross_val_no))

#####################################################
#This for loop runs cross-validation iterated on 'l' random 60/40 data partitions
#####################################################

pb_f <- winProgressBar(title="Progress bar - cross-validations", label="% tests done", min=0, max=100, initial=0)

for (l in 1:cross_val_no) {
  
  #first split the model data 60/40 (known/unknown)
  
  train_ind <- sample(seq_len(nrow(model_data_alpha)), size = smp_size)
  
  model_data <- model_data_alpha[train_ind, ]
  data <- model_data_alpha[-train_ind, ]
  
  known_IDs <- data.frame(data$ind,data$gr)
  names(known_IDs)[1] <- "ind"
  names(known_IDs)[2] <- "gr"
  data <- data.frame(data[,-1:-2])
  data_w_blanks <- data.frame(rep("unknown",nrow(data)),rep("unknown",nrow(data)),data)
  names(data_w_blanks)[1] <- "ind"
  names(data_w_blanks)[2] <- "gr"
  
  #then recombine known and unknown splits into one dataset
  total_dataset <- rbind(data_w_blanks,model_data)
  
  #gotta avoid groups with <3 observations for CVA to work
  SINGLETONchecktable <- data.frame(table(model_data$gr))
  SINGLETONchecktest <-  SINGLETONchecktable$Freq
  
  
  while(any(SINGLETONchecktest < 2)) {
    
    train_ind <- sample(seq_len(nrow(model_data_alpha)), size = smp_size)
    
    model_data <- model_data_alpha[train_ind, ]
    data <- model_data_alpha[-train_ind, ]
    
    known_IDs <- data.frame(data$ind,data$gr)
    names(known_IDs)[1] <- "ind"
    names(known_IDs)[2] <- "gr"
    data <- data.frame(data[,-1:-2])
    
    data_w_blanks <- data.frame(rep("unknown",nrow(data)),rep("unknown",nrow(data)),data)
    names(data_w_blanks)[1] <- "ind"
    names(data_w_blanks)[2] <- "gr"
    total_dataset <- rbind(data_w_blanks,model_data)
    
    
    SINGLETONchecktable <- data.frame(table(model_data$gr))
    SINGLETONchecktest <-  SINGLETONchecktable$Freq
    
  }
  
  MODEL_DERIV <- data.frame(model_data[,-1:-2])
  
  data_to_ID <- data.frame(rep("unknown",nrow(data)),data)
  names(data_to_ID)[1] <- "gr"
  
  
  f = file()
  sink(file=f) ## silence upcoming output using anonymous file connection
  
  model_CVA <- CVA(MODEL_DERIV,model_data$gr)
  
  sink() ## undo silencing
  close(f)
  
  result_IDs <- vector('numeric',length(row.names(data)))
  
  
  groups_for_test <- unique(model_data$gr)
  
  #######################################################################
  #####THIS FOR LOOP 'fz' iterates per unknown observation
  #######################################################################
  
  for (fz in 1:length(row.names(data))) {
    
    model_test_data <- data.frame(model_data$gr,MODEL_DERIV)
    names(model_test_data)[1] <- "gr"
    #add unknown observation to model data
    model_test_data <- rbind(model_test_data, data_to_ID[fz,])
    
    test_mean <- vector("numeric", length(groups_for_test)) 
    ratio_test_within <- vector("numeric", length(groups_for_test))
    ratio_test_between <- vector("numeric", length(groups_for_test))
    
    #####################################################################################
    #####THIS FOR LOOP 'ff' iterates per possible group ID for the one unknown observation defined above ('fz')
    ########################################################################################
    
    for (ff in 1:length(groups_for_test)) {
      
      model_test_data_iteration <- model_test_data
      model_test_data_iteration[which(model_test_data_iteration$gr == "unknown", arr.ind=TRUE),1] <- groups_for_test[ff]
      
      f = file()
      sink(file=f) ## silence upcoming output using anonymous file connection
      
      test_CVA <- CVA(model_test_data[,-1:-2],model_test_data_iteration$gr)
      
      sink() ## undo silencing
      close(f)
      
      #RATIO TEST 
      
      ratio_test_CV_scores <- data.frame(test_CVA[["CVscores"]])
      ratio_test_CV_scores_ind_fordist <- cbind(model_test_data_iteration$gr,data.frame(ratio_test_CV_scores, row.names = row.names(model_test_data_iteration)))
      
      row.names(ratio_test_CV_scores_ind_fordist) <- paste(row.names(ratio_test_CV_scores_ind_fordist),ratio_test_CV_scores_ind_fordist[,1], sep = "&")
      
      ratio_test_ind_Mdist_matrix <- dist(ratio_test_CV_scores_ind_fordist) 
      ratio_test_ind_Mdist_pairwise <- matrixConvert(ratio_test_ind_Mdist_matrix)
      ratio_test_strSPLITsp1 <- str_split_fixed(ratio_test_ind_Mdist_pairwise$sp1, "&", 2)
      ratio_test_strSPLITsp2 <- str_split_fixed(ratio_test_ind_Mdist_pairwise$sp2, "&", 2)
      
      ratio_test_ind_Mdist_pairwise <- data.frame(ratio_test_strSPLITsp1[,2],ratio_test_strSPLITsp2[,2],ratio_test_ind_Mdist_pairwise$dist)
      names(ratio_test_ind_Mdist_pairwise)[1] <- "sp1"
      names(ratio_test_ind_Mdist_pairwise)[2] <- "sp2"
      names(ratio_test_ind_Mdist_pairwise)[3] <- "dist"
      
      ratio_test_within_grp_comps <- subset(ratio_test_ind_Mdist_pairwise, sp1 == sp2)
      ratio_test_between_grp_comps <- subset(ratio_test_ind_Mdist_pairwise, sp1 != sp2)
      
      ratio_test_mean_within_grp_dist <- mean(ratio_test_within_grp_comps$dist)
      ratio_test_mean_between_grp_dist <- mean(ratio_test_between_grp_comps$dist)
      
      ratio_test_within[ff] <- ratio_test_mean_within_grp_dist
      ratio_test_between[ff] <- ratio_test_mean_between_grp_dist
      test_mean[ff] <- ratio_test_mean_between_grp_dist/ratio_test_mean_within_grp_dist
      
      
    }
    
    test_mean_ordered <- data.frame(test_mean)
    test_mean_ordered_num <- data.frame(as.numeric(row.names(test_mean_ordered)), test_mean_ordered)
    test_mean_ordered <- arrange(test_mean_ordered_num, -test_mean)
    
    iter_of_true_ID <- test_mean_ordered[1,1]
    
    true_ID <- groups_for_test[iter_of_true_ID]
    
    result_IDs[fz] <- true_ID
    
    
  }
  
  data_to_ID_final <- data.frame(row.names(data_to_ID),result_IDs,data)
  names(data_to_ID_final)[1] <- "ind"
  names(data_to_ID_final)[2] <- "group"
  
  f = file()
  sink(file=f) ## silence upcoming output using anonymous file connection
  
  FINAL_CVA <- CVA(data_to_ID_final[,-1:-2],data_to_ID_final$group)
  
  sink() ## undo silencing
  close(f)
  
  FINAL_CVA_mean <- min(FINAL_CVA[["Dist"]][["GroupdistMaha"]])
  
  FINAL_IDS_DATA <- arrange(data_to_ID_final, ind)
  known_IDs <- arrange(known_IDs, ind)
  
  #calculate accurate IDs for each cross-validation
  accuracy <- data.frame(ifelse(known_IDs$gr==as.character(FINAL_IDS_DATA$group),"Yes","No"))
  names(accuracy)[1] <- "correct ID"
  
  accuracy_val[l] <- (100/nrow(accuracy))*sum(accuracy$`correct ID` == "Yes")
  
  #repeat 60/40% cross_validation 'l' times
  
  info_f <- sprintf("%f%% cross-validations done",ceiling(100/cross_val_no*l))
  setWinProgressBar(pb_f, ceiling(100/cross_val_no*l), label=info_f)
  
}
#calcuate mean and 95% C.I. for all cross-validations to evaluate model
Boot_CI <- data.frame(CI(accuracy_val,0.95))

write.table(Boot_CI, file = "Bootstrap CI.txt", sep="\t")


