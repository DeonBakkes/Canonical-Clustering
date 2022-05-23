require(dplyr)
require(Morpho)
require(microbenchmark)
require(usedist)
require(stringr)
require(ggplot2)
require(svDialogs)
require(otuSummary)
require(fs)

####  ####  ####  ####
# # This scripts test a fundamental assumption about the data - that the RATIO TEST distances (between and within group ~ Mahalanobis)
# # increase as the dataset is progressively classified more accurately - if met (skewed curve increasing to the right), this method is reliable
# # This phenomenon is related to how the CVA 'chooses' data features (axes) at oblique angles to each other (not orthogonal like PCA) in order 
# # to maximise between-group distance. Heterogenous groups (containing some incorrect IDs) will limit the effectiveness of distinct data 
# # features standing out to be detected by the CVA because feature variation is bi-modal (two peaks) - forcing different groups to appear 
# # more similar. Homogenous groups will not have these 'clashes' in feature distinctiveness as there is no competition (bi-modal) for the 
# # most distinct feature (given the observations are the same type ['species'] and variation between features is uni-modal (normal distribution)
####  ####  ####  ####


#reads MODEL data .txt file exported from ID unknown script
model_data <- read.table(choose.files(caption = "Choose data file (.txt) with MODEL", multi = FALSE),sep="\t",header = T)
names(model_data)[1] <- "ind"
names(model_data)[2] <- "gr"
model_data <- data.frame(row.names = model_data$ind,model_data)

model_data_alpha <- model_data

#randomise group labels to make most of them 'wrong'
model_data_wrong <- model_data 
model_data_wrong$gr <- sample(model_data$gr)

#create empty tables for final scores
accuracy_val <- vector("numeric", nrow(model_data))

#####################################################
#This for loop runs iterations on 'l' corrections of wrong labels
#####################################################

test_mean <- vector("numeric", nrow(model_data)) 
ratio_test_within <- vector("numeric", nrow(model_data))
ratio_test_between <- vector("numeric", nrow(model_data))

pb_f <- winProgressBar(title="Progress bar - labels corrected", label="% tests done", min=0, max=100, initial=0)

for (l in 1:nrow(model_data)) {
  
  wrong_labels <- model_data_wrong[1:nrow(model_data_wrong),]
  
  MODEL_wrong_labels <- data.frame(model_data_wrong[,-1:-2])
  
  f = file()
  sink(file=f) ## silence upcoming output using anonymous file connection
  
  model_CVA <- CVA(MODEL_wrong_labels,model_data_wrong$gr)
  
  sink() ## undo silencing
  close(f)
  
  #RATIO TEST
  
  ratio_test_CV_scores <- data.frame(model_CVA[["CVscores"]])
  ratio_test_CV_scores_ind_fordist <- cbind(model_data_wrong$gr,data.frame(ratio_test_CV_scores, row.names = row.names(model_data_wrong)))
  
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
  
  ratio_test_within[l] <- ratio_test_mean_within_grp_dist
  ratio_test_between[l] <- ratio_test_mean_between_grp_dist
  test_mean[l] <- ratio_test_mean_between_grp_dist/ratio_test_mean_within_grp_dist
  
  #then add +1 correct ID to wrong IDs dataset
  
  model_data_wrong$gr[l] <- model_data_alpha$gr[l]
  
  info_f <- sprintf("%f%% labels corrected",ceiling(100/nrow(model_data)*l))
  setWinProgressBar(pb_f, ceiling(100/nrow(model_data)*l), label=info_f)
  
}

#plot the RATIO TEST for each iteration of corrected labels
test_mean_table <- data.frame(test_mean)

test_mean_table <- data.frame(as.numeric(row.names(test_mean_table)), test_mean_table$test_mean)
names(test_mean_table)[1] <- "accurate labels added"
names(test_mean_table)[2] <- "mean Mdist"

pdf("Mdist-accuracy regression.pdf", width = 10, height = 10)

plot(test_mean_table$`accurate labels added`,test_mean_table$`mean Mdist`,
     xlab=paste("accurate labels added (x-y correlation below)"),
     ylab=paste("mean Mdist"), main = "Label accuracy by mean Mdist", sub = cor(test_mean_table$`accurate labels added`,test_mean_table$`mean Mdist`))

dev.off()

