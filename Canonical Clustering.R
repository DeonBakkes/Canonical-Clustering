require(Morpho)
require(microbenchmark)
require(usedist)
require(stringr)
require(dplyr)
require(ggplot2)
require(svDialogs)
require(otuSummary)
require(fs)

#reads data .txt file in the vector form: observation ID (1st column) --> variables

data <- read.table(choose.files(caption = "Choose data file (.txt)", multi = FALSE),sep="\t",header = T, row.names = 1)

#Sorts data, observation order is important, so we start with an alphabetical system and stick with it

sorted_data <- data
sorted_data <- data.frame(row.names(data),data)
sorted_data <- arrange(sorted_data, row.names(sorted_data))
sorted_data <- data.frame(row.names = sorted_data[,1], sorted_data[,-1])

######################################################
### run elbow plot to estimate best a priori group number (k) - best to do this a few times
######################################################

kmean_withinss <- function(k) {
  cluster <- kmeans(sorted_data, k)
  return (cluster$tot.withinss)
}

# Set maximum possible groups (more than 10 is a bit crazy for k-means) 
max_k <-10

# Run algorithm over a range of k 
wss <- sapply(2:max_k, kmean_withinss)
elbow <-data.frame(2:max_k, wss)

# Plot the graph with gglop
ggplot(elbow, aes(x = X2.max_k, y = wss)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = seq(1, 20, by = 1))

ggsave("elbow plot.pdf",scale = 1)

file_show(path(getwd(), "elbow plot.pdf"))

#input info from elbow plot
iteration_groups <- dlg_input("Please specify the number of apriori groups (k) (check 'elbow plot.pdf' saved to working directory): ", Sys.info()["iteration_groups"])$res
iteration_groups <- as.numeric(iteration_groups)

#######################################################
### Begin with initial K-means search
#######################################################

groups_vector <- sample(runif(n = 50, min = 0, max = iteration_groups), size = nrow(sorted_data), replace = TRUE)
groups_vector <- ceiling(groups_vector)                          

k_n <- dlg_input("Please specify the number of iterations for model search (10000 or more is usually best): ", Sys.info()["k_n"])$res
k_n <- as.numeric(k_n)
#set up vectors to fill with for loop
Group_Maha <- vector('numeric',k_n)
Kmeans_IDs_stored <- list()

pb_k <- winProgressBar(title="Progress bar - model search", label="% iterations done", min=0, max=100, initial=0)

###########################################################
###this for loop 'explores' possibilities in the k-means space
#########################################################

for (k in 1:k_n) {
  set.seed(get_nanotime()/1000000) 
  kmeans <- kmeans(sorted_data,iteration_groups, iter.max = 1000)
  
  Kmeans_IDs <- data.frame(row.names = row.names(sorted_data), kmeans[["cluster"]])
  names(Kmeans_IDs)[1] <- "gr"
  
  #CVA doesn't work when there are less than 3 samples in a group, so we have to exclude these from label randomisations
  #K-means normally doesn't produce <3 observation groups, but you never know...
  
  SINGLETONchecktable <- data.frame(table(Kmeans_IDs$gr))
  SINGLETONchecktest <-  SINGLETONchecktable$Freq
  
  while(any(SINGLETONchecktest < 2)) {
    
    kmeans <- kmeans(sorted_data,iteration_groups, iter.max = 1000)
    
    Kmeans_IDs <- data.frame(row.names = row.names(sorted_data), kmeans[["cluster"]])
    names(Kmeans_IDs)[1] <- "gr"
    
    SINGLETONchecktable <- data.frame(table(Kmeans_IDs$gr))
    SINGLETONchecktest <-  SINGLETONchecktable$Freq
    
  }
  
  Kmeans_IDs_for_CVA <- data.frame(row.names = row.names(Kmeans_IDs), Kmeans_IDs$gr,sorted_data)
  
  f = file()
  sink(file=f) ## silence upcoming output using anonymous file connection
  
  #run a CVA based on the k-means result
  Kmeans_IDs_CVA <- CVA(Kmeans_IDs_for_CVA[,-1],Kmeans_IDs_for_CVA[,1])
  
  sink() ## undo silencing
  close(f)  
  
  #calculate Mahalanobis distances ('spread') between groups (higher tends to have less misclassifications when used in an iterative context 
  #see 'test assumption' script)
  Kmeans_calc_spread <- min(Kmeans_IDs_CVA[["Dist"]][["GroupdistMaha"]])
  
  Kmeans_calc <- Kmeans_calc_spread
  
  #repeat this k times to comprehensively 'search' all k-means space possibilities
  Group_Maha[k] <- Kmeans_calc
  
  Kmeans_IDs_stored[[k]] <- Kmeans_IDs_for_CVA[,1]
  
  info_k <- sprintf("%f%% steps done",ceiling(100/k_n*k))
  setWinProgressBar(pb_k, ceiling(100/k_n*k), label=info_k)
}

#compile Mahalanobis distances from k-means searches and select the highest! which has max liklihood of being closer to the truth 
#(see 'fundamental assumption')

fit_table <- data.frame(Group_Maha)
fit_table_num <- data.frame(as.numeric(row.names(fit_table)),Group_Maha)
fit_table <- arrange(fit_table_num, -Group_Maha)
max_likelihood <- fit_table[1,2]
max_likelihood_iteration <- fit_table[1,1]
max_likelihood_IDs <- Kmeans_IDs_stored[[max_likelihood_iteration]]
best_Kmeans <- data.frame(max_likelihood_IDs,max_likelihood)
best_Kmeans_IDs <- data.frame(row.names = row.names(sorted_data), best_Kmeans$max_likelihood_IDs)
names(best_Kmeans_IDs)[1] <- "gr"

#run CVA on 'best' k-means result
step_IDs_for_CVA <- data.frame(best_Kmeans_IDs, sorted_data)

f = file()
sink(file=f) ## silence upcoming output using anonymous file connection

step_IDs_CVA <- CVA(step_IDs_for_CVA[,-1],step_IDs_for_CVA[,1])

sink() ## undo silencing
close(f) 

k_means_mean <- min(step_IDs_CVA[["Dist"]][["GroupdistMaha"]])

pdf("k-means trained.pdf", width = 10, height = 10)

plot(step_IDs_CVA$CVscores, col=step_IDs_for_CVA[,1], pch=as.numeric(step_IDs_for_CVA[,1]), typ="n",asp=1,
     xlab=paste("1st canonical axis", paste(round(step_IDs_CVA$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(step_IDs_CVA$Var[2,2],1),"%")), main = "k-means trained", sub = k_means_mean)
text(step_IDs_CVA$CVscores, as.character(row.names(sorted_data)), col=as.numeric(step_IDs_for_CVA[,1]), cex=.7)

legend("topleft", 
       legend = unique(step_IDs_for_CVA[,1]), 
       cex = 1, 
       text.col = unique(step_IDs_for_CVA[,1]), fill = unique(step_IDs_for_CVA[,1]), 
       horiz = T)

dev.off()

#finalise IDs from K-means search and compile with original data
KMEANS_IDS_DATA <- data.frame(step_IDs_for_CVA)
KMEANS_IDS_DATA <- arrange(KMEANS_IDS_DATA, gr)

write.table(KMEANS_IDS_DATA, file = "Kmeans Grouped data.txt", sep="\t")

########################################################################################
######### NOW MAKE THE MODEL FROM THE BEST (most central) INDIVIDUALS FOR EACH GROUP
########################################################################################

CVA_ind <- data.frame(step_IDs_CVA[["CVscores"]])
CVA_ind_fordist <- cbind(step_IDs_for_CVA$gr,data.frame(CVA_ind, row.names = row.names(step_IDs_for_CVA)))

row.names(CVA_ind_fordist) <- paste(row.names(CVA_ind_fordist),CVA_ind_fordist[,1], sep = "&")

ind_Mdist_matrix <- dist(CVA_ind_fordist) 

ind_centroid_distance <- dist_to_centroids(ind_Mdist_matrix,as.factor(step_IDs_for_CVA$gr))

ind_centroid_distance_SPLIT <- str_split_fixed(ind_centroid_distance$Item, "&", 2) #splits strings apart
ind_centroid_distance <- data.frame(ind_centroid_distance_SPLIT,ind_centroid_distance[,-1], stringsAsFactors = FALSE) #places distance again
names(ind_centroid_distance)[2] <- "group"
names(ind_centroid_distance)[1] <- "ind"

ind_centroid_distance <- subset.data.frame(ind_centroid_distance, as.character(group) == as.character(CentroidGroup))

model_size_runs_mean <- vector('numeric',min(table(unlist(step_IDs_for_CVA$gr)))-2)

within <- vector('numeric',min(table(unlist(step_IDs_for_CVA$gr)))-2)
between <- vector('numeric',min(table(unlist(step_IDs_for_CVA$gr)))-2)
search_best_model <- vector('numeric',min(table(unlist(step_IDs_for_CVA$gr)))-2)
a_w_kmeans_groups <- data.frame(row.names(sorted_data),step_IDs_for_CVA$gr,sorted_data)

pb_m_s <- winProgressBar(title="Progress bar - model testing", label="% models done", min=0, max=100, initial=0)
pb_f <- winProgressBar(title="Progress bar - sample clustering", label="% samples done", min=0, max=100, initial=0)

################################################################################
############# This for loop searches for best model - iteratively includes one extra individual (all groups) and tests the between group distances
#################################################################################

#-2 to account for starting at three observations per group 
for (ms in 1:(min(table(unlist(step_IDs_for_CVA$gr)))-2)) {
  
  model_size <- ms+2 #+2 because we need at least 3 observations per groups for CVA to work
  
  #collect individual observations in each k-means group by distance from its centroid
  inds_for_model <- ind_centroid_distance %>% 
    mutate(gr_pair = paste(pmin(group, CentroidGroup), pmax(group, CentroidGroup), sep = ",")) %>%
    group_by(gr_pair) %>% 
    top_n(-model_size, CentroidDistance) 
  
  model_data <- subset(a_w_kmeans_groups,!is.na(match(row.names(a_w_kmeans_groups),inds_for_model$ind)))
  names(model_data)[1] <- "ind"
  names(model_data)[2] <- "gr"
  
  #bring in the unknown data
  data_to_ID <- subset(sorted_data,is.na(match(row.names(sorted_data),inds_for_model$ind)))
  data_to_ID <- data.frame(row.names(data_to_ID),rep("unknown",nrow(data_to_ID)),data_to_ID)
  names(data_to_ID)[1] <- "ind"
  names(data_to_ID)[2] <- "gr"
  
  result_IDs <- vector('numeric',length(data_to_ID$ind))
  
  groups_for_test <- unique(model_data$gr)
  
  #######################################################################
  #####THIS FOR LOOP 'fz' iterates per unknown observation
  #######################################################################
  
  for (fz in 1:length(data_to_ID$ind)) {
    
    model_test_data <- model_data
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
      model_test_data_iteration[which(model_test_data_iteration$gr == "unknown", arr.ind=TRUE),2] <- groups_for_test[ff]
      
      f = file()
      sink(file=f) ## silence upcoming output using anonymous file connection
      
      test_CVA <- CVA(model_test_data[,-1:-2],model_test_data_iteration$gr)
      
      sink() ## undo silencing
      close(f)
      
      #RATIO TEST OF WITHIN AND BETWEEN GROUP DISTANCES ('the test statistic' to evaluate most likely group ID for that observation)
      
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
    
    #evaluate which group got the best 'test statistic'
    test_mean_ordered <- data.frame(test_mean)
    test_mean_ordered_num <- data.frame(as.numeric(row.names(test_mean_ordered)), test_mean_ordered)
    test_mean_ordered <- arrange(test_mean_ordered_num, -test_mean)
    
    iter_of_true_ID <- test_mean_ordered[1,1]
    
    true_ID <- groups_for_test[iter_of_true_ID]
    
    result_IDs[fz] <- true_ID
    
    test_mean_ordered <- data.frame()
    model_test_data_iteration <- data.frame()
    model_test_data <- data.frame()
    
    #repeat for each unknown observation 
    
    info_f <- sprintf("%f%% steps done",ceiling(100/length(data_to_ID$ind)*fz))
    setWinProgressBar(pb_f, ceiling(100/length(data_to_ID$ind)*fz), label=info_f)
    
  }
  #evaluate the total set of tested IDs from current model size ('ms' iterator) with CVA
  data_to_ID_final <- data.frame(data_to_ID$ind,result_IDs)
  names(data_to_ID_final)[1] <- "ind"
  names(data_to_ID_final)[2] <- "group"
  
  FINAL_IDS <- rbind(data_to_ID_final[,1:2],inds_for_model[,1:2])
  FINAL_IDS$ind <- as.character(FINAL_IDS$ind)
  
  FINAL_IDS <- arrange(FINAL_IDS, ind)
  
  FINAL_IDS_DATA <- data.frame(FINAL_IDS,sorted_data)
  
  f = file()
  sink(file=f) ## silence upcoming output using anonymous file connection
  
  FINAL_CVA <- CVA(FINAL_IDS_DATA[,-1:-2],FINAL_IDS_DATA$group)
  
  sink() ## undo silencing
  close(f)
  
  #RATIO TEST
  
  search_best_model_CVA_ind <- data.frame(FINAL_CVA[["CVscores"]])
  search_best_model_CVA_ind_fordist <- cbind(FINAL_IDS_DATA$group,data.frame(search_best_model_CVA_ind, row.names = row.names(FINAL_IDS_DATA)))
  
  row.names(search_best_model_CVA_ind_fordist) <- paste(row.names(search_best_model_CVA_ind_fordist),search_best_model_CVA_ind_fordist[,1], sep = "&")
  
  search_best_model_ind_Mdist_matrix <- dist(search_best_model_CVA_ind_fordist) 
  search_best_model_ind_Mdist_pairwise <- matrixConvert(search_best_model_ind_Mdist_matrix)
  strSPLITsp1 <- str_split_fixed(search_best_model_ind_Mdist_pairwise$sp1, "&", 2)
  strSPLITsp2 <- str_split_fixed(search_best_model_ind_Mdist_pairwise$sp2, "&", 2)
  
  search_best_model_ind_Mdist_pairwise <- data.frame(strSPLITsp1[,2],strSPLITsp2[,2],search_best_model_ind_Mdist_pairwise$dist)
  names(search_best_model_ind_Mdist_pairwise)[1] <- "sp1"
  names(search_best_model_ind_Mdist_pairwise)[2] <- "sp2"
  names(search_best_model_ind_Mdist_pairwise)[3] <- "dist"
  
  within_grp_comps <- subset(search_best_model_ind_Mdist_pairwise, sp1 == sp2)
  between_grp_comps <- subset(search_best_model_ind_Mdist_pairwise, sp1 != sp2)
  
  mean_within_grp_dist <- mean(within_grp_comps$dist)
  mean_between_grp_dist <- mean(between_grp_comps$dist)
  
  #add test statistic for current model size to list for comparison later
  within[ms] <- mean_within_grp_dist
  between[ms] <- mean_between_grp_dist
  search_best_model[ms] <- mean_between_grp_dist/mean_within_grp_dist
  
  ####
  FINAL_CVA_mean <- min(FINAL_CVA[["Dist"]][["GroupdistMaha"]])
  
  model_size_runs_mean[ms] <- FINAL_CVA_mean
  
  #clean vectors to repeat test with model size + 1
  inds_for_model <- data.frame()
  model_data <- data.frame()
  FINAL_IDS <- data.frame()
  data_to_ID <- data.frame()
  
  info_m_s <- sprintf("%f%% steps done",ceiling(100/(min(table(unlist(step_IDs_for_CVA$gr)))-2)*ms))
  setWinProgressBar(pb_m_s, ceiling(100/(min(table(unlist(step_IDs_for_CVA$gr)))-2)*ms), label=info_m_s)
  
}

#search through all model size tests for the point where highest RATIO TEST distance arises and stabilises (informativeness saturation point)

#calculate iterative absolute differences between increased model sizes
d1 <- diff(search_best_model)
abs_d1 <- abs(d1)

#select iterations with no change (shows stability)
stability <- which(abs_d1 / max(abs_d1) < min(abs_d1[abs_d1 > 0]))

#select iteration with no change (stability) AND which has highest RATIO TEST distance
search_best_model_stability <- data.frame(stability, search_best_model[stability])
saturation_point <- 3+(search_best_model_stability$stability[which.max(search_best_model_stability$search_best_model.stability.)])
#+3 because we skip single and two-point datagroups and diffs substract 1

#if no saturation point, then just select highest
isEmpty <- function(x) {
  return(length(x)==0)
}


if (isEmpty(saturation_point)) {
  
  saturation_point <- 2+which.max(search_best_model)
  
}


pdf("model saturation curve.pdf", width = 10, height = 10)

plot(3:(length(search_best_model)+2),search_best_model, main = "Model saturation", sub = saturation_point)
abline(v=saturation_point)

dev.off()

#########################################################################################################
#####Use best model size based on highest informativeness (saturation point) to finally ID observations
########### this is the same as above fz and ff loops above (unknowns & group iterations), but we now know which model is best for this data
########################################################################################################

best_model_size <- saturation_point
best_model_size <- as.numeric(best_model_size)

inds_for_model <- ind_centroid_distance %>% 
  mutate(gr_pair = paste(pmin(group, CentroidGroup), pmax(group, CentroidGroup), sep = ",")) %>%
  group_by(gr_pair) %>% 
  top_n(-best_model_size, CentroidDistance) 

inds_selected_for_model <- data.frame(inds_for_model$ind,inds_for_model$group)

write.table(inds_selected_for_model, file = "inds selected for model.txt", sep="\t")

model_data <- subset(a_w_kmeans_groups,!is.na(match(row.names(a_w_kmeans_groups),inds_for_model$ind)))
names(model_data)[1] <- "ind"
names(model_data)[2] <- "gr"

data_to_ID <- subset(sorted_data,is.na(match(row.names(sorted_data),inds_for_model$ind)))
data_to_ID <- data.frame(row.names(data_to_ID),rep("unknown",nrow(data_to_ID)),data_to_ID)
names(data_to_ID)[1] <- "ind"
names(data_to_ID)[2] <- "gr"


f = file()
sink(file=f) ## silence upcoming output using anonymous file connection

model_CVA <- CVA(model_data[,-1:-2],model_data$gr)

sink() ## undo silencing
close(f)

model_mean <- min(model_CVA[["Dist"]][["GroupdistMaha"]])

pdf("model CVA.pdf", width = 10, height = 10)

plot(model_CVA$CVscores, col=model_data$gr, pch=as.numeric(model_data$gr), typ="n",asp=1,
     xlab=paste("1st canonical axis", paste(round(model_CVA$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(model_CVA$Var[2,2],1),"%")), main = "model CVA", sub = model_mean)
text(model_CVA$CVscores, as.character(model_data$ind), col=as.numeric(model_data$gr), cex=.7)

dev.off()

result_IDs <- vector('numeric',length(data_to_ID$ind))

groups_for_test <- unique(model_data$gr)

#####FOR LOOP WILL START HERE use 'fz' as iterator

for (fz in 1:length(data_to_ID$ind)) {
  
  model_test_data <- model_data
  
  model_test_data <- rbind(model_test_data, data_to_ID[fz,])
  
  test_mean <- vector("numeric", length(groups_for_test)) 
  ratio_test_within <- vector("numeric", length(groups_for_test))
  ratio_test_between <- vector("numeric", length(groups_for_test))
  
  for (ff in 1:length(groups_for_test)) {
    
    model_test_data_iteration <- model_test_data
    model_test_data_iteration[which(model_test_data_iteration$gr == "unknown", arr.ind=TRUE),2] <- groups_for_test[ff]
    
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
  
  
  info_f <- sprintf("%f%% steps done",ceiling(100/length(data_to_ID$ind)*fz))
  setWinProgressBar(pb_f, ceiling(100/length(data_to_ID$ind)*fz), label=info_f)
  
}

data_to_ID_final <- data.frame(data_to_ID$ind,result_IDs)
names(data_to_ID_final)[1] <- "ind"
names(data_to_ID_final)[2] <- "group"

FINAL_IDS <- rbind(data_to_ID_final[,1:2],inds_for_model[,1:2])
FINAL_IDS$ind <- as.character(FINAL_IDS$ind)

FINAL_IDS <- arrange(FINAL_IDS, ind)

FINAL_IDS_DATA <- data.frame(FINAL_IDS,sorted_data)

f = file()
sink(file=f) ## silence upcoming output using anonymous file connection

FINAL_CVA <- CVA(FINAL_IDS_DATA[,-1:-2],FINAL_IDS_DATA$group)

sink() ## undo silencing
close(f)

FINAL_CVA_mean <- min(FINAL_CVA[["Dist"]][["GroupdistMaha"]])


pdf("Machine Learning Clusters.pdf", width = 10, height = 10)

plot(FINAL_CVA$CVscores, col=FINAL_IDS_DATA$group, pch=as.numeric(FINAL_IDS_DATA$group), typ="n",asp=1,
     xlab=paste("1st canonical axis", paste(round(FINAL_CVA$Var[1,2],1),"%")),
     ylab=paste("2nd canonical axis", paste(round(FINAL_CVA$Var[2,2],1),"%")), main = "Machine Learning Clusters", sub = FINAL_CVA_mean)
text(FINAL_CVA$CVscores, as.character(FINAL_IDS_DATA$ind), col=as.numeric(FINAL_IDS_DATA$group), cex=.7)

legend("topleft", 
       legend = unique(FINAL_IDS_DATA$group), 
       cex = 1, 
       text.col = unique(FINAL_IDS_DATA$group), fill = unique(FINAL_IDS_DATA$group), 
       horiz = T)

dev.off()

mean_per_var_per_group <- aggregate(FINAL_IDS_DATA[, -1:-2], list(FINAL_IDS_DATA$group), mean)
names(mean_per_var_per_group)[1] <- "group"

write.table(mean_per_var_per_group, file = "means per variable per group.txt", sep="\t")
write.table(model_data, file = "MODEL data.txt", sep="\t")

FINAL_IDS_DATA <- arrange(FINAL_IDS_DATA, group)

write.table(FINAL_IDS_DATA, file = "Grouped data.txt", sep="\t")

file_show(path(getwd(), "Machine Learning Clusters.pdf"))

