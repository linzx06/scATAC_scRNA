# scATAC_scRNA
The R package will be released soon. Source code for the manuscript "Model-based approach to the joint analysis of single-cell data on chromatin accessibility and gene expression".

# run it on test data
library(mixtools)

library(label.switching)

load("test_data.RData")

source("gibbs_functions.R")

result_test <- getClusterGibbs(data_acc=test_data$data_acc, data_exp=test_data$data_exp, 
                               overlap_seq_acc=test_data$overlap_seq_acc, overlap_seq_exp=test_data$overlap_seq_exp, 
                               nCluster=2, niter=1000)
                               
# compare the result with the true cell type label
table((test_data$acc_true_cluster)[which(result_test$cluster_acc==1)]) 

table((test_data$acc_true_cluster)[which(result_test$cluster_acc==2)])

table((test_data$exp_true_cluster)[which(result_test$cluster_exp==1)]) 

table((test_data$exp_true_cluster)[which(result_test$cluster_exp==2)])

