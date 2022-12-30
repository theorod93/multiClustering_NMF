#' In this script we will test the functions on simulations found in the script:
#' #### 'rest_multi_NMTF_functions.R' ####
#'

#libraries
library(MASS)
library(data.table)
library(clue)
library(fossil)
library(mclust)
library(aricode)


# Initialisation

n_v <- 2 # Number of views
X_simulated <- vector("list", length = n_v)
F_simulated <- vector("list", length = n_v)
S_simulated <- vector("list", length = n_v)
G_simulated <- vector("list", length = n_v)

# Input parameters
true_row_clusterings <- vector("list", length = n_v)
true_row_clusterings_original <- vector("list", length = n_v)
true_col_clusterings <- vector("list", length = n_v)
true_col_clusterings_original <- vector("list", length = n_v)
NN_input <- vector("list", length = n_v)
PP_input <- vector("list", length = n_v)
NN_input[[1]] <- c(50,50)
NN_input[[2]] <- c(50,50)
PP_input[[1]] <- c(50,50,50)
PP_input[[2]] <- c(50,50,50)
# Generate data
for (v in 1:n_v){
  nn_input <- NN_input[[v]]
  pp_input <- PP_input[[v]]
  temp <- simulate_multiClusters(nn = nn_input,
                                 pp = pp_input)
  X_simulated[[v]]<- temp$X
  F_simulated[[v]]<- temp$F
  S_simulated[[v]]<- temp$S
  G_simulated[[v]]<- temp$G
  true_row_clusterings[[v]] <- temp$row_clusters
  true_col_clusterings[[v]] <- temp$col_clusters
  true_row_clusterings_original[[v]] <- temp$row_clusters_original
  true_col_clusterings_original[[v]] <- temp$col_clusters_original
}
# Restricions : F^(1) = F^(2)
X_simulated[[1]] <- F_simulated[[1]] %*% S_simulated[[1]] %*% t(G_simulated[[1]])
X_simulated[[2]] <- F_simulated[[2]] %*% S_simulated[[2]] %*% t(G_simulated[[1]])
true_col_clusterings[[2]] <- true_col_clusterings[[1]]
true_col_clusterings_original[[2]] <- true_col_clusterings_original[[1]]

# Apply Restrictive Multi-NMTF
# Initialisation
Finit <- vector("list", length = length(X_simulated))
Sinit <- vector("list", length = length(X_simulated))
Ginit <- vector("list", length = length(X_simulated))
for (i in 1:length(X_simulated)){
  ss <- svd(X_simulated[[i]])
  Finit[[i]] <- ss$u
  Sinit[[i]] <- diag(ss$d)
  Ginit[[i]] <- ss$v
}
KK <- lapply(NN_input, length)
LL <- lapply(PP_input, length)
for (i in 1:length(X_simulated)){
  L <- LL[[i]]
  K <- KK[[i]]
  Sinit[[i]] <- mvrnorm(n = K, mu = rep(0,L), Sigma = diag(L))
}
# Select top K for row-clusters and L for column clusters

for (i in 1:length(X_simulated)){
  L <- LL[[i]]
  K <- KK[[i]]
  Finit[[i]] <- abs(Finit[[i]][,1:K])
  Sinit[[i]] <- abs(Sinit[[i]][1:K,1:L])
  Ginit[[i]] <- abs(Ginit[[i]][,1:L])
}

# Tuning parameters
phi <- matrix(0, nrow = length(X_simulated), ncol = length(X_simulated))
xi <- matrix(0, nrow = length(X_simulated), ncol = length(X_simulated))
psi <- matrix(0, nrow = length(X_simulated), ncol = length(X_simulated))
nIter <- 5000
#test_psi <- c(0,0.1,0.3,0.5,0.8,1.0,2.0,5.0,10.0,20.0,50.0)
#test_psi <- c(0, 0.01, 0.1, 0.5, 0.8, 1.0, 1.2, 5.0, 10.0, 100.0)
test_psi <- c(0, 0.01, 0.1, 0.5, 0.8, 1.0, 1.2,10.0)
tuningAccuracy <- vector("list", length = length(test_psi))
tuningAdjRandValue <- vector("list", length = length(test_psi))
tuningNMIvalue <- vector("list", length = length(test_psi))
tuningError <- vector("list", length = length(test_psi))
for (k in 1:length(test_psi)){
  psi[1,2] <- test_psi[k]
  X_simulated_NMTF <- restMultiNMTF_run(Xinput = X_simulated,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)
  eval_measures_NMTF <- evaluate_simulation(X_nmtf = X_simulated_NMTF,
                                            true_row_clustering = true_row_clusterings,
                                            true_col_clustering = true_col_clusterings)
  tuningAccuracy[[k]] <- eval_measures_NMTF$accuracy
  tuningAdjRandValue[[k]] <- eval_measures_NMTF$ARI
  tuningNMIvalue[[k]] <- eval_measures_NMTF$NMI
  tuningError[[k]] <- X_simulated_NMTF$Error
  if (test_psi[k] == 0){
    finalAccuracy_single <- tuningAccuracy[[k]]
    finalAdjRandValue_single <- tuningAdjRandValue[[k]]
    finalNMIvalue_single <- tuningNMIvalue[[k]]
    finalError_single <- tuningError[[k]]
  }
}
# Run once more with optimised error
#weights_cv <- c(0.9,0.1)
#cv_metrics <- cv_eval_metric(accList = tuningAccuracy,errList = tuningError, testValues = test_psi, weights = weights_cv)
#psi[1,2] <- cv_metrics$minValue
min_tuning <- which.max(lapply(tuningAccuracy, function(x) mean(x[2,])))
psi[1,2] <- test_psi[min_tuning]
X_simulated_NMTF <- restMultiNMTF_run(Xinput = X_simulated,
                                      Finput = Finit,
                                      Sinput = Sinit,
                                      Ginput = Ginit,
                                      phi = phi,
                                      xi = xi,
                                      psi = psi,
                                      nIter = nIter)
eval_measures_NMTF <- evaluate_simulation(X_nmtf = X_simulated_NMTF,
                                          true_row_clustering = true_row_clusterings,
                                          true_col_clustering = true_col_clusterings)
finalAccuracy_multiple <- eval_measures_NMTF$accuracy
finalAdjRandValue_multiple <- eval_measures_NMTF$ARI
finalNMIvalue_multiple <- eval_measures_NMTF$NMI
selectedTuning <- cv_metrics$minValue
finalError_multiple <- X_simulated_NMTF$Error

# Store evaluation measures for each value in test_psi
plotRowAccuracy_Xa <- c()
plotRowAccuracy_Xb <- c()
plotColAccuracy_Xa <- c()
plotColAccuracy_Xb <- c()
for (k in 1:length(test_psi)){
  plotRowAccuracy_Xa <- c(plotRowAccuracy_Xa, tuningAccuracy[[k]][1,1])
  plotRowAccuracy_Xb <- c(plotRowAccuracy_Xb, tuningAccuracy[[k]][1,2])
  plotColAccuracy_Xa <- c(plotColAccuracy_Xa, tuningAccuracy[[k]][2,1])
  plotColAccuracy_Xb <- c(plotColAccuracy_Xb, tuningAccuracy[[k]][2,2])
}



## Implement HSIC:
# Based on Euclidean function as a kernel function, i.e. measureing distance between samples
hsic_values <- hsic_independence_test(lapply(X_simulated,t))
psi[1,2] <- hsic_values[1,2]
X_simulated_NMTF <- restMultiNMTF_run(Xinput = X_simulated,
                                      Finput = Finit,
                                      Sinput = Sinit,
                                      Ginput = Ginit,
                                      phi = phi,
                                      xi = xi,
                                      psi = psi,
                                      nIter = nIter)
eval_measures_NMTF <- evaluate_simulation(X_nmtf = X_simulated_NMTF,
                                          true_row_clustering = true_row_clusterings,
                                          true_col_clustering = true_col_clusterings)
finalAccuracy_hsic <- eval_measures_NMTF$accuracy
finalAdjRandValue_hsic <- eval_measures_NMTF$ARI
finalNMIvalue_hsic <- eval_measures_NMTF$NMI
finalError_hsic <- X_simulated_NMTF$Error


# Save
save.image("scenAc_test.RData")
write.table(finalError_single, "finalError_scenAc_single.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalAdjRandValue_single, "finalAdjRandValue_scenAc_single.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalNMIvalue_single, "finalNMIvalue_scenAc_single.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalAccuracy_single, "finalAccuracy_scenAc_single.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(finalError_multiple, "finalError_scenAc_multiple.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalAdjRandValue_multiple, "finalAdjRandValue_scenAc_multiple.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalNMIvalue_multiple, "finalNMIvalue_scenAc_multiple.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalAccuracy_multiple, "finalAccuracy_scenAc_multiple.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(selectedTuning, "selectedTuning_scenAc.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(finalError_hsic, "finalError_scenAc_hsic.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalAdjRandValue_hsic, "finalAdjRandValue_scenAc_hsic.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalNMIvalue_hsic, "finalNMIvalue_scenAc_hsic.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalAccuracy_hsic, "finalAccuracy_scenAc_hsic.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(t(plotRowAccuracy_Xa), "plotRowAccuracy_Xa_scenAc.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(t(plotRowAccuracy_Xb), "plotRowAccuracy_Xb_scenAc.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(t(plotColAccuracy_Xa), "plotColAccuracy_Xa_scenAc.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(t(plotColAccuracy_Xb), "plotColAccuracy_Xb_scenAc.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(hsic_values[1,2], "hsic_12_scenAc.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

