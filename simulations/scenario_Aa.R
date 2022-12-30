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
NN_input[[1]] <- c(50,50,50)
NN_input[[2]] <- c(50,50,50)
PP_input[[1]] <- c(50,50)
PP_input[[2]] <- c(50,50)
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
X_simulated[[2]] <- F_simulated[[1]] %*% S_simulated[[2]] %*% t(G_simulated[[2]])

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
nIter <- 2000
test_phi <- c(0,0.1,0.3,0.5,0.8,1.0,2.0,5.0,10.0,20.0,50.0)
test_phi <- c(0, 0.1, 0.5, 1.0, 10.0)
tuningAccuracy <- vector("list", length = length(test_phi))
tuningAdjRandValue <- vector("list", length = length(test_phi))
tuningNMIvalue <- vector("list", length = length(test_phi))
tuningError <- vector("list", length = length(test_phi))
for (k in 1:length(test_phi)){
  phi[1,2] <- test_phi[k]
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
  if (test_phi[k] == 0){
    finalAccuracy_single <- tuningAccuracy[[k]]
    finalAdjRandValue_single <- tuningAdjRandValue[[k]]
    finalNMIvalue_single <- tuningNMIvalue[[k]]
    finalError_single <- tuningError[[k]]
  }
}
# Run once more with optimised error
weights_metrics <- c(0.9,0.1) # First to error, second to accuracy
cv_metrics <- cv_eval_metric(accList = tuningAccuracy,errList = tuningError, testValues = test_phi, weights = weights_metrics)
phi[1,2] <- cv_metrics$minValue
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
selectedTuning <- test_phi[min_tuning]
finalError_multiple <- X_simulated_NMTF$Error

# Save
save.image("scenAa_test.RData")
write.table(finalError_single, "finalError_scenAa_single.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalAdjRandValue_single, "finalAdjRandValue_scenAa_single.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalNMIvalue_single, "finalNMIvalue_scenAa_single.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalAccuracy_single, "finalAccuracy_scenAa_single.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(finalError_multiple, "finalError_scenAa_multiple.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalAdjRandValue_multiple, "finalAdjRandValue_scenAa_multiple.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalNMIvalue_multiple, "finalNMIvalue_scenAa_multiple.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(finalAccuracy_multiple, "finalAccuracy_scenAa_multiple.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(selectedTuning, "selectedTuning_scenAa.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

