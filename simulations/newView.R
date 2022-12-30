#' In this script we will test the functions on simulations found in the script:
#' #### 'rest_multi_NMTF_functions.R' ####
#'
#' In this simulation, we will assume that a new data-view is introduced, in which the samples remain the same, 
#' but the features are not previously seen. However, the features lie on the exact same clusters as other data-views
#' 
#' The aim of this simulation is to assess whether through restrictive Mutli-NMTF, one can cluster unseen features
#' based on the information obtained from analysing other data-views from the same samples
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
PP_input[[1]] <- c(60,40)
PP_input[[2]] <- c(40,60)
# Generate data
for (v in 1:n_v){
  nn_input <- NN_input[[v]]
  pp_input <- PP_input[[v]]
  temp <- simulate_biclusters(nn = nn_input,
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
X_simulated[[2]] <- F_simulated[[1]] %*% S_simulated[[1]] %*% t(G_simulated[[2]])
true_row_clusterings[[2]] <- true_row_clusterings[[1]]
true_row_clusterings_original[[2]] <- true_row_clusterings_original[[1]]
# Simulate new unseen data-view
temp <- simulate_biclusters(nn = c(50,50,50),
                            pp = c(50,50))
X_test<- temp$X
F_test<- temp$F
S_test<- temp$S
G_test<- temp$G
true_row_clusterings_test <- temp$row_clusters
true_col_clusterings_test <- temp$col_clusters
true_row_clusterings_original_test <- temp$row_clusters_original
true_col_clusterings_original_test <- temp$col_clusters_original
X_test <- F_simulated[[1]] %*% S_simulated[[1]] %*% t(G_test)


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
nIter <- 1000
test_phi <- c(0, 0.01, 0.1, 0.5, 0.8, 1.0, 1.2, 5.0, 10.0, 100.0)
#test_phi <- c(0,0.1,0.3,0.5,0.8,1.0,2.0,5.0,10.0,20.0,50.0)
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

weights_cv <- c(0.25,0.75)
cv_metrics <- cv_eval_metric(accList = tuningAccuracy,errList = tuningError, testValues = test_phi, weights = weights_cv)
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



# Compute the new G_hat
library(MASS)
G_hat <- t( ginv(X_simulated_NMTF$Soutput[[1]]) %*% ginv(X_simulated_NMTF$Foutput[[1]]) %*% X_test )
apply(G_hat, 1, which.max)
table(apply(G_hat, 1, which.max), apply(G_test, 1, which.max))
apply(abs(G_hat), 1, which.max)
G_hat <- t( t(X_simulated_NMTF$Soutput[[2]]) %*% t(X_simulated_NMTF$Foutput[[2]]) %*% X_test )
table(apply(G_hat, 1, which.max), apply(G_test, 1, which.max))
