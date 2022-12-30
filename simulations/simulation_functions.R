#' In this script we will write functions to generate data to simulate
#' Restrictive Multi-NMTF (for which the functions are the file mentioned below)
#' #### 'rest_multi_NMTF_functions.R' ####
#' 
#' The idea is 
#' 1. generate 'true' matrices F^v, S^v, G^v
#' 2. let X^v = F^v % S^v % t(G^v)
#' Then to test them:
#' 3. initialise Fhat^v, Shat^v, Ghat^v
#' 4. run function
#' 5. compute Xhat^v = Fhat^v % Shat^v % t(Ghat^v)
#' 6. compare with truth
#' 
#' NOTE: Store labels as well as data

# Directory
setwd("C:/Users/theod/OneDrive - Imperial College London/Documents/Imperial/PhD/Multi-Clustering/Restrictive Multi-NMTF")

# Libraries
library(MASS)
library(fossil)
library(data.table)

# Functions

simulate_biclusters <- function(nn, pp){
  #' N: Number of samples
  #' P: Number of features
  #' K: Number of sample-clusters (row-clusters)
  #' L: Number of feature-clusters (column-clusters)
  #' nn: Vector with each element correpsonding to the number of samples in each row-cluster
  #' pp: Vector with each element correpsonding to the number of samples in each col-cluster
  #' i.e. sum(nn) = N
  #'      sum(pp) = P

  N <- sum(nn)
  P <- sum(pp)
  K <- length(nn)
  L <- length(pp)
  ## Simulate F, G and S separately:
  # Define row-clusters and column-clusters
  row_clusters <- c()
  for (i in 1:K){
    row_clusters <- c(row_clusters, rep(i,nn[i]))
  }
  col_clusters <- c()
  for (i in 1:L){
    col_clusters <- c(col_clusters, rep(i,pp[i]))
  }
  
  # Simulate F
  F_noise <- abs(mvrnorm(n = N, mu = seq(1, K,1), Sigma = diag(K)))
  F_clust <- matrix(0, nrow = nrow(F_noise), ncol = ncol(F_noise))
  s <- 1
  e <- 0
  for (j in 1:K){
    e <- e + nn[j]
    F_clust[s:e,j] <- 1 + mean(F_noise[s:e,])
    s <- e + 1
  }
  F_original <- abs(scale(F_noise) + F_clust)
  row_clusters_original <- apply(F_original,1,which.max)
  
  # Simulate G
  G_noise <- abs(mvrnorm(n = P, mu = seq(1, L,1), Sigma = diag(L)))
  G_clust <- matrix(0, nrow = nrow(G_noise), ncol = ncol(G_noise))
  s <- 1
  e <- 0
  for (j in 1:L){
    e <- e + pp[j]
    G_clust[s:e,j] <- 1 + mean(G_noise[s:e,])
    s <- e + 1
  }
  G_original <- abs(scale(G_noise) + G_clust)
  #table(apply(abs(G_original),1,which.max), col_clusters)
  col_clusters_original <- apply(G_original,1,which.max)
  
  # Simulate S
  S_original <- abs(mvrnorm(n = K, mu = rep(0,L), Sigma = diag(L)))
  
  # Compute X
  X_original <- F_original %*% S_original %*% t(G_original)
  
  return(list("X" = X_original,
              "F" = F_original,
              "S" = S_original,
              "G" = G_original,
              "row_clusters" = row_clusters,
              "row_clusters_original" = row_clusters_original,
              "col_clusters" = col_clusters,
              "col_clusters_original" = col_clusters_original))
}


simulate_multiClusters <- function(nn, pp){
  #' N: Number of samples
  #' P: Number of features
  #' K: Number of sample-clusters (row-clusters)
  #' L: Number of feature-clusters (column-clusters)
  #' nn: Vector with each element correpsonding to the number of samples in each row-cluster
  #' pp: Vector with each element correpsonding to the number of samples in each col-cluster
  #' i.e. sum(nn) = N
  #'      sum(pp) = P
  #'      
  #' F -> Can be represented as a matrix of squares (KxK):
  #'      with K squares in rows and K squares in columns
  #'      One square per clustering
  #' G -> Same as F
  #' S -> Scale multivariate normally distributed matrix (?)
  
  N <- sum(nn)
  P <- sum(pp)
  K <- length(nn)
  L <- length(pp)
  ## Simulate F, G and S separately:
  # Define row-clusters and column-clusters
  row_clusters <- c()
  for (i in 1:K){
    row_clusters <- c(row_clusters, rep(i,nn[i]))
  }
  col_clusters <- c()
  for (i in 1:L){
    col_clusters <- c(col_clusters, rep(i,pp[i]))
  }
  
  # Simulate F
  # Generate data, cluster-by-cluster
  F_simulated <- matrix(0, nrow = N, ncol = K)
  row_end <- 0
  mean_matrix <- diag(rep(10,K))
  for (i in 1:K){
    row_start <- row_end + 1
    row_end <- row_end + nn[i]
    col_end <- 0
    for (j in 1:K){
      col_start <- col_end + 1
      col_end  <- col_end + pp[j]
      F_temp <- rnorm(n = nn[i], mean = mean_matrix[i,j], sd = 1)
      F_simulated[row_start:row_end, j] <- F_temp
    }
  }
  F_noise <- abs(mvrnorm(n = N, mu = rep(0, K), Sigma = diag(K)))
  F_original <- abs(scale(F_noise + F_simulated))
  row_clusters_original <- apply(F_original,1,which.max)
  
  # Simulate G
  G_simulated <- matrix(0, nrow = P, ncol = L)
  row_end <- 0
  mean_matrix <- diag(rep(10,L))
  for (i in 1:L){
    row_start <- row_end + 1
    row_end <- row_end + pp[i]
    col_end <- 0
    for (j in 1:L){
      col_start <- col_end + 1
      col_end  <- col_end + pp[j]
      G_temp <- rnorm(n = pp[i], mean = mean_matrix[i,j], sd = 1)
      G_simulated[row_start:row_end, j] <- G_temp
    }
  }
  G_noise <- abs(mvrnorm(n = P, mu = rep(0, L), Sigma = diag(L)))
  G_original <- abs(scale(G_noise + G_simulated))
  col_clusters_original <- apply(G_original,1,which.max)
  
  # Simulate S
  S_original <- abs(mvrnorm(n = K, mu = rep(0,L), Sigma = diag(L)))
  
  # Compute X
  X_original <- F_original %*% S_original %*% t(G_original)
  
  return(list("X" = X_original,
              "F" = F_original,
              "S" = S_original,
              "G" = G_original,
              "row_clusters" = row_clusters,
              "row_clusters_original" = row_clusters_original,
              "col_clusters" = col_clusters,
              "col_clusters_original" = col_clusters_original))
}

####################################################################################################
########  Test  ###########################################################################
####################################################################################################



n_v <- 3 # Number of views
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
NN_input[[1]] <- c(100,100,100,100)
NN_input[[2]] <- c(100,100,100,100,100)
NN_input[[3]] <- c(100,100,100)
PP_input[[1]] <- c(50,50,50,50,50)
PP_input[[2]] <- c(50,50,50,50,50)
PP_input[[3]] <- c(50,50,50,50,50)
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

# Assume that G^1 = G^2 = G^3
# i.e.:
X_simulated[[2]] <- F_simulated[[1]] %*% S_simulated[[2]] %*% t(G_simulated[[2]])
#X_simulated[[3]] <- F_simulated[[1]] %*% S_simulated[[3]] %*% t(G_simulated[[3]])
X_simulated[[3]] <- F_simulated[[3]] %*% S_simulated[[3]] %*% t(G_simulated[[1]])
true_row_clusterings[[2]] <- true_row_clusterings[[1]]
true_col_clusterings[[3]] <- true_col_clusterings[[1]]
true_row_clusterings_original[[2]] <- true_row_clusterings_original[[1]]
true_col_clusterings_original[[3]] <- true_col_clusterings_original[[1]]
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
phi <- matrix(0, nrow = length(X_simulated), ncol = length(X_simulated))
xi <- matrix(0, nrow = length(X_simulated), ncol = length(X_simulated))
psi <- matrix(0, nrow = length(X_simulated), ncol = length(X_simulated))

# Start with general: R = all-ones
Rinput <- matrix(1, nrow = 3, ncol = length(X_simulated))
psi[1,2] <- 0.5
psi[1,3] <- 1
psi[2,3] <- 2
#xi[2,3] <- 1
Rinput[1,1] <- 0
Rinput[1,2] <- 0
#Rinput[1,3] <- 0
Rinput[3,2] <- 0
Rinput[3,3] <- 0
#Rinput[2,2] <- 0
#Rinput[2,3] <- 0
X_simulated_NMTF <- restMultiNMTF_run(Xinput = X_simulated,
                                      Finput = Finit,
                                      Sinput = Sinit,
                                      Ginput = Ginit,
                                      phi = phi,
                                      xi = xi,
                                      psi = psi,
                                      nIter = nIter)
plot(X_simulated_NMTF$Error, type = "o")
# Get clusterings
eval_measures_NMTF <- evaluate_simulation(X_nmtf = X_simulated_NMTF,
                                          true_row_clustering = true_row_clusterings,
                                          true_col_clustering = true_col_clusterings)
eval_measures_NMTF$accuracy

# Test tuning parameter
test_phi <- c(0,0.1,0.3,0.5,0.8,1.0,2.0,5.0,10.0,20.0,50.0)
tuningAccuracy <- vector("list", length = length(test_phi))
tuningAdjRandValue <- vector("list", length = length(test_phi))
tuningNMIvalue <- vector("list", length = length(test_phi))
tuningError <- matrix(0, nrow = length(test_phi), ncol = 1000)
for (k in 1:length(test_phi)){
  phi[1,2] <- test_phi[k]
  X_simulated_NMTF <- restrictiveMultiNMTF(X = X_simulated,
                                           R = Rinput,
                                           Finit = Finit,
                                           Sinit = Sinit,
                                           Ginit = Ginit,
                                           phi = phi,
                                           xi = xi,
                                           psi = psi,
                                           nIter = 1000)
  # Get clusterings
  row_clustering <- vector("list", length = length(X_simulated))
  column_clustering <- vector("list", length = length(X_simulated))
  for (i in 1:length(X_simulated)){
    row_clustering[[i]] <- apply(X_simulated_NMTF$Foutput[[i]], 1, which.max)
    column_clustering[[i]] <- apply(X_simulated_NMTF$Goutput[[i]], 1, which.max)
  }
  column_table <- vector("list", length = length(X_simulated))
  row_table <- vector("list", length = length(X_simulated))
  accuracy <- matrix(0, nrow = 2, ncol = length(X_simulated))
  adjRandValue <- matrix(0, nrow = 2, ncol = length(X_simulated))
  nmiValue <- matrix(0, nrow = 2, ncol = length(X_simulated))
  rownames(accuracy) <- c("Row-clustering", "Column-clustering")
  rownames(adjRandValue) <- c("Row-clustering", "Column-clustering")
  rownames(nmiValue) <- c("Row-clustering", "Column-clustering")
  #column_clustering[[2]][203:223] <- 4
  for (i in 1:length(X_simulated)){
    row_table[[i]] <- table(row_clustering[[i]], true_row_clusterings[[i]])
    if (nrow(row_table[[i]]) < ncol(row_table[[i]])){
      for (g in 1:(ncol(row_table[[i]]) - nrow(row_table[[i]]))){
        row_table[[i]] <- rbind(row_table[[i]], rep(0,ncol(row_table[[i]])))
      }
    }else if (nrow(row_table[[i]]) > ncol(row_table[[i]])){
      for (g in 1:(nrow(row_table[[i]]) - ncol(row_table[[i]]))){
        row_table[[i]] <- cbind(row_table[[i]], rep(0,nrow(row_table[[i]])))
      }
    }
    row_table[[i]] <- fixTablingClustering(row_table[[i]])
    column_table[[i]] <- table(column_clustering[[i]], true_col_clusterings[[i]])
    if (nrow(column_table[[i]]) < ncol(column_table[[i]])){
      for (g in 1:(ncol(column_table[[i]]) - nrow(column_table[[i]]))){
        column_table[[i]] <- rbind(column_table[[i]], rep(0,ncol(column_table[[i]])))
      }
    }else if (nrow(column_table[[i]]) > ncol(column_table[[i]])){
      for (g in 1:(nrow(column_table[[i]]) - ncol(column_table))){
        column_table[[i]] <- cbind(column_table[[i]], rep(0,nrow(column_table[[i]])))
      }
    }
    column_table[[i]] <- fixTablingClustering(column_table[[i]])
    accuracy[1,i] <- accuracyTable(row_table[[i]])
    accuracy[2,i] <- accuracyTable(column_table[[i]])
    adjRandValue[1,i] <- adjustedRandIndex(row_clustering[[i]], true_row_clusterings[[i]])
    adjRandValue[2,i] <- adjustedRandIndex(column_clustering[[i]], true_col_clusterings[[i]])
    nmiValue[1,i] <- NMI(row_clustering[[i]], true_row_clusterings[[i]])
    nmiValue[2,i] <- NMI(column_clustering[[i]], true_col_clusterings[[i]])
  }
  tuningAccuracy[[k]] <- accuracy
  tuningAdjRandValue[[k]] <- adjRandValue
  tuningNMIvalue[[k]] <- nmiValue
  tuningError[k,] <- X_simulated_NMTF$Error
}
