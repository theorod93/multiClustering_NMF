#' In this script, we will create simulations by generating block-diagonal matrices as X matrices
#' 
#' That is instead of generating separately F, S and G and then taking X = FSG^T
#' 


# Libraries
library(MASS)
library(Matrix)

## Scratch notes ## 
# Number of row-clusters and column-clusters
n_rowClusters <- 3
n_colClusters <- 3
# Introduce simulated data-views 
X_trial <- vector("list", length = 2)
true_row_clusterings <- vector("list", length = length(X_trial))
true_col_clusterings <- vector("list", length = length(X_trial))
# First data-view
# Dimensions of each block
# vector of length = n_rowClusters
row_dim <- c(100, 200, 300)
col_dim <- c(50, 100, 250)
# Create a list as input for block-diagonal data generation
inputList <- vector("list", length = max(length(row_dim), length(col_dim)))
inputList[[1]] <- mvrnorm(n =row_dim[1], mu = rep(5,col_dim[1]), Sigma = diag(col_dim[1]))
inputList[[2]] <- mvrnorm(n =row_dim[2], mu = rep(5,col_dim[2]), Sigma = diag(col_dim[2]))
inputList[[3]] <- mvrnorm(n =row_dim[3], mu = rep(5,col_dim[3]), Sigma = diag(col_dim[3]))
X_trial[[1]] <- as.matrix(bdiag(inputList))
X_noise <- mvrnorm(n = nrow(X_trial[[1]]), mu = rep(0, ncol(X_trial[[1]])), Sigma = diag(10,ncol(X_trial[[1]])))
X_trial[[1]] <- X_trial[[1]] + X_noise
true_row_clusterings[[1]] <- c(rep(1,100),
                               rep(2,200),
                               rep(3,300))
true_col_clusterings[[1]] <- c(rep(1,50),
                               rep(2,100),
                               rep(3,250))
# Repeat for second data-view (with different column dimensions)
# vector of length = n_rowClusters
row_dim <- c(100, 200, 300)
col_dim <- c(100, 100, 100)
# Create a list as input for block-diagonal data generation
inputList <- vector("list", length = max(length(row_dim), length(col_dim)))
inputList[[1]] <- mvrnorm(n =row_dim[1], mu = rep(6,col_dim[1]), Sigma = diag(col_dim[1]))
inputList[[2]] <- mvrnorm(n =row_dim[2], mu = rep(6,col_dim[2]), Sigma = diag(col_dim[2]))
inputList[[3]] <- mvrnorm(n =row_dim[3], mu = rep(6,col_dim[3]), Sigma = diag(col_dim[3]))
X_trial[[2]] <- as.matrix(bdiag(inputList))
X_noise <- mvrnorm(n = nrow(X_trial[[2]]), mu = rep(0, ncol(X_trial[[2]])), Sigma = diag(10,ncol(X_trial[[2]])))
X_trial[[2]] <- X_trial[[2]] + X_noise
true_row_clusterings[[2]] <- c(rep(1,100),
                               rep(2,200),
                               rep(3,300))
true_col_clusterings[[2]] <- c(rep(1,100),
                               rep(2,100),
                               rep(3,100))
#### Run Restrictive Multi-NMTF ####

# Initialisation
# Initialisation
Finit <- vector("list", length = length(X_trial))
Sinit <- vector("list", length = length(X_trial))
Ginit <- vector("list", length = length(X_trial))
for (i in 1:length(X_trial)){
  ss <- svd(X_trial[[i]])
  Finit[[i]] <- ss$u
  Sinit[[i]] <- diag(ss$d)
  Ginit[[i]] <- ss$v
}
KK <- c(3,3,3)
LL <- c(3,3,3)
for (i in 1:length(X_trial)){
  L <- LL[[i]]
  K <- KK[[i]]
  Sinit[[i]] <- mvrnorm(n = K, mu = rep(0,L), Sigma = diag(L))
}
# Select top K for row-clusters and L for column clusters

for (i in 1:length(X_trial)){
  L <- LL[[i]]
  K <- KK[[i]]
  Finit[[i]] <- abs(Finit[[i]][,1:K])
  Sinit[[i]] <- abs(Sinit[[i]][1:K,1:L])
  Ginit[[i]] <- abs(Ginit[[i]][,1:L])
}

# Tuning parameters
phi <- matrix(0, nrow = length(X_trial), ncol = length(X_trial))
xi <- matrix(0, nrow = length(X_trial), ncol = length(X_trial))
psi <- matrix(0, nrow = length(X_trial), ncol = length(X_trial))
nIter <- 1000
test_phi <- c(0, 0.1, 0.5, 1.0, 10.0)
tuningAccuracy <- vector("list", length = length(test_phi))
tuningAdjRandValue <- vector("list", length = length(test_phi))
tuningNMIvalue <- vector("list", length = length(test_phi))
tuningError <- vector("list", length = length(test_phi))
for (k in 1:length(test_phi)){
  phi[1,2] <- test_phi[k]
  X_trial_NMTF <- restMultiNMTF_run(Xinput = X_trial,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)
  eval_measures_NMTF <- evaluate_simulation(X_nmtf = X_trial_NMTF,
                                            true_row_clustering = true_row_clusterings,
                                            true_col_clustering = true_col_clusterings)
  tuningAccuracy[[k]] <- eval_measures_NMTF$accuracy
  tuningAdjRandValue[[k]] <- eval_measures_NMTF$ARI
  tuningNMIvalue[[k]] <- eval_measures_NMTF$NMI
  tuningError[[k]] <- X_trial_NMTF$Error
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
X_trial_NMTF <- restMultiNMTF_run(Xinput = X_trial,
                                      Finput = Finit,
                                      Sinput = Sinit,
                                      Ginput = Ginit,
                                      phi = phi,
                                      xi = xi,
                                      psi = psi,
                                      nIter = nIter)
eval_measures_NMTF <- evaluate_simulation(X_nmtf = X_trial_NMTF,
                                          true_row_clustering = true_row_clusterings,
                                          true_col_clustering = true_col_clusterings)
finalAccuracy_multiple <- eval_measures_NMTF$accuracy

# Based on Euclidean function as a kernel function, i.e. measureing distance between samples
hsic_values <- hsic_independence_test(X_trial)
phi[1,2] <- abs(scale(hsic_values))[1,2]
X_trial_NMTF <- restMultiNMTF_run(Xinput = X_trial,
                                      Finput = Finit,
                                      Sinput = Sinit,
                                      Ginput = Ginit,
                                      phi = phi,
                                      xi = xi,
                                      psi = psi,
                                      nIter = nIter)
eval_measures_NMTF_hsic <- evaluate_simulation(X_nmtf = X_trial_NMTF,
                                               true_row_clustering = true_row_clusterings,
                                               true_col_clustering = true_col_clusterings)
eval_measures_NMTF_hsic$accuracy


#### Create functions

generate_blockDiagonalMatrices_sameNclusters <- function(n_views,row_dim,col_dim){
  #' This is a function to simulate Block-diagonal matrices
  #' NOTE: In this function, we assume that the number of row-clusters and column-clusters is the same!!
  #' Input:
  #' n_views -> Number of data-views
  #' row_dim -> List of length = n_views: - Each element is a vector with the number of samples in each row-cluster.
  #'                                      - Of length = num of row-clusters
  #' col_dim -> List of length = n_views: - Each element is a vector with the number of samples in each col-cluster.
  #'                                      - Of length = num of col-clusters
  #' 

  # Introduce simulated data-views 
  X_trial <- vector("list", length = n_views)
  true_row_clusterings <- vector("list", length = n_views)
  true_col_clusterings <- vector("list", length = n_views)
  # Generate each view :
  for (v in 1:n_views){
    # Generate the separate blocks first
    inputList <- vector("list", length = length(row_dim[[v]]))
    mean_value <- round(runif(1,1,5))
    for (i in 1:length(row_dim[[v]])){
      inputList[[i]] <- mvrnorm(n =row_dim[[v]][i], mu = rep(mean_value,col_dim[[v]][i]), Sigma = diag(col_dim[[v]][i]))
      true_row_clusterings[[v]] <- c(true_row_clusterings[[v]], rep(i, row_dim[[v]][i]))
      true_col_clusterings[[v]] <- c(true_col_clusterings[[v]], rep(i, col_dim[[v]][i]))
    }
    X_trial_temp <- as.matrix(bdiag(inputList))
    sd_value_noise <- round(runif(1,7,15))
    mean_value_noise <- round(runif(1,7,15))
    X_noise <- mvrnorm(n = nrow(X_trial_temp),
                       mu = rep(mean_value_noise, ncol(X_trial_temp)), 
                       Sigma = diag(sd_value_noise,ncol(X_trial_temp)))
    X_trial[[v]] <- X_trial_temp + X_noise
  }
  return(list("X_trial" = X_trial,
              "true_row_clustering" = true_row_clusterings,
              "true_col_clustering" = true_col_clusterings))
}

generate_blockDiagonalMatrices_sameNclusters_uniform <- function(n_views,row_dim,col_dim){
  #' This is a function to simulate Block-diagonal matrices
  #' NOTE: In this function, we assume that the number of row-clusters and column-clusters is the same!!
  #' Input:
  #' n_views -> Number of data-views
  #' row_dim -> List of length = n_views: - Each element is a vector with the number of samples in each row-cluster.
  #'                                      - Of length = num of row-clusters
  #' col_dim -> List of length = n_views: - Each element is a vector with the number of samples in each col-cluster.
  #'                                      - Of length = num of col-clusters
  #' 
  
  # Introduce simulated data-views 
  X_trial <- vector("list", length = n_views)
  true_row_clusterings <- vector("list", length = n_views)
  true_col_clusterings <- vector("list", length = n_views)
  # Generate each view :
  for (v in 1:n_views){
    # Generate the separate blocks first
    inputList <- vector("list", length = length(row_dim[[v]]))
    mean_value <- round(runif(1,1,5))
    for (i in 1:length(row_dim[[v]])){
      inputList[[i]] <- matrix(runif(n = (row_dim[[v]][i]*col_dim[[v]][i])))
      #inputList[[i]] <- mvrnorm(n =row_dim[[v]][i], mu = rep(mean_value,col_dim[[v]][i]), Sigma = diag(col_dim[[v]][i]))
      true_row_clusterings[[v]] <- c(true_row_clusterings[[v]], rep(i, row_dim[[v]][i]))
      true_col_clusterings[[v]] <- c(true_col_clusterings[[v]], rep(i, col_dim[[v]][i]))
    }
    X_trial_temp <- as.matrix(bdiag(inputList))
    sd_value_noise <- round(runif(1,7,15))
    mean_value_noise <- round(runif(1,7,15))
    X_noise <- mvrnorm(n = nrow(X_trial_temp),
                       mu = rep(mean_value_noise, ncol(X_trial_temp)), 
                       Sigma = diag(sd_value_noise,ncol(X_trial_temp)))
    X_trial[[v]] <- X_trial_temp + X_noise
  }
  return(list("X_trial" = X_trial,
              "true_row_clustering" = true_row_clusterings,
              "true_col_clustering" = true_col_clusterings))
}


#### Test the above function:
input_n_views <- 2
input_row_dim <- vector("list", length = input_n_views)
input_row_dim[[1]] <- c(50,100,150,200)
input_row_dim[[2]] <- c(50,100,150,200)
#input_row_dim[[3]] <- c(100,100,100,100,100,100)
input_col_dim <- vector("list", length = input_n_views)
input_col_dim[[1]] <- c(500,50,100,250)
input_col_dim[[2]] <- c(200,50,120,220)
#input_col_dim[[3]] <- c(450,250,50,5,5,250)
data_block_diagonal <- generate_blockDiagonalMatrices_sameNclusters(n_views = input_n_views,
                                                                    row_dim = input_row_dim,
                                                                    col_dim = input_col_dim)
X_trial <- data_block_diagonal$X_simulated
true_row_clusterings <- data_block_diagonal$true_row_clustering
true_col_clusterings <- data_block_diagonal$true_col_clustering

## UNIFORM DISTRIBUTION
data_block_diagonal <- generate_blockDiagonalMatrices_sameNclusters_uniform(n_views = input_n_views,
                                                                    row_dim = input_row_dim,
                                                                    col_dim = input_col_dim)
X_trial <- data_block_diagonal$X_simulated
true_row_clusterings <- data_block_diagonal$true_row_clustering
true_col_clusterings <- data_block_diagonal$true_col_clustering
