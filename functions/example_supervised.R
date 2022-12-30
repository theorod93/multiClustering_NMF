#' An example of running Multi-view restrictive Non-negative Matrix Tri-Factorisation
#' 
#' Use dataExample as input matrix
#' 

# Directory
#setwd("C:/Users/theod/OneDrive - Imperial College London/Documents/Imperial/PhD/Multi-Clustering/Restrictive Multi-NMTF/Data/nutriMouse/")

# Libraries
library(data.table)
library(MASS)
library(clue)
library(fossil)
library(mclust)
library(aricode)
library(whitening)

# Data
# dataExample: List 
# labelsExample

# Split train and test
dataExample_train <- vector("list")
dataExample_test <- vector("list")
sam <- sample(nrow(dataExample_data$genes), floor(0.7*nrow(dataExample_data$genes)))
dataExample_train <- dataExample[sam,]
dataExample_test <- dataExample[-sam,]
labelsExample_train <- vector("list")
labelsExample_train <- labelsExample[sam]
labelsExample_test <- vector("list")
labelsExample_test <- labelsExample[-sam]

#### RUN Restrictive Multi-NMTF ####
# True clustering

true_row_clusterings <- vector("list", length = length(dataExample_train))
true_col_clusterings <- vector("list", length = length(dataExample_train))
names(true_row_clusterings) <- names(dataExample_train)
for (i in 1:length(dataExample_train)){
  true_row_clusterings[[i]] <- labelsExample_train
  true_col_clusterings[[i]] <- rep(1, ncol(dataExample_train[[i]]))
}

# Initialization
Finit <- vector("list", length = length(dataExample_train))
Sinit <- vector("list", length = length(dataExample_train))
Ginit <- vector("list", length = length(dataExample_train))
for (i in 1:length(dataExample_train)){
  ss <- svd(dataExample_train[[i]])
  Finit[[i]] <- ss$u
  Sinit[[i]] <- diag(ss$d)
  Ginit[[i]] <- ss$v
}
KK <- rep(length(unique(labelsExample)),length(dataExample_train))
LL <- rep(2, length(dataExample_train))
for (i in 1:length(dataExample_train)){
  L <- LL[[i]]
  K <- KK[[i]]
  Sinit[[i]] <- mvrnorm(n = K, mu = rep(0,L), Sigma = diag(L))
}
# Select top K for row-clusters and L for column clusters

for (i in 1:length(dataExample_train)){
  L <- LL[[i]]
  K <- KK[[i]]
  Finit[[i]] <- abs(Finit[[i]][,1:K])
  Sinit[[i]] <- abs(Sinit[[i]][1:K,1:L])
  Ginit[[i]] <- abs(Ginit[[i]][,1:L])
}

# Tuning parameters
phi <- matrix(0, nrow = length(dataExample_train), ncol = length(dataExample_train))
xi <- matrix(0, nrow = length(dataExample_train), ncol = length(dataExample_train))
psi <- matrix(0, nrow = length(dataExample_train), ncol = length(dataExample_train))
nIter <- 1000

test_phi <- c(0, 0.01, 0.05,0.1, 0.2, 0.5, 0.65, 0.8, 1.0, 5.0, 10.0, 20.0, 50.0)
tuningAccuracy <- vector("list", length = length(test_phi))
tuningAdjRandValue <- vector("list", length = length(test_phi))
tuningNMIvalue <- vector("list", length = length(test_phi))
tuningError <- vector("list", length = length(test_phi))

# Assume n_v = 2 data-views
k <- 1
for (k in 1:length(test_phi)){
  phi[1,2] <- test_phi[k]
  dataExample_NMTF <- restMultiNMTF_run(Xinput = dataExample_train,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)
  eval_measures_dataExample_NMTF <- evaluate_simulation(X_nmtf = dataExample_NMTF,
                                                        true_row_clustering = true_row_clusterings,
                                                        true_col_clustering = true_col_clusterings)
  tuningAccuracy[[k]] <- eval_measures_dataExample_NMTF$accuracy
  tuningAdjRandValue[[k]] <- eval_measures_dataExample_NMTF$ARI
  tuningNMIvalue[[k]] <- eval_measures_dataExample_NMTF$NMI
  tuningError[[k]] <- dataExample_NMTF$Error
  if (k == 1){
    finalAccuracy_single <- tuningAccuracy[[k]]
    finalAdjRandValue_single <- tuningAdjRandValue[[k]]
    finalNMIvalue_single <- tuningNMIvalue[[k]]
    finalError_single <- tuningError[[k]]
  }
  if (test_phi[k] == 1){
    finalAccuracy_allOnes <- tuningAccuracy[[k]]
    finalAdjRandValue_allOnes <- tuningAdjRandValue[[k]]
    finalNMIvalue_allOnes <- tuningNMIvalue[[k]]
    finalError_allOnes <- tuningError[[k]]
  }
}

# Find optimal parameter:
meanACC <- lapply(tuningAccuracy, rowMeans)
chosen_metric <- 1
temp_mean <- meanACC[[1]][1]
for (i in 2:length(tuningAccuracy)){
  if (meanACC[[i]][1] > temp_mean){
    temp_mean <- meanACC[[i]][1]
    chosen_metric <- i
  }
}

phi[1,2] <- test_phi[chosen_metric]
dataExample_NMTF <- restMultiNMTF_run(Xinput = dataExample_train,
                                      Finput = Finit,
                                      Sinput = Sinit,
                                      Ginput = Ginit,
                                      phi = phi,
                                      xi = xi,
                                      psi = psi,
                                      nIter = nIter)
eval_measures_dataExample_NMTF <- evaluate_simulation(X_nmtf = dataExample_NMTF,
                                                      true_row_clustering = true_row_clusterings,
                                                      true_col_clustering = true_col_clusterings)

finalAccuracy_multiple <- eval_measures_dataExample_NMTF$accuracy
finalAdjRandValue_multiple <- eval_measures_dataExample_NMTF$ARI
finalNMIvalue_multiple <- eval_measures_dataExample_NMTF$NMI
finalError_multiple <- dataExample_NMTF$Error

## HSIC
hsic_values <- hsic_independence_test(dataExample_train)
phi <- hsic_values
phi[lower.tri(phi)] <- 0
diag(phi) <- 0
dataExample_NMTF <- restMultiNMTF_run(Xinput = dataExample_train,
                                      Finput = Finit,
                                      Sinput = Sinit,
                                      Ginput = Ginit,
                                      phi = phi,
                                      xi = xi,
                                      psi = psi,
                                      nIter = nIter)
eval_measures_dataExample_NMTF <- evaluate_simulation(X_nmtf = dataExample_NMTF,
                                                      true_row_clustering = true_row_clusterings,
                                                      true_col_clustering = true_col_clusterings)
finalAccuracy_hsic <- eval_measures_dataExample_NMTF$accuracy
finalAdjRandValue_hsic <- eval_measures_dataExample_NMTF$ARI
finalNMIvalue_hsic <- eval_measures_dataExample_NMTF$NMI
finalError_hsic <- dataExample_NMTF$Error

## Test dataset
Ftest <- dataExample_test %*% ginv(t(dataExample_NMTF$Goutput[[1]])) %*% ginv(dataExample_NMTF$Soutput[[1]])
F_clusters <- apply(Ftest, 1, which.max)


