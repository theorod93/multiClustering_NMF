#' In this script we will write function with the updates of the restrictive Multi-NMTF
#' 
#' Functions of additional solutions will be presented for comparison
#' 
#' General problem: min{||X-FSG^T||^2}
#' |-> Need to find updates for F, S, G
#' 
#' X (NxP)
#' F (NxK)
#' S (KxL)
#' G (PxL)

# Directory
setwd("C:/Users/theod/OneDrive - Imperial College London/Documents/Imperial/PhD/Multi-Clustering/Restrictive Multi-NMTF")

# Libraries
library(MASS)

# General solution
update_F <- function(Xinput, Finput, Sinput, Ginput, phi, k){
  #' X: Input matrix
  #' F: row-clustering -- Entire list as input of length n_v
  #' S: connection between row-clustering and column-clustering
  #' G: column-clustering
  #' phi: weight on restrictions for F -> matrix of size (n_v x n_v)
  #' k: which view to update
  #' Output: An update for Finput[[k]] 
  
  # Find numerator
  currentF <- Finput[[k]]
  numerator_matrix <- Xinput %*% Ginput %*% t(Sinput)
  denominator_matrix <- currentF %*% Sinput %*% t(Ginput) %*% Ginput %*% t(Sinput)
  #denominator_matrix <- currentF %*% t(currentF) %*% Xinput %*% Ginput %*% t(Sinput)
  numerator <- matrix(0, nrow = nrow(currentF), ncol = ncol(currentF))
  denominator <- matrix(1, nrow = nrow(currentF), ncol = ncol(currentF))
  outputF <- currentF
  temp_den <- 0
  temp_num <- 0
  for (i in 1:nrow(currentF)){
    for (j in 1:ncol(currentF)){
      temp_1 <- 0
      temp_2 <- 0
      temp_3 <- 0
      for (u in 1:length(Finput)){
        if (phi[u,k] !=0){
          temp_1 <- temp_1 + phi[u,k]*Finput[[u]]
          temp_3 <- temp_3 + phi[u,k]
        }
        if (phi[k,u] !=0){
          temp_2 <- temp_2 + phi[k,u]*Finput[[u]]
          temp_3 <- temp_3 + phi[k,u]
        }
      }
      temp_den <- temp_3*Finput[[k]]
      temp_num <- temp_1 + temp_2
      if (is.na(sum(temp_num))){
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else if (sum(temp_num)== 0){
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else{
        numerator[i,j] <- numerator_matrix[i,j] + temp_num[i,j] 
        denominator[i,j] <- denominator_matrix[i,j] + temp_den[i,j]
      }
      outputF[i,j] <- currentF[i,j] * (numerator[i,j]/denominator[i,j])
    }
  }
  return(outputF)
}
update_S <- function(Xinput, Finput, Sinput, Ginput, xi, k){
  #' X: Input matrix
  #' F: row-clustering
  #' S: connection between row-clustering and column-clustering -- Entire list as input of length n_v
  #' G: column-clustering
  #' xi: weight on restrictions for S -> matrix of size (n_v x n_v)
  #' k: which view to update
  #' Output: An update for Finput[[k]] 

  # Find numerator
  currentS <- Sinput[[k]]
  numerator_matrix <- t(Finput) %*% Xinput %*% Ginput
  denominator_matrix <- t(Finput) %*% Finput %*% currentS %*% t(Ginput) %*% Ginput
  numerator <- matrix(0, nrow = nrow(currentS), ncol = ncol(currentS))
  denominator <- matrix(1, nrow = nrow(currentS), ncol = ncol(currentS))
  outputS <- currentS
  
  for (i in 1:nrow(currentS)){
    for (j in 1:ncol(currentS)){
      temp_1 <- 0
      temp_2 <- 0
      temp_3 <- 0
      for (u in 1:length(Sinput)){
        if (xi[u,k] !=0){
          temp_1 <- temp_1 + xi[u,k]*Sinput[[u]]
          temp_3 <- temp_3 + xi[u,k]
        }
        if (xi[k,u] !=0){
          temp_2 <- temp_2 + xi[k,u]*Sinput[[u]]
          temp_3 <- temp_3 + xi[k,u]
        }
      }
      temp_den <- temp_3*Sinput[[k]]
      temp_num <- temp_1 + temp_2
      if (is.na(sum(temp_num))){
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else if (sum(temp_num) == 0){
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else{
        numerator[i,j] <- numerator_matrix[i,j] + temp_num[i,j] 
        denominator[i,j] <- denominator_matrix[i,j] + temp_den[i,j]
      }
      outputS[i,j] <- currentS[i,j] * (numerator[i,j]/denominator[i,j])
    }
  }
  return(outputS)
}
update_G <- function(Xinput, Finput, Sinput, Ginput, psi, k){
  #' X: Input matrix
  #' F: row-clustering
  #' S: connection between row-clustering and column-clustering
  #' G: column-clustering -- Entire list as input of length n_v
  #' psi: weight on restrictions for G -> matrix of size (n_v x n_v)
  #' k: which view to update
  #' Output: An update for Finput[[k]] 
  
  # Find numerator
  currentG <- Ginput[[k]]
  numerator_matrix <- t(Xinput) %*% Finput %*% Sinput
  denominator_matrix <- currentG %*% t(Sinput) %*% t(Finput) %*% Finput %*% Sinput
  #denominator_matrix <- currentG %*% t(currentG) %*% t(Xinput) %*% Finput %*% Sinput
  numerator <- matrix(0, nrow = nrow(currentG), ncol = ncol(currentG))
  denominator <- matrix(1, nrow = nrow(currentG), ncol = ncol(currentG))
  outputG <- currentG
  temp_den <- 0
  temp_num <- 0
  for (i in 1:nrow(currentG)){
    for (j in 1:ncol(currentG)){
      temp_1 <- 0
      temp_2 <- 0
      temp_3 <- 0
      for (u in 1:length(Ginput)){
        if (psi[u,k] !=0){
          temp_1 <- temp_1 + psi[u,k]*Ginput[[u]]
          temp_3 <- temp_3 + psi[u,k]
        }
        if (psi[k,u] !=0){
          temp_2 <- temp_2 + psi[k,u]*Ginput[[u]]
          temp_3 <- temp_3 + psi[k,u]
        }
      }
      temp_den <- temp_3*Ginput[[k]]
      temp_num <- temp_1 + temp_2
      if (is.na(sum(temp_num))){
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else if (sum(temp_num) == 0){
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else{
        numerator[i,j] <- numerator_matrix[i,j] + temp_num[i,j] 
        denominator[i,j] <- denominator_matrix[i,j] + temp_den[i,j]
      }
      outputG[i,j] <- currentG[i,j] * (numerator[i,j]/denominator[i,j])
    }
  }
  return(outputG)
}

# Main function
restrictiveMultiNMTF_algo <- function(X, R, Finit, Sinit, Ginit, phi, xi, psi){
  #' X: Input matrix
  #' R: Matrix of size (3 x n_v) -> restrictions matrix
  #' phi: weight on restrictions for F -> matrix of size (n_v x n_v)
  #' psi: weight on restrictions for S -> matrix of size (n_v x n_v)
  #' xi: weight on restrictions for G -> matrix of size (n_v x n_v)
  #' Finit: Inital F matrix
  #' Sinit: Inital S matrix
  #' Ginit: Inital G matrix
  #' Output: Foutput, Soutput, Goutput
  
  # Update view-by-view
  n_v <- length(X)
  
#  # Nomralise data:
#  currentF <- vector("list", length = length(Finit))
#  currentS <- vector("list", length = length(Sinit))
#  currentG <- vector("list", length = length(Ginit))
#  for (i in 1:n_v){
#    X[[i]] <- scale(X[[i]])
#    currentF[[i]] <- scale(Finit[[i]])
#    currentS[[i]] <- scale(Sinit[[i]])
#    currentG[[i]] <- scale(Ginit[[i]])
#  }
  currentF <- Finit
  currentS <- Sinit
  currentG <- Ginit
  # Normalise:
#  normMatrices <- normalisation_l1(Finput = currentF,
#                                   Sinput = currentS, 
#                                   Ginput = currentG)
#  currentF <- normMatrices$Foutput
#  currentS <- normMatrices$Soutput
#  currentG <- normMatrices$Goutput
  for (v in 1:n_v){
    # Update F
    currentF[[v]] <- update_F(Xinput = X[[v]],
                              Finput = currentF,
                              Sinput = currentS[[v]],
                              Ginput = currentG[[v]],
                              phi = phi, k = v)
    # Normalise F
#    currentF[[v]] <- scale(currentF[[v]])
    # Update S
    currentS[[v]] <- update_S(Xinput = X[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS,
                              Ginput = currentG[[v]],
                              xi = xi, k = v)
    # Normalise S
#    currentS[[v]] <- scale(currentS[[v]])
    # Update G
    currentG[[v]] <- update_G(Xinput = X[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS[[v]],
                              Ginput = currentG,
                              psi = psi, k = v)
    # Normalise G
#    currentF[[v]] <- scale(currentF[[v]])
#    currentS[[v]] <- scale(currentS[[v]])
#    currentG[[v]] <- scale(currentG[[v]])
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  return(list("Foutput" = Foutput,"Soutput" = Soutput,"Goutput" = Goutput))
}

# Run the algorithm of Restrictive Multi-NMTF until convergence:
restrictiveMultiNMTF <- function(X, R, Finit, Sinit, Ginit, phi, xi, psi, nIter = NULL){
  #' X: Input matrix
  #' R: Matrix of size (3 x n_v) -> restrictions matrix
  #' phi: weight on restrictions for F -> matrix of size (n_v x n_v)
  #' psi: weight on restrictions for S -> matrix of size (n_v x n_v)
  #' xi: weight on restrictions for G -> matrix of size (n_v x n_v)
  #' Finit: Inital F matrix
  #' Sinit: Inital S matrix
  #' Ginit: Inital G matrix
  #' Output: Foutput, Soutput, Goutput
  #' 
  
  # Nomralise data:
#  normMatrices <- normalisation_l1(Finput = Finit,
#                                   Sinput = Sinit, 
#                                   Ginput = Ginit)
#  currentF <- normMatrices$Foutput
#  currentS <- normMatrices$Soutput
#  currentG <- normMatrices$Goutput
#  currentF <- vector("list", length = length(Finit))
#  currentS <- vector("list", length = length(Sinit))
#  currentG <- vector("list", length = length(Ginit))  
#  for (i in 1:n_v){
#    X[[i]] <- scale(X[[i]])
#    currentF[[i]] <- scale(Finit[[i]])
#    currentS[[i]] <- scale(Sinit[[i]])
#    currentG[[i]] <- scale(Ginit[[i]])
#  }
  # Initialisation
  currentF <- Finit
  currentS <- Sinit
  currentG <- Ginit
  Xhat <- vector("list", length = length(currentF))
  total_err <- c()
  if (is.null(nIter)){
    # Run while-loop until convergence
    mean_err <- 1
    while (mean_err > 0.01){
      err <- numeric(length = length(currentF))
      new_parameters <- restrictiveMultiNMTF_algo(X = X, 
                                R = R, 
                                Finit = currentF, 
                                Sinit = currentS,
                                Ginit = currentG,
                                phi = phi,
                                xi = xi, 
                                psi = psi)  
      currentF <- new_parameters$Foutput
      currentS <- new_parameters$Soutput
      currentG <- new_parameters$Goutput
      
      for (i in 1:length(currentF)){
        Xhat[[i]] <- currentF[[i]] %*% currentS[[i]] %*% t(currentG[[i]])
        err[i] <- sum(abs(X[[i]]-Xhat[[i]]))
      }
      mean_err <- mean(err)
      total_err <- c(total_err, mean_err)
    }
  } else {
    for (t in 1:nIter){
      err <- numeric(length = length(currentF))
      new_parameters <- restrictiveMultiNMTF_algo(X = X, 
                                                  R = R, 
                                                  Finit = currentF, 
                                                  Sinit = currentS,
                                                  Ginit = currentG,
                                                  phi = phi,
                                                  xi = xi, 
                                                  psi = psi)  
      currentF <- new_parameters$Foutput
      currentS <- new_parameters$Soutput
      currentG <- new_parameters$Goutput
      for (i in 1:length(currentF)){
        Xhat[[i]] <- currentF[[i]] %*% currentS[[i]] %*% t(currentG[[i]])
        err[i] <- sum(abs(X[[i]]-Xhat[[i]]))
      }
      mean_err <- mean(err)
      total_err <- c(total_err, mean_err)
    }
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  return(list("Foutput" = Foutput,"Soutput" = Soutput,
              "Goutput" = Goutput, "Error" = total_err))
}

# L1-Normalisation
normalisation_l1 <- function(Finput, Sinput, Ginput){
  #' Finput: List 
  #' Sinput: List 
  #' Ginput: List 
  #' 
  #' Output -> F,S,G :Lists
  
  # Compute normalization matrices:
  n_v <- length(Finput)
  Q_F <- vector("list", length = n_v)
  Q_G <- vector("list", length = n_v)
  for (v in 1:n_v){
    Q_F[[v]] <- diag(colSums(Finput[[v]]))
    Q_G[[v]] <- diag(colSums(Ginput[[v]]))
  }
  # Normalise:
  Foutput <- vector("list", length = n_v)
  Soutput <- vector("list", length = n_v)
  Goutput <- vector("list", length = n_v)
  for (v in 1:n_v){
    Foutput[[v]] <- Finput[[v]] %*% solve(Q_F[[v]])
    Soutput[[v]] <- Q_F[[v]] %*% Sinput[[v]] %*% Q_G[[v]]
    Goutput[[v]] <- Ginput[[v]] %*% t(solve(Q_G[[v]]))
  }
  return(list("Foutput" = Foutput, "Soutput" = Soutput, "Goutput" = Goutput))
}

## Additional functions
fixTablingClustering <- function(tableInput, max = TRUE){
  if (nrow(tableInput) != ncol(tableInput)){
    stop("PLease enter a square matrix")
  }
  lengthTable <- nrow(tableInput)
  orderTemp <- solve_LSAP(tableInput, maximum = max)
  tableTemp <- tableInput
  for (i in 1:lengthTable){
    tableTemp[i,] <- tableInput[which(orderTemp==i),]
  }
  return(tableTemp)
}
# Measuer Accuracy
accuracyTable <- function(tableInput){
  sumDiag <- sum(diag(tableInput))
  sumAll <- sum(tableInput)
  acc <- sumDiag/sumAll
  return(acc)
}


# Check if a matrix is singular
f <- function(m) class(try(solve(m),silent=T))=="matrix"
# if f(x) = FALSE --> x is  a singular matrix


#' New approach on the final function:
#' Follow this structure:
#' 1. Normalise X^v s.t ||X||_1 = 1
#' 2. Initialise F, S, G
#' 3. Repeat for each view u
#'    a. Fix ALL, update F^u
#'    b. Normalise F^u and S^u
#'    c. Fix ALL, update G^u
#'    c. Normalise G^u and S^u
#'    d. Fix ALL, update S^u
#' until convergence

single_l1_normalisation <- function(Xmatrix){
  #newMatrix <- matrix(0, nrow = nrow(Xmatrix), ncol = ncol(Xmatrix))
  newMatrix <- apply(Xmatrix,1, function(x) x/sum(x))
  return(newMatrix)
}

single_alt_l1_normalisation <- function(Xmatrix){
  #newMatrix <- matrix(0, nrow = nrow(Xmatrix), ncol = ncol(Xmatrix))
  # Introduce Q matrix
  Q <- diag(colSums(Xmatrix))
  newMatrix <- Xmatrix %*% solve(Q)
  return(list("Q" = Q, "newMatrix" = newMatrix))
}

restMultiNMTF_test <- function(Xinput, Finput, Sinput, Ginput, phi, xi, psi, nIter){
  #' Run Restrictive-Multi NMTF, following the above algorithm
  n_v <- length(Xinput)
  # Normalise Xinput 
  for (v in 1:n_v){
    Xinput[[v]] <- single_alt_l1_normalisation(Xinput[[v]])$newMatrix
  }
  
  # Take Finit, Sinit, Ginit as the initialised latent representations
  currentF <- Finput
  currentS <- Sinput
  currentG <- Ginput
  # Initialising additional parameters
  Xhat <- vector("list", length = n_v)
  total_err <- c()
  # Update until convergence, or for nIter times
  if (is.null(nIter)){
    # Run while-loop until convergence
    err_diff <- 1
    err_temp <- 0
    while (err_diff > 1.0e-4){
      err <- numeric(length = n_v)
      for (v in 1:n_v){
        # Update F
        currentF[[v]] <- update_F(Xinput = Xinput[[v]],
                                  Finput = currentF,
                                  Sinput = currentS[[v]],
                                  Ginput = currentG[[v]],
                                  phi = phi, k = v)
        # Normalise F and S
        Q_F <- single_alt_l1_normalisation(currentF[[v]])$Q
        currentF[[v]] <- currentF[[v]] %*% solve(Q_F)
        currentS[[v]] <- Q_F %*% currentS[[v]] 
        # Update G
        currentG[[v]] <- update_G(Xinput = Xinput[[v]],
                                  Finput = currentF[[v]],
                                  Sinput = currentS[[v]],
                                  Ginput = currentG,
                                  psi = psi, k = v)
        # Normalise F and S
        Q_G <- single_alt_l1_normalisation(currentG[[v]])$Q
        currentG[[v]] <- currentG[[v]] %*% solve(Q_G)
        currentS[[v]] <- currentS[[v]] %*% t(Q_G)
        # Update S
        currentS[[v]] <- update_S(Xinput = Xinput[[v]],
                                  Finput = currentF[[v]],
                                  Sinput = currentS,
                                  Ginput = currentG[[v]],
                                  xi = xi, k = v)
      }
      for (v in 1:n_v){
        Xhat[[v]] <- currentF[[v]] %*% currentS[[v]] %*% t(currentG[[v]])
        err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]]))
      }
      mean_err <- mean(err)
#      err_diff <- abs(mean_err - err_temp)
      total_err <- c(total_err, mean_err)
      err_diff <- abs(mean_err-err_temp)
      err_temp <- mean_err
    }
  } else {
    for (t in 1:nIter){
      err <- numeric(length = length(currentF))
      for (v in 1:n_v){
        # Update F
        currentF[[v]] <- update_F(Xinput = Xinput[[v]],
                                  Finput = currentF,
                                  Sinput = currentS[[v]],
                                  Ginput = currentG[[v]],
                                  phi = phi, k = v)
        # Normalise F and S
        Q_F <- single_alt_l1_normalisation(currentF[[v]])$Q
        currentF[[v]] <- currentF[[v]] %*% solve(Q_F)
        currentS[[v]] <- Q_F %*% currentS[[v]] 
        # Update G
        currentG[[v]] <- update_G(Xinput = Xinput[[v]],
                                  Finput = currentF[[v]],
                                  Sinput = currentS[[v]],
                                  Ginput = currentG,
                                  psi = psi, k = v)
        # Normalise F and S
        Q_G <- single_alt_l1_normalisation(currentG[[v]])$Q
        currentG[[v]] <- currentG[[v]] %*% solve(Q_G)
        currentS[[v]] <- currentS[[v]] %*% t(Q_G)
        # Update S
        currentS[[v]] <- update_S(Xinput = Xinput[[v]],
                                  Finput = currentF[[v]],
                                  Sinput = currentS,
                                  Ginput = currentG[[v]],
                                  xi = xi, k = v)
      }
      for (v in 1:n_v){
        Xhat[[v]] <- currentF[[v]] %*% currentS[[v]] %*% t(currentG[[v]])
        err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]]))
      }
      mean_err <- mean(err)
      total_err <- c(total_err, mean_err)
    }
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  return(list("Foutput" = Foutput, "Soutput" = Soutput,
              "Goutput" = Goutput, "Error" = total_err))
}

restMultiNMTF_algo <- function(Xinput, Finput, Sinput, Ginput, phi, xi, psi, nIter){
  #' Follow this structure:
  #' 1. Normalise X^v s.t ||X||_1 = 1
  #' 2. Initialise F, S, G
  #' 3. Repeat for each view u
  #'    a. Fix ALL, update F^u
  #'    b. Normalise F^u and S^u
  #'    c. Fix ALL, update G^u
  #'    c. Normalise G^u and S^u
  #'    d. Fix ALL, update S^u
  #' until convergence
  #' 
  #' X: Input matrix
  #' phi: weight on restrictions for F -> matrix of size (n_v x n_v)
  #' psi: weight on restrictions for S -> matrix of size (n_v x n_v)
  #' xi: weight on restrictions for G -> matrix of size (n_v x n_v)
  #' Finit: Inital F matrix
  #' Sinit: Inital S matrix
  #' Ginit: Inital G matrix
  #' Output: Foutput, Soutput, Goutput
  
  # Update view-by-view
  n_v <- length(Xinput)
  currentF <- Finput
  currentS <- Sinput
  currentG <- Ginput
  for (v in 1:n_v){
    # Update F
    currentF[[v]] <- update_F(Xinput = Xinput[[v]],
                              Finput = currentF,
                              Sinput = currentS[[v]],
                              Ginput = currentG[[v]],
                              phi = phi, k = v)
    # Normalise F and S
    # Check singular
    if (f(currentF[[v]])){
      Q_F <- single_alt_l1_normalisation(currentF[[v]])$Q
      currentF[[v]] <- currentF[[v]] %*% solve(Q_F)
      currentS[[v]] <- Q_F %*% currentS[[v]]   
    } else {
      svd_F <- svd(currentF[[v]])
      Q_F <- single_alt_l1_normalisation(svd_F$u)$Q
      currentF[[v]] <- currentF[[v]] %*% solve(Q_F)
      currentS[[v]] <- Q_F %*% currentS[[v]] 
    }
    # Update G
    currentG[[v]] <- update_G(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS[[v]],
                              Ginput = currentG,
                              psi = psi, k = v)
    # Normalise F and S
    if (f(currentG[[v]])){
      Q_G <- single_alt_l1_normalisation(currentG[[v]])$Q
      currentG[[v]] <- currentG[[v]] %*% solve(Q_G)
      currentS[[v]] <- currentS[[v]] %*% t(Q_G)
    } else {
      svd_G <- svd(currentG[[v]])
      Q_G <- single_alt_l1_normalisation(svd_G$u)$Q
      currentG[[v]] <- currentG[[v]] %*% solve(Q_G)
      currentS[[v]] <- currentS[[v]] %*% t(Q_G)
    }
    # Update S
    currentS[[v]] <- update_S(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS,
                              Ginput = currentG[[v]],
                              xi = xi, k = v)
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  return(list("Foutput" = Foutput,"Soutput" = Soutput,"Goutput" = Goutput))
}

# Main function
restMultiNMTF_run <- function(Xinput, Finput, Sinput, Ginput, phi, xi, psi, nIter){
  #' Run Restrictive-Multi NMTF, following the above algorithm
  n_v <- length(Xinput)
  # Normalise Xinput 
  for (v in 1:n_v){
    Xinput[[v]] <- single_alt_l1_normalisation(Xinput[[v]])$newMatrix
  }
  
  # Take Finit, Sinit, Ginit as the initialised latent representations
  currentF <- Finput
  currentS <- Sinput
  currentG <- Ginput
  # Initialising additional parameters
  Xhat <- vector("list", length = n_v)
  total_err <- c()
  # Update until convergence, or for nIter times
  if (is.null(nIter)){
    # Run while-loop until convergence
    err_diff <- 1
    err_temp <- 0
    while (err_diff > 1.0e-4){
      err <- numeric(length = n_v)
      new_parameters <- restMultiNMTF_algo(X = Xinput, 
                                           Finput = currentF, 
                                           Sinput = currentS,
                                           Ginput = currentG,
                                           phi = phi,
                                           xi = xi, 
                                           psi = psi)  
      currentF <- new_parameters$Foutput
      currentS <- new_parameters$Soutput
      currentG <- new_parameters$Goutput
      for (v in 1:n_v){
        Xhat[[v]] <- currentF[[v]] %*% currentS[[v]] %*% t(currentG[[v]])
        err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]]))
      }
      mean_err <- mean(err)
      #      err_diff <- abs(mean_err - err_temp)
      total_err <- c(total_err, mean_err)
      err_diff <- abs(mean_err-err_temp)
      err_temp <- mean_err
    }
  } else {
    for (t in 1:nIter){
      err <- numeric(length = length(currentF))
      new_parameters <- restMultiNMTF_algo(X = Xinput, 
                                           Finput = currentF, 
                                           Sinput = currentS,
                                           Ginput = currentG,
                                           phi = phi,
                                           xi = xi, 
                                           psi = psi)  
      currentF <- new_parameters$Foutput
      currentS <- new_parameters$Soutput
      currentG <- new_parameters$Goutput
      for (v in 1:n_v){
        Xhat[[v]] <- currentF[[v]] %*% currentS[[v]] %*% t(currentG[[v]])
        err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]]))
      }
      mean_err <- mean(err)
      total_err <- c(total_err, mean_err)
    }
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  return(list("Foutput" = Foutput, "Soutput" = Soutput,
              "Goutput" = Goutput, "Error" = total_err))
}


# Cross-validation Evalutation metric  
cv_eval_metric <- function(accList, errList, testValues, weights){
  #' This function is used to determine the optimal  tuning parameters
  #' Based on Accuracy and Error
  #' 
  #' accList: A list of with measurements on accuracy with each tuning parameter
  #' errList: A list of with measurements on error with each tuning parameter
  #' testValues: A vector with the tuning parameter values that were tested
  #' Length of lists = # of testing tuning parameters
  #' Length of testValues = Length of lists
  #' weights: vector of length 2, first element is the weight to error, and second to accuracy
  #' 
  newAccList <- lapply(accList, function(x) (1-x))
  metricList <- vector("list", length = length(accList))
  for (i in 1:length(metricList)){
    metricList[[i]] <- weights[1]*min(errList[[i]]) + weights[1]*mean(newAccList[[i]])
  }
  minCase <- which.min(metricList)
  minValue <- testValues[minCase]
  return(list("metricList" = metricList, "minCase" = minCase, "minValue" = minValue))
}

evaluate_simulation <- function(X_nmtf, true_row_clustering, true_col_clustering){
  #' X_nmtf: Output of restMultiNMTF_run
  #' true_row/col_clustering: What the name states
  #' 
  row_clustering <- vector("list", length = length(X_nmtf$Foutput))
  column_clustering <- vector("list", length = length(X_nmtf$Goutput))
  for (i in 1:length(X_nmtf$Foutput)){
    row_clustering[[i]] <- apply(X_nmtf$Foutput[[i]], 1, which.max)
    column_clustering[[i]] <- apply(X_nmtf$Goutput[[i]], 1, which.max)
  }
  column_table <- vector("list", length = length(X_nmtf$Foutput))
  row_table <- vector("list", length = length(X_nmtf$Foutput))
  accuracy <- matrix(0, nrow = 2, ncol = length(X_nmtf$Foutput))
  rownames(accuracy) <- c("Row-clustering", "Column-clustering")
  adjRandValue <- matrix(0, nrow = 2, ncol = length(X_nmtf$Foutput))
  rownames(adjRandValue) <- c("Row-clustering", "Column-clustering")
  nmiValue <- matrix(0, nrow = 2, ncol = length(X_nmtf$Foutput))
  rownames(nmiValue) <- c("Row-clustering", "Column-clustering")
  for (i in 1:length(X_nmtf$Foutput)){
    row_table[[i]] <- table(row_clustering[[i]], true_row_clustering[[i]])
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
    column_table[[i]] <- table(column_clustering[[i]], true_col_clustering[[i]])
    if (nrow(column_table[[i]]) < ncol(column_table[[i]])){
      for (g in 1:(ncol(column_table[[i]]) - nrow(column_table[[i]]))){
        column_table[[i]] <- rbind(column_table[[i]], rep(0,ncol(column_table[[i]])))
      }
    }else if (nrow(column_table[[i]]) > ncol(column_table[[i]])){
      for (g in 1:(nrow(column_table[[i]]) - ncol(column_table[[i]]))){
        column_table[[i]] <- cbind(column_table[[i]], rep(0,nrow(column_table[[i]])))
      }
    }
    column_table[[i]] <- fixTablingClustering(column_table[[i]])
    accuracy[1,i] <- accuracyTable(row_table[[i]])
    accuracy[2,i] <- accuracyTable(column_table[[i]])
    adjRandValue[1,i] <- adjustedRandIndex(row_clustering[[i]], true_row_clustering[[i]])
    adjRandValue[2,i] <- adjustedRandIndex(column_clustering[[i]], true_col_clustering[[i]])
    nmiValue[1,i] <- NMI(row_clustering[[i]], true_row_clustering[[i]])
    nmiValue[2,i] <- NMI(column_clustering[[i]], true_col_clustering[[i]])
  }
  return(list("accuracy" = accuracy, "ARI" = adjRandValue, "NMI" = nmiValue))
}

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
  
  # Return
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
  mean_matrix <- diag(rep(3,K))
  for (i in 1:K){
    row_start <- row_end + 1
    row_end <- row_end + nn[i]
    col_end <- 0
    for (j in 1:K){
#      col_start <- col_end + 1
#      col_end  <- col_end + pp[j]
      F_temp <- rnorm(n = nn[i], mean = mean_matrix[i,j], sd = 1)
      F_simulated[row_start:row_end, j] <- F_temp
    }
  }
  #F_simulated <- abs(F_simulated)
  #F_noise <- mvrnorm(n = N, mu = rep(0, K), Sigma = diag(K))
  #F_original <- abs(F_noise + F_simulated)
  F_original <- abs(F_simulated)
  row_clusters_original <- apply(F_original,1,which.max)
  
  # Simulate G
  G_simulated <- matrix(0, nrow = P, ncol = L)
  row_end <- 0
  mean_matrix <- diag(rep(3,L))
  for (i in 1:L){
    row_start <- row_end + 1
    row_end <- row_end + pp[i]
    col_end <- 0
    for (j in 1:L){
#      col_start <- col_end + 1
#      col_end  <- col_end + pp[j]
      G_temp <- rnorm(n = pp[i], mean = mean_matrix[i,j], sd = 1)
      G_simulated[row_start:row_end, j] <- G_temp
    }
  }
  #G_simulated <- abs(G_simulated)
  #G_noise <- mvrnorm(n = P, mu = rep(0, L), Sigma = diag(L))
  #G_original <- abs(G_noise + G_simulated)
  G_original <- abs(G_simulated)
  col_clusters_original <- apply(G_original,1,which.max)
  
  # Simulate S
  S_original <- abs(mvrnorm(n = K, mu = rep(0,L), Sigma = diag(L)))
  
  # Compute X
  X_simulated <- F_original %*% S_original %*% t(G_original)
  X_noise <- mvrnorm(n = nrow(X_simulated), mu = rep(1,ncol(X_simulated)), Sigma = diag(ncol(X_simulated)))
  X_original <- abs(X_noise + X_simulated)
  
  return(list("X" = X_original,
              "F" = F_original,
              "S" = S_original,
              "G" = G_original,
              "row_clusters" = row_clusters,
              "row_clusters_original" = row_clusters_original,
              "col_clusters" = col_clusters,
              "col_clusters_original" = col_clusters_original))
}

#########################################################################################################
########## HSIC #########################################################################################
#########################################################################################################

#' Measure HSIC:
#' |-> A promising dependence test
#' Potentially to be used as an indication to the tuning parameters in Restrictive Multi-NMTF

# Existing libraries:
library(kernelPSI)
library(kpcalg)

#' Hilbert-Schmidt Independence Criterion (HSIC)

hsic_measure <- function(K,L){
  #' Input:
  #' Two matrices, K and L, for which HSIC is to be computed
  #' NOTE: K and L are usually kernel functions
  #' Output:
  #' HSIC measure
  #' 
  #' NOTES: K and L must share the same dimensions and so does H
  N <- nrow(K)
  if (nrow(K) != nrow(L)){
    stop("Dimensions of the two input matrices must be the same")
  }
  if (ncol(K) != ncol(L)){
    stop("Dimensions of the two input matrices must be the same")
  }
  # Compute H matrix
  H <- matrix(0, nrow = nrow(K), ncol = ncol(K))
  for (i in 1:nrow(H)){
    for (j in 1:ncol(H)){
      if (i==j){
        H[i,j] <- 1 - (1/N^2)
      } else {
        H[i,j] <- 0 - (1/N^2)
        
      }
    }
  }
  
  hsic_denominator <- (N-1)^2
  hsic_numerator <- sum(diag(K%*%H%*%L%*%H))
  hsic_measure <- hsic_numerator / hsic_denominator
  return(hsic_measure)
}

gaussian_kernel <- function(X){
  sigma <- sd(X)
  
}

hsic_independence_test <- function(Xinput, kernel_choice = "euclidean"){
  #' Input:
  #' Xinput: A list of input matrices. HSIC of all possible combinations are to be computed
  #' kernel_choice: Compte the kernel function under the given function
  #' |-> Options: "gaussian", "euclidean"
  #' Analysis:
  #' K: Kernel matrix of first matrix
  #' L: Kernel matrix of second matrix
  #' HSIC: Use existing OR developed function
  #' Output: HSIC of (K,L)
  
  n_v <- length(Xinput)
  # Compute Kernels
  K <- vector("list", length = n_v)
  if (kernel_choice == "euclidean"){
    for (v in 1:n_v){
      K[[v]] <- scale(Xinput[[v]]) %*% t(scale(Xinput[[v]])) / ncol(Xinput[[v]])
    }
  } else if(kernel_choice == "gaussian"){
    for (v in 1:n_v){
      K[[v]] <- rbfkernel(scale(Xinput[[v]]))
    }
  }
  
  # Compute HSIC
  hsic_output <- matrix(0, nrow = n_v, ncol = n_v)
  for (i in 1:n_v){
    for (j in 1:n_v){
      hsic_output[i,j] <- hsic_measure(K[[i]], K[[j]])
    }
  }
  return(hsic_output)
}




#########################################################################################################
########## NO Normalisation #############################################################################
#########################################################################################################

# Main function
restMultiNMTF_run_noScale <- function(Xinput, Finput, Sinput, Ginput, phi, xi, psi, nIter){
  #' Run Restrictive-Multi NMTF, following the above algorithm
  n_v <- length(Xinput)
  # Normalise Xinput 
  #  for (v in 1:n_v){
  #    Xinput[[v]] <- single_alt_l1_normalisation(Xinput[[v]])$newMatrix
  #  }
  
  # Take Finit, Sinit, Ginit as the initialised latent representations
  currentF <- Finput
  currentS <- Sinput
  currentG <- Ginput
  # Initialising additional parameters
  Xhat <- vector("list", length = n_v)
  total_err <- c()
  # Update until convergence, or for nIter times
  if (is.null(nIter)){
    # Run while-loop until convergence
    err_diff <- 1
    err_temp <- 0
    while (err_diff > 1.0e-4){
      err <- numeric(length = n_v)
      new_parameters <- restMultiNMTF_algo(X = Xinput, 
                                           Finput = currentF, 
                                           Sinput = currentS,
                                           Ginput = currentG,
                                           phi = phi,
                                           xi = xi, 
                                           psi = psi)  
      currentF <- new_parameters$Foutput
      currentS <- new_parameters$Soutput
      currentG <- new_parameters$Goutput
      for (v in 1:n_v){
        Xhat[[v]] <- currentF[[v]] %*% currentS[[v]] %*% t(currentG[[v]])
        err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]]))
      }
      mean_err <- mean(err)
      #      err_diff <- abs(mean_err - err_temp)
      total_err <- c(total_err, mean_err)
      err_diff <- abs(mean_err-err_temp)
      err_temp <- mean_err
    }
  } else {
    for (t in 1:nIter){
      err <- numeric(length = length(currentF))
      new_parameters <- restMultiNMTF_algo(X = Xinput, 
                                           Finput = currentF, 
                                           Sinput = currentS,
                                           Ginput = currentG,
                                           phi = phi,
                                           xi = xi, 
                                           psi = psi)  
      currentF <- new_parameters$Foutput
      currentS <- new_parameters$Soutput
      currentG <- new_parameters$Goutput
      for (v in 1:n_v){
        Xhat[[v]] <- currentF[[v]] %*% currentS[[v]] %*% t(currentG[[v]])
        err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]]))
      }
      mean_err <- mean(err)
      total_err <- c(total_err, mean_err)
    }
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  return(list("Foutput" = Foutput, "Soutput" = Soutput,
              "Goutput" = Goutput, "Error" = total_err))
}
