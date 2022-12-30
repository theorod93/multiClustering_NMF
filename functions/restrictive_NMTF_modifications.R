#' Modifications on the rest_multi_NMTF_functions.R
#' 

# Check if a matrix is singular
f <- function(m) class(try(solve(m),silent=T))=="matrix"
# if f(x) = FALSE --> x is  a singular matrix

# Algorithm
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
#    Q_F <- single_alt_l1_normalisation(currentF[[v]])$Q
#    currentF[[v]] <- currentF[[v]] %*% solve(Q_F)
#    currentS[[v]] <- Q_F %*% currentS[[v]] 
    # Update G
    currentG[[v]] <- update_G(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS[[v]],
                              Ginput = currentG,
                              psi = psi, k = v)
    # Normalise F and S
#    Q_G <- single_alt_l1_normalisation(currentG[[v]])$Q
#    currentG[[v]] <- currentG[[v]] %*% solve(Q_G)
#    currentS[[v]] <- currentS[[v]] %*% t(Q_G)
    # Update S
    currentS[[v]] <- update_S(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS,
                              Ginput = currentG[[v]],
                              xi = xi, k = v)
    # Normalise F and S
    StGt <- t(currentS[[v]]%*% t(currentG[[v]]))
    # Replace NaN with 0 -> To be eligible for SVD
    StGt[which(is.na(StGt))] <- 0
    # Check if StGt, which is the input of the normalisation process, is singular
    if (f(StGt)){ # True means non-singular and we can proceed.
      Q <- single_alt_l1_normalisation(StGt)$Q
      currentF[[v]] <- currentF[[v]] %*% (Q)
      currentS[[v]] <- solve(Q) %*% currentS[[v]]  
    } else {
      svd_StGt <- svd(StGt)
      Q <- single_alt_l1_normalisation(svd_StGt$u)$Q
      currentF[[v]] <- currentF[[v]] %*% (Q)
      currentS[[v]] <- solve(Q) %*% currentS[[v]]  
    }
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  return(list("Foutput" = Foutput,"Soutput" = Soutput,"Goutput" = Goutput))
}


# Main function
restMultiNMTF_run_mod <- function(Xinput, Finput, Sinput, Ginput, phi, xi, psi, nIter){
  #' Run Restrictive-Multi NMTF, following the above algorithm
  n_v <- length(Xinput)
  # Normalise Xinput 
  for (v in 1:n_v){
    Xinput[[v]] <- single_alt_l1_normalisation(Xinput[[v]])$newMatrix
  }
#    for (v in 1:n_v){
#      Xinput[[v]] <- scale(Xinput[[v]])
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


# Initialisation function
initialisation_NMTF <- function(Xinput, K, L){
  #' Function to initialise F, S, G using the input matrix X
  #' Xinput : List with data-views as its elements (length = n_v)
  #' K : vector indicating the number of row-clusters in each data-view (length = n_v)
  #' L : vector indicating the number of column-clusters in each data-view (length = n_v)
  #' 
  #' 1. Apply k-means on columns -> G = C_c + 0.2, C_c:cluster centroids
  #' 2. Apply k-means on rows -> F = C_r + 0.2, C_r:cluster centroids
  #' 3. S = F^T X G
  #' ** Repeat the abovefor each view v 
  n_v <- length(Xinput)
  Finit <- vector("list", length = n_v)
  Sinit <- vector("list", length = n_v)
  Ginit <- vector("list", length = n_v)
  for (v in 1:n_v){
    kmeans_col <- kmeans(x = Xinput[[v]], centers = L[v])
    kmeans_row <- kmeans(x = t(Xinput[[v]]), centers = K[v])
    Finit[[v]] <- t(kmeans_row$centers) + 0.2
    Ginit[[v]] <- t(kmeans_col$centers) + 0.2
    Sinit[[v]] <- t(Finit[[v]]) %*% Xinput[[v]] %*% Ginit[[v]]
  }
  return(list("Finit" = Finit,
              "Sinit" = Sinit,
              "Ginit" = Ginit))
  #' Output: Finit, Sinit, Ginit -- All lists with data-views as its elements (length = n_v)
}


############## NEW UPDATES ###################
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
#  denominator_matrix <- currentF %*% Sinput %*% t(Ginput) %*% Ginput %*% t(Sinput)
  denominator_matrix <- currentF %*% t(currentF) %*% Xinput %*% Ginput %*% t(Sinput)
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
#  denominator_matrix <- currentG %*% t(Sinput) %*% t(Finput) %*% Finput %*% Sinput
  denominator_matrix <- currentG %*% t(currentG) %*% t(Xinput) %*% Finput %*% Sinput
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
