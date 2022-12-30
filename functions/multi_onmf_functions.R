#' Create function to implement NMF
#' 
#' Following the results from
#' Ding et al (2006) Orthogonal Nonnegative Matrix Tri-factorizations for Clustering
#' 
#' and expanding them to multi-view data, V:num of views
#' Updates:
#' |-> S_v - forall V
#' |-> G_v - forall V
#' |-> F would be the same forall v
#' 

# Directory
setwd("C:/Users/theod/OneDrive - Imperial College London/Documents/Imperial/PhD/Multi-Clustering/Functions")

# Libraries
library(clue)
library(aricode)
library(fossil)
library(mclust)
library(data.table)
library(MASS)
library(stringr)

# Functions:

update_F_onmf_fTf <- function(Finput, Ginput, Xinput){
  # Updating F with the constraint F^TF=I
  # Input
  # F - row_clusters, nxk
  # Return F
  # Update given by:
  new_numer <- Xinput %*% Ginput
  new_denom <- (Finput %*% t(Finput)) %*% (Xinput %*% Ginput)
  # element-by-element
  Ftemp <- Finput
  for (i in nrow(Finput)){
    for (k in ncol(Finput)){
      Ftemp[i,k] <- Finput[i,k]*sqrt(new_numer[i,k]/new_denom[i,k])
    }
  }
  # Output
  return(Ftemp)
}
update_G_onmf_fTf <- function(Finput, Ginput, Xinput){
  # Updating G with the constraint F^TF=I
  # Input
  # G - column_cluster, pxk
  # Return G
  # Update given by:
  new_numer <- t(Xinput) %*% Finput
  new_denom <- Ginput %*% (t(Finput) %*% Finput)
  # element-by-element
  Gtemp <- Ginput
  for (j in nrow(Ginput)){
    for (k in ncol(Ginput)){
      Gtemp[j,k] <- Ginput[j,k]*(new_numer[j,k]/new_denom[j,k])
    }
  }
  # Output
  return(Gtemp)
}
update_F_onmf_gTg <- function(Finput, Ginput, Xinput){
  # Updating F with the constraint G^TG=I
  # Input
  # F - row_clusters, nxk
  # Return F
  # Update given by:
  new_numer <- Xinput %*% Ginput
  new_denom <- Finput %*% (t(Ginput) %*% Ginput)
  # element-by-element
  Ftemp <- Finput
  for (i in nrow(Finput)){
    for (k in ncol(Finput)){
      Ftemp[i,k] <- Finput[i,k]*(new_numer[i,k]/new_denom[i,k])
    }
  }
  # Output
  return(Ftemp)
}
update_G_onmf_gTg <- function(Finput, Ginput, Xinput){
  # Updating G with the constraint G^TG=I
  # Input
  # G - column_cluster, pxk
  # Return G
  # Update given by:
  new_numer <- t(Xinput) %*% Finput
  new_denom <- (Ginput %*% t(Ginput)) %*% (t(Xinput) %*% Finput)
  # element-by-element
  Gtemp <- Ginput
  for (j in nrow(Ginput)){
    for (k in ncol(Ginput)){
      Gtemp[j,k] <- Ginput[j,k]*sqrt(new_numer[j,k]/new_denom[j,k])
    }
  }
  # Output
  return(Gtemp)
}


onmf_single <- function(X, Finit, Ginit, num_clusters, nIter){
  # Input
  # X: Data matrix, nxp
  # Finit: Initial matrix of F, o.w. simulate standard multi-Normal
  # Ginit: Initial matrix of G, o.w. simulate standard multi-Normal
  # num_clusters: Number of clusters
  # nIIter: Max number of iterations
  
  if (is.null(Finit)){
    Finit <- mvrnorm(n = nrow(X), mu = rep(0,num_clusters), Sigma = diag(num_clusters))
  }
  if (is.null(Ginit)){
    Ginit <- mvrnorm(n = ncol(X), mu = rep(0,num_clusters), Sigma = diag(num_clusters))
  }
  
  # Return both F_orthogonal and G_orthohonal
  F_Forthogonal <- Finit
  G_Forthogonal <- Ginit
  F_Gorthogonal <- Finit
  G_Gorthogonal <- Ginit
  for (n in 1:nIter){
    # Update F:
    F_Forthogonal <- update_F_onmf_fTf(Finput = F_Forthogonal, 
                                       Ginput = G_Forthogonal, 
                                       Xinput = X)
    F_Gorthogonal <- update_F_onmf_gTg(Finput = F_Gorthogonal, 
                                       Ginput = G_Gorthogonal, 
                                       Xinput = X)
    # Update G:
    G_Forthogonal <- update_G_onmf_fTf(Finput = F_Forthogonal, 
                                       Ginput = G_Forthogonal, 
                                       Xinput = X)
    G_Gorthogonal <- update_G_onmf_gTg(Finput = F_Gorthogonal, 
                                       Ginput = G_Gorthogonal, 
                                       Xinput = X)
  }
  # Return list:
  # F-orthogonal
  output_Forthogonal <- vector("list")
  output_Forthogonal$Foutput <- F_Forthogonal
  output_Forthogonal$Goutput <- G_Forthogonal
  # G-orthogonal
  output_Gorthogonal <- vector("list")
  output_Gorthogonal$Foutput <- F_Gorthogonal
  output_Gorthogonal$Goutput <- G_Gorthogonal
  # Final Output

  return(list("X" = X,
              "F_orthogonal" = output_Forthogonal,
              "G_orthogonal" = output_Gorthogonal))
}

# Store
save.image("onmf_functions.RData")


# Test functions:
# 'simulations_multi_onmf.R'
