#This code appears to implement Stage 1 calculations for an EM algorithm to 
#estimate latent positions in a network model. It includes functions to convert between 
#vector and matrix representations of lower triangular matrices, as well as the main EM algorithm function itself.

func.mtv<-function(mat){
  lower=mat
  lower[upper.tri(lower,diag=T)]<-NA
  c.lower=c(lower)
  c.lower[!is.na(c.lower)]
}

vector_to_lower_tri <- function(vec, matrix_size = NULL, 
                                         fill_upper = 0, fill_diag = 0,
                                         symmetric = TRUE, validate = TRUE) {
  k <- length(vec)
  
  if (!is.null(matrix_size)) {
    # Use provided matrix size
    n <- matrix_size
    expected_k <- n * (n - 1) / 2
    
    if (k != expected_k && validate) {
      stop(paste("Vector length", k, "doesn't match expected", expected_k, 
                 "for", n, "x", n, "matrix"))
    }
    
    # Pad or truncate vector if needed
    if (k < expected_k) {
      vec <- c(vec, rep(0, expected_k - k))
      warning(paste("Vector padded with", expected_k - k, "zeros"))
    } else if (k > expected_k) {
      vec <- vec[1:expected_k]
      warning(paste("Vector truncated to", expected_k, "elements"))
    }
    
  } else {
    # Calculate matrix size from vector length
    n <- (1 + sqrt(1 + 8*k)) / 2
    
    if (abs(n - round(n)) > 1e-10 && validate) {
      valid_lengths <- sapply(2:20, function(x) x*(x-1)/2)
      stop(paste("Invalid vector length:", k, 
                 "\nValid lengths:", paste(valid_lengths, collapse = ", ")))
    }
    
    n <- round(n)
  }
  
  # Create matrix
  result_matrix <- matrix(fill_upper, n, n)
  diag(result_matrix) <- fill_diag
  
  # Fill lower triangle
  result_matrix[lower.tri(result_matrix)] <- vec
  
  # Make symmetric if requested
  if (symmetric) {
    result_matrix[upper.tri(result_matrix)] <- t(result_matrix)[upper.tri(result_matrix)]
  }
  
  return(result_matrix)
}

func.vtm<-function(vec){
  v = length(c(vec))
  n = (1+sqrt(8*v+1))/2
  lower=matrix(rep(0,n^2),n)
  lower[lower.tri(lower, diag=FALSE)]=vec
  lower+t(lower)
}


EM_toU <- function(G, d=3, u1=0.5, sig.a=2, sig.u=2, tol=10^(-10), maxiter=200){
  library(igraph);library(MASS);
  # G=adj_matrix[[i]]; d=3; u1=0.5; sig.a=2; sig.u=2; tol=10^(-10); maxiter=200
  G = c(G[lower.tri(G, diag=FALSE)]);v = length(G); p = (1+sqrt(8*v+1))/2
  sG = func.vtm(vec = G) 
  g = graph_from_adjacency_matrix(sG)
  D = distances(g, mode='all')
  #D = mean_distance(g)
  D[D==Inf] = max(D[D!=Inf])+mean_distance(g)
  #z0 = cmdscale(D,(d-1))
  z0 <- D[, 1:(d-1), drop = FALSE]
  # Normalize
  z0 <- scale(z0)
  
  u1=0.5
  U = cbind(u1,matrix(z0,ncol=(d-1)))
  # a = -log(1/mean(G)-1)
  a = 0
  step = 0; diff = 1
  while(step<=maxiter & diff>tol){
    # E step
    dotU = func.mtv(U%*%t(U))
    phi = a + dotU
    W = (exp(phi)-1)/(2*phi*(1+exp(phi)))
    mW = func.vtm(vec = W)
    
    # M step
    # max a
    old.a = a
    a = sum(G-W*dotU-0.5)/(1/sig.a + sum(W))
    
    old.U = U
    
    # max U[k,-1]
    for(k in 1:p){
      tmp1 = t(U[-k,-1])%*%diag(mW[-k,k])%*%U[-k,-1] + (1/sig.u)*diag(rep(1,(d-1)))
      tmp2 = matrix(sG[-k,k]-mW[-k,k]*(a+u1^2)-0.5,nrow=1)%*%U[-k,-1]
      U[k,-1] = tmp2%*%ginv(tmp1)
    }
    
    diff.U = max(abs(c(old.U-U)))
    diff.a = abs(old.a-a)
    diff = max(diff.U, diff.a)
    step = step + 1
    print(step)
  }
  list(U=U,a=a)
}

EM_toU_ver0 <- function(G, d=3, tol=10^(-10), maxiter=200){
  library(igraph)
  G = c(G);v = length(G); p = (1+sqrt(8*v+1))/2
  sG = func.vtm(G)
  g = graph_from_adjacency_matrix(sG)
  D = distances(g, mode='all')
  D[D==Inf] = max(D[D!=Inf])+1
  z0 = cmdscale(D,(d-1))
  
  u1=0.5
  U = cbind(u1,matrix(z0,ncol=(d-1)))
  a = -log(1/mean(G)-1)
  step = 0; diff = 1
  while(step<=maxiter & diff>tol){
    # E step
    dotU = func.mtv(U%*%t(U))
    phi = a + dotU
    W = (exp(phi)-1)/(2*phi*(1+exp(phi)))
    mW = func.vtm(W)
    
    # M step
    # max a
    old.a = a
    a = sum(G-W*dotU-0.5)/sum(W)
    
    old.U = U
    
    # max U[k,-1]
    for(k in 1:p){
      tmp1 = t(U[-k,-1])%*%diag(mW[-k,k])%*%U[-k,-1] + 0.5*diag(rep(1,(d-1)))
      tmp2 = matrix(sG[-k,k]-mW[-k,k]*(a+u1^2)-0.5,nrow=1)%*%U[-k,-1]
      U[k,-1] = tmp2%*%solve(tmp1)
    }
    
    diff.U = norm(old.U-U,"F")
    diff.a = abs(old.a-a)
    diff = max(diff.U, diff.a)
    step = step + 1
  }
  list(U=U,a=a)
}
