#--------------------------------------------------------------------------------
#' Construct a constrained covariance matrix from the adjacency matrix
#--------------------------------------------------------------------------------
#'
#' This function constructs a covariance function from the graph's adjacency
#' matrix. The covariance function may be squared exponential,
#' rational quadratic, Matern or the matrix exponential.
#' It includes a constraint, where a linear combination of the parameters
#' can be fixed.
#'
#'
#' @param adj.matrix The graph adjacency matrix
#' @param type The type of covariance function used. One of "sqexp", "ratquad", "matern" or "matrix". Note: only matern with nu = 5/2 is supported.
#' @param hyperparameters A vector containing the covariance function hyperparameters. 
#' For the squared exponential and matern, the vector should contain the variance 
#' and length scale, for the rational quadratic, the vector should contain the variance, l
#' ength scale and scaling parameters
#' @param linear.combination A matrix which defines the linear combination of 
#' the parameter vector lambda = (lambda_1, ..., lambda_N)^T. The linear 
#' combination is a vector of coefficients such that 
#' linear.combination %*% lambda = linear.constraint.
#' @param linear.constraint The value the linear constraint takes. 
#' Defaults to 0.
#' @return The mean vector and covariance matrix
#'
#' @seealso For more information about covariance functions see \url{https://www.cs.toronto.edu/~duvenaud/cookbook/} or \url{http://gaussianprocess.org/gpml/chapters/RW4.pdf}
#'
#' @examples
#' #Construct covariance matrix of Dar es Salaam, Tanzania, using network metric
#' data(dar.adj.matrix, package = "BSBT") #load dar es salaam adjacency matrix
#' k <- constrained_adjacency_covariance_function(dar.adj.matrix, type = "sqexp",
#'        hyperparameters = c(1, 1), rep(1, dim(dar.adj.matrix)[1]), 0)
#'        #Covariance registered by sum of objects is 0 using rational quadratic function

#--------------------------------------
# Additive decomposition
#--------------------------------------

# https://www.cs.toronto.edu/~duvenaud/cookbook/
# One nice thing about constructing an additive function is that you can decompose 
# your posterior over functions into additive parts. That is, if the kernel is a sum
# ksum(x,x′)=k1(x,x′)+k2(x,x′)+⋯+kD(x,x′)
# then your posterior can also be decomposed into a sum of Gaussian processes, each with mean and 
# covariance given by the corresponding part of the kernel:
# Mean(fd(x⋆))=kd(x⋆,X)Ksum(X,X)−1f(X)
# and variance
# Cov(fd(x⋆),fd(x⋆))=kd(x⋆,x⋆)−kd(x⋆,X)Ksum(X,X)−1kd(X,x⋆)


constrained_adjacency_covariance_function<- function (adj.matrix, type, hyperparameters, 
                                                      linear.combination,linear.constraint = 0) {
  set.seed(1234)
  if (type != "matrix") { 
    # directly get the graph here as adj_mtx
    # adj.matrix = adj_matrix
    # # if you directly feed the graph object as adj.mtx
    object.network <- igraph::graph_from_adjacency_matrix(adj.matrix)
    shortest.path.matrix <- igraph::distances(object.network, algorithm = 'dijkstra')
    #shortest.paths(object.network, algorithm = "dijkstra")
  }
  if (type == "sqexp" & length(hyperparameters) != 2) 
    stop("Insufficient hyperparameters. Squared Exponential requires 2 values.")
  if (type == "ratquad" & length(hyperparameters) != 3) 
    stop("Insufficient hyperparameters. Rational Quadratic requires 3 values.")
  if (type == "matern" & length(hyperparameters) != 2) 
    stop("Insufficient hyperparameters. Matern requires 2 values.")
  if (type == "matrix" & length(hyperparameters) != 2) 
    stop("Insufficient hyperparameters. Matrix exponential requires 2 values.")
  if (dim(adj.matrix)[1] != length(linear.combination)) 
    stop("Could not constrain distirbution. Linear constraint dimensions does not match number of objects.")
  
  #-------------------------
  # KERNELS
  #-------------------------
  if (type == "sqexp") { # Squared Exponential Kernel (SE)
    k <- hyperparameters[1]^2 * exp(-shortest.path.matrix^2/hyperparameters[2]^2)
  }
  else if (type == "ratquad") { # Rational Quadratic Kernel (RQ)
    k <- hyperparameters[1]^2 * (1 + shortest.path.matrix^2/
                                   (2 *(hyperparameters[2]^2) * hyperparameters[3]))^(-hyperparameters[3])
  }
  else if (type == "matern") { # this is matern 5/2 Kernel
    k <- (1 + sqrt(5)/hyperparameters[2] * shortest.path.matrix + 
            5/(3 * hyperparameters[2]^2) * shortest.path.matrix^2) * 
      exp(-sqrt(5)/hyperparameters[2] * shortest.path.matrix)
  }
  else if (type == "matrix") { # Matrix Exponential Kernel
    k <- expm::expm(hyperparameters[2] * adj.matrix)
    k <- hyperparameters[1]^2 * diag(diag(k)^-0.5) %*% k %*% 
      diag(diag(k)^-0.5)
  }
  else {
    stop("Could not construct covariance matrix. Unrecognised covariance type.")
  }
  if (length(linear.constraint) > 1) 
    stop("Currently only scalar constraints are supported")
  #---------------------------------------
  prior.mean <- rep(0, dim(k)[1])
  # Linear Constraint Handling enforcing a linear constraint on
  # the Gaussian process. Given the covariance matrix K, 
  # and we define a constrained covariance:
  constrained.k <- k - k %*% linear.combination %*% t(k %*% 
                                                        linear.combination) * as.numeric((1/(linear.combination %*% 
                                                                                               k %*% linear.combination)))
  constrained.mean <- prior.mean + as.numeric((linear.constraint - 
                                                 t(prior.mean) %*% linear.combination)/(t(linear.combination) %*% 
                                                                                          k %*% linear.combination)) * k %*% linear.combination
  # Eigen Decomposition and Inversion
  spectral.decomp <- eigen(constrained.k)
  # make sure the non negative eigen values are made 0
  spectral.decomp$values[spectral.decomp$values < 0] <- 0
  constrained.decomp <- spectral.decomp$vectors %*% diag(sqrt(spectral.decomp$values))
  constrained.inv <- MASS::ginv(constrained.k)
  colnames(constrained.decomp)<- colnames(constrained.inv) <- colnames(constrained.k)
  rownames(constrained.decomp)<- rownames(constrained.inv) <- rownames(constrained.k)
  return(list(mean = constrained.mean, covariance = constrained.k, 
              decomp = constrained.decomp, inv = constrained.inv))
}

