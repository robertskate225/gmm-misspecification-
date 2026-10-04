#libraries
library(MASS)
#this library will be used to generate the skew-normal data
library(sn)
library(mvtnorm)
#this library is used to fit the skew normal mixture model
library(mixsmsn)



set.seed(2026)

#generating gmm data
gmm_generate <- function(n, pi, mu, sigma){
  
  #finding parameters
  k <- as.numeric(length(pi))
  p <- ncol(mu)
  
  #generate latent components
  #this will also allow us to see how well the model actually does asign the points
  Z_true <- sample(1:k, size = n, replace = TRUE, prob = pi)
  
  #initialise X matrix
  #will be using for loop to generate data
  X <- matrix(0, nrow = n, ncol = p)
  
  for(i in 1:k) {
    label_x <- which(Z_true == i)
    n_k <- length(label_x)
    
    if(n_k > 0){
      X[label_x, ] <- MASS::mvrnorm(n = n_k, mu = mu[i, ], Sigma = sigma[[i]])
    }
  }
  
  return(list(X=X, Z_true = Z_true))
}


#generate skew normal
#NOTE - the labels will still be included here (to be used for the label switching)
skew_norm_generate <- function(n, pi, mu, sigma, alpha){
  
  k <- as.numeric(length(pi))
  p <- ncol(mu)
  
  Z_true <- sample(1:k, size = n, replace = TRUE, prob = pi)
  
  X <- matrix(0, ncol = p, nrow = n)
  for(i in 1:k){
    label_x <- which(Z_true == i)
    n_k <- length(label_x)
    if(n_k > 0){
      X[label_x, ] <- sn::rmsn(n = n_k, xi = mu[i,], Omega = sigma[[i]], alpha = alpha[[i]])
    }
  }
  
  #store the true component labels in the same order
  return(list(X = X, Z_true = Z_true))
}



#em algorithm
#done for the estimation of the gmm parameters
em_gmm <- function(X, k, tol, max_iter = 1000){
  
  n <- nrow(X)
  p <- ncol(X)
  
  #1 - initial values
  #STEP 2 - K-means (initialisation for EM)
  kmean <- kmeans(X, centers = k, nstart = 25)
  
  #parameters - do not want them separately (want them "bind")
  mu0 <- kmean$centers
  #see which obs is assigned to which cluster (like an estimed Z)
  cluster <- kmean$cluster
  pi0 <- as.numeric(table(factor(cluster, levels = 1:k))) / n
  #initialise vector for sigma
  sigma0 <- vector("list", k)
  
  for(i in 1:k){
    X_k <- X[cluster == i, , drop = FALSE]
    sigma0[[i]] <- if (nrow(X_k) > p) {cov(X_k)} else {diag(p)}
  }
  
  
  #iteration counter
  iter <- 0
  
  #convergence indicator
  converged <- FALSE
  
  loglik <- -Inf
  
  while(!converged && iter < max_iter){
    
    #increased iteration counter
    iter <- iter+1
    
    
    #2 - E step
    
    #initialise matrix
    tau <- matrix(0, nrow = n, ncol = k)
    
    for(i in 1:k){
      tau[,i] <- pi0[i]*mvtnorm::dmvnorm(X, mean = mu0[i,], sigma = sigma0[[i]])
    }
    
    #normalise
    tau <- tau/rowSums(tau)
    
    
    #3 - M step
    #effective number of obs assigned to each component
    nk <- colSums(tau)
    
    #updated mixing proportions
    pi_new <- nk/n
    
    #updated means
    mu_new <- matrix(0, nrow = k, ncol = p)
    for(i in 1:k){
      mu_new[i,] <- colSums(tau[,i]*X)/nk[i]
    }
    
    #updated cov matrices
    sigma_new <- vector("list", k)
    for(i in 1:k){
      #difference obs and mean
      diff <- sweep(X, 2, mu_new[i,], "-")
      #weighted cov matrix
      sigma_new[[i]] <- t(diff) %*% (diff*tau[,i])/nk[i]
    }
    
    #log-likelihood
    density_sum <- matrix(0, nrow = n, ncol = k)
    
    for(i in 1:k){
      density_sum[,i] <- pi_new[i]*mvtnorm::dmvnorm(X, mean=mu_new[i,], sigma=sigma_new[[i]])
    }
    
    loglik_new <- sum(log(rowSums(density_sum)))
    
    
    #check convergence 
    loglik_diff <- abs(loglik_new-loglik)
    
    if(loglik_diff < tol){
      converged <- TRUE
    }
    
    #update parameters
    pi0 <- pi_new
    mu0 <- mu_new
    sigma0 <- sigma_new
    loglik <- loglik_new
  }
  
  return(list(pi = pi_new,mu =  mu_new, sigma = sigma_new, loglik = loglik_new))
  
}



#mapping funtion
#NB - creating a function to assist us with the label switching (using one of Yao's methods)
#this is now the generalized method (taking any p variables into account)
mapping <- function(k){
  if(k==1){
    return(matrix(1, 1, 1))
  }
  
  vec <- 1:k
  map_matrix <- do.call(rbind, lapply(vec, function(i){
    sub_map <- mapping(k-1)
    sub_map[sub_map >= i] <- sub_map[sub_map >= i] + 1
    cbind(i, sub_map)
  }))
  return(map_matrix)
}



#relabel yao (fix for label switching)
#general 
#NOTE - if k=1, this is not necessary!!!
relabel_yao <- function(X, Z_true, par){
  
  k <- length(par$pi)
  k_true <- length(unique(Z_true))
  
  #Constraint condition (do not to relabel when there is only 1 cluster)
  if(k == 1){
    return(par)
  }
  
  if(k == k_true){
  n <- nrow(X)
  
  #evaluate the densities under the current estimates
  dens <- matrix(0, nrow = n, ncol = k)
  for (j in 1:k){
    dens[, j] <- par$pi[j] * mvtnorm::dmvnorm(X, par$mu[j, ], par$sigma[[j]])
  }
  
  #using the general mapping function to take all possible mappings into account (permutations)
  map <- mapping(k)
  best_loglik <- -Inf 
  best_map <- 1:k
  
  for(j in 1:nrow(map)){
    current_map <- map[j, ]
    mapped_assign <- current_map[Z_true]
    current_loglik <- sum(log(dens[cbind(1:n, mapped_assign)]))
    
    if(current_loglik > best_loglik){
      best_loglik <- current_loglik
      best_map <- current_map
    }
  }
  
  #reorder the estimated parameters
  par$pi <- par$pi[best_map]
  par$mu <- par$mu[best_map, , drop = FALSE]
  par$sigma <- par$sigma[best_map]
  }
  #this will be done for the cluster misspecification
  else{
    #simply order
    order <- order(par$mu[,1])
    par$pi <- par$pi[order]
    par$mu <- par$mu[order, , drop = FALSE]
    par$sigma <- par$sigma[order]
  }
  return(par)
}

#calculating true density 
#this will be necessary for the MSE calculation
calc_true_dens <- function(X, pi, mu, sigma, alpha, type){
  
  n <- nrow(X)
  k_true <- as.numeric(length(pi))
  dens_matrix <- matrix(0, nrow = n, ncol = k_true)
  
  for(j in 1:k_true){
    if(type == "skew_normal"){
      dens_matrix[, j] <- pi[j]*sn::dmsn(X, xi = mu[j,], Omega = sigma[[j]], alpha = alpha[[j]])
    } else {
      dens_matrix[, j] <- pi[j]*mvtnorm::dmvnorm(X, mean = mu[j,], sigma = sigma[[j]])
    }
  }
  return(rowSums(dens_matrix))
}

#will also be computing the density MSE 
#this is another model evaluation
#NOTE - the fitted values will be input into this function
dens_mse <- function(X, true_dens_vec, pi, mu, sigma){
  
  n <- nrow(X)
  k_fit <- length(pi)
  
  #fitted density
  #this will be for each iteration in the MC sim
  fit_dens <- matrix(0, nrow = n, ncol = k_fit)
  for(j in 1:k_fit){
    sigma_j <- if(is.list(sigma)) {sigma[[j]]} else {sigma[j, , ]}
    fit_dens[, j] <- pi[j]*mvtnorm::dmvnorm(X, mean = mu[j, ], sigma = sigma_j)
  }
  fit_dens_vec <- rowSums(fit_dens)
  
  #MSE (Mean Squared Error)
  mse <- mean((fit_dens_vec - true_dens_vec)^2)
  return(mse)
}



#monte carlo simulation (as function)
#NOTE - creating a function for this as well - then we can simply use it for the different cases
#shortening model
mc <- function(sim, n, k, pi, mu, sigma, alpha, type, tol = 1e-6, max_iter = 1000){
  k_true <- length(pi)
  p <- ncol(mu)
  
  #initialise matrices
  pi_res <- matrix(0, nrow = sim, ncol = k)
  mu_res <- array(0, dim = c(sim, k, p))
  sigma_res <- array(0, dim = c(sim, k, p, p))
  loglik_res <- numeric(sim)
  AIC_res <- numeric(sim)
  BIC_res <- numeric(sim)
  dens_MSE_res <- numeric(sim)
  
  #number of parameters (used to calculate the AIC and BIC)
  num_par <- (k-1) + k*p + k*p*(p+1)/2
  
  #for loop
  for(i in 1:sim){
    if(type == "gaussian"){
      data <- gmm_generate(n, pi, mu, sigma)
    }
    else if (type == "skew_normal"){
      data <- skew_norm_generate(n, pi, mu, sigma, alpha)
    }
    
    #find parameters using the em algorithm
    par <- em_gmm(data$X, k = k, tol = tol, max_iter = max_iter)
    
    #relabel data (for the label switching)
    par <- relabel_yao(data$X, data$Z_true, par)
    
    
    #save the estimates
    pi_res[i, ] <- par$pi
    mu_res[i, , ] <- par$mu
    for(j in 1:k){
      sigma_res[i, j, , ] <- par$sigma[[j]]
    }
    
    #save the model evaluations
    loglik_res[i] <- par$loglik
    AIC_res[i] <- -2*par$loglik + 2*num_par
    BIC_res[i] <- -2*par$loglik + num_par*log(n)
    
    #calculate true density (and then Density MSE)
    true_dens <- calc_true_dens(data$X, pi, mu, sigma, alpha, type)
    dens_MSE_res[i] <- dens_mse(data$X, true_dens, par$pi, par$mu, par$sigma)
  }
  
  #need to return the estimates
  #want the mean, and variance of each estimate
  return(list( data = data,
    #first the raw estimates
    raw_res = list(pi = pi_res, mu = mu_res, sigma = sigma_res, loglik = loglik_res, AIC = AIC_res, BIC = BIC_res, Density_MSE = dens_MSE_res),
    Mean = list(pi_hat = colMeans(pi_res), mu_hat = apply(mu_res, c(2, 3), mean), sigma_hat = apply(sigma_res, c(2, 3, 4), mean), 
                Evaluation = data.frame(loglik = mean(loglik_res), AIC = mean(AIC_res), BIC = mean(BIC_res), Density_MSE = mean(dens_MSE_res))),
    Variance = list(pi_hat_var = apply(pi_res, 2, var), mu_hat_var = apply(mu_res, c(2, 3), var), sigma_hat_var = apply(sigma_res, c(2, 3, 4), var),
                    Evaluation = data.frame(loglik_var = var(loglik_res), AIC_var = var(AIC_res), BIC_var = var(BIC_res), Density_MSE_var = var(dens_MSE_res)))
    
  ))
  
}





#implementing everything:
#NOTE - this is simply an example (want to ensure that the code is working correctly)

############################################################################################
#  CORRECT MODEL
############################################################################################

#parameters
k <- 3
p <- 5
n <- 500
sim <- 500

pi <- c(0.5, 0.3, 0.2)
mu <- matrix(c(0, 0, 0, 0, 0,
               4, 4, 4, 4, 4,
               8, 8, 8, 8, 8), nrow = k, ncol = p, byrow = TRUE)
sigma <- list(diag(p), diag(p), diag(p))

alpha <- list(c(3, 3, 3, 3, 3), c(-3, -3, -3, -3, -3), c(3, 3, 3, 3, 3))

#run monte carlo
mc_gmm <- mc(sim = sim, n = n, k = k, pi = pi, mu = mu, sigma = sigma, type = "gaussian", alpha = alpha)



############################################################################################
#  NON-GAUSSIAN DATA
############################################################################################

mc_skew <- mc(sim = sim, n = n, k = k, pi = pi, mu = mu, sigma = sigma, alpha = alpha, type = "skew_normal")





############################################################################################
#  CLUSTER MISSPECIFICATION
############################################################################################

k_new <- 4

mc_cluster <- mc(sim = sim, n = n, k = k_new, pi = pi, mu = mu, sigma = sigma, alpha = alpha, type = "gaussian")



############################################################################################
#  PLOTS - ONLY FOR THE BIVARIATE CASE
############################################################################################

plot_gmm_2d <- function(X, k, k_new, pi, mu, sigma, alpha, pi_est, mu_est, sigma_est, type){
  
  p <- ncol(X)
  
  if(p == 2){
    #grid
    #NOTE - the grid will be used to evaluate the 2 densities (allows for smooth contours)
    x_grid <- seq(min(X[,1])-1, max(X[,1])+1, length.out = 200)
    y_grid <- seq(min(X[,2])-1, max(X[,2])+1, length.out = 200)
    grid <- expand.grid(X1 = x_grid, X2 = y_grid)
    
    k_true <- length(pi)
    k <- length(pi_est)
    
    #true density
    true_density <- numeric(nrow(grid))
  for(i in 1:k_true){
    
    if(type == "skew_normal"){
    true_density <- true_density + pi[i]*sn::dmsn(grid, xi = mu[i,], Omega = sigma[[i]], alpha = alpha[[i]])
  }
    else {  
    true_density <- true_density + pi[i]*mvtnorm::dmvnorm(grid, mean = mu[i,], sigma = sigma[[i]])
    }}

  #estimated gmm density
  #will be done using k (therefore, will be possible for ALL cases)
  est_density <- numeric(nrow(grid))
  for(j in 1:k){
    #NOTE - must extract the sigma correctly
    sigma_j <- if(is.list(sigma_est)) sigma_est[[j]] else sigma_est[j, , ]
    est_density <- est_density + pi_est[j]*mvtnorm::dmvnorm(grid, mean = mu_est[j,], sigma = sigma_j)
  }
  
  #plot
  plot(X, col = "blue", pch = 19, 
       xlab = "X1", ylab = "X2")
  contour(x = x_grid, y = y_grid, z = matrix(true_density, nrow = length(x_grid)),
          add = TRUE, drawlabels = FALSE, lwd = 2, lty = 2)
  contour(x = x_grid, y = y_grid, z = matrix(est_density, nrow = length(x_grid)),
          add = TRUE, drawlabels = FALSE, lwd = 2, lty = 1, col = "orange")
  #legend("bottomright", legend = c("Randomly simulated dataset", "True Density", "Monte Carlo Estimates"),
  #       col = c("blue", "black", "orange"), pch = c(19, NA, NA),
  #       lty = c(NA, 2, 1 ), cex = 0.5)
    }
  else{message("Plotting skipped: p not equal to 2.")}
}


#creating the plots for all 3 cases
#full model
plot_gmm_2d(mc_gmm$data$X, k = k, k_new = k_new, pi = pi, mu = mu, sigma = sigma, alpha = alpha,
            pi_est = mc_gmm$Mean$pi_hat, mu_est = mc_gmm$Mean$mu_hat, sigma_est = mc_gmm$Mean$sigma_hat, type = "guassian")

#cluster misspecification  
plot_gmm_2d(mc_cluster$data$X, k = k, k_new = k_new, pi = pi, mu = mu, sigma = sigma, alpha = alpha,
            pi_est = mc_cluster$Mean$pi_hat, mu_est = mc_cluster$Mean$mu_hat, sigma_est = mc_cluster$Mean$sigma_hat, type = "guassian")

#skew normal data
plot_gmm_2d(mc_skew$data$X, k = k, k_new = k_new, pi = pi, mu = mu, sigma = sigma, alpha = alpha,
            pi_est = mc_skew$Mean$pi_hat, mu_est = mc_skew$Mean$mu_hat, sigma_est = mc_skew$Mean$sigma_hat, type = "skew_normal")


#NOTE - this function does not account for label switching!
#However, since the focus for this is finding the loglik, AIC and BIC this does not matter
#Since: label switching does not effect the density! (and therefore won't effect these metrics)

#will now fit the skew normal mixture model onto the skew normal data that was generated
#this is to compare the AIC and BIC with what was found when fitting the GMM on this data
fit_sn <- smsn.mmix(y = mc_skew$data$X, g = k, family = "Skew.normal",
                   get.init = TRUE, group = TRUE, criteria = TRUE)


#extracting fit metrics
sn_loglik <- fit_sn$logLik
sn_AIC <- fit_sn$aic
sn_BIC <- fit_sn$bic


