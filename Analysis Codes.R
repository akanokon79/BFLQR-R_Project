library(MASS)

# Reproducibility
set.seed(2026)

# =========================================
# BFLQR Gibbs Sampler
# =========================================

bflqr_gibbs <- function(y, X, tau,
                        niter=6000,
                        burnin=2000,
                        tau_beta=5,
                        tau_fuse=5){
  
  n <- length(y)
  p <- ncol(X)
  
  theta <- (1-2*tau)/(tau*(1-tau))
  
  beta <- rep(0,p)
  w <- rep(1,n)
  
  lambda_beta <- rep(1,p)
  lambda_fuse <- rep(1,p-1)
  
  D <- diff(diag(p))
  
  beta_store <- matrix(0,niter-burnin,p)
  
  for(iter in 1:niter){
    
    resid <- as.numeric(y - X %*% beta)
    
    w <- rexp(n, rate=1/(abs(resid)+1e-6))
    
    z <- theta*w
    
    W <- diag(1/(w+1e-6))
    
    V_inv <- t(X)%*%W%*%X +
      diag(1/(tau_beta*lambda_beta+1e-6)) +
      t(D)%*%diag(1/(tau_fuse*lambda_fuse+1e-6))%*%D +
      diag(1e-6,p)
    
    V <- solve(V_inv)
    
    mu <- V %*% t(X)%*%W%*%(y-z)
    
    beta <- as.numeric(mvrnorm(1,mu,V))
    
    lambda_beta <- 1/rgamma(
      p,1,abs(beta)+1e-6
    )
    
    if(p>1){
      
      diff_beta <- diff(beta)
      
      lambda_fuse <- 1/rgamma(
        p-1,1,
        abs(diff_beta)+1e-6
      )
    }
    
    if(iter > burnin)
      beta_store[iter-burnin,] <- beta
  }
  
  list(
    beta_mean = colMeans(beta_store),
    
    beta_ci = t(apply(
      beta_store,2,
      quantile,c(0.05,0.95)
    )),
    
    pip = colMeans(
      abs(beta_store)>0.2
    )
  )
}

# =========================================
# Load Data
# =========================================

ak <- na.omit(read.csv(
  "C:/Users/AKANINYENE O. EKONG/OneDrive/Desktop/emissionsdata.csv"
))

ak <- data.frame(lapply(
  ak[,c("Y",
        "X1","X2","X3","X4","X5",
        "X6","X7","X8","X9","X10")],
  as.numeric
))

y <- ak$Y

X <- scale(as.matrix(
  ak[,c("X1","X2","X3","X4","X5",
        "X6","X7","X8","X9","X10")]
))

taus <- c(0.25,0.50,0.75)

# =========================================
# Main Analysis
# =========================================

for(tau in taus){
  
  cat("\n============================\n")
  cat("Quantile =",tau,"\n")
  cat("============================\n")
  
  fit <- bflqr_gibbs(
    y,X,tau,
    tau_beta=5,
    tau_fuse=5
  )
  
  res <- data.frame(
    Variable = paste0("X",1:ncol(X)),
    beta_mean = round(fit$beta_mean,3),
    pip = round(fit$pip,3),
    ci05 = round(fit$beta_ci[,1],3),
    ci95 = round(fit$beta_ci[,2],3)
  )
  
  res$selected <- ifelse(
    res$pip > 0.7,1,0
  )
  
  print(res)
  
  # Prediction Performance
  
  yhat <- X %*% fit$beta_mean
  
  residuals <- y - yhat
  
  mse <- mean(residuals^2)
  
  rmse <- sqrt(mse)
  
  cat("\nMSE  =",round(mse,5))
  cat("\nRMSE =",round(rmse,5),"\n")
}