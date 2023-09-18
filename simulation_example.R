#### simulation ####
# 1. two states
# 2. three subtypes
# 3. Each state has 3*3 modular structure

rm(list=ls())

#### helper functions ####
mat_to_vec_i <- function(mat){
  # convert matrix to vector with lower triangle (no diagonal)
  # input
  # mat: p*p array
  p <- dim(mat)[1]
  vec <- matrix(0, ncol = 1, nrow = p*(p-1)/2)
  t <- 0
  for(j in 2:p){
    for(k in seq_len(j-1)){
      t <- t+1
      vec[t, 1] <- mat[j, k]
    }
  }
  
  return(vec)
}

cbase_block_to_node <- function(node_group, cbase){
  mat <- matrix(NA, nrow = length(node_group), ncol = length(node_group))
  for(i in 1:length(node_group)){
    for(j in 1:length(node_group)){
      mat[i, j] <- cbase[node_group[i], node_group[j]]
    }
  }
  return(mat)
}


#### load function to be tested ####
setwd("/home/tc732/project/ABCD/github/")
source("./MMBeans_function.R")

#### generate data ####
n <- 100 # number of sample size
p <- 60 # number of nodes
K <- 3 # number of modules
H <- 3 # number of subtypes
M <- 2 # number of states

#### true values
# noise level
mu0 <- 0
sig0 <- sqrt(100)

# signal level
mu_seq <- c(-5, 0, 5, 
            0, 5, -5, 
            5, -5, 0, 
            5, -5, 0,
            -5, 0, 5, 
            0, 5, -5) + 2

sig_seq <- c(7, 5, 3,
             3, 7, 5,
             5, 3, 7,
             7, 5, 3,
             3, 7, 5,
             5, 3, 7)

# Z: subject clusterG1
set.seed(98765); Z <- sample(seq(1, H, 1), size = n, replace = T)

# G: block membership
G1 <- c(rep(1, 0.25*p), rep(2, 0.4*p), rep(3, 0.35*p))
G2 <- c(rep(1, 0.30*p), rep(2, 0.30*p), rep(3, 0.40*p))


#### simulated values
# gamma: informative block selection
mat_gamma1 <- matrix(NA, nrow = p, ncol = p)
mat_gamma2 <- matrix(NA, nrow = p, ncol = p)

set.seed(1837502)
# state 1
X1 <- array(NA, c(p, p, n))
# (1, 1)-th block
for(j in which(G1==1)){
  for(k in which(G1==1)){
    X1[j,k,] <- rnorm(n, mean = mu0, sd = sig0)
    mat_gamma1[j, k] <- 0
  }
}

# (1, 2)-th block
for(j in which(G1==1)){
  for(k in which(G1==2)){
    X1[k,j,Z==1] <- X1[j,k,Z==1] <- rnorm(sum(Z==1), mean = mu_seq[1], sd = sig_seq[1])
    X1[k,j,Z==2] <- X1[j,k,Z==2] <- rnorm(sum(Z==2), mean = mu_seq[2], sd = sig_seq[2])
    X1[k,j,Z==3] <- X1[j,k,Z==3] <- rnorm(sum(Z==3), mean = mu_seq[3], sd = sig_seq[3])
    mat_gamma1[k, j] <- mat_gamma1[j, k] <- 1
  }
}

# (1, 3)-th block
for(j in which(G1==1)){
  for(k in which(G1==3)){
    X1[k,j,] <- X1[j,k,] <- rnorm(n, mean = mu0, sd = sig0)
    mat_gamma1[k, j] <- mat_gamma1[j, k] <- 0
  }
}

# (2, 2)-th block
for(j in which(G1==2)){
  for(k in which(G1==2)){
    X1[j,k,Z==1] <- rnorm(sum(Z==1), mean = mu_seq[4], sd = sig_seq[4])
    X1[j,k,Z==2] <- rnorm(sum(Z==2), mean = mu_seq[5], sd = sig_seq[5])
    X1[j,k,Z==3] <- rnorm(sum(Z==3), mean = mu_seq[6], sd = sig_seq[6])
    mat_gamma1[j, k] <- 1
  }
}

# (2, 3)-th block
for(j in which(G1==2)){
  for(k in which(G1==3)){
    X1[k,j,] <- X1[j,k,] <- rnorm(n, mean = mu0, sd = sig0)
    mat_gamma1[k, j] <- mat_gamma1[j, k] <- 0
  }
}

# (3, 3)-th block
for(j in which(G1==3)){
  for(k in which(G1==3)){
    X1[j,k,Z==1] <- rnorm(sum(Z==1), mean = mu_seq[7], sd = sig_seq[7])
    X1[j,k,Z==2] <- rnorm(sum(Z==2), mean = mu_seq[8], sd = sig_seq[8])
    X1[j,k,Z==3] <- rnorm(sum(Z==3), mean = mu_seq[9], sd = sig_seq[9])
    mat_gamma1[j, k] <- 1
  }
}

for(i in 1:n){
  X1[,,i][lower.tri(X1[,,i], diag = TRUE)] <- 0
  X1[,,i] <- X1[,,i] + t(X1[,,i])
  if(!isSymmetric(X1[, , i])){"Stop! Matrix is not symmetric"}
  
  mat_gamma1[lower.tri(mat_gamma1, diag = TRUE)] <- 0
  mat_gamma1 <- mat_gamma1 + t(mat_gamma1)
  if(!isSymmetric(mat_gamma1)){"Stop! mat_gamma1 Matrix is not symmetric"}
} 

# state 2:
X2 <- array(NA, c(p, p, n))
# (1, 1)-th block
for(j in which(G2==1)){
  for(k in which(G2==1)){
    X2[j,k,] <- rnorm(n, mean = mu0, sd = sig0)
    mat_gamma2[j, k] <- 0
  }
}

# (1, 2)-th block
for(j in which(G2==1)){
  for(k in which(G2==2)){
    X2[k,j,Z==1] <- X2[j,k,Z==1] <- rnorm(sum(Z==1), mean = mu_seq[10], sd = sig_seq[10])
    X2[k,j,Z==2] <- X2[j,k,Z==2] <- rnorm(sum(Z==2), mean = mu_seq[11], sd = sig_seq[11])
    X2[k,j,Z==3] <- X2[j,k,Z==3] <- rnorm(sum(Z==3), mean = mu_seq[12], sd = sig_seq[12])
    mat_gamma2[k, j] <- mat_gamma2[j, k] <- 1
  }
}

# (1, 3)-th block
for(j in which(G2==1)){
  for(k in which(G2==3)){
    X2[k,j,Z==1] <- X2[j,k,Z==1] <- rnorm(sum(Z==1), mean = mu_seq[13], sd = sig_seq[13])
    X2[k,j,Z==2] <- X2[j,k,Z==2] <- rnorm(sum(Z==2), mean = mu_seq[14], sd = sig_seq[14])
    X2[k,j,Z==3] <- X2[j,k,Z==3] <- rnorm(sum(Z==3), mean = mu_seq[15], sd = sig_seq[15])
    mat_gamma2[k, j] <- mat_gamma2[j, k] <- 1
  }
}

# (2, 2)-th block
for(j in which(G2==2)){
  for(k in which(G2==2)){
    X2[j,k,] <- rnorm(n, mean = mu0, sd = sig0)
    mat_gamma2[j, k] <- 0
  }
}

# (2, 3)-th block
for(j in which(G2==2)){
  for(k in which(G2==3)){
    X2[k,j,Z==1] <- X2[j,k,Z==1] <- rnorm(sum(Z==1), mean = mu_seq[16], sd = sig_seq[16])
    X2[k,j,Z==2] <- X2[j,k,Z==2] <- rnorm(sum(Z==2), mean = mu_seq[17], sd = sig_seq[17])
    X2[k,j,Z==3] <- X2[j,k,Z==3] <- rnorm(sum(Z==3), mean = mu_seq[18], sd = sig_seq[18])
    mat_gamma2[k, j] <- mat_gamma2[j, k] <- 1
  }
}

# (3, 3)-th block
for(j in which(G2==3)){
  for(k in which(G2==3)){
    X2[j,k,] <- rnorm(n, mean = mu0, sd = sig0)
    mat_gamma2[j, k] <- 0
  }
}

for(i in 1:n){
  X2[,,i][lower.tri(X2[,,i], diag = TRUE)] <- 0
  X2[,,i] <- X2[,,i] + t(X2[,,i])
  if(!isSymmetric(X2[, , i])){"Stop! Matrix is not symmetric"}
  
  mat_gamma2[lower.tri(mat_gamma2, diag = TRUE)] <- 0
  mat_gamma2 <- mat_gamma2 + t(mat_gamma2)
  if(!isSymmetric(mat_gamma2)){"Stop! mat_gamma2 Matrix is not symmetric"}
} 

vec_gamma1 <- mat_to_vec_i(mat_gamma1)
vec_gamma2 <- mat_to_vec_i(mat_gamma2)

multi_X <- array(0, dim = c(p, p, n, M))
multi_X[, , , 1] <- X1
multi_X[, , , 2] <- X2


#### MMBeans ####
cat("MMBeans starts!\n")

# initial: EI
rZ <- sample(seq(1, H+2), size = n, replace = T)
rH <- length(unique(rZ))
b <- matrix(rnorm(n*rH, 0, 0.1), n, rH)
for(h in 1:rH){
  b[rZ==h, h] <- 0.6
}
eb <- exp(b)
EI <- eb/apply(eb,1,sum)

# initial: PP
VBIC_vec <- c()
rslt_list <- list()
i <- 0
for(PP1_inc in c(1, 2)){
  for(PP2_inc in c(1, 2)){
    i <- i+1
    rG1 <- sample(seq(1, K+PP1_inc), size = n, replace = T)
    rG2 <- sample(seq(1, K+PP2_inc), size = n, replace = T)
    rK <-  max(max(rG1), max(rG2))
    
    PP <- array(rnorm(rK*p*M, 0, 0.1), dim = c(p, rK, M))
    for(j in 1:p){
      PP[j, rG1[j], 1] <- 0.6
      PP[j, rG2[j], 2] <- 0.6
    }
    
    for(mi in 1:M){
      ep <- exp(PP[, , mi])
      PP[, , mi] <- ep/apply(ep, 1, sum)
    }
    
    
    rslt_tmp <- MMBeans(multi_X, PP, H=rH, 
                    alpha = 1, lam = 1,
                    mu0 = 0, sig0 = var(multi_X), 
                    dir_phi = rep(1, rK),
                    alpha_sigma = 0.02, beta_sigma = 0.02, 
                    EI, b, W = 1, maxit = 100)
    
    VBIC_vec <- c(VBIC_vec, rslt_tmp$VBIC)
    rslt_list[[i]] <- rslt_tmp
  }
}
rslt <- rslt_list[[which.max(VBIC_vec)]]

#### simulation result ####
Z_MMBeans <- apply(rslt$EI, 1, which.max)
G1_MMBeans <- apply(rslt$PP[, , 1], 1, which.max)
G2_MMBeans <- apply(rslt$PP[, , 2], 1, which.max)


vec_cbase1 <- mat_to_vec_i(cbase_block_to_node(node_group=G1_MMBeans, 
                                               cbase = rslt$cbase[, , 1]))
vec_cbase2 <- mat_to_vec_i(cbase_block_to_node(node_group=G2_MMBeans, 
                                               cbase = rslt$cbase[, , 2]))

gamma_true <- factor(c(vec_gamma1, vec_gamma2), levels = c(1, 0))

gamma_MMBeans <- factor(c(as.numeric(vec_cbase1 > median(unique(vec_cbase1))), 
                          as.numeric(vec_cbase2 > median(unique(vec_cbase2)))), 
                        levels = c(1, 0))
tab_gamma <- caret::confusionMatrix(table(gamma_MMBeans, gamma_true))

roc_cbase <- pROC::roc(c(vec_gamma1, vec_gamma2), c(vec_cbase1, vec_cbase2), direction = "<")


cat("Subtyping ARI: ", mclust::adjustedRandIndex(Z_MMBeans, Z), "\n")
cat("Modular structure: \n")
cat("M1 ARI: ", mclust::adjustedRandIndex(G1_MMBeans, G1), "\n")
cat("M2 ARI: ", mclust::adjustedRandIndex(G2_MMBeans, G2), "\n")
cat("Feature selection: \n")
cat("Sensitivity: ", tab_gamma$byClass["Sensitivity"], "\n")
cat("Specificity: ", tab_gamma$byClass["Specificity"], "\n")
cat("Y-index: ", tab_gamma$byClass["Sensitivity"] + tab_gamma$byClass["Specificity"] - 1, "\n")
cat("AUC: ", pROC::auc(roc_cbase), "\n")