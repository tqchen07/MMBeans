MMBeans <- function(X, PP, H, alpha, lam, mu0, sig0, dir_phi,
                    alpha_sigma, beta_sigma, EI, b, 
                    W=1, maxit=100, ELBO_tol = 1e-4){

  ## Args:
  # X : c(p, p, n, M) quantitative array
  # EI: c(n, H) matrix, initial start for subject cluster
  # PP: c(p, K, M) matrix, initial start for node membership
  # alpha: hyper-parameter for stick breaking w_l
  # mu0: mean for noisy features, usually set to 0
  # sig0: variance for noisy features
  # dir_phi: K-vector hyper-parameter for tau's dirichlet distribution
  # alpha_sigma: hyperparameter of inverse-gamma distribution for variance
  # beta_sigma: hyperparameter of inverse-gamma distribution for variance
  # W: optional weights for subjects, default is 1
  # maxit: maximum iterations
  # ELBO_tol: 
  
  t0 <- Sys.time()
  M <- dim(X)[4]
  n <- dim(X)[3]
  p <- dim(X)[2]
  K <- dim(PP)[2]
  
  Elgau <- array(0, c(n, K, K, H, M)) 
  Elgau0 <- array(0, c(n, K, K, H, M)) 
  
  now_cbase <- array(0, c(K, K, H, M))
  cbase_tmp <- array(0, c(K, K, H, M)) # same value across H 
  Egam <- 1/(1 + exp(-now_cbase)) # Eq[gamma_hts]
  
  d <- array(0, c(K, K, H, M))
  v <- array(0, c(K, K, H, M))
  m <- array(0, c(K, K, H, M))
  r <- array(0, c(K, K, H, M))
  
  Xsq <- X^2
  X_vec  <- array(X, dim = c(p*p, n, M)) 
  Xsq_vec <- array(Xsq, dim = c(p*p, n, M))
  mu0_sig0 <- mu0*mu0/sig0 + log(sig0)
  
  ELBO_iter0 <- 1
  ELBO_list <- c()
  ElogP_list <- c()
  ElogQ_list <- c()
  ElogP_X_Xi_list <- c()
  iter <- 0
  while (iter < maxit){
    iter = iter + 1
    cat("iter:", iter, " ")
    
    # preparation
    PPPP <- array(0, c(p*p, K, K, M))
    for(mi in 1:M){
      for(s in 1:K){
        for(t in 1:K){
          tmp <- PP[, s, mi] %*% t(PP[, t, mi])
          diag(tmp) <- 0
          PPPP[, s, t, mi] <- matrix(tmp, ncol = 1)
        }
      }
    }
    
    XsqPPsum <- XPPsum <- array(NA, dim = c(K, K, n, M))
    for(mi in 1:M){
      for(s in 1:K){
        for(t in 1:K){
          XPPsum[s, t, , mi] <- t(X_vec[, , mi]) %*% PPPP[, s, t, mi]
          
          XsqPPsum[s, t, , mi] <- t(Xsq_vec[, , mi]) %*% PPPP[, s, t, mi]
        }
      }
    }
    
    
    # f, g: beta distribution
    sEI <- apply(EI, 2, sum) # sum_1^n E[I(Zi=h)]
    
    f <- 1 + sEI # f_h: H-vector
    tmp <- diffinv(-sEI)
    g <- alpha + tmp[-1] - tmp[H + 1] # g_h: H-vector
    
    tmp <- digamma(f + g)
    Elw <- digamma(f) - tmp # Eq[log(Wh)]
    El1w <- digamma(g) - tmp # Eq[log(1-Wh)]
    Elw[H] <- 0 # for h>=H, wh=1, thus E[log(wh)] = 0
    El1w[H] <- -Inf # E[log(1-wh)] = -inf
    cumEl1w <- diffinv(El1w) #length = H+1; sum_{l=1}^{h-1}Eq[log(1-wl)], cumEl1w[1]=0 
    
    EIW <- EI * W
    sEIW <- apply(EIW, 2, sum) # sum_i=1^n E[I(Zi=h)]*Wi
    
    # d, v, m, r: Normal-inverse-gamma
    forc1 <- array(0, dim = c(K, K, H, M))
    forc2 <- array(0, dim = c(K, K, H, M))
    forc3 <- array(0, dim = c(K, K, H, M))
    
    tmpsum1 <- 0
    tmpsum2 <- 0
    tmpsum3 <- 0
    for(mi in 1:M){
      for(h in 1:H){
        for(s in 1:K){
          for(t in 1:K){
            tmpsum1 <- sum(PPPP[, s, t, mi])
            tmpsum2 <- XPPsum[s, t, , mi] %*%  EIW[, h]
            tmpsum3 <- XsqPPsum[s, t, , mi] %*%  EIW[, h]
            
            forc1[s, t, h, mi] <- tmpsum1 * sEIW[h] 
            forc2[s, t, h, mi] <- tmpsum2
            forc3[s, t, h, mi] <- tmpsum3
          }
        }
      }
    }
    d <- alpha_sigma + forc1*Egam
    v <- lam + forc1*Egam
    m <- forc2/v * Egam
    r <- beta_sigma - v*m*m + forc3*Egam
    Elsig <- log(r/2) - digamma(d/2)
    
    
    # cbase: feature selection
    Elgau  <- array(0, c(n, K, K, H, M)) 
    Elgau0 <- array(0, c(n, K, K, H, M)) 
    
    d_r <- d/r
    dm_r <- d*m/r
    dm2_r_v_Elsig <- d*m*m/r + 1/v + Elsig
    
    for(mi in 1:M){
      for(h in 1:H){
        for(s in 1:K){
          for(t in 1:K){
            tmp1sum <- - 0.5*XsqPPsum[s, t, , mi] * d_r[s, t, h, mi] + XPPsum[s, t, , mi] * dm_r[s, t, h, mi] # n
            tmp2sum <- sum(- 0.5*PPPP[, s, t, mi] * dm2_r_v_Elsig[s, t, h, mi]) # dim = p*p => 1
            Elgau[ , s, t, h, mi] <- tmp1sum + tmp2sum
            
            tmp3sum <- 0.5*XsqPPsum[s, t, , mi] * (1/sig0) - XPPsum[s, t, , mi] * mu0/sig0 # dim = n
            tmp4sum <- sum(0.5*PPPP[, s, t, mi] * mu0_sig0) # dim = p*p => 1
            Elgau0[, s, t, h, mi] <- tmp3sum + tmp4sum
          }# for t
        }# for s
      }# for h
    }# for mi
    
    for(mi in 1:M){
      for(s in 1:K){
        for(t in 1:K){
          for(h in 1:H ){
            cbase_tmp[s, t, h, mi] <- t(Elgau[, s, t, h, mi]) %*% EIW[, h]
            cbase_tmp[s, t, h, mi] <- cbase_tmp[s, t, h, mi] + t(Elgau0[,s,t,h,mi]) %*% EIW[, h]
          }
        }
      }
    }
    
    cbase <- apply(cbase_tmp, c(1, 2, 4), sum)
    cbase_arr <- array(0, dim = c(K, K, H, M))
    for(h in 1:H){
      cbase_arr[, , h, ] <- cbase
    }
    
    old_cbase <- now_cbase
    now_cbase <- cbase_arr
    Egam <- 1/(1 + exp(-cbase_arr))
    
    # b: subtyping
    old_b <- b
    for(h in 1:H ){
      b[, h] <- matrix(Elgau[, , , h, ], nrow = n) %*% c(Egam[, , h, ])
      b[, h] <- b[,h] - matrix(Elgau0[ , , ,h, ], nrow=n) %*% c(1-Egam[ , , h, ])
      b[, h] <- b[,h] + Elw[h] + cumEl1w[h]
    }
    b <- b - apply(b, 1, max) # log-sum-exp trick
    eb <- exp(b - apply(b, 1, max))
    EI <- eb/apply(eb, 1, sum)
    
    EIW_tmp <- EI * W
    sEIW_tmp <- apply(EIW_tmp, 2, sum)
    
    # PP: modular structure
    dir_phi_ti <- matrix(NA, nrow = K, ncol = M)
    Eltao <- matrix(NA, nrow = K, ncol = M)
    for(mi in 1:M){
      dir_phi_ti[, mi] <- apply(PP[, , mi], 2, sum) + dir_phi
      Eltao[, mi] <- digamma(dir_phi_ti[, mi]) - digamma(sum(dir_phi_ti[, mi]) ) # E[log(tao_t)]
    }
    
    old_PP <- PP
    new_PP <- array(NA, dim = c(p, K, M))
    
    Egam_dm2_r_v_Elsig <- Egam * dm2_r_v_Elsig
    Egam_d_r <- Egam * d_r
    Egam_dm_r <- Egam * dm_r
    
    for(mi in 1:M){
      for(t in 1:K){
        tmp1 <- PP[, , mi] %*% Egam_dm2_r_v_Elsig[, t, , mi] %*% sEIW_tmp #p*1
        tmp4 <- PP[, , mi] %*% (1- Egam[, t, , mi]) %*% sEIW_tmp * mu0_sig0
        for(j in 1:p){
          tmp2 <- PP[, , mi] * Xsq[j, , , mi] %*% EIW_tmp %*% t(Egam_d_r[, t, , mi])
          tmp3 <- -2 * PP[, , mi] * X[j, , , mi] %*% EIW_tmp %*% t(Egam_dm_r[, t, , mi])
          
          tmp1sum <- sum(tmp1) - sum(tmp1[j])
          tmp2sum <- sum(tmp2) - sum(tmp2[j, ]) 
          tmp3sum <- sum(tmp3) - sum(tmp3[j, ])
          
          tmp5 <- PP[, , mi] * Xsq[j, , , mi] %*% EIW_tmp %*% t(1- Egam[, t, , mi])/sig0
          tmp6 <- -2 * PP[, , mi] * X[j, , , mi] %*% EIW_tmp %*% t(1- Egam[, t, , mi])*mu0/sig0
          
          tmp4sum <- sum(tmp4) - sum(tmp4[j])
          tmp5sum <- sum(tmp5) - sum(tmp5[j, ])
          tmp6sum <- sum(tmp6) - sum(tmp6[j, ])
          new_PP[j, t, mi] <- -0.5*(tmp1sum + tmp2sum + tmp3sum)-0.5*(tmp4sum + tmp5sum + tmp6sum)
        }
        new_PP[, t, mi] <- new_PP[, t, mi] + Eltao[t, mi]
      }
    }
    
    for(mi in 1:M){
      ePPmi <- exp(new_PP[, , mi] - apply(new_PP[, , mi], 1, max))
      PP[, , mi] <- ePPmi/apply(ePPmi, 1, sum)
    }
    
    #---------------------------# 
    # ElogP calculation
    Elp_PPPP <- array(0, c(p*p, K, K, M))
    for(mi in 1:M){
      for(s in 1:K){
        for(t in 1:K){
          tmp <- PP[, s, mi] %*% t(PP[, t, mi])
          diag(tmp) <- 0
          Elp_PPPP[, s, t, mi] <- matrix(tmp, ncol = 1)
        }
      }
    }
    
    Elp_XsqPPsum <- Elp_XPPsum <- array(NA, dim = c(K, K, n, M))
    for(mi in 1:M){
      for(s in 1:K){
        for(t in 1:K){
          Elp_XPPsum[s, t, , mi]  <- t(X_vec[, , mi]) %*% Elp_PPPP[, s, t, mi]
          
          Elp_XsqPPsum[s, t, , mi] <- t(Xsq_vec[, , mi]) %*% Elp_PPPP[, s, t, mi]
        }
      }
    }
    
    EIW <- EI*W
    ElogP_X <- 0 # Eq(logP(X|Xi)), X is data
    for(mi in 1:M){
      # (n, H) %*% (H, K*K) %*% (K*K, 1)
      tmp1 <- sum(EIW %*% t(matrix( -dm2_r_v_Elsig[, , , mi]/2, ncol = H)) %*% matrix(colSums(Elp_PPPP[, , , mi])*Egam[, ,1,mi], ncol = 1))
      tmp2 <- sum(EIW * sum(colSums(Elp_PPPP[, , , mi])*(1-Egam[, , 1, mi])) * (-log(sig0)/2) - mu0^2/sig0/2)
      
      # (1, K*K) %*% (K*K, n) %*% (n, H)
      tmp3 <- sum(t(matrix((1-Egam[, , 1, mi]), ncol = 1)) %*% matrix(-Elp_XsqPPsum[, , , mi]/sig0/2, ncol = n) %*% EIW)
      tmp4 <- sum(t(matrix((1-Egam[, , 1, mi]), ncol = 1)) %*% matrix(Elp_XPPsum[, , , mi]*mu0/sig0, ncol = n) %*% EIW)
      
      # (K*K, n) %*% (n, H) * (K*K, H) * (K*K, H)
      tmp5 <- sum(matrix(Elp_XsqPPsum[, , , mi], ncol = n) %*% EIW * (matrix((-d[, , , mi]/r[, , , mi]/2), ncol = H)) * matrix(Egam[, , 1, mi], ncol = H, nrow = K*K))
      # (K*K, n) %*% (n, H) * (K*K, H) * (K*K, H)
      tmp6 <- sum(matrix(Elp_XPPsum[, , , mi], ncol = n) %*% EIW * (matrix((d[, , , mi]*m[, , , mi]/r[, , , mi]), ncol = H)) * matrix(Egam[, , 1, mi], ncol = H, nrow = K*K))

      tmp7 <- sum(EIW * sum(Elp_PPPP[, , , mi]) * log(1/sqrt(2*pi)))
      
      ElogP_X <- ElogP_X + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6 + tmp7
    }
    ElogP_Xi_Z <- sum(EIW %*% (Elw + cumEl1w[-H-1])) # Eq[logP(Z|w)]
    ElogP_Xi_G <- 0
    for(j in 1:p){ElogP_Xi_G <- ElogP_Xi_G + sum(PP[j, , ]*Eltao)} # Eq[logP(G|tau)]
    ElogP_Xi_musig <- sum(-(3/2+alpha_sigma/2)*Elsig-beta_sigma/2*d/r-lam/2*(1/v + m^2*d/r)) # Eq[logPmusig]
    ElogP_Xi_w <- sum((alpha - 1) * El1w[-H]) # Eq[logP(w)]
    ElogP_Xi_phi <- sum((dir_phi-1) * Eltao) # Eq[logP(phi|tau)]
    
    ElogP_Xi <- ElogP_Xi_Z + ElogP_Xi_G + ElogP_Xi_musig + ElogP_Xi_w + ElogP_Xi_phi
    ElogP_X_Xi <- ElogP_X + ElogP_Xi
    
    ElogQ_Z <- sum(EIW * sweep(b, 1, log(apply(exp(b), 1, sum)), `-`))
    ElogQ_G <- sum(PP * log(PP), na.rm = T)
    ElogQ_musig <- sum(1/2*log(v) + d/2*log(r/2)-lgamma(d/2)-(3/2 + d/2)*Elsig - d/2-1/2)
    ElogQ_w <- sum((f-1)*Elw + (g-1)*El1w - lgamma(f) - lgamma(g) + lgamma(f+g), na.rm = T)
    ElogQ_phi <- 0
    for(mi in 1:M){ElogQ_phi <- ElogQ_phi + sum((dir_phi_ti[,mi] - 1)*Eltao[,mi])*(-sum(lgamma(dir_phi_ti[,mi])) + lgamma(sum(dir_phi_ti[,mi])))}
    ElogQ_gamma <- 0
    if(any(exp(cbase) == Inf)){
      ElogQ_gamma <- sum(Egam[, , 1, ]*cbase - cbase)
    }else{
      ElogQ_gamma <- sum(Egam[, , 1, ]*cbase - log(1 + exp(cbase)))
    }
    ElogQ <- ElogQ_Z + ElogQ_G + ElogQ_musig + ElogQ_w + ElogQ_phi + ElogQ_gamma
    
    ElogP_list <- c(ElogP_list, ElogP_X)
    ElogP_X_Xi_list <- c(ElogP_X_Xi_list, ElogP_X_Xi)
    ElogQ_list <- c(ElogQ_list, ElogQ)
    ELBO_list <- c(ELBO_list, ElogP_X_Xi - ElogQ)
    
    cat("ELBO:", ElogP_X_Xi - ElogQ, "\n")
    
    # convergence criteria
    ## criteria 1: relative absolute difference is less than 1e-4
    abs_lift <- abs((c(ELBO_iter0, ELBO_list)[iter+1] - c(ELBO_iter0, ELBO_list)[iter])/c(ELBO_iter0, ELBO_list)[iter]) 
    if ((abs_lift <= 1e-4) & (abs_lift > 0)) {
      cat("Convergence criteria meets (abs lift < 1e-4)!\n")
      break 
    }
    ## criteria 2: absolute difference is not changing for 3 iterations.
    if (iter > 3){
      osc_seq <- tail(abs(diff(c(ELBO_iter0, ELBO_list))), 3)
      if (all(diff(osc_seq) == 0)) { 
        cat("Convergence criteria meets (stable abs differences)!\n")
        break 
      }
    }
    
  } # end while true 
  if (iter == maxit) {cat("Max iteration is reached!\n")}
  
  comp_time <- difftime(Sys.time(), t0, units = "secs")
  
  VBIC <- -2*ElogP_X + 2*ElogQ
  
  #---------------------------# 
  list(iter = iter, EI = EI, b = b, PP = PP,
       ElogP = ElogP_X, ElogQ = ElogQ, VBIC = VBIC,
       ElogP_list = ElogP_list, ElogP_X_Xi_list = ElogP_X_Xi_list,
       ElogQ_list = ElogQ_list, ELBO_list = ELBO_list,
       cbase = cbase, Egam = Egam[, , 1, ], W = W,
       f = f, g = g, Elw = Elw, El1w = El1w, d = d, m = m,  r = r, v = v, 
       Elsig = Elsig, Eltao = Eltao, dir_phi_ti = dir_phi_ti,
       comp_time = comp_time)
  
} #end function
