KF_BSd_3F_uKD <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){
    
    r_1 <- log(r[1])
    r_2 <- log(r[2])
    r_c <- log(r[3])
    
    n_factors <- length(kappa)  
    
    n_ages <- nrow(mu_bar)   # - Number of ages
    n_years <- ncol(mu_bar)  # - Number of years
    
    delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix
    
    v_ti <- mu_bar
    F_ti <- mu_bar
    
    # - Def. variables
    ## - State variables
    X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
    X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)
    
    S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
    S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors
    
    ## - Factor loading matrices
    A_tT <- matrix(0, n_ages, 1)
    B_tT <- matrix(NA, n_ages, n_factors)
    
    # - Initial values of states and of covariance
    X_t_c[,1] <- x0
    X_t[,1] <- x0
    
    # - Initialize X and Sigma
    x_ti <- x0 #init_X
    P_ti <- diag(1, n_factors) * 1e-10
    
    S_t_c[,(1:n_factors)] <- diag(1, n_factors) * 1e-10
    S_t[,(1:n_factors)] <- diag(1, n_factors) * 1e-10
    
    R <- matrix(0, n_factors, n_factors) # - Factor covariance
    
    Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
    
    # - Build diffusion process
    ## - Build lower cholesky factor
    dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
    odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
    Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
    diag(Low_chol) <- exp(dg_l_Sigma_chol)
    
    # - Get Sigma (covariance matrix of the diffusion process)
    Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)
    
    # - Get R (covariance of the state variable)
    for(row in 1:n_factors){
        for(col in 1:n_factors){
            R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
        }
    }
    
    for(age in 1:n_ages){    # - scroll over the ages
        A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)  
        B_tT[age,] <- B_BSd_3F(age, delta_matrix)  ###### B_ind(age,delta)  
    }
    
    for(t in 1:n_years){
        
        # - First observation
        x_ti <- Phi %*% x_ti    # - x_{1,t}
        P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
        v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
        F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
        #    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + max(1e-8, r_c + exp(r_1) * exp(exp(r_2))) #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
        
        X_t_c[,t+1] <- x_ti
        S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti
        
        for(i in 2:n_ages){
            x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}
            
            # - Joseph formula univariate
            K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
            P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)
            
            # - log-likelihood values from Koopman and Durbin
            F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i 
            v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti
            
        }
        
        X_t[,t+1] <- x_ti
        S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti
    }
    
    return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}
# - Fill lower triangular matrix (by row) from a vector: this turns out useful when coding delta and Sigma
# - in the dependent AFNS and Blackburn-Sherris model.
low_trg_fill <- function(par_vector){
    n_factors <- (-1 + sqrt(1 + 8 * length(par_vector))) * 0.5
    par_matrix <- matrix(0, n_factors, n_factors) # - lower triangular matrix of delta coefficients
    count_vec <- 1
    for(row in 1:n_factors){
        for(col in 1:row){
            par_matrix[row,col] <- par_vector[count_vec]
            count_vec <- count_vec + 1
        }
    }
    return(par_matrix)
}
cov2par <- function(Sigma_el){
    n_factors <- (-1 + sqrt(1 + 8 * length(Sigma_el))) * 0.5
    
    Sigma_diffusion <- matrix(diag(Sigma_el[1:n_factors]), n_factors, n_factors)
    Sigma_diffusion <- Sigma_diffusion + low_trg_fill_0diag(Sigma_el[(n_factors+1):length(Sigma_el)]) + t(low_trg_fill_0diag(Sigma_el[(n_factors+1):length(Sigma_el)]))
    
    Low_Sigma <- t(chol(Sigma_diffusion))
    odg_Sigma_chol <- rep(0, length(Sigma_el) - n_factors)
    
    count_vec <- 1
    for(row in 2:n_factors){
        for(col in 1:(row-1)){
            odg_Sigma_chol[count_vec] <- Low_Sigma[row,col]
            count_vec <- count_vec + 1
        }
    }
    
    return(list(dg_l_Sigma_chol=log(diag(Low_Sigma)), odg_Sigma_chol=odg_Sigma_chol))
}
# - Fill lower triangular matrix with 0s on the main diagonal (this is also useful when coding Sigma)
low_trg_fill_0diag <- function(par_vector){
    n_factors <- (1 + sqrt(1 + 8 * length(par_vector))) * 0.5
    par_matrix <- matrix(0, n_factors, n_factors) # - lower triangular matrix of delta coefficients
    count_vec <- 1
    for(row in 2:n_factors){
        for(col in 1:(row-1)){
            par_matrix[row,col] <- par_vector[count_vec]
            count_vec <- count_vec + 1
        }
    }
    return(par_matrix)
}
#### -  A(t,T) function for dependent model (already divided by T-t, and coded as if t=0)
A_BSd_3F <- function(T, sigma_mat, delta_mat){   # - function of matrices delta and sigma which are lower diagonal
    
    # D terms as in Huang et. al. 2020
    D1 <- delta_mat[2,1] / (delta_mat[1,1] - delta_mat[2,2])
    D2 <- delta_mat[3,2] / (delta_mat[2,2] - delta_mat[3,3])
    D3 <- delta_mat[2,1] / (delta_mat[1,1] - delta_mat[3,3])
    D4 <- delta_mat[3,1] / (delta_mat[1,1] - delta_mat[3,3])
    D5 <- delta_mat[3,2] / (delta_mat[1,1] - delta_mat[3,3])
    
    # E terms as in Huang et. al. 2020
    E1 <- 1 + D1 + D1 * D5 + D4
    E2 <- D1 * (1 + D2)
    E3 <- D2 * D3 - D4
    
    # F terms as in Huang et. al. 2020
    F1 <- (sigma_mat[1,1] ^ 2) * (E1 ^2)
    F2 <- (sigma_mat[1,1] ^ 2) * (E2 ^2) - 2 * sigma_mat[1,1] * sigma_mat[2,1] * E2 * (1 + D2) + ((sigma_mat[2,1] ^ 2) + (sigma_mat[2,2] ^ 2)) * ((1 + D2) ^ 2)  
    F3 <- (sigma_mat[1,1] ^ 2) * (E3 ^2) - 2 * sigma_mat[1,1] * sigma_mat[2,1] * E3 * D2 + ((sigma_mat[2,1] ^ 2) + (sigma_mat[2,2] ^ 2)) * (D2 ^ 2) - 2 * (sigma_mat[2,1] * sigma_mat[3,1] + sigma_mat[2,2] * sigma_mat[3,2]) * D2 + 2 * sigma_mat[1,1] * sigma_mat[3,1] * E3 + ((sigma_mat[3,1] ^ 2) + (sigma_mat[3,2] ^ 2) + (sigma_mat[3,3] ^ 2))
    F4 <- - 2 * (sigma_mat[1,1] ^ 2) * E1 * E2 + 2 * sigma_mat[1,1] * sigma_mat[2,1] * E1 * (1 + D2)
    F5 <-   2 * (sigma_mat[1,1] ^ 2) * E1 * E3 - 2 * sigma_mat[1,1] * sigma_mat[2,1] * E1 * D2 + 2 * sigma_mat[1,1] * sigma_mat[3,1] * E1
    F6 <- - 2 * (sigma_mat[1,1] ^ 2) * E2 * E3 + 2 * sigma_mat[1,1] * sigma_mat[2,1] * (E3 * (1 + D2) + E2 * D2) - 2 * sigma_mat[1,1] * sigma_mat[3,1] * E2 + 2 * (sigma_mat[2,1] * sigma_mat[3,1] + sigma_mat[2,2] * sigma_mat[3,2]) * (1 + D2) - 2 * ((sigma_mat[2,1] ^ 2) + (sigma_mat[2,2] ^ 2)) * D2 * (1 + D2)
    
    value <- - 0.5 * ( (F1 / (delta_mat[1,1]^3)) * ( 0.5 * (1 - exp( - 2 * T * delta_mat[1,1])) - 2 * (1 - exp( - T * delta_mat[1,1])) + delta_mat[1,1] * T) + 
                           (F2 / (delta_mat[2,2]^3)) * ( 0.5 * (1 - exp( - 2 * T * delta_mat[2,2])) - 2 * (1 - exp( - T * delta_mat[2,2])) + delta_mat[2,2] * T) + 
                           (F3 / (delta_mat[3,3]^3)) * ( 0.5 * (1 - exp( - 2 * T * delta_mat[3,3])) - 2 * (1 - exp( - T * delta_mat[3,3])) + delta_mat[3,3] * T) + 
                           (F4 / (delta_mat[1,1] * delta_mat[2,2])) * (T - (1 - exp(- T * delta_mat[1,1])) / delta_mat[1,1] - (1 - exp(- T * delta_mat[2,2])) / delta_mat[2,2] + (1 - exp(- T * (delta_mat[1,1] + delta_mat[2,2]))) / (delta_mat[1,1] + delta_mat[2,2])) + 
                           (F5 / (delta_mat[1,1] * delta_mat[3,3])) * (T - (1 - exp(- T * delta_mat[1,1])) / delta_mat[1,1] - (1 - exp(- T * delta_mat[3,3])) / delta_mat[3,3] + (1 - exp(- T * (delta_mat[1,1] + delta_mat[3,3]))) / (delta_mat[1,1] + delta_mat[3,3])) + 
                           (F6 / (delta_mat[2,2] * delta_mat[3,3])) * (T - (1 - exp(- T * delta_mat[2,2])) / delta_mat[2,2] - (1 - exp(- T * delta_mat[3,3])) / delta_mat[3,3] + (1 - exp(- T * (delta_mat[2,2] + delta_mat[3,3]))) / (delta_mat[2,2] + delta_mat[3,3])) ) / T
    
    return(value)
}
### - 3 factors
#### - B(t,T) function for dependent model (already divided by T-t, and coded as if t=0)
B_BSd_3F <- function(T,delta_mat){   # - function of a matrix delta which is lower diagonal
    D1 <- delta_mat[2,1] / (delta_mat[1,1] - delta_mat[2,2])
    D2 <- delta_mat[3,2] / (delta_mat[2,2] - delta_mat[3,3])
    D3 <- delta_mat[2,1] / (delta_mat[1,1] - delta_mat[3,3])
    D4 <- delta_mat[3,1] / (delta_mat[1,1] - delta_mat[3,3])
    D5 <- delta_mat[3,2] / (delta_mat[1,1] - delta_mat[3,3])
    E1 <- 1 + D1 + D1 * D5 + D4
    E2 <- D1 * (1 + D2)
    E3 <- D2 * D3 - D4
    
    B1 <- - E1 * (1 - exp(- delta_mat[1,1] * T)) / (delta_mat[1,1] * T) + E2 * (1 - exp(- delta_mat[2,2] * T)) / (delta_mat[2,2] * T) - E3 * (1 - exp(- delta_mat[3,3] * T)) / (delta_mat[3,3] * T)
    B2 <- - (1 + D2) * (1 - exp(- delta_mat[2,2] * T)) / (delta_mat[2,2] * T) + D2 * (1 - exp(- delta_mat[3,3] * T)) / (delta_mat[3,3] * T)
    B3 <- - (1 - exp(- delta_mat[3,3] * T)) / (delta_mat[3,3] * T)
    
    return(-c(B1, B2, B3))
}

#### -  A(t,T) function for dependent model (already divided by T-t, and coded as if t=0)
generate_BSd_3F_proj <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, total_proj_years,number_gene){
    
    n_factors <- length(kappa)
    n_ages <- nrow(mu_bar)
    n_years <- ncol(mu_bar)
    
    A_tT <- matrix(0, n_ages, 1)
    B_tT <- matrix(NA, n_ages, n_factors)
    X_t_last <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]
    #print(X_t_last)
    delta_matrix <- low_trg_fill(delta)
    
    #E_X_t1 <- exp(-kappa * proj_years) * X_t_last
     
        
    # - Build diffusion process
    ## - Build lower cholesky factor
    dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
    odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
    Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
    diag(Low_chol) <- exp(dg_l_Sigma_chol)
    #my correction
    R<-(diag(n_factors)-diag(exp(-kappa),n_factors))%*%Low_chol%*%t(Low_chol)%*%t(diag(n_factors)-diag(exp(-kappa),n_factors))
    library(MASS)
    eta<-mvrnorm(n=number_gene,mu=rep(0,n_factors),Sigma=R)
    
    
    omega_2 <- c()
    for(age in 1:n_ages){
        om_2 <- r[3]+r[1]*sum(exp(r[2]*(1:age)))/age
        omega_2 <- append(omega_2,om_2)
        
    }
    omega_2_mat <- diag(omega_2)
    epsilon <- mvrnorm(n=number_gene,mu=rep(0,n_ages),Sigma=omega_2_mat)
    
    # list for iteration
    Itr_list <- list()
    
    for (iter in 1:number_gene){
        S_prj <- matrix(NA, n_ages, total_proj_years)
        #mu_hat_prj <- matrix(NA, n_ages, number_gene)
        
        
        for (proj_years in 1:total_proj_years){
            for(age in 1:n_ages){    # - scroll over the ages
                A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)
                B_tT[age,] <- B_BSd_3F(age, delta_matrix) 
                
                X_t1 <- exp(-kappa * proj_years) * X_t_last + eta[iter,]
                
                
                
                S_prj[age, proj_years] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% X_t1 + epsilon[iter,age])) 
                
            }
            
        }
        S_young<-matrix(rep(Sx_young,total_proj_years),nrow=40)
        Itr_list[[iter]] <- rbind(S_young,S_prj)
        
    }
    
    
    
    return(Itr_list)
}


fm_parameters <- fm_bsd_3f$fit$par_est
m_parameters <- fm_bsd_3f$fit$par_est
t_parameters <- t_bsd_3f$fit$par_est

Generated_t <-generate_BSd_3F_proj(x0 = t_parameters$x0, 
                                    delta = t_parameters$delta, kappa = t_parameters$kappa, 
                                    sigma_dg = t_parameters$sigma_dg, Sigma_cov = t_parameters$Sigma_cov, 
                                    r = c(t_parameters$r1, t_parameters$r2, t_parameters$rc), 
                                    mu_bar_tot, total_proj_years = 100,number_gene = 200)

# 70 year fixed
for (i in 1:200){
    Generated_t_30[[i]]<-cbind(Generated_t_30[[i]],matrix(rep(Generated_t_30[[i]][,30],70),ncol=70))
    
}

Generated_fm <-generate_BSd_3F_proj(x0 = fm_parameters$x0, 
                     delta = fm_parameters$delta, kappa = fm_parameters$kappa, 
                     sigma_dg = fm_parameters$sigma_dg, Sigma_cov = fm_parameters$Sigma_cov, 
                     r = c(fm_parameters$r1, fm_parameters$r2, fm_parameters$rc), 
                     mu_bar_fem, total_proj_years = 50,number_gene = 100)
Generated_m <-generate_BSd_3F_proj(x0 = m_parameters$x0, 
                                    delta = m_parameters$delta, kappa = m_parameters$kappa, 
                                    sigma_dg = m_parameters$sigma_dg, Sigma_cov = m_parameters$Sigma_cov, 
                                    r = c(m_parameters$r1, m_parameters$r2, m_parameters$rc), 
                                    mu_bar_man, total_proj_years = 100,number_gene = 100)

plot(rowMeans(Generated_t[[1]]))

# life expectancy increase rate
lifeexp_mat<-matrix(nrow=100,ncol=(100-1))
for (i in 1:100){
    Itr1_lifeexp<- sapply(1:100, function(j) sum(Generated_fm[[i]][,j]))
    lifeexp_mat[i,]<-diff(Itr1_lifeexp)/head(Itr1_lifeexp,-1)
}














