

data{
    int <lower=1> N; // Number of individuals
    vector [N] d13C_ind; // 13 Carbon signature for each Individual
    vector [N] d15N_ind; // 15 Nitrogen signature for each Individual
    
    int <lower=0> src_no; // Number of posible basal resources
    vector [src_no] src_C; // Mean value of 13 Carbon in each basal resources
    vector [src_no] SD_src_C; // SD value of 13 Carbon in each basal resources
    vector [src_no] src_N; // Mean value of 15 Nitrogen in each basal resources
    vector [src_no] SD_src_N; // SD value of 15 Nitrogen in each basal resources
    }
    
parameters{
    
    real <lower=0> DeltaC;
    real <lower=0> DeltaN;
    real <lower=0> L;   // Trophic Level
    
    ordered [src_no] src_C_mean;
    ordered [src_no] src_N_mean;
    simplex [src_no] Theta;
    vector <lower=0> [src_no] e_C; // residual error Jackson et al. 2009
    vector <lower=0> [src_no] e_N; // residual error Jackson et al. 2009
    vector [src_no] sigma_C;
    vector [src_no] sigma_N;
    }
    
transformed parameters{
    
    vector [src_no] mu_C;
    vector [src_no] mu_N;
    vector [src_no] Theta2;
    vector [src_no] Theta_mean;
    vector [src_no] ilr_global;
    
    // 13C Stable Isotopes
    for (k in 1:src_no){
    mu_C[k] <- src_C_mean[k] + DeltaC * L;
    }
    
    // 15N Stable Isotopes
    for (k in 1:src_no){
    mu_N[k] <- src_N_mean[k] + DeltaN * L;
    }
    
    // to integrate proportion on sigma estimation
    for(k in 1:src_no){
    Theta2[k] <- Theta[k]^2;
    }
    
    // ilr transform of global proportions (Theta) same MIXSIAR and Egozcue 2003 (pages 296)
    for(k in 1:(src_no - 1)){
    Theta_mean <- rows_dot_self(Theta[k]) ^ (1/k);
    ilr_global[k] <- sqrt(k/(k+1)) * log(Theta_mean[k]/Theta[k+1]);
    }
    
    // Dont generate individual deviates from the global mean, but keep same model structure
    for(i in 1:N){
    for(src in 1:(n.sources - 1)){
        ilr_ind[i,src] <- 0;
        ilr_total[i,src] <- ilr_global[src] + ilr_ind[i,src]; // add all effects togeter for each individual
        }
    }
}

model{
vector[src_no] contributions;

// priors
    Theta ~ dirichlet(rep_vector(1, src_no));
    for(k in 1:src_no){
    src_C_mean[k] ~ normal (src_C, SD_src_C);
    }
    
    for(k in 1:src_no){
    src_N_mean[k] ~ normal (src_N, SD_src_N);
    }
    
    DeltaC ~ normal(0.39, 1.14);
    DeltaN ~ normal(3.4, 0.99);
    
    L ~ uniform (0,6)       // Trophic level ranges from 0 to 6
    
// likelihood

    for(k in 1:src_no){
    sigma_C[k] ~ normal(SD_src_C[K] * Theta2[k], e_C);
    }
    
        for(k in 1:src_no){
    sigma_N[k] ~ normal(SD_src_N[K] * Theta2[k], e_N);
        }
    
    for(i in 1:N){
    for(k in 1:src_no){
    contributions[k] = log(Theta[k]) + normal_lpdf (d13C_ind| mu_C, sigma_C[k]);
    }
        }
    
        for(i in 1:N){
    for(k in 1:src_no){
    contributions[k] = log(Theta[k]) + normal_lpdf (d15N_ind| mu_N, sigma_N[k]);
    }
        }
        
    target += log_sum_exp(contributions);
        
    }
        

