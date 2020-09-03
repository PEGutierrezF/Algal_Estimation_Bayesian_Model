

data{
    int <lower=1> N; // Number of individuals
    vector [N] d13C_ind; // 13Carbon signature for each Individual
    vector [N] d15N_ind; // 15Nitrogen signature for each Individual
    
    int <lower=0> src_no; // Number of posible basal resources
    vector [src_no] src_C; // Mean value of 13 Carbon in each basal resources
    vector [src_no] SD_src_C; // Standard deviation of 13 Carbon in each basal resources
    vector [src_no] src_N; // Mean value of 15 Nitrogen in each basal resources
    vector [src_no] SD_src_N; // Standard deviation of 15 Nitrogen in each basal resources
    }
    
parameters{
    
    real <lower=0> DeltaC;
    real <lower=0> DeltaN;
    real <lower=0> L;   // Trophic Level
    
    ordered [src_no] src_C_mean; // An ordered vector type in Stan represents a vector whose entries are sorted in ascending order
    ordered [src_no] src_N_mean;
    simplex [src_no] Theta; // Is a vector with non-negative values whose entries sum to 1
    
    // residual error (e_) representing additional unquantified
    // variation between individual; normally distributed with
    // mean = 0 and standard deviation= s. (Jackson et al. 2009)
    vector <lower=0> [src_no] e_C;  
    vector <lower=0> [src_no] e_N;
    
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
    mu_C[k] = src_C_mean[k] + DeltaC * L;
    }
    
    // 15N Stable Isotopes
    for (k in 1:src_no){
    mu_N[k] = src_N_mean[k] + DeltaN * L;
    }
    
    // to integrate proportion on sigma estimation
    for(k in 1:src_no){
    Theta2[k] = Theta[k]^2;
    }
    
    // draw p.global (global proportion means) from an uninformative Dirichlet,
    // then ilr_global is the ILR-transform of p.global. Egozcue 2003 (pages 296)
    // Pavel: gmean = Theta, k = src.  
    for(k in 1:(src_no - 1)){
    Theta_mean = rows_dot_self(Theta[k])^(1/k);
    ilr_global[k] = sqrt(k/(k+1)) * log(Theta_mean[k] / Theta[k+1]); // page 296, Egozcue 2003
    }
    
    // DON'T generate individual deviates from the global/region/pack mean (but keep same model structure)
    for(i in 1:N){
    for(src in 1:(n.sources-1)) {
        ilr_ind[i,src] = 0;
        ilr_total[i,src] = ilr_global[src] + ilr_ind[i,src]; // add all effects together for each individual (in ilr-space)
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
    
    DeltaC ~ normal(-0.41, 1.14);   // (Vander Zanden and Rasmussen (2001)
    DeltaN ~ normal(0.6, 1.7);      // Bunn, Leigh, & Jardine, (2013)
    L ~ uniform (0,5);               // Trophic level ranges from 0 to 5 
    
// likelihood

    for(k in 1:src_no){
    sigma_C[k] ~ normal(SD_src_C[K] * Theta2[k], e_C);
    }
    
        for(k in 1:src_no){
    sigma_N[k] ~ normal(SD_src_N[K] * Theta2[k], e_N);
        }
    
    for(i in 1:N){
    for(k in 1:src_no){
    contributions[k] = log(Theta[k]) + normal_lpdf (d13C_ind | mu_C, sigma_C[k]);
    }
        }
    
        for(i in 1:N){
    for(k in 1:src_no){
    contributions[k] = log(Theta[k]) + normal_lpdf (d15N_ind | mu_N, sigma_N[k]);
    }
        }
    
    // The log density increment statement (target +=) is used to add log_sum_exp(contributions)
    // to the log density defined by the rest of the program    
    // log density Jacobian adjustment
    
    target += log_sum_exp(contributions);  
        
    }
        

