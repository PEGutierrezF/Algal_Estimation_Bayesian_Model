

data{
    int <lower=1> N; // Number of individuals
    vector [N] d13C_ind; // 13Carbon signature for each Individual
    vector [N] d15N_ind; // 15Nitrogen signature for each Individual
    
    int src_no; // Number of posible basal resources
    vector [src_no] src_C; // Mean value of 13 Carbon in each basal resources
    vector [src_no] SD_src_C; // Standard deviation of 13 Carbon in each basal resources
    vector [src_no] src_N; // Mean value of 15 Nitrogen in each basal resources
    vector [src_no] SD_src_N; // Standard deviation of 15 Nitrogen in each basal resources
    }
    
parameters{
    
    real Delta_C;
    real Delta_N;
    real <lower=0> L;   // Trophic Level
    
    ordered [src_no] src_C_mean; // An ordered vector type in Stan represents a vector whose entries are sorted in ascending order
    ordered [src_no] src_N_mean;
    simplex [src_no] Theta; // Is a vector with non-negative values whose entries sum to 1

    vector <lower = 0>[src_no] sigma_C;
    vector <lower = 0>[src_no] sigma_N;
    }

model{
    vector[src_no] contributions;
  // priors
  for (k in 1:src_no){
  src_C_mean[k] ~ normal(src_C[k], SD_src_C[k]); //values from field Data samples of sources
  }
  
  for (k in 1:src_no){
  src_N_mean[k] ~ normal(src_N[k], SD_src_N[k]); //values from field Data samples of sources
  }
  
    Theta ~ dirichlet(rep_vector(1, src_no)); // Uniform prior
    
    Delta_C ~ normal(0.41, 1.14);   // Vander Zanden and Rasmussen (2001)
    Delta_N ~ normal(0.6, 1.7);      // Bunn, Leigh, & Jardine, (2013)
    L ~ uniform (0,10);               // Trophic level ranges from 0 to 10 
    
    sigma_C ~ normal(0,10);
    sigma_N ~ normal(0,10);
  
// likelihood

  for(i in 1:N) {
    for(k in 1:src_no) {
      contributions[k] = log(Theta[k]) + normal_lpdf(d13C_ind[i] | (src_C_mean[k]+(Delta_C*L)), sigma_C[k]);
    }
    target += log_sum_exp(contributions);
  }
    
   for(i in 1:N) {
    for(k in 1:src_no) {
      contributions[k] = log(Theta[k]) + normal_lpdf(d15N_ind[i] | (src_N_mean[k]+(Delta_N*L)), sigma_N[k]);
    }
    target += log_sum_exp(contributions);
   }
    
}
   
    
