
    data{
    int <lower = 1> N; // number of data points
    vector [N] d13C_P; // Periphyton isotopic signal
    vector [N] d15N_P; // Periphyton isotopic signal 
    
    vector [N] d13C_T; // Terrestrial isotopic signal
    vector [N] d15N_T; // Terrestrial isotopic signal
    
    int  <lower=0> Stream_N [N]; // November 2017
    int  <lower=0> Stream_N_no;  // number of stream type
    
    int  <lower=0> Stream_J [N]; // June 2018
    int  <lower=0> Stream_J_no;  // number of stream type
    
    int  <lower=0> Stream_F [N]; // February 2019
    int  <lower=0> Stream_F_no;  // number of stream type
    
    }
    
    parameters {
    
    real d13C_A;
    real d15N_A;
    real <lower=0, upper=1> F_T;
    real <lower=0, upper=1> alpha1;
    real <lower=0, upper=1> alpha2;
    real <lower=0, upper=1> alpha3;
    real <lower=0> sigma_C;
    real <lower=0> sigma_N;
    
    real beta1;
    real beta2;
    
    real beta3;
    real beta4;
    
    real beta5;
    real beta6;
    
    real beta7;
    real beta8;
    
    }
   
    model{
    //priors
    d13C_A ~ normal(-23.72, 4.10); //Mean and sd of d13C of fillamentous algae, cyanobacteria and isolated diatoms on tropical streams.
    
    d15N_A ~ normal (4.10, 3.13); //Mean and sd of d13C of fillamentous algae, cyanobacteria and isolated diatoms on tropical streams.
    
    beta1 ~normal (0,1);
    beta2 ~normal (0,1);
    
    beta3 ~normal (0,1);
    beta4 ~normal (0,1);
    
    beta5 ~normal (0,1);
    beta6 ~normal (0,1);
    
    beta7 ~normal (0,1);
    beta8 ~normal (0,1);
    
    //likelihood
    
    d13C_P ~ normal ((d13C_A + beta1*Stream_N_no + beta2*Stream_J_no + beta3*Stream_F_no)*(1-(F_T- alpha1*Stream_N_no - alpha2*Stream_J_no - alpha3*Stream_F_no)) + d13C_T* (F_T- alpha1*Stream_N_no - alpha2*Stream_J_no - alpha3*Stream_F_no), sigma_C);
    d15N_P ~ normal ((d15N_A + beta4*Stream_N_no + beta5*Stream_J_no + beta6*Stream_F_no )*(1- (F_T- alpha1*Stream_N_no - alpha2*Stream_J_no - alpha3*Stream_F_no)) + d15N_T*(F_T- alpha1*Stream_N_no - alpha2*Stream_J_no - alpha3*Stream_F_no), sigma_N);
    
    }
    
    
