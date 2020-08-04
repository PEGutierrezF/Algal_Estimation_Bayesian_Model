
    data{
    int <lower = 1> N; // number of data points
    vector [N] d13C_P; // Periphyton isotopic signal
    vector [N] d15N_P; // Periphyton isotopic signal 
    
    vector [N] d13C_T; // Terrestrial isotopic signal
    vector [N] d15N_T; // Terrestrial isotopic signal
    
    int  <lower=0> Date [N]; // Date
    int  <lower=0> Date_no;  // number of date
    }
    
    parameters { //parametros calculados y a ser calculados
    
    vector [Date_no] C;             // Date para Carbon
    vector [Date_no] Ni;             // Date para Nitrogen
    real d13C_A;                    // Carbon Algae
    real d15N_A;                    // Nitrogen Algae
    
    real <lower=0, upper=1> F_T;    // Fracion F:Terrestrial

    real beta;                      //
    
    }
   
    model{
    //priors
    d13C_A ~ normal(-23.72, 4.10); // Mean and sd of d13C of fillamentous algae, cyanobacteria and isolated diatoms on tropical streams.
    d15N_A ~ normal (4.10, 3.13); // Mean and sd of d13C of fillamentous algae, cyanobacteria and isolated diatoms on tropical streams.
    
    beta ~normal (0,1); // mean=0, Sd=1

    //likelihood
for(i in 1:N){
    d13C_P[i] ~ normal ((d13C_A + beta*C[Date[i]]) * (1- (F_T- C[Date[i]])) + d13C_T * (F_T- C[Date[i]]));
    d15N_P[i] ~ normal ((d15N_A + beta*Ni[Date[i]]) * (1- (F_T- Ni[Date[i]])) + d15N_T *(F_T- Ni[Date[i]]));
}
    }
    
    
