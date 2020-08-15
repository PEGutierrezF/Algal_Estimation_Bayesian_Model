
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
   
    vector [Date_no] d13C_A;             // vector de valores a estimar d13C para algas en las distintas fechas
    vector [Date_no] d15N_A;             // vector de valores a estimar d15N para algas en las distintas fechas
    vector <lower=0, upper=1> [Date_no] F_T;    // vector de la fracción de material aloctono en el perifiton por fecha
   
    //partial pooling
    real mu_d13C_A;
    real sd_d13C_A;
   
    real mu_d15N_A;
    real sd_d15N_A;
   
    real mu_F_T;
    real sd_F_T;
   
    real <lower=0> sigma_C;         // Distribution
    real <lower=0> sigma_N;         // Distribution
   
   
    }
   
    model{
    //priors
    mu_d13C_A ~ normal(-23.72, 4.10); // Mean and sd of d13C of fillamentous algae, cyanobacteria and isolated diatoms on tropical streams.
    sd_d13C_A~ normal(0, 10);
    d13C_A ~ normal(mu_d13C_A, sd_d13C_A);
   
    mu_d15N_A ~ normal (4.10, 3.13); // Mean and sd of d13C of fillamentous algae, cyanobacteria and isolated diatoms on tropical streams.
    sd_d15N_A~ normal(0, 10);
    d15N_A ~ normal(mu_d15N_A, sd_d15N_A);
   
    mu_F_T ~ normal(0.5, 0.5); // no estoy muy seguro de este valor. Estoy dandole uno bastante amplio si va de 0 a 1
    sd_F_T ~ normal (0, 0.1);
    F_T ~ normal(mu_F_T, sd_F_T);

    //likelihood
    for(i in 1:N){
    d13C_P[i] ~ normal ((d13C_A[Date[i]]) * (1- (F_T[Date[i]])) +
                            d13C_T[i] * (F_T[Date[i]]), sigma_C);
    }
   
    for (i in 1:N){
    d15N_P[i] ~ normal ((d15N_A[Date[i]]) * (1- (F_T[Date[i]])) +
                            d15N_T[i] * (F_T[Date[i]]), sigma_N);
    }
   

    }
   
    
