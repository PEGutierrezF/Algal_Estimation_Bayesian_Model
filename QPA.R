

options(scipen=999) ## Cancel the scientific notation

#Loading required packages
library(ggplot2)
#library (ggedit)
library(plyr)
library(rstan)
rstan_options(auto_write = TRUE)##To avoid recompilation of unchanged Stan programs
options(mc.cores = parallel::detectCores()) ## Le dice que use los nucleos disponibles para el analisis
Sys.setenv(LOCAL_CPPFLAGS = '-march=native') ## Rstan recommend this for improved execution time  

sources <- read.csv("sourcesQP.csv")
sources

sink("QP_PR.stan")
cat("
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
    real d13C_A;                    // Carbon Algae, to estimate
    real d15N_A;                    // Nitrogen Algae, to estimate
    
    real <lower=0> sigma_C;         // Distribution
    real <lower=0> sigma_N;         // Distribution
    
    real <lower=0, upper=1> F_T;    // Contribution of the terrestrial fraction

    real beta;                      // ?
    
    }
   
    model{
    //priors
    d13C_A ~ normal(-23.72, 4.10); // Mean and sd of d13C of fillamentous algae, cyanobacteria and isolated diatoms on tropical streams.
    d15N_A ~ normal (4.10, 3.13); // Mean and sd of d13C of fillamentous algae, cyanobacteria and isolated diatoms on tropical streams.
    


    //likelihood
for(i in 1:N){
    d13C_P[i] ~ normal ((d13C_A + C[Date[i]]) * (1- (F_T- C[Date[i]])) + 
    d13C_T * (F_T- C[Date[i]]), sigma_C);
    
    d15N_P[i] ~ normal ((d15N_A + Ni[Date[i]]) * (1- (F_T- Ni[Date[i]])) +
    d15N_T *(F_T- Ni[Date[i]]), sigma_N);
}
    }
    
    ",
    
    fill=TRUE)
sink()


# Model -------------------------------------------------------------------
# delta13C_P =  Delta Periphyton (Biofilm), same for Nitrogen
# delta13C_T = Delta Terrestrial, same for Nitrogen

QP <- list(d13C_P=(sources$delta13C_P), d15N_P=(sources$delta15N_P), 
                       d13C_T = (sources$delta13C_T), d15N_T = (sources$delta15N_T), 
                       N = length(sources$delta13C_P), Date =(sources$Date_no), Date_no = 4)


QPCA <- rstan::stan(file = "QP_PR.stan", data = QP,
                    # control= list(adapt_delta = 0.99),
                    chains = 4, iter = 5000) # warmup= 30000,

QPCA


# Plot results ------------------------------------------------------------

traceplot(QPCA)
print(QPCA)


# Sources -----------------------------------------------------------------

mod_sources_corr <- extract(QPCA)
sources_cor <- mod_sources_corr
head(sources_cor) 

sources_QPA <- data.frame (
    
    February.2017 = c (
        mean(sources_cor$d13C_A), 
        mean (sources_cor$d15N_A), 
        sd(sources_cor$d13C_A), 
        sd(sources_cor$d15N_A), 
        
        mean(sources$delta13C_P[sources$Use=="1"]), 
        mean(sources$delta15N_P[sources$Use=="1"]), 
        sd(sources$delta13C_P[sources$Use=="1"]), 
        sd(sources$delta15N_P[sources$Use=="1"]),
        
        mean(sources$delta13C_T[sources$Use=="1"]), 
        mean(sources$delta15N_T[sources$Use=="1"]), 
        sd(sources$delta13C_T[sources$Use=="1"]), 
        sd(sources$delta15N_T[sources$Use=="1"])),
    
    November2017 = c(
        mean(sources_cor$d13C_A + sources_cor$beta1), 
        mean(sources_cor$d15N_A + sources_cor$beta4),
        sd(sources_cor$d13C_A + sources_cor$beta1),
        sd(sources_cor$d15N_A + sources_cor$beta4),
        
        mean(sources$delta13C_P[sources$Use=="2"]), 
        mean(sources$delta15N_P[sources$Use=="2"]), 
        sd(sources$delta13C_P[sources$Use=="2"]), 
        sd(sources$delta15N_P[sources$Use=="2"]),
        
        mean(sources$delta13C_T[sources$Use=="2"]), 
        mean(sources$delta15N_T[sources$Use=="2"]), 
        sd(sources$delta13C_T[sources$Use=="2"]), 
        sd(sources$delta15N_T[sources$Use=="2"])),
   
    June2018= c(
        mean(sources_cor$d13C_A + sources_cor$beta2),
        mean(sources_cor$d15N_A + sources_cor$beta5),
        sd(sources_cor$d13C_A + sources_cor$beta2),
        sd(sources_cor$d15N_A + sources_cor$beta5),
        
        mean(sources$delta13C_P[sources$Use=="3"]), 
        mean(sources$delta15N_P[sources$Use=="3"]), 
        sd(sources$delta13C_P[sources$Use=="3"]), 
        sd(sources$delta15N_P[sources$Use=="3"]),
        
        mean(sources$delta13C_T[sources$Use=="3"]), 
        mean(sources$delta15N_T[sources$Use=="3"]), 
        sd(sources$delta13C_T[sources$Use=="3"]), 
        sd(sources$delta15N_T[sources$Use=="3"])),
    
    February2019= c(
        mean(sources_cor$d13C_A + sources_cor$beta3),
        mean(sources_cor$d15N_A + sources_cor$beta6),
        sd(sources_cor$d13C_A + sources_cor$beta3),
        sd(sources_cor$d15N_A + sources_cor$beta6),
        
        mean(sources$delta13C_P[sources$Use=="4"]), 
        mean(sources$delta15N_P[sources$Use=="4"]), 
        sd(sources$delta13C_P[sources$Use=="4"]), 
        sd(sources$delta15N_P[sources$Use=="4"]),
        
        mean(sources$delta13C_T[sources$Use=="4"]), 
        mean(sources$delta15N_T[sources$Use=="4"]), 
        sd(sources$delta13C_T[sources$Use=="4"]), 
        sd(sources$delta15N_T[sources$Use=="4"]))
)

rownames(sources_QPA) <- c("d13C_A","d15N_A", "sd_d13C_A", "sd_d15N_A", "d13C_P",
                         "d15N_P", "sd_d13C_P", "sd_d15N_P", "d13C_T", "d15N_T",
                         "sd_d13C_T", "sd_d15N_T")
print(sources_QPA, digits=2)


