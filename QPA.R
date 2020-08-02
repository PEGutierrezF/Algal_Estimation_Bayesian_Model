

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
    
    ",
    
    fill=TRUE)
sink()


# Model -------------------------------------------------------------------
# delta13C_P =  Delta Periphyton (Biofilm), same for Nitrogen
# delta13C_T = Delta Terrestrial, same for Nitrogen

QP <- list(d13C_P=(sources$delta13C_P), d15N_P=(sources$delta15N_P), 
                       d13C_T = (sources$delta13C_T), d15N_T = (sources$delta15N_T), 
                       N = length(sources$delta13C_P), Stream =(sources$Use), Stream_no = 4, 
                        Stream_N = (sources$November), Stream_N_no = 1, 
                        Stream_J = (sources$June), Stream_J_no = 1,
                        Stream_F = (sources$February), Stream_F_no = 1)


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

sources_A <- data.frame (
    
    Forest = c (
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
        sd(sources$delta15N_T[sources$Use=="1"]),
        
        
        mean(sources_2$d13C_C3[sources_2$use=="1"], na.rm=T), 
        mean(sources_2$d15N_C3[sources_2$use=="1"], na.rm=T), 
        sd(sources_2$d13C_C3[sources_2$use=="1"], na.rm=T), 
        sd(sources_2$d15N_C3[sources_2$use=="1"], na.rm=T),
        
        mean(sources_2$d13C_C4[sources_2$use=="1"], na.rm=T), 
        mean(sources_2$d15N_C4[sources_2$use=="1"], na.rm=T), 
        sd(sources_2$d13C_C4[sources_2$use=="1"], na.rm=T), 
        sd(sources_2$d15N_C4[sources_2$use=="1"], na.rm=T),
        
        mean(sources_2$d13C_H[sources_2$use=="1"], na.rm=T), 
        mean(sources_2$d15N_H[sources_2$use=="1"], na.rm=T), 
        sd(sources_2$d13C_H[sources_2$use=="1"], na.rm=T), 
        sd(sources_2$d15N_H[sources_2$use=="1"], na.rm=T)),
    
    Coffee = c(
        mean(sources_cor$d13C_A + sources_cor$beta1), 
        mean(sources_cor$d15N_A + sources_cor$beta3),
        sd(sources_cor$d13C_A + sources_cor$beta1),
        sd(sources_cor$d15N_A + sources_cor$beta3),
        
        mean(sources$delta13C_P[sources$Use=="2"]), 
        mean(sources$delta15N_P[sources$Use=="2"]), 
        sd(sources$delta13C_P[sources$Use=="2"]), 
        sd(sources$delta15N_P[sources$Use=="2"]),
        
        mean(sources$delta13C_T[sources$Use=="2"]), 
        mean(sources$delta15N_T[sources$Use=="2"]), 
        sd(sources$delta13C_T[sources$Use=="2"]), 
        sd(sources$delta15N_T[sources$Use=="2"]),
        
        mean(sources_2$d13C_C3[sources_2$use=="2"], na.rm=T), 
        mean(sources_2$d15N_C3[sources_2$use=="2"], na.rm=T), 
        sd(sources_2$d13C_C3[sources_2$use=="2"], na.rm=T), 
        sd(sources_2$d15N_C3[sources_2$use=="2"], na.rm=T),
        
        mean(sources_2$d13C_C4[sources_2$use=="2"], na.rm=T), 
        mean(sources_2$d15N_C4[sources_2$use=="2"], na.rm=T), 
        sd(sources_2$d13C_C4[sources_2$use=="2"], na.rm=T), 
        sd(sources_2$d15N_C4[sources_2$use=="2"], na.rm=T),
        
        mean(sources_2$d13C_H[sources_2$use=="2"], na.rm=T), 
        mean(sources_2$d15N_H[sources_2$use=="2"], na.rm=T), 
        sd(sources_2$d13C_H[sources_2$use=="2"], na.rm=T), 
        sd(sources_2$d15N_H[sources_2$use=="2"], na.rm=T)),
    
    Pasture= c(
        mean(sources_cor$d13C_A + sources_cor$beta2),
        mean(sources_cor$d15N_A + sources_cor$beta4),
        sd(sources_cor$d13C_A + sources_cor$beta2),
        sd(sources_cor$d15N_A + sources_cor$beta4),
        
        mean(sources$delta13C_P[sources$Use=="3"]), 
        mean(sources$delta15N_P[sources$Use=="3"]), 
        sd(sources$delta13C_P[sources$Use=="3"]), 
        sd(sources$delta15N_P[sources$Use=="3"]),
        
        mean(sources$delta13C_T[sources$Use=="3"]), 
        mean(sources$delta15N_T[sources$Use=="3"]), 
        sd(sources$delta13C_T[sources$Use=="3"]), 
        sd(sources$delta15N_T[sources$Use=="3"]),
        
        mean(sources_2$d13C_C3[sources_2$use=="3"], na.rm=T), 
        mean(sources_2$d15N_C3[sources_2$use=="3"], na.rm=T), 
        sd(sources_2$d13C_C3[sources_2$use=="3"], na.rm=T), 
        sd(sources_2$d15N_C3[sources_2$use=="3"], na.rm=T),
        
        mean(sources_2$d13C_C4[sources_2$use=="3"], na.rm=T), 
        mean(sources_2$d15N_C4[sources_2$use=="3"], na.rm=T), 
        sd(sources_2$d13C_C4[sources_2$use=="3"], na.rm=T), 
        sd(sources_2$d15N_C4[sources_2$use=="3"], na.rm=T),
        
        mean(sources_2$d13C_H[sources_2$use=="3"], na.rm=T), 
        mean(sources_2$d15N_H[sources_2$use=="3"], na.rm=T), 
        sd(sources_2$d13C_H[sources_2$use=="3"], na.rm=T), 
        sd(sources_2$d15N_H[sources_2$use=="3"], na.rm=T))
    
)


