

options(scipen=999) ## Cancel the scientific notation

#Loading required packages
library(ggplot2)
#library (ggedit)
library(plyr)
library(rstan)
rstan_options(auto_write = TRUE)##To avoid recompilation of unchanged Stan programs
options(mc.cores = parallel::detectCores()) ## Le dice que use los nucleos disponibles para el analisis
Sys.setenv(LOCAL_CPPFLAGS = '-march=native') ## Rstan recommend this for improved execution time  

setwd("D:/LTER/Manuscript 2019 Stable Isotopes/Estimations/Pavel/")
sources <- read.csv("D:/LTER/Manuscript 2019 Stable Isotopes/Estimations/Pavel/sourcesQP.csv")
sources

sink("QP_PR.stan")
cat("
    data{
    int <lower = 1> N; // number of data points
    vector [N] d13C_P; // Periphyton isotopic signal
    vector [N] d15N_P; // Periphyton isotopic signal 
    
    vector [N] d13C_T; // Terrestrial isotopic signal
    vector [N] d15N_T; // Terrestrial isotopic signal
    
    int  <lower=0> Stream_N [N]; // Stream Coffee
    int  <lower=0> Stream_N_no;  // number of stream type
    
    int  <lower=0> Stream_J [N]; // Stream Coffee
    int  <lower=0> Stream_J_no;  // number of stream type
    
    int  <lower=0> Stream_F [N]; // Stream Coffee
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




