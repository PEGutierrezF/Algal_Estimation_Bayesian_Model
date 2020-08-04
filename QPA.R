

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
                    chains = 4, iter = 30000) # warmup= 30000,

QPCA


# Plot results ------------------------------------------------------------

traceplot(QPCA)
print(QPCA)


# Sources -----------------------------------------------------------------

mod_sources_corr <- extract(QPCA)
sources_cor <- mod_sources_corr
head(sources_cor) 

sources_QPA <- data.frame (
    
    February.2017 = c(
            mean(sources_cor$d13C_A),
            mean(sources_cor$d15N_A),
            sd(sources_cor$d13C_A),
            sd(sources_cor$d15N_A),
            
            mean(sources_cor$d13C_A),
            mean(sources_cor$d15N_A),
            sd(sources_cor$d13C_A),
            sd(sources_cor$d15N_A),
            
            mean(sources$delta13C_P[sources$Date_no=="1"]), 
            mean(sources$delta15N_P[sources$Date_no=="1"]), 
            sd(sources$delta13C_P[sources$Date_no=="1"]), 
            sd(sources$delta15N_P[sources$Date_no=="1"]),
            
            mean(sources$delta13C_T[sources$Date_no=="1"]), 
            mean(sources$delta15N_T[sources$Date_no=="1"]), 
            sd(sources$delta13C_T[sources$Date_no=="1"]), 
            sd(sources$delta15N_T[sources$Date_no=="1"])),
    
    November.2017 = c(
        mean(sources_cor$d13C_A), 
        mean(sources_cor$d15N_A),
        sd(sources_cor$d13C_A),
        sd(sources_cor$d15N_A),
        
        mean(sources$delta13C_P[sources$Date_no=="2"]), 
        mean(sources$delta15N_P[sources$Date_no=="2"]), 
        sd(sources$delta13C_P[sources$Date_no=="2"]), 
        sd(sources$delta15N_P[sources$Date_no=="2"]),
        
        mean(sources$delta13C_T[sources$Date_no=="2"]), 
        mean(sources$delta15N_T[sources$Date_no=="2"]), 
        sd(sources$delta13C_T[sources$Date_no=="2"]), 
        sd(sources$delta15N_T[sources$Date_no=="2"])),
   
    June.2018= c(
        mean(sources_cor$d13C_A),
        mean(sources_cor$d15N_A),
        sd(sources_cor$d13C_A),
        sd(sources_cor$d15N_A),
        
        mean(sources$delta13C_P[sources$Date_no=="3"]), 
        mean(sources$delta15N_P[sources$Date_no=="3"]), 
        sd(sources$delta13C_P[sources$Date_no=="3"]), 
        sd(sources$delta15N_P[sources$Date_no=="3"]),
        
        mean(sources$delta13C_T[sources$Date_no=="3"]), 
        mean(sources$delta15N_T[sources$Date_no=="3"]), 
        sd(sources$delta13C_T[sources$Date_no=="3"]), 
        sd(sources$delta15N_T[sources$Date_no=="3"])),
    
    February.2019= c(
        mean(sources_cor$d13C_A),
        mean(sources_cor$d15N_A),
        sd(sources_cor$d13C_A),
        sd(sources_cor$d15N_A),
        
        mean(sources$delta13C_P[sources$Date_no=="4"]), 
        mean(sources$delta15N_P[sources$Date_no=="4"]), 
        sd(sources$delta13C_P[sources$Date_no=="4"]), 
        sd(sources$delta15N_P[sources$Date_no=="4"]),
        
        mean(sources$delta13C_T[sources$Date_no=="4"]), 
        mean(sources$delta15N_T[sources$Date_no=="4"]), 
        sd(sources$delta13C_T[sources$Date_no=="4"]), 
        sd(sources$delta15N_T[sources$Date_no=="4"]))
)

rownames(sources_QPA) <- c("d13C_A","d15N_A", "sd_d13C_A", "sd_d15N_A", "d13C_P",
                         "d15N_P", "sd_d13C_P", "sd_d15N_P", "d13C_T", "d15N_T",
                         "sd_d13C_T", "sd_d15N_T")

print(sources_QPA, digits=2)


