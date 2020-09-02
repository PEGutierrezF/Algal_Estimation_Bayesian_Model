



#--------------------------------------------
# Algae contribution from Biofilm -  Terrestrial sources
# 15 Aug 2020
#PEGF
#--------------------------------------------
#

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
                    control= list(adapt_delta = 0.99, max_treedepth=12),
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
    
    February.2017 = c(
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


###########################################################################
# Model to TEF and Food webs ----------------------------------------------
###########################################################################

# src= Sources

sink("QPA_FW.stan")
cat("

data{
    int <lower=1> N; // Number of individuals
    vector [N] d13C_ind; // 13Carbon signature for each Individual
    vector [N] d15N_ind; // 15Nitrogen signature for each Individual
    
    int <lower=0> src_no; // Number of posible basal resources
    vector [src_no] src_C; // Mean value of 13 Carbon in each basal resources
    vector [src_no] SD_src_C; // Standar deviation value of 13 Carbon in each basal resources
    vector [src_no] src_N; // Mean value of 15 Nitrogen in each basal resources
    vector [src_no] SD_src_N; // Standar deviation value value of 15 Nitrogen in each basal resources
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
    vector [src_no] ilr.global;
    
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
    
    // draw p.global (global proportion means) from an uninformative Dirichlet,
    // then ilr.global is the ILR-transform of p.global. Egozcue 2003 (pages 296)
    // Pavel -> Theta.  
    for(k in 1:(src_no - 1)){
    Theta_mean <- rows_dot_self(Theta[k]) ^ (1/k);
    ilr.global[k] <- sqrt(k/(k+1)) * log(Theta_mean[k]/Theta[k+1]);
    }
    
    // Dont generate individual deviates from the global mean, but keep same model structure
    for(i in 1:N){
    for(src in 1:(n.sources - 1)){
        ilr.ind[i,src] <- 0;
        ilr.total[i,src] <- ilr.global[src] + ilr.ind[i,src]; // add all effects togeter for each individual
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
    L ~ uniform (0,5)               // Trophic level ranges from 0 to 5 
    
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
        
"

    ,fill=TRUE)
sink()

QPA_Feb <- read.csv("Stan_QPA/QPA_Feb.csv")
QPA_Feb 

sources_QPA_Feb <- read.csv("Stan_QPA/sources_QPA_Feb.csv")
sources_QPA_Feb 


QPAlist <- list(d13C_ind=(QPA_Feb$d13C_ind), d15N_ind=(QPA_Feb$d15N_ind), 
                N= length(QPA_Feb$d13C_ind), src_no=3,
                src_C=(sources_QPA_Feb$meand13CPl), SD_src_C=(sources_QPA_Feb$SDd13C), 
                src_N = (sources_QPA_Feb$meand15NPl), SD_src_N=(sources_QPA_Feb$SDd15N))

QPA_FW <- stan(file='QPA_FW.stan', data= QPAlist,
               warmup=98000,
               chains=4, iter= 1000)    
    
    
    
    

