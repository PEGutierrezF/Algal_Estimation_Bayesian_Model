



#--------------------------------------------
# Algae contribution from Biofilm -  Terrestrial sources
# 15 Aug 2020
#PEGF
#--------------------------------------------
#

#Loading required packages
library(ggplot2)
#library (ggedit)
library(plyr)
library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)##To avoid recompilation of unchanged Stan programs
options(mc.cores = parallel::detectCores()) ## Le dice que use los nucleos disponibles para el analisis
Sys.setenv(LOCAL_CPPFLAGS = '-march=native') ## Rstan recommend this for improved execution time  

sources<- read.csv("sourcesQPB.csv")
head(sources)

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
    real <lower=0> sd_d13C_A;
   
    real mu_d15N_A;
    real <lower=0> sd_d15N_A;
   
    //real  <lower=0, upper=1> mu_F_T;
    //real <lower=0> sd_F_T;
   
    real <lower=0> sigma_C;         // Distribution
    real <lower=0> sigma_N;         // Distribution
   
   
    }
   
    model{
    //priors
    mu_d13C_A ~ normal(-22.6352941176471, 3.9344887178118); // Mean and sd of d13C of fillamentous algae, cyanobacteria and isolated diatoms on tropical streams.
    sd_d13C_A ~ normal(0, 10);
    d13C_A ~ normal(mu_d13C_A, sd_d13C_A);
   
    mu_d15N_A ~ normal (3.67588235294118, 2.46396798585008); // Mean and sd of d13C of fillamentous algae, cyanobacteria and isolated diatoms on tropical streams.
    sd_d15N_A ~ normal(0, 10);
    d15N_A ~ normal(mu_d15N_A, sd_d15N_A);
   
    //mu_F_T ~ normal(0.5, 0.5); // no estoy muy seguro de este valor. Estoy dandole uno bastante amplio si va de 0 a 1
    //sd_F_T ~ normal (0, 1);
    F_T ~ beta (1, 1); // uniform prior for every value between 0 - 1

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

datalist_QP <- list(d13C_P=(sources$delta13C_P), d15N_P=(sources$delta15N_P), 
                       d13C_T = (sources$delta13C_T), d15N_T = (sources$delta15N_T), 
                       N = length(sources$delta13C_P), Date =(sources$Use), Date_no = 4)

datalist_QP
QPCA <- stan(file = "QP_PR.stan", data = datalist_QP,
                    #control= list(adapt_delta = 0.999, max_treedepth=12),
                    warmup= 48000, chains = 4, iter = 50000)

# Plot results ------------------------------------------------------------

traceplot(QPCA)
print(QPCA)


# Sources -----------------------------------------------------------------

sources_cor <- extract(QPCA)
head(sources_cor) 

sources_QPB <- data.frame (
    February_2017 = c (
        mean(sources_cor$d13C_A[,1]),
        sd(sources_cor$d13C_A[,1]),
        mean(sources_cor$d15N_A[,1]),
        sd(sources_cor$d15N_A[,1]),
        
        mean(sources$delta13C_P[sources$Use=="1"]), 
        sd(sources$delta13C_P[sources$Use=="1"]),
        mean(sources$delta15N_P[sources$Use=="1"]), 
        sd(sources$delta15N_P[sources$Use=="1"]), 
        
        mean(sources$delta13C_T[sources$Use=="1"]), 
        sd(sources$delta13C_T[sources$Use=="1"]),
        mean(sources$delta15N_T[sources$Use=="1"]), 
        sd(sources$delta15N_T[sources$Use=="1"])
        ),
    
    November_2017 = c (
        mean(sources_cor$d13C_A[,2]),
        sd(sources_cor$d13C_A[,2]),
        mean(sources_cor$d15N_A[,2]),
        sd(sources_cor$d15N_A[,2]),
        
        mean(sources$delta13C_P[sources$Use=="2"]), 
        sd(sources$delta13C_P[sources$Use=="2"]),
        mean(sources$delta15N_P[sources$Use=="2"]), 
        sd(sources$delta15N_P[sources$Use=="2"]), 
        
        mean(sources$delta13C_T[sources$Use=="2"]), 
        sd(sources$delta13C_T[sources$Use=="2"]),
        mean(sources$delta15N_T[sources$Use=="2"]), 
        sd(sources$delta15N_T[sources$Use=="2"])
        ),
    June_2018 = c(
        mean(sources_cor$d13C_A[,3]),
        sd(sources_cor$d13C_A[,3]),
        mean(sources_cor$d15N_A[,3]),
        sd(sources_cor$d15N_A[,3]),
        
        mean(sources$delta13C_P[sources$Use=="3"]), 
        sd(sources$delta13C_P[sources$Use=="3"]),
        mean(sources$delta15N_P[sources$Use=="3"]), 
        sd(sources$delta15N_P[sources$Use=="3"]), 
        
        mean(sources$delta13C_T[sources$Use=="3"]), 
        sd(sources$delta13C_T[sources$Use=="3"]),
        mean(sources$delta15N_T[sources$Use=="3"]), 
        sd(sources$delta15N_T[sources$Use=="3"])
        ),

    February_2019= c(
        mean(sources_cor$d13C_A[,4]),
        sd(sources_cor$d13C_A[,4]),
        mean(sources_cor$d15N_A[,4]),
        sd(sources_cor$d15N_A[,4]),
    
        mean(sources$delta13C_P[sources$Use=="4"]), 
        sd(sources$delta13C_P[sources$Use=="4"]),
        mean(sources$delta15N_P[sources$Use=="4"]), 
        sd(sources$delta15N_P[sources$Use=="4"]), 
    
        mean(sources$delta13C_T[sources$Use=="4"]), 
        sd(sources$delta13C_T[sources$Use=="4"]),
        mean(sources$delta15N_T[sources$Use=="4"]), 
        sd(sources$delta15N_T[sources$Use=="4"])
     )
    )


rownames(sources_QPB) <- c("d13C_A","sd_d13C_A","d15N_A", "sd_d15N_A",
                           "d13C_P","sd_d13C_P", "d15N_P", "sd_d15N_P", 
                           "d13C_T", "sd_d13C_T", "d15N_T","sd_d15N_T")
sources_QPB <-as.data.frame(t(as.matrix(sources_QPB)))
print(sources_QPB, digits=2)
write.csv(sources_QPB, file= "Sources_QPB_results.csv") ## export as csv for further results

