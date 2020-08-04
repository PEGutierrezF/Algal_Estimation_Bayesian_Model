
    data {
    int <lower=1> N; //number of data points
    vector [N] vf;
    vector [N] no3;
    int  <lower=0> region [N];
    int  <lower=0> region_no;
    }
    
    parameters {
    vector [region_no] a;
    vector [region_no] b;
    real <lower=0> sigma;
    real mub;
    real mua;    
    real <lower=0> sdb;
    real <lower=0> sda;
    }
    
    model { 
    //priors. 
    mub~ normal(0,1);
    mua~ normal(0,1);
    sdb~ normal(0,1);
    sda~ normal(0,1);
    a ~ normal(mua, sda); //prior for intercepts is weighted by the mean intercept
    b~ normal(mub,sdb);   // prior for slope is weighted by the mean slope 

    //likelihood    	
    for(i in 1:N){
    vf[i] ~ normal(a[region[i]] + b[region[i]]* no3[i] , sigma);
    }
    }
    
