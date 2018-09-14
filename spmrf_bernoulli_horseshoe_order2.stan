
data {
    int<lower=0> N;              // number of observations
    int<lower=0> J;              // number of grid cells
    vector [N] xvar1;            // locations for observations
    vector [J-1] duxvar1;        // distances between unique locations
    int<lower=0> xrank1[N];      // rank order of location for each obs
    int <lower=0,upper=1> y[N];  // response for obs j
}
transformed data {
  	vector <lower=0> [J-2] drat;
    vector <lower=0> [J-2] sdrat;
    real muy;
    real sdy;
    real <lower=0,upper=1> pp[N];
	real logp[N];

	for (i in 1:N){
		pp[i] = y[i]+0.0;
		if (pp[i]==0.0)
			pp[i] = 0.005;
		if (pp[i]==1.0)
			pp[i] = 0.995;
		logp[i] = logit(pp[i]);
	}

    muy = logit(mean(pp));
    sdy = sd(logp);
    for (k in 1:(J-2)){
	    drat[k] = duxvar1[k+1]/duxvar1[k];
	    sdrat[k] = sqrt(0.5*square(duxvar1[k+1])*(duxvar1[k] + duxvar1[k+1]));
    }
}
parameters {
    real zdelta[J-1];
    real ztheta1;
    real <lower=0, upper=1> ztau[J-2];
    real <lower=0, upper=1> zgam;
    real <lower=0, upper=1> zptau2;
}
transformed parameters{
    vector[J] theta;
    real <lower=0> gam;
    real <lower=0> ptau2;
    vector[J-2] tau;

    gam = 0.01*tan(zgam*pi()/2);
    ptau2 = gam*(1/sqrt(3.0))*tan(zptau2*pi()/2);
    theta[1] = 5*sdy*ztheta1 + muy;
    theta[2] = ptau2*zdelta[1] + theta[1];
    for (j in 1:(J-2)){
        tau[j] = gam*tan(ztau[j]*pi()/2);
        theta[j+2] = zdelta[j+1]*tau[j]*sdrat[j] + (1+drat[j])*theta[j+1]-drat[j]*theta[j];
    }
}
model
{
    zgam ~ uniform(0, 1);
    zptau2 ~ uniform(0, 1);
    ztau ~ uniform(0, 1);
    ztheta1 ~ normal(0, 1);
    zdelta ~ normal(0, 1);

    for (i in 1:N)
    {
        y[i] ~ bernoulli_logit(theta[xrank1[i]]);
    }
}

generated quantities
{
  vector[N] log_lik;

  for( i in 1:N )
  {
    log_lik[i] = bernoulli_logit_lpmf(y[i] | theta[xrank1[i]]);
  }
}
