# based on EWJ and also Xuejia Ke, who made a demo R package: https://github.com/thebestecho/demo-RNAFracQuant
# I rewrote Edward's functions to be flexible both for the number of columns to mix and the names of the input columns
# see below for example usage
make_stan_data <- function(ctdata,column_list) {
  # make the data into the right format for stan
  stan_dat <- list(NRNA=nrow(ctdata),
                   tot_obs=as.integer(round(ctdata$Total)))
  for (col in column_list) {
    t <- structure(list(as.integer(round(ctdata[[col]]))),
                   names=list(paste0(col,"_obs")))
    stan_dat <- append(stan_dat,t)
  }
  return(stan_dat)
}

make_stan_model <- function(column_list) {
  model <- paste0(
    '// -*- mode: C -*-
  data {
    // Number of RNAs
    int<lower=1> NRNA;     
    
    // Note: These are all integers
    // columns t, s, p
    int<lower=0> tot_obs[NRNA];\n',
    paste(paste0("\tint<lower=0> ",column_list,"_obs[NRNA];\n"),collapse=""),
    '}
  parameters {\n',
    paste(paste0("\treal<lower=0> mixing_",column_list,";\n"),collapse=""),
    '  
    // dispersion parameter for counts
    real<lower=0> phi;
  }
  model{
    // mixing ratios\n',
    paste(paste0("\tmixing_",column_list," ~gamma(1,1);\n"),collapse=""),
    '
    // Cauchy prior for negbin dispersion parameter
    phi ~ cauchy(0,3);
    
    for(idx in 1:NRNA){ 
      // count distn negative binomial with specified means
      // Total
      tot_obs[idx] ~ neg_binomial_2(\n',
    paste("\t",paste0("mixing_",column_list[1:length(column_list)-1]," * ",
                      column_list[1:length(column_list)-1],"_obs[idx] + "),collapse=""),
    paste0("mixing_",column_list[length(column_list)]," * ",
           column_list[length(column_list)],"_obs[idx], phi);\n"),
    ' }
  
  }')
  return(model)
}

initialize_stan_mixing_model <- function(ctdata,column_list,model=NULL) {
  # first generate the matrix for the specific columns to mix
  stan_dat <- make_stan_data(ctdata,column_list)
  
  # now generate the text for the stan model from the columns
  if (is.null(model)) {model <- make_stan_model(column_list)}
  
  #initialize the stan model
  stan(model_code=model,data=stan_dat,chains = 1,iter = 100)
}

fit_mixing_model <- function(ctdata,column_list,n_iter=1000,n_chains=4,
                             control=list(adapt_delta=0.85),
                             stan_mixing=NULL,...) {
  stan_dat <- make_stan_data(ctdata,column_list)
  if (is.null(stan_mixing)) {
    stan_mixing <- initialize_stan_mixing_model(ctdata,column_list)
  }
  stan_mixing_fit <- stan(fit = stan_mixing,data = stan_dat,
                          iter = n_iter,chains = n_chains,...)
  return(stan_mixing_fit %>% summary())
}

get_mixing_params <- function(item){
  # extracts all the parameters from the stan model
  out <- tibble(Variable=rownames(item$summary),
                mean=item$summary[,"mean"],
                se_mean=item$summary[,"se_mean"],
                sd=item$summary[,"sd"],
                "x2.5" = item$summary[,"2.5%"],
                x25 = item$summary[,"25%"],
                x50 = item$summary[,"50%"],
                x75 = item$summary[,"75%"],
                "x97.5" = item$summary[,"97.5%"],
                n_eff = item$summary[,"n_eff"],
                Rhat = item$summary[,"Rhat"]) %>%
    pivot_longer(cols=c(mean:Rhat),names_to="Term",values_to="Value")
}
