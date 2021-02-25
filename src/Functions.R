#### Defining functions ####

###  Determines the squared distance between the preferred and observed R2 and prevalence.
dist_R2_prev <- function(par, pref_R2, pref_prev){
  # par is a vector with initial guesses for both
  # the intercept and the regression coefficients.
  # However since regression have been restricted
  # only one value is necessary.
  
  # Providing the intercept and regression coefficients as specified in the par object. 
  dgm_par <-
    c(par[1], 
      rep(par[2] * 3, round(0.3 * n_pred)),  # strong
      rep(par[2], round(0.5 *  n_pred)),  # weaker
      rep(par[2] * 0, round(0.2 * n_pred))) # noise
  
  # Obtain values for y based on Bernoulli distribution, with input p
  p <- 1/(1+exp(-dm %*% dgm_par))
  y <- as.numeric(p > runif(length(p))) 
  
  # Obtain observed values of c-statistic and 
  # average predicted probability of an event
  
  obs_R2 <- pseudo_Rsqrs(p = p, y = y) # obtain R2cs based on p and y
  obs_prev <- mean(y) # "Observed" prevalence
  
  # Sum of squared differences of both values:
  (obs_R2-pref_R2)^2 + (obs_prev-pref_prev)^2
  
}

#### Checks whether the optimized parameters look like we want them to:
checking <- function(par){
  ## What do the observed prevalence and c-statistic look like?
  # Providing the  and regression coefficients as specified in the par object. 
  
  # Defining candidate predictors
  dgm_par <-
    c(par[1], 
      rep(par[2] * 3, round(0.3 * n_pred)), # Strong
      rep(par[2], round(0.5 *  n_pred)),  # Weaker
      rep(par[2] * 0, round(0.2 * n_pred))) # Noise
  
  # Obtain values for y based on uniform distribtuin, with input p
  p <- 1/(1+exp(-dm %*% dgm_par))
  y <- as.numeric(p > runif(length(p))) 
  
  # Obtain observed values
  obs_cstat <- fastAUC(p = p, y = y)
  obs_prev <- mean(y) 
  c("cstat" = obs_cstat, "prev" = obs_prev)
}

### Same idea as the checking function, however, using the validation data:
checking_val <- function(par){
  
  # Providing the intercept and regression coefficients as specified in the par object. 
  dgm_par_val <-
    c(par[1], 
      rep(par[2] * 3, round(0.3 * n_pred)), 
      rep(par[2], round(0.5 *  n_pred)), 
      rep(par[2] * 0, round(0.2 * n_pred)))
  
  
  # Obtain values for y based on uniform distribution, with input p
  p_val <- 1/(1+exp(-dm_val %*% dgm_par_val))
  y_val <- as.numeric(p_val > runif(length(p_val)))
  
  # Obtain observed values
  #obs_cstat <- c_stat2(preds = p_val, outcome = y_val) # obtain c-statistic based on p and y
  obs_cstat <- fastAUC(p = p_val, y = y_val)
  obs_prev <- mean(y_val)
  
  # return results
  c("cstat" = obs_cstat, "prev" = obs_prev)
}

## R^2 Cox-Snell
pseudo_Rsqrs <- function(p, y){ 
  
  .LL <- function(p, y){
    sum(y*log(p)+(1-y)*log(1-p))
  }
  
  LL_fit  <- .LL(p=p, y=y) 
  LL_null <- .LL(p=mean(y), y=y)
  
  cox <- 1-exp(-(LL_fit-LL_null)*2/length(y)) 
  cox_max <- 1 - exp(2 * length(y) ^ (-1) * LL_null)
  c("cox"=cox)
  
}

## AUC (C-statistic)
fastAUC <- function(p, y) {
  x1 = p[y==1]; n1 = length(x1); 
  x2 = p[y==0]; n2 = length(x2);
  r = rank(c(x1,x2))  
  auc = (sum(r[1:n1]) - n1*(n1+1)/2) / n1 / n2
  return(auc)
}
