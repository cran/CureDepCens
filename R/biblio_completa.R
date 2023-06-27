
###--------------------------------------------------------------------------------------------------------------------------------------###
### Function to generate the frailty terms (w_k) of the complete log-conditional, where k denote the group index.                        ###
###--------------------------------------------------------------------------------------------------------------------------------------###
log_cond_wi <- function(w_k_aux,risco_a_T, risco_a_C, delta_T,delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta){

  w_k <- log(w_k_aux/(1-w_k_aux))

  if (ncol(t(t(X_cure))) == 1){
    pred_cure <- exp((t(t(X_cure))*beta_cure)+w_k)
  }
  if (ncol(t(t(X_cure))) != 1){
    pred_cure <- exp((X_cure%*%beta_cure)+w_k)
  }
  if (ncol(t(t(X_C))) == 1){
    exp_C <- exp((t(t(X_C))*beta_C)+alpha*w_k)
  }
  if (ncol(t(t(X_C))) != 1){
    exp_C <- exp((X_C%*%beta_C)+alpha*w_k)
  }

  log_vero_w <- sum(delta_T*(w_k) - pred_cure*(1- exp(-risco_a_T)) + delta_C*alpha*w_k - risco_a_C*exp_C)-((w_k^2)/(2*theta)) - log(w_k_aux*(1-w_k_aux))

  return(log_vero_w)
}
###--------------------------------------------------------------------------------------------------------------------------------------
support_wi <-  function(w_k_aux,risco_a_T, risco_a_C, delta_T,delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta){(w_k_aux>0)*(w_k_aux<1)}
###--------------------------------------------------------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------------------------------------------------------###
###                                                                                                                                      ###
### Piecewise exponential model                                                                                                          ###
###                                                                                                                                      ###
###--------------------------------------------------------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------------------------------------------------------###
### Function with the log-likelihood function, failure part, used in optim to estimate the parameters of the piecewise exponential       ###
### distribution and the cure fraction, considering covariates and frailty in the cure fraction                                          ###
###--------------------------------------------------------------------------------------------------------------------------------------###
modelo_param_T_MEP <-  function(param_T,delta_t, xi_T,X_cure,num_int_T,bi){  #lambda_T= vector of param_T
  colcura <- ncol(X_cure)
  n <- length(X_cure[,1])
  Delta <- delta_t%*%t(rep(1,num_int_T))
  lambda <- t((param_T[1:num_int_T])%*%t(rep(1,n)))  #estimate the vector of lambdas_j

  if (ncol(t(t(X_cure))) == 1){
    x_cure_b <- t(t(X_cure))*param_T[(num_int_T+1):(num_int_T+colcura)]   #estimate the vector of betas
  }
  if (ncol(t(t(X_cure))) != 1){
    x_cure_b <- X_cure%*%param_T[(num_int_T+1):(num_int_T+colcura)]   #estimate the vector of betas
  }

  prop <- exp(x_cure_b)*rowMeans(exp(bi))
  aux <- apply(exp(lambda)*xi_T,1,sum)
  soma <- aux%*%t(rep(1,num_int_T))
  aux2 <- as.vector(prop)*(1- exp(-apply(exp(lambda)*xi_T,1,sum)))
  aux3 <- aux2%*%t(rep(1,num_int_T))
  mean_bi <- rowMeans(bi)%*%t(rep(1,num_int_T))

  U <- -sum(Delta*(as.vector(x_cure_b) + lambda  + mean_bi - soma) - aux3)
}
###--------------------------------------------------------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------------------------------------------------------###
### Function with the log-likelihood function used in multiroot() to estimate vector of beta_C and alpha parameters of the dependent     ###
### censoring times.                                                                                                                     ###
###--------------------------------------------------------------------------------------------------------------------------------------###
modelo_C_MEP <-  function(beta_C, X_C,delta_c,risco_a_C,bi,n_intMC){
  q <- ncol(X_C)
  n <- nrow(X_C)
  if(ncol(t(t(X_C))) == 1){
    w_kl_beta_C <- risco_a_C*exp((X_C[,]*beta_C[1:q]))*(rowSums(exp(beta_C[q+1]*bi[,]))/n_intMC)
    w_kl_beta_C_alpha_num <- risco_a_C*exp(X_C[,]*beta_C[1:q])*(rowSums(bi[,]*exp(beta_C[q+1]*bi[,]))/n_intMC)
  }
  else{
    w_kl_beta_C <- risco_a_C*exp((X_C[,]%*%beta_C[1:q]))*(rowSums(exp(beta_C[q+1]*bi[,]))/n_intMC)
    w_kl_beta_C_alpha_num <- risco_a_C*exp(X_C[,]%*%beta_C[1:q])*(rowSums(bi[,]*exp(beta_C[q+1]*bi[,]))/n_intMC)
  }
  w_kl_beta_C_num <- cbind(matrix(w_kl_beta_C, nrow = n, ncol = q))*X_C
  U_C_1 <- colSums((X_C*delta_c - w_kl_beta_C_num))

  U_C_alpha <- sum(delta_c*(rowSums(bi[,])/n_intMC) - w_kl_beta_C_alpha_num)

  c(U_C_1 = U_C_1,U_C_alpha=U_C_alpha)
}
###--------------------------------------------------------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------------------------------------------------------###
### Function that partitions the time axis into intervals, also returns the interval indicator (observations of each interval)           ###
###--------------------------------------------------------------------------------------------------------------------------------------###
time_grid <- function(time, event, n_int=NULL)
{
  o <- order(time)
  time <- time[o]
  event <- event[o]
  time_aux <- unique(time[event==1])
  if(is.null(n_int))
  {
    n_int <- length(time_aux)
  }

  m <- length(time_aux)
  if(n_int > m)
  {
    a <- c(0,unique(time[event==1]))
    a[length(a)] <- Inf
  }
  else
  {
    b <- min(m,n_int)
    k1 <- trunc(m/b)
    r <- m-b*k1
    k2 <- k1+1
    idf1 <- seq(k1,(b-r)*k1, k1)
    idf2 <- sort(seq(m,max(idf1),-k2))
    idf <- unique(c(idf1,idf2))
    a_inf <- c(0,time_aux[idf])
    a_inf[length(a_inf)] <- Inf
    a_s_inf  <- c(0,time_aux[idf])
  }
  saida <- list(a_inf,a_s_inf)

  return(saida)
}
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Function with the calculation of first order derivatives for the piecewise exponential model.                                        ###
### Returns the expected array of first-order derivatives.                                                                               ###
### Calculation is done as in the article by Louis et al (1982)                                                                          ###
###--------------------------------------------------------------------------------------------------------------------------------------###
Esp_Deriv_Prim_Ordem_MEP <-  function( X_cure, X_C,delta_T,delta_C,id_T, beta_cure, beta_C, alpha, theta,nu_ind_j_T,nu_ind_j_C, lambda_T_j,
                                       lambda_C_j, deltaij_T,deltaij_C, w_k_grupo, ident){

  nc <- nrow(X_cure)
  p <- ncol(X_cure)
  q <- ncol(X_C)
  lambdalengthJ <- length(lambda_T_j)
  n <- nrow(X_C)
  wk = w_k_grupo[ident,]
  m <- nrow(w_k_grupo)
  deriv1 <- matrix(NA,(p+q+2+length(lambda_T_j)+length(lambda_C_j)),ncol(wk))  #matrix that with the first order derivatives

  if(ncol(t(t(X_cure))) == 1){
    pred_cure <- as.vector(exp(t(t(X_cure))*beta_cure))*exp(wk)
  }
  else{
    pred_cure <- as.vector(exp(X_cure%*%beta_cure))*exp(wk)
  }
  if(ncol(t(t(X_C))) == 1){
    pred_C <- as.vector(exp(t(t(X_C))*beta_C))*exp(alpha*wk)
  }
  else{
    pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  }

  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
  lambda_t <- t(lambda_T_j%*%t(rep(1,n)))
  lambda_c <- t(lambda_C_j%*%t(rep(1,n)))
  risco_a_T <- apply(t(as.vector(lambda_T_j)*t(deltaij_T)),1,sum)
  risco_a_C <- apply(t(as.vector(lambda_C_j)*t(deltaij_C)),1,sum)

  for (i in 1:p){
    deriv1[i,] <- colSums(X_cure[,i]*(Delta_t - pred_cure*(1-exp(-risco_a_T))))  #derived of beta_T's
  }

  for ( j in 1:length(lambda_T_j)){
    deriv1[(j+p),] <- nu_ind_j_T[j]*(lambda_T_j[j]^(-1)) - sum(delta_T[id_T==j]*deltaij_T[id_T==j,j]) - colSums(pred_cure*exp(-risco_a_T)*deltaij_T[,j]) #derived of lambda_T_j's
  }

  for (i in 1:q){
    deriv1[(length(lambda_T_j)+p+i),] <- colSums(X_C[,i]*(Delta_c - risco_a_C*pred_C))  #derived of beta_C's
  }

  deriv1[(length(lambda_T_j)+p+q+1),] <- colSums(wk*(Delta_c - risco_a_C*pred_C))  #derived of alpha

  for ( j in 1:length(lambda_C_j)){
    deriv1[(length(lambda_T_j)+p+q+j+1),] <- nu_ind_j_C[j]*(lambda_C_j[j]^(-1)) - colSums(deltaij_C[,j]*pred_C)    #derived of lambda_C_j's
  }

  deriv1[length(lambda_C_j)+length(lambda_T_j)+p+q+2,] <- -0.5*(m*theta^(-1) - (theta^(-2))*colSums(w_k_grupo^(2)))  #derivada de tau

  aux <- deriv1[,1]%*%t(deriv1[,1])

  for( i in 2:ncol(wk)){
    aux <- aux + deriv1[,i]%*%t(deriv1[,i])
  }

  return(aux/ncol(wk))
}
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Function with the calculation of second order derivatives for the piecewise exponential model.                                       ###
### Returns the expected array of second-order derivatives.                                                                              ###
### Calculation is done as in the article by Louis et al (1982).                                                                         ###
###--------------------------------------------------------------------------------------------------------------------------------------###
Esp_DerivParciais_MEP  <-  function( X_cure, X_C,delta_T,delta_C, beta_cure, beta_C, alpha, theta,nu_ind_j_T,nu_ind_j_C, lambda_T_j,lambda_C_j,
                                     deltaij_T,deltaij_C, w_k_grupo, ident){
  nc <- nrow(X_cure)
  p <- ncol(X_cure)
  q <- ncol(X_C)

  n <- nrow(X_C)
  wk = w_k_grupo[ident,]
  m <- nrow(w_k_grupo)
  deriv2 <- matrix(0,(p+q+2+length(lambda_T_j)+length(lambda_C_j)),(p+q+2+length(lambda_T_j)+length(lambda_C_j)))

  if (ncol(t(t(X_cure))) == 1){
    pred_cure <- as.vector(exp(t(t(X_cure))*beta_cure))*exp(wk)
  }
  if (ncol(t(t(X_cure))) != 1){
    pred_cure <- as.vector(exp(X_cure%*%beta_cure))*exp(wk)
  }
  if (ncol(t(t(X_C))) == 1){
    pred_C <- as.vector(exp(t(t(X_C))*beta_C))*exp(alpha*wk)
  }
  if (ncol(t(t(X_C))) != 1){
    pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  }

  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
  lambda_t <- t(lambda_T_j%*%t(rep(1,n)))
  lambda_c <- t(lambda_C_j%*%t(rep(1,n)))
  risco_a_T <- apply(t(as.vector(lambda_T_j)*t(deltaij_T)),1,sum)
  risco_a_C <- apply(t(as.vector(lambda_C_j)*t(deltaij_C)),1,sum)

  for (i in 1:p){
    for (u in 1:p) {
      deriv2[u,i] <- mean(colSums( - X_cure[,u]*X_cure[,i]*pred_cure*(1-exp(-risco_a_T))))   #derivative of beta_T's from the derivative of beta_T's
      deriv2[i,u] <- deriv2[u,i]
    }
  }

  for (i in 1:p) {
    for (j in 1:length(lambda_T_j)){
      deriv2[i,(j+p)] <- - mean(colSums(X_cure[,i]*pred_cure*exp(-risco_a_T)*deltaij_T[,j]))  #derivative of vector lambda_T_j from the derivative of beta_T
      deriv2[(j+p),i] <- deriv2[i,(j+p)]
    }
  }

  for (j in 1:length(lambda_T_j)){
    deriv2[(j+p),(j+p)] <- - nu_ind_j_T[j]*(lambda_T_j[j]^(-2)) + mean(colSums(pred_cure*exp(-risco_a_T)*deltaij_T[,j]*deltaij_T[,j]))   #derivative of vector lambda_T_j from the derivative of lambda_T_j
  }

  for (i in 1:q){
    for (u in 1:q) {
      deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+u)] <- - mean(colSums(X_C[,i]*X_C[,u]*risco_a_C*pred_C))  # derivative of beta_C's from the derivative of beta_C's
      deriv2[(length(lambda_T_j)+p+u),(length(lambda_T_j)+p+i)] <- deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+u)]
    }
  }

  for (i in 1:q){
    deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+q+1)]  <- - mean(colSums(X_C[,i]*wk*risco_a_C*pred_C)) # derivative of alpha from the derivative of beta_C's
    deriv2[(length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+q+i)]  <- deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+q+1)]
  }

  deriv2[(length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+q+1)] <- - mean(colSums(wk*wk*risco_a_C*pred_C))  # derivative of alpha from the derivative of alpha

  for (j in 1:length(lambda_C_j)){
    deriv2[(j+length(lambda_T_j)+p+q+1),(j+length(lambda_T_j)+p+q+1)] <- - nu_ind_j_C[j]*(lambda_C_j[j]^(-2))  # derivative of vector lambda_C_j from the derivative of vector lambda_C_j
  }

  for (i in 1:q) {
    for (j in 1:length(lambda_C_j)){
      deriv2[(length(lambda_T_j)+p+i),(j+length(lambda_T_j)+p+q+1)] <- - mean(colSums(deltaij_C[,j]*pred_C*X_C[,i]))
      deriv2[(j+length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+i)] <- deriv2[(length(lambda_T_j)+p+i),(j+length(lambda_T_j)+p+q+1)]   # derivative of beta_C's from the derivative of vector lambda_C_j
    }
  }

  for (j in 1:length(lambda_C_j)){
    deriv2[(length(lambda_T_j)+p+q+1),(j+length(lambda_T_j)+p+q+1)] <- - mean(colSums(deltaij_C[,j]*pred_C*wk))   # derivative of alpha from the derivative of lambda_C_j
    deriv2[(j+length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+q+1)] <- deriv2[(length(lambda_T_j)+p+q+1),(j+length(lambda_T_j)+p+q+1)]
  }

  deriv2[(length(lambda_T_j)+length(lambda_C_j)+p+q+2),(length(lambda_T_j)+length(lambda_C_j)+p+q+2)] <-  -0.5*(-m*(theta^(-2)) + (theta^(-3))*2*mean(colSums(w_k_grupo^(2))))   # derivative of theta from the derivative of theta

  return((as.matrix(deriv2)))
}
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Function with the calculation of model selection criteria, considering a piecewise exponential model                                 ###                                                                            ###
###--------------------------------------------------------------------------------------------------------------------------------------###
crit_MEP <- function(X_cure,X_C,bi,n,beta_cure,lambda_T_j,beta_C,alpha,lambda_C_j,delta_t,delta_c,risco_a_T,risco_a_C,id_T,id_C,m,ident) {
  w <- bi  # frailty matrix saved in the last iteration
  L <- ncol(w)
  num_param <- length(c(beta_cure,lambda_T_j,beta_C,alpha,lambda_C_j)) # number of parameters, parameters saved in fit (fit <- c((out[s,]+out[s-1,]+out[s-2,])/3))
  log_vero <- matrix(NA,n,L)  # log-likelihood, L is the number of Monte Carlo replicates of frailty w

  for (l in 1:L){
    if (ncol(t(t(X_cure))) == 1 && ncol(t(t(X_C))) == 1){

      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_cure*beta_cure + w[,l] - risco_a_T) - exp(X_cure*beta_cure + w[,l])*(1-exp(-risco_a_T))
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C*beta_C + alpha*w[,l]) - risco_a_C*exp(X_C*beta_C + alpha*w[,l])
    }

    if (ncol(t(t(X_cure))) == 1 && ncol(t(t(X_C))) != 1){

      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_cure*beta_cure + w[,l] - risco_a_T) - exp(X_cure*beta_cure + w[,l])*(1-exp(-risco_a_T))
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C%*%beta_C + alpha*w[,l]) - risco_a_C*exp(X_C%*%beta_C + alpha*w[,l])
    }

    if (ncol(t(t(X_cure))) != 1 && ncol(t(t(X_C))) == 1){

      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_cure%*%beta_cure + w[,l] - risco_a_T) - exp(X_cure%*%beta_cure + w[,l])*(1-exp(-risco_a_T))
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C*beta_C + alpha*w[,l]) - risco_a_C*exp(X_C*beta_C + alpha*w[,l])

    }
    if (ncol(t(t(X_cure))) != 1 && ncol(t(t(X_C))) != 1){

      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_cure%*%beta_cure + w[,l] - risco_a_T) - exp(X_cure%*%beta_cure + w[,l])*(1-exp(-risco_a_T))
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C%*%beta_C + alpha*w[,l]) - risco_a_C*exp(X_C%*%beta_C + alpha*w[,l])
    }
  }

  vero <- exp(log_vero)
  vero_grupo <- matrix(NA,m,L)  #m is the number of groups/cluster

  for (i in 1:m){
    vero_grupo[i,] <- matrixStats::colProds(vero, rows = which(ident==i))  # likelihood for each group
  }

  log_lik <- sum(log(rowSums(vero_grupo)/L))

  AIC <- 2*(-log_lik + num_param)
  BIC <- 2*(-log_lik + 0.5*log(n)*num_param)
  HQ <- 2*(-log_lik + log(log(n))*num_param)

  return(cbind(AIC,BIC,HQ))
}
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Function with all steps of the piecewise exponential model maximization, considering EM-Monte Carlo algorithm in the fit             ###
###--------------------------------------------------------------------------------------------------------------------------------------###
model_MEP_dep <-  function(formula, data, delta_t, delta_c, ident, Num_intervals){

  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- stats::model.matrix(formula, data = mf, rhs = 2)
  X_cura <- Z[,-1]
  X_cura <- as.matrix(X_cura)
  X_cura <- cbind(1,X_cura)
  X_C <- X[,-1]
  X_C <- as.matrix(X_C)
  Xlabels <- colnames(X_cura)
  Zlabels <- colnames(X_C)
  time <- stats::model.response(mf)
  t <- stats::model.response(mf)
  q <- ncol(X_C)
  p <- ncol(X_cura)
  n <- nrow(data)
  m <- max(ident)

  ###--------------------------------------------------------------------------------------
  #initial value for parameters
  beta_cura <- rep(0.1,p)
  beta_C <- rep(0.1,q)
  theta <- 1
  alpha <- 0
  param <- c(alpha, beta_cura,beta_C,theta)
  risco_a_T <- rep(1,n)
  risco_a_C <- rep(1,n)
  bmax <- 5
  lambda_T_j <- rep(0.1,bmax)
  ###----------------------------------------------------------------------------------------------------
  # EMMC Algorithm Specifications
  maxit <- 150
  eps1= rep(1e-7, length(param))
  eps2= rep(1e-8, length(param))
  n_intMCs = c(rep(10,20),rep(20,30),rep(30,50),rep(50,30),rep(75,10),rep(100,10))

  ## Initializing the objects used
  out =  matrix(NA,maxit+1,length(param))
  dif =  matrix(NA,maxit+1,length(param))
  final = length(param)
  count = rep(0,maxit+1)
  s=1
  continue=TRUE

  while (continue == TRUE) {

    count = rep(0,maxit+1)
    out[s,] =  c(alpha, beta_cura,beta_C,theta)
    n_intMC = n_intMCs[s]
    m <- max(ident)
    w_chapeu_grupo <- matrix(NA, m, n_intMC)

    for ( k in 1:m){  #k eh grupo
      w_trans <- arms(0.5, myldens=log_cond_wi, indFunc=support_wi,n.sample=n_intMC,risco_a_T=risco_a_T[ident==k], risco_a_C=risco_a_C[ident==k], delta_T=delta_t[ident==k],delta_C=delta_c[ident==k], X_cure=X_cura[ident==k,],X_C=X_C[ident==k,], beta_cure=beta_cura,beta_C=beta_C,alpha=alpha, theta=theta)
      w_auxi <- log(w_trans)-log(1-w_trans)
      w_chapeu_grupo[k,] <- w_auxi
    }
    bi = w_chapeu_grupo[ident,]
    theta <- mean(w_chapeu_grupo^2)
    ###----------------------------------------------------------------------------------------
    a <- time_grid(t,delta_t,bmax); a <- a[[1]]; #failure times grid
    c <- time_grid(t,delta_c,bmax); c <- c[[1]]; #dependent censoring times grid
    b <- length(a)-1
    d <- length(c)-1
    a[length(a)] <- max(t)
    c[length(c)] <- max(t)
    num_int_T <- bmax
    num_int_C <- bmax

    id_T <- as.numeric(cut(t,a))
    nu_T <- tapply(delta_t,id_T ,sum)

    xi.falha <- matrix(0,n,b)
    xi.cens <- matrix(0,n,d)

    for(i in 1:n){
      for(j in 1:b){
        xi.falha[i,j] <- (min(t[i], a[j+1])-a[j])*((t[i]-a[j])>0)
      }
    }
    id.falha <- as.numeric(cut(t,a))
    id.cens  <- as.numeric(cut(t,c))
    for(i in 1:n){
      for(j in 1:d){
        xi.cens[i,j] <- (min(t[i], c[j+1])-c[j])*((t[i]-c[j])>0)
      }
    }
    param_T <- optim(par = c(lambda_T_j, beta_cura), fn = modelo_param_T_MEP, control = list(maxit = 5000), method = c( "BFGS"), delta_t=delta_t, xi_T=xi.falha,X_cure=X_cura,num_int_T=bmax, bi=bi)
    param_Ts <- param_T$par
    lambda_T_j <- exp(param_Ts[1:num_int_T])
    beta_cura  <- param_Ts[(num_int_T+1):(num_int_T+p)]
    risco_a_T <- apply(t(lambda_T_j*t(xi.falha)),1,sum)
    ###--------------------------------------------------------------------------------
    pred_linear_C <- exp(X_C%*%beta_C)*rowMeans(exp(alpha*bi[,]))
    id_C <- as.numeric(cut(t,c))
    nu_C <- tapply(delta_c,id_C ,sum)

    xi.cens_pred <- matrix(0,n,d)
    for(i in 1:n){
      for(j in 1:d){
        xi.cens_pred[i,j] <- (min(t[i], c[j+1])-c[j])*((t[i]-c[j])>0)*pred_linear_C[i]
      }
    }
    lambda_C_j <- nu_C/apply(xi.cens_pred,2,sum)
    risco_a_C <- apply(t(as.vector(lambda_C_j)*t(xi.cens)),1,sum)

    S_C <- multiroot(f = modelo_C_MEP, start = c(beta_C,alpha), X_C=X_C,delta_c=delta_c,risco_a_C=risco_a_C,bi=bi,n_intMC=n_intMC)
    betas_C <- S_C$root
    beta_C <- betas_C[1:q]
    alpha <- betas_C[q+1]
    ###-----------------------------------------------------------------------------------
    # stopping criterion
    out[s+1,]= c(alpha, beta_cura,beta_C,theta)
    #print(out[s+1,])
    if (s>2) {
      dif[s,] <- (abs(out[s+1,]-out[s,]))/(abs(out[s,])-eps1)

      for (z in 1: length(param)){
        if (dif[s,z]<eps2[z]) {
          count[s] = count[s] + 1
        }
      }
    }
    s=s+1
    if (s>3) {
      continue=((s <= maxit) & (count[s] < length(param))&(count[s-1] < length(param))&(count[s-2] < length(param)))
    }
  } ## end of EMMC
  param_est <- c((out[s,]+out[s-1,]+out[s-2,])/3)

  ###-----------------------------------------------------------------------------------
  # calculation of the expectation of the derivatives for the calculation of the Standard Error, according to Louis et al (1982)
  Esp_deriv_ordem1 <- Esp_Deriv_Prim_Ordem_MEP( X_cure=X_cura, X_C=X_C,delta_T=delta_t,delta_C=delta_c,id_T=id_T, beta_cure=beta_cura, beta_C=beta_C, alpha=alpha,
                                                theta=theta,nu_ind_j_T=nu_T,nu_ind_j_C=nu_C, lambda_T_j=lambda_T_j,lambda_C_j=lambda_C_j,
                                                deltaij_T=xi.falha,deltaij_C=xi.cens, w_k_grupo=w_chapeu_grupo, ident=ident)

  Esp_deriv_ordem2 <- Esp_DerivParciais_MEP( X_cure=X_cura, X_C=X_C,delta_T=delta_t,delta_C=delta_c, beta_cure=beta_cura, beta_C=beta_C, alpha=alpha,
                                             theta=theta,nu_ind_j_T=nu_T,nu_ind_j_C=nu_C, lambda_T_j=lambda_T_j,lambda_C_j=lambda_C_j,
                                             deltaij_T=xi.falha,deltaij_C=xi.cens, w_k_grupo=w_chapeu_grupo, ident=ident)

  InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
  Var <- solve(InfFisher)
  ErroPadrao <- sqrt(diag(Var))
  if (any(is.na(ErroPadrao))) warning("The algorithm did not converge. It might converge if you run the function again.", call. = FALSE)
  ###-----------------------------------------------------------------------------------
  # p-value
  p_value_alpha <- 2*stats::pnorm(-abs(param_est[1]/ErroPadrao[(2*p+q+2)]))
  p_value_t <- vector()
  for (i in 1:p){
    p_value_t[i] <- 2*stats::pnorm(-abs((param_est[c(i)])/(ErroPadrao[c(1+i)])))
  }
  p_value_c <- vector()
  for (j in 1:q){
    p_value_c[j] <- 2*stats::pnorm(-abs((param_est[c(2+p+j)])/(ErroPadrao[c(2*p+1+j)])))
  }


  # calculation of model selection criteria
  criterios <- crit_MEP(X_cure=X_cura,X_C,bi,n,beta_cure=beta_cura,lambda_T_j,beta_C,alpha,lambda_C_j,delta_t=delta_t,delta_c=delta_c,risco_a_T,risco_a_C,id_T,id_C,m,ident)

  # Adjusting the class/list
  fit <- param_est
  fit <- list(fit=fit)
  fit$stde <- ErroPadrao
  fit$crit <- criterios
  fit$pvalue <- c(p_value_alpha, p_value_t, p_value_c)
  fit$n <- n
  fit$p <- p
  fit$q <- q
  fit$call <- match.call()
  fit$formula <- stats::formula(Terms)
  fit$terms <- stats::terms.formula(formula)
  fit$labels1 <- Zlabels
  fit$labels2 <- Xlabels
  fit$bmax <- bmax
  fit$risco_a_T <- risco_a_T
  fit$risco_a_C <- risco_a_C
  fit$bi <- bi
  fit$X_cura <- X_cura
  fit$X_C <- X_C
  fit$time <- time
  fit$t <- t
  class(fit) <- "dcensoring"
  return(fit)
}
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
###                                                                                                                                      ###
### Weibull model                                                                                                                        ###
###                                                                                                                                      ###
###--------------------------------------------------------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------------------------------------------------------###
### Function with the log-likelihood function, failure part, used in optim to estimate the parameters of the Weibull distribution        ###
### and the cure fraction, considering covariates and frailty in the cure fraction                                                       ###
###--------------------------------------------------------------------------------------------------------------------------------------###
modelo_param_T_Weib <-  function(par,t,delta_T, X_cure, bi){
  p <- ncol(X_cure)
  n <- nrow(X_cure)
  if(ncol(t(t(X_cure))) == 1){
    pred <- X_cure*par[3]
    pred_cure <- exp(pred)*rowMeans(exp(bi))
  }
  else{
    pred <- X_cure%*%par[3:(p+2)]
    pred_cure <- exp(pred)*rowMeans(exp(bi))
  }
  a <- par[1]
  b <- exp(par[2])
  exp_T <- (t^a)*b

  U <- -sum(delta_T*(pred + rowMeans(bi) + log(a) + (a-1)*log(t) + log(b) - exp_T)- pred_cure*(1-exp(-exp_T)))
}
###--------------------------------------------------------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------------------------------------------------------###
### Function with the log-likelihood function, dependent censoring part, used in optim to estimate the parameters of the Weibull         ###
### distribution, considering covariates in the hazard function of dependent censoring time C.                                           ###
###--------------------------------------------------------------------------------------------------------------------------------------###
modelo_param_C <-  function(par,t,delta_C, X_C, bi){

  q <- ncol(X_C)
  n <- nrow(X_C)
  if(ncol(t(t(X_C))) == 1){
    pred_C <- X_C*par[4]
  }
  else{
    pred_C <- X_C%*%par[4:(q+3)]
  }
  a <- par[1]
  b <- exp(par[2])
  exp_C <- (t^a)*b*exp(pred_C)

  U <- -sum(delta_C*(log(a) + (a-1)*log(t) + log(b) + pred_C + par[3]*rowMeans(bi))- exp_C*rowMeans(exp(par[3]*bi)))
}
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Functions with first order derivatives for the Weibull model                                                                         ###
###--------------------------------------------------------------------------------------------------------------------------------------###
Esp_Deriv_Prim_Ordem_Weib <-  function(t,delta_T, delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta, alpha_T, alpha_C, lambda_T,
                                       lambda_C, w_k_grupo, ident){

  p <- ncol(X_cure)
  q <- ncol(X_C)
  wk = w_k_grupo[ident,]
  m <- nrow(w_k_grupo)
  num_param <- length(beta_cure)+2+length(beta_C)+2+length(alpha)+length(theta)
  deriv1 <- matrix(NA,num_param,ncol(wk))
  if (ncol(t(t(X_cure))) == 1){
    pred_cure <- as.vector(exp(t(t(X_cure))*beta_cure))*exp(wk)
  }
  if (ncol(t(t(X_cure))) != 1){
    pred_cure <- as.vector(exp(X_cure%*%beta_cure))*exp(wk)
  }
  if (ncol(t(t(X_C))) == 1){
    pred_C <- as.vector(exp(t(t(X_C))*beta_C))*exp(alpha*wk)
  }
  if (ncol(t(t(X_C))) != 1){
    pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  }
  risco_a_T <- (t^alpha_T)*lambda_T
  risco_a_C <- (t^alpha_C)*lambda_C
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))

  for (i in 1:p){
    deriv1[i,] <- colSums(X_cure[,i]*(Delta_t - pred_cure*(1-exp(-risco_a_T)))) # derivative of vector beta_T's
  }

  deriv1[1+p,] <- colSums(Delta_t*(alpha_T^(-1) + log(t) - risco_a_T*log(t))- pred_cure*exp(-risco_a_T)*risco_a_T*log(t))   # derivative of alpha.t

  deriv1[2+p,] <- colSums(Delta_t*(lambda_T^(-1) - t^alpha_T) - pred_cure*exp(-risco_a_T)*(t^alpha_T))   # derivative of lambda.t

  for (i in 1:q){
    deriv1[2+p+i,] <- colSums(X_C[,i]*(Delta_c - risco_a_C*pred_C))  # derivative of vector beta_C's
  }

  deriv1[2+p+q+1,] <- colSums(wk*(Delta_c - risco_a_C*pred_C))  # derivative of alpha

  deriv1[2+p+q+2,] <- colSums(Delta_c*(alpha_C^(-1) + log(t)) - risco_a_C*log(t)*pred_C)   # derivative of alpha.c

  deriv1[2+p+q+3,] <- colSums(Delta_c*(lambda_C^(-1)) - (t^alpha_C)*pred_C)   # derivative of lambda.c

  deriv1[2+p+q+4,] <- -0.5*(m*theta^(-1) - (theta^(-2))*colSums(w_k_grupo^(2)))   # derivative of theta

  aux <- deriv1[,1]%*%t(deriv1[,1])

  for( i in 2:ncol(wk)){
    aux <- aux + deriv1[,i]%*%t(deriv1[,i])
  }

  return(aux/ncol(wk))
}
###--------------------------------------------------------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------------------------------------------------------###
### Function with second order partial derivatives for the Weibull model                                                                 ###
###--------------------------------------------------------------------------------------------------------------------------------------###
Esp_DerivParciais_Weib <-  function(t,delta_T, delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta, alpha_T, alpha_C, lambda_T, lambda_C,
                                    w_k_grupo, ident){

  p <- ncol(X_cure)
  q <- ncol(X_C)
  num_param <- length(beta_cure)+2+length(beta_C)+2+length(alpha)+length(theta)
  deriv2 <- matrix(0,num_param,num_param)
  wk = w_k_grupo[ident,]
  m <- nrow(w_k_grupo)
  if (ncol(t(t(X_cure))) == 1){
    pred_cure <- as.vector(exp(t(t(X_cure))*beta_cure))*exp(wk)
  }
  if (ncol(t(t(X_cure))) != 1){
    pred_cure <- as.vector(exp(X_cure%*%beta_cure))*exp(wk)
  }
  if (ncol(t(t(X_C))) == 1){
    pred_C <- as.vector(exp(t(t(X_C))*beta_C))*exp(alpha*wk)
  }
  if (ncol(t(t(X_C))) != 1){
    pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  }
  risco_a_T <- (t^alpha_T)*lambda_T
  risco_a_C <- (t^alpha_C)*lambda_C
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))

  for (i in 1:p){
    for (u in 1:p) {
      deriv2[u,i] <- mean(colSums( - X_cure[,u]*X_cure[,i]*pred_cure*(1-exp(-risco_a_T))))   #derivative of vector beta_cura from the derivative of vector beta_cura
      deriv2[i,u] <- deriv2[u,i]
    }
  }

  for (i in 1:p){
    deriv2[i,p+1] <- mean(colSums( - X_cure[,i]*pred_cure*exp(-risco_a_T)*log(t)*risco_a_T))   #derivative of alpha.t from derivative of vector beta_cura
    deriv2[p+1,i] <- deriv2[i,p+1]
  }

  for (i in 1:p){
    deriv2[i,p+2] <- mean(colSums( - X_cure[,i]*pred_cure*exp(-risco_a_T)*(t^alpha_T)))    #derivative of lambda.t from derivative of vector beta_cura
    deriv2[p+2,i] <- deriv2[i,p+2]
  }

  deriv2[p+1,p+1] <- mean(colSums(Delta_t*(-alpha_T^(-2) - risco_a_T*log(t)*log(t)) + (pred_cure*exp(-risco_a_T)*(-risco_a_T)*log(t)*log(t)*(-risco_a_T))))    #derivative of alpha.t from the derivative of alpha.t

  deriv2[p+1,p+2] <- mean(colSums(Delta_t*(-(t^alpha_T)*log(t)) + (pred_cure*exp(-risco_a_T)*(-risco_a_T)*(-t^alpha_T)*log(t))))  #derivative of lambda.t from the derivative of alpha.t
  deriv2[p+2,p+1] <- deriv2[p+1,p+2]

  deriv2[p+2,p+2] <- mean(colSums(Delta_t*(-lambda_T^(-2)) + pred_cure*exp(-risco_a_T)*(t^alpha_T)*(t^alpha_T))) #derivative of lambda.t from the derivative of lambda.t

  for (i in 1:q){
    for (u in 1:q) {
      deriv2[p+2+u,p+2+i] <- mean(colSums(- X_C[,u]*X_C[,i]*risco_a_C*pred_C))   #derivative of vector beta_C's from the derivative of beta_C's
      deriv2[p+2+i,p+2+u] <- deriv2[p+2+u,p+2+i]
    }
  }

  for (i in 1:q){
    deriv2[p+2+i,p+q+3] <- mean(colSums(- X_C[,i]*wk*risco_a_C*pred_C))   #derivative of alpha from the derivative of beta_C's
    deriv2[p+q+3,p+2+i] <- deriv2[p+2+i,p+q+3]
  }

  for (i in 1:q){
    deriv2[p+2+i,p+q+4] <- mean(colSums(- X_C[,i]*risco_a_C*log(t)*pred_C))   #derivative of alpha.c from the derivative of beta_C's
    deriv2[p+q+4,p+2+i] <- deriv2[p+2+i,p+q+4]
  }

  for (i in 1:q){
    deriv2[p+2+i,p+q+5] <- mean(colSums(- X_C[,i]*(t^alpha_C)*pred_C))   #derivative of lambda.c from the derivative of beta_C's
    deriv2[p+q+5,p+2+i] <- deriv2[p+2+i,p+q+5]
  }

  deriv2[p+q+3,p+q+3] <- mean(colSums(- risco_a_C*pred_C*wk*wk))  #derivative of alpha from the derivative of alpha

  deriv2[p+q+3,p+q+4] <- mean(colSums( - risco_a_C*log(t)*pred_C*wk))  #derivative of alpha.c from the derivative of alpha
  deriv2[p+q+4,p+q+3] <-  deriv2[p+q+3,p+q+4]

  deriv2[p+q+3,p+q+5] <- mean(colSums(- (t^alpha_C)*pred_C*wk))  #derivative of lambda.c from the derivative of alpha
  deriv2[p+q+5,p+q+3] <-  deriv2[p+q+3,p+q+5]

  deriv2[p+q+4,p+q+4] <- mean(colSums(Delta_c*(-alpha_C^(-2)) - risco_a_C*pred_C*log(t)*log(t)))   #derivative of alpha.c from the derivative of alpha.c

  deriv2[p+q+4,p+q+5] <- mean(colSums( - (t^alpha_C)*log(t)*pred_C))   #derivative of lambda.c from the derivative of alpha.c
  deriv2[p+q+5,p+q+4] <- deriv2[p+q+4,p+q+5]

  deriv2[p+q+5,p+q+5] <- mean(colSums(Delta_c*(-lambda_C^(-2))))  #derivative of lambda.c from the derivative of lambda.c

  deriv2[p+q+6,p+q+6] <- mean(-0.5*(-m*(theta^(-2)) + (theta^(-3))*2*colSums(w_k_grupo^(2))))  #derivative of theta from the derivative of theta

  return((as.matrix(deriv2)))
}
###--------------------------------------------------------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------------------------------------------------------###
### Function for calculating model selection criteria, considering fit with the Weibull model                                            ###
###--------------------------------------------------------------------------------------------------------------------------------------###
crit_weibull <- function(X_cure,X_C,bi,n,alpha_T,lambda_T,beta_cura,alpha_C,lambda_C,beta_C,alpha,delta_t,delta_c,time,m,ident) {
  w <- bi      #frailty matrix saved in the last iteration
  L <- ncol(w)
  num_param <- length(c(alpha,beta_cura,beta_C,alpha_T,alpha_C,lambda_T,lambda_C)) # number of parameters, parameters saved in fit (fit <- c((out[s,]+out[s-1,]+out[s-2,])/3))
  log_vero <- matrix(NA,n,L)  #log-likelihood, L number of Monte Carlo replicates of frailty w

  pred <- X_cure%*%beta_cura
  a_T <- alpha_T
  b_T <- lambda_T
  exp_T <- (time^a_T)*b_T

  pred_C <- X_C%*%beta_C
  a_C <- alpha_C
  b_C <- lambda_C
  exp_C <- (time^a_C)*b_C*exp(pred_C)

  for ( l in 1:L){
    log_vero[,l] <- delta_t*(pred + w[,l] + log(a_T) + (a_T-1)*log(time) + log(b_T) - exp_T) - exp(pred)*exp(w[,l])*(1-exp(-exp_T))
    + delta_c*(log(a_C) + (a_C-1)*log(time) + log(b_C) + pred_C + alpha*w[,l])- exp_C*exp(alpha*w[,l])
  }

  vero <- exp(log_vero)
  vero_grupo <- matrix(NA,m,L)  #m is the number of groups/cluster

  for (i in 1:m){
    vero_grupo[i,] <- matrixStats::colProds(vero, rows = which(ident==i))  # likelihood for each group
  }

  log_lik <- sum(log(rowSums(vero_grupo)/L))

  AIC <- 2*(-log_lik + num_param)

  BIC <- 2*(-log_lik + 0.5*log(n)*num_param)

  HQ <- 2*(-log_lik + log(log(n))*num_param)

  return(cbind(AIC,BIC,HQ))
}
###--------------------------------------------------------------------------------------------------------------------------------------###



###--------------------------------------------------------------------------------------------------------------------------------------###
### Function with all steps of the Weibull model maximization, considering EM-Monte Carlo algorithm in the fit                           ###
###--------------------------------------------------------------------------------------------------------------------------------------###
model_Weibull_dep <-  function(formula, data, delta_t, delta_c, ident){

  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- stats::model.matrix(formula, data = mf, rhs = 2)
  X_cura <- Z[,-1]
  X_cura <- as.matrix(X_cura)
  X_cura <- cbind(1,X_cura)
  X_C <- X[,-1]
  X_C <- as.matrix(X_C)
  Xlabels <- colnames(X_cura)
  Zlabels <- colnames(X_C)
  time <- stats::model.response(mf)
  t <- stats::model.response(mf)
  q <- ncol(X_C)
  p <- ncol(X_cura)
  n <- nrow(data)
  m <- max(ident)
  ###--------------------------------------------------------------------------------------
  # initial value
  beta_cura <- rep(0.1,p)
  beta_C <- rep(0.1,q)
  alpha_T <- 1
  lambda_T <- 1
  alpha_C <- 1
  lambda_C <- 1
  theta <- 1
  alpha <- 0
  param <- c(alpha, beta_cura,beta_C,alpha_T,alpha_C,lambda_T,lambda_C,theta)
  risco_a_T <- rep(1,n)
  risco_a_C <- rep(1,n)
  ###----------------------------------------------------------------------------------------------------
  #EMMC Algorithm Specifications
  maxit <- 150 #maximum number of iterations
  eps1= rep(1e-7, length(param))
  eps2= rep(1e-8, length(param))
  n_intMCs = c(rep(10,20),rep(20,30),rep(30,50),rep(50,30),rep(75,10),rep(100,10))

  out =  matrix(NA,maxit+1,length(param))
  dif =  matrix(NA,maxit+1,length(param))
  final = length(param)
  count = rep(0,maxit+1)
  s=1
  continue=TRUE

  while (continue == TRUE) {

    count = rep(0,maxit+1)
    out[s,] =  c(alpha, beta_cura,beta_C,alpha_T,alpha_C,lambda_T,lambda_C,theta)
    n_intMC = n_intMCs[s]
    m <- max(ident)
    w_chapeu_grupo <- matrix(NA, m, n_intMC)

    for ( k in 1:m){  #k eh grupo
      w_trans <- arms(0.4, myldens=log_cond_wi, indFunc=support_wi,n.sample=n_intMC,risco_a_T=risco_a_T[ident==k], risco_a_C=risco_a_C[ident==k], delta_T=delta_t[ident==k],delta_C=delta_c[ident==k], X_cure=X_cura[ident==k,],X_C=X_C[ident==k,], beta_cure=beta_cura,beta_C=beta_C,alpha=alpha, theta=theta)
      w_auxi <- log(w_trans)-log(1-w_trans)
      w_chapeu_grupo[k,] <- w_auxi
    }
    bi = w_chapeu_grupo[ident,]
    theta <- mean(w_chapeu_grupo^2)
    ###----------------------------------------------------------------------------------------
    param_T <- optim(par = c(alpha_T,lambda_T, beta_cura), fn = modelo_param_T_Weib, control = list(maxit = 5000), method = c( "Nelder-Mead"),t=t,delta_T=delta_t, X_cure=X_cura, bi=bi)
    par_T <- param_T$par
    alpha_T <- par_T[1]
    lambda_T <- exp(par_T[2])
    beta_cura <- par_T[3:(2+p)]
    risco_a_T <- (t^alpha_T)*(lambda_T)
    ###--------------------------------------------------------------------------------
    param_C <- optim(par = c(alpha_C,lambda_C,alpha,beta_C), fn = modelo_param_C, control = list(maxit = 5000), method = c( "Nelder-Mead"),t=t,delta_C=delta_c, X_C=X_C, bi=bi)
    par_C <- param_C$par
    alpha_C <- par_C[1]
    lambda_C <- exp(par_C[2])
    alpha <- par_C[3]
    beta_C <- par_C[4:(3+q)]
    risco_a_C <- (t^alpha_C)*(lambda_C)
    ###-----------------------------------------------------------------------------------
    # stopping criterion
    out[s+1,]= c(alpha, beta_cura,beta_C,alpha_T,alpha_C,lambda_T,lambda_C,theta)
    #print(out[s+1,])
    if (s>2) {
      dif[s,] <- (abs(out[s+1,]-out[s,]))/(abs(out[s,])-eps1)

      for (z in 1: length(param)){
        if (dif[s,z]<eps2[z]) {
          count[s] = count[s] + 1
        }
      }
    }
    s=s+1
    if (s>3) {
      continue=((s <= maxit) & (count[s] < length(param))&(count[s-1] < length(param))&(count[s-2] < length(param)))
    }
  } ## end of EMMC
  param_est <- c((out[s,]+out[s-1,]+out[s-2,])/3)

  # standard error calculation
  Esp_deriv_ordem1 <- Esp_Deriv_Prim_Ordem_Weib(t,delta_T=delta_t, delta_C=delta_c, X_cure=X_cura, X_C=X_C,
                                                beta_cure=beta_cura, beta_C=beta_C, alpha=alpha, theta=theta,
                                                alpha_T=alpha_T, alpha_C=alpha_C, lambda_T=lambda_T, lambda_C=lambda_C,
                                                w_k_grupo=w_chapeu_grupo, ident=ident)

  Esp_deriv_ordem2 <- Esp_DerivParciais_Weib(t,delta_T=delta_t, delta_C=delta_c, X_cure=X_cura, X_C=X_C,
                                             beta_cure=beta_cura, beta_C=beta_C, alpha=alpha, theta=theta,
                                             alpha_T=alpha_T, alpha_C=alpha_C, lambda_T=lambda_T, lambda_C=lambda_C,
                                             w_k_grupo=w_chapeu_grupo, ident=ident)

  InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
  Var <- solve(InfFisher)
  ErroPadrao <- sqrt(diag(Var))
  if (any(is.na(ErroPadrao))) warning("The algorithm did not converge. It might converge if you run the function again.", call. = FALSE)

  p_value_alpha <- 2*stats::pnorm(-abs(param_est[1]/ErroPadrao[(3+p+q)]))
  p_value_t <- vector()
  for (i in 1:p){
    p_value_t[i] <- 2*stats::pnorm(-abs((param_est[c(1+i)])/(ErroPadrao[c(i)])))
  }
  p_value_c <- vector()
  for (j in 1:q){
    p_value_c[j] <- 2*stats::pnorm(-abs((param_est[c(2+p+j)])/(ErroPadrao[c(2+p+j)])))
  }

  # criteria calculation
  criterios <- crit_weibull(X_cure=X_cura,X_C,bi,n,alpha_T,lambda_T,beta_cura,alpha_C,lambda_C,beta_C,alpha,delta_t=delta_t,delta_c=delta_c,time=t,m,ident)

  # Adjusting the class/list
  fit <- param_est
  fit <- list(fit=fit)
  fit$stde <- ErroPadrao
  fit$crit <- criterios
  fit$pvalue <- c(p_value_alpha, p_value_t, p_value_c)
  fit$n <- n
  fit$p <- p
  fit$q <- q
  fit$call <- match.call()
  fit$formula <- stats::formula(Terms)
  fit$terms <- stats::terms.formula(formula)
  fit$labels1 <- Zlabels
  fit$labels2 <- Xlabels
  fit$risco_a_T <- risco_a_T
  fit$risco_a_C <- risco_a_C
  fit$bi <- bi
  fit$X_cura <- X_cura
  fit$X_C <- X_C
  fit$time <- time
  fit$t <- t
  class(fit) <- "dcensoring"
  return(fit)
}
