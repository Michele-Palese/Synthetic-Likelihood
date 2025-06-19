
sigmoid <- function(x) 1/(1 + exp(-x))
softmax <- function(hidden) exp(hidden)/rowSums(exp(hidden)) 
dim_input <- 1
k <- 3
hidden_size <- 5
# Verosimiglianza ed intervalli 
logL <- function(theta, dati) -sum(dgamma(dati, theta,
                                          1, log = T))
simulatore <- function(theta, n_sim = 20)rgamma(n_sim, theta,1)
studio_simulazione <- function(numero_simulazioni = 500)
{
  TRUE_THETA <- numeric(numero_simulazioni)
  Y_OSS <- matrix(NA,nrow = numero_simulazioni, ncol = 20)
  out <- numeric(numero_simulazioni)
  intervalli <- matrix(NA, nrow = numero_simulazioni, ncol = 2)
  for (i in 1:numero_simulazioni)
  { set.seed(i)   # fisso il seme per la riproducibilità
    true_theta <- runif(1,0.5,4)
    TRUE_THETA[i] <- true_theta
    y_oss <- simulatore(theta = true_theta, n_sim = 20)
    Y_OSS[i,] <- y_oss
    mle <- nlminb(start = 1, function(u) logL(u,dati = y_oss))
    out[i] <- mle$par
    IC_l <- uniroot(function(u) - logL(u, dati = y_oss) 
                    + mle$objective
                    + qchisq(0.95,1)/2, c(1e-07,mle$par))$root
    IC_u <- uniroot(function(u) - logL(u, dati = y_oss) 
                    +  mle$objective
                    + qchisq(0.95,1)/2, c(mle$par,10))$root
    intervalli[i,] <- c(IC_l,IC_u) 
  }
  list(TRUE_THETA = TRUE_THETA,Y_OSS = Y_OSS, 
       result =cbind(intervalli, out) )
}
intervalli_W <- studio_simulazione()
Y_OSS <- intervalli_W$Y_OSS

# write.csv(as.data.frame(Y_OSS), 'campioni_simulati.csv')
# grafici in pyton 
TRUE_THETA <-intervalli_W$TRUE_THETA
set.seed(299)  ; theta_star <- rgamma(60,3,0.5)
set.seed(9389) ; x_star     <- unlist(lapply(theta_star, 
                                             function(x) simulatore(x)))
set.seed(982)  ; initial_params <- runif(40)


IC <- intervalli_W$result[,1:2]


library(modeest)
# Funzione che trasorma l'input theta nei
# parametri della mistura p(x|theta) (surrogato)
# a partire da valori causali dei pesi W (sono 50 in tutto)

# costruzione della rete 
SYNT_logL <- function(input, param, target)
{
  w1      <- matrix(param[1:5],
                    nrow = dim_input, ncol = hidden_size)
  w_mu    <- matrix(param[6:20], nrow = hidden_size, ncol = k)
  w_sigma <- matrix(param[21:25],
                    nrow = hidden_size, ncol = 1)
  w_pi    <- matrix(param[26:40],
                    nrow = hidden_size, ncol = k)
  h       <- input %*% w1
  media   <- h %*% w_mu
  sigma   <- exp(h %*% w_sigma)
  pi      <- softmax(h %*% w_pi)
  # verosimiglianza della mistura
  sum(log(dnorm(target, media[1], sigma)*pi[1] +
            dnorm(target, media[2], sigma)*pi[2]+
            dnorm(target, media[3], sigma)*pi[3]))
}

# Funzione che fa MCMC sulla u(\theta) 
# aggiornata ottenuta ad ogni passo di SNL
MCMC <- function(n_sim = 3000, burnin = 1000, 
                 eps = 2.5, th_0 = runif(1,1,5), param, target)
{
  acc <- 0
  th_values <- numeric(n_sim)
  th <- th_0
  th_values[1] <- th_0
  for (j in 2:n_sim)
  {
    th_star <- th + runif(1,-eps,eps)
    if (th_star <= 0) {alpha <- 0}
    else
    {
      alpha   <- exp(SYNT_logL(param = param, 
                               input = th_star, target = target))/
        (exp(SYNT_logL(param = param, 
                       input = th, target = target)))
    }
    # print(alpha)
    
    if (runif(1) < alpha)
    {
      th <- th_star
      acc <- acc + 1
    }
    th_values[j] <- th
  }
  list(th_dist = th_values[500:n_sim], acc = acc/n_sim)
}

studio_simulazione_surrogato <- function(n_simulazioni)
{
  mode <- numeric(n_simulazioni)
  for (i in 1:n_simulazioni)
  {
    obs <- Y_OSS[i,]
    set.seed(299)  ; theta_star <- rgamma(50,3,0.5)
    set.seed(9389) ; x_star     <- unlist(lapply(theta_star, 
                                                 function(x) simulatore(x)))
    set.seed(982)  ; initial_params <- runif(40)
    opt <- optim(initial_params, function(x) -SYNT_logL(x, 
                                                        input = theta_star, target = x_star),
                 method = 'BFGS', control = list(trace = 0,
                                                 maxit = 2000))
    
    par_update <- opt$par
    theta_dist <- MCMC(n_sim = 3000, burnin = 1500, eps = 1, 
                       th_0 = 5.60, param = par_update, target = obs)
    
    for (r in 1:8)
    {
      new_theta <- sample(theta_dist$th_dist, 20, replace = F)
      theta_star <- c(theta_star,new_theta)
      x_star     <- c(x_star,unlist(lapply(new_theta, 
                                           function(x) simulatore(x))))
      
      
      # aggiorno i pesi della reta
      
      opt <- optim(initial_params, function(x) -SYNT_logL(x, 
                                                          input = theta_star, target = x_star),
                   method = 'BFGS', control = list(trace = 0, maxit = 2000))
      par_update <- opt$par
      # campiono dalla posteriori aggiornata
      
      theta_dist <- MCMC(n_sim = 2000, eps= 1, burnin = 500, 
                         param = par_update, target = obs)
      cat('=')
    }
    mode[i] <- mlv(theta_dist$th_dist, method = 'shorth')
    print(IC[i,1] <=  mode[i] & IC[i,2] >= mode[i] )
    print(c(IC[i,1],mode[i] ,IC[i,2]))
    theta_dist <- NULL
  }
  out <- mode
}
stime_surrogato <- studio_simulazione_surrogato(500)
stime_surrogato

studio_simulazione_risultati <- cbind(inizio_surrogato, IC)

sum(inizio_surrogato >= IC[,1] & inizio_surrogato <= IC[,2])/500

