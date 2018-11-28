# ---------------------------------------------------------------------------------

# This R code is to propogate uncertainties in linear regression and isochron dating for the synthetic dataset for the paper "Monte Carlo sampling for error propagation in linear regression and applications in isochron geochronology"


# ---------------------------------------------------------------------------------

library(MASS)
library(ggplot2)
library(dplyr)
library(rjags)
library(openxlsx)

# need to change this to your own directory!!!
setwd("/Users/snoopy/Documents/Research/Projects/Other_People_Projects/Yang Li/Isochron/Code/R-code/Local code")

# Ludwig ------------------------------------------------------------------


LudwigRegression <- function(df_data, decay){
  
  x <- df_data$x
  dx <- df_data$dx * 0.5 # now always play with 1 sigma
  y <- df_data$y
  dy <- df_data$dy * 0.5
  rho <- df_data$rho
  num <- nrow(df_data) # total sample number
  
  # Note: a is intercept and b is slope
  
  
  N_itermax <- 1000 # max N of iterations, <50 times more than enough?
  tol <- 1e-15 # termination criterion
  
  # Initiallize weight parameter for ISOPLOT
  
  wx <- dx^-2 # x weight
  wy <- dy^-2 # y weight
  alpha <- (wx*wy)^0.5 # correlated weight of x,y
  
  # The initial guess of b can be anything, so just start as 1 is fine
  b0 <- rep(1, N_itermax+1)
  
  # Loop to get the best b
  
  for (i in 1:N_itermax){
    w <- wx*wy/(wx + b0[i]^2*wy - 2*b0[i]*rho*alpha) # weight of each sample
    w_sum <- sum(w)
    x_bar <- sum(w*x)/w_sum # all x's weigthed mean
    y_bar <- sum(w*y)/w_sum # all y's weigthed mean
    u <- x - x_bar # derivation of x[i] to x_bar
    v <- y - y_bar # derivation of y[i] to y_bar
    beta <- w*(u/wy + b0[i]*v/wx - (b0[i]*u + v)*rho/alpha) # residual
    b0[i+1] <- sum(w*beta*v)/sum(w*beta*u) # b's mathmatic expression
    
    if (abs((b0[i+1] - b0[i])/b0[i+1]) < tol) {
      b <- b0[i+1]
      break #yes, then stop here and the best b=b0(i+1)
    }
  }
  
  if (i == N_itermax){
    b <- b0[i+1]
  }
  
  # Use the best b to recalculate the parameter below
  
  w <- wx*wy/(wx + b^2*wy - 2*b*rho*alpha) # weight of each sample
  w_sum <- sum(w)
  x_bar <- sum(w*x)/w_sum # all x's weigthed mean
  y_bar <- sum(w*y)/w_sum # all y's weigthed mean
  u <- x - x_bar # derivation of x[i] to x_bar
  v <- y - y_bar # derivation of y[i] to y_bar
  beta <- w*(u/wy + b*v/wx - (b*u + v)*rho/alpha) # residual
  
  
  # Best y-intercept and uncertainties calculation
  
  a <- y_bar - b*x_bar # math expression of a=y-b*x
  x_adj <- x_bar + beta # adjusted x, or the expected value of x
  x_adj_bar <- sum(w*x_adj)/w_sum # weighted mean of adjusted x
  u <- x_adj - x_adj_bar # derivation of x[i] to x_bar
  sigma_b <- sum(w*u^2)^-0.5 # 1 sigma priori error of b
  sigma_a <- (w_sum^-1 + (x_adj_bar*sigma_b)^2)^0.5 # 1 sigma priori error of a
  
  # Parameters for Model 1
  
  delta_y <- y - (b*x + a) # residual on Y
  SSE <- delta_y^2*w # sum of squared error
  MSWD <- sum(SSE)/(num-2)
  prob <- 1 - pchisq(sum(SSE), num-2)
  b1 <- b
  a1 <- a
  
  # Expand or not
  
  if (prob >= 0.05){
    sigma_b1 <- sigma_b*2
    sigma_a1 <- sigma_a*2
  }else{
    sigma_b1 <- sigma_b*(MSWD^0.5)*qt(1 - (1 - 0.95)*0.5, num-2)
    sigma_a1 <- sigma_a*(MSWD^0.5)*qt(1 - (1 - 0.95)*0.5, num-2)
  }
  
  
  # Start model 3 solution --------------------------------------------------
  
  sigma_a3 <- sigma_a*(MSWD^0.5)
  b2 <- rep(1, N_itermax)
  tmp <- dy^2
  
  rho_y <- rho*(tmp/(tmp + sigma_a3^2))^0.5
  
  wt_x <- dx^-2
  wt_y <- (tmp + sigma_a3^2)^-1
  wt_mean <- (wt_x*wt_y)^0.5
  
  for (j in 1:N_itermax){
    wt <- wt_x*wt_y/(b2[j]^2*wt_y + wt_x - 2*b2[j]*rho_y*wt_mean)
    wt_sum <- sum(wt)
    x_bar <- sum(wt*x)/wt_sum
    y_bar <- sum(wt*y)/wt_sum
    wt_sqr <- wt^2
    u <- x - x_bar
    v <- y - y_bar
    uu <- u^2
    vv <- v^2
    uv <- u*v
    c <- (uu/wt_y - vv/wt_x)*wt_sqr
    d <- (uv/wt_x - rho_y*uu/wt_mean)*wt_sqr
    e <- (uv/wt_y - rho_y*vv/wt_mean)*wt_sqr
    Test <- sum(c)^2 + 4*sum(d)*sum(e)
    b2[j+1] <- (Test^0.5 - sum(c))/(2*sum(d))
    
    if (abs((b2[j+1] - b2[j])/b2[j+1]) < tol){
      b3 <- b2[j+1]
      break
    }
  }
  
  
  if (j == N_itermax){
    b3 <- b2[j+1]
  }
  
  # Use the best b3 to recalculate the parameter below
  
  wt <- wt_x*wt_y/(b3^2*wt_y + wt_x - 2*b3*rho_y*wt_mean)
  wt_sum <- sum(wt)
  x_bar <- sum(wt*x)/wt_sum
  y_bar <- sum(wt*y)/wt_sum
  wt_sqr <- wt^2
  
  
  # Best slop,intercept and uncertainties calculation
  
  a3 <- y_bar - b3*x_bar
  yf_resid <- a3 + b3*x - y
  wtd_resid <- yf_resid^2*wt
  sums <- sum(wtd_resid)
  MSWD3 <- sums/(num-2)
  Number <- wt*yf_resid*(rho_y*wt_mean - wt_y*b3)
  True_x <- x + Number/(wt_mean^2)
  Sum_xz <- sum(True_x*wt)
  Sum_x2z <- sum(True_x^2*wt)
  Denom <- Sum_x2z*wt_sum - Sum_xz^2
  VarSlApr <- wt_sum/Denom
  VarIntApr <- Sum_x2z/Denom
  CovInterSlope <- -Sum_xz/Denom
  RhoInterSlope <- CovInterSlope/(VarSlApr*VarIntApr)^0.5
  ErrSlApr <- VarSlApr^0.5
  ErrIntApr <- VarIntApr^0.5
  
  sigma_b3 <- ErrSlApr*(MSWD3^0.5)*qt(1 - (1 - 0.95)*0.5, num-2)
  sigma_a3 <- ErrIntApr*(MSWD3^0.5)*qt(1 - (1 - 0.95)*0.5, num-2)
  
  
  age1 <- log(b1+1)/decay/1e6
  age3 <- log(b3+1)/decay/1e6
  
  
  sigma_age1 <- sigma_b1/(decay*(b1+1)*1e6) # this is 2 sigma
  sigma_age3 <- sigma_b3/(decay*(b3+1)*1e6) # this is 2 sigma
  
  
  # Model 1 or model 3 ------------------------------------------------------
  
  if (prob >= 0.15){
    a <- a1
    sigma_a <- sigma_a1
    age <- age1
    sigma_age <- sigma_age1
    b <- b1
    sigma_b <- sigma_b1
    model_name <- 'Model 1'
    
  }else{
    a <- a3
    sigma_a <- sigma_a3  
    age <- age3
    sigma_age <- sigma_age3
    b <- b3
    sigma_b <- sigma_b3
    model_name <- 'Model 3'
  }
  
  return(list(lud_a=a, lud_a_2sigma=sigma_a, lud_age=age, lud_age_2sigma=sigma_age, model_name=model_name, n=num, MSWD=MSWD, prob=prob))
}


# MonteCarlo --------------------------------------------------------------


MonteCarloRegression <- function(df_data, decay){
  
  x <- df_data$x
  dx <- df_data$dx * 0.5 # now always play with 1 sigma
  y <- df_data$y
  dy <- df_data$dy * 0.5
  rho <- df_data$rho
  num <- nrow(df_data) # total sample number
  
  
  # Default input for simulation
  
  N_pick <- (100*2.58/1*1)^2 # maxmium times for simulation 
  
  N_pick <- as.integer(N_pick)
  
  N_newdata <- N_pick*10 # NO.of data generated in Gaussian distribution, better to be larger than N_pick so no overlap
  nnn <- 100 # For model uncertainties, smooth the final countor plot
  
  
  # Convert error ellipse to Gaussian distribution
  
  sigmaxy <- dx*dy*rho
  
  # construct newdata matrixs
  newdata <- matrix(NA, nrow = N_newdata, ncol = 2*num)
  
  for (m in 1:num){
    C <- matrix(c(dx[m]^2, sigmaxy[m], sigmaxy[m], dy[m]^2), 2, 2)
    # set.seed(121)
    R <- mvrnorm(n = N_newdata, c(x[m], y[m]), C)
    newdata[,2*m-1] <- R[,1]
    newdata[,2*m] <- R[,2]
  }
  
  # Analytical uncertainties
  
  b4sm <- numeric(N_pick)
  a4sm <- numeric(N_pick)
  r4sm <- numeric(N_pick)
  db4sm <- numeric(N_pick)
  da4sm <- numeric(N_pick)
  
  # Random pick data from above and fit line for N_pick times
  total <- num*N_pick
  
  randomNo <- matrix(sample.int(N_newdata, total, replace = T), nrow = num) # Get it out could speed the process a bit
  
  new_x <- double(num)
  new_y <- double(num)
  
  r <- 0.317311 # 68% confidence interval and this is fixed
  t <- qt(1-r/2, num-2)
  
  
  for (n in 1:N_pick){
    for (nn in 1:num){
      new_x[nn] <- newdata[randomNo[nn, n], 2*nn-1]
      new_y[nn] <- newdata[randomNo[nn, n], 2*nn]
    }
    
    Sx <- sum(new_x)
    Sy <- sum(new_y)
    Sxx <- sum(new_x^2)
    Syy <- sum(new_y^2)
    Sxy <- sum(new_x*new_y)
    
    beta <- (num*Sxy - Sx*Sy)/(num*Sxx - Sx^2)
    alpha <- Sy/num - beta*Sx/num
    se_sqr <- 1/num/(num-2)*(num*Syy - Sy^2 - beta^2*(num*Sxx - Sx^2)) # error variance (from wiki)
    sbeta_sqr <- num*se_sqr/(num*Sxx - Sx^2) # This is same to the ISLR book
    salpha_sqr <- sbeta_sqr*Sxx/num # This is same to the ISLR book
    sbeta <- sbeta_sqr^0.5*t
    salpha <- salpha_sqr^0.5*t
    
    
    b4sm[n] <- beta
    a4sm[n] <- alpha
    r4sm[n] <- -Sx/((num^0.5)*Sxx^0.5)
    db4sm[n] <- sbeta
    da4sm[n] <- salpha
    
  }
  
  # This is a trick to calculate the model uncertainty
  newa  = r4sm * da4sm
  
  spread  = matrix(rnorm(nnn*N_pick), nrow = N_pick, ncol = nnn)
  temp = matrix(rnorm(nnn*N_pick), nrow = N_pick, ncol = nnn)
  
  newa  = as.numeric(spread*newa + temp*sqrt(da4sm^2 - newa^2) + a4sm)
  newb  = as.numeric(b4sm + spread*db4sm)
  
  rm(spread, temp)
  
  loc <- ! (is.na(newa) | is.na(newb))
  
  newa <- newa[loc] # get rid of na value. coming from salpha_sqr negative
  newb <- newb[loc] # get rid of na value. coming from sbeta_sqr negative
  
  
  
  
  # Calculate age and uncertainty -------------------------------------------
  
  # analytical
  
  age4sm <- log(b4sm+1)/decay/1e6 # analytical
  
  
  age5 <- mean(age4sm)
  sigma_age5 <- 2*sd(age4sm)
  
  
  # analytical and model
  
  newage <- log(newb+1)/decay/1e6 # analytical and model
  
  

  age55 <- mean(newage)
  sigma_age55 <- 2*sd(newage)
  
  # intercept
  a5 <- mean(a4sm) # analytical
  sigma_a5 <- 2*sd(a4sm)
  
  a55 <- mean(newa) # analytical and model
  sigma_a55 <- 2*sd(newa)
  
  # cor
  
  r5=cor(age4sm, a4sm)
  
  r55 <- cor(newage, newa)
  
  
  return(list(ana_age=age5, ana_age_2sigma=sigma_age5, ana_a=a5, ana_a_2sigma=sigma_a5, ana_corr=r5, ana_model_age=age55, ana_model_age_2sigma=sigma_age55, ana_model_a=a55, ana_model_a_2sigma=sigma_a55, ana_model_corr=r55))
  
}



# Create the dataset (no geological scatter)


total_output <- data.frame(
  lud_age = numeric(0), 
  lud_age_2sigma = numeric(0), 
  lud_a= numeric(0),
  lud_a_2sigma = numeric(0),
  model_name = character(0),
  n = integer(0),
  MSWD = numeric(0),
  prob = numeric(0),
  
  ana_age = numeric(0),
  ana_age_2sigma = numeric(0),
  ana_a= numeric(0),
  ana_a_2sigma = numeric(0),
  ana_corr = numeric(0),
  
  ana_model_age = numeric(0),
  ana_model_age_2sigma = numeric(0),
  ana_model_a= numeric(0),
  ana_model_a_2sigma = numeric(0),
  ana_model_corr = numeric(0),
  
  bayes_ana_age = numeric(0),
  bayes_ana_age_2sigma = numeric(0),
  bayes_ana_a= numeric(0),
  bayes_ana_a_2sigma = numeric(0),
  bayes_ana_corr = numeric(0),
  
  bayes_ana_model_age = numeric(0),
  bayes_ana_model_age_2sigma = numeric(0),
  bayes_ana_model_a= numeric(0),
  bayes_ana_model_a_2sigma = numeric(0),
  bayes_ana_model_corr = numeric(0)
  
)


# Total number of resampling

num_sampling = 10

size_low <- 5
size_high <- 30


time_low <- 100e6
time_high <- 4500e6


initial_low <- 0.2 # note that changes from 0.3 to 0.05
initial_high <- 1.2


x_low <- 100
x_high <- 1000

decay <- 1.666e-11

dx_low <- 0.002
dx_high <- 0.01

dy_low <- 0.002
dy_high <- 0.01

rho_low <- 0.4
rho_high <- 0.999

N <- sample(size_low:size_high, size = num_sampling, replace = T)

# Begin data generation

count <- 0

start <- Sys.time()

for (size_sample in N){
  
  # create the dataset
  x_low_high <- runif(n = 2, min = 100, max = 1000)
  
  x_low <- min(x_low_high)
  x_high <- max(x_low_high)
  
  x <- runif(n = size_sample, min = x_low, max = x_high)
  dx <- x * runif(n = size_sample, min = dx_low, max = dx_high)
  
  initial <- runif(1, min = initial_low, max = initial_high)
  
  time <- runif(1, min = time_low, max = time_high)
  
  y <- x * (exp(decay * time) - 1) + initial
  
  
  z <- y
  
  # z <- y * (1 + runif(size_sample, min = 0.08, max = 0.12) * sample(x = c(-1,1), size = size_sample, replace = T)) # old one
  
  # z <- y * (1 + runif(size_sample, min = 0.02, max = 0.06) * sample(x = c(-1,1), size = size_sample, replace = T)) # before jacky
  
  # z <- y * (1 + runif(size_sample, min = 0.00036, max = 0.0036) * sample(x = c(-1,1), size = size_sample, replace = T))
  
  # z <- y * (1 + runif(size_sample, min = 0.006, max = 0.013) * sample(x = c(-1,1), size = size_sample, replace = T)) # good for generating 0 to 60
  
  # z <- y * (1 + runif(size_sample, min = 0.0025, max = 0.007) * sample(x = c(-1,1), size = size_sample, replace = T)) # good for generating 2 to 5
  
  z <- y * (1 + runif(size_sample, min = 0.002, max = 0.005) * sample(x = c(-1,1), size = size_sample, replace = T)) # good for generating 2 to 5
  
  dz <- z * runif(n = size_sample, min = dy_low, max = dy_high)
  
  rho <- runif(n = size_sample, min = rho_low, max = rho_high)
  
  # run Ludwig model
  df_data <- data.frame(x = x, y = z, dx = dx, dy = dz, rho = rho)
  
  Ludwig_out <- LudwigRegression(df_data, decay)
  print('Ludwig finished')
  print(Sys.time() - start)
  
  
  Monte_out <- MonteCarloRegression(df_data, decay)
  print('Monte finished')
  print(Sys.time() - start)
  
  
  # Bayes_out <- BayesRegression(df_data, decay)
  # print('Bayes finished')
  # print(Sys.time() - start)
  
  
  total_output <- rbind(total_output, c(Ludwig_out, Monte_out))
  
  count <- count + 1
  
  cat('\n', count, 'iteration finished', '\n')
  print(Sys.time() - start)
  cat('\n', '---------------------------------')
  cat('\n')
  
}

print('program finished')
print(Sys.time() - start)


now <- Sys.time()

file_name <- paste0(format(now, "%Y%m%d_%H%M%S_"), "total_output.xlsx")

write.xlsx(total_output, file = file_name, colNames = TRUE)


