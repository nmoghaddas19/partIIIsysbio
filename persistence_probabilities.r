# measuring the effect of varying the landscape, independent of R_0,
# on the number of times the pathogen persisted in 100 simulations

N <- 500
r <- 1
th <- 0.25
B <- 0.09
w <- 0.25
g <- 0.25
a <- 0.25

# preparing matrix where each column is a particular event 
# and the entry in each row will become the rate of 
# the event occuring in a particular field
rates <- c(rep(1/(r-th),N),
           rep(1/(r-th), N),
           rep((1/th), N),
           rep(B,N),
           rep(1/th, N),
           rep(g, N))
events <- matrix(rates, byrow = FALSE, nrow = N)

# Gillespie's method

# initialising vectors
Z <- 15000 # the number of time steps to run for
J <- 100 # the number of times to run the whole simulation
t <- matrix(0, J, Z)
S <- matrix(0,Z,N)
I <- matrix(0,Z,N)
E <- matrix(0,Z,N)
C <- matrix(0,Z,N)
Su <- matrix(0,J,Z)
In <- matrix(0,J,Z)
Em <- matrix(0,J,Z)
Co <- matrix(0,J,Z)
event_tracker <- matrix(0,J,Z)
field_tracker <- matrix(0,J,Z)
test_tracker <- matrix(0,J,Z)
r_0 <- numeric(20)
eigenvalues <- numeric(20)
p_invasion <- numeric(20)

A <- matrix(0, N, 3)
A[,1] <- sample(1:250, N, TRUE)
A[,2] <- sample(1:250, N, TRUE)
A[,3] <- c(rep(1, N))

area <- matrix(0,N,N)
dist <- matrix(0,N,N)
for (j in seq(1, N, 1)) {
  for (i in seq(1, N, 1)) {
    area[j,i] <- A[j,3] * A[i,3]
    dist[j,i] <- sqrt((A[j,1] - A[i,1])^2 + (A[j,2] - A[i,2])^2)
  }
}
for (k in 1:20) {
  r_0[k] <- B*th^2/r * (g*(r-th)+1)/(g*(r-th)+1-a*w)
  eta <- seq(0.007,0.017,0.0005)[k]
  # creating matrix of dispersion kernels
  kernel <- exp(dist * eta * -1)
  M <- kernel * area
  
  for (i in seq(1, N, 1)) {
    M[i,i] <- 0
  }
  eigenvalues[k] <- eigen(M)$value[1] # storing each time the leading eignevalue
  
  for (j in 1:J) {
    time <- 0
    S[1,] <- c(rep(0,20), rep(1,N-20)) # initial conditions. First 10 fields infected, remainder susceptible. None empty or colonised.
    I[1,] <- c(rep(1,20),rep(0,N-20))
    E[1,] <- c(rep(0,N))
    C[1,] <- c(rep(0,N))
    
    for (i in seq(1, Z-1, 1)) {
      new_events <- matrix(0, N, 6)
      new_events[,1] <- E[i,] * events[,1] 
      new_events[,2] <- C[i,] * events[,2]
      new_events[,3] <- S[i,] * events[,3]
      
      betas <- c(M %*% I[i,]) 
      new_events[,4] <- betas * S[i,] * events[,4] 
      
      new_events[,5] <- I[i,] * events[,5]
      new_events[,6] <- C[i,] * events[,6]
      
      # time step
      delta_t <- rexp(1, sum(new_events))
      time <- time + delta_t
      t[j,i+1] <- time + delta_t
      
      # choose which event occured
      col_sum <- apply(new_events, 2, FUN=sum)
      col_sums <- col_sum / sum(new_events)
      event <- sample(c(1,2,3,4,5,6), size = 1, replace = TRUE, prob = col_sums)
      event_tracker[j,i+1] <- event
      
      # choose which field it occured in
      S_delta_t <- S[i,]
      I_delta_t <- I[i,]
      E_delta_t <- E[i,]
      C_delta_t <- C[i,]
      
      if (event == 1) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,1])
        E_delta_t[field] <- 0
        S_delta_t[field] <- 1
      }
      if (event == 2) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,2])
        C_delta_t[field] <- 0
        test <- runif(1)
        test_tracker[j,i+1] <- test
        if (test <= a) {
          I_delta_t[field] <- 1
        } 
        else { 
          S_delta_t[field] <- 1
        }
      }
      if (event == 3) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,3])
        S_delta_t[field] <- 0
        E_delta_t[field] <- 1
      }
      if (event == 4) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,4])
        S_delta_t[field] <- 0
        I_delta_t[field] <- 1
      }
      if (event == 5) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,5])
        I_delta_t[field] <- 0
        test <- runif(1)
        test_tracker[j,i+1] <- test
        if (test <= w) {
          C_delta_t[field] <- 1
        } 
        else { 
          E_delta_t[field] <- 1
        }
      }
      if (event == 6) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,6])
        C_delta_t[field] <- 0
        E_delta_t[field] <- 1
      }
      field_tracker[j,i+1] <- field
      
      S[i+1,] <- S_delta_t
      I[i+1,] <- I_delta_t
      E[i+1,] <- E_delta_t
      C[i+1,] <- C_delta_t
    }
    Su[j,] <- apply(S, 1, FUN = sum)
    In[j,] <- apply(I, 1, FUN = sum)
    Em[j,] <- apply(E, 1, FUN = sum)
    Co[j,] <- apply(C, 1, FUN = sum)
  }
  p_invasion[k] <- sum((In[,Z] + Co[,Z]) >= 1)/J # tracking whether or not the pathogen persisted at the end of the simulation
}


# measuring the effect of varying R_0, independent of the landscape,
# on the number of times the pathogen persisted in 100 simulations

N <- 500
r <- 1
th <- 0.25
B <- 0.10
w <- 0.25
g <- 0.25
a <- 0.25

# Gillespie's method

# initialising vectors
Z <- 15000 # the number of time steps to run for
J <- 100 # the number of times to run the whole simulation
t <- matrix(0, J, Z)
S <- matrix(0,Z,N)
I <- matrix(0,Z,N)
E <- matrix(0,Z,N)
C <- matrix(0,Z,N)
Su <- matrix(0,J,Z)
In <- matrix(0,J,Z)
Em <- matrix(0,J,Z)
Co <- matrix(0,J,Z)
event_tracker <- matrix(0,J,Z)
field_tracker <- matrix(0,J,Z)
test_tracker <- matrix(0,J,Z)
r_0 <- numeric(20)
p_invasion <- numeric(20)
for (k in 1:20) {
  B <- seq(0.04,0.20,0.008)[k]
  r_0[k] <- B*th^2/r * (g*(r-th)+1)/(g*(r-th)+1-a*w)
  
  # preparing matrix where each column is a particular event 
  # and the entry in each row will become the rate of 
  # the event occuring in a particular field
  rates <- c(rep(1/(r-th),N),
             rep(1/(r-th), N),
             rep((1/th), N),
             rep(B,N),
             rep(1/th, N),
             rep(g, N))
  events <- matrix(rates, byrow = FALSE, nrow = N)
  
  for (j in 1:J) {
    time <- 0
    S[1,] <- c(rep(0,20), rep(1,N-20)) # initial conditions. First 10 fields infected, remainder susceptible. None empty or colonised.
    I[1,] <- c(rep(1,20), rep(0,N-20))
    E[1,] <- c(rep(0,N))
    C[1,] <- c(rep(0,N))
    
    for (i in seq(1, Z-1, 1)) {
      new_events <- matrix(0, N, 6)
      new_events[,1] <- E[i,] * events[,1] 
      new_events[,2] <- C[i,] * events[,2]
      new_events[,3] <- S[i,] * events[,3]
      
      betas <- c(M %*% I[i,]) 
      new_events[,4] <- betas * S[i,] * events[,4] 
      
      new_events[,5] <- I[i,] * events[,5]
      new_events[,6] <- C[i,] * events[,6]
      
      # time step
      delta_t <- rexp(1, sum(new_events))
      time <- time + delta_t
      t[j,i+1] <- time + delta_t
      
      # choose which event occured
      col_sum <- apply(new_events, 2, FUN=sum)
      col_sums <- col_sum / sum(new_events)
      event <- sample(c(1,2,3,4,5,6), size = 1, replace = TRUE, prob = col_sums)
      event_tracker[j,i+1] <- event
      
      # choose which field it occured in
      S_delta_t <- S[i,]
      I_delta_t <- I[i,]
      E_delta_t <- E[i,]
      C_delta_t <- C[i,]
      
      if (event == 1) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,1])
        E_delta_t[field] <- 0
        S_delta_t[field] <- 1
      }
      if (event == 2) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,2])
        C_delta_t[field] <- 0
        test <- runif(1)
        test_tracker[j,i+1] <- test
        if (test <= a) {
          I_delta_t[field] <- 1
        } 
        else { 
          S_delta_t[field] <- 1
        }
      }
      if (event == 3) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,3])
        S_delta_t[field] <- 0
        E_delta_t[field] <- 1
      }
      if (event == 4) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,4])
        S_delta_t[field] <- 0
        I_delta_t[field] <- 1
      }
      if (event == 5) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,5])
        I_delta_t[field] <- 0
        test <- runif(1)
        test_tracker[j,i+1] <- test
        if (test <= w) {
          C_delta_t[field] <- 1
        } 
        else { 
          E_delta_t[field] <- 1
        }
      }
      if (event == 6) {
        field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,6])
        C_delta_t[field] <- 0
        E_delta_t[field] <- 1
      }
      field_tracker[j,i+1] <- field
      
      S[i+1,] <- S_delta_t
      I[i+1,] <- I_delta_t
      E[i+1,] <- E_delta_t
      C[i+1,] <- C_delta_t
    }
    Su[j,] <- apply(S, 1, FUN = sum)
    In[j,] <- apply(I, 1, FUN = sum)
    Em[j,] <- apply(E, 1, FUN = sum)
    Co[j,] <- apply(C, 1, FUN = sum)
  }
  p_invasion[k] <- sum((In[,Z] + Co[,Z]) >= 1)/J
}
