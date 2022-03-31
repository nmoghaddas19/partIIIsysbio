install.packages("deSolve")
require(deSolve)

# non spatial gillespie simulation ----

for (j in seq(1, 100, 1)) {
  
  # setting initial conditions and parameters ----
  S <- 490 
  I <- 50
  E <- 0
  C <- 0
  
  r <- 1
  th <- 0.25
  B <- 0.10
  w <- 0.25
  g <- 0.25
  a <- 0.25
  
  ts <- c(0)
  time <- 0
  susceptibles <- c(S)
  infected <- c(I)
  empty <- c(E)
  colonised <- c(C)
  
  for (i in seq(2, 10000, 1)) {
    
    # defining rates of transitions between states ----
    R1 <- 1/(r-th) * E
    R2 <- 1/(r-th) * C
    R3 <- (1/th) * S
    R4 <- B * S * I
    R5 <- 1/th * I
    R6 <- g * C
    
    # total rate ----
    T <- R1 + R2 + R3 + R4 + R5 + R6
    
    # choosing time interval
    time_jump <- rexp(1, T)
    time <- time + time_jump
    ts[i] <- time
    
    # choosing which event occurred ----
    test <- runif(1) 
    if (test <= R1/T) {
      S <- S + 1
      E <- E - 1
    }
    
    if (R1/T < test & test <= (R1 + R2)/T) {
      C <- C - 1
      if (runif(1) <= a) {
        I <- I + 1
      } else {
        S <- S + 1
      } 
    }
    
    if ((R1 + R2)/T < test & test <= (R1 + R2 + R3)/T) {
      S <- S - 1
      E <- E + 1
    }
    
    if ((R1 + R2 + R3)/T < test & test <= (R1 + R2 + R3 + R4)/T) {
      S <- S - 1
      I <- I + 1
    }
    
    if ((R1 + R2 + R3 + R4)/T < test & test <= (R1 + R2 + R3 + R4 + R5)/T) {
      I <- I - 1
      if (runif(1) <= w) {
        C <- C + 1
      } else {
        E <- E + 1
      }
    }
    
    if ((R1 + R2 + R3 + R4 + R5)/T < test & test <= (R1 + R2 + R3 + R4 + R5 + R6)/T) {
      E <- E + 1
      C <- C - 1
    }
    
    susceptibles[[i]] <- S
    infected[[i]] <- I
    colonised[[i]] <- C
    empty[[i]] <- E
  }
  
  if (j == 1) {
    Sus <- susceptibles
    Infe <- infected
    Col <- colonised
    Emp <- empty
    Tim <- ts
  }
  else {
    Sus <- rbind(Sus, susceptibles)
    Infe <- rbind(Infe, infected)
    Col <- rbind(Col, colonised)
    Emp <- rbind(Emp, empty)
    Tim <- rbind(Tim, ts)
  }
}
Sus_means <- apply(X=Sus, 2, FUN=mean)
Sus_sd <- apply(X=Sus, 2, FUN=sd)
# X is the matrix,
# 2 means we are applying it column-wise,
# FUN is the operation to be peformed
Inf_means <- apply(X=Infe, 2, FUN=mean)
Inf_sd <- apply(X=Infe, 2, FUN=sd)
Col_means <- apply(X=Col, 2, FUN=mean)
Col_sd <- apply(X=Col, 2, FUN=sd)
Emp_means <- apply(X=Emp, 2, FUN=mean)
Emp_sd <- apply(X=Emp, 2, FUN=sd)
Time_means <- apply(X=Tim, 2, FUN=mean)

# non spatial deterministic solution ----
test.RHS <- function(t, s, p) {
  dS <- 1 / (p[1] - p[2]) * s[3] + (1 - p[6]) / (p[1] - p[2]) * s[4]- (1 / p[2]) * s[1]- p[3] * s[1] * s[2]
  
  dI <- p[6] / (p[1] - p[2]) * s[4]+ p[3] * s[1] * s[2]- 1/p[2] * s[2]
  
  dE <- 1 / p[2] * s[1]+ (1 - p[5]) / p[2] * s[2]+ p[4] * s[4]- 1 / (p[1]-p[2]) * s[3]
  
  dC <- p[5] / p[2] * s[2]- p[4] * s[4]- 1 / (p[1] - p[2]) * s[4]
  
  return(list(c(dS, dI, dE, dC)))
} 

p <- c(1, 0.25, 0.1, 0.25, 0.25, 0.25)
ic <- c(490, 50, 0, 0)
times <- seq(0, 6.5, 0.001)

SI.soln <- ode(y=ic, times=times, func=test.RHS, parms=p)

# landscape matrix ----
N <- 540
C <- matrix(0, N, 3)
C[,1] <- sample(1:250, 540, TRUE)
C[,2] <- sample(1:250, 540, TRUE)
C[,3] <- c(rep(1, N))

eta <- 0.000001 # very small landscape parameter removes effect of space
area <- matrix(0,N,N)
dist <- matrix(0,N,N)
for (j in seq(1, N, 1)) {
  for (i in seq(1, N, 1)) {
    area[j,i] <- C[j,3] * C[i,3]
    dist[j,i] <- sqrt((C[j,1] - C[i,1])^2 + (C[j,2] - C[i,2])^2)
  }
}

# creating matrix of dispersion kernels
kernel <- exp(dist * eta * -1)
M <- kernel * area

for (i in seq(1, N, 1)) {
  M[i,i] <- 0
}


# spatially explicit deterministic solution 
para <- c(r=1, th=0.25, B=0.1, g=0.25, w=0.25, a=0.25)

S.init <- c(rep(0,50), rep(1,N-50))
I.init <- c(rep(1,50),rep(0,N-50))
E.init <- c(rep(0,N))
C.init <- c(rep(0,N))
ic <- c(S.init, I.init, E.init, C.init)

times <- seq(0, 6.5, 0.01)

SIECmodel <- function(t, y, p) {
  S <- y[1:N]
  I <- y[(N+1):(2*N)]
  E <- y[(2*N+1):(3*N)]
  C <- y[(3*N+1):(4*N)]
  
  output <- numeric(4*N)
  
  betas <- numeric(N)
  for (k in 1:N) {
    betas[k] <- sum(M[,k] * I) * p[[3]]
  }
  
  # maybe need to times by area and beta here
  
  for (k in 1:N) {
    output[k] <- 1/(p[1]-p[2]) * E[k] + (1-p[6])/(p[1]-p[2]) * C[k] - 1/p[2] * S[k] - betas[k] * S[k]
    output[k+N] <- p[6]/(p[1]-p[2]) * C[k] + betas[k] * S[k]-1/p[2] * I[k]
    output[k+2*N] <- 1/p[2] * S[k] + (1-p[5])/p[2] * I[k]+ p[4] * C[k] - 1/(p[1]-p[2]) * E[k]
    output[k+3*N] <- p[5]/p[2] * I[k] -p[4] * C[k] -1/(p[1]-p[2]) * C[k]
  }
  return(list(output))
}
SIEC.soln <- ode(y=ic, times=times, func=SIECmodel, parms=para)

SIEC.soln <- SIEC.soln[,-1]
S.out <- SIEC.soln[,c(1:N)]
I.out <- SIEC.soln[,c((N+1):(2*N))]
E.out <- SIEC.soln[,c((2*N+1):(3*N))]
C.out <- SIEC.soln[,c((3*N+1):(4*N))]
S.out <- apply(S.out, MARGIN=1, FUN=sum)
I.out <- apply(I.out, MARGIN=1, FUN=sum)
E.out <- apply(E.out, MARGIN=1, FUN=sum)
C.out <- apply(C.out, MARGIN=1, FUN=sum)

# spatially explicit stochastic simulation
rates <- c(rep(1/(r-th),N),
           rep(1/(r-th), N),
           rep((1/th), N),
           rep(B,N),
           rep(1/th, N),
           rep(g, N))
events <- matrix(rates, byrow = FALSE, nrow = N)

Z <- 10000
J <- 20
t <- matrix(0,J,Z)
S <- matrix(0,Z,N)
I <- matrix(0,Z,N)
E <- matrix(0,Z,N)
C <- matrix(0,Z,N)
Su <- matrix(0,J,Z)
In <- matrix(0,J,Z)
Em <- matrix(0,J,Z)
Co <- matrix(0,J,Z)

for (j in 1:J) {
  time <- 0
  S[1,] <- c(rep(0,50), rep(1,N-50))
  I[1,] <- c(rep(1,50),rep(0,N-50))
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
    
    # choose event 
    col_sum <- apply(new_events, 2, FUN=sum)
    col_sums <- col_sum / sum(new_events)
    event <- sample(c(1,2,3,4,5,6), size = 1, replace = TRUE, prob = col_sums)
    #event_tracker[j,i+1] <- event
    
    S_delta_t <- S[i,]
    I_delta_t <- I[i,]
    E_delta_t <- E[i,]
    C_delta_t <- C[i,]
    
    # choose field
    if (event == 1) {
      field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,1])
      E_delta_t[field] <- 0
      S_delta_t[field] <- 1
    }
    if (event == 2) {
      field <- sample(1:N, size = 1, replace = TRUE, prob = new_events[,2])
      C_delta_t[field] <- 0
      test <- runif(1)
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

S_n <- apply(Su, 2, FUN=mean)
S_sd <- apply(Su, 2, FUN=sd)
I_n <- apply(In, 2, FUN=mean)
I_sd <- apply(In, 2, FUN=sd)
E_n <- apply(Em, 2, FUN=mean)
E_sd <- apply(Em, 2, FUN=sd)
C_n <- apply(Co, 2, FUN=mean)
C_sd <- apply(Co, 2, FUN=sd)
t_n <- apply(t, 2, FUN=mean)

# plotting ----
plot(SI.soln[,"time"], SI.soln[,"2"], type = "l", xlab = "Time", ylab = "Number of fields", main = "Infected", ylim=c(0,400))
lines(Time_means, Inf_means, col = "red")
arrows(x0=Time_means, y0=Inf_means-Inf_sd*1.96/sqrt(100), x1=Time_means, y1=Inf_means+Inf_sd*1.96/sqrt(100), code=0, col='#FF000008', lwd=1)
lines(t_n, I_n, col = "blue")
arrows(x0=t_n, y0=I_n-I_sd*1.96/sqrt(20), x1=Time_means, y1=I_n+I_sd*1.96/sqrt(20), code=0, col="#0000FF08", lwd=1)
lines(times, I.out, col = "orange")
legend(x='topright', legend=c('Non spatial deterministic', 'Spatial deterministic', 'Non spatial simulation (100 runs)', 'Spatial simulation (20 runs)'), lty=1, col=c('black', 'orange','red', 'blue'), cex=1.1, lwd=3, bty='n')

plot(SI.soln[,"time"], SI.soln[,"1"], type = "l", xlab = "Time", ylab = "Number of fields", main = "Susceptible")
lines(Time_means, Sus_means, col = "red")
arrows(x0=Time_means, y0=Sus_means-Sus_sd*1.96/sqrt(100), x1=Time_means, y1=Sus_means+Sus_sd*1.96/sqrt(100), code=0, col='#FF000008', lwd=1)
lines(t_n, S_n, col = "blue")
arrows(x0=t_n, y0=S_n-S_sd*1.96/sqrt(20), x1=Time_means, y1=S_n+S_sd*1.96/sqrt(20), code=0, col="#0000FF08", lwd=1)
lines(times, S.out, col = "orange")
legend(x='topright', legend=c('Non spatial deterministic', 'Spatial deterministic', 'Non spatial simulation (100 runs)', 'Spatial simulation (20 runs)'), lty=1, col=c('black', 'orange','red', 'blue'), cex=1.1, lwd=3, bty='n')

plot(SI.soln[,"time"], SI.soln[,"3"], type = "l", xlab = "Time", ylab = "Number of fields", main = "Empty")
lines(Time_means, Emp_means, col = "red")
arrows(x0=Time_means, y0=Emp_means-Emp_sd*1.96/sqrt(100), x1=Time_means, y1=Emp_means+Emp_sd*1.96/sqrt(100), code=0, col='#FF000008', lwd=1)
lines(t_n, E_n, col = "blue")
arrows(x0=t_n, y0=E_n-E_sd*1.96/sqrt(20), x1=Time_means, y1=E_n+E_sd*1.96/sqrt(20), code=0, col="#0000FF08", lwd=1)
lines(times, E.out, col = "orange")
legend(x='bottomright', legend=c('Non spatial deterministic', 'Spatial deterministic', 'Non spatial simulation (100 runs)', 'Spatial simulation (20 runs)'), lty=1, col=c('black', 'orange','red', 'blue'), cex=1.1, lwd=3, bty='n')

plot(SI.soln[,"time"], SI.soln[,"4"], type = "l", xlab = "Time", ylab = "Number of fields", main = "Colonised")
lines(Time_means, Col_means, col = "red")
arrows(x0=Time_means, y0=Col_means-Col_sd*1.96/sqrt(100), x1=Time_means, y1=Col_means+Col_sd*1.96/sqrt(100), code=0, col='#FF000008', lwd=1)
lines(t_n, C_n, col = "blue")
arrows(x0=t_n, y0=C_n-C_sd*1.96/sqrt(20), x1=Time_means, y1=C_n+C_sd*1.96/sqrt(20), code=0, col="#0000FF08", lwd=1)
lines(times, C.out, col = "orange")
legend(x='bottomright', legend=c('Non spatial deterministic', 'Spatial deterministic', 'Non spatial simulation (100 runs)', 'Spatial simulation (20 runs)'), lty=1, col=c('black', 'orange','red', 'blue'), cex=1.1, lwd=3, bty='n')

