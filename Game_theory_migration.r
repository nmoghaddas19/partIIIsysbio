# defining survival and body size functions and their products for use later
FPm <- function(s) { 
  (0.48 + 0.21 / (1 +exp(15-s))) * exp(-0.5*(s-15)^2/1^2)/(1*sqrt(2*pi)) # survival rate of migrants times normal distribution of body weights
}
FPn <- function(s) {
  (0.10 + 0.8 / (1 +exp(15-s))) * exp(-0.5*(s-15)^2/1^2)/(1*sqrt(2*pi)) # survival rate of nonmigrants times normal distribution of body weights
}
P_survive_migrant <- function(s) {
  0.48 + 0.21 / (1 +exp(15-s))
}
P_survive_nonmigrant <- function(s) {
  0.10 + 0.8 / (1 +exp(15-s))
}
Fs <- function(s) {
  exp(-0.5*(s-15)^2/1^2)/(1*sqrt(2*pi))
}
s <- seq(0,30,0.05)
plot(s, P_survive_nonmigrant(s), pch=20, cex=0.7, ylim= c(0,1), ylab='Survival Probability', xlab='Body Size', frame.plot=FALSE)
points(s , P_survive_migrant(s), pch=20, cex=0.7, col='blue')
legend(x='topleft', legend=c('Migrant','Non migrant'), lty=c(1,1), col=c('blue','black'), cex=1.1, lwd=4, bty='n')

# plotting the effect of varying habitat quality on migration threshold and equilibrium population and proportions
best_territory <- seq(0.65,1.75,0.002)
s <- seq(0,30,length.out=length(best_territory))
Nm_soln <- numeric(length(s))
Nn_soln <- numeric(length(s))
Rn <- numeric(length(s))
Rm <- numeric(length(s))
fitness_m <- numeric(length(s))
fitness_n <- numeric(length(s))
threshold <- numeric(length(s))
Nm_therealdeal <- numeric(length(best_territory))
Nn_therealdeal <- numeric(length(best_territory))

for (j in 1:length(best_territory)) {
  a <- best_territory[j]
  b <- a/500000
  for (i in 1:length(s)) {
    S <- s[i]
    FPm_integral <- integrate(FPm, 0, S)
    FPn_integral <- integrate(FPn, S, 30)
    A <- ((FPm_integral$value + FPn_integral$value) * (2 + 2*a + b) - 2)/ (b *(FPm_integral$value + FPn_integral$value)^2)
    Nm_soln[i] <- A * FPm_integral$value
    Nn_soln[i] <- A * FPn_integral$value
    Rn[i] <- a - b*(Nn_soln[i]-1)/2
    Rm[i] <- a - b*Nn_soln[i] - b*(Nm_soln[i]-1)/2 
    fitness_m[i] <- P_survive_migrant(S) * (Rm[i] + 1)
    fitness_n[i] <- P_survive_nonmigrant(S) * (Rn[i] + 1)
  }
  threshold[j] <- approx(fitness_m-fitness_n, s, xout=0)$y
  therealFPm_integral <- integrate(FPm, 0, threshold[j])
  therealFPn_integral <- integrate(FPn, threshold[j], 30)
  A <- ((therealFPm_integral$value + therealFPn_integral$value) * (2 + 2*a + b) - 2)/ (b *(therealFPm_integral$value + therealFPn_integral$value)^2)
  Nm_therealdeal[j] <- A * therealFPm_integral$value
  Nn_therealdeal[j] <- A * therealFPn_integral$value
}
plot(best_territory*3, threshold, pch=20, cex = 0.7, cex.lab=1.2, xlab='Maximum Clutch Size', ylab='Migration Threshold (Body Size)', frame.plot=FALSE)

plot(best_territory*3, Nn_therealdeal/100000, ylim = c(min(Nm_therealdeal/100000),max(Nn_therealdeal/100000)), xlab='Maximum Clutch Size', ylab='Population (x100,000)', pch=20, cex = 0.7, cex.lab=1.2, frame.plot=FALSE) # xlim=c(2.5, max(best_territory*3)), ylim = c(0,max(Nn_therealdeal)),
points(best_territory*3, Nm_therealdeal/100000, col='blue', pch=20, cex = 0.7)

N <- Nm_therealdeal + Nn_therealdeal
prop_Nm <- Nm_therealdeal/N
prop_Nn <- Nn_therealdeal/N
plot(best_territory*3, prop_Nn, ylim=c(0,1), xlab='Maximum Clutch Size', ylab='Proportion', pch=20, cex = 0.7, cex.lab=1.2, frame.plot=FALSE, bg="NA")
points(best_territory*3, prop_Nm, col='blue', pch=20, cex = 0.7)
plot(best_territory*3, N/100000, xlab='Maximum Clutch Size', ylab='Total Population (x100,000)', pch=20, cex=0.7, cex.lab=1.2, frame.plot=FALSE)

# the effect of varying the habitat quality/clutch size on fitness difference between migrants and non-migrants
s <- seq(0,30,0.01)
Nm_soln2 <- numeric(length(s))
Nn_soln2 <- numeric(length(s))
Rn2 <- numeric(length(s))
Rm2 <- numeric(length(s))
fitness_m2 <- numeric(length(s))
fitness_n2 <- numeric(length(s))
threshold2 <- numeric(length(s))
#for (j in 1:length(best_territory)) {
for (j in 1:6) {
  #a <- c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8)[j]
  a <- c(0.2,0.7,1.2,1.7,2.3, 2.7)[j]
  b <- 0.0001
  for (i in 1:length(s)) {
    S <- s[i]
    FPm_integral <- integrate(FPm, min(s), S)
    FPn_integral <- integrate(FPn, S, max(s))
    A <- ((FPm_integral$value + FPn_integral$value) * (2 + 2*a + b) - 2)/ (b *(FPm_integral$value + FPn_integral$value)^2)
    Nm_soln2[i] <- A * FPm_integral$value
    Nn_soln2[i] <- A * FPn_integral$value
    Rn2[i] <- a - b*(Nn_soln2[i]-1)/2
    Rm2[i] <- a - b*Nn_soln2[i] - b*(Nm_soln2[i]-1)/2
    fitness_m2[i] <- P_survive_migrant(S) * (Rm2[i] + 1)
    fitness_n2[i] <- P_survive_nonmigrant(S) * (Rn2[i] + 1)
  }
  threshold <- approx(fitness_m2-fitness_n2, s, xout=0)$y
  approx(fitness_m2-fitness_n2, s, xout=0)$y
  if (j == 1) {
    plot(s, fitness_m2-fitness_n2, ylim=c(-2.3,1.2), pch=20, cex=0.7, col=j, frame.plot=FALSE, xlab='Body Size', ylab='Fitness Difference')
  } else {
    points(s, fitness_m2-fitness_n2, pch=20, cex=0.7, col=j)
  } 
}
abline(h=0,lty='dashed')
legend(x='topright', legend=c('a=0.6','a=2.1', 'a=3.6', 'a=4.1', 'a=6.6', 'a=7.1'), lty=c(1,1), col=c(1,2,3,4,5,6), cex=1.1, lwd=4, bty='n')

# the effect of habitat number on equilibrium populations and proportions 
best_territory <- seq(0.65,1.75,0.002)
s <- seq(0,30,length.out=length(best_territory))
Nm_soln <- numeric(length(s))
Nn_soln <- numeric(length(s))
Rn <- numeric(length(s))
Rm <- numeric(length(s))
fitness_m <- numeric(length(s))
fitness_n <- numeric(length(s))
threshold <- numeric(length(s))
Nm_therealdeal <- numeric(length(best_territory))
Nn_therealdeal <- numeric(length(best_territory))
total_habitats <- seq(100000,1000000, 1000)

for (j in 1:length(total_habitats)) {
  a <- 1.1
  b <- a/total_habitats[j]
  for (i in 1:length(s)) {
    S <- s[i]
    FPm_integral <- integrate(FPm, 0, S)
    FPn_integral <- integrate(FPn, S, 30)
    A <- ((FPm_integral$value + FPn_integral$value) * (2 + 2*a + b) - 2)/ (b *(FPm_integral$value + FPn_integral$value)^2)
    Nm_soln[i] <- A * FPm_integral$value
    Nn_soln[i] <- A * FPn_integral$value
    Rn[i] <- a - b*(Nn_soln[i]-1)/2
    Rm[i] <- a - b*Nn_soln[i] - b*(Nm_soln[i]-1)/2 
    fitness_m[i] <- P_survive_migrant(S) * (Rm[i] + 1)
    fitness_n[i] <- P_survive_nonmigrant(S) * (Rn[i] + 1)
  }
  threshold[j] <- approx(fitness_m-fitness_n, s, xout=0)$y
  therealFPm_integral <- integrate(FPm, 0, threshold[j])
  therealFPn_integral <- integrate(FPn, threshold[j], 30)
  A <- ((therealFPm_integral$value + therealFPn_integral$value) * (2 + 2*a + b) - 2)/ (b *(therealFPm_integral$value + therealFPn_integral$value)^2)
  Nm_therealdeal[j] <- A * therealFPm_integral$value
  Nn_therealdeal[j] <- A * therealFPn_integral$value
}
plot(total_habitats/100000, threshold, pch=20, cex = 0.7, xlab='Total Number of Habitats (x100000)', ylab='Migration Threshold (Body Size)', frame.plot=FALSE, cex.lab=1.2,)
plot(total_habitats/100000, Nn_therealdeal/100000, ylim = c(min(Nm_therealdeal/100000),max(Nn_therealdeal/100000)), xlab='Total Number of Habitats (x100000)', ylab='Population (x100,000)', pch=20, cex = 0.7, cex.lab=1.2, frame.plot=FALSE) # xlim=c(2.5, max(best_territory*3)), ylim = c(0,max(Nn_therealdeal)),
points(total_habitats/100000, Nm_therealdeal/100000, col='blue', pch=20, cex = 0.7)
legend(x=1, y=3.5, legend=c('Migrant','Non migrant'), lty=c(1,1), col=c('blue','black'), cex=1.2, lwd=4, bty='n')
N <- Nm_therealdeal + Nn_therealdeal
prop_Nm <- Nm_therealdeal/N
prop_Nn <- Nn_therealdeal/N
plot(total_habitats/100000, prop_Nn, ylim=c(0,1), xlab='Total Number of Habitats (x100000)', ylab='Proportion', pch=20, cex = 0.7, cex.lab=1.2, frame.plot=FALSE)
points(total_habitats/100000, prop_Nm, col='blue', pch=20, cex = 0.7)

plot(total_habitats/100000, N/100000, xlab='Total Number of Habitats (x100000)', ylab='Total Population (x100,000)', pch=20, cex=0.7, cex.lab=1.2, frame.plot=FALSE)
legend(x='topleft', legend=c('Migrant','Non migrant'), lty=c(1,1), col=c('blue','black'), cex=1.2, lwd=4, bty='n')

