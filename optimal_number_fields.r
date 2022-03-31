install.packages("RSpectra")
require(RSpectra)

N <- 300
runs <- 5
areas <- matrix(0,runs,N)
sums <- matrix(0,runs,N)
rand_fields <- matrix(0,runs,N)
values_area <- matrix(0,runs,N)
values_sums <- matrix(0,runs,N)
values_rand <- matrix(0,runs,N)
removed_fields_matrix <- matrix(0,runs,N)
values_step_matrix <- matrix(0,runs,N)

for(i in 1:runs) {
  # creating a clustered landscape and its corresponding landscape matrix
  n_clusters <- 10
  cluster_centres_x <- runif(n_clusters, 0, 250)
  cluster_centres_y <- runif(n_clusters, 0, 250)
  fields_per_cluster <- round(N/n_clusters)
  N <- fields_per_cluster * n_clusters
  C <- matrix(0, N, 3)
  C[,1] <- rnorm(N, cluster_centres_x, 5) # x coordinates
  C[,2] <- rnorm(N, cluster_centres_y, 5) # y coordinates
  C[,3] <- rnorm(N, 2, 0.5) # areas
  #plot(C[,1], C[,2], pch = 0, cex = 0.7)
  
  eta <- 0.015
  area <- matrix(0,N,N)
  dist <- matrix(0,N,N)
  for (j in seq(1, N, 1)) {
    for (k in seq(1, N, 1)) {
      area[j,k] <- C[j,3] * C[k,3]
      dist[j,k] <- sqrt((C[j,1] - C[k,1])^2 + (C[j,2] - C[k,2])^2)
    }
  }
  
  # creating matrix of dispersion kernels
  kernel <- exp(dist * eta * -1)
  M <- kernel * area
  
  # replacing diagonals with 0
  for (l in seq(1, N, 1)) {
    M[l,l] <- 0
  }
  
  M_areas <- M
  M_sums <- M
  M_rand <- M
  M_step <- M
  
  # ordering fields by area and row sum for later
  areas[i,] <- sort(C[,3], decreasing = TRUE, index.return = TRUE)$ix[1:N] 
  row_sums <- apply(M, MARGIN=1, FUN=sum)
  sums[i,] <- sort(row_sums, decreasing = TRUE, index.return = TRUE)$ix[1:N]
  
  # random control
  rand_fields[i,] <- sample(1:N, replace = FALSE, N) 
  
  # choosing fields based on the greatest effect on the leading eigenvalue
  M_temp <- M
  temp_eigs <- numeric(N)
  removed_fields <- c()
  values_step <- numeric(N)
  for (n in 1:N) {
    for (o in 1:N) {
      temp_removed_fields <- c(removed_fields, o)
      temp_landscape <- M_temp
      temp_landscape[o,] <- numeric(N)
      temp_landscape[,o] <- numeric(N)
      temp_eigs[o] <- eigs_sym(temp_landscape, 1, which="LM")$values
    }
    removed_fields[n] <- sort(temp_eigs, index.return=TRUE)$ix[1]
    M_temp[removed_fields[n], ] <- numeric(N)
    M_temp[ ,removed_fields[n]] <- numeric(N)
    values_step[n] <- eigs_sym(M_temp, 1, which = "LM")$values
  }
  removed_fields_matrix[i,] <- removed_fields
  values_step_matrix[i,] <- values_step
  
  for(m in 1:N) {
    M_areas[areas[i,m],] <- numeric(N)
    M_areas[,areas[i,m]] <- numeric(N)
    values_area[i,m] <- eigs_sym(M_areas, 1, which = "LM")$values
    M_sums[sums[i,m],] <- numeric(N)
    M_sums[,sums[i,m]] <- numeric(N)
    values_sums[i,m] <- eigs_sym(M_sums, 1, which = "LM")$values
    M_rand[rand_fields[i,m],] <- numeric(N)
    M_rand[,rand_fields[i,m]] <- numeric(N)
    values_rand[i,m] <- eigs_sym(M_rand, 1, which = "LM")$values
  }
}
values_area_mean <- apply(values_area, MARGIN=2, FUN=mean)
values_area_sd <- apply(values_area, MARGIN=2, FUN=sd)
values_sums_mean <- apply(values_sums, MARGIN=2, FUN=mean)
values_sums_sd <- apply(values_sums, MARGIN=2, FUN=sd)
values_rand_mean <- apply(values_rand, MARGIN=2, FUN=mean)
values_rand_sd <- apply(values_rand, MARGIN=2, FUN=sd)

values_step_mean <- apply(values_step_matrix, MARGIN = 2, FUN=mean)
values_step_sd <- apply(values_step_matrix, MARGIN = 2, FUN=sd)

plot(values_step_mean, pch=20, col="red", frame.plot = FALSE, xlab="Number of fields", ylab="Leading eigenvalue")
points(values_rand_mean, pch=20)
points(values_sums_mean, pch=20, col="blue")
points(values_area_mean, pch=20, col="purple")
legend(x='topright', legend=c('Eigenvalue', 'Area', 'Row sum', 'Random'), lty=1, col=c('red', 'purple','blue', 'black',), cex=1.1, lwd=3, bty='n')


plot(values_rand_mean-values_step_mean, pch=20, col="red", frame.plot = FALSE, xlab="Number of fields protected", ylab="Leading eigenvalue (relative to random)")
points(values_rand_mean-values_rand_mean, pch=20)
points(values_rand_mean-values_sums_mean, pch=20, col="blue")
points(values_rand_mean-values_area_mean, pch=20, col="purple")
legend(x='topright', legend=c('Eigenvalue', 'Area', 'Row sum', 'Random'), lty=1, col=c('red', 'purple','blue', 'black',), cex=1.1, lwd=3, bty='n')

plot((values_rand_mean-values_step_mean)/values_rand_mean, pch=20, col="red", frame.plot = FALSE, xlab="Number of fields", ylab="Leading eigenvalue (fractional change relative to random)")
points((values_rand_mean-values_rand_mean)/values_rand_mean, pch=20)
points((values_rand_mean-values_sums_mean)/values_rand_mean, pch=20, col="blue")
points((values_rand_mean-values_area_mean)/values_rand_mean, pch=20, col="purple")
legend(x='topleft', legend=c('Eigenvalue', 'Area', 'Row sum', 'Random'), lty=1, col=c('red', 'purple','blue', 'black'), cex=1.1, lwd=3, bty='n')