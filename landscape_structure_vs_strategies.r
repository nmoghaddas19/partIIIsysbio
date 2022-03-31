# this creates 300 small landscapes of 16 clusters and varies the average 
# between cluster distance to show the effect of protecting 10 fields according 
# to different strategies: 1.simulated annealing optmisation algorithm 2. step-wise
# selection of the field with greatest effect on the eigenvalue 3. by largest
# area 4. by largest rowsum in the landscape matrix 5. random selection 6. exhaustive
# search of all possible combinations of fields

N <- 64
clusters <- 16
step_matrix <- matrix(0, 300, 20)
rand_matrix <- matrix(0, 300, 20)
area_matrix <- matrix(0, 300, 20)
sums_matrix <- matrix(0, 300, 20)
annealing_matrix <- matrix(0, 300, 20)
betclusters_matrix <- matrix(0, 300, 20)
withbet_matrix <- matrix(0, 300, 20)
for(m in 1:300) {
  values_exhaustive_min <- numeric(20)
  values_step_min <- numeric(20)
  values_annealing_min <- numeric(20)
  values_area_min <- numeric(20)
  values_sums_min <- numeric(20)
  values_rand_min <- numeric(20)
  within_clusters <- numeric(20)
  #st_devs <- c(1,2,3,5,10)
  scale <- c(0.005, 0.01, 0.05,0.1,0.2,0.4,0.6,0.8,1,1.1,1.2,1.3,1.4,1.6,1.8,2.0,2.5,3,3.5,4)
  with_bet_clusters <- numeric(20)
  bet_clusters <- numeric(20)
  centre_x <- runif(clusters,0,350)
  centre_y <- runif(clusters,0,350)
  
  for(i in 1:20) {
    n_clusters <- clusters
    cluster_centres_x <- centre_x * scale[i]
    cluster_centres_y <- centre_y * scale[i]
    fields_per_cluster <- round(N/n_clusters)
    N <- fields_per_cluster * n_clusters
    C <- matrix(0, N, 3)
    C[,1] <- rnorm(N, cluster_centres_x, 2)
    C[,2] <- rnorm(N, cluster_centres_y, 2)
    C[,3] <- rnorm(N, 2, 0.5)
    #plot(C[,1], C[,2], pch = 0, cex = 0.7)
    
    tf_vector <- rep(FALSE, 16)
    this_mean <- 0
    for(h in 1:16) {
      x_this <- tf_vector
      x_this[h] <- T
      this_mean <- this_mean + mean(dist(C[x_this,1:2]))
    }
    within_clusters[i] <- this_mean/16
    
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
    
    for (l in seq(1, N, 1)) {
      M[l,l] <- 0
    }
    
    nF <- 10 # number of fields we can protect
    
    # true solution: testing every single possibility - only feasible for small N
    
    #fields <- combn(1:N, nF)
    #values_exhaustive <- numeric(ncol(fields))
    #for(j in 1:ncol(fields)) {
    #  M_exhaustive <- M
    #  this_fields <- fields[,j]
    #  M_exhaustive[this_fields,] <- numeric(N)
    #  M_exhaustive[,this_fields] <- numeric(N)
    #  values_exhaustive[j] <- eigs_sym(M_exhaustive, 1, which="LM")$values
    #}
    #values_exhaustive_min[i] <- min(values_exhaustive)
    
    # step-wise selestion of the field with the greatest effect on the eigenvalue
    M_temp <- M
    temp_eigs <- numeric(5)
    removed_fields <- c()
    values_step <- numeric(5)
    for (n in 1:10) {
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
    values_step_min[i] <- min(values_step)
  
    f1 <- sample(1:N, replace = FALSE, nF)
    bad_jumps <- numeric(2500)
    eigen_i <- numeric(2500)
    Es <- numeric(2500)
    Temp <- seq(0.005, 0.000001, -0.000002)
    
    # simulated annealing
    for (k in 1:2500) {
      reduced_M1 <- M
      reduced_M1[f1,] <- numeric(N)
      reduced_M1[,f1] <- numeric(N)
      eigen_M1 <- eigs_sym(reduced_M1, 1, which = "LM")$values
      f2 <- f1
      new_field <- sample((1:N)[-f1], 1)
      switched_field <- sample(1:nF,1)
      f2[switched_field] <- new_field
      reduced_M2 <- M
      reduced_M2[f2,] <- numeric(N)
      reduced_M2[,f2] <- numeric(N)
      eigen_M2 <- eigs_sym(reduced_M2, 1, which = "LM")$values
      E <- eigen_M2 - eigen_M1
      Es[i]<- E
      if (E < 0) {
        f1 <- f2
        eigen_i[k] <- eigen_M2
      } else if (exp(-E/Temp[i]) > runif(1)) {
        f1 <- f2
        eigen_i[k] <- eigen_M2
        bad_jumps[k] <- 1
      } else {
        eigen_i[k] <- eigen_M1
      }
    }
    values_annealing_min[i] <- min(eigen_i)
    
    M_areas <- M
    M_sums <- M
    M_rand1 <- M
    M_rand2 <- M
    M_rand3 <- M
    M_rand4 <- M
    M_rand5 <- M
    M_step <- M
    
    values_area <- numeric(nF)
    values_sums <- numeric(nF)
    values_rand1 <- numeric(nF)
    values_rand2 <- numeric(nF)
    values_rand3 <- numeric(nF)
    values_rand4 <- numeric(nF)
    values_rand5 <- numeric(nF)
    
    areas <- sort(C[,3], decreasing = TRUE, index.return = TRUE)$ix[1:nF] 
    row_sums <- apply(M, MARGIN=1, FUN=sum)
    sums <- sort(row_sums, decreasing = TRUE, index.return = TRUE)$ix[1:nF]
    rand_fields1 <- sample(1:N, replace = FALSE, nF) 
    rand_fields2 <- sample(1:N, replace = FALSE, nF) 
    rand_fields3 <- sample(1:N, replace = FALSE, nF) 
    rand_fields4 <- sample(1:N, replace = FALSE, nF) 
    rand_fields5 <- sample(1:N, replace = FALSE, nF) 
    
    M_areas[areas,] <- numeric(N)
    M_areas[,areas] <- numeric(N)
    values_area <- eigs_sym(M_areas, 1, which = "LM")$values
    M_sums[sums,] <- numeric(N)
    M_sums[,sums] <- numeric(N)
    values_sums <- eigs_sym(M_sums, 1, which = "LM")$values
    M_rand1[rand_fields1,] <- numeric(N)
    M_rand1[,rand_fields1] <- numeric(N)
    values_rand1 <- eigs_sym(M_rand1, 1, which = "LM")$values
    M_rand2[rand_fields2,] <- numeric(N)
    M_rand2[,rand_fields2] <- numeric(N)
    values_rand2 <- eigs_sym(M_rand2, 1, which = "LM")$values
    M_rand3[rand_fields3,] <- numeric(N)
    M_rand3[,rand_fields3] <- numeric(N)
    values_rand3 <- eigs_sym(M_rand3, 1, which = "LM")$values
    M_rand4[rand_fields4,] <- numeric(N)
    M_rand4[,rand_fields4] <- numeric(N)
    values_rand4 <- eigs_sym(M_rand4, 1, which = "LM")$values
    M_rand5[rand_fields5,] <- numeric(N)
    M_rand5[,rand_fields5] <- numeric(N)
    values_rand5 <- eigs_sym(M_rand5, 1, which = "LM")$values
    
    values_area_min[i] <- min(values_area)
    values_sums_min[i] <- min(values_sums)
    values_rand_min[i] <- (values_rand1+values_rand2+values_rand3+values_rand4+values_rand5)/5
    bet_clusters[i] <- mean(dist(cbind(cluster_centres_x,cluster_centres_y)))
    with_bet_clusters[i] <- within_clusters[i] / bet_clusters[i]
  }
  
  step_matrix[m,] <- values_step_min
  rand_matrix[m,] <- values_rand_min
  area_matrix[m,] <- values_area_min
  sums_matrix[m,] <- values_sums_min
  annealing_matrix[m,] <- values_annealing_min
  betclusters_matrix[m,] <- bet_clusters
  withbet_matrix[m,] <- with_bet_clusters
}

values_step_min <- apply(step_matrix, MARGIN=2,FUN=mean)
values_step_sd <- apply(step_matrix, MARGIN=2,FUN=sd)
values_rand_min <- apply(rand_matrix, MARGIN=2,FUN=mean)
values_rand_sd <- apply(rand_matrix, MARGIN=2,FUN=sd)
values_area_min <- apply(area_matrix, MARGIN=2,FUN=mean)
values_area_sd <- apply(area_matrix, MARGIN=2,FUN=sd)
values_sums_min <- apply(sums_matrix, MARGIN=2,FUN=mean)
values_sums_sd <- apply(sums_matrix, MARGIN=2,FUN=sd)
values_annealing_min <- apply(annealing_matrix, MARGIN=2,FUN=mean)
values_annealing_sd <- apply(annealing_matrix, MARGIN=2,FUN=sd)
bet_clusters <- apply(betclusters_matrix, MARGIN=2,FUN=mean)
with_bet_clusters <- apply(withbet_matrix, MARGIN=2,FUN=mean)
 
# this creates 200 small landscapes of 64 fields and varies the number of 
# clusters to show the effect of protecting 10 fields according 
# to different strategies: 1.simulated annealing optmisation algorithm 2. step-wise
# selection of the field with greatest effect on the eigenvalue 3. by largest
# area 4. by largest rowsum in the landscape matrix 5. random selection 6. exhaustive
# search of all possible combinations of fields

N <- 64
clusters <- c(1,2,4,8,16,32,64)
step_matrix <- matrix(0, 200, 7)
rand_matrix <- matrix(0, 200, 7)
area_matrix <- matrix(0, 200, 7)
sums_matrix <- matrix(0, 200, 7)
annealing_matrix <- matrix(0, 200, 7)
betclusters_matrix <- matrix(0, 200, 7)
withbet_matrix <- matrix(0, 200, 7)
for(m in 1:200) {
  values_exhaustive_min <- numeric(7)
  values_step_min <- numeric(7)
  values_annealing_min <- numeric(7)
  values_area_min <- numeric(7)
  values_sums_min <- numeric(7)
  values_rand_min <- numeric(7)
  within_clusters <- numeric(7)
  #st_devs <- c(1,2,3,5,10)
  #scale <- c(0.005, 0.01, 0.05,0.1,0.2,0.4,0.6, 0.8, 1,1.1, 1.2,1.3, 1.4,1.6,1.8,2.0,2.5,3,3.5,4)
  with_bet_clusters <- numeric(7)
  bet_clusters <- numeric(7)
  
  for(i in 1:7) {
    n_clusters <- clusters[i]
    centre_x <- runif(n_clusters,0,450)
    centre_y <- runif(n_clusters,0,450)
    cluster_centres_x <- centre_x 
    cluster_centres_y <- centre_y 
    
    
    C <- matrix(0, N, 3)
    C[,1] <- rnorm(N, cluster_centres_x, 5)
    C[,2] <- rnorm(N, cluster_centres_y, 5)
    C[,3] <- runif(N, 0.1, 2)
    #plot(C[,1], C[,2], pch = 0, cex = 0.7)
    
    tf_ting <- rep(FALSE, n_clusters)
    this_mean <- 0
    for(h in 1:n_clusters) {
      x_this <- tf_ting
      x_this[h] <- T
      this_mean <- this_mean + mean(dist(C[x_this,1:2]))
    }
    within_clusters[i] <- this_mean/n_clusters
    
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
    
    for (l in seq(1, N, 1)) {
      M[l,l] <- 0
    }
    
    #fields <- combn(1:N, nF)
    #values_exhaustive <- numeric(ncol(fields))
    #for(j in 1:ncol(fields)) {
    #  M_exhaustive <- M
    #  this_fields <- fields[,j]
    #  M_exhaustive[this_fields,] <- numeric(N)
    #  M_exhaustive[,this_fields] <- numeric(N)
    #  values_exhaustive[j] <- eigs_sym(M_exhaustive, 1, which="LM")$values
    #}
    #values_exhaustive_min[i] <- min(values_exhaustive)
    
    nF <- 10
    M_temp <- M
    temp_eigs <- numeric(nF)
    removed_fields <- c()
    values_step <- numeric(nF)
    for (n in 1:nF) {
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
    values_step_min[i] <- min(values_step)
    
    
    f1 <- sample(1:N, replace = FALSE, nF)
    bad_jumps <- numeric(2500)
    eigen_i <- numeric(2500)
    Es <- numeric(2500)
    Temp <- seq(0.005, 0.000001, -0.000002)
    
    for (k in 1:2500) {
      reduced_M1 <- M
      reduced_M1[f1,] <- numeric(N)
      reduced_M1[,f1] <- numeric(N)
      eigen_M1 <- eigs_sym(reduced_M1, 1, which = "LM")$values
      f2 <- f1
      new_field <- sample((1:N)[-f1], 1)
      switched_field <- sample(1:nF,1)
      f2[switched_field] <- new_field
      reduced_M2 <- M
      reduced_M2[f2,] <- numeric(N)
      reduced_M2[,f2] <- numeric(N)
      eigen_M2 <- eigs_sym(reduced_M2, 1, which = "LM")$values
      E <- eigen_M2 - eigen_M1
      Es[i]<- E
      if (E < 0) {
        f1 <- f2
        eigen_i[k] <- eigen_M2
      } else if (exp(-E/Temp[i]) > runif(1)) {
        f1 <- f2
        eigen_i[k] <- eigen_M2
        bad_jumps[k] <- 1
      } else {
        eigen_i[k] <- eigen_M1
      }
    }
    values_annealing_min[i] <- min(eigen_i)
    
    M_areas <- M
    M_sums <- M
    M_rand1 <- M
    M_rand2 <- M
    M_rand3 <- M
    M_rand4 <- M
    M_rand5 <- M
    M_step <- M
    
    values_area <- numeric(nF)
    values_sums <- numeric(nF)
    values_rand1 <- numeric(nF)
    values_rand2 <- numeric(nF)
    values_rand3 <- numeric(nF)
    values_rand4 <- numeric(nF)
    values_rand5 <- numeric(nF)
    
    areas <- sort(C[,3], decreasing = TRUE, index.return = TRUE)$ix[1:nF] 
    row_sums <- apply(M, MARGIN=1, FUN=sum)
    sums <- sort(row_sums, decreasing = TRUE, index.return = TRUE)$ix[1:nF]
    rand_fields1 <- sample(1:N, replace = FALSE, nF) 
    rand_fields2 <- sample(1:N, replace = FALSE, nF) 
    rand_fields3 <- sample(1:N, replace = FALSE, nF) 
    rand_fields4 <- sample(1:N, replace = FALSE, nF) 
    rand_fields5 <- sample(1:N, replace = FALSE, nF) 
    
    M_areas[areas,] <- numeric(N)
    M_areas[,areas] <- numeric(N)
    values_area <- eigs_sym(M_areas, 1, which = "LM")$values
    M_sums[sums,] <- numeric(N)
    M_sums[,sums] <- numeric(N)
    values_sums <- eigs_sym(M_sums, 1, which = "LM")$values
    M_rand1[rand_fields1,] <- numeric(N)
    M_rand1[,rand_fields1] <- numeric(N)
    values_rand1 <- eigs_sym(M_rand1, 1, which = "LM")$values
    M_rand2[rand_fields2,] <- numeric(N)
    M_rand2[,rand_fields2] <- numeric(N)
    values_rand2 <- eigs_sym(M_rand2, 1, which = "LM")$values
    M_rand3[rand_fields3,] <- numeric(N)
    M_rand3[,rand_fields3] <- numeric(N)
    values_rand3 <- eigs_sym(M_rand3, 1, which = "LM")$values
    M_rand4[rand_fields4,] <- numeric(N)
    M_rand4[,rand_fields4] <- numeric(N)
    values_rand4 <- eigs_sym(M_rand4, 1, which = "LM")$values
    M_rand5[rand_fields5,] <- numeric(N)
    M_rand5[,rand_fields5] <- numeric(N)
    values_rand5 <- eigs_sym(M_rand5, 1, which = "LM")$values
    
    values_area_min[i] <- min(values_area)
    values_sums_min[i] <- min(values_sums)
    values_rand_min[i] <- (values_rand1+values_rand2+values_rand3+values_rand4+values_rand5)/5
    bet_clusters[i] <- mean(dist(cbind(cluster_centres_x,cluster_centres_y)))
    with_bet_clusters[i] <- within_clusters[i] / bet_clusters[i]
  }
  
  step_matrix[m,] <- values_step_min
  rand_matrix[m,] <- values_rand_min
  area_matrix[m,] <- values_area_min
  sums_matrix[m,] <- values_sums_min
  annealing_matrix[m,] <- values_annealing_min
  betclusters_matrix[m,] <- bet_clusters
  withbet_matrix[m,] <- with_bet_clusters
}

values_step_min <- apply(step_matrix, MARGIN=2,FUN=mean)
values_step_sd <- apply(step_matrix, MARGIN=2,FUN=sd)
values_rand_min <- apply(rand_matrix, MARGIN=2,FUN=mean)
values_rand_sd <- apply(rand_matrix, MARGIN=2,FUN=sd)
values_area_min <- apply(area_matrix, MARGIN=2,FUN=mean)
values_area_sd <- apply(area_matrix, MARGIN=2,FUN=sd)
values_sums_min <- apply(sums_matrix, MARGIN=2,FUN=mean)
values_sums_sd <- apply(sums_matrix, MARGIN=2,FUN=sd)
values_annealing_min <- apply(annealing_matrix, MARGIN=2,FUN=mean)
values_annealing_sd <- apply(annealing_matrix, MARGIN=2,FUN=sd)
bet_clusters <- apply(betclusters_matrix, MARGIN=2,FUN=mean)
with_bet_clusters <- apply(withbet_matrix, MARGIN=2,FUN=mean)
