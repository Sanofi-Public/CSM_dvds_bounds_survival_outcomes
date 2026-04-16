# This code allows for generating n.mc datasets

source("src/data_simulation_fun.R")

version <- 1

# Number of Monte-Carlo samples
n.mc <- 20

data.list <- list()

for (i in 1:n.mc) {

  #######################
  ### Data generation ###
  #######################

  p.X <- 2
  p.U <- 2
  n <- 1000

  gamma.data <- 3
  mixture.coeff <- 0.9
  shape <- 6

  simul.data.list <- genData(p.X=p.X, p.U=p.U, n=n, alpha=mixture.coeff,
                             gamma.data=gamma.data, shape=shape)

  data.list[[i]] <- simul.data.list

}

# Save list of data sets
file.name <- paste0("./data/simul/data_list_v", version, ".rds")
saveRDS(data.list, file=file.name)
