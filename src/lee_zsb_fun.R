# Functions used by the method of Lee et al. (ZSB-like method)


#' Function to compute the bounds from Lee et al. (2024) for the survival function on Monte-Carlo samples. Works under parallel computation.
#'
#' @param data.list A list containing the Monte-Carlo samples
#' @param data.name Either "simul", "rhc", or "gbsg"
#' @param t.vec A vector of times for which to compute the bounds of the survival function
#' @param a If a=1, computes the bounds for the survival function among the treated and, if a=0, among the control
#' @param gamma Sensitivity parameter/Confounding strength
#' @param nb.folds Number of folds for cross-fitting
#' @param RMST If TRUE, computes the bounds for the RMST, else, computes the bounds for the survival functions
#'
#' @return A list of the bounds for each Monte-Carlo sample at each time point in t.vec
#' @export
#'
#' @examples
survLeeMc <- function(data.list, data.name="simul",
                      t.vec=seq(from=0.1, to=9, length.out=10), a=1,
                      gamma=1, nb.folds=5, RMST=FALSE) {
  
  # Number of Monte-Carlo samples
  n.mc <- length(data.list)
  t.vec.length <- length(t.vec)
  tau <- gamma / (1 + gamma)
  gamma.len <- length(gamma)
  
  exec.times.list <- list()
  bounds.list <- list()
  
  for (i in 1:n.mc) {
    
    message("Monte-Carlo dataset: ", i)
    
    # Retrieve dataset i
    simul.data.list <- data.list[[i]]
    p.X <- ncol(simul.data.list$X)
    n <- length(simul.data.list$A)
    
    if (data.name == "simul") {
      Delta <- 1 * (simul.data.list$T01 <= simul.data.list$C)
    } else if (data.name %in% c("rhc", "gbsg")) {
      Delta <- simul.data.list$Delta
    }
    
    sim.data.df <- data.frame(X=simul.data.list$X, U=simul.data.list$U,
                              A=simul.data.list$A, T.obs=simul.data.list$T.obs,
                              Delta=Delta)
    
    sim.data.df <- sim.data.df[order(sim.data.df$T.obs), ]
    
    # Quantile for administrative censoring
    q.adm.cens <- quantile(sim.data.df$T.obs, prob=0.95)
    sim.data.df[sim.data.df$T.obs > q.adm.cens, "T.obs"] <- q.adm.cens
    sim.data.df[sim.data.df$T.obs > q.adm.cens, "Delta"] <- 0
    
    if (a == 0) {
      # Invert treated and control group
      sim.data.df$A <- 1 - sim.data.df$A
    }
    
    if (nb.folds > 1) {
      # Index for each fold
      folds.ind <- split(sample(1:n, replace=FALSE),
                         cut(x=1:n, breaks=nb.folds, labels=FALSE))
    } else {
      folds.ind <- NULL
    }
    
    tk.unique.a.delta <- unique(sim.data.df$T.obs[sim.data.df$A == 1 & sim.data.df$Delta == 1])
    
    
    ###################################
    ### Propensity score estimation ###
    ###################################
    e.predict <- propScore(sim.data.df=sim.data.df, p.X=p.X,
                           nb.folds=nb.folds, folds.ind=folds.ind)
    weights <- (1 - e.predict) / e.predict
    
    # Weights for A=a
    weights.a.delta1 <- weights[sim.data.df$A == 1 & sim.data.df$Delta == 1]
    weights.a.delta0 <- weights[sim.data.df$A == 1 & sim.data.df$Delta == 0]
    
    T.obs.a.delta1 <- sim.data.df$T.obs[sim.data.df$A == 1 & sim.data.df$Delta == 1]
    T.obs.a.delta0 <- sim.data.df$T.obs[sim.data.df$A == 1 & sim.data.df$Delta == 0]
    
    nb.group.a.delta1 <- nrow(sim.data.df[sim.data.df$A == 1 & sim.data.df$Delta == 1, ])
    
    init <- rep(1, nb.group.a.delta1)
    
    
    # To compute the bounds under ignorability (simple Kaplan-Meier estimator)
    ign.KM.fun <- function(t) {
      
      num <- sum((1 + weights.a.delta1)[T.obs.a.delta1 == t])
      den <- sum((1 + weights.a.delta1)[T.obs.a.delta1 >= t]) + sum((1 + weights.a.delta0)[T.obs.a.delta0 >= t])
      
      c(num=num, den=den)
    }
    
    
    lb.a <- rep(NA, t.vec.length)
    ub.a <- rep(NA, t.vec.length)
    ign.a <- rep(NA, t.vec.length)
    
    start.time <- Sys.time()
    
    for (j in 1:t.vec.length) {
      
      message("Time index: ", j)
      
      t.j <- t.vec[j]
      
      message("Computing lower bound")
      # For the lower bound
      # fnscale=1 because it is a minimization problem
      lb.t.a.optim <- optim(par=init, fn=objectiveFun, time=t.j,
                            tk.unique.a.delta=tk.unique.a.delta,
                            sim.data.df=sim.data.df,
                            optim.prob="min", gamma=gamma,
                            weights.a.delta1=weights.a.delta1,
                            weights.a.delta0=weights.a.delta0,
                            T.obs.a.delta1=T.obs.a.delta1,
                            T.obs.a.delta0=T.obs.a.delta0,
                            RMST=RMST,
                            lower=1/gamma, upper=gamma, method="L-BFGS-B",
                            control=list(fnscale=1, maxit=100))
      
      message("Computing upper bound")
      # For the upper bound
      # fnscale=-1 because it is a maximization problem
      ub.t.a.optim <- optim(par=init, fn=objectiveFun, time=t.j,
                            tk.unique.a.delta=tk.unique.a.delta,
                            sim.data.df=sim.data.df,
                            optim.prob="max", gamma=gamma,
                            weights.a.delta1=weights.a.delta1,
                            weights.a.delta0=weights.a.delta0,
                            T.obs.a.delta1=T.obs.a.delta1,
                            T.obs.a.delta0=T.obs.a.delta0,
                            RMST=RMST,
                            lower=1/gamma, upper=gamma, method="L-BFGS-B",
                            control=list(fnscale=-1, maxit=100))
      
      message("Computing Kaplan-Meier estimator")
      tk.vec <- tk.unique.a.delta[tk.unique.a.delta <= t.j]
      
      if (sum(tk.unique.a.delta <= t.j) == 0) {
        # If the desired time point is lower than the minimum observed time
        if (RMST) {
          # The area is 0
          ign.a[j] <- 0
        } else {
          # The probability is one 1
          ign.a[j] <- 1
        }
        
      } else {
        
        num.den <- sapply(tk.vec, ign.KM.fun)
        num <- num.den[1, ]
        den <- num.den[2, ]
        
        if (RMST) {
          
          S.t <- cumprod(1 - num / den)
          # Sort by decreasing order
          S.t <- sort(S.t, decreasing=TRUE)
          # Midpoint rule (faster convergence rate than left or right rectangle method)
          S.t <- (c(1, S.t[-length(S.t)]) + S.t) / 2
          # Compute the time difference
          time.diff <- diff(c(0, tk.vec))
          # Compute the area by the rectangle method
          ign.a[j] <- sum(S.t * time.diff)
          
        } else {
          ign.a[j] <- prod(1 - num / den)
        }
        
      }
      
      lb.a[j] <- lb.t.a.optim$value
      ub.a[j] <- ub.t.a.optim$value
      
    }
    
    end.time <- Sys.time()
    elapsed <- end.time - start.time
    message(elapsed)
    exec.times.list[[i]] <- list(start.time=start.time, end.time=end.time, elapsed=elapsed)
    
    bounds.df <- data.frame(lb=lb.a, ub=ub.a, ign=ign.a, time=t.vec)
    
    if (RMST) {
      # All areas are null for RMST at time 0
      bounds.list[[i]] <- rbind(c(0, 0, 0, 0), bounds.df)
    } else {
      # All probabilities are 1 for survival functions at time 0
      bounds.list[[i]] <- rbind(c(1, 1, 1, 0), bounds.df)
    }
    
  }
  
  return(list(bounds.list=bounds.list, exec.times.list=exec.times.list))
}


#' Objective function to minimize or maximize
#'
#' @param init Vector to initialize the optimization
#' @param time The time point of interest
#' @param tk.unique.a.delta Vector of time points in the data (treated and uncensored units)
#' @param sim.data.df A dataframe containing the data
#' @param optim.prob Either "min" for a minimization problem, or "max" for a maximization problem
#' @param gamma The confounding strength
#' @param weights.a.delta1 Weights for treated and uncensored units
#' @param weights.a.delta0 Weights for treated and censored units
#' @param T.obs.a.delta1 Observed times for treated and uncensored units
#' @param T.obs.a.delta0 Observed times for treated and censored units
#' @param RMST If TRUE, computes the objective function for the RMST, else, computes the objective for the survival functions
#'
#' @return The survival function estimate
#' @keywords internal
#'
#' @examples
objectiveFun <- function(init, time, tk.unique.a.delta, sim.data.df,
                         optim.prob, gamma,
                         weights.a.delta1, weights.a.delta0,
                         T.obs.a.delta1, T.obs.a.delta0, RMST=FALSE) {
  
  z <- init
  
  if (optim.prob == "min") {
    gamma.fact <- 1 / gamma
  } else if (optim.prob == "max") {
    gamma.fact <- gamma
  } else {
    stop("optim.prob must be 'min' or 'max'")
  }
  
  KM.fun <- function(t) {
    
    num <- sum((1 + z * weights.a.delta1)[T.obs.a.delta1 == t])
    den <- sum((1 + z * weights.a.delta1)[T.obs.a.delta1 >= t]) + sum((1 + gamma.fact * weights.a.delta0)[T.obs.a.delta0 >= t])
    
    c(num=num, den=den)
  }
  
  tk.vec <- tk.unique.a.delta[tk.unique.a.delta <= time]
  
  # If the desired time point is lower than the minimum observed time
  if (sum(tk.unique.a.delta <= time) == 0) {
    if (RMST) {
      
      # The area under the curve is 0
      return(0)
      
    } else {
      
      # The probability is equal to 1
      return(1)
      
    }
  }
  
  num.den <- sapply(tk.vec, KM.fun)
  num <- num.den[1, ]
  den <- num.den[2, ]
  
  if (RMST) {
    
    # Compute the probability at each time point lower than "time"
    S.t <- cumprod(1 - num / den)
    # Sort by decreasing order
    S.t <- sort(S.t, decreasing=TRUE)
    # Midpoint rule (faster convergence rate than left or right rectangle method)
    S.t <- (c(1, S.t[-length(S.t)]) + S.t) / 2
    # Compute the time difference
    time.diff <- diff(c(0, tk.vec))
    # Compute the area by the rectangle method
    area <- sum(S.t * time.diff)
    
    return(area)
    
  } else {
    
    # Compute the probability at time point "time"
    S.t <- prod(1 - num / den)
    return(S.t)
  }
}


#' Function to launch the experiment with the bounds of Lee et al. (2024)
#'
#' @param version Experiment number
#' @param data.version Simulated data version number
#' @param data.name Either "simul", "rhc", or "gbsg"
#' @param gamma The sensitivity parameter/confounding strength
#' @param nb.folds Number of folds for cross-fitting
#' @param RMST Boolean. If TRUE, bounds for the RMST are computed. Else, bounds for the survival function are computed.
#' @param bootstrap Boolean. If TRUE, the bounds are computed for different bootstrap samples and confidence intervals are returned.
#' @param B Number of bootstrap samples
#' @param alpha Significance level of the confidence intervals
#'
#' @return
#' @export
leeMcExp <- function(version=1, data.version=2, data.name="simul", gamma=1,
                     nb.folds=5, RMST=FALSE, bootstrap=FALSE, B=100,
                     alpha=0.05) {
  
  if (data.name == "simul") {
    
    # Read list of data sets
    file.name <- paste0("./data/simul/data_list_v", data.version, ".rds")
    data.list <- readRDS(file.name)
    
    tk.unique <- seq(from=0.1, to=9, length.out=10)

  } else if (data.name == "rhc") {
    
    file.name <- paste0("./data/RHC/preprocessed_data.rds")
    data.list <- readRDS(file.name)
    
    tk.unique <- seq(from=min(data.list[[1]]$T.obs), to=50, length.out=5)
    tk.unique <- sort(c(tk.unique, 30))  # Add 30 days
    
    if (bootstrap) {
      data.list <- bootstrapFun(data.list=data.list, B=B)
    }
    
  } else if (data.name == "gbsg") {
    
    file.name <- paste0("./data/GBSG/preprocessed_data.rds")
    data.list <- readRDS(file.name)
    
    tk.unique <- seq(from=min(data.list[[1]]$T.obs), to=1826.25, length.out=5)

    if (bootstrap) {
      data.list <- bootstrapFun(data.list=data.list, B=B)
    }
    
  }
  
  start.time <- Sys.time()
  
  treated.bounds.lee.list <- survLeeMc(data.list=data.list, data.name=data.name,
                                       t.vec=tk.unique,
                                       a=1, gamma=gamma, nb.folds=nb.folds,
                                       RMST=RMST)
  
  control.bounds.lee.list <- survLeeMc(data.list=data.list, data.name=data.name,
                                       t.vec=tk.unique,
                                       a=0, gamma=gamma, nb.folds=nb.folds,
                                       RMST=RMST)
  
  end.time <- Sys.time()
  elapsed <- end.time - start.time
  message(elapsed)
  exec.time.list <- list(start=start.time, end=end.time, elapsed=elapsed)
  
  if (RMST) {
    measure <- "rmst"
  } else {
    measure <- "surv"
  }
  
  saveRDS(exec.time.list, file=paste0("results/", data.name, "_lee_", measure, "_tot_exec_time_v", version, ".rds"))
  saveRDS(control.bounds.lee.list$exec.times.list, file=paste0("results/", data.name, "_lee_", measure, "_control_exec_time_mc_v", version, ".rds"))
  saveRDS(treated.bounds.lee.list$exec.times.list, file=paste0("results/", data.name, "_lee_", measure, "_treated_exec_time_mc_v", version, ".rds"))
  
  # Transform list into dataframe
  treated.bounds.lee <- do.call(rbind.data.frame, treated.bounds.lee.list$bounds.list)
  control.bounds.lee <- do.call(rbind.data.frame, control.bounds.lee.list$bounds.list)
  
  # Save results
  saveRDS(treated.bounds.lee, file=paste0("results/", data.name, "_lee_bounds_", measure, "_mc_treated_v", version, ".rds"))
  saveRDS(control.bounds.lee, file=paste0("results/", data.name, "_lee_bounds_", measure, "_mc_control_v", version, ".rds"))
  
  # For confidence intervals
  col.names <- c("ub", "lb", "time")
  
  treated.bounds.ci <- foreach(k=1:length(tk.unique)+1) %dopar% {
    apply(X=treated.bounds.lee[treated.bounds.lee$time == c(0, tk.unique)[k],
                               col.names],
          MARGIN=2,
          FUN=quantile,
          probs=c(alpha/2, 1-alpha/2))
  }
  
  control.bounds.ci <- foreach(k=1:length(tk.unique)+1) %dopar% {
    apply(X=control.bounds.lee[control.bounds.lee$time == c(0, tk.unique)[k],
                               col.names],
          MARGIN=2,
          FUN=quantile,
          probs=c(alpha/2, 1-alpha/2))
  }
  
  saveRDS(treated.bounds.ci, file=paste0("results/", data.name, "_lee_bounds_ci_", measure, "_mc_treated_v", version, ".rds"))
  saveRDS(control.bounds.ci, file=paste0("results/", data.name, "_lee_bounds_ci_", measure, "_mc_control_v", version, ".rds"))
}
