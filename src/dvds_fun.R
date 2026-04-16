# Functions specific to the DVDS bounds


#' Function to compute the DVDS bounds for the survival function on Monte-Carlo samples. Works under parallel computation.
#'
#' @param data.list A list containing the Monte-Carlo samples
#' @param data.name Either "simul" or "rhc"
#' @param tk.unique A vector of times for which to compute the bounds of the survival function
#' @param a If a=1, computes the bounds for the survival function among the treated and, if a=0, among the control
#' @param gamma Sensitivity parameter/Confounding strength
#' @param nb.folds Number of folds for cross-fitting
#' @param use.U Whether to use the (unobserved) values of U
#' @param num.trees Number of trees for random forests
#' @param cq.method Method to compute the conditional quantiles. Default is "fn".
#' @param RMST If TRUE, computes the bounds for the RMST, else, computes the bounds for the survival functions
#'
#' @return A list of the bounds for each Monte-Carlo sample at each time point in tk.unique
#' @export
#'
#' @examples
survDvdsMc <- function(data.list, data.name="simul",
                       tk.unique=seq(from=0.1, to=9, length.out=10),
                       a=1, gamma=1, nb.folds=5, use.U=FALSE, num.trees=500,
                       cq.method="fn", RMST=FALSE) {
  
  require(foreach)
  
  # Number of Monte-Carlo samples
  n.mc <- length(data.list)
  # Number of unique time-points
  tk.unique.length <- length(tk.unique)
  tau <- gamma / (1 + gamma)
  # Initialize empty list that will contain the estimations
  bounds.list <- list()
  
  for (j in 1:n.mc) {
    
    # Retrieve dataset i
    simul.data.list <- data.list[[j]]
    p.X <- ncol(simul.data.list$X)
    p.U <- ncol(simul.data.list$U)
    n <- length(simul.data.list$A)
    
    if (data.name == "simul") {
      Delta <- 1 * (simul.data.list$T01 <= simul.data.list$C)
    } else if (data.name == "rhc") {
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
    
    
    #####################################
    ### Censoring function estimation ###
    #####################################
    
    if (a == 0) {
      # Invert treated and control group
      sim.data.df$A <- 1 - sim.data.df$A
    }
    
    # # Estimation via a Cox model for A=1 (treated)
    G.predict1 <- censoringFun(sim.data.df=sim.data.df, p.X=p.X, p.U=p.U,
                               nb.folds=nb.folds, use.U=use.U, t.predict=NULL)
    
    if (RMST) {
      unique.T.obs <- unique(sim.data.df$T.obs)
      unique.T.obs.inf.max.tk <- unique.T.obs[unique.T.obs <= max(tk.unique)]
      G.predict2 <- censoringFun(sim.data.df=sim.data.df, p.X=p.X, p.U=p.U,
                                 nb.folds=nb.folds, use.U=use.U,
                                 t.predict=unique.T.obs.inf.max.tk)
    } else {
      G.predict2 <- censoringFun(sim.data.df=sim.data.df, p.X=p.X, p.U=p.U,
                                 nb.folds=nb.folds, use.U=use.U,
                                 t.predict=tk.unique)
    }

    
    ###################################
    ### Propensity score estimation ###
    ###################################
    
    e.predict <- propScore(sim.data.df=sim.data.df, p.X=p.X, nb.folds=nb.folds)
    weights <- (1 - e.predict) / e.predict
    
    if (RMST) {
      
      bounds.df <- data.frame(ub1=rep(NA, tk.unique.length),
                              lb1=rep(NA, tk.unique.length),
                              ign.vec1=rep(NA, tk.unique.length),
                              ub2=rep(NA, tk.unique.length),
                              lb2=rep(NA, tk.unique.length),
                              ign.vec2=rep(NA, tk.unique.length),
                              time=rep(NA, tk.unique.length))
      
      # For each time for which we want the RMST
      for (i in 1:tk.unique.length) {
        
        # Get all the (unique) observed times that are lower than tk
        unique.T.obs.inf.tk <- unique.T.obs[unique.T.obs <= tk.unique[i]]
        
        if (length(unique.T.obs.inf.tk) == 0) {
          
          bounds.df[i, ] <- c(0, 0, 0, 0, 0, 0, tk.unique[i])
          
        } else {
          
          # Get the bounds for all these time points
          bounds.rect.df <- boundsFun(sim.data.df=sim.data.df, n=n, p.X=p.X,
                                      tk.unique=unique.T.obs.inf.tk,
                                      gamma=gamma,
                                      G.predict1=G.predict1,
                                      G.predict2=G.predict2,
                                      weights=weights, num.trees=num.trees,
                                      nb.folds=nb.folds, cq.method=cq.method)
          
          # Add bounds for time 0 (the probability is 1)
          bounds.rect.df <- rbind(c(1, 1, 1, 1, 1, 1, 0), bounds.rect.df)
          
          # Get time difference for rectangle method
          time.diff <- diff(c(0, unique.T.obs.inf.tk))
          
          midpoint.lb1 <- (bounds.rect.df$lb1[-length(bounds.rect.df$lb1)] + bounds.rect.df$lb1[-1]) / 2
          lb1 <- sum(midpoint.lb1 * time.diff)
          midpoint.ub1 <- (bounds.rect.df$ub1[-length(bounds.rect.df$ub1)] + bounds.rect.df$ub1[-1]) / 2
          ub1 <- sum(midpoint.ub1 * time.diff)
          midpoint.ign.vec1 <- (bounds.rect.df$ign.vec1[-length(bounds.rect.df$ign.vec1)] + bounds.rect.df$ign.vec1[-1]) / 2
          ign.vec1 <- sum(midpoint.ign.vec1 * time.diff)
          
          midpoint.lb2 <- (bounds.rect.df$lb2[-length(bounds.rect.df$lb2)] + bounds.rect.df$lb2[-1]) / 2
          lb2 <- sum(midpoint.lb2 * time.diff)
          midpoint.ub2 <- (bounds.rect.df$ub2[-length(bounds.rect.df$ub2)] + bounds.rect.df$ub2[-1]) / 2
          ub2 <- sum(midpoint.ub2 * time.diff)
          midpoint.ign.vec2 <- (bounds.rect.df$ign.vec2[-length(bounds.rect.df$ign.vec2)] + bounds.rect.df$ign.vec2[-1]) / 2
          ign.vec2 <- sum(midpoint.ign.vec2 * time.diff)
          
          bounds.df[i, ] <- c(ub1, lb1, ign.vec1,
                              ub2, lb2, ign.vec2,
                              tk.unique[i])
          
        }
      }
      
      bounds.list[[j]] <- rbind(c(0, 0, 0, 0, 0, 0, 0), bounds.df)
      
    } else {
      
      bounds.df <- boundsFun(sim.data.df=sim.data.df, n=n, p.X=p.X, tk.unique=tk.unique,
                             gamma=gamma,
                             G.predict1=G.predict1, G.predict2=G.predict2,
                             weights=weights, num.trees=num.trees,
                             nb.folds=nb.folds, cq.method=cq.method)
      
      bounds.list[[j]] <- rbind(c(1, 1, 1, 1, 1, 1, 0), bounds.df)
    }
  }
  
  return(bounds.list)
}


#' Function that computes the DVDS bounds
#'
#' @param sim.data.df A data.frame object containing the data
#' @param n The sample size
#' @param p.X The dimension of X
#' @param gamma The sensitivity parameter
#' @param G.predict1 Censoring function predictions for Estimator I
#' @param G.predict2 Censoring function predictions for Estimator II
#' @param weights Weights involving the propensity score in our bounds
#' @param num.trees Number of trees for random forests
#' @param tk.unique.vec Vector of unique times for which to compute the intervals
#' @param nb.folds Number of folds for cross-fitting
#' @param cq.method Method to compute the conditional quantiles
#'
#' @return A dataframe with the bounds
#' @keywords internal
boundsFun <- function(sim.data.df, n, p.X, tk.unique.vec, gamma, G.predict1,
                      G.predict2, weights, num.trees, nb.folds, cq.method) {
  
  tk.unique.vec.length <- length(tk.unique.vec)
  tau <- gamma / (1 + gamma)
  
  # For the progress bar
  pb <- txtProgressBar(max=tk.unique.vec.length, style=3)
  progress <- function(N) {
    if (N != tk.unique.vec.length) {
      setTxtProgressBar(pb, N)
    } else {
      setTxtProgressBar(pb, N)
      close(pb)
    }
  }
  opts <- list(progress=progress)
  
  bounds.df <- foreach(k=1:tk.unique.vec.length,
                       .combine="rbind",
                       .options.snow=opts,
                       .export=c("outcomeRegression", "condQuantile", "quantileRegression", "randomForestFun")) %dopar% {
                         
                         indic.tk <- sim.data.df$T.obs > tk.unique.vec[k]
                         # Equivalent of Y for estimator 1
                         sim.data.df$Y.est1 <- sim.data.df$Delta * (1 - indic.tk) / G.predict1
                         # Equivalent of Y for estimator 2
                         sim.data.df$Y.est2 <- indic.tk / G.predict2[k, ]
                         
                         X.formula <- paste("X", 1:p.X, sep=".", collapse=" + ")
                         Q.formula1 <- paste("Y.est1 ~", X.formula, "+ A")
                         Q.formula2 <- paste("Y.est2 ~", X.formula, "+ A")
                         
                         Q.predict1 <- condQuantile(sim.data.df=sim.data.df, p.X=p.X,
                                                    Q.formula=Q.formula1, tau=tau,
                                                    nb.folds=nb.folds, method=cq.method)
                         
                         Q.predict2 <- condQuantile(sim.data.df=sim.data.df, p.X=p.X,
                                                    Q.formula=Q.formula2, tau=tau,
                                                    nb.folds=nb.folds, method=cq.method)
                         
                         # To avoid dimension error when Gamma=1
                         if (gamma == 1) {
                           
                           Q.lb1 <- Q.ub1 <- Q.predict1
                           Q.lb2 <- Q.ub2 <- Q.predict2
                           
                         } else {
                           
                           Q.lb1 <- Q.predict1[, 1]
                           Q.ub1 <- Q.predict1[, 2]
                           
                           Q.lb2 <- Q.predict2[, 1]
                           Q.ub2 <- Q.predict2[, 2]
                           
                         }
                         
                         # Estimator I
                         sim.data.df$Y.rho.lb1 <- sim.data.df$Y.est1 / gamma + (1 - 1 / gamma) * (Q.lb1 + pmin(rep(0, n), sim.data.df$Y.est1 - Q.lb1) / (1 - tau))
                         rho.lb.formula1 <- as.formula(paste("Y.rho.lb1 ~", paste("X", 1:p.X, sep=".", collapse=" + "), "+ A"))

                         rho.lb1 <- outcomeRegression(sim.data.df=sim.data.df,
                                                      rho.formula=rho.lb.formula1,
                                                      nb.folds=nb.folds,
                                                      num.trees=num.trees)
                         
                         sim.data.df$Y.rho.ub1 <- sim.data.df$Y.est1 / gamma + (1 - 1 / gamma) * (Q.ub1 + pmax(rep(0, n), sim.data.df$Y.est1 - Q.ub1) / (1 - tau))
                         rho.ub.formula1 <- as.formula(paste("Y.rho.ub1 ~", paste("X", 1:p.X, sep=".", collapse=" + "), "+ A"))

                         rho.ub1 <- outcomeRegression(sim.data.df=sim.data.df,
                                                      rho.formula=rho.ub.formula1,
                                                      nb.folds=nb.folds,
                                                      num.trees=num.trees)
                         
                         # Estimator II
                         sim.data.df$Y.rho.lb2 <- sim.data.df$Y.est2 / gamma + (1 - 1 / gamma) * (Q.lb2 + pmin(rep(0, n), sim.data.df$Y.est2 - Q.lb2) / (1 - tau))
                         rho.lb.formula2 <- as.formula(paste("Y.rho.lb2 ~", paste("X", 1:p.X, sep=".", collapse=" + "), "+ A"))

                         rho.lb2 <- outcomeRegression(sim.data.df=sim.data.df,
                                                      rho.formula=rho.lb.formula2,
                                                      nb.folds=nb.folds,
                                                      num.trees=num.trees)
                         
                         sim.data.df$Y.rho.ub2 <- sim.data.df$Y.est2 / gamma + (1 - 1 / gamma) * (Q.ub2 + pmax(rep(0, n), sim.data.df$Y.est2 - Q.ub2) / (1 - tau))
                         rho.ub.formula2 <- as.formula(paste("Y.rho.ub2 ~", paste("X", 1:p.X, sep=".", collapse=" + "), "+ A"))

                         rho.ub2 <- outcomeRegression(sim.data.df=sim.data.df,
                                                      rho.formula=rho.ub.formula2,
                                                      nb.folds=nb.folds,
                                                      num.trees=num.trees)
                         
                         # Estimator I (invert lb and ub because of the "-" sign)
                         ub1 <- 1 - mean(sim.data.df$A * sim.data.df$Y.est1 + (1-sim.data.df$A) * rho.lb1 + sim.data.df$A * weights * (Q.lb1 + gamma**(- sign(sim.data.df$Y.est1 - Q.lb1)) * (sim.data.df$Y.est1 - Q.lb1) - rho.lb1))
                         lb1 <- 1 - mean(sim.data.df$A * sim.data.df$Y.est1 + (1-sim.data.df$A) * rho.ub1 + sim.data.df$A * weights * (Q.ub1 + gamma**(sign(sim.data.df$Y.est1 - Q.ub1)) * (sim.data.df$Y.est1 - Q.ub1) - rho.ub1))
                         
                         # Compute estimator I under ignorability
                         # Compute the modified outcome regression
                         rho.formula1 <- as.formula(paste("Y.est1 ~", paste("X", 1:p.X, sep=".", collapse=" + "), "+ A"))
                         
                         rho1 <- outcomeRegression(sim.data.df=sim.data.df,
                                                   rho.formula=rho.formula1,
                                                   nb.folds=nb.folds,
                                                   num.trees=num.trees)
                         
                         ign.vec1 <- 1 - mean(sim.data.df$A * sim.data.df$Y.est1 + (1-sim.data.df$A) * rho1 + sim.data.df$A * weights * (sim.data.df$Y.est1 - rho1))
                         
                         
                         # Estimator II
                         ub2 <- mean(sim.data.df$A * sim.data.df$Y.est2 + (1-sim.data.df$A) * rho.ub2 + sim.data.df$A * weights * (Q.ub2 + gamma**(sign(sim.data.df$Y.est2 - Q.ub2)) * (sim.data.df$Y.est2 - Q.ub2) - rho.ub2))
                         lb2 <- mean(sim.data.df$A * sim.data.df$Y.est2 + (1-sim.data.df$A) * rho.lb2 + sim.data.df$A * weights * (Q.lb2 + gamma**(- sign(sim.data.df$Y.est2 - Q.lb2)) * (sim.data.df$Y.est2 - Q.lb2) - rho.lb2))
                         
                         # Compute estimator II under ignorability
                         # Compute the modified outcome regression
                         rho.formula2 <- as.formula(paste("Y.est2 ~", paste("X", 1:p.X, sep=".", collapse=" + "), "+ A"))
                         
                         rho2 <- outcomeRegression(sim.data.df=sim.data.df,
                                                   rho.formula=rho.formula2,
                                                   nb.folds=nb.folds,
                                                   num.trees=num.trees)
                         
                         ign.vec2 <- mean(sim.data.df$A * sim.data.df$Y.est2 + (1-sim.data.df$A) * rho2 + sim.data.df$A * weights * (sim.data.df$Y.est2 - rho2))
                         
                         bounds.df <- data.frame(ub1=ub1, lb1=lb1, ign.vec1=ign.vec1,
                                                 ub2=ub2, lb2=lb2, ign.vec2=ign.vec2,
                                                 time=tk.unique.vec[k])
                         
                         return(bounds.df)
                       }
  
  return(bounds.df)
}


#' Function to launch the experiment with the DVDS bounds
#'
#' @param version Experiment number
#' @param data.version Simulated data version number
#' @param data.name Either "simul" or "rhc"
#' @param gamma The sensitivity parameter/confounding strength
#' @param nb.folds Number of folds for cross-fitting
#' @param num.trees Number of trees in the random forest
#' @param RMST Boolean. If TRUE, bounds for the RMST are computed. Else, bounds for the survival function are computed.
#'
#' @return
#' @export
dvdsMcExp <- function(version=1, data.version=2, data.name="simul", gamma=1,
                      nb.folds=5, num.trees=500, RMST=FALSE) {
  
  use.U <- FALSE
  
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
    
  }
  
  start.time <- Sys.time()
  
  treated.bounds.list <- survDvdsMc(data.list=data.list, data.name=data.name, tk.unique=tk.unique,
                                    a=1, gamma=gamma, nb.folds=nb.folds,
                                    use.U=use.U, num.trees=num.trees,
                                    RMST=RMST)
  
  control.bounds.list <- survDvdsMc(data.list=data.list, data.name=data.name, tk.unique=tk.unique,
                                    a=0, gamma=gamma, nb.folds=nb.folds,
                                    use.U=use.U, num.trees=num.trees,
                                    RMST=RMST)
  
  end.time <- Sys.time()
  elapsed <- end.time - start.time
  message(elapsed)
  exec.time.list <- list(start=start.time, end=end.time, elapsed=elapsed)
  
  # Transform list into dataframe
  treated.bounds <- do.call(rbind.data.frame, treated.bounds.list)
  control.bounds <- do.call(rbind.data.frame, control.bounds.list)
  
  if (RMST) {
    measure <- "rmst"
  } else {
    measure <- "surv"
  }
  
  saveRDS(treated.bounds, file=paste0("results/", data.name, "_dvds_bounds_", measure, "_mc_treated_v", version, ".rds"))
  saveRDS(control.bounds, file=paste0("results/", data.name, "_dvds_bounds_", measure, "_mc_control_v", version, ".rds"))
  saveRDS(exec.time.list, file=paste0("results/", data.name, "_dvds_", measure, "_tot_exec_time_v", version, ".rds"))
  
}
