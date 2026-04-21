# Different functions used to estimate the bounds on S^{(a)}(t), plot or initialize parallel computing


#' Function to initialize a cluster for parallel computation
#'
#' @return A cluster
#' @export
initCluster <- function() {
  
  nb.cores <- parallel::detectCores() - 1  # Get number of cores and subtract 1 (DO NOT USE ALL THE CORES)
  cl <- parallel::makeCluster(nb.cores)  # Create the cluster
  invisible(parallel::clusterExport(cl, c("quantileRegression", "coxFun", "randomForestFun")))  # Export functions
  # invisible(parallel::clusterEvalQ(cl, {library(libraryName)}))  # Uncomment to register a library in the cluster
  
  #########################################################################################################
  # THE FOLLOWING LINES OF CODE WERE MODIFIED FOR LICENSING REASONS                                       #
  # CHANGE THEM TO ENABLE PARALLEL COMPUTING                                                              #
  # FOR EXAMPLE, USE X::registerX(cl), WHERE X IS THE NAME OF THE PACKAGE USED TO REGISTER THE CLUSTER cl #
  # WHEN RUNNING CODE, IF A FUNCTION FROM A LIBRARY libraryName IS NOT RECOGNIZED, REGISTER IT WITH       #
  # invisible(parallel::clusterEvalQ(cl, {library(libraryName)}))                                         #
  #########################################################################################################
  
  # ======> CHANGE CODE FROM HERE...
  message("IMPORTANT: To enable parallel computing, see instructions in initCluster() function in src/utils.R.")
  
  foreach::registerDoSEQ()
  # ======> ... TO HERE
  
  return(cl)
}


#' Function to close a cluster and stop parallel computation
#'
#' @param cl A cluster to close
#'
#' @export
closeCluster <- function(cl) {
  parallel::stopCluster(cl)
}


#' A function that fits a linear quantile regression
#'
#' @param formula A formula object or string of the form Y ~ X1 + X2 + ..., where Y is the response and X1, X2, ... are the explanatory variables
#' @param quantiles A vector of quantiles to be estimated (numbers strictly between 0 and 1)
#' @param data A data.frame object containing the variables named in the formula
#' @param method A string indicating the quantile regression algorithm. Default is "fn".
#'
#' @return A quantile regression fit object
#' @export
quantileRegression <- function(formula, quantiles, data, method="fn") {
  
  #########################################################################################################
  # THE FOLLOWING LINES OF CODE WERE MODIFIED FOR LICENSING REASONS                                       #
  # CHANGE THEM TO ENABLE QUANTILE REGRESSION                                                             #
  # CHANGE THE NAME OF THE LIBRARY IN THE require() FUNCTION TO IMPORT THE QUANTILE REGRESSION FITTING    #
  # AND PREDICTION FUNCTIONS OF YOUR CHOICE                                                               #
  # CHANGE THE NAME OF THE FUNCTION qrFun AND THE NAME OF THE DIFFERENT ARGUMENTS WHEN CALLING qrFun      #
  # (formula, quantiles, data, and method) ACCORDING TO YOUR CHOICE                                       #
  # FOR EXAMPLE, YOU CAN CHOOSE A LINEAR QUANTILE REGRESSION FUNCTION                                     #
  # ONCE MODIFICATIONS ARE DONE, CHANGE stop() INTO message() OR REMOVE IT                                #
  #########################################################################################################
  
  # ======> CHANGE CODE FROM HERE...
  stop("IMPORTANT: A (linear) quantile regression fitting function must be provided by the user.\nSee function 'quantileRegression' in src/utils.R.")
  
  require(quantileRegressionLibrary)
  
  qr.fit <- qrFun(formula=formula,
                  quantiles=quantiles,
                  data=data,
                  method=method)
  # ======> ... TO HERE
  
  return(qr.fit)
}


#' A function that fits a Cox proportional hazards model
#'
#' @param formula A formula object or string of the form Surv(time, status) ~ X1 + X2 + ..., where time is the observed time, status is the censoring indicator, and X1, X2, ... are the explanatory variables
#' @param data A data.frame object containing the variables named in the formula
#' @param ties A method to handle ties (multiple individuals having an event at the same time)
#'
#' @return A survival model fit object
#' @export
coxFun <- function(formula, data, ties="breslow") {
  
  #########################################################################################################
  # THE FOLLOWING LINES OF CODE WERE MODIFIED FOR LICENSING REASONS                                       #
  # CHANGE THEM TO ENABLE THE ESTIMATION OF THE COX MODEL                                                 #
  # CHANGE THE NAME OF THE LIBRARY IN THE require() FUNCTION TO IMPORT THE COX MODEL                      #
  # AND PREDICTION FUNCTIONS OF YOUR CHOICE                                                               #
  # CHANGE THE NAME OF THE FUNCTION cmFun AND THE NAME OF THE DIFFERENT ARGUMENTS WHEN CALLING cmFun      #
  # (formula and data) ACCORDING TO YOUR CHOICE                                                           #
  # ONCE MODIFICATIONS ARE DONE, CHANGE stop() INTO message() OR REMOVE IT                                #
  #########################################################################################################
  
  # ======> CHANGE CODE FROM HERE...
  stop("IMPORTANT: A Cox proportional hazards model must be provided by the user.\nSee function 'coxFun' in src/utils.R.")
  
  require(coxModelLibrary)
  
  cm.fit <- cmFun(formula=formula,
                  data=data,
                  ties=ties)
  # ======> ... TO HERE
  
  return(cm.fit)
}


#' A function that fits a random forest model
#'
#' @param formula A formula object or string of the form Y ~ X1 + X2 + ..., where Y is the response and X1, X2, ... are the explanatory variables
#' @param data A data.frame object containing the variables named in the formula
#' @param num.trees The number of trees in the random forest
#'
#' @return A random forest fit object
#' @export
randomForestFun <- function(formula, data, num.trees) {
  
  #########################################################################################################
  # THE FOLLOWING LINES OF CODE WERE MODIFIED FOR LICENSING REASONS                                       #
  # CHANGE THEM TO ENABLE THE ESTIMATION OF THE RANDOM FOREST MODEL                                       #
  # CHANGE THE NAME OF THE LIBRARY IN THE require() FUNCTION TO IMPORT THE RANDOM FOREST MODEL            #
  # AND PREDICTION FUNCTIONS OF YOUR CHOICE                                                               #
  # CHANGE THE NAME OF THE FUNCTION rfFun AND THE NAME OF THE DIFFERENT ARGUMENTS WHEN CALLING rfFun      #
  # (formula, data, and n.trees) ACCORDING TO YOUR CHOICE                                                 #
  # ONCE MODIFICATIONS ARE DONE, CHANGE stop() INTO message() OR REMOVE IT                                #
  #########################################################################################################
  
  # ======> CHANGE CODE FROM HERE...
  stop("IMPORTANT: A random forest model must be provided by the user.\nSee function 'randomForestFun' in src/utils.R.")
  
  require(randomForestLibrary)
  
  rf.fit <- rfFun(formula=formula,
                  data=data,
                  n.trees=n.trees)
  # ======> ... TO HERE
  
  return(rf.fit)
}


#' Function to estimate the censoring function G by cross-fitting
#'
#' @param sim.data.df A dataframe containing the data
#' @param p.X Dimension of X
#' @param p.U Dimension of U
#' @param nb.folds Number of folds for cross-fitting
#' @param folds.ind Folds indices for K-fold cross-fitting
#' @param use.U Boolean. If TRUE, use U to predict the censoring function in the simulations.
#' @param t.predict Vector of times for prediction
#' @param censoring Boolean. If TRUE, computes the survival function. Else, computes the cumulative distribution function.
#' @param fast.RMST Boolean. If TRUE, the censoring function is computed for fast RMST.
#'
#' @return A vector of length n, the sample size of the data
#' @keywords internal
censoringFun <- function(sim.data.df, p.X, p.U, nb.folds, folds.ind=NULL, use.U=FALSE,
                         t.predict=NULL, censoring=TRUE, fast.RMST=FALSE) {
  
  require(survival)
  
  n <- nrow(sim.data.df)
  
  predict.data <- sim.data.df
  
  if (!is.null(t.predict)) {
    t.length <- length(t.predict)
  }
  
  if (censoring) {
    surv.string <- "Surv(T.obs, 1-Delta) ~"
  } else {
    surv.string <- "Surv(T.obs, Delta) ~"
  }
  
  if (use.U) {
    
    G.formula <- paste(surv.string, paste("X.", 1:p.X, sep="", collapse=" + "),
                       "+", paste("U.", 1:p.U, sep="", collapse=" + "),
                       "+", "A")

  } else {
    
    G.formula <- paste(surv.string, paste("X.", 1:p.X, sep="", collapse=" + "),
                       "+", "A")
    
  }
  
  if (nb.folds == 1) {
    
    if (is.null(t.predict)) {
      
      G.model <- coxFun(as.formula(G.formula), data=sim.data.df)
      G.predict <- predict(G.model, predict.data, type="survival")
      
    } else {
      
      G.predict <- matrix(NA, nrow=t.length, ncol=n)
      G.model <- coxFun(as.formula(G.formula), data=sim.data.df)
      
      for (j in 1:t.length) {
        
        if (fast.RMST) {
          predict.data$T.obs <- pmin(sim.data.df$T.obs, t.predict[j])
        } else {
          predict.data$T.obs <- t.predict[j]
        }
        
        G.predict[j, ] <- predict(G.model, predict.data, type="survival")
        
      }
      
    }
    
  } else { # K-fold cross-fitting
    
    if (is.null(t.predict)) {
      
      G.predict <- rep(NA, n)

      for (i in 1:nb.folds) {
        
        # Fit a Cox model and predict the values for each fold
        G.model <- coxFun(as.formula(G.formula), data=sim.data.df[-folds.ind[[i]], ])
        G.predict[folds.ind[[i]]] <- predict(G.model, predict.data[folds.ind[[i]], ], type="survival")
        
      }
      
    } else {
      
      G.predict <- matrix(NA, nrow=t.length, ncol=n)
      G.models <- list()
      
      for (i in 1:nb.folds) {
        
        # Fit a Cox model and predict the values for each fold
        G.model <- coxFun(as.formula(G.formula), data=sim.data.df[-folds.ind[[i]], ])
        G.models[[i]] <- G.model
        
      }
      
      for (j in 1:t.length) {
        
        if (fast.RMST) {
          predict.data$T.obs <- pmin(sim.data.df$T.obs, t.predict[j])
        } else {
          predict.data$T.obs <- t.predict[j]
        }
        
        for (i in 1:nb.folds) {
          
          G.predict[j, folds.ind[[i]]] <- predict(G.models[[i]], predict.data[folds.ind[[i]], ], type="survival")
          
        }
        
      }
      
    }
    
  }
  
  return(G.predict)
}


#' Function to estimate the propensity score
#'
#' @param sim.data.df A dataframe containing the data
#' @param p.X Dimension of X
#' @param nb.folds Number of folds for cross-fitting
#' @param folds.ind Folds indices for K-fold cross-fitting
#'
#' @return A vector of length n, the sample size of the data
#' @keywords internal
propScore <- function(sim.data.df, p.X, nb.folds, folds.ind=NULL) {
  
  n <- nrow(sim.data.df)
  
  # Estimate the propensity score with all observed covariates X
  e.X.formula <- paste("A ~", paste("X", 1:p.X, sep=".", collapse=" + "))
  
  if (nb.folds == 1) {
    
    # Fit a logistic model and predict the values
    e.X.predict <- unname(glm(e.X.formula, family=binomial(link='logit'),
                              data=sim.data.df)$fitted.values)
    
  } else {
    
    e.X.predict <- rep(NA, n)
    
    for (i in 1:nb.folds) {
      
      # Fit a logistic model and predict the values for each fold
      e.X.model <- glm(e.X.formula, family=binomial(link='logit'),
                       data=sim.data.df[-folds.ind[[i]], ])
      e.X.predict[folds.ind[[i]]] <- unname(predict(e.X.model,
                                                    sim.data.df[folds.ind[[i]], ],
                                                    type="response"))
      
    }
    
  }
  
  return(e.X.predict)
}


#' Function to estimate the conditional quantiles
#'
#' @param sim.data.df A dataframe containing the data
#' @param p.X Dimension of X
#' @param Q.formula Formula for the quantile regression model
#' @param tau 1 / (1 + Gamma)
#' @param nb.folds Number of folds for cross-fitting
#' @param method Method to compute the conditional quantiles. Default is "fn".
#' @param folds.ind Folds indices for K-fold cross-fitting
#' @param cq.type Type of conditional quantiles. Either "linear" or "rf" (random forest)
#'
#' @return Either a vector (if tau=1/2) or a matrix
#' @keywords internal
condQuantile <- function(sim.data.df, p.X, Q.formula, tau,
                         nb.folds, folds.ind=NULL,
                         method="fn", cq.type=c("linear", "rf")) {
  
  n <- nrow(sim.data.df)
  
  if (nb.folds == 1) {
    
    if (cq.type == "linear") {
      
      # Fit a quantile regression model and predict the values
      Q.model <- quantileRegression(formula=Q.formula,
                                    quantiles=c(tau, 1-tau),
                                    data=sim.data.df,
                                    method=method)
      
      Q.predict <- predict(Q.model,
                           data.frame(sim.data.df[, c(paste("X", 1:p.X, sep="."), "A")]))
      
    } else if (cq.type == "rf") {
      
      stop('"rf" method for cq.type not implemented')
      
    } else {
      stop('cq.type must be "linear" or "rf"')
    }
    
  } else {
    
    if (cq.type == "linear") {
      
      if (tau == 1/2) {
        # When gamma=1, we only have 1 value/column
        Q.predict <- rep(NA, n)
      } else {
        # Otherwise, we have 2 values/columns
        Q.predict <- matrix(NA, nrow=n, ncol=2)
      }
      
      for (i in 1:nb.folds) {
        
        # Fit a quantile regression model on each fold except the i-th and predict the values on the i-th fold
        Q.model <- quantileRegression(formula=Q.formula,
                                      quantiles=c(tau, 1-tau),
                                      data=sim.data.df[-folds.ind[[i]], ],
                                      method=method)
        
        if (tau == 1/2) {
          Q.predict[folds.ind[[i]]] <- predict(Q.model,
                                               data.frame(sim.data.df[folds.ind[[i]],
                                                                      c(paste("X", 1:p.X, sep="."), "A")]))
        } else {
          Q.predict[folds.ind[[i]], ] <- predict(Q.model,
                                                 data.frame(sim.data.df[folds.ind[[i]],
                                                                        c(paste("X", 1:p.X, sep="."), "A")]))
        }
        
      }
      
    } else if (cq.type == "rf") {
      
      stop('"rf" method for cq.type not implemented')
      
    } else {
      stop('cq.type must be "linear" or "rf"')
    }
    
  }
  
  return(Q.predict)
}


#' Function to estimate the modified outcome regression
#'
#' @param sim.data.df A dataframe containing the data
#' @param rho.formula Formula for the random forest model
#' @param nb.folds Number of folds for cross-fitting
#' @param folds.ind Folds indices for K-fold cross-fitting
#' @param num.trees Number of trees for random forests
#'
#' @return A vector of predictions
#' @keywords internal
outcomeRegression <- function(sim.data.df, rho.formula, nb.folds, folds.ind=NULL, num.trees) {
  
  n <- nrow(sim.data.df)
  
  if (nb.folds == 1) {
    
    # Fit a random forest model and predict the values
    rho.model <- randomForestFun(rho.formula, data=sim.data.df, num.trees=num.trees)
    
    # For at least rho, I need to predict on A==1 because there is a factor (1-A) before this term
    rho.pred.data <- sim.data.df
    rho.pred.data$A <- 1
    
    rho.predict <- predict(rho.model, rho.pred.data)$predictions
    
  } else {
    
    rho.predict <- rep(NA, n)
    
    for (i in 1:nb.folds) {
      
      # Fit a random forest model on each fold except the i-th and predict the values on the i-th fold
      rho.model <- randomForestFun(rho.formula, data=sim.data.df[-folds.ind[[i]], ],
                                   num.trees=num.trees)
      
      rho.pred.data <- sim.data.df[folds.ind[[i]], ]
      rho.pred.data$A <- 1
      
      rho.predict[folds.ind[[i]]] <- predict(rho.model, rho.pred.data)$predictions
      
    }
  }
  
  return(rho.predict)
}


#' Function to create bootstrap resamples for confidence interval computation
#'
#' @param data.list The list of data to resample
#' @param B Number of bootstrap resamples
#'
#' @returns A list of B bootstrap resamples
#' @export
bootstrapFun <- function(data.list, B=100) {
  
  # Retrieve list of data
  data.list1 <- data.list[[1]]
  
  # Get number of individuals
  n.obs <- length(data.list1$A)
  
  btstrp.list <- foreach(b=1:B) %dopar% {
    
    bootstrap.ind <- sample(1:n.obs, n.obs, replace=TRUE)
    
    if (is.null(dim(data.list1$U))) {
      U <- data.list1$U[bootstrap.ind]
    } else {
      U <- data.list1$U[bootstrap.ind, ]
    }
    
    return(list(X=data.list1$X[bootstrap.ind, ],
                U=U,
                T1=data.list1$T1[bootstrap.ind],
                T0=data.list1$T0[bootstrap.ind],
                T1.distrib.fun=NULL,
                T0.distrib.fun=NULL,
                T01=data.list1$T01[bootstrap.ind],
                C=data.list1$C[bootstrap.ind],
                T.obs=data.list1$T.obs[bootstrap.ind],
                A=data.list1$A[bootstrap.ind],
                Delta=data.list1$Delta[bootstrap.ind]))
  }
  
  return(btstrp.list)
}


#' Function to save the plots for the different experiments
#'
#' @param dvds.version Experiment number for the DVDS bounds
#' @param lee.version Experiment number for the bounds of Lee et al. (2024)
#' @param dvds.boot.version Experiment number for the DVDS bounds with bootstrap samples
#' @param lee.boot.version Experiment number for the bounds of Lee et al. (2024) with bootstrap samples
#' @param data.name Either "simul" or "rhc"
#' @param data.version Simulated data version number
#' @param gamma The sensitivity parameter/confounding strength
#' @param RMST Boolean. If TRUE, bounds for the RMST were computed. Else, bounds for the survival function were computed.
#' @param fast.RMST If TRUE, the fast version of the RMST was used (no rectangle integration)
#' @param bootstrap Boolean. If TRUE, the bounds are computed for different bootstrap samples and confidence intervals are returned.
#'
#' @return
#' @export
plotFun <- function(dvds.version=1, lee.version=1,
                    dvds.boot.version=1, lee.boot.version=1, data.name="simul",
                    data.version=2, gamma=1, RMST=FALSE, fast.RMST=FALSE,
                    bootstrap=FALSE) {
  
  width <- 6.79
  height <- 5.38
  
  est1.bnds.col <- "blue"
  est1.ign.col <- alpha("blue", 0.5)
  est2.bnds.col <- "green3"
  est2.ign.col <- alpha("green3", 0.5)
  lee.bnds.col <- "red"
  lee.ign.col <- alpha("red", 0.5)
  true.color <- "magenta"
  
  est1.label <- "Estimator I"
  est2.label <- "Estimator II"
  lee.label <- "Lee et al. (2024)"
  true.label <- "Ground truth"
  
  colors <- c("Estimator I"=est1.bnds.col, "Estimator II"=est2.bnds.col,
              "Lee et al. (2024)"=lee.bnds.col, "Ground truth"=true.color)
  
  if (RMST && !fast.RMST) {
    measure <- "rmst"
    fast <- ""
  } else if (RMST && fast.RMST) {
    measure <- "rmst"
    fast <- "fast_"
  } else {
    measure <- "surv"
    fast <- ""
  }
  
  # Read list of data sets
  file.name <- paste0("./data/simul/data_list_v", data.version, ".rds")
  data.list <- readRDS(file.name)
  
  treated.bounds.lee <- readRDS(paste0("results/", data.name, "_lee_bounds_", measure, "_mc_treated_v", lee.version, ".rds"))
  control.bounds.lee <- readRDS(paste0("results/", data.name, "_lee_bounds_", measure, "_mc_control_v", lee.version, ".rds"))
  exec.time.lee <- readRDS(paste0("results/", data.name, "_lee_", measure, "_tot_exec_time_v", lee.version, ".rds"))

  treated.bounds.dvds <- readRDS(paste0("results/", data.name, "_dvds_bounds_", fast, measure, "_mc_treated_v", dvds.version, ".rds"))
  control.bounds.dvds <- readRDS(paste0("results/", data.name, "_dvds_bounds_", fast, measure, "_mc_control_v", dvds.version, ".rds"))
  exec.time.dvds <- readRDS(paste0("results/", data.name, "_dvds_", fast, measure, "_tot_exec_time_v", dvds.version, ".rds"))
  
  message("Lee et al. (2024): Total execution time = ", as.numeric(exec.time.lee$elapsed, units="hours"), " hours")
  message("DVDS (Estimators I and II): Total execution time = ", as.numeric(exec.time.dvds$elapsed, units="hours"), " hours")
  
  if (bootstrap) {
    treated.bounds.ci.lee <- readRDS(paste0("results/", data.name, "_lee_bounds_ci_", measure, "_mc_treated_v", lee.boot.version, ".rds"))
    control.bounds.ci.lee <- readRDS(paste0("results/", data.name, "_lee_bounds_ci_", measure, "_mc_control_v", lee.boot.version, ".rds"))
    
    treated.bounds.ci.dvds <- readRDS(paste0("results/", data.name, "_dvds_bounds_ci_", fast, measure, "_mc_treated_v", dvds.boot.version, ".rds"))
    control.bounds.ci.dvds <- readRDS(paste0("results/", data.name, "_dvds_bounds_ci_", fast, measure, "_mc_control_v", dvds.boot.version, ".rds"))
    
    n.bounds <- nrow(do.call(rbind.data.frame, treated.bounds.ci.lee))
    treated.lb.ci.lee <- do.call(rbind.data.frame, treated.bounds.ci.lee)[seq(1, n.bounds-1, 2), "lb"]
    treated.ub.ci.lee <- do.call(rbind.data.frame, treated.bounds.ci.lee)[seq(2, n.bounds, 2), "ub"]
    control.lb.ci.lee <- do.call(rbind.data.frame, control.bounds.ci.lee)[seq(1, n.bounds-1, 2), "lb"]
    control.ub.ci.lee <- do.call(rbind.data.frame, control.bounds.ci.lee)[seq(2, n.bounds, 2), "ub"]
    
    treated.lb.ci.dvds <- do.call(rbind.data.frame, treated.bounds.ci.dvds)[seq(1, n.bounds-1, 2), c("lb1", "lb2")]
    treated.ub.ci.dvds <- do.call(rbind.data.frame, treated.bounds.ci.dvds)[seq(2, n.bounds, 2), c("ub1", "ub2")]
    control.lb.ci.dvds <- do.call(rbind.data.frame, control.bounds.ci.dvds)[seq(1, n.bounds-1, 2), c("lb1", "lb2")]
    control.ub.ci.dvds <- do.call(rbind.data.frame, control.bounds.ci.dvds)[seq(2, n.bounds, 2), c("ub1", "ub2")]
  }
  
  if (data.name == "simul") {
    
    # Monte-Carlo dataset index to get the true survival curves
    ind <- 1
    # Retrieve dataset ind
    simul.data.list <- data.list[[ind]]
    T1.distrib <- simul.data.list$T1.distrib.fun
    T0.distrib <- simul.data.list$T0.distrib.fun
    
    unique.times <- sort(unique(treated.bounds.dvds$time))
    if (RMST) {
      # Number of time points for integration by rectangle method
      nb.timepts <- 100
      true.RMST.treated <- rep(NA, length(unique.times))
      true.RMST.control <- rep(NA, length(unique.times))
      
      for (i in 1:length(unique.times)) {
        time.vec <- seq(from=0, to=unique.times[i], length.out=nb.timepts)
        S.mean.treated <- colMeans(sapply(time.vec, T1.distrib))
        S.mean.control <- colMeans(sapply(time.vec, T0.distrib))
        
        S.mean.treated.midpoint <- (S.mean.treated[-length(S.mean.treated)] + S.mean.treated[-1]) / 2
        true.RMST.treated[i] <- sum(S.mean.treated.midpoint * diff(time.vec))
        
        S.mean.control.midpoint <- (S.mean.control[-length(S.mean.control)] + S.mean.control[-1]) / 2
        true.RMST.control[i] <- sum(S.mean.control.midpoint * diff(time.vec))
      }
      
      # Easier for me to replace like this, but it is the RMST, not the survival function
      S.mean.treated <- true.RMST.treated
      S.mean.control <- true.RMST.control
      
    } else {
      S.mean.treated <- colMeans(sapply(unique.times, T1.distrib))
      S.mean.control <- colMeans(sapply(unique.times, T0.distrib))
    }
    
    diff.surv.lb1 <- treated.bounds.dvds$lb1 - control.bounds.dvds$ub1
    diff.surv.ub1 <- treated.bounds.dvds$ub1 - control.bounds.dvds$lb1
    diff.surv.ign1 <- treated.bounds.dvds$ign.vec1 - control.bounds.dvds$ign.vec1
    
    diff.surv.lb2 <- treated.bounds.dvds$lb2 - control.bounds.dvds$ub2
    diff.surv.ub2 <- treated.bounds.dvds$ub2 - control.bounds.dvds$lb2
    diff.surv.ign2 <- treated.bounds.dvds$ign.vec2 - control.bounds.dvds$ign.vec2
    
    lee.diff.surv.lb <- treated.bounds.lee$lb - control.bounds.lee$ub
    lee.diff.surv.ub <- treated.bounds.lee$ub - control.bounds.lee$lb
    lee.diff.surv.ign <- treated.bounds.lee$ign - control.bounds.lee$ign
    
    # Compute the true difference of survival functions
    true.diff <- S.mean.treated - S.mean.control
    
    
    # Plots
    surv.plot.df <- data.frame(treated.lb1=treated.bounds.dvds$lb1, treated.ub1=treated.bounds.dvds$ub1,
                               control.lb1=control.bounds.dvds$lb1, control.ub1=control.bounds.dvds$ub1,
                               treated.ign1=treated.bounds.dvds$ign.vec1, control.ign1=control.bounds.dvds$ign.vec1,
                               treated.lb2=treated.bounds.dvds$lb2, treated.ub2=treated.bounds.dvds$ub2,
                               control.lb2=control.bounds.dvds$lb2, control.ub2=control.bounds.dvds$ub2,
                               treated.ign2=treated.bounds.dvds$ign.vec2, control.ign2=control.bounds.dvds$ign.vec2,
                               treated.lb.lee=treated.bounds.lee$lb, treated.ub.lee=treated.bounds.lee$ub,
                               control.lb.lee=control.bounds.lee$lb, control.ub.lee=control.bounds.lee$ub,
                               treated.ign.lee=treated.bounds.lee$ign, control.ign.lee=control.bounds.lee$ign,
                               S.mean.treated=S.mean.treated,
                               S.mean.control=S.mean.control,
                               time=treated.bounds.dvds$time)
    
    diff.plot.df <- data.frame(diff.surv.lb1=diff.surv.lb1, diff.surv.ub1=diff.surv.ub1,
                               diff.surv.ign1=diff.surv.ign1,
                               diff.surv.lb2=diff.surv.lb2, diff.surv.ub2=diff.surv.ub2,
                               diff.surv.ign2=diff.surv.ign2,
                               lee.diff.surv.lb=lee.diff.surv.lb, lee.diff.surv.ub=lee.diff.surv.ub,
                               lee.diff.surv.ign=lee.diff.surv.ign,
                               true.diff=true.diff,
                               time=treated.bounds.dvds$time)
    
    if (bootstrap) {
      
      if (RMST) {
        # Initialize with null RMST at time 0 for the 12 columns
        init.ci <- c(rep(0, 12), 0)
      } else {
        # Initialize with survival probability 1 at time 0 for the 12 columns
        init.ci <- c(rep(1, 12), 0)
      }
      
      ci.df <- rbind(init.ci,
                     data.frame(treated.lb1.ci=treated.lb.ci.dvds$lb1, treated.ub1.ci=treated.ub.ci.dvds$ub1,
                                control.lb1.ci=control.lb.ci.dvds$lb1, control.ub1.ci=control.ub.ci.dvds$ub1,
                                treated.lb2.ci=treated.lb.ci.dvds$lb2, treated.ub2.ci=treated.ub.ci.dvds$ub2,
                                control.lb2.ci=control.lb.ci.dvds$lb2, control.ub2.ci=control.ub.ci.dvds$ub2,
                                treated.lb.lee.ci=treated.lb.ci.lee, treated.ub.lee.ci=treated.ub.ci.lee,
                                control.lb.lee.ci=control.lb.ci.lee, control.ub.lee.ci=control.ub.ci.lee,
                                time=unique(treated.bounds.dvds$time)[-1]))
      
      diff.surv.lb1.boot <- ci.df$treated.lb1.ci - ci.df$control.ub1.ci
      diff.surv.ub1.boot <- ci.df$treated.ub1.ci - ci.df$control.lb1.ci
      
      diff.surv.lb2.boot <- ci.df$treated.lb2.ci - ci.df$control.ub2.ci
      diff.surv.ub2.boot <- ci.df$treated.ub2.ci - ci.df$control.lb2.ci
      
      lee.diff.surv.lb.boot <- ci.df$treated.lb.lee.ci - ci.df$control.ub.lee.ci
      lee.diff.surv.ub.boot <- ci.df$treated.ub.lee.ci - ci.df$control.lb.lee.ci
      
      diff.ci.df <- data.frame(diff.surv.lb1.boot=diff.surv.lb1.boot, diff.surv.ub1.boot=diff.surv.ub1.boot,
                               diff.surv.lb2.boot=diff.surv.lb2.boot, diff.surv.ub2.boot=diff.surv.ub2.boot,
                               lee.diff.surv.lb.boot=lee.diff.surv.lb.boot,
                               lee.diff.surv.ub.boot=lee.diff.surv.ub.boot,
                               time=ci.df$time)
    }
    
    if (RMST) {
      
      y.legend.treated <- "RMST for the treated and bounds"
      y.legend.control <- "RMST for the control and bounds"
      
      surv.plot.treated <- ggplot(surv.plot.df) +
        geom_line(aes(x=time, y=S.mean.treated, color=true.label), linetype="dashed", key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ign.lee, group=time), width=0.2, color=lee.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.lb.lee, group=time, color=lee.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ub.lee, group=time, color=lee.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ign2, group=time), width=0.2, color=est2.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.lb2, group=time, color=est2.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ub2, group=time, color=est2.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ign1, group=time), width=0.2, color=est1.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.lb1, group=time, color=est1.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ub1, group=time, color=est1.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        labs(x="Time", y=y.legend.treated, color="Legend") +
        theme_linedraw() +
        theme(legend.position=c(0.85, 0.18),
              legend.background=element_rect(linewidth=0.2, linetype="solid",
                                             colour="black")) +
        scale_color_manual(values=colors)
      
      surv.plot.control <- ggplot(surv.plot.df) +
        geom_line(aes(x=time, y=S.mean.control, color=true.label), linetype="dashed", key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ign.lee, group=time), width=0.2, color=lee.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.lb.lee, group=time, color=lee.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ub.lee, group=time, color=lee.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ign2, group=time), width=0.2, color=est2.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.lb2, group=time, color=est2.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ub2, group=time, color=est2.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ign1, group=time), width=0.2, color=est1.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.lb1, group=time, color=est1.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ub1, group=time, color=est1.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        labs(x="Time", y=y.legend.control, color="Legend") +
        theme_linedraw() +
        theme(legend.position=c(0.85, 0.18),
              legend.background=element_rect(linewidth=0.2, linetype="solid",
                                             colour="black")) +
        scale_color_manual(values=colors)
      
    } else {
      
      y.legend.treated <- "Survival function for the treated and bounds"
      y.legend.control <- "Survival function for the control and bounds"
      
      surv.plot.treated <- ggplot(surv.plot.df) +
        geom_line(aes(x=time, y=S.mean.treated, color=true.label), linetype="dashed", key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ign.lee, group=time), width=0.2, color=lee.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.lb.lee, group=time, color=lee.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ub.lee, group=time, color=lee.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ign2, group=time), width=0.2, color=est2.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.lb2, group=time, color=est2.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ub2, group=time, color=est2.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ign1, group=time), width=0.2, color=est1.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.lb1, group=time, color=est1.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=treated.ub1, group=time, color=est1.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        labs(x="Time", y=y.legend.treated, color="Legend") +
        ylim(0, 1) +
        theme_linedraw() +
        theme(legend.position=c(0.15, 0.18),
              legend.background=element_rect(linewidth=0.2, linetype="solid",
                                             colour="black")) +
        scale_color_manual(values=colors)
      
      surv.plot.control <- ggplot(surv.plot.df) +
        geom_line(aes(x=time, y=S.mean.control, color=true.label), linetype="dashed", key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ign.lee, group=time), width=0.2, color=lee.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.lb.lee, group=time, color=lee.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ub.lee, group=time, color=lee.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ign2, group=time), width=0.2, color=est2.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.lb2, group=time, color=est2.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ub2, group=time, color=est2.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ign1, group=time), width=0.2, color=est1.ign.col, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.lb1, group=time, color=est1.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        geom_boxplot(aes(x=time, y=control.ub1, group=time, color=est1.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
        labs(x="Time", y=y.legend.control, color="Legend") +
        ylim(0, 1) +
        theme_linedraw() +
        theme(legend.position=c(0.15, 0.18),
              legend.background=element_rect(linewidth=0.2, linetype="solid",
                                             colour="black")) +
        scale_color_manual(values=colors)
      
    }
    
    surv.plot.treated.file.name <- paste0("figures/", data.name, "_", fast, measure, "_plot_treated_lee_dvds_estim_I_and_II_v", lee.version, "_v", dvds.version, ".pdf")
    pdf(surv.plot.treated.file.name, width=width, height=height)
    print(surv.plot.treated)
    dev.off()
    
    surv.plot.control.file.name <- paste0("figures/", data.name, "_", fast, measure, "_plot_control_lee_dvds_estim_I_and_II_v", lee.version, "_v", dvds.version, ".pdf")
    pdf(surv.plot.control.file.name, width=width, height=height)
    print(surv.plot.control)
    dev.off()
    
    if (RMST) {
      
      y.legend.diff12 <- "Difference of RMST"
      
    } else {
      
      y.legend.diff12 <- "Difference of survival functions"
      
    }
    
    diff.surv.plot12 <- ggplot(diff.plot.df) +
      geom_hline(yintercept=0, color="black", linetype="dotted") +
      geom_line(aes(x=time, y=true.diff, color=true.label), linetype="dashed", key_glyph=draw_key_rect) +
      geom_boxplot(aes(x=time, y=lee.diff.surv.ign, group=time), width=0.2, color=lee.ign.col, alpha=0, key_glyph=draw_key_rect) +
      geom_boxplot(aes(x=time, y=lee.diff.surv.lb, group=time, color=lee.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
      geom_boxplot(aes(x=time, y=lee.diff.surv.ub, group=time, color=lee.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
      geom_boxplot(aes(x=time, y=diff.surv.ign2, group=time), width=0.2, color=est2.ign.col, alpha=0, key_glyph=draw_key_rect) +
      geom_boxplot(aes(x=time, y=diff.surv.lb2, group=time, color=est2.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
      geom_boxplot(aes(x=time, y=diff.surv.ub2, group=time, color=est2.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
      geom_boxplot(aes(x=time, y=diff.surv.ign1, group=time), width=0.2, color=est1.ign.col, alpha=0, key_glyph=draw_key_rect) +
      geom_boxplot(aes(x=time, y=diff.surv.lb1, group=time, color=est1.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
      geom_boxplot(aes(x=time, y=diff.surv.ub1, group=time, color=est1.label), width=0.5, alpha=0, key_glyph=draw_key_rect) +
      labs(x="Time", y=y.legend.diff12, color="Legend") +
      theme_linedraw() +
      theme(legend.position=c(0.15, 0.18),
            legend.background=element_rect(linewidth=0.2, linetype="solid",
                                           colour="black")) +
      scale_color_manual(values=colors)
    
    diff.surv.plot12.file.name <- paste0("figures/", data.name, "_diff_", fast, measure, "_plot_lee_dvds_estim_I_and_II_v", lee.version, "_v", dvds.version, ".pdf")
    pdf(diff.surv.plot12.file.name, width=width, height=height)
    print(diff.surv.plot12)
    dev.off()
    
  } else if (data.name %in% c("rhc", "gbsg")) {
    
    diff.surv.lb1 <- treated.bounds.dvds$lb1 - control.bounds.dvds$ub1
    diff.surv.ub1 <- treated.bounds.dvds$ub1 - control.bounds.dvds$lb1
    diff.surv.ign1 <- treated.bounds.dvds$ign.vec1 - control.bounds.dvds$ign.vec1
    
    diff.surv.lb2 <- treated.bounds.dvds$lb2 - control.bounds.dvds$ub2
    diff.surv.ub2 <- treated.bounds.dvds$ub2 - control.bounds.dvds$lb2
    diff.surv.ign2 <- treated.bounds.dvds$ign.vec2 - control.bounds.dvds$ign.vec2
    
    lee.diff.surv.lb <- treated.bounds.lee$lb - control.bounds.lee$ub
    lee.diff.surv.ub <- treated.bounds.lee$ub - control.bounds.lee$lb
    lee.diff.surv.ign <- treated.bounds.lee$ign - control.bounds.lee$ign
    
    
    # Plots
    surv.plot.df <- data.frame(treated.lb1=treated.bounds.dvds$lb1, treated.ub1=treated.bounds.dvds$ub1,
                               control.lb1=control.bounds.dvds$lb1, control.ub1=control.bounds.dvds$ub1,
                               treated.ign1=treated.bounds.dvds$ign.vec1, control.ign1=control.bounds.dvds$ign.vec1,
                               treated.lb2=treated.bounds.dvds$lb2, treated.ub2=treated.bounds.dvds$ub2,
                               control.lb2=control.bounds.dvds$lb2, control.ub2=control.bounds.dvds$ub2,
                               treated.ign2=treated.bounds.dvds$ign.vec2, control.ign2=control.bounds.dvds$ign.vec2,
                               treated.lb.lee=treated.bounds.lee$lb, treated.ub.lee=treated.bounds.lee$ub,
                               control.lb.lee=control.bounds.lee$lb, control.ub.lee=control.bounds.lee$ub,
                               treated.ign.lee=treated.bounds.lee$ign, control.ign.lee=control.bounds.lee$ign,
                               time=treated.bounds.dvds$time)
    
    diff.plot.df <- data.frame(diff.surv.lb1=diff.surv.lb1, diff.surv.ub1=diff.surv.ub1,
                               diff.surv.ign1=diff.surv.ign1,
                               diff.surv.lb2=diff.surv.lb2, diff.surv.ub2=diff.surv.ub2,
                               diff.surv.ign2=diff.surv.ign2,
                               lee.diff.surv.lb=lee.diff.surv.lb, lee.diff.surv.ub=lee.diff.surv.ub,
                               lee.diff.surv.ign=lee.diff.surv.ign,
                               time=treated.bounds.dvds$time)
    
    if (bootstrap) {

      if (RMST) {
        # Initialize with null RMST at time 0 for the 12 columns
        init.ci <- c(rep(0, 12), 0)
      } else {
        # Initialize with survival probability 1 at time 0 for the 12 columns
        init.ci <- c(rep(1, 12), 0)
      }
      
      ci.df <- rbind(init.ci,
                     data.frame(treated.lb1.ci=treated.lb.ci.dvds$lb1, treated.ub1.ci=treated.ub.ci.dvds$ub1,
                                control.lb1.ci=control.lb.ci.dvds$lb1, control.ub1.ci=control.ub.ci.dvds$ub1,
                                treated.lb2.ci=treated.lb.ci.dvds$lb2, treated.ub2.ci=treated.ub.ci.dvds$ub2,
                                control.lb2.ci=control.lb.ci.dvds$lb2, control.ub2.ci=control.ub.ci.dvds$ub2,
                                treated.lb.lee.ci=treated.lb.ci.lee, treated.ub.lee.ci=treated.ub.ci.lee,
                                control.lb.lee.ci=control.lb.ci.lee, control.ub.lee.ci=control.ub.ci.lee,
                                time=unique(treated.bounds.dvds$time)[-1]))
      
      diff.surv.lb1.boot <- ci.df$treated.lb1.ci - ci.df$control.ub1.ci
      diff.surv.ub1.boot <- ci.df$treated.ub1.ci - ci.df$control.lb1.ci
      
      diff.surv.lb2.boot <- ci.df$treated.lb2.ci - ci.df$control.ub2.ci
      diff.surv.ub2.boot <- ci.df$treated.ub2.ci - ci.df$control.lb2.ci
      
      lee.diff.surv.lb.boot <- ci.df$treated.lb.lee.ci - ci.df$control.ub.lee.ci
      lee.diff.surv.ub.boot <- ci.df$treated.ub.lee.ci - ci.df$control.lb.lee.ci
      
      diff.ci.df <- data.frame(diff.surv.lb1.boot=diff.surv.lb1.boot, diff.surv.ub1.boot=diff.surv.ub1.boot,
                               diff.surv.lb2.boot=diff.surv.lb2.boot, diff.surv.ub2.boot=diff.surv.ub2.boot,
                               lee.diff.surv.lb.boot=lee.diff.surv.lb.boot,
                               lee.diff.surv.ub.boot=lee.diff.surv.ub.boot,
                               time=ci.df$time)
    }
    
    if (RMST) {
      
      y.legend.treated <- "RMST for the treated and bounds"
      y.legend.control <- "RMST for the control and bounds"
      
      surv.plot.treated <- ggplot(surv.plot.df) +
        geom_line(aes(x=time, y=treated.ign.lee), color=lee.ign.col, linetype="dashed") +
        geom_line(aes(x=time, y=treated.lb.lee), color=lee.bnds.col) +
        geom_line(aes(x=time, y=treated.ub.lee), color=lee.bnds.col) +
        geom_line(aes(x=time, y=treated.ign2), color=est2.ign.col, linetype="dotted") +
        geom_line(aes(x=time, y=treated.lb2), color=est2.bnds.col) +
        geom_line(aes(x=time, y=treated.ub2), color=est2.bnds.col) +
        geom_line(aes(x=time, y=treated.ign1), color=est1.ign.col, linetype="dotted") +
        geom_line(aes(x=time, y=treated.lb1), color=est1.bnds.col) +
        geom_line(aes(x=time, y=treated.ub1), color=est1.bnds.col) +
        geom_point(aes(x=time, y=treated.ign.lee), color=lee.ign.col, shape=1) +
        geom_point(aes(x=time, y=treated.lb.lee), color=lee.bnds.col) +
        geom_point(aes(x=time, y=treated.ub.lee), color=lee.bnds.col) +
        geom_point(aes(x=time, y=treated.ign2), color=est2.ign.col, shape=1) +
        geom_point(aes(x=time, y=treated.lb2), color=est2.bnds.col) +
        geom_point(aes(x=time, y=treated.ub2), color=est2.bnds.col) +
        geom_point(aes(x=time, y=treated.ign1), color=est1.ign.col, shape=1) +
        geom_point(aes(x=time, y=treated.lb1), color=est1.bnds.col) +
        geom_point(aes(x=time, y=treated.ub1), color=est1.bnds.col) +
        labs(x="Time", y=y.legend.treated) +
        theme_linedraw()
      
      surv.plot.control <- ggplot(surv.plot.df) +
        geom_line(aes(x=time, y=control.ign.lee), color=lee.ign.col, linetype="dashed") +
        geom_line(aes(x=time, y=control.lb.lee), color=lee.bnds.col) +
        geom_line(aes(x=time, y=control.ub.lee), color=lee.bnds.col) +
        geom_line(aes(x=time, y=control.ign2), color=est2.ign.col, linetype="dotted") +
        geom_line(aes(x=time, y=control.lb2), color=est2.bnds.col) +
        geom_line(aes(x=time, y=control.ub2), color=est2.bnds.col) +
        geom_line(aes(x=time, y=control.ign1), color=est1.ign.col, linetype="dotted") +
        geom_line(aes(x=time, y=control.lb1), color=est1.bnds.col) +
        geom_line(aes(x=time, y=control.ub1), color=est1.bnds.col) +
        geom_point(aes(x=time, y=control.ign.lee), color=lee.ign.col, shape=1) +
        geom_point(aes(x=time, y=control.lb.lee), color=lee.bnds.col) +
        geom_point(aes(x=time, y=control.ub.lee), color=lee.bnds.col) +
        geom_point(aes(x=time, y=control.ign2), color=est2.ign.col, shape=1) +
        geom_point(aes(x=time, y=control.lb2), color=est2.bnds.col) +
        geom_point(aes(x=time, y=control.ub2), color=est2.bnds.col) +
        geom_point(aes(x=time, y=control.ign1), color=est1.ign.col, shape=1) +
        geom_point(aes(x=time, y=control.lb1), color=est1.bnds.col) +
        geom_point(aes(x=time, y=control.ub1), color=est1.bnds.col) +
        labs(x="Time", y=y.legend.control) +
        theme_linedraw()
      
    } else {
      
      y.legend.treated <- "Survival function for the treated and bounds"
      y.legend.control <- "Survival function for the control and bounds"
      
      surv.plot.treated <- ggplot(surv.plot.df) +
        geom_line(aes(x=time, y=treated.ign.lee), color=lee.ign.col, linetype="dashed") +
        geom_line(aes(x=time, y=treated.lb.lee), color=lee.bnds.col) +
        geom_line(aes(x=time, y=treated.ub.lee), color=lee.bnds.col) +
        geom_line(aes(x=time, y=treated.ign2), color=est2.ign.col, linetype="dotted") +
        geom_line(aes(x=time, y=treated.lb2), color=est2.bnds.col) +
        geom_line(aes(x=time, y=treated.ub2), color=est2.bnds.col) +
        geom_line(aes(x=time, y=treated.ign1), color=est1.ign.col, linetype="dotted") +
        geom_line(aes(x=time, y=treated.lb1), color=est1.bnds.col) +
        geom_line(aes(x=time, y=treated.ub1), color=est1.bnds.col) +
        geom_point(aes(x=time, y=treated.ign.lee), color=lee.ign.col, shape=1) +
        geom_point(aes(x=time, y=treated.lb.lee), color=lee.bnds.col) +
        geom_point(aes(x=time, y=treated.ub.lee), color=lee.bnds.col) +
        geom_point(aes(x=time, y=treated.ign2), color=est2.ign.col, shape=1) +
        geom_point(aes(x=time, y=treated.lb2), color=est2.bnds.col) +
        geom_point(aes(x=time, y=treated.ub2), color=est2.bnds.col) +
        geom_point(aes(x=time, y=treated.ign1), color=est1.ign.col, shape=1) +
        geom_point(aes(x=time, y=treated.lb1), color=est1.bnds.col) +
        geom_point(aes(x=time, y=treated.ub1), color=est1.bnds.col) +
        labs(x="Time", y=y.legend.treated) +
        ylim(0, 1) +
        theme_linedraw()
      
      surv.plot.control <- ggplot(surv.plot.df) +
        geom_line(aes(x=time, y=control.ign.lee), color=lee.ign.col, linetype="dashed") +
        geom_line(aes(x=time, y=control.lb.lee), color=lee.bnds.col) +
        geom_line(aes(x=time, y=control.ub.lee), color=lee.bnds.col) +
        geom_line(aes(x=time, y=control.ign2), color=est2.ign.col, linetype="dotted") +
        geom_line(aes(x=time, y=control.lb2), color=est2.bnds.col) +
        geom_line(aes(x=time, y=control.ub2), color=est2.bnds.col) +
        geom_line(aes(x=time, y=control.ign1), color=est1.ign.col, linetype="dotted") +
        geom_line(aes(x=time, y=control.lb1), color=est1.bnds.col) +
        geom_line(aes(x=time, y=control.ub1), color=est1.bnds.col) +
        geom_point(aes(x=time, y=control.ign.lee), color=lee.ign.col, shape=1) +
        geom_point(aes(x=time, y=control.lb.lee), color=lee.bnds.col) +
        geom_point(aes(x=time, y=control.ub.lee), color=lee.bnds.col) +
        geom_point(aes(x=time, y=control.ign2), color=est2.ign.col, shape=1) +
        geom_point(aes(x=time, y=control.lb2), color=est2.bnds.col) +
        geom_point(aes(x=time, y=control.ub2), color=est2.bnds.col) +
        geom_point(aes(x=time, y=control.ign1), color=est1.ign.col, shape=1) +
        geom_point(aes(x=time, y=control.lb1), color=est1.bnds.col) +
        geom_point(aes(x=time, y=control.ub1), color=est1.bnds.col) +
        labs(x="Time", y=y.legend.control) +
        ylim(0, 1) +
        theme_linedraw()
      
    }
    
    surv.plot.treated.file.name <- paste0("figures/", data.name, "_", fast, measure, "_plot_treated_lee_dvds_estim_I_and_II_v", lee.version, "_v", dvds.version, ".pdf")
    pdf(surv.plot.treated.file.name, width=width, height=height)
    print(surv.plot.treated)
    dev.off()
    
    surv.plot.control.file.name <- paste0("figures/", data.name, "_", fast, measure, "_plot_control_lee_dvds_estim_I_and_II_v", lee.version, "_v", dvds.version, ".pdf")
    pdf(surv.plot.control.file.name, width=width, height=height)
    print(surv.plot.control)
    dev.off()
    
    if (RMST) {

      y.legend.diff <- "Difference of RMST"
      
    } else {
      
      y.legend.diff <- "Difference of survival functions"
      
    }
    
    linewidth <- 0.7
    
    if (RMST) {
      if (data.name == "rhc") {
        lower.ylim <- -10
        upper.ylim <- 5
        x.intercept <- 30
      } else if (data.name == "gbsg") {
        lower.ylim <- -300
        upper.ylim <- 550
        x.intercept <- 1826.25
      }
    } else {
      if (data.name == "rhc") {
        lower.ylim <- -0.3
        upper.ylim <- 0.15
        x.intercept <- 30
      } else if (data.name == "gbsg") {
        lower.ylim <- -0.3
        upper.ylim <- 0.6
        x.intercept <- 1826.25
      }
    }
    
    diff.surv.plot12 <- ggplot(diff.plot.df) +
      geom_hline(yintercept=0, color="black", linetype="solid", linewidth=0.5) +
      geom_vline(xintercept=x.intercept, color="black", linetype="dotted", linewidth=0.6) +
      geom_line(aes(x=time, y=lee.diff.surv.ign), color=lee.ign.col, linetype="dashed") +
      geom_line(aes(x=time, y=lee.diff.surv.lb, color=lee.label), linewidth=linewidth) +
      geom_line(aes(x=time, y=lee.diff.surv.ub, color=lee.label), linewidth=linewidth) +
      geom_line(aes(x=time, y=diff.surv.ign2), color=est2.ign.col, linetype="dashed") +
      geom_line(aes(x=time, y=diff.surv.lb2, color=est2.label), linewidth=linewidth) +
      geom_line(aes(x=time, y=diff.surv.ub2, color=est2.label), linewidth=linewidth) +
      geom_line(aes(x=time, y=diff.surv.ign1), color=est1.ign.col, linetype="dashed") +
      geom_line(aes(x=time, y=diff.surv.lb1, color=est1.label), linewidth=linewidth) +
      geom_line(aes(x=time, y=diff.surv.ub1, color=est1.label), linewidth=linewidth) +
      geom_point(aes(x=time, y=lee.diff.surv.ign), color=lee.ign.col, shape=1) +
      geom_point(aes(x=time, y=lee.diff.surv.lb, color=lee.label)) +
      geom_point(aes(x=time, y=lee.diff.surv.ub, color=lee.label)) +
      geom_point(aes(x=time, y=diff.surv.ign2), color=est2.ign.col, shape=1) +
      geom_point(aes(x=time, y=diff.surv.lb2, color=est2.label)) +
      geom_point(aes(x=time, y=diff.surv.ub2, color=est2.label)) +
      geom_point(aes(x=time, y=diff.surv.ign1), color=est1.ign.col, shape=1) +
      geom_point(aes(x=time, y=diff.surv.lb1, color=est1.label)) +
      geom_point(aes(x=time, y=diff.surv.ub1, color=est1.label)) +
      labs(x="Time", y=y.legend.diff, color="Legend") +
      ylim(lower.ylim, upper.ylim) +
      theme_linedraw() +
      theme(legend.position=c(0.15, 0.15),
            legend.background=element_rect(linewidth=0.2, linetype="solid", 
                                           colour="black")) +
      scale_color_manual(values=colors)
    
    diff.surv.plot12.file.name <- paste0("figures/", data.name, "_diff_", fast, measure, "_plot_lee_dvds_estim_I_and_II_v", lee.version, "_v", dvds.version, ".pdf")
    pdf(diff.surv.plot12.file.name, width=width, height=height)
    print(diff.surv.plot12)
    dev.off()
    
    if (bootstrap) {
      
      diff.surv.plot12.ci <- ggplot(diff.ci.df) +
        geom_hline(yintercept=0, color="black", linetype="solid", linewidth=0.5) +
        geom_vline(xintercept=x.intercept, color="black", linetype="dotted", linewidth=0.6) +
        geom_line(aes(x=time, y=lee.diff.surv.ign), data=diff.plot.df, color=lee.ign.col, linetype="dashed") +
        geom_line(aes(x=time, y=lee.diff.surv.lb.boot, color=lee.label), linewidth=linewidth) +
        geom_line(aes(x=time, y=lee.diff.surv.ub.boot, color=lee.label), linewidth=linewidth) +
        geom_line(aes(x=time, y=diff.surv.ign2), data=diff.plot.df, color=est2.ign.col, linetype="dashed") +
        geom_line(aes(x=time, y=diff.surv.lb2.boot, color=est2.label), linewidth=linewidth) +
        geom_line(aes(x=time, y=diff.surv.ub2.boot, color=est2.label), linewidth=linewidth) +
        geom_line(aes(x=time, y=diff.surv.ign1), data=diff.plot.df, color=est1.ign.col, linetype="dashed") +
        geom_line(aes(x=time, y=diff.surv.lb1.boot, color=est1.label), linewidth=linewidth) +
        geom_line(aes(x=time, y=diff.surv.ub1.boot, color=est1.label), linewidth=linewidth) +
        geom_point(aes(x=time, y=lee.diff.surv.ign), data=diff.plot.df, color=lee.ign.col, shape=1) +
        geom_point(aes(x=time, y=lee.diff.surv.lb.boot, color=lee.label)) +
        geom_point(aes(x=time, y=lee.diff.surv.ub.boot, color=lee.label)) +
        geom_point(aes(x=time, y=diff.surv.ign2), data=diff.plot.df, color=est2.ign.col, shape=1) +
        geom_point(aes(x=time, y=diff.surv.lb2.boot, color=est2.label)) +
        geom_point(aes(x=time, y=diff.surv.ub2.boot, color=est2.label)) +
        geom_point(aes(x=time, y=diff.surv.ign1), data=diff.plot.df, color=est1.ign.col, shape=1) +
        geom_point(aes(x=time, y=diff.surv.lb1.boot, color=est1.label)) +
        geom_point(aes(x=time, y=diff.surv.ub1.boot, color=est1.label)) +
        labs(x="Time", y=y.legend.diff, color="Legend") +
        ylim(lower.ylim, upper.ylim) +
        theme_linedraw() +
        theme(legend.position=c(0.15, 0.15),
              legend.background=element_rect(linewidth=0.2, linetype="solid", 
                                             colour="black")) +
        scale_color_manual(values=colors)
      
      diff.surv.plot12.ci.file.name <- paste0("figures/", data.name, "_diff_", fast, measure, "_plot_ci_lee_dvds_estim_I_and_II_v", lee.boot.version, "_v", dvds.boot.version, ".pdf")
      pdf(diff.surv.plot12.ci.file.name, width=width, height=height)
      print(diff.surv.plot12.ci)
      dev.off()
    }
    
  }
  
}
