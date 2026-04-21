library(ggplot2)
library(foreach)

source("src/utils.R")
source("src/dvds_fun.R")
source("src/lee_zsb_fun.R")


# Set random seed for reproducibility
seed <- 1
set.seed(seed)

# Cluster creation for parallel computation
cl <- initCluster()


#####################
#  Data generation  #
#####################

source("src/mc_data_generation.R")


###################################
#  Experiments on simulated data  #
###################################

data.name <- "simul"
nb.folds <- 5
num.trees <- 500
cq.type <- "linear"


##############################################################################
#                          Exp. 1: Survival function                         #
##############################################################################

set.seed(seed)

exp.nb <- 1
data.version <- 1
gamma <- 3
RMST <- FALSE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, gamma=gamma, RMST=RMST)


##############################################################################
#                               Exp. 2: RMST                                 #
##############################################################################

set.seed(seed)

exp.nb <- 2
data.version <- 1
gamma <- 3
RMST <- TRUE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, gamma=gamma, RMST=RMST, fast.RMST=TRUE)


#############################
#  Experiments on RHC data  #
#############################

data.name <- "rhc"


##############################################################################
#                   Exp. 3: Survival function, Gamma = 1.15                  #
##############################################################################

set.seed(seed)

exp.nb <- 3
gamma <- 1.15
RMST <- FALSE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, gamma=gamma, RMST=RMST)


##############################################################################
#                   Exp. 4: Survival function, Gamma = 1.2                   #
##############################################################################

set.seed(seed)

exp.nb <- 4
gamma <- 1.2
RMST <- FALSE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, gamma=gamma, RMST=RMST)


##############################################################################
#                   Exp. 5: Survival function, Gamma = 1.5                   #
##############################################################################

set.seed(seed)

exp.nb <- 5
gamma <- 1.5
RMST <- FALSE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, gamma=gamma, RMST=RMST)


##############################################################################
#                         Exp. 6: RMST, Gamma = 1.15                         #
##############################################################################

set.seed(seed)

exp.nb <- 6
gamma <- 1.15
RMST <- TRUE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, gamma=gamma, RMST=RMST)


##############################################################################
#                          Exp. 7: RMST, Gamma = 1.2                         #
##############################################################################

set.seed(seed)

exp.nb <- 7
gamma <- 1.2
RMST <- TRUE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, gamma=gamma, RMST=RMST)


##############################################################################
#                          Exp. 8: RMST, Gamma = 1.5                         #
##############################################################################

set.seed(seed)

exp.nb <- 8
gamma <- 1.5
RMST <- TRUE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, gamma=gamma, RMST=RMST)


#############################
#  Experiments on GBSG data  #
#############################

data.name <- "gbsg"


##############################################################################
#                   Exp. 9: Survival function, Gamma = 1.15                  #
##############################################################################

set.seed(seed)

exp.nb <- 9
gamma <- 1.15
RMST <- FALSE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                   Exp. 10: Survival function, Gamma = 1.2                  #
##############################################################################

set.seed(seed)

exp.nb <- 10
gamma <- 1.2
RMST <- FALSE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                   Exp. 11: Survival function, Gamma = 1.5                  #
##############################################################################

set.seed(seed)

exp.nb <- 11
gamma <- 1.5
RMST <- FALSE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                         Exp. 12: RMST, Gamma = 1.15                        #
##############################################################################

set.seed(seed)

exp.nb <- 12
gamma <- 1.15
RMST <- TRUE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                         Exp. 13: RMST, Gamma = 1.2                         #
##############################################################################

set.seed(seed)

exp.nb <- 13
gamma <- 1.2
RMST <- TRUE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                         Exp. 14: RMST, Gamma = 1.5                         #
##############################################################################

set.seed(seed)

exp.nb <- 14
gamma <- 1.5
RMST <- TRUE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                 Exp. 15: Survival function, Gamma = 1.15, CI               #
##############################################################################

set.seed(seed)

exp.nb <- 15
gamma <- 1.15
RMST <- FALSE

bootstrap <- TRUE
B <- 200
alpha <- 0.05

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type, bootstrap=bootstrap, B=B, alpha=alpha)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST, bootstrap=bootstrap,
         B=B, alpha=alpha)

# Save plots
plotFun(dvds.version=9, lee.version=9, # To plot the curves under ignorability
        dvds.boot.version=exp.nb, lee.boot.version=exp.nb,
        data.name=data.name, data.version=data.version,
        RMST=RMST, bootstrap=bootstrap)

##############################################################################
#                 Exp. 16: Survival function, Gamma = 1.08, CI               #
##############################################################################

set.seed(seed)

exp.nb <- 16
gamma <- 1.08
RMST <- FALSE

bootstrap <- TRUE
B <- 200
alpha <- 0.05

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type, bootstrap=bootstrap, B=B, alpha=alpha)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST, bootstrap=bootstrap,
         B=B, alpha=alpha)

# Save plots
plotFun(dvds.version=9, lee.version=9, # To plot the curves under ignorability
        dvds.boot.version=exp.nb, lee.boot.version=exp.nb,
        data.name=data.name, data.version=data.version,
        RMST=RMST, bootstrap=bootstrap)

##############################################################################
#                   Exp. 17: Survival function, Gamma = 1.08                 #
##############################################################################

set.seed(seed)

exp.nb <- 17
gamma <- 1.08
RMST <- FALSE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          cq.type=cq.type)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)

##############################################################################
#                              Exp. 18: fast RMST                            #
##############################################################################

set.seed(seed)

data.name <- "simul"

exp.nb <- 18
data.version <- 1
gamma <- 3
RMST <- TRUE

# DVDS bounds for fast RMST
dvdsMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST,
          fast.RMST=TRUE, cq.type=cq.type)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=2, data.name=data.name,
        data.version=data.version, gamma=gamma, RMST=RMST, fast.RMST=TRUE)


# Stop the cluster for parallel computing
closeCluster(cl)
