library(foreach)
library(ggplot2)

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
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)

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
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


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
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                   Exp. 4: Survival function, Gamma = 1.2                   #
##############################################################################

set.seed(seed)

exp.nb <- 4
gamma <- 1.2
RMST <- FALSE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                   Exp. 5: Survival function, Gamma = 1.5                   #
##############################################################################

set.seed(seed)

exp.nb <- 5
gamma <- 1.5
RMST <- FALSE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                         Exp. 6: RMST, Gamma = 1.15                         #
##############################################################################

set.seed(seed)

exp.nb <- 6
gamma <- 1.15
RMST <- TRUE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                          Exp. 7: RMST, Gamma = 1.2                         #
##############################################################################

set.seed(seed)

exp.nb <- 7
gamma <- 1.2
RMST <- TRUE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


##############################################################################
#                          Exp. 8: RMST, Gamma = 1.5                         #
##############################################################################

set.seed(seed)

exp.nb <- 8
gamma <- 1.5
RMST <- TRUE

# DVDS bounds
dvdsMcExp(version=exp.nb, data.name=data.name,
          gamma=gamma, nb.folds=nb.folds, num.trees=num.trees, RMST=RMST)

# Bounds of Lee et al. (2024)
leeMcExp(version=exp.nb, data.version=data.version, data.name=data.name,
         gamma=gamma, nb.folds=nb.folds, RMST=RMST)

# Save plots
plotFun(dvds.version=exp.nb, lee.version=exp.nb, data.name=data.name,
        data.version=data.version, RMST=RMST)


# Stop the cluster for parallel computing
closeCluster(cl)
