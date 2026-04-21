library(fastDummies)

# Read .csv file
gbsg.data <- read.csv("data/GBSG/gbsg.csv")
# Remove first column (index)
gbsg.data <- gbsg.data[-1]

summary(gbsg.data)

# Status
# 1 recurrence or death, 0 alive without recurrence
Delta <- gbsg.data$status

# Observed time
# Recurrence free survival time; days to first recurrence, death or last contact
T.obs <- gbsg.data$rfstime

# Treatment is 1 if tamoxifen is used, 0 otherwise
A <- gbsg.data$hormon

# Create observed covariates columns
X.cont.labels <- c("age", "size", "nodes", "pgr", "er")

X.bin.labels <- c("meno")

X.cat.labels <- c("grade")

X.cont <- gbsg.data[, X.cont.labels]
X.bin <- gbsg.data[, X.bin.labels]
X.cat <- dummy_cols(gbsg.data[, X.cat.labels], remove_selected_columns=TRUE, remove_first_dummy=TRUE)
X <- cbind(X.cont, X.bin, X.cat)
names(X) <- paste0("X.", 1:ncol(X))

# Unname covariates to avoid errors (because it could be renamed X.X.1 instead of X.1)
X <- unname(X)

n.obs <- nrow(gbsg.data)

# I put the list in a list. Easier after to combine with the code for simulated data.
preprocessed.data.list <- list(list(X=X, U=rep(NA, n.obs), T1=rep(NA, n.obs), T0=rep(NA, n.obs),
                                    T1.distrib.fun=NULL,
                                    T0.distrib.fun=NULL,
                                    T01=rep(NA, n.obs), C=rep(NA, n.obs), T.obs=T.obs, A=A,
                                    Delta=Delta))

saveRDS(preprocessed.data.list, file="./data/GBSG/preprocessed_data.rds")
