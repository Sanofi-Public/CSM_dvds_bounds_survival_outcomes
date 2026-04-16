library(fastDummies)

# Read .csv file
rhc.data <- read.csv("data/RHC/rhc.csv")
# Remove first column (index)
rhc.data <- rhc.data[-1]

# 1 if death and 0 if censored because of discharge or date of last contact (Delta).
Delta <- 1 * !is.na(rhc.data$dthdte)

# Observed time
# If death is observed, date of death - admission date
# If death is not observed, max(date of last contact, discharge date) - admission date
T.obs <- ifelse(Delta == 1,
                rhc.data$dthdte - rhc.data$sadmdte,
                pmax(rhc.data$lstctdte, rhc.data$dschdte) - rhc.data$sadmdte)

# Administrative censoring at 180 days
Delta <- ifelse(T.obs > 180, 0, Delta)
T.obs <- ifelse(T.obs > 180, 180, T.obs)

# Treatment is 1 if RHC is used, 0 otherwise
A <- ifelse(rhc.data$swang1 == "RHC", 1, 0)

summary(rhc.data)

# Create observed covariates columns
X.cont.labels <- c("age", "edu", "aps1", "scoma1", "meanbp1", "wblc1", "hrt1",
                   "resp1", "temp1", "pafi1", "alb1", "hema1", "bili1", "crea1",
                   "sod1", "pot1", "paco21", "ph1", "wtkilo1",
                   "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx",
                   "renalhx", "liverhx", "gibledhx", "malighx", "immunhx",
                   "transhx", "amihx")

X.cat.labels <- c("sex", "ca", "dnr1", "ninsclas", "resp", "card", "neuro",
                  "gastr", "renal", "meta", "hema", "seps", "trauma", "ortho",
                  "race", "income")

X.cont <- rhc.data[, X.cont.labels]
X.cat <- dummy_cols(rhc.data[, X.cat.labels], remove_selected_columns=TRUE, remove_first_dummy=TRUE)
X <- cbind(X.cont, X.cat)
names(X) <- paste0("X.", 1:ncol(X))

# Unname covariates to avoid errors (because it could be renamed X.X.1 instead of X.1)
X <- unname(X)

n.obs <- nrow(rhc.data)

# Check for NAN values
if (sum(is.na(X)) + sum(is.na(T.obs)) + sum(is.na(A)) + sum(is.na(Delta)) != 0) {
  stop("NAN values in the data")
}

# I put the list in a list. Easier after to combine with the code for simulated data.
preprocessed.data.list <- list(list(X=X, U=rep(NA, n.obs), T1=rep(NA, n.obs), T0=rep(NA, n.obs),
                                    T1.distrib.fun=NULL,
                                    T0.distrib.fun=NULL,
                                    T01=rep(NA, n.obs), C=rep(NA, n.obs), T.obs=T.obs, A=A,
                                    Delta=Delta))

saveRDS(preprocessed.data.list, file="./data/RHC/preprocessed_data.rds")
