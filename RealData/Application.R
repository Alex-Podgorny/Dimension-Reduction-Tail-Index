
# Load all R scripts from the 'Methods' directory
# These scripts define the estimation methods used in the analysis
methods_func <- list.files(path = "Methods", pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
sapply(methods_func, source)

# Load the dataset
load("RealData/Data/Data.RData")

# Normalize the covariates (X) to [0, 1] for each column
X <- apply(DATA[,-1], 2, function(x) (x - min(x)) / (max(x) - min(x)))

# Transform the response variable (y) by scaling and applying an exponential function
y <- c(DATA[,1])
y <- exp(y / sd(y))

# Set basic variables
n <- length(y)  # Number of observations
p <- ncol(X)    # Number of covariates

# Visualize the tail behavior with a QQ-plot
Y_sort <- sort(y)  # Sort the response variable
k <- 100           # Number of largest observations considered
x <- log(k / 1:k)  # X-axis values for the QQ-plot
z <- log(Y_sort[n - 1:k + 1]) - log(Y_sort[n - k])  # Y-axis values
plot(x = x, y = z, xlab = "log(k / i)", ylab = "log differences")
abline(0, (lm(z ~ x + 0)$coefficients), col = "red")  # Add a regression line

# Estimation setup
set.seed(1)  # Set seed for reproducibility

# Build the compact set X0 by selecting observations within a radius
r0 <- quantile(as.numeric(apply(X, 1, function(x) norm(x, "2"))), 0.9)
X0 <- which(as.numeric(apply(X, 1, function(x) norm(x, "2"))) < r0)

# Set hyperparameters for the estimation
alpha <- n^(-0.3)  # Tail level
b <- 0.2           # Bandwidth scaling exponent

# Initialization for CTI subspace estimation
Bhat_CTI_q <- list()  # List to store Bhat estimates for different q values
c_q <- c()            # Vector to store the average tail index for each q

# Estimate the CTI subspace dimension q
for (q in 1:p) {
  h <- n^(-b / q) / 2                       # Bandwidth for dimension q
  n0 <- ceiling(n * h^q * alpha * log(n)^1.1)  # Number of observations for the estimation
  Bhat_CTI <- CTI(X, y, X0[1:n0], q, alpha, h)  # Estimate the CTI subspace
  Bhat_CTI_q[[q]] <- sign(Bhat_CTI[1,1])*Bhat_CTI                  
  c_q[q] <- mean(local_Hill(X, y, X[X0,], Bhat_CTI, alpha, h), na.rm = TRUE)
  
  # Stop if the tail index increases for the current dimension
  if (q > 1 && c_q[q - 1] < c_q[q]) {
    q_hat <- q - 1  # Optimal dimension
    break
  }
}

# Save the results
Results <- list(Bhat_CTI_q = Bhat_CTI_q, c_q = c_q)
save(Results, file = "RealData/Results.RData")

# Estimate the tail index gamma(x) for observations in X0
gamma_hat <- local_Hill(X, y, X[X0,], Bhat_CTI_q[[q_hat]], alpha, n^(-b) / 2, mink = 1)
IndNA <- which(!is.na(gamma_hat))  # Indices of non-NA gamma estimates

# Plot gamma_hat(x) vs. the projection on the estimated CTI subspace
plot(x = c(t(Bhat_CTI_q[[q_hat]]) %*% t(X[X0,]))[IndNA],
     y = na.omit(gamma_hat),
     xlab = "Estimated CTI Subspace",
     ylab = "Estimated Tail Index",
     type = "p",
     pch = 3)

# Plot the largest ozone observations vs. the CTI subspace projection
Z <- c(DATA[,1])               # Original response variable
Z_sort <- sort(Z)              # Sort response variable
large_obs <- Z > Z_sort[n - alpha * n]  # Select largest observations

plot(x = c(t(Bhat_CTI_q[[q_hat]]) %*% t(X))[large_obs],
     y = Z[large_obs],
     xlab = "Estimated CTI Subspace",
     ylab = "Ozone",
     type = "p",
     pch = 1)

# Fit a quadratic model to the largest observations
BX1 <- c(t(Bhat_CTI_q[[q_hat]]) %*% t(X))[large_obs]
Y_log <- Z[large_obs]
lm_model <- lm(Y_log ~ BX1 + I(BX1^2))  # Quadratic regression
fit <- c(cbind(rep(1, length(BX1)), BX1, BX1^2) %*% lm_model$coefficients)

# Add the fitted curve to the plot
lines(BX1[order(fit)], fit[order(fit)], col = "red")
