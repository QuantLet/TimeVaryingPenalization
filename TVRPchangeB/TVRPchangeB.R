# ----------------------------------------------------------------------------------------
# Simulations of changes in a linear model affecting LASSO parameter lambda
# ----------------------------------------------------------------------------------------

# Clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# Set working directory
setwd("")

# Install and load packages
libraries = c("MASS", "lars", "scales", "doParallel", "rRAP")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Call functions to fit the model
source("TVRPfuncsimB.r")

# Simulation setup
n.obs    = 400        # Number of observations
n.param  = 20         # Number of parameters
n.sim    = 100        # Number of simulations
w        = 50         # Length of moving windows
w.rap    = 40         # Burn-in period for RAP algorithm
seed1    = 20150206   # Seed to simulate design matrix X
seed2    = 20150602   # Seed to simulate error terms
sd       = 1          # Standard deviation of error term 
r        = 0.1        # Correlation coefficient for design
m.type   = "BIC"      # Method of choosing lambda (BIC or GCV)

# Define observation with a change point
if(n.obs %% 2 == 1) n.obs = n.obs + 1 ;
n.cp  = n.obs/2

# Initiate cluster for parallel computing
n.cores = detectCores()   # Number of cores to be used
cl      = makeCluster(n.cores)
registerDoParallel(cl)
getDoParWorkers()

# Define number of scenarios in every core
n.simcores = rep((n.sim %/% n.cores), n.cores)
h.simcores = n.sim %% n.cores
if (h.simcores != 0){
  n.simcores[1:h.simcores] = n.simcores[1:h.simcores] + 1
}
n.parallel = sum(n.simcores != 0)

# LARS algorithm: Lasso estimation with moving windows of length w 
out_lars    = list()
Y           = Y_sim(sd, r) 
out_tmp     = foreach(ic = 1:n.parallel, 
                      .packages = c('lars')) %dopar% par.lassosim(ic)
out_lars    = res.sim(n.parallel, out_tmp) 
lars.lambda = out_lars$mean.lb.fit

# RAP algorithm
out_rap    = list()
Y          = Y_sim(sd, r) 
out_tmp    = foreach(ic = 1:n.parallel, 
                     .packages = c('lars', 'rRAP')) %dopar% par.rapsim(ic)
out_rap    = res.sim(n.parallel, out_tmp) 
rap.lambda = out_rap$mean.lb.fit[(w - w.rap):(n.obs - w.rap)]

# Close cluster
stopCluster(cl)

# Normalize values of lambda to [0, 1] interval
norm.lars = norm.0to1(lars.lambda)
norm.rap  = norm.0to1(rap.lambda)

# Plot LARS and RAP lambdas together
par(mar = c(3, 5, 1, 1))
plot(norm.lars, type = "l",  col = "black", axes = FALSE, 
     frame = TRUE, cex.main = 1.5, ylab = expression(paste("Average ", lambda)),
     xlim = c(-(w + 10), (n.obs - w + 10)), ylim = c(0, 1.1), 
     cex.axis = 1, cex.lab = 1.2, lwd = 1, xlab = "")
axis(1, at = c(-w, n.cp - w, n.obs - w), 
     labels = c("0", paste(expression("t ="), n.cp), n.obs), cex.axis = 1.2)
axis(2, cex.axis = 1.2)
abline(v = (n.cp - w), lty = 3)
lines(norm.rap, col = "blue", lty = 2, lwd = 1.5)
legend("topright", c("BIC", "RAP"), col = c("black", "blue"),
       ncol = 2, cex = 0.9, lwd = 1.5, lty = c(1, 4, 2))

