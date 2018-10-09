# Clear all variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# Set working directory
setwd("")

# Install and load packages
libraries = c("MASS", "scales", "foreach", "doParallel", "lars", "rRAP")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("TVRPfuncapp.r")

# Initiate cluster for parallel computing
n.cores = detectCores()   # Number of cores to be used
cl      = makeCluster(n.cores)
registerDoParallel(cl)
getDoParWorkers()

# Load fMRI data
load("EMO.RData")
n.objects   = 8
data.fear   = read.table("EMO_fear_LR.txt", header = FALSE)
data.neut   = read.table("EMO_neut_LR.txt", header = FALSE)

# Computation setup for LARS
n.nodes     = 15                # Number of brain regions
w           = 30                # Length of moving windows
w.rap       = 20                # Burn-in period for RAP algorithm
m.type      = "BIC"             # Method of choosing lambda (BIC or GCV)
r.factor    = 0.95              # Forgetting factor for RAP
s.size      = 0.025             # Step-size parameter for RAP
lambda.init = 0.3               # Inital lambda value for RAP

# Define number of nodes to regress in every core
n.appcores = rep((n.nodes %/% n.cores), n.cores)
h.appcores = n.nodes %% n.cores
if (h.appcores != 0){
  n.appcores[1:h.appcores] = n.appcores[1:h.appcores] + 1
}
n.parallel = sum(n.appcores != 0)

# LARS: Lasso estimation with moving windows of length w 
out_lars = list()
for (idat in 1:n.objects){
  Sys.time()
  data      = EMO.data[[idat]]
  out_tmp   = foreach(icor = 1:n.parallel,            # Parallel computing
                      .packages = c("lars")) %dopar% par.lassoapp(icor) 
  Sys.time()
  out_lars[[idat]] = res.app(n.parallel, out_tmp)
  print(idat)
}

lambda.lars = numeric(0)
for (k3 in 1:n.objects){
  lambda.lars = cbind(lambda.lars, apply(out_lars[[k3]]$lambda.fit, 1, mean))
}
mean.lars = apply(lambda.lars, 1, mean)

# RAP: Lasso estimation
out_rap = list()
for (idat in 1:n.objects){
  Sys.time()
  data      = EMO.data[[idat]]
  out_tmp   = foreach(icor = 1:n.parallel,            # Parallel computing
                      .packages = c("lars", "rRAP")) %dopar% par.rapapp(icor)   
  Sys.time()
  out_rap[[idat]] = res.app(n.parallel, out_tmp)
  print(idat)
}
lambda.rap      = numeric(0)
for (k3 in 1:n.objects){
  lambda.rap = cbind(lambda.rap, apply(out_rap[[k3]]$lambda.fit, 1, mean))
}
mean.rap = apply(lambda.rap, 1, mean)

# Close cluster
stopCluster(cl)

# Normalize values of lambda to [0, 1] interval
norm.lars = norm.0to1(mean.lars)
norm.rap  = norm.0to1(mean.rap[-seq(1:(w - w.rap))])

# Plot settings
par(mfrow = c(1, 1))
par(mar = c(5, 5, 1, 1))

# Plot time series of lambda both for LARS and RAP
plot(norm.lars, type = "l",  col = "black",
     xlab = "Time", frame = TRUE, cex.main = 1.5, 
     ylab = expression(paste("Average ", lambda)), cex.lab = 1.2, lwd = 1)
lines(norm.rap, col = "blue", lty = 2, lwd = 1.5)
polygon(c((48 - w), max(0, (23 - w)), max(0, (23 - w)), (48 - w)), 
        c(0, 0, 1, 1), col = alpha("red2", 0.2),
        border = FALSE)
polygon(c((106 - w), (81 - w), (81 - w), (106 - w)), 
        c(0, 0, 1, 1), col = alpha("red2", 0.2),
        border = FALSE)
polygon(c((165 - w), (140 - w), (140 - w), (165 - w)), 
        c(0, 0, 1, 1), col = alpha("red2", 0.2),
        border = FALSE)
polygon(c((78 - w), (51 - w), (51 - w), (78 - w)), 
        c(0, 0, 1, 1), col = alpha("blue", 0.2),
        border = FALSE)
polygon(c((136 - w), (111 - w), (111 - w), (136 - w)), 
        c(0, 0, 1, 1), col = alpha("blue", 0.2),
        border = FALSE)
polygon(c((175 - w), (169 - w), (169 - w), (175 - w)), 
        c(0, 0, 1, 1), col = alpha("blue", 0.2),
        border = FALSE)
legend("bottomright", c("BIC", "RAP"), col = c("black", "blue"),
       ncol = 2, cex = 0.9, lwd = 1.5, lty = c(1, 4, 2))

