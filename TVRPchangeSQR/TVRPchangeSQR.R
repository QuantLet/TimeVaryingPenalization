# ----------------------------------------------------------------------------------------
# Simulations of changes in a linear model affecting LASSO parameter lambda
# ----------------------------------------------------------------------------------------

# Clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# Set working directory
setwd("H:/DNAS")

# Install and load packages
libraries = c("MASS", "lars", "scales", "doParallel", "rRAP")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("TVRPfuncsim.r")

# Simulation setup
n.obs    = 400        # Number of observations
n.param  = 20         # Number of parameters
n.sim    = 100        # Number of simulations
w        = 50         # Length of moving windows
seed1    = 20150206   # Seed to simulate design matrix X
seed2    = 20150602   # Seed to simulate error terms
sd.start = 1          # Standard deviation of error term for i <= n.cp
sd.end   = 1          # Standard deviation of error term for i > n.cp
q.start  = 5          # Number of nonzero parameters for i <= n.cp
q.end    = 5         # Number of nonzero parameters for i > n.cp
r.start  = 0.1        # Correlation coefficient for design for i <= n.cp
r.end    = seq(0.2, 0.9, 0.1)        # Correlation coefficient for design for i > n.cp
m.type   = "GCV"      # Method of choosing lambda (BIC or GCV - leads to overfit)

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
out_lars  = list()
dellam.lars = list()
for (i.sd in 1:length(sd.end)){
  delta.lam.q   = list()
  for (i.q in 1:length(q.end)){
    delta.lam.r = numeric(0)
    for (i.r in 1:length(r.end)){
      Y         = Y_sim(sd.start, sd.end[i.sd], q.start, q.end[i.q], r.start, r.end[i.r]) 
      out_tmp   = foreach(ic = 1:n.parallel, .packages = c('lars')) %dopar% par.lassosim(ic)
      out_lars  = res.sim(n.parallel, out_tmp) 
      delta.lam.r[i.r] = (mean(out_lars$mean.lb.fit[n.cp:(n.obs - w)])/
                            mean(out_lars$mean.lb.fit[1:(n.cp - w - 1)]))
    }
    delta.lam.q[[i.q]] = delta.lam.r
  }
  dellam.lars[[i.sd]]  = delta.lam.q
}
lars.lambda = unlist(dellam.lars)

if(m.type == "BIC"){dellam.bic = lars.lambda
} else {dellam.gcv = lars.lambda}

# rRAP algorithm
w.rap      = 40
out_rap    = list()
dellam.rap = list()
for (i.sd in 1:length(sd.end)){
  delta.lam.q   = list()
  for (i.q in 1:length(q.end)){
    delta.lam.r = numeric(0)
    for (i.r in 1:length(r.end)){
      Y         = Y_sim(sd.start, sd.end[i.sd], q.start, q.end[i.q], r.start, r.end[i.r]) 
      out_tmp   = foreach(ic = 1:n.parallel, .packages = c('lars', 'rRAP')) %dopar% par.rapsim(ic)
      out_rap   = res.sim(n.parallel, out_tmp) 
      delta.lam.r[i.r] = (mean(out_rap$mean.lb.fit[(n.cp + w - w.rap):(n.obs - w.rap)])/
                            mean(out_rap$mean.lb.fit[(w - w.rap + 1):(n.cp - w.rap - 1)]))
    }
    delta.lam.q[[i.q]] = delta.lam.r
  }
  dellam.rap[[i.sd]] = delta.lam.q
}
rap.lambda = unlist(dellam.rap)

# Close cluster
stopCluster(cl)

# Normalize values of lambda to [0, 1] interval
norm.lars = norm.0to1(out_lars$mean.lb.fit)
norm.rap  = norm.0to1(out_rap$mean.lb.fit[-seq(1:(w - w.rap))])

# Plot LARS and RAP lambdas
par(mar = c(3, 5, 1, 1))
plot(norm.lars, type = "l",  col = "black", axes = FALSE, 
     frame = TRUE, cex.main = 1.5, ylab = expression(paste("Average ", lambda)),
     xlim = c(-(w + 10), (n.obs - w + 10)), ylim = c(0, 1.1), 
     cex.axis = 1, cex.lab = 1.2, lwd = 1, xlab = "")
axis(1, at = c(-w, n.cp - w, n.obs - w), labels = c("0", paste(expression("t ="), n.cp), n.obs),
     cex.axis = 1.2)
axis(2, cex.axis = 1.2)
abline(v = (n.cp - w), lty = 3)
lines(norm.rap, col = "blue", lty = 2, lwd = 1.5)
legend("topright", c("BIC", "RAP"), col = c("black", "blue"),
       ncol = 2, cex = 0.9, lwd = 1.5, lty = c(1, 4, 2))

delta.sd  = as.vector(sd.end/sd.start)
delta.q   = as.vector(q.end/q.start)
delta.r   = as.vector(r.end/r.start)

# Plot relative changes of lambda wrt relative changes of sigma
par(mar = c(5, 5, 1, 1))
plot(dellam.bic ~ delta.sd, type = "l", ylim = c(1, 2.4), xlim = c(1, 2),
     xlab = expression(paste(sigma [2]/sigma [1])), ylab = expression(paste(lambda [2]/lambda [1])),
     cex.axis = 1, cex.lab = 1.2, lwd = 2)
lines(dellam.gcv ~ delta.sd, col = "red", lwd = 2, lty = 4)
lines(rap.lambda ~ delta.sd, col = "blue", lwd = 2, lty = 2)
legend("topright", c("BIC", "GCV", "RAP"), col = c("black", "red", "blue"),
       ncol = 3, cex = 0.6, lwd = 1.5, lty = c(1, 4, 2))

# Plot relative changes of lambda wrt relative changes of q
plot(dellam.bic ~ delta.q, type = "l", xlim = c(1, 3), ylim = c(0.2, 1),
     xlab = expression(paste("q" [2]/ "q" [1])), ylab =  expression(paste(lambda [2]/lambda [1])),
     cex.axis = 1, cex.lab = 1.2, lwd = 2)
lines(dellam.gcv ~ delta.q, col = "red", lwd = 2, lty = 4)
lines(rap.lambda ~ delta.q, col = "blue", lwd = 2, lty = 2)
legend("topright", c("BIC", "GCV", "RAP"), col = c("black", "red", "blue"),
       ncol = 3, cex = 0.6, lwd = 1.5, lty = c(1, 4, 2))

# Plot relative changes of lambda wrt relative changes of rho
plot(dellam.bic ~ delta.r, type = "l", ylim = c(0.85, 1.23),
     xlab = expression(paste(rho [2]/rho [1])), ylab = expression(paste(lambda [2]/lambda [1])),
     cex.axis = 1, cex.lab = 1.2, lwd = 2)
lines(dellam.gcv ~ delta.r, col = "red", lwd = 2, lty = 4)
lines(rap.lambda ~ delta.r, col = "blue", lwd = 2, lty = 2)
legend("topright", c("BIC", "GCV", "RAP"), col = c("black", "red", "blue"),
       ncol = 3, cex = 0.6, lwd = 1.5, lty = c(1, 4, 2))


