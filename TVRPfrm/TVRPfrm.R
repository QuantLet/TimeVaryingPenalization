# Clear all variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# Set working directory
# setwd("/Users/Lenka/Documents/IRTG 1792/TVP/TVRPfrm")

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

# Computation setup
n.firm      = 100    # Number of companies to include in computation (maximum is 200)
w           = 63     # Length of moving windows
w.rap       = 53     # Burn-in period for RAP algorithm
m.type      = "BIC"  # Method of choosing lambda (BIC or GCV - leads to overfit)
r.factor    = 0.95   # Forgetting factor for RAP
s.size      = 0.025  # Step-size parameter for RAP
lambda.init = 0.3    # Inital lambda value for RAP

# Load data: 100 companies
tmpdata = read.csv("100_firms_returns_and_scaled_macro_2018-08-13.csv",sep=",") # FRM data
data    = subset(tmpdata, select = c(2:(n.firm + 1)))
dates   = tmpdata[, 1]

# Define number of companies to regress in every core
n.appcores = rep((n.firm %/% n.cores), n.cores)
h.appcores = n.firm %% n.cores
if (h.appcores != 0){
  n.appcores[1:h.appcores] = n.appcores[1:h.appcores] + 1
}
n.parallel = sum(n.appcores != 0)

# LARS: Lasso estimation with moving windows of length w 
Sys.time()
out_tmp    = foreach(i = 1:n.parallel, .packages = c("lars")) %dopar% par.lassoapp(i)  
Sys.time()
out_lars   = res.app(n.parallel, out_tmp)      # Collect results from the cores

# RAP: Lasso estimation 
Sys.time()
out_tmp   = foreach(i = 1:n.parallel, .packages = c("lars", "rRAP")) %dopar% par.rapapp(i)   
Sys.time()
out_rap   = res.app(n.parallel, out_tmp)       # Collect results from the cores

# Close cluster
stopCluster(cl)

# Normalize values of lambda to [0, 1] interval
norm.lars = norm.0to1(out_lars$mean.lb.fit)
norm.rap  = norm.0to1(out_rap$mean.lb.fit[-seq(1:(w - w.rap))])

# Plot settings
par(mfrow = c(1, 1))
par(mar = c(5, 5, 1, 1))
at.tmp = c(grep("2008", dates)[1] - w, grep("2009", dates)[1] - w, 
           grep("2010", dates)[1] - w, grep("2011", dates)[1] - w, 
           grep("2012", dates)[1] - w, grep("2013", dates)[1] - w, 
           grep("2014", dates)[1] - w, grep("2015", dates)[1] - w, 
           grep("2016", dates)[1] - w, grep("2017", dates)[1] - w,
           grep("2018", dates)[1] - w)

# Plot time series of lambda both for LARS and RAP
plot(norm.lars, type = "l",  col = "black", axes = FALSE, 
     xlab = "Year", frame = TRUE, cex.main = 1.5, 
     ylab = expression(paste("Average ", lambda)), cex.lab = 1.2, lwd = 1)
axis(1, cex.axis = 1, labels = c(2008:2018), at = at.tmp)
axis(2, cex.axis = 1)
lines(norm.rap, col = "blue", lty = 2, lwd = 1.5)
legend("topright", c("BIC", "RAP"), col = c("black", "blue"),
       ncol = 2, cex = 0.9, lwd = 1.5, lty = c(1, 4, 2))

# Plot daily close stock prices 
tmpdata.stock = read.csv("100_firms_stocks.csv",sep=",") # FRM stock prices
col.data      = colorRampPalette(c("darkblue", "darkblue", "red2"))(n.firm)
stock.tmp     = tmpdata.stock[, -1]
stock.new     = stock.tmp
stock.order   = order(stock.tmp[1, ])
for (i in 1:n.firm){
  k              = stock.order[i]
  stock.new[, i] = stock.tmp[, k]
}

# All of the stock prices for 100 companies
plot(stock.new[, 1], type = "l",  col = alpha(col.data[1], 0.5), axes = FALSE, 
     xlab = "Year", frame = TRUE, cex.main = 1.5, ylab = "Stock price in USD",
     cex.lab = 1.5, ylim = c(0, 650))
axis(1, cex.axis = 1.5, labels = c(2008:2018), at = at.tmp)
axis(2, cex.axis = 1.5)
for (i in 2:100) {
  lines(stock.new[, i], col = alpha(col.data[i], 0.5))
}
lines(stock.new[, 100], col = alpha(col.data[100], 0.5))


