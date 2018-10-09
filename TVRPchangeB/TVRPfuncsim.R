# LARS -----------------------------------------------------------------------------------
# Lasso from LARS with BIC or GCV as a stopping rule
lasso.fctsim = function(x, y, win, m.type = c("BIC", "GCV")) {
  
  k1            = ifelse(group.nr == 1, 1, sum(n.simcores[1:(group.nr - 1)],1))
  k2            = sum(n.simcores[1:group.nr])
  m             = n.param
  n             = n.obs
  res.norm      = numeric(0)
  coeff.norm    = numeric(0) 
  lambda.fit    = numeric(0)
  act.set       = numeric(0)
  cond.num      = numeric(0)
  max.eigen     = numeric(0)
  
  for (j in k1:k2) {
    
    res.norm1    = numeric(0)
    coeff.norm1  = numeric(0)
    lambda.fit1  = numeric(0)
    act.set1     = numeric(0)
    cond.num1    = numeric(0)
    max.eigen1   = numeric(0)
    
    y.tmp        = y[[j]]
    x.tmp        = x[[j]]
    
    for (i in 1:(n - win + 1)) {
      
      # Normalization of columns of x
      ywin          = y.tmp[i:(i + win - 1)]
      xwin          = x.tmp[i:(i + win - 1), ]
      nwin          = nrow(xwin)
      
      # OLS fit to standardized x and y
      out.ls        = lm(ywin ~ xwin)
      beta.ols      = out.ls$coeff[2:(m + 1)]
      object        = lars(xwin, ywin, type = "lasso", normalize = FALSE)
      sig2f         = summary(out.ls)$sigma^2
      
      # Get min BIC or GCV
      if (m.type == "BIC"){
        bic         = (log(nwin) * object$df) + (log(as.vector(object$RSS)/nwin) * nwin)
        step        = which.min(bic)
      } else {
        gcv         = (as.vector(object$RSS)/(1 - (object$df/nwin))^2)/nwin
        step        = which.min(gcv)
      }
      
      lambda        = c(object$lambda, 0)[step]  # Lambda minimizing BIC
      
      fit           = predict.lars(object, xwin, s = step, type = "fit",  
                                   mode = "step")$fit
      coeff         = predict.lars(object, xwin, s = step, type = "coef", 
                                   mode = "step")$coefficients
      st            = sum(coeff != 0)          # Number of nonzero coefficients
      mse           = sum((ywin - fit)^2)/(nwin - st - 1)
      
      xbtmp         = xwin %*% coeff
      restmp        = ywin - xbtmp
      lambda.tmp    = (t(restmp) %*% xbtmp) / (sqrt(nwin) * sum(abs(coeff)))
      res.norm1     = c(res.norm1, sqrt(sum(restmp^2)))
      coeff.norm1   = c(coeff.norm1, sum(abs(coeff)))
      lambda.fit1   = c(lambda.fit1, lambda)
      act.set1      = c(act.set1, st)
      cond.num1     = c(cond.num1, kappa(xwin))
      max.eigen1    = c(max.eigen1, max(eigen(t(xwin) %*% xwin)$values))  
    }
    
    res.norm        = cbind(res.norm, res.norm1)
    coeff.norm      = cbind(coeff.norm, coeff.norm1)
    lambda.fit      = cbind(lambda.fit, lambda.fit1)
    act.set         = cbind(act.set, act.set1)
    cond.num        = cbind(cond.num, cond.num1)
    max.eigen       = cbind(max.eigen, max.eigen1)
    
    print(j)
  } 
  
  mean.rn           = apply(res.norm, 1, mean)
  mean.cn           = apply(coeff.norm, 1, mean)
  mean.lb           = apply(lambda.fit, 1, mean)
  mean.as           = apply(act.set, 1, mean)
  mean.k            = apply(cond.num, 1, mean)
  mean.eigen        = apply(max.eigen, 1, mean)
  
  values            = list(lambda.fit, act.set, res.norm, coeff.norm, cond.num, 
                           max.eigen, mean.rn, mean.cn, mean.lb,mean.as, mean.k,
                           mean.eigen)
  names(values)     = c("lambda.fit", "act.set", "res.norm", "coeff.norm", 
                        "cond.num", "max.eigen", "mean.rn", "mean.cn", "mean.lb", 
                        "mean.as", "mean.k", "mean.eigen")
  return(values)
}

# Function to use lasso.fctsim parallelly
par.lassosim = function(gr.nr){ 
  
  group.nr  <<- gr.nr   # Define group of simulations to be computed on the core
  out_cores = lasso.fctsim(X, Y, w, m.type) 
  return(out_cores)
}
# ----------------------------------------------------------------------------------------

# rRAP -----------------------------------------------------------------------------------
# Lasso from rRAP with BIC stopping rule
rap.fctsim = function(x, y, win, m.type = c("BIC", "GCV")) {
  
  updateFun      = getFromNamespace('update.RAP', 'rRAP')
  
  k1             = ifelse(group.nr == 1, 1, sum(n.simcores[1:(group.nr - 1)],1))
  k2             = sum(n.simcores[1:group.nr])
  m              = n.param
  n              = n.obs
  xbeta          = numeric(0)
  res            = numeric(0)
  res.norm       = numeric(0)
  coeff.norm     = numeric(0) 
  lambda.fit     = numeric(0)
  act.set        = numeric(0)
  
  for (j in k1:k2) {
    
    y.tmp            = y[[j]]
    x.tmp            = x[[j]]
    
    xbeta1           = numeric(0)
    res1             = numeric(0)
    res.norm1        = numeric(0)
    coeff.norm1      = numeric(0)
    act.set1         = numeric(0)
    
    # Select burn-in period for the initial estimation of lambda
    ywin.tmp         = y.tmp[1:win]
    xwin.tmp         = x.tmp[1:win, ]
    nwin.tmp         = nrow(xwin.tmp)
    
    # OLS fit to standardized x and y
    out.ls           = lm(ywin.tmp ~ xwin.tmp)
    beta.ols         = out.ls$coeff[2:(m + 1)]
    object           = lars(xwin.tmp, ywin.tmp, type = "lasso", normalize = FALSE)
    sig2f            = summary(out.ls)$sigma^2
    
    # Get min BIC or GCV
    if (m.type == "BIC"){
      bic         = ((log(nwin.tmp) * object$df) + 
                     (log(as.vector(object$RSS)/nwin.tmp) * nwin.tmp))
      step        = which.min(bic)
    } else {
      gcv         = (as.vector(object$RSS)/(1 - (object$df/nwin.tmp))^2)/nwin.tmp
      step        = which.min(gcv)
    }
    
    # Find initial value for RAP algorithm
    lambda.init   = (c(object$lambda, 0)[step])/(1 * nwin.tmp)  # Lambda minimizing BIC or GCV
    
    # Start RAP
    Rtmp          = RAP(X = matrix(xwin.tmp, nrow = win), y = ywin.tmp, r = 0.95, eps = 0.01, 
                      l0 = lambda.init)
    beta.rap      = Rtmp$beta
    
    for (irap in (win + 1):nrow(x.tmp)){
      Rtmp        = updateFun(Rtmp, Ynew = y.tmp[irap], 
                              Xnew = matrix(x.tmp[irap, ], nrow = 1))
      
      beta.tmp    = Rtmp$beta 
      beta.rap    = rbind(beta.rap, beta.tmp)
      q.rap       = sum(Rtmp$beta != 0)
      
      # Collect results for the fit with minimal BIC
      xbtmp       = xwin.tmp %*% as.vector(Rtmp$beta)
      restmp      = ywin.tmp - xbtmp
      xbeta1      = c(xbeta1, xbtmp)
      res1        = c(res1, restmp)
      res.norm1   = c(res.norm1, sqrt(sum(restmp^2)))
      coeff.norm1 = c(coeff.norm1, sum(abs(beta.tmp)))
      act.set1    = c(act.set1, q.rap)
    }
    
    # Collect results for all simulations
    xbeta         = cbind(xbeta, xbeta1)
    res           = cbind(res, res1)
    res.norm      = cbind(res.norm, res.norm1)
    coeff.norm    = cbind(coeff.norm, coeff.norm1)
    lambda.fit    = cbind(lambda.fit, Rtmp$l1Track)
    act.set       = cbind(act.set, act.set1)
    
    print(j)
  } 
  
  values              = list(lambda.fit, act.set, res.norm, coeff.norm)
  names(values)       = c("lambda.fit", "act.set", "res.norm", "coeff.norm")
  
  return(values)
}

# Function to use rap.fctsim parallelly
par.rapsim = function(gr.nr){
  group.nr  <<- gr.nr
  out_cores = rap.fctsim(X, Y, w.rap, m.type) 
  return(out_cores)
}

# ----------------------------------------------------------------------------------------

# Collect results ------------------------------------------------------------------------
# Function collecting results from all the cores
res.sim = function(n.cores, input){
  
  res.norm   = numeric(0)
  coeff.norm = numeric(0) 
  lambda.fit = numeric(0)
  act.set    = numeric(0)
  cond.num   = numeric(0)
  max.eigen  = numeric(0)
  
  # Collect results from all the cores
  for (i in 1:n.cores){
    res.norm    = cbind(res.norm, input[[i]]$res.norm)
    coeff.norm  = cbind(coeff.norm, input[[i]]$coeff.norm)
    lambda.fit  = cbind(lambda.fit, input[[i]]$lambda.fit)
    act.set     = cbind(act.set, input[[i]]$act.set)
    cond.num    = cbind(cond.num, input[[i]]$cond.num)
    max.eigen   = cbind(max.eigen, input[[i]]$max.eigen)
  }
  
  # Compute means and medians over all simulations
  mean.rn       = apply(res.norm, 1, mean)
  mean.cn       = apply(coeff.norm, 1, mean)
  mean.lb.fit   = apply(lambda.fit, 1, mean) 
  mean.as       = apply(act.set, 1, mean)
  mean.k        = apply(cond.num, 1, mean)
  mean.eigen    = apply(max.eigen, 1, mean)
  
  med.rn        = apply(res.norm, 1, median)
  med.cn        = apply(coeff.norm, 1, median)
  med.lb.fit    = apply(lambda.fit, 1, median) 
  med.as        = apply(act.set, 1, median)
  med.k         = apply(cond.num, 1, median)
  med.eigen     = apply(max.eigen, 1, median)
  
  values        = list(res.norm, coeff.norm, lambda.fit, act.set, cond.num,
                       max.eigen, mean.rn, mean.cn, mean.lb.fit, mean.as,
                       mean.k, mean.eigen, med.rn, med.cn, med.lb.fit,  
                       med.as, med.k, med.eigen)
  
  names(values) = c("res.norm", "coeff.norm", "lambda.fit", "act.set", "cond.num",
                    "max.eigen", "mean.rn", "mean.cn", "mean.lb.fit", "mean.as",
                    "mean.k", "mean.eigen", "med.rn", "med.cn", "med.lb.fit",  
                    "med.as", "med.k", "med.eigen")
  
  return(values)  
}
# ----------------------------------------------------------------------------------------

# Scenarios ------------------------------------------------------------------------------
# Function to simulate true coefficients (beta) with a change of q after i = n.cp
beta_sim = function(q.start, q.end){
  
  b      = list()
  
  tmp1   = rep(1, q.start)
  tmp2   = rep(0, n.param - length(tmp1))
  b[[1]] = c(tmp1, tmp2)   # Coefficients for i <= n.cp
  
  tmp3   = rep(1, q.end)
  tmp4   = rep(0, n.param - length(tmp3))
  b[[2]] = c(tmp3, tmp4)   # Coefficients for i > n.cp
  
  return(b)
}

# Function to simulate design matrix (X) with a change of rho after i = n.cp
design_sim = function(r.start, r.end){
  
  mu     = rep(0, n.param) # Mean is set to be 0
  Sigma1 = matrix(0, nrow = n.param, ncol = n.param)   # Covariance matrix for i <= n.cp
  Sigma2 = matrix(0, nrow = n.param, ncol = n.param)   # Covariance matrix for i > n.cp
  
  for (i in 1:n.param) { 
    for (j in 1:n.param) {
      if (i == j){
        Sigma1[i, j] = 1
      } else {
        Sigma1[i, j] = r.start^abs(i - j)
      }
    }
  }
  
  for (i in 1:n.param) { 
    for (j in 1:n.param) {
      if (i == j){
        Sigma2[i, j] = 1
      } else {
        Sigma2[i, j] = r.end^abs(i - j)
      }
    }
  }
  
  X = list()
  set.seed(seed1)
  for (i in 1:n.sim){
    X1     = mvrnorm(n = n.cp, mu, Sigma1)             # Design matrix for i <= n.cp
    X2     = mvrnorm(n = (n.obs - n.cp), mu, Sigma2)   # Design matrix for i > n.cp
    X[[i]] = rbind(X1, X2)
  } 
  
  return(X)  
}

# Function to simulate error term with a change of variance after i = n.cp
eps_sim = function(sd.start, sd.end) {  
  
  eps1  = list()   # Error term for i <= n.cp
  set.seed(seed2)
  for (i in 1:n.sim){
    eps1[[i]] = rnorm(n.cp, mean = 0, sd = sd.start)
  } 
  
  eps2  = list()   # Error term for i > n.cp
  set.seed(seed2)
  for (i in 1:n.sim){
    eps2[[i]] = rnorm((n.obs - n.cp), mean = 0, sd = sd.end)
  }
  
  eps   = list()
  for (i in 1:n.sim){
    eps[[i]] = c(eps1[[i]],eps2[[i]])
  }
  
  return(eps)  
}

# Function to compute all n.obs observations of Y
Y_sim = function(sd.start, sd.end, q.start, q.end, r.start, r.end){
  
  Y     = list()
  eps   = eps_sim(sd.start, sd.end)       # Simulate n.sim vectors of error term
  b     = beta_sim(q.start, q.end)        # Simulate n.sim vectors of beta
  X     <<- design_sim(r.start, r.end)    # Simulate n.sim designs (define X globally)
  
  Y1    = list()                          # Y for i <= n.cp
  for (i in 1:n.sim){
    Y.tmp1 = numeric(0)
    for (j in 1:n.cp){
      Y.tmp1 = c(Y.tmp1, b[[1]] %*% X[[i]][j, ] + eps[[i]][j])
    }
    Y1[[i]] = Y.tmp1
  }
  
  Y2    = list()                          # Y for i > n.cp
  for (i in 1:n.sim){
    Y.tmp2 = numeric(0)
    for (j in (n.cp + 1):n.obs){
      Y.tmp2 = c(Y.tmp2, b[[2]] %*% X[[i]][j, ] + eps[[i]][j])
    }
    Y2[[i]] = Y.tmp2
  }
  
  Y.tmp3 = list()
  for (i in 1:n.sim){
    Y.tmp3[[i]] = c(Y1[[i]],Y2[[i]])
  }
  
  Y = Y.tmp3                              # n.sim vectors of Y of length n.obs
  return(Y)
}

# Function to normalize data to [0, 1] interval
norm.0to1 = function(dat){
  new.dat = (dat - min(dat))/(max(dat) - min(dat))
  new.dat
}
# ----------------------------------------------------------------------------------------