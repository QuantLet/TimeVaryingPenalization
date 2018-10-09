# LARS -----------------------------------------------------------------------------------
# Lasso from LARS with BIC or GCV as a stopping rule
lasso.fctapp = function(data, win, m.type = c("BIC", "GCV")) {
  
  k1            = ifelse(group.nr == 1, 1, sum(n.appcores[1:(group.nr - 1)],1))
  k2            = sum(n.appcores[1:group.nr])
  m             = ncol(data) - 1
  n             = nrow(data)
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
    
    y.tmp        = as.vector(data[, j])
    x.tmp        = as.matrix(data[, -j])
    
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
par.lassoapp = function(gr.nr){ 
  
  group.nr  <<- gr.nr   # Define group of simulations to be computed on the core
  out_cores = lasso.fctapp(data, w, m.type) 
  return(out_cores)
}
# ----------------------------------------------------------------------------------------

# rRAP -----------------------------------------------------------------------------------
# Lasso from rRAP with BIC stopping rule
rap.fctapp = function(data, win, r.factor, s.size, lambda.init) {
  
  updateFun      = getFromNamespace('update.RAP', 'rRAP')
  
  k1             = ifelse(group.nr == 1, 1, sum(n.appcores[1:(group.nr - 1)],1))
  k2             = sum(n.appcores[1:group.nr])
  m              = ncol(data) - 1
  n              = nrow(data)
  xbeta          = numeric(0)
  res            = numeric(0)
  res.norm       = numeric(0)
  coeff.norm     = numeric(0) 
  lambda.fit     = numeric(0)
  act.set        = numeric(0)
  
  for (j in k1:k2) {
    
    y.tmp            = as.vector(data[, j])
    x.tmp            = as.matrix(data[, -j])
    
    xbeta1           = numeric(0)
    res1             = numeric(0)
    res.norm1        = numeric(0)
    coeff.norm1      = numeric(0)
    act.set1         = numeric(0)
    
    # Select burn-in period for the initial estimation of lambda
    ywin.tmp         = y.tmp[1:win]
    xwin.tmp         = x.tmp[1:win, ]
    nwin.tmp         = nrow(xwin.tmp)
    
    # # OLS fit to standardized x and y
    # out.ls           = lm(ywin.tmp ~ xwin.tmp)
    # beta.ols         = out.ls$coeff[2:(m + 1)]
    # object           = lars(xwin.tmp, ywin.tmp, type = "lasso", normalize = FALSE)
    # sig2f            = summary(out.ls)$sigma^2
    # 
    # # Get min BIC or GCV
    # if (m.type == "BIC"){
    #   bic         = ((log(nwin.tmp) * object$df) + 
    #                  (log(as.vector(object$RSS)/nwin.tmp) * nwin.tmp))
    #   step        = which.min(bic)
    # } else {
    #   gcv         = (as.vector(object$RSS)/(1 - (object$df/nwin.tmp))^2)/nwin.tmp
    #   step        = which.min(gcv)
    # }
    # 
    # # Find initial value for RAP algorithm
    # lambda.init   = (c(object$lambda, 0)[step])/(1 * nwin.tmp)  # Lambda minimizing BIC or GCV
    
    # Start RAP
    Rtmp          = RAP(X = matrix(xwin.tmp, nrow = win), y = ywin.tmp, r = r.factor,
                        eps = s.size, l0 = lambda.init)
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
par.rapapp = function(gr.nr){
  group.nr  <<- gr.nr
  out_cores = rap.fctapp(data, w.rap, r.factor, s.size, lambda.init) 
  return(out_cores)
}
# ----------------------------------------------------------------------------------------

# Collect results ------------------------------------------------------------------------
# Function collecting results from all the cores
res.app = function(n.cores, input){
  
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

# Function to normalize data to [0, 1] interval
norm.0to1 = function(dat){
  new.dat = (dat - min(dat))/(max(dat) - min(dat))
  new.dat
}
# ----------------------------------------------------------------------------------------