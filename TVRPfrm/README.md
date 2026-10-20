<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of QuantLet: TVRPfrm

Published in: submitted to N/A

Description: Performs the Lasso regression with two distinct algorithms. The first one uses moving window method and the Bayesian information criterion (BIC) or the generalized cross-validation (GCV) to calibrate the penalty parameter (lambda), and the second is called the real-time adaptive penalization (RAP). The input data are daily logarithmic stock returns of 100 largest U.S. financial companies listed on NASDAQ. In the fitting procedure, each of the selected companies is taken as a dependent variable and is regressed by means of the Lasso regression on the remaining ones. The code plots the time series of the average penalization parameter, lambda. Visualisation of daily stock prices of the financial companies from the dataset is also available.’

Keywords: Lasso, shrinkage, L1-norm penalty, change point, bic, gcv, adaptive penalization, regression, plot, moving window, time-series, stock return

See also: TVRPchangeSQR, TVRPchangeB, TVRPfmri, XFGTVP_BetaChange, XFGTVP_FRM, XFGTVP_LambdaSim, TVPactiveset, TVPbetanorm, TVPdesign, TVPvariance

Author: Lenka Zboňáková

Submitted: 9 October 2018 by Lenka Zboňáková

Datafile: 100_firms_returns_and_scaled_macro_2018-08-13.csv, 100_firms_stocks.csv

Input: 
- n.firm:      Number of companies (maximum is 100)
- w :          Size of the moving window
- w.rap :      Burn-in period for the RAP algorithm
- m.type:      Type of the method to calibrate lambda (“BIC” or “GCV”)
- r.factor:    Forgetting factor for the RAP
- s.size:      Step-size parameter for the RAP
- lambda.init: Initial lambda value for the RAP

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/TimeVaryingPenalization/master/TVRPfrm/TVRPfrm-1.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/TimeVaryingPenalization/master/TVRPfrm/TVRPfrm-2.png" alt="Image" />
</div>

