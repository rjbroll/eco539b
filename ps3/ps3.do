**************** ECO 539b: PS3 *********************
* CD
cd "/Users/ryanbroll/Documents/GitHub/class/eco539b/ps3"


******************* Q1 *********************
* Import data
use famine.dta, clear

* Construct x1 and x2
gen x1 = lgrain_pred * famine
gen x2 = lgrain_pred * (1-famine)

* Regression 1, (a) (i) - (iii)
reg ldeaths x1 x2 ltotpop lurbpop i.year, robust
reg ldeaths x1 x2 ltotpop lurbpop i.year, vce(hc2, dfadjust)
matrix list e(adj_df)

* Regression 1, (b) (i) - (iii)
reg ldeaths x1 x2 ltotpop lurbpop i.year, vce(cluster prov)
reg ldeaths x1 x2 ltotpop lurbpop i.year, vce(hc2 prov, dfadjust)
matrix list e(adj_df)

* Regression 2, (a) (i) - (iii)
reg ldeaths x1 x2 ltotpop lurbpop i.year if inrange(year,1953,1965), robust
reg ldeaths x1 x2 ltotpop lurbpop i.year if inrange(year,1953,1965), vce(hc2, dfadjust)
matrix list e(adj_df)

* Regression 2, (b) (i) - (iii)
reg ldeaths x1 x2 ltotpop lurbpop i.year if inrange(year,1953,1965), vce(cluster prov)
reg ldeaths x1 x2 ltotpop lurbpop i.year if inrange(year,1953,1965), vce(hc2 prov, dfadjust)
matrix list e(adj_df)

* Regression 1 - Construct dataset to export to Matlab for bootstrapping
preserve
keep prov year ldeaths x1 x2 ltotpop lurbpop
order prov year ldeaths x1 x2 ltotpop lurbpop
tab year, g(ydum)
drop ydum1 // omitted year dummy
export delimited reg1.csv, replace
restore

* Regression 2 - Construct dataset to export to Matlab for bootstrapping
keep if inrange(year,1953,1965)
keep prov year ldeaths x1 x2 ltotpop lurbpop
order prov year ldeaths x1 x2 ltotpop lurbpop
tab year, g(ydum)
drop ydum1 // omitted year dummy
export delimited reg2.csv, replace



******************* Q3 - Pacific *********************
* Import data
use ak91.dta, clear

* sample selection
keep if census == 1980 & cohort == 2 & division == 9

* create instrument
gen q1 = (age == floor(age))


* Run IV, first stage, and reduced form regressions to get parameters
	* IV regression
ivregress 2sls lwage (educ = q1), robust
local bhat=_b[educ]                              
local se_bhat=_se[educ] 
predict uhat, resid

	* First stage regression
reg educ q1, robust
predict vhat, resid
local F = (_b[q1]/_se[q1])^2
egen mean_Z = mean(q1)
gen Z_tilde = q1 - mean_Z
gen Zuhat = Z_tilde*uhat
gen Zvhat = Z_tilde*vhat
correlate Zuhat Zvhat
local rho_tilde_hat = r(rho)

	* Reduced form regression
reg lwage q1, robust
predict epshat, resid
gen Zepshat = Z_tilde*epshat

* Compute Anderson-Rubin confidence interval
local q = (invnorm(0.975))^2
local m_arL = `bhat' + ((-(`q'*`rho_tilde_hat'/sqrt(`F')) - sqrt(`q')* sqrt(1 - ((`q'*(1 - `rho_tilde_hat'^2))/`F')))/(1 - (`q'/`F')))*`se_bhat'
local m_arR = `bhat' + ((-(`q'*`rho_tilde_hat'/sqrt(`F')) + sqrt(`q')* sqrt(1 - ((`q'*(1 - `rho_tilde_hat'^2))/`F')))/(1 - (`q'/`F')))*`se_bhat'
display `m_arL' 
display `m_arR'

* Estimate R:= V(Zv)/V(Zeps) and C:= Corr(Zv,Zeps)
corr Zvhat Zepshat, cov
local R = sqrt(r(Var_1)/r(Var_2))
corr Zvhat Zepshat
local rho_rf = r(rho)
local beta = 0
display (`rho_rf' - `beta'*`R')/sqrt(1+(`beta'*`R')^2 - 2*`rho_rf'*`beta'*`R')
local beta = .25
display (`rho_rf' - `beta'*`R')/sqrt(1+(`beta'*`R')^2 - 2*`rho_rf'*`beta'*`R')

******************* Q3 - mid-Atlantic *********************
* Import data
use ak91.dta, clear

* sample selection
keep if census == 1980 & cohort == 2 & division == 2

* create instrument
gen q1 = (age == floor(age))


* Run IV, first stage, and reduced form regressions to get parameters
	* IV regression
ivregress 2sls lwage (educ = q1), robust
local bhat=_b[educ]                              
local se_bhat=_se[educ] 
predict uhat, resid

	* First stage regression
reg educ q1, robust
predict vhat, resid
local F = (_b[q1]/_se[q1])^2
egen mean_Z = mean(q1)
gen Z_tilde = q1 - mean_Z
gen Zuhat = Z_tilde*uhat
gen Zvhat = Z_tilde*vhat
correlate Zuhat Zvhat
local rho_tilde_hat = r(rho)

	* Reduced form regression
reg lwage q1, robust
predict epshat, resid
gen Zepshat = Z_tilde*epshat

* Compute Anderson-Rubin confidence interval
local q = (invnorm(0.975))^2
local m_arL = `bhat' + ((-(`q'*`rho_tilde_hat'/sqrt(`F')) - sqrt(`q')* sqrt(1 - ((`q'*(1 - `rho_tilde_hat'^2))/`F')))/(1 - (`q'/`F')))*`se_bhat'
local m_arR = `bhat' + ((-(`q'*`rho_tilde_hat'/sqrt(`F')) + sqrt(`q')* sqrt(1 - ((`q'*(1 - `rho_tilde_hat'^2))/`F')))/(1 - (`q'/`F')))*`se_bhat'
display `m_arL' 
display `m_arR'

* Estimate R:= V(Zv)/V(Zeps) and C:= Corr(Zv,Zeps)
corr Zvhat Zepshat, cov
local R = sqrt(r(Var_1)/r(Var_2))
corr Zvhat Zepshat
local rho_rf = r(rho)
local beta = 0
display (`rho_rf' - `beta'*`R')/sqrt(1+(`beta'*`R')^2 - 2*`rho_rf'*`beta'*`R')
local beta = .25
display (`rho_rf' - `beta'*`R')/sqrt(1+(`beta'*`R')^2 - 2*`rho_rf'*`beta'*`R')
