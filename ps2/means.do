
* Simulate outcome process
cd "/Users/ryanbroll/Documents/GitHub/research/flexible_did/figures/"

clear all 
set seed 1616


* Set parameters
	* Panel parameters
local N = 10000
local T = 15

	* Outcome process parameters
local rho = .65
local beta = 1

* Set seed and observations
local obs = `N'*`T'
set obs `obs'

* Form panel
egen i = seq(), block(`T')
bysort i : gen t = _n 
xtset i t 


* Generate outcome process: y_{it} = rho * y_{i,t-1} + beta * D_{it} + e_{it}, e_{it} i.i.d. N(0,1)
	* Base process for period before treatment
gen e = rnormal()
gen y = i/(1-`rho')  if t == 1
replace y = i+ `rho'*L.y + e if t > 1 & t <= 9
reg y L.y
by i: gen ytilde = y - y[1]
reg ytilde L.ytilde

	* Selectively assign treatment in the 10th period to high values of y in the 9th period
gen r = runiform() if t == 10
gen D = (r > .2) if (t == 10) & (L.y >= 1) // selection
replace D = (r > .8) if (t == 10) & (L.y < 1) // selection
// gen D = (r>.5) if t == 10 // no selection
drop r
replace D = L.D if t >=11
replace D = 0 if D == .
replace y = `rho'*L.y + D*`beta' + e if t >=10

* Generate dummies and run event study regression
by i: gen evertreated = (D[10] == 1)
gen dynamict = t if evertreated == 1
replace dynamict = 0 if dynamict == .
xtreg y ib9.dynamict i.t, fe vce(robust)

* Plot event study
coefplot, keep(*.dynamict) vertical base ///
rename(1.dynamict="-9" 2.dynamict="-8" 3.dynamict="-7" 4.dynamict="-6" 5.dynamict="-5" ///
	   6.dynamict="-4" 7.dynamict="-3" 8.dynamict="-2" 9.dynamict = "-1" 10.dynamict="0"  ///
	   11.dynamict="1" 12.dynamict="2" 13.dynamict="3" 14.dynamict="4" 15.dynamict="5") ///
	   ytitle("coefficient on treatment") xtitle("event time") title("An Event Study") ///
	   yline(0)

	   
graph export fake_eventstudy.pdf,replace


* estimate correctly
reg y L.y ib9.dynamict , robust

// // * Estimate dynamic panel AR(1) using pooled regression?
// reg y L.y D0
//
//
// reg y D0
