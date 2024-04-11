**************** ECO 539b: PS3 *********************
* CD
cd "/Users/ryanbroll/Documents/GitHub/class/eco539b/ps3"

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
