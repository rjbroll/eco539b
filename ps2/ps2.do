********************** ECO 539b: Problem Set 2 *************************
cd "/Users/ryanbroll/Documents/GitHub/class/eco539b/ps2/"

* import data
use ak91.dta, clear

* keep relevant sample
keep if census == 1980 & cohort == 2
tab census cohort

* (i) wage and education
reg lwage educ, robust
reg lwage educ, cluster(SOB)

* (ii) growing up in NJ
gen nj = (SOB == 34)
reg lwage nj, robust
reg lwage nj, cluster(SOB)
