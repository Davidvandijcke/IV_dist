* =============================================================================
* Replicate MP First Stage: Within-Group Quantile Regression
* =============================================================================
*
* For each county-trimester-race group j and quantile tau:
*   qreg bweight sex mom_age mom_age_sq legit, quantile(tau)
* Extract the intercept → this is the QR-adjusted group quantile function.
*
* Then save the intercepts as group-level data for the second stage in R.
*
* =============================================================================

clear all
set more off
set maxvar 10000

local projroot "."
local natdir   "`projroot'/data/in/nber_natality"
local almdir   "`projroot'/data/in/almond_replication/Data-Programs-All"
local outdir   "`projroot'/data/out"

local qtiles "5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95"

* =============================================================================
* PHASE 1: Load individual data (same pipeline as before, but keep covariates)
* =============================================================================

tempfile master
local first = 1

forvalues yr = 1968/1977 {
    di "--- `yr' ---"
    use "`natdir'/natl`yr'.dta", clear
    keep stateres cntyres cityres dbirwt csex dmage mrace birmon dlegit

    rename dbirwt bweight
    rename csex sex
    rename dmage mom_age
    rename birmon month
    rename dlegit legit

    gen year = `yr'

    replace bweight = . if bweight == 9999 | bweight == 0
    replace mrace = . if mrace == 9
    replace mrace = 9 if mrace != 1 & mrace != 2 & mrace != .
    drop if mrace == 9 | mrace == .
    drop if bweight == .
    drop if mom_age == . | mom_age == 0

    * Individual covariates for QR (matching MP Section 6)
    gen mom_age_sq = mom_age^2
    * Recode legitimacy: 1=legit, 0=illegit/unknown
    gen legit_bin = (legit == 1) if legit != .
    replace legit_bin = 0 if legit_bin == .

    * County crosswalk
    drop if cntyres == "ZZZ" | cntyres == "99999" | cntyres == ""
    gen nstate = real(stateres)
    gen ncnty = real(substr(cntyres, 3, 3))
    gen ncity = real(cityres)
    drop if nstate == . | ncnty == .
    replace ncity = 999 if ncity == .
    drop stateres cntyres cityres

    local yy = `yr' - 1900
    preserve
    use "`almdir'/nat_cwalk.dta", clear
    keep stfips countyfips nstate ncnty`yy' ncity`yy'
    drop if ncnty`yy' == -99
    rename ncnty`yy' ncnty
    rename ncity`yy' ncity
    duplicates drop
    tempfile cwalk
    save `cwalk'
    restore

    merge m:1 nstate ncnty ncity using `cwalk'
    keep if _merge == 3
    drop _merge nstate ncnty ncity
    drop if countyfips == . | countyfips == -99

    gen trimester = ceil(month / 3)

    * FSP
    sort stfips countyfips
    merge m:1 stfips countyfips using "`almdir'/foodstamps.dta", keepusing(fs_month fs_year)
    drop if _merge == 2
    drop _merge
    gen time_fs = (fs_year - 1959) * 12 + fs_month
    gen time_birth = (year - 1959) * 12 + month
    gen fsp = (time_fs <= time_birth - 3) if time_fs != .
    replace fsp = 0 if fsp == .
    drop time_fs time_birth fs_month fs_year

    gen time = (year - 1968) * 4 + trimester
    gen state_year = stfips * 10000 + year
    egen county_id = group(stfips countyfips)

    if `first' {
        save `master', replace
        local first = 0
    }
    else {
        append using `master'
        save `master', replace
    }
}

di "Total individuals: " _N

* Group identifier
egen group_id = group(stfips countyfips year trimester mrace)

* Group size — drop small groups (MP threshold: K1 + 1 + 25 = 30)
bysort group_id: gen nbirths = _N
drop if nbirths < 30

di "Individuals after dropping small groups: " _N
distinct group_id
di "Groups: " r(ndistinct)


* =============================================================================
* PHASE 2: Run within-group QR for each quantile (using statsby)
* =============================================================================

* Keep only what we need for QR
keep bweight sex mom_age mom_age_sq legit_bin ///
     stfips countyfips year trimester mrace month ///
     fsp time state_year county_id group_id nbirths

* For efficiency, process BLACK mothers only first (smaller sample)
* Comment out this line to process all races:
* keep if mrace == 2

di ""
di "============================================="
di "Running within-group quantile regressions"
di "============================================="

* Create a dataset to store QR intercepts
* For each group, we need the intercept at each quantile level

tempfile indiv_data
save `indiv_data'

* Run statsby for each quantile
foreach tau of local qtiles {
    di ""
    di "--- Quantile tau = `tau' ---"

    use `indiv_data', clear

    * statsby runs qreg for each group and saves coefficients
    * We only need the intercept (_b[_cons])
    statsby qr_intercept_`tau'=_b[_cons], by(group_id) ///
        saving("`outdir'/qr_tau`tau'", replace) ///
        noisily: ///
        qreg bweight sex mom_age mom_age_sq legit_bin, quantile(`tau')
}


* =============================================================================
* PHASE 3: Merge all quantile intercepts into one dataset
* =============================================================================

di ""
di "--- Merging QR intercepts ---"

use "`outdir'/qr_tau5", clear

foreach tau in 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 {
    merge 1:1 group_id using "`outdir'/qr_tau`tau'", nogen
}

* Rename for consistency with D-IV pipeline
foreach tau of local qtiles {
    rename qr_intercept_`tau' qr_`tau'
}

* Merge back group identifiers
merge 1:1 group_id using `indiv_data', keepusing(stfips countyfips year trimester mrace ///
    fsp time state_year county_id nbirths month) keep(3) nogen

* Keep one row per group (statsby already collapsed)
* But we need group-level controls too

* Merge controls
sort stfips countyfips year
merge m:1 stfips countyfips year using "`almdir'/reistran.dta", nogen keep(1 3)
sort stfips countyfips year
merge m:1 stfips countyfips year using "`almdir'/reisinc.dta", nogen keep(1 3)
sort stfips countyfips
merge m:1 stfips countyfips using "`almdir'/fscbdata_short.dta", nogen keep(1 3)
capture drop inc3k60 rural60 age560 age6560 employagpct60
capture gen lnpop60 = ln(pop60)
foreach control in black60 urban60 farmlandpct60 lnpop60 {
    capture gen `control'_time = `control' * time
}

compress

tab mrace
sum nbirths

save "`outdir'/mp_qr_intercepts.dta", replace

* Clean up temporary quantile files
foreach tau of local qtiles {
    capture erase "`outdir'/qr_tau`tau'.dta"
}

di ""
di "Saved: `outdir'/mp_qr_intercepts.dta"
di "Done."
