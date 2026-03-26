* =============================================================================
* Save individual-level data for mdqr (black mothers only, to save memory)
* =============================================================================

clear all
set more off
set maxvar 10000

local projroot "."
local natdir   "`projroot'/data/in/nber_natality"
local almdir   "`projroot'/data/in/almond_replication/Data-Programs-All"
local outdir   "`projroot'/data/out"

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
    * Keep only black mothers
    keep if mrace == 2
    drop if bweight == .
    drop if mom_age == . | mom_age == 0

    * Individual covariates
    gen mom_age_sq = mom_age^2
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

    * Group ID
    egen group_id = group(stfips countyfips year trimester)

    if `first' {
        save `master', replace
        local first = 0
    }
    else {
        append using `master'
        save `master', replace
    }
}

* Re-create group_id consistently across all years
drop group_id
egen group_id = group(stfips countyfips year trimester)

di "Black mothers: " _N
bysort group_id: gen nbirths = _N

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

* Keep essential variables only
keep bweight sex mom_age mom_age_sq legit_bin ///
     fsp stfips countyfips year trimester ///
     time state_year county_id group_id nbirths ///
     tranpcret tranpcmed tranpcvet ripc ///
     black60_time urban60_time farmlandpct60_time lnpop60_time

compress
save "`outdir'/mp_individual_black.dta", replace
di "Saved: `outdir'/mp_individual_black.dta"
di "Obs: " _N
