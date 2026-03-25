* =============================================================================
* Build enhanced density sample with continuous individual covariates
* =============================================================================
* Quick extraction: read 1970-1974 natality (5 years), merge crosswalks,
* keep a 5% random sample with bweight + mom_age + sex + legit + group info.
* =============================================================================

clear all
set more off

local projroot "."
local natdir   "`projroot'/data/in/nber_natality"
local almdir   "`projroot'/data/in/almond_replication/Data-Programs-All"
local outdir   "`projroot'/data/out"

tempfile master
local first = 1

forvalues yr = 1969/1975 {
    di "--- `yr' ---"
    use "`natdir'/natl`yr'.dta", clear
    keep stateres cntyres cityres dbirwt csex dmage mrace birmon dlegit

    rename dbirwt bweight
    rename csex sex
    rename dmage mom_age
    rename birmon month
    rename dlegit legit

    gen year = `yr'

    * Clean
    replace bweight = . if bweight == 9999 | bweight == 0
    replace mrace = . if mrace == 9
    replace mrace = 9 if mrace != 1 & mrace != 2 & mrace != .
    drop if mrace == 9 | mrace == .
    drop if bweight == .
    drop if mom_age == . | mom_age == 0

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

    * County ID
    egen county_id = group(stfips countyfips)

    * 5% random sample
    set seed `yr'
    gen rand = runiform()
    keep if rand < 0.05
    drop rand

    keep stfips countyfips year trimester mrace bweight mom_age sex legit fsp county_id

    if `first' {
        save `master', replace
        local first = 0
    }
    else {
        append using `master'
        save `master', replace
    }
}

compress
save "`outdir'/mp_density_enhanced.dta", replace
di "Enhanced density sample: " _N " obs"
tab mrace
sum mom_age bweight
