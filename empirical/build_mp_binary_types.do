* =============================================================================
* Build binary-type groups: Young (<24) vs Old (>=24)
* =============================================================================
* Simpler than 4-way split — most groups have both types present.
* =============================================================================

clear all
set more off
set maxvar 10000

local projroot "."
local natdir   "`projroot'/data/in/nber_natality"
local almdir   "`projroot'/data/in/almond_replication/Data-Programs-All"
local outdir   "`projroot'/data/out"

local qtiles "5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95"

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

    * Binary type: 1 = old (>=24), 2 = young (<24)
    gen byte btype = 1 if mom_age >= 24 & mom_age != .
    replace btype = 2 if mom_age < 24 & mom_age != .
    drop if btype == .

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

* --- Compute type-specific group QFs ---
egen fine_group_id = group(stfips countyfips year trimester mrace btype)
foreach p of local qtiles {
    bysort fine_group_id: egen q_`p' = pctile(bweight), p(`p')
}
bysort fine_group_id: gen nbirths_type = _N

* Original group size
egen coarse_group_id = group(stfips countyfips year trimester mrace)
bysort coarse_group_id: gen nbirths_total = _N
gen type_share = nbirths_type / nbirths_total

* Collapse to type-specific group level
bysort fine_group_id: keep if _n == 1
keep stfips countyfips year trimester month mrace btype ///
     nbirths_type nbirths_total type_share ///
     fsp county_id state_year time q_*

replace month = (trimester - 1) * 3 + 2

* Drop small type-cells
di "Type groups before filter: " _N
drop if nbirths_type < 25
di "Type groups after filter: " _N

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

tab mrace btype
tab btype if mrace == 2

save "`outdir'/mp_binary_type_groups.dta", replace
di "Done."
