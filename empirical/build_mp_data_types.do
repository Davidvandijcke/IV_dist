* =============================================================================
* Build Type-Specific Groups + Composition Data for Decomposition
* =============================================================================
*
* Extends build_mp_data.do to produce:
*   1. mp_type_groups.dta     — group QFs at county × trimester × race × type level
*   2. mp_type_composition.dta — type shares within each original group
*   3. mp_type_density.dta     — individual bweight + type for kernel density estimation
*
* Types = mom_age_group (< 24 / >= 24) × legitimacy (legit / illegit) = 4 types
*
* Usage:
*   cd IV_dist
*   /Applications/Stata/StataSE.app/Contents/MacOS/stataSE -b do empirical/build_mp_data_types.do
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

di "============================================="
di "Building type-specific datasets"
di "============================================="


* =============================================================================
* PHASE 1: Read natality, merge crosswalks, create types — year by year
* =============================================================================

tempfile master
local first = 1

forvalues yr = 1968/1977 {

    di ""
    di "--- Processing `yr' ---"

    use "`natdir'/natl`yr'.dta", clear
    keep stateres cntyres cityres dbirwt csex dmage mrace birmon dlegit dlivord datayear

    rename dbirwt   bweight
    rename csex     sex
    rename dmage    mom_age
    rename birmon   month
    rename dlegit   legit
    rename dlivord  nchildltot

    gen year = `yr'
    drop datayear

    * Clean
    replace bweight = . if bweight == 9999 | bweight == 0
    replace mrace = . if mrace == 9
    replace mrace = 9 if mrace != 1 & mrace != 2 & mrace != .
    drop if mrace == 9 | mrace == .
    drop if bweight == .

    * County codes
    drop if cntyres == "ZZZ" | cntyres == "99999" | cntyres == ""
    drop if stateres == "" | stateres == "ZZ"
    gen nstate = real(stateres)
    gen ncnty  = real(substr(cntyres, 3, 3))
    gen ncity  = real(cityres)
    drop if nstate == . | ncnty == .
    replace ncity = 999 if ncity == .
    drop stateres cntyres cityres

    * Crosswalk merge
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

    * Trimester
    gen trimester = ceil(month / 3)

    * --- Create type variable ---
    * Type = mom_age_group × legitimacy (2 × 2 = 4 types)
    gen mom_young = (mom_age < 24) if mom_age != .
    gen illegit = (legit == 2) if legit != . & legit != 8 & legit != 9
    * Drop observations with missing type info
    drop if mom_young == . | illegit == .

    gen type_id = 1 + mom_young + 2 * illegit
    label define typelbl 1 "Old+Legit" 2 "Young+Legit" 3 "Old+Illegit" 4 "Young+Illegit"
    label values type_id typelbl

    di "  Obs after type creation: " _N

    * Append
    if `first' {
        save `master', replace
        local first = 0
    }
    else {
        append using `master'
        save `master', replace
    }
}

di ""
di "============================================="
di "Total individual observations: " _N
di "============================================="


* =============================================================================
* PHASE 2: Merge food stamps + controls (at individual level, before collapse)
* =============================================================================

* Food stamp dates
sort stfips countyfips
merge m:1 stfips countyfips using "`almdir'/foodstamps.dta", keepusing(fs_month fs_year)
drop if _merge == 2
drop _merge

* FSP treatment
gen time_fs = (fs_year - 1959) * 12 + fs_month
gen time_birth = (year - 1959) * 12 + month
gen fsp = (time_fs <= time_birth - 3) if time_fs != .
replace fsp = 0 if fsp == .
drop time_fs time_birth fs_month fs_year

* Time variable
gen time = (year - 1968) * 4 + trimester

* State-year
gen state_year = stfips * 10000 + year

* County ID
egen county_id = group(stfips countyfips)

di "Individual obs with FSP: " _N


* =============================================================================
* PHASE 3: Save Dataset C — individual bweight + type for density estimation
* =============================================================================

di ""
di "--- Saving density data (sample for KDE) ---"

* Save a random 10% sample of individual data for density estimation
* (full dataset is too large for R)
set seed 42
gen rand = runiform()
preserve
keep if rand < 0.10
keep stfips countyfips year trimester mrace type_id bweight fsp county_id
compress
save "`outdir'/mp_type_density.dta", replace
di "  Density sample: " _N " observations"
restore
drop rand


* =============================================================================
* PHASE 4: Compute type-specific group quantile functions (Dataset A)
* =============================================================================

di ""
di "--- Computing type-specific group quantile functions ---"

* Fine group = county × trimester × year × race × type
egen fine_group_id = group(stfips countyfips year trimester mrace type_id)

* Quantile functions within each fine group
foreach p of local qtiles {
    bysort fine_group_id: egen q_`p' = pctile(bweight), p(`p')
}

* Group size
bysort fine_group_id: gen nbirths_type = _N

* Also compute original (coarse) group size
egen coarse_group_id = group(stfips countyfips year trimester mrace)
bysort coarse_group_id: gen nbirths_total = _N

* Type share within original group
gen type_share = nbirths_type / nbirths_total


* =============================================================================
* PHASE 5: Save Dataset B — type composition per original group
* =============================================================================

di ""
di "--- Saving type composition data ---"

preserve
* Collapse to original group × type level
* Use mean for fsp since it can vary slightly within groups (~1% of cases)
collapse (first) type_share nbirths_type nbirths_total (mean) fsp, ///
    by(stfips countyfips year trimester mrace type_id county_id state_year time)

* For reshape, all non-reshape vars must be constant within i-groups.
* county_id, state_year, time are constant within (stfips countyfips year trimester mrace).
* fsp may vary slightly by type — take group-level mean.
bysort stfips countyfips year trimester mrace: egen fsp_group = mean(fsp)
drop fsp
rename fsp_group fsp

* Reshape wide: one row per original group, columns for each type share
reshape wide type_share nbirths_type, ///
    i(stfips countyfips year trimester mrace county_id state_year time fsp nbirths_total) ///
    j(type_id)

* Fill missing type shares with 0
forvalues t = 1/4 {
    capture replace type_share`t' = 0 if type_share`t' == .
    capture replace nbirths_type`t' = 0 if nbirths_type`t' == .
}

* Merge controls
sort stfips countyfips year
merge m:1 stfips countyfips year using "`almdir'/reistran.dta", nogen keep(1 3)
sort stfips countyfips year
merge m:1 stfips countyfips year using "`almdir'/reisinc.dta", nogen keep(1 3)
sort stfips countyfips
merge m:1 stfips countyfips using "`almdir'/fscbdata_short.dta", nogen keep(1 3)
capture drop inc3k60 rural60 age560 age6560 employagpct60

* 1960 chars × time
capture gen lnpop60 = ln(pop60)
foreach control in black60 urban60 farmlandpct60 lnpop60 {
    capture gen `control'_time = `control' * time
}

compress
save "`outdir'/mp_type_composition.dta", replace
di "  Composition groups: " _N
restore


* =============================================================================
* PHASE 6: Collapse type-specific groups (Dataset A) and save
* =============================================================================

di ""
di "--- Collapsing to type-specific group level ---"

* Collapse to one row per fine group
bysort fine_group_id: keep if _n == 1
keep stfips countyfips year trimester month mrace type_id ///
     nbirths_type nbirths_total type_share ///
     fsp county_id state_year time ///
     q_* fine_group_id coarse_group_id

replace month = (trimester - 1) * 3 + 2

* Drop small groups
di "  Fine groups before filter: " _N
drop if nbirths_type < 25
di "  Fine groups after dropping < 25 births: " _N

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

* Summary
tab mrace type_id
tab type_id if mrace == 2
tab type_id if mrace == 1

save "`outdir'/mp_type_groups.dta", replace

di ""
di "Saved all datasets to `outdir'/"
di "Done."
