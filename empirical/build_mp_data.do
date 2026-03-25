* =============================================================================
* Build MP (2025) Analysis Dataset for D-IV
* =============================================================================
*
* Reads NBER natality .dta files (1968-1977), merges with Almond et al. (2011)
* auxiliary data, and produces a group-level dataset with within-group quantile
* functions of birth weight.
*
* Groups = county × trimester × year × race (same as MP Section 6)
*
* Input:
*   - data/in/nber_natality/natl{year}.dta        (NBER pre-built natality)
*   - data/in/almond_replication/Data-Programs-All/  (Almond auxiliary files)
*
* Output:
*   - data/out/mp_analysis_data.dta               (group-level QFs + controls)
*
* Usage:
*   cd IV_dist
*   /Applications/Stata/StataSE.app/Contents/MacOS/stataSE -b do empirical/build_mp_data.do
*
* =============================================================================

clear all
set more off
set maxvar 10000

* --- Paths ---
local projroot "."
local natdir   "`projroot'/data/in/nber_natality"
local almdir   "`projroot'/data/in/almond_replication/Data-Programs-All"
local outdir   "`projroot'/data/out"

* --- Quantile grid (matching D-IV: 19 points from 0.05 to 0.95) ---
local qtiles "5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95"

di "============================================="
di "Building MP analysis dataset"
di "============================================="


* =============================================================================
* PHASE 1: Read and clean natality microdata, year by year
* =============================================================================

tempfile master
local first = 1

forvalues yr = 1968/1977 {

    di ""
    di "--- Processing `yr' ---"

    * Load NBER natality file
    use "`natdir'/natl`yr'.dta", clear

    * Keep only needed variables
    keep stateres cntyres cityres dbirwt csex dmage mrace birmon dlegit dlivord datayear

    * --- Rename to common names ---
    rename dbirwt   bweight
    rename csex     sex
    rename dmage    mom_age
    rename birmon   month
    rename dlegit   legit
    rename dlivord  nchildltot

    * Year: use the loop variable (NBER datayear codes are inconsistent)
    gen year = `yr'
    drop datayear

    * --- Clean variables ---

    * Birth weight: drop missing/invalid
    replace bweight = . if bweight == 9999 | bweight == 0

    * Race: recode to Almond convention (1=white, 2=black, 9=other)
    replace mrace = . if mrace == 9
    replace mrace = 9 if mrace != 1 & mrace != 2 & mrace != .

    * Drop "other" race (following Almond/MP)
    drop if mrace == 9 | mrace == .

    * Drop missing birth weight (can't compute quantile functions)
    drop if bweight == .

    * --- Parse county codes ---
    * NBER stateres = 2-char NCHS state code (string)
    * NBER cntyres  = 5-char NCHS state+county (string, format SSCCC)
    * NBER cityres  = 3-char NCHS city code (string)

    * Drop foreign residents
    drop if cntyres == "ZZZ" | cntyres == "99999" | cntyres == ""
    drop if stateres == "" | stateres == "ZZ"

    * Convert to numeric for crosswalk merge
    gen nstate = real(stateres)
    gen ncnty  = real(substr(cntyres, 3, 3))
    gen ncity  = real(cityres)

    drop if nstate == . | ncnty == .
    replace ncity = 999 if ncity == .

    drop stateres cntyres cityres

    * --- Merge with county FIPS crosswalk ---
    * The crosswalk maps NCHS state/county codes to standard FIPS codes.
    * Year-specific county code columns: ncnty68, ncnty69, etc.
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
    tab _merge
    keep if _merge == 3
    drop _merge nstate ncnty ncity

    * Drop missing county FIPS
    drop if countyfips == . | countyfips == -99
    replace countyfips = . if countyfips == -99

    * --- Create trimester ---
    gen trimester = ceil(month / 3)
    label define trimlbl 1 "Q1 (Jan-Mar)" 2 "Q2 (Apr-Jun)" 3 "Q3 (Jul-Sep)" 4 "Q4 (Oct-Dec)"
    label values trimester trimlbl

    * --- Create group identifier ---
    * Group = county × trimester × year × race (same as MP)
    egen group_id = group(stfips countyfips year trimester mrace)

    di "  Observations: " _N
    di "  Groups: " r(max)

    * --- Compute within-group quantile functions ---
    * For each quantile tau, compute pctile of bweight within each group
    foreach p of local qtiles {
        bysort group_id: egen q_`p' = pctile(bweight), p(`p')
    }

    * --- Count births per group ---
    bysort group_id: gen nbirths = _N

    * --- Collapse to group level ---
    * Keep one row per group with quantile functions + identifiers
    bysort group_id: keep if _n == 1

    keep stfips countyfips year trimester month mrace nbirths group_id q_*

    * Rename month to representative month (first month of trimester)
    replace month = (trimester - 1) * 3 + 2  // middle month of trimester

    di "  Groups after collapse: " _N

    * --- Append to master ---
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
di "Total groups before merges: " _N
di "============================================="


* =============================================================================
* PHASE 2: Merge auxiliary data (food stamps, controls)
* =============================================================================

* --- Merge food stamp start dates ---
di "Merging food stamp dates..."
sort stfips countyfips
merge m:1 stfips countyfips using "`almdir'/foodstamps.dta", ///
    keepusing(fs_month fs_year)
tab _merge
drop if _merge == 2
drop _merge

* --- Compute FSP treatment variable ---
* fsp = 1 if food stamps were in place 3+ months before birth
* time_fs: month when FSP started (months since Jan 1959)
* time_birth: birth month (months since Jan 1959)
gen time_fs = (fs_year - 1959) * 12 + fs_month
gen time_birth = (year - 1959) * 12 + month

* FSP must be in place by 3rd trimester (3 months before birth)
gen fsp = (time_fs <= time_birth - 3) if time_fs != .
replace fsp = 0 if fsp == .

tab fsp
tab year fsp

drop time_fs time_birth fs_month fs_year

* --- Merge REIS transfer data ---
di "Merging REIS transfers..."
sort stfips countyfips year
merge m:1 stfips countyfips year using "`almdir'/reistran.dta"
tab _merge if year >= 1968 & year <= 1977
drop if _merge == 2
drop _merge

* --- Merge REIS income data ---
di "Merging REIS income..."
sort stfips countyfips year
merge m:1 stfips countyfips year using "`almdir'/reisinc.dta"
tab _merge if year >= 1968 & year <= 1977
drop if _merge == 2
drop _merge

* --- Merge 1960 county characteristics ---
di "Merging 1960 county book..."
sort stfips countyfips
merge m:1 stfips countyfips using "`almdir'/fscbdata_short.dta"
tab _merge
drop if _merge == 2
drop _merge

* Drop unnecessary 1960 vars (following Almond's regprep.do)
capture drop inc3k60 rural60 age560 age6560 employagpct60


* =============================================================================
* PHASE 3: Generate fixed effects and controls
* =============================================================================

* --- Drop small groups ---
di "Dropping groups with < 25 births..."
count if nbirths < 25
drop if nbirths < 25
di "Groups remaining: " _N

* --- Time variable (trimester index for FE) ---
gen time = (year - 1968) * 4 + trimester

* --- County identifier for FE ---
gen double stcntyfips = stfips + 0.001 * countyfips

* --- State-year identifier for FE ---
gen state_year = stfips * 10000 + year

* --- Southern states indicator ---
gen south = 0
foreach s in 1 5 10 11 12 13 21 22 24 28 37 40 45 47 48 51 54 {
    replace south = 1 if stfips == `s'
}

* --- Log population ---
capture gen lnpop60 = ln(pop60)

* --- 1960 characteristics × time trend (following MP/Almond) ---
foreach control in black60 urban60 farmlandpct60 lnpop60 {
    capture gen `control'_time = `control' * time
}

* --- Label race ---
label define racelbl 1 "white" 2 "black"
label values mrace racelbl

* --- Unique county ID for clustering ---
egen county_id = group(stfips countyfips)


* =============================================================================
* PHASE 4: Save
* =============================================================================

compress

di ""
di "============================================="
di "Final dataset summary"
di "============================================="

tab mrace
sum nbirths, detail
sum fsp

di "Groups (black): "
count if mrace == 2
di "Groups (white): "
count if mrace == 1

* Save
sort stfips countyfips year trimester mrace
save "`outdir'/mp_analysis_data.dta", replace

di ""
di "Saved: `outdir'/mp_analysis_data.dta"
di "Done."
