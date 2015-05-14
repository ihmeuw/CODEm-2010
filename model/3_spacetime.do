/*
Created:	1 May 2011
Updated:	18 August 2011
Purpose: 	Run the spatiotemporal portion of the model
*/

// setup the program
	timer on 1
	set mem 2g
	set more off
	set matsize 11000
	set maxvar 32000
	log using "${base_dir}/logs/3_spacetime_${super_region}_${super_region_chunk}_${test}.smcl", replace
	display "Running spacetime for super region $super_region (part $super_region_chunk) on node `r(o1)'"

// load in the linear model results
	use "$temp_dir/1_linear_region_1_${test}_${super_region}.dta", clear
	forvalues i = 2 / $number_submodel_chunks {
		merge 1:1 n using "$temp_dir/1_linear_region_`i'_${test}_${super_region}.dta", nogen
	}

    qui count
    local numobs `r(N)'
    local inname "$temp_dir/3_drop_${super_region}_${super_region_chunk}_${test}.csv"
    local outname "$temp_dir/3temp_${super_region}_${super_region_chunk}_${test}.csv"
    format ln_rate residual_* linear_* %12.0g
    outsheet using `inname', comma replace
    local fast_location "${code_dir}/model/fast_spacetime"
    di "`fast_location' `numobs' $number_submodels $start_age $end_age $start_year $end_data_year $lambda $lambda_no_data $zeta $zeta_no_data $omega $test `inname' `outname'"
    !`fast_location' $number_submodels $start_age $end_age $start_year $end_data_year $lambda $lambda_no_data $zeta $zeta_no_data $omega $test `inname' `outname'
    tempfile current_data
    save `current_data'

    capture insheet using `outname', comma names clear
    merge 1:1 _n using `current_data'

    // if no data
    count if has_data == 1
    if `r(N)' == 0 {
        capture replace simple_spacetime_ln_rate ="0"
        capture destring simple_spacetime_ln_rate, replace
        di in red _rc
        capture replace simple_spacetime_ln_rate = .
        di in red _rc
    }

// generate the final spacetime predictions
	forvalues i = 1 / $number_submodels {
		quietly generate spacetime_`i' = linear_`i' + pred_residual_`i'
	}

// put a ceiling on the spacetime rates to never exceed the envelope
	forvalues i = 1 / $number_submodels {
		local dv: word `i' of $dv_list
		if "`dv'" == "ln_rate" quietly replace spacetime_`i' = ln(envelope/pop) if spacetime_`i' > ln(envelope/pop) & spacetime_`i' != .
	}

// keep just the relevant isos
	generate keepme = 0
	foreach c of global isos {
		quietly replace keepme = 1 if iso3 == "`c'"
	}
	keep if keepme == 1

// save the results
    if "$good_data_c_years_path" != "" {
        keep n iso3 year age lt_* ln_* super_region national region test_* spacetime_* envelope pop simple_spacetime_ln_rate predictme_* good_data
    }
    else {
	keep n iso3 year age lt_* ln_* super_region national region test_* spacetime_* envelope pop simple_spacetime_ln_rate predictme_*
    }
	save "$temp_dir/3_spacetime_${super_region}_${super_region_chunk}_${test}.dta", replace

// close the logs
	timer off 1
	timer list 1
	log close
