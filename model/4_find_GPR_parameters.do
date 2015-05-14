/*
Created:	1 May 2011
Updated:	18 August 2011
Purpose: 	Find the global MAD estimates from spacetime, to be used as the amplitude/variance in GPR
*/

// setup the program
	timer on 1
	clear
	set mem 2g
	set more off
	set matsize 11000
	set maxvar 32000
	log using "${base_dir}/logs/4_find_gpr_params_${test}.smcl", replace
	display "Compiling GPR parameters on node `r(o1)'"

// load in the spacetime results
	local num_sr: word count $super_regions
	forvalues i = 1/`num_sr' {
		local this_sr: word `i' of $super_regions
		local num_chunks: word `i' of $chunks_per_sr
		forvalues this_chunk = 1/`num_chunks' {
			append using "$temp_dir/3_spacetime_`this_sr'_`this_chunk'_${test}.dta"
		}
	}

// group data by VR only w/ 10+ years vs other
	quietly {
		preserve
		duplicates drop n, force
		merge 1:1 n using "$temp_dir/input_database.dta"
		bysort iso3 age: egen tmp = count(lt_cf) if source_type == "VR" & national == 1
		by iso3 age: egen num_vr = mean(tmp)
		by iso3 age: egen num_tot = count(lt_cf)
		capture drop mad_nsv_type
		generate mad_nsv_type = "VR" if num_vr == num_tot & num_vr >= 10
		generate tmp_national = (lt_cf != . & national == 1)
		by iso3 age: egen has_national = max(tmp_national)
		generate tmp_subnational = (lt_cf != . & national != 1)
		by iso3 age: egen has_subnational = max(tmp_subnational)
		generate has_both = (has_national == 1 & has_subnational == 1)
		replace mad_nsv_type = "subnational" if national != 1 & has_both == 1
		replace mad_nsv_type = "other" if mad_nsv_type == ""
		generate mad_amp_type = mad_nsv_type
		replace mad_amp_type = "other" if mad_amp_type == "subnational"
		keep iso3 age mad_nsv_type mad_amp_type
		duplicates drop iso3 age, force
        tempfile mad
		save `mad'
		restore
		merge m:1 iso3 age using `mad', nogen keep(match master)

	}

// find the MAD estimators
	quietly {
		generate residual_universal_ln_rate = (simple_spacetime_ln_rate - ln_rate) if test_${test} == 0
		generate residual_universal_lt_cf = (logit(exp(simple_spacetime_ln_rate) * (pop / envelope)) - lt_cf) if test_${test} == 0
		bysort age mad_nsv_type: egen mad_nsv_universal_ln_rate = mad(residual_universal_ln_rate)
		bysort age mad_nsv_type: egen mad_nsv_universal_lt_cf = mad(residual_universal_lt_cf)
	}

// loop through each submodel
	forvalues i = 1 / $number_submodels {

// identify the dependent variable for the submodel
		local dv: word `i' of $dv_list

// find the spacetime error
		quietly generate spacetime_error_`i' = spacetime_`i' - `dv' if test_${test} == 0

// find MAD of the residuals for each age group
		bysort age mad_amp_type: egen spacetime_mad_amp_`i' = mad(spacetime_error_`i')
	}

// save results
	duplicates drop n, force
	keep n mad* *_mad_* *_type age
	save "$temp_dir/4_gpr_parameters_${test}.dta", replace

// close the logs
	timer off 1
	timer list 1
	log close
