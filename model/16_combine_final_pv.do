/*
Created:	1 May 2011
Updated:	26 April 2012
Purpose: 	Combine the PV metrics from the final ensemble, including coverage
*/

// setup the program
	timer on 1
	set mem 1g
	set more off
	log using "${base_dir}/logs/16_final_pv.smcl", replace
	display "Combining final predictive validity metrics on node `r(o1)'"

// load in the pv from each test
	foreach t of global tests {
		append using "$temp_dir/15_ensemble_coverage_`t'.dta"
	}

// create composite PV metrics by taking the mean across holdouts
	preserve
		keep if substr(test, 1, 7) == "holdout"
		collapse (mean) *same_direction* *rmse* *coverage*
		quietly generate test = "combined out-of-sample"
		tempfile composite
		save `composite', replace
	restore
	append using `composite'

// save the results
	quietly replace test = "in-sample" if test == "insample"
	quietly generate order = real(substr(test,9,.))
	quietly replace order = -1 if test == "combined out-of-sample"
	quietly replace order = 0 if test == "in-sample"
	sort order
	drop order
	order test rmse_t2 same_direction_t2 coverage_lower_t2 coverage_upper_t2 top_submodel_rmse_t2 top_submodel_same_direction_t2 top_submodel_coverage_lower_t2 top_submodel_coverage_upper_t2
	save "${base_dir}/results/final_pv.dta", replace
	outsheet using "${base_dir}/results/final_pv.csv", comma replace

** calculate coverage metrics for top submodel in a given family as specified by the user

** top spacetime cause fraction model
	if $top_st_cf == 1 {
	// load in the pv from each test
		foreach t of global tests {
			append using "$temp_dir/15_ensemble_coverage_`t'_top_st_cf_model.dta"
		}

	// create composite PV metrics by taking the mean across holdouts
		preserve
			keep if substr(test, 1, 7) == "holdout"
			collapse (mean) *same_direction* *rmse* *coverage*
			quietly generate test = "combined out-of-sample"
			tempfile composite
		save `composite', replace
	restore
	append using `composite'
	// save the results
		quietly replace test = "in-sample" if test == "insample"
		quietly generate order = real(substr(test,9,.))
		quietly replace order = -1 if test == "combined out-of-sample"
		quietly replace order = 0 if test == "in-sample"
		sort order
		drop order
		order test rmse_t2 same_direction_t2 coverage_lower_t2 coverage_upper_t2 top_submodel_rmse_t2 top_submodel_same_direction_t2 top_submodel_coverage_lower_t2 top_submodel_coverage_upper_t2
		save "${base_dir}/results/final_pv_top_st_cf_model.dta", replace
		outsheet using "${base_dir}/results/final_pv_top_st_cf_model.csv", comma replace
	}

** top spacetime rate model
	if $top_st_rate == 1 {
	// load in the pv from each test
		foreach t of global tests {
			append using "$temp_dir/15_ensemble_coverage_`t'_top_st_rate_model.dta"
		}

	// create composite PV metrics by taking the mean across holdouts
		preserve
			keep if substr(test, 1, 7) == "holdout"
			collapse (mean) *same_direction* *rmse* *coverage*
			quietly generate test = "combined out-of-sample"
			tempfile composite
			save `composite', replace
		restore
		append using `composite'

	// save the results
		quietly replace test = "in-sample" if test == "insample"
		quietly generate order = real(substr(test,9,.))
		quietly replace order = -1 if test == "combined out-of-sample"
		quietly replace order = 0 if test == "in-sample"
		sort order
		drop order
		order test rmse_t2 same_direction_t2 coverage_lower_t2 coverage_upper_t2 top_submodel_rmse_t2 top_submodel_same_direction_t2 top_submodel_coverage_lower_t2 top_submodel_coverage_upper_t2
		save "${base_dir}/results/final_pv_top_st_cf_model.dta", replace
		outsheet using "${base_dir}/results/final_pv_top_st_cf_model.csv", comma replace
	}

** top linear cause fraction model
	if $top_lin_cf == 1 {
	// load in the pv from each test
		foreach t of global tests {
			append using "$temp_dir/15_ensemble_coverage_`t'_top_lin_cf_model.dta"
		}

	// create composite PV metrics by taking the mean across holdouts
		preserve
			keep if substr(test, 1, 7) == "holdout"
			collapse (mean) *same_direction* *rmse* *coverage*
			quietly generate test = "combined out-of-sample"
			tempfile composite
			save `composite', replace
		restore
		append using `composite'

	// save the results
		quietly replace test = "in-sample" if test == "insample"
		quietly generate order = real(substr(test,9,.))
		quietly replace order = -1 if test == "combined out-of-sample"
		quietly replace order = 0 if test == "in-sample"
		sort order
		drop order
		order test rmse_t2 same_direction_t2 coverage_lower_t2 coverage_upper_t2 top_submodel_rmse_t2 top_submodel_same_direction_t2 top_submodel_coverage_lower_t2 top_submodel_coverage_upper_t2
		save "${base_dir}/results/final_pv_top_lin_cf_model.dta", replace
		outsheet using "${base_dir}/results/final_pv_top_lin_cf_model.csv", comma replace
	}

** top linear rate model
	if $top_lin_rate == 1 {
	// load in the pv from each test
		foreach t of global tests {
			append using "$temp_dir/15_ensemble_coverage_`t'_top_st_rate_model.dta"
		}

	// create composite PV metrics by taking the mean across holdouts
		preserve
			keep if substr(test, 1, 7) == "holdout"
			collapse (mean) *same_direction* *rmse* *coverage*
			quietly generate test = "combined out-of-sample"
			tempfile composite
			save `composite', replace
		restore
		append using `composite'

	// save the results
		quietly replace test = "in-sample" if test == "insample"
		quietly generate order = real(substr(test,9,.))
		quietly replace order = -1 if test == "combined out-of-sample"
		quietly replace order = 0 if test == "in-sample"
		sort order
		drop order
		order test rmse_t2 same_direction_t2 coverage_lower_t2 coverage_upper_t2 top_submodel_rmse_t2 top_submodel_same_direction_t2 top_submodel_coverage_lower_t2 top_submodel_coverage_upper_t2
		save "${base_dir}/results/final_pv_top_st_rate_model.dta", replace
		outsheet using "${base_dir}/results/final_pv_top_st_rate_model.csv", comma replace
	}

// close the logs
	timer off 1
	timer list 1
	log close
