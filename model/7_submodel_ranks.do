/*
Created:	1 May 2011
Updated:	26 April 2012
Purpose: 	Combine all the PV tests to calculate submodel ranks
*/

// setup the program
	timer on 1
	set mem 1g
	log using "${base_dir}/logs/7_submodel_ranks.smcl", replace
	display "Calculating submodel weights on node `r(o1)'"

if "$psi_min" == "" global psi_min = 1
if "$psi_int" == "" global psi_int = 0.01
global psi_7_int = $psi_int/2

// load in the error from each test
	foreach t of global tests {
		append using "$temp_dir/6_submodel_pv_`t'.dta"
	}


// create a composite weight by using the median PV tests from the holdouts to create more stable rankings
	preserve
		keep if substr(test, 1, 7) == "holdout"
		collapse (median) same_direction* rmse_t1 rmse_t2, by(submodel)
		quietly {
			egen rmse_rank_t1 = rank(rmse_t1), track
			egen rmse_rank_t2 = rank(rmse_t2), track
			if "$trend_method" == "RMSE" {
				egen trend_rank_t1 = rank(same_direction_t1), track
				egen trend_rank_t2 = rank(same_direction_t2), track
			}
			else {
				egen trend_rank_t1 = rank(same_direction_t1), field
				egen trend_rank_t2 = rank(same_direction_t2), field
			}


			generate tmp = trend_rank_t1 + rmse_rank_t1
			egen total_rank_t1 = rank(tmp), unique
			drop tmp
			generate tmp = trend_rank_t2 + rmse_rank_t2
			egen total_rank_t2 = rank(tmp), unique
			drop tmp

			// output ranks
			save "$base_dir/results/submodel_median_ranks.dta", replace

		}
		generate test = "combined out-of-sample"
		forvalues p = $psi_min($psi_7_int)$psi_max {
			local ps = subinstr(string(`p'), ".", "x", .)
			quietly generate double points_psi_`ps' = `p'^(_N - total_rank_t1)
			summarize points_psi_`ps', meanonly
			quietly generate double weight_psi_`ps' = points_psi_`ps' / `r(sum)'
		}

// convert to numbers of draws from each submodel
		forvalues p = $psi_min($psi_7_int)$psi_max {
			local ps = subinstr(string(`p'), ".", "x", .)
			quietly generate double draws_psi_`ps' = round(1000 * weight_psi_`ps') if test == "combined out-of-sample"
		}
		tempfile composite
		save `composite', replace
	restore
	append using `composite'

// add in submodel info
	quietly generate spacetime_or_linear = "spacetime" if substr(submodel,1,9) == "spacetime"
	quietly replace spacetime_or_linear = "linear" if spacetime_or_linear == ""
	quietly generate submodel_number = real(substr(submodel, length(spacetime_or_linear)+2, .))
	drop submodel
	rename submodel_number submodel
	quietly generate type = ""
	quietly generate name = ""
	quietly generate covariates = ""
	quietly generate dependent_variable = ""
	forvalues i = 1 / $number_submodels {
		mata st_local("type_`i'", st_global("type_`i'"))
		quietly replace type = "`type_`i''" if submodel == `i'
		mata st_local("name_`i'", st_global("name_`i'"))
		quietly replace name = "`name_`i''" if submodel == `i'
		mata st_local("dv_`i'", st_global("dv_`i'"))
		quietly replace dependent_variable = "`dv_`i''" if submodel == `i'
		mata st_local("covariates_`i'", st_global("covariates_`i'"))
		quietly replace covariates = "`covariates_`i''" if submodel == `i'
	}
	drop points* *weight*
	order submodel name spacetime_or_linear dependent_variable type covariates test rmse_t* rmse_rank* same_direction* trend_rank* total_rank* draws*
	quietly replace test = "in-sample" if test == "insample"
	quietly generate order = real(substr(test,9,.))
	quietly replace order = -1 if test == "combined out-of-sample"
	quietly replace order = 0 if test == "in-sample"
	sort order test submodel spacetime_or_linear
	drop order

// save the weights
	save "$temp_dir/7_submodel_ranks.dta", replace
	save "${base_dir}/results/submodel_ranks.dta", replace
	outsheet using "${base_dir}/results/submodel_ranks.csv", comma replace

// close the logs
	timer off 1
	timer list 1
	log close
