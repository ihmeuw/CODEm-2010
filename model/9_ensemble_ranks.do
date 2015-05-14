/*
Created:	1 May 2011
Updated:	11 April 2012
Purpose: 	combine all the PV tests from the ensemble models to determine the optimal value of psi
*/

// setup the program
	timer on 1
	clear
	set mem 1g
	set more off
	log using "${base_dir}/logs/9_ensemble_ranks.smcl", replace
	display "Combining PV metrics for ensemble models on node `r(o1)'"

// load in the error from each test
	foreach t of global tests {
		append using "$temp_dir/8_ensemble_pv_`t'.dta"
	}

// create a composite weight by taking the sum of the points from the 3 holdouts
	preserve
		keep if substr(test, 1, 7) == "holdout"
		collapse (sum) *rank* (median) same_direction* rmse_t1, by(psi)
		quietly egen tmp = rank(total_rank_t1), track
		drop total_rank_t1
		rename tmp total_rank_t1

		// output trend ranks
		save "${base_dir}/results/pool_median_ranks.dta", replace

		quietly generate test = "combined out-of-sample"
		tempfile composite
		save `composite', replace
	restore
	append using `composite'
	quietly generate model = "ensemble, psi=" + string(psi)

// add in results from the top submodel
	quietly {
		preserve
		use "${base_dir}/results/submodel_ranks.dta", clear
		replace name = name + " " + spacetime_or_linear
		levelsof name if total_rank_t1 == 1 & test == "combined out-of-sample", l(sm) c
		keep if name == "`sm'"
		keep test rmse_t1 same_direction_t1
		generate model = "top submodel"
		tempfile sub
		save `sub', replace
		restore
		append using `sub'
	}

// make a pretty table of the results
	quietly {
		replace test = "in-sample" if test == "insample"
		generate order = real(substr(test,9,.))
		replace order = -1 if test == "combined out-of-sample"
		replace order = 0 if test == "in-sample"
		sort order total_rank_t1
		keep test model rmse_t1 same_direction_t1 psi total_rank_t1
		order model test psi total_rank_t1 rmse_t1 same_direction_t1
		save "$temp_dir/9_ensemble_pv.dta", replace
		save "${base_dir}/results/ensemble_pv.dta", replace
		outsheet using "${base_dir}/results/ensemble_pv.csv", replace comma
	}

// close the logs
	timer off 1
	timer list 1
	log close

