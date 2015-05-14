/*
Created:	1 May 2011
Updated:	6 June 2012
Purpose: 	Calculate predictive validity for the ensemble models
*/

// setup the program
	timer on 1
	clear
	set mem 5g
	set maxvar 32000
	set more off
	log using "${base_dir}/logs/8_ensemble_pv_${test}.smcl", replace
	display "Finding pooled model PV for test $test on node `r(o1)'"

** safeguard to prevent errors for people using old runfiles
if "$psi_min" == "" global psi_min = 1
if "$psi_int" == "" global psi_int = 0.01

// load in the submodel weights
	use "$temp_dir/7_submodel_ranks.dta", clear
	forvalues p = $psi_min($psi_int)$psi_max {
		local ps = subinstr(string(`p'), ".", "x", .)
		summarize draws_psi_`ps', meanonly
		local total_draws_`ps' = `r(sum)'
		forvalues i = 1 / $number_submodels {
			quietly levelsof draws_psi_`ps' if submodel == `i' & test == "combined out-of-sample" & spacetime_or_linear == "spacetime", l(s_d_`ps'_`i') c
			if missing("`s_d_`ps'_`i''") local s_d_`ps'_`i' 0
			local s_w_`ps'_`i' = `s_d_`ps'_`i'' / `total_draws_`ps''
			quietly levelsof draws_psi_`ps' if submodel == `i' & test == "combined out-of-sample" & spacetime_or_linear == "linear", l(l_d_`ps'_`i') c
			if missing("`l_d_`ps'_`i''") local l_d_`ps'_`i' 0
			local l_w_`ps'_`i' = `l_d_`ps'_`i'' / `total_draws_`ps''
		}
	}

// load in the pooled results
	local num_sr: word count $super_regions
	forvalues i = 1/`num_sr' {
		local this_sr: word `i' of $super_regions
		local num_chunks: word `i' of $chunks_per_sr
		forvalues this_chunk = 1/`num_chunks' {
			** adding capture for not so that models with very little data do not break because GPR can't run in places with no data
			cap append using "$temp_dir/5_gpr_`this_sr'_`this_chunk'_${test}.dta"
		}
	}
	compress

// merge on the linear predictions
	drop if gpr_1_spacetime_mean == .
	duplicates drop n, force
	forvalues i = 1 / $number_submodel_chunks {
		quietly merge 1:1 n using "$temp_dir/2_linear_country_`i'_${test}.dta", keepusing(linear_*) keep(match master) nogen
	}

// keep just the test data
	if "${test}" == "insample" replace test_insample = 1
	quietly generate true_ln_rate = ln_rate
	drop if test_${test} == 0 | ln_rate == .
	compress

// generate draws of the log rates based on the proportions found above
	forvalues p = $psi_min($psi_int)$psi_max {
		local ps = subinstr(string(`p'), ".", "x", .)
		quietly generate pool_psi_`ps' = 0
		forvalues i = 1 / $number_submodels {
			local dv: word `i' of $dv_list
			if "`dv'" == "ln_rate" quietly replace pool_psi_`ps' = pool_psi_`ps' + (`l_w_`ps'_`i'' * linear_`i') + (`s_w_`ps'_`i'' * gpr_`i'_spacetime_mean)
			else if "`dv'" == "lt_cf" quietly replace pool_psi_`ps' = pool_psi_`ps' + (`l_w_`ps'_`i'' * ln(invlogit(linear_`i') * (envelope / pop))) + (`s_w_`ps'_`i'' * ln(invlogit(gpr_`i'_spacetime_mean) * (envelope / pop)))
		}
	}

** find first differences in predictions/data
	keep iso3 year age true_ln_rate pool_* ln_rate_sd test_${test}
	preserve
	collapse (mean) pool_* true_ln_rate, by(iso3 year age test_${test})
	egen ca = group(iso3 age)
	tsset ca year
if "$trend_window" == "" global trend_window 1

	forvalues j = $trend_window_min / $trend_window {
			generate data_`j'_fd = true_ln_rate - L`j'.true_ln_rate
	}

	forvalues p = $psi_min($psi_int)$psi_max {
		forvalues j = $trend_window_min / ${trend_window} {

			local ps = subinstr(string(`p'), ".", "x", .)
			generate pool_psi_`ps'_`j'_fd = pool_psi_`ps' - L`j'.pool_psi_`ps'
		}
	}



** count how often the first difference of the prediction lines up with the data (no eta)
	forvalues p = $psi_min($psi_int)$psi_max {
		local ps = subinstr(string(`p'), ".", "x", .)

		** quietly {
			generate same_direction_t1_psi_`ps' = 0
			generate same_direction_t2_psi_`ps' = 0

			local denominator_pool_t1_`ps' 0
			local denominator_pool_t2_`ps' 0
		** }

		forvalues j = $trend_window_min / ${trend_window} {
			** quietly {
				if "$trend_method" == "RMSE" {
					replace same_direction_t1_psi_`ps' = same_direction_t1_psi_`ps' + ((data_`j'_fd - pool_psi_`ps'_`j'_fd)/`j')^2 if test_${test} == 1 & pool_psi_`ps'_`j'_fd != . & data_`j'_fd != .
					replace same_direction_t2_psi_`ps' = same_direction_t2_psi_`ps' + ((data_`j'_fd - pool_psi_`ps'_`j'_fd)/`j')^2 if test_${test} == 2 & pool_psi_`ps'_`j'_fd != . & data_`j'_fd != .
				}
				else {
					replace same_direction_t1_psi_`ps' = same_direction_t1_psi_`ps' + (data_`j'_fd / pool_psi_`ps'_`j'_fd > 0) if test_${test} == 1 & pool_psi_`ps'_`j'_fd != . & data_`j'_fd != . & pool_psi_`ps'_`j'_fd != 0 & data_`j'_fd != 0
					replace same_direction_t2_psi_`ps' = same_direction_t2_psi_`ps' + (data_`j'_fd / pool_psi_`ps'_`j'_fd > 0) if test_${test} == 2 & pool_psi_`ps'_`j'_fd != . & data_`j'_fd != . & pool_psi_`ps'_`j'_fd != 0 & data_`j'_fd != 0

					replace same_direction_t1_psi_`ps' = same_direction_t1_psi_`ps' + .5*(pool_psi_`ps'_`j'_fd == 0) if test_${test} == 1 & pool_psi_`ps'_`j'_fd != . & data_`j'_fd != . & data_`j'_fd != 0
					replace same_direction_t2_psi_`ps' = same_direction_t2_psi_`ps' + .5*(pool_psi_`ps'_`j'_fd == 0) if test_${test} == 2 & pool_psi_`ps'_`j'_fd != . & data_`j'_fd != . & data_`j'_fd != 0

					replace same_direction_t1_psi_`ps' = same_direction_t1_psi_`ps' + (pool_psi_`ps'_`j'_fd == data_`j'_fd) if test_${test} == 1 & pool_psi_`ps'_`j'_fd != . & data_`j'_fd != . & pool_psi_`ps'_`j'_fd == 0 & data_`j'_fd == 0
					replace same_direction_t2_psi_`ps' = same_direction_t2_psi_`ps' + (pool_psi_`ps'_`j'_fd == data_`j'_fd) if test_${test} == 2 & pool_psi_`ps'_`j'_fd != . & data_`j'_fd != . & pool_psi_`ps'_`j'_fd == 0 & data_`j'_fd == 0
				}



				count if pool_psi_`ps'_`j'_fd != . & test_${test} == 1
				local denominator_pool_t1_`ps' `denominator_pool_t1_`ps'' + `r(N)'
				count if pool_psi_`ps'_`j'_fd != . & test_${test} == 2
				local denominator_pool_t2_`ps' `denominator_pool_t2_`ps'' + `r(N)'
			** }
		}
	}
	collapse (sum) same_direction*


	forvalues p = $psi_min($psi_int)$psi_max {
		** quietly {
			local ps = subinstr(string(`p'), ".", "x", .)
			replace same_direction_t1_psi_`ps' = same_direction_t1_psi_`ps'/((`denominator_pool_t1_`ps''))
			replace same_direction_t2_psi_`ps' = same_direction_t2_psi_`ps'/((`denominator_pool_t2_`ps''))
		** }
	}

	if "$trend_method" == "RMSE" {
		forvalues p = $psi_min($psi_int)$psi_max {
			local ps = subinstr(string(`p'), ".", "x", .)
				replace same_direction_t1_psi_`ps' = sqrt(same_direction_t1_psi_`ps')
				replace same_direction_t2_psi_`ps' = sqrt(same_direction_t2_psi_`ps')
		}
	}

// format the first difference results
	generate test = "$test"
	tempfile tmp_fd
	save `tmp_fd', replace
	restore


// find the RMSE
	keep pool_psi_* true_ln_rate test_${test}
	forvalues p = $psi_min($psi_int)$psi_max {
		local ps = subinstr(string(`p'), ".", "x", .)
		quietly generate squared_err_t1_psi_`ps' = (pool_psi_`ps' - true_ln_rate)^2 if test_${test} == 1
		quietly generate squared_err_t2_psi_`ps' = (pool_psi_`ps' - true_ln_rate)^2 if test_${test} == 2
	}
	collapse (mean) squared_err*
	forvalues p = $psi_min($psi_int)$psi_max {
		local ps = subinstr(string(`p'), ".", "x", .)
		quietly generate rmse_t1_psi_`ps' = sqrt(squared_err_t1_psi_`ps')
		quietly generate rmse_t2_psi_`ps' = sqrt(squared_err_t2_psi_`ps')
	}
	generate test = "$test"

// add the first difference results back on
	merge 1:1 test using `tmp_fd', nogen


// reshape so that values of psi are long
	reshape long rmse_t1_psi_ rmse_t2_psi_ same_direction_t1_psi_ same_direction_t2_psi_, i(test) j(psi_str) str
	quietly generate psi = real(subinstr(psi_str,"x",".",.))
	drop psi_str
	drop squared_err*
	rename rmse_t1_psi_ rmse_t1
	rename rmse_t2_psi_ rmse_t2
	rename same_direction_t1_psi_ same_direction_t1
	rename same_direction_t2_psi_ same_direction_t2

// rank by each test
	quietly egen rmse_rank_t1 = rank(rmse_t1), track
	quietly egen rmse_rank_t2 = rank(rmse_t2), track
	if "$trend_method" == "RMSE" {
		quietly egen trend_rank_t1 = rank(same_direction_t1), track
		quietly egen trend_rank_t2 = rank(same_direction_t2), track
	}
	else {
		quietly egen trend_rank_t1 = rank(same_direction_t1), field
		quietly egen trend_rank_t2 = rank(same_direction_t2), field
	}


// calculate the overall rank for each model
	quietly generate tmp = trend_rank_t1 + rmse_rank_t1
	quietly egen total_rank_t1 = rank(tmp), track
	drop tmp
	quietly generate tmp = trend_rank_t2 + rmse_rank_t2
	quietly egen total_rank_t2 = rank(tmp), track
	drop tmp

    //debug
	save "${base_dir}/results/pool_ranks_${test}.dta", replace

// output PV for this test
	save "$temp_dir/8_ensemble_pv_${test}.dta", replace

// close the logs
	timer off 1
	timer list 1
	log close
