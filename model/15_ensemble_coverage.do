/*
Created:	1 May 2011
Updated:	26 April 2012
Purpose: 	Calculate predictive validity and coverage for the ensemble model
*/

// setup the program (copy globals into locals, etc)
	timer on 1
	clear
	set mem 6g
	set maxvar 32000
	set more off
	log using "${base_dir}/logs/15_ensemble_coverage_${test}.smcl", replace
	display "Finding ensemble model PV for test $test on node `r(o1)'"

** adding bit to prevent errors for currently running models
if "$psi_int" == "" global psi_int 0.01

// find how many draws were made in total
	use "$temp_dir/9_ensemble_pv.dta", clear
	summarize psi if test == "combined out-of-sample" & total_rank_t1 == 1, meanonly
	if substr(string(`r(mean)'), 5,1) != "" {
		di in red "This model has a 3 way tie between nonadjacent psi values.  We rounded the average psi to the nearest values of $psi_int."
		local optimal_psi = round(`r(mean)', $psi_int)
		local optimal_psi_str =  subinstr(string(`optimal_psi'), ".", "x", .)
	}
	else {
    local optimal_psi = `r(mean)'
    local optimal_psi_str = subinstr(string(`r(mean)'), ".", "x", .)
    }

    // save optimal psi in a csv for fun
    drop *
    set obs 1
    gen optimal_psi = `optimal_psi'
    outsheet using "${base_dir}/optimal_psi.csv", comma replace

	use "$temp_dir/7_submodel_ranks.dta", clear
	summarize draws_psi_`optimal_psi_str' if test == "combined out-of-sample", meanonly
	local total_draws = `r(sum)'

// load in the draws
	local num_sr: word count $super_regions
	forvalues i = 1/`num_sr' {
		local this_sr: word `i' of $super_regions
		local num_chunks: word `i' of $chunks_per_sr
		clear
		forvalues this_chunk = 1/`num_chunks' {
            append using "$temp_dir/10_spacetime_draws_`this_sr'_`this_chunk'_${test}.dta"
		}
		quietly merge 1:1 iso3 year age using "$temp_dir/11_linear_draws_`this_sr'_${test}.dta", keep(match master) nogen
		quietly compress
		tempfile tmp_`i'
        save `tmp_`i''

	}
	clear
	forvalues i = 1/`num_sr' {
		quietly append using `tmp_`i''
	}
	keep iso3 year age ensemble_* top_submodel_*
	tempfile draws
    save `draws'

// load in the data
	use "$temp_dir/input_database.dta", clear

// add on the NSV and draws
	merge m:1 iso3 year age using `draws', keep(match master) nogen

	merge 1:1 n using "$temp_dir/4_gpr_parameters_${test}.dta", keep(match master) nogen

// keep just the rows with test data
	if "${test}" == "insample" replace test_insample = 2
	keep if ln_rate != . & test_${test} == 2

// convert all the draws to log rates
	forvalues i = 1/`total_draws' {
		quietly generate ensemble_ln_rate_d`i' = ln(ensemble_d`i' / pop)
		drop ensemble_d`i'
	}
	quietly egen ensemble_mean = rowmean(ensemble_ln_rate_d*)
	quietly egen ensemble_sd = rowsd(ensemble_ln_rate_d*)
	forvalues i = 1/100 {
		quietly generate top_submodel_ln_rate_d`i' = ln(top_submodel_d`i' / pop)
		drop top_submodel_d`i'
	}
	quietly egen top_submodel_mean = rowmean(top_submodel_ln_rate_d*)
	quietly egen top_submodel_sd = rowsd(top_submodel_ln_rate_d*)

// find the error in the test data
	quietly generate true_ln_rate = ln_rate
	quietly generate squared_error = (ensemble_mean - true_ln_rate)^2
	quietly generate top_submodel_squared_error = (top_submodel_mean - true_ln_rate)^2

// find coverage of each datapoint
	quietly generate coverage_lower_t2 = inrange(true_ln_rate, ensemble_mean - (1.96 * sqrt(ln_rate_sd^2 + ensemble_sd^2)), ensemble_mean + (1.96 * sqrt(ln_rate_sd^2 + ensemble_sd^2)))
	quietly generate top_submodel_coverage_lower_t2 = inrange(true_ln_rate, top_submodel_mean - (1.96 * sqrt(ln_rate_sd^2 + top_submodel_sd^2)), top_submodel_mean + (1.96 * sqrt(ln_rate_sd^2 + top_submodel_sd^2)))
	quietly generate coverage_upper_t2 = inrange(true_ln_rate, ensemble_mean - (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + ensemble_sd^2)), ensemble_mean + (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + ensemble_sd^2)))
	quietly generate top_submodel_coverage_upper_t2 = inrange(true_ln_rate, top_submodel_mean - (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + top_submodel_sd^2)), top_submodel_mean + (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + top_submodel_sd^2)))

// find aggregate coverage/RMSE for this holdout
	preserve
	collapse (mean) *coverage* squared_error top_submodel_squared_error
	quietly generate rmse_t2 = sqrt(squared_error)
	quietly generate top_submodel_rmse_t2 = sqrt(top_submodel_squared_error)
	quietly generate test = "$test"
	tempfile rmse_cov
    save `rmse_cov'
	restore

// find whether the trend is in the correct direction
	collapse (mean) ensemble_mean true_ln_rate top_submodel_mean, by(iso3 year age)
	egen ca = group(iso3 age)
	tsset ca year
    if "$trend_window" == "" global trend_window 1
	forvalues j = $trend_window_min / $trend_window {
			 generate data_`j'_fd = true_ln_rate - L`j'.true_ln_rate
			 generate ensemble_`j'_fd = ensemble_mean - L`j'.ensemble_mean
			 generate top_submodel_`j'_fd = top_submodel_mean - L`j'.top_submodel_mean

	}

	local denominator_ensemble_coverage 0
	generate same_direction_t2 = 0
	generate top_submodel_same_direction_t2 = 0
	forvalues j = $trend_window_min / $trend_window {

		if "$trend_method" == "RMSE" {
			replace same_direction_t2 = same_direction_t2 + (data_`j'_fd - ensemble_`j'_fd)^2 if data_`j'_fd != . & ensemble_`j'_fd !=.
			replace top_submodel_same_direction_t2 = top_submodel_same_direction_t2 + (data_`j'_fd - top_submodel_`j'_fd)^2 if data_`j'_fd != . & ensemble_`j'_fd !=.
		}
		else {
			replace same_direction_t2 = same_direction_t2 + (data_`j'_fd / ensemble_`j'_fd >= 0) if data_`j'_fd != . & ensemble_`j'_fd !=.
			replace top_submodel_same_direction_t2 = top_submodel_same_direction_t2 + (data_`j'_fd / top_submodel_`j'_fd >= 0) if data_`j'_fd != . & ensemble_`j'_fd !=.
			replace same_direction_t2 = same_direction_t2 + .5*(ensemble_`j'_fd == 0) if data_`j'_fd != . & ensemble_`j'_fd !=. & data_`j'_fd != 0
			replace top_submodel_same_direction_t2 = top_submodel_same_direction_t2 + .5*(top_submodel_`j'_fd == 0) if data_`j'_fd != . & ensemble_`j'_fd !=. & data_`j'_fd != 0
			replace same_direction_t2 = same_direction_t2 + (data_`j'_fd == ensemble_`j'_fd) if data_`j'_fd == 0 & ensemble_`j'_fd == 0
			replace top_submodel_same_direction_t2 = top_submodel_same_direction_t2 + (data_`j'_fd == top_submodel_`j'_fd) if data_`j'_fd == 0 & ensemble_`j'_fd == 0

		}


		count if ensemble_`j'_fd != .
		local denominator_ensemble_coverage `denominator_ensemble_coverage' + `r(N)'
	}

	collapse (sum) same_direction_t2 top_submodel_same_direction_t2
	replace same_direction_t2 = same_direction_t2/(`denominator_ensemble_coverage')
	replace top_submodel_same_direction_t2 = top_submodel_same_direction_t2/(`denominator_ensemble_coverage')
	quietly generate test = "$test"
// join the results together
	merge 1:1 test using `rmse_cov', keep(match master) nogen

	keep test rmse_t2 same_direction_t2 coverage_lower_t2 coverage_upper_t2 top_submodel_rmse_t2 top_submodel_same_direction_t2 top_submodel_coverage_lower_t2 top_submodel_coverage_upper_t2

// output PV for this test
	save "$temp_dir/15_ensemble_coverage_${test}.dta", replace

** calculate coverage for top submodel in a given family as specified by the user

** top spacetime cause fraction model
	if $top_st_cf == 1 {
	** bring all of the data together and keep only the test data
		local num_sr: word count $super_regions
		forvalues i = 1/`num_sr' {
			local this_sr: word `i' of $super_regions
			local num_chunks: word `i' of $chunks_per_sr
			forvalues this_chunk = 1/`num_chunks' {
				quietly append using "$temp_dir/10_spacetime_draws_`this_sr'_`this_chunk'_top_cf_model_${test}.dta"
			}
		}
		use "/home/j/Project/Causes of Death/codem/models/${cause}/${model_name}/results/death_draws_top_st_cf_model.dta", clear
		merge 1:m iso3 year age using "$temp_dir/input_database.dta", keep(match using) nogen
		merge 1:1 n using "$temp_dir/4_gpr_parameters_${test}.dta", keep(match master) nogen
		if "${test}" == "insample" replace test_insample = 2
		keep if ln_rate != . & test_${test} == 2

	// convert all the draws to log rates
		forvalues i = 1/1000 {
			quietly generate ln_rate_d`i' = ln(draw_`i' / pop)
			drop draw_`i'
		}
		quietly egen draw_mean = rowmean(ln_rate_d*)
		quietly egen draw_sd = rowsd(ln_rate_d*)

	// find the error in the test data
		quietly generate true_ln_rate = ln_rate
		quietly generate squared_error = (draw_mean - true_ln_rate)^2

	// find coverage of each datapoint
		quietly generate coverage_lower_t2 = inrange(true_ln_rate, draw_mean - (1.96 * sqrt(ln_rate_sd^2 + draw_sd^2)), draw_mean + (1.96 * sqrt(ln_rate_sd^2 + draw_sd^2)))
		quietly generate coverage_upper_t2 = inrange(true_ln_rate, draw_mean - (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + draw_sd^2)), draw_mean + (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + draw_sd^2)))

	// find aggregate coverage/RMSE for this holdout
		preserve
		collapse (mean) *coverage* squared_error
		quietly generate rmse_t2 = sqrt(squared_error)
		quietly generate test = "$test"
		tempfile rmse_cov
        save `rmse_cov'
		restore

	// find whether the trend is in the correct direction
        collapse (mean) draw_mean true_ln_rate, by(iso3 year age)
        egen ca = group(iso3 age)
        tsset ca year
        forvalues j = $trend_window_min / $trend_window {
                 generate data_`j'_fd = true_ln_rate - L`j'.true_ln_rate
                 generate draw_`j'_fd = draw_mean - L`j'.draw_mean
        }

        local denominator_ensemble_coverage 0
        generate same_direction_t2 = 0
        generate top_submodel_same_direction_t2 = 0
        forvalues j = $trend_window_min / $trend_window {

            if "$trend_method" == "RMSE" {
                replace same_direction_t2 = same_direction_t2 + ((data_`j'_fd - draw_`j'_fd)/`j')^2 if data_`j'_fd != . & draw_`j'_fd !=.
            }
            else {
                replace same_direction_t2 = same_direction_t2 + (data_`j'_fd / draw_`j'_fd >= 0) if data_`j'_fd != . & draw_`j'_fd !=.
            }


            count if draw_`j'_fd != .
            local denominator_draw_coverage `denominator_draw_coverage' + `r(N)'
        }



        collapse (sum) same_direction_t2
        replace same_direction_t2 = same_direction_t2/(`denominator_draw_coverage')
        if "$trend_method" == "RMSE" {
            replace same_direction_t2 = sqrt(same_direction_t2)
        }

        quietly generate test = "$test"

// join the results together
		merge 1:1 test using `rmse_cov', keep(match master) nogen

		keep test rmse_t2 same_direction_t2 coverage_lower_t2 coverage_upper_t2

	// output PV for this test
		save "$temp_dir/15_ensemble_coverage_${test}_top_st_cf_model.dta", replace
	}

** top spacetime rate model
	if $top_st_rate == 1 {
		local num_sr: word count $super_regions
		forvalues i = 1/`num_sr' {
			local this_sr: word `i' of $super_regions
			local num_chunks: word `i' of $chunks_per_sr
			forvalues this_chunk = 1/`num_chunks' {
				quietly append using "$temp_dir/10_spacetime_draws_`this_sr'_`this_chunk'_top_rate_model_${test}.dta"
			}
		}
		merge 1:m iso3 year age using "$temp_dir/input_database.dta", keep(match using) nogen
		merge 1:1 n using "$temp_dir/4_gpr_parameters_${test}.dta", keep(match master) nogen
		if "${test}" == "insample" replace test_insample = 2
		keep if ln_rate != . & test_${test} == 2

	// convert all the draws to log rates
		forvalues i = 1/1000 {
			quietly generate ln_rate_d`i' = ln(draw_`i' / pop)
			drop draw_`i'
		}
		quietly egen draw_mean = rowmean(ln_rate_d*)
		quietly egen draw_sd = rowsd(ln_rate_d*)

	// find the error in the test data
		quietly generate true_ln_rate = ln_rate
		quietly generate squared_error = (draw_mean - true_ln_rate)^2

	// find coverage of each datapoint
		quietly generate coverage_lower_t2 = inrange(true_ln_rate, draw_mean - (1.96 * sqrt(ln_rate_sd^2 + draw_sd^2)), draw_mean + (1.96 * sqrt(ln_rate_sd^2 + draw_sd^2)))
		quietly generate coverage_upper_t2 = inrange(true_ln_rate, draw_mean - (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + draw_sd^2)), draw_mean + (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + draw_sd^2)))

	// find aggregate coverage/RMSE for this holdout
		preserve
		collapse (mean) *coverage* squared_error
		quietly generate rmse_t2 = sqrt(squared_error)
		quietly generate test = "$test"
		tempfile rmse_cov
        save `rmse_cov'
		restore

	// find whether the trend is in the correct direction
        collapse (mean) draw_mean true_ln_rate, by(iso3 year age)
        egen ca = group(iso3 age)
        tsset ca year
        forvalues j = $trend_window_min / $trend_window {
                 generate data_`j'_fd = true_ln_rate - L`j'.true_ln_rate
                 generate draw_`j'_fd = draw_mean - L`j'.draw_mean
        }

        local denominator_ensemble_coverage 0
        generate same_direction_t2 = 0
        generate top_submodel_same_direction_t2 = 0
        forvalues j = $trend_window_min / $trend_window {

            if "$trend_method" == "RMSE" {
                replace same_direction_t2 = same_direction_t2 + ((data_`j'_fd - draw_`j'_fd)/`j')^2 if data_`j'_fd != . & draw_`j'_fd !=.
            }
            else {
                replace same_direction_t2 = same_direction_t2 + (data_`j'_fd / draw_`j'_fd >= 0) if data_`j'_fd != . & draw_`j'_fd !=.
            }


            count if draw_`j'_fd != .
            local denominator_draw_coverage `denominator_draw_coverage' + `r(N)'
        }



        collapse (sum) same_direction_t2
        replace same_direction_t2 = same_direction_t2/(`denominator_draw_coverage')
        if "$trend_method" == "RMSE" {
            replace same_direction_t2 = sqrt(same_direction_t2)
        }

        quietly generate test = "$test"

// join the results together
		merge 1:1 test using `rmse_cov', keep(match master) nogen

		keep test rmse_t2 same_direction_t2 coverage_lower_t2 coverage_upper_t2

	// output PV for this test
		save "$temp_dir/15_ensemble_coverage_${test}_top_st_rate_model.dta", replace
	}

** top linear cause fraction model
	if $top_lin_cf == 1 {
		local num_sr: word count $super_regions
		forvalues i = 1/`num_sr' {
			local this_sr: word `i' of $super_regions
			quietly append using "$temp_dir/11_top_lin_cf_model_draws_`this_sr'_${test}.dta"
		}
		merge 1:m iso3 year age using "$temp_dir/input_database.dta", keep(match using) nogen
		merge 1:1 n using "$temp_dir/4_gpr_parameters_${test}.dta", keep(match master) nogen
		if "${test}" == "insample" replace test_insample = 2
		keep if ln_rate != . & test_${test} == 2

	// convert all the draws to log rates
		forvalues i = 1/1000 {
			quietly generate ln_rate_d`i' = ln(draw_`i' / pop)
			drop draw_`i'
		}
		quietly egen draw_mean = rowmean(ln_rate_d*)
		quietly egen draw_sd = rowsd(ln_rate_d*)

	// find the error in the test data
		quietly generate true_ln_rate = ln_rate
		quietly generate squared_error = (draw_mean - true_ln_rate)^2

	// find coverage of each datapoint
		quietly generate coverage_lower_t2 = inrange(true_ln_rate, draw_mean - (1.96 * sqrt(ln_rate_sd^2 + draw_sd^2)), draw_mean + (1.96 * sqrt(ln_rate_sd^2 + draw_sd^2)))
		quietly generate coverage_upper_t2 = inrange(true_ln_rate, draw_mean - (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + draw_sd^2)), draw_mean + (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + draw_sd^2)))

	// find aggregate coverage/RMSE for this holdout
		preserve
		collapse (mean) *coverage* squared_error
		quietly generate rmse_t2 = sqrt(squared_error)
		quietly generate test = "$test"
		tempfile rmse_cov
        save `rmse_cov'
		restore

	// find whether the trend is in the correct direction
        collapse (mean) draw_mean true_ln_rate, by(iso3 year age)
        egen ca = group(iso3 age)
        tsset ca year
        forvalues j = $trend_window_min / $trend_window {
                 generate data_`j'_fd = true_ln_rate - L`j'.true_ln_rate
                 generate draw_`j'_fd = draw_mean - L`j'.draw_mean
        }

        local denominator_ensemble_coverage 0
        generate same_direction_t2 = 0
        generate top_submodel_same_direction_t2 = 0
        forvalues j = $trend_window_min / $trend_window {

            if "$trend_method" == "RMSE" {
                replace same_direction_t2 = same_direction_t2 + ((data_`j'_fd - draw_`j'_fd)/`j')^2 if data_`j'_fd != . & draw_`j'_fd !=.
            }
            else {
                replace same_direction_t2 = same_direction_t2 + (data_`j'_fd / draw_`j'_fd >= 0) if data_`j'_fd != . & draw_`j'_fd !=.
            }


            count if draw_`j'_fd != .
            local denominator_draw_coverage `denominator_draw_coverage' + `r(N)'
        }



        collapse (sum) same_direction_t2
        replace same_direction_t2 = same_direction_t2/(`denominator_draw_coverage')
        if "$trend_method" == "RMSE" {
            replace same_direction_t2 = sqrt(same_direction_t2)
        }

        quietly generate test = "$test"

// join the results together
		merge 1:1 test using `rmse_cov', keep(match master) nogen

		keep test rmse_t2 same_direction_t2 coverage_lower_t2 coverage_upper_t2

	// output PV for this test
		save "$temp_dir/15_ensemble_coverage_${test}_top_lin_cf_model.dta", replace
	}

** top linear rate model
	if $top_lin_rate == 1 {
		local num_sr: word count $super_regions
		forvalues i = 1/`num_sr' {
			local this_sr: word `i' of $super_regions
			quietly append using "$temp_dir/11_top_lin_rate_model_draws_`this_sr'_${test}.dta"
		}
		merge 1:m iso3 year age using "$temp_dir/input_database.dta", keep(match using) nogen
		merge 1:1 n using "$temp_dir/4_gpr_parameters_${test}.dta", keep(match master) nogen

	// convert all the draws to log rates
		forvalues i = 1/1000 {
			quietly generate ln_rate_d`i' = ln(draw_`i' / pop)
			drop draw_`i'
		}
		quietly egen draw_mean = rowmean(ln_rate_d*)
		quietly egen draw_sd = rowsd(ln_rate_d*)

	// find the error in the test data
		quietly generate true_ln_rate = ln_rate
		quietly generate squared_error = (draw_mean - true_ln_rate)^2

	// find coverage of each datapoint
		quietly generate coverage_lower_t2 = inrange(true_ln_rate, draw_mean - (1.96 * sqrt(ln_rate_sd^2 + draw_sd^2)), draw_mean + (1.96 * sqrt(ln_rate_sd^2 + draw_sd^2)))
		quietly generate coverage_upper_t2 = inrange(true_ln_rate, draw_mean - (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + draw_sd^2)), draw_mean + (1.96 * sqrt(ln_rate_sd^2 + (1.4826*mad_nsv_universal_ln_rate)^2 + draw_sd^2)))

	// find aggregate coverage/RMSE for this holdout
		preserve
		collapse (mean) *coverage* squared_error
		quietly generate rmse_t2 = sqrt(squared_error)
		quietly generate test = "$test"
		tempfile rmse_cov
        save `rmse_cov'

		restore

	// find whether the trend is in the correct direction
        collapse (mean) draw_mean true_ln_rate, by(iso3 year age)
        egen ca = group(iso3 age)
        tsset ca year
        forvalues j = $trend_window_min / $trend_window {
                 generate data_`j'_fd = true_ln_rate - L`j'.true_ln_rate
                 generate draw_`j'_fd = draw_mean - L`j'.draw_mean
        }

        local denominator_ensemble_coverage 0
        generate same_direction_t2 = 0
        generate top_submodel_same_direction_t2 = 0
        forvalues j = $trend_window_min / $trend_window {

            if "$trend_method" == "RMSE" {
                replace same_direction_t2 = same_direction_t2 + ((data_`j'_fd - draw_`j'_fd)/`j')^2 if data_`j'_fd != . & draw_`j'_fd !=.
            }
            else {
                replace same_direction_t2 = same_direction_t2 + (data_`j'_fd / draw_`j'_fd >= 0) if data_`j'_fd != . & draw_`j'_fd !=.
            }


            count if draw_`j'_fd != .
            local denominator_draw_coverage `denominator_draw_coverage' + `r(N)'
        }



        collapse (sum) same_direction_t2
        replace same_direction_t2 = same_direction_t2/(`denominator_draw_coverage')
        if "$trend_method" == "RMSE" {
            replace same_direction_t2 = sqrt(same_direction_t2)
        }

        quietly generate test = "$test"

// join the results together
		merge 1:1 test using `rmse_cov', keep(match master) nogen

		keep test rmse_t2 same_direction_t2 coverage_lower_t2 coverage_upper_t2

	// output PV for this test
		save "$temp_dir/15_ensemble_coverage_${test}_top_lin_rate_model_.dta", replace
	}

// close the logs
	timer off 1
	timer list 1
	log close
