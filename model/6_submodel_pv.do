/*
Created:	1 May 2011
Updated:	6 June 2012
Purpose: 	Calculate predictive validity for each submodel
*/

// setup the program
	clear all
	cap log close
	timer on 1
	set mem 5g
	set more off
	set matsize 11000
	set maxvar 32000
	log using "${base_dir}/logs/6_submodel_pv_${test}.smcl", replace
	display "Finding submodel PV for test $test on node `r(o1)'"

// load in the GPR results
	local num_sr: word count $super_regions
	forvalues i = 1/`num_sr' {
		local this_sr: word `i' of $super_regions
		local num_chunks: word `i' of $chunks_per_sr
		forvalues this_chunk = 1/`num_chunks' {
			** adding capture for not so that models with very little data do not break because GPR can't run in places with no data
			cap append using "$temp_dir/5_gpr_`this_sr'_`this_chunk'_${test}.dta", keep(iso3 year age spacetime_* gpr_*_mean lt_cf ln_rate ln_rate_sd envelope pop test_${test} n simple_spacetime_ln_rate)
		}
	}

// merge on the linear predictions
	forvalues i = 1 / $number_submodel_chunks {
		quietly merge 1:1 n using "$temp_dir/2_linear_country_`i'_${test}.dta", keepusing(linear_*) keep(match master) nogen
	}

// keep just the test data
	if "${test}" == "insample" replace test_insample = 1
	quietly generate true_ln_rate = ln_rate
	drop if test_${test} == 0 | ln_rate == .
	compress

// change output to log death rates for each submodel in order to rank them (previously created new variables but some models exceeded Stata's maximum variable limit)
	forvalues i = 1 / $number_submodels {
		rename linear_`i' ln_rate_`i'_linear
		rename gpr_`i'_spacetime_mean ln_rate_`i'_spacetime
		local dv: word `i' of $dv_list
		quietly {
			if "`dv'" == "lt_cf" {
				quietly replace ln_rate_`i'_linear = ln(invlogit(ln_rate_`i'_linear) * (envelope / pop))
				quietly replace ln_rate_`i'_spacetime = ln(invlogit(ln_rate_`i'_spacetime) * (envelope / pop))
			}
		}
	}

** find first differences in predictions and data
	preserve
	collapse (mean) ln_rate_*_linear ln_rate_*_spacetime true_ln_rate, by(iso3 year age test_${test})
	egen ca = group(iso3 age test_${test})
	tsset ca year

	if "$trend_window" == "" global trend_window 1
	forvalues j = $trend_window_min / $trend_window {
			 generate data_`j'_fd = true_ln_rate - L`j'.true_ln_rate
	}


	forvalues i = 1 / $number_submodels {
		forvalues j =  $trend_window_min / ${trend_window} {

			 generate linear_`i'_`j'_fd = ln_rate_`i'_linear - L`j'.ln_rate_`i'_linear
			 generate spacetime_`i'_`j'_fd = ln_rate_`i'_spacetime - L`j'.ln_rate_`i'_spacetime
		}
	}



** count how often the first difference of the prediction lines up with the data
	forvalues i = 1 / $number_submodels {
		quietly {
			generate same_direction_t1linear_`i' = 0
			generate same_direction_t2linear_`i' = 0
			generate same_direction_t1spacetime_`i' = 0
			generate same_direction_t2spacetime_`i' = 0

			local denominator_t1linear_`i' 0
			local denominator_t1spacetime_`i' 0
			local denominator_t2linear_`i' 0
			local denominator_t2spacetime_`i' 0

            if "$good_data_c_years_path" != "" {
                local good_denominator_t1linear_`i' 0
                local good_denominator_t1spacetime_`i' 0
                local good_denominator_t2linear_`i' 0
                local good_denominator_t2spacetime_`i' 0
            }

			}
		forvalues j =  $trend_window_min / ${trend_window} {
			if "$trend_method" == "RMSE" {
				replace same_direction_t1linear_`i' = same_direction_t1linear_`i' + ((data_`j'_fd - linear_`i'_`j'_fd)/`j')^2 if test_${test} == 1 & linear_`i'_`j'_fd != . & data_`j'_fd != .
				replace same_direction_t2linear_`i' = same_direction_t2linear_`i' + ((data_`j'_fd - linear_`i'_`j'_fd)/`j')^2 if test_${test} == 2 & linear_`i'_`j'_fd != . & data_`j'_fd != .


				replace same_direction_t1spacetime_`i' = same_direction_t1spacetime_`i' + ((data_`j'_fd - spacetime_`i'_`j'_fd)/`j')^2 if test_${test} == 1 & spacetime_`i'_`j'_fd !=. & data_`j'_fd != .
				replace same_direction_t2spacetime_`i' = same_direction_t2spacetime_`i' + ((data_`j'_fd - spacetime_`i'_`j'_fd)/`j')^2 if test_${test} == 2 & spacetime_`i'_`j'_fd !=. & data_`j'_fd != .
			}
			else {
					replace same_direction_t1linear_`i' = same_direction_t1linear_`i' + (data_`j'_fd / linear_`i'_`j'_fd > 0) if test_${test} == 1 & linear_`i'_`j'_fd != . & data_`j'_fd != . & linear_`i'_`j'_fd != 0 & data_`j'_fd != 0
					replace same_direction_t2linear_`i' = same_direction_t2linear_`i' + (data_`j'_fd / linear_`i'_`j'_fd > 0) if test_${test} == 2 & linear_`i'_`j'_fd != . & data_`j'_fd != . & linear_`i'_`j'_fd != 0 & data_`j'_fd != 0
					replace same_direction_t1linear_`i' = same_direction_t1linear_`i' + .5*(linear_`i'_`j'_fd ==0) if test_${test} == 1 & data_`j'_fd != . & data_`j'_fd != 0
					replace same_direction_t2linear_`i' = same_direction_t2linear_`i' + .5*(linear_`i'_`j'_fd ==0) if test_${test} == 2 & data_`j'_fd != . & data_`j'_fd != 0
					replace same_direction_t1linear_`i' = same_direction_t1linear_`i' + (linear_`i'_`j'_fd == data_`j'_fd) if test_${test} == 1 & data_`j'_fd != . & data_`j'_fd == 0 & linear_`i'_`j'_fd == 0
					replace same_direction_t2linear_`i' = same_direction_t2linear_`i' + (linear_`i'_`j'_fd == data_`j'_fd) if test_${test} == 2 & data_`j'_fd != . & data_`j'_fd == 0 & linear_`i'_`j'_fd == 0


					replace same_direction_t1spacetime_`i' = same_direction_t1spacetime_`i' + (data_`j'_fd / spacetime_`i'_`j'_fd > 0) if test_${test} == 1 & spacetime_`i'_`j'_fd !=. & data_`j'_fd != . & spacetime_`i'_`j'_fd != 0
					replace same_direction_t2spacetime_`i' = same_direction_t2spacetime_`i' + (data_`j'_fd / spacetime_`i'_`j'_fd > 0) if test_${test} == 2 & spacetime_`i'_`j'_fd !=. & data_`j'_fd != . & spacetime_`i'_`j'_fd != 0
					replace same_direction_t1spacetime_`i' = same_direction_t1spacetime_`i' + .5*(spacetime_`i'_`j'_fd ==0) if test_${test} == 1 & data_`j'_fd != . & data_`j'_fd != 0
					replace same_direction_t2spacetime_`i' = same_direction_t2spacetime_`i' + .5*(spacetime_`i'_`j'_fd ==0) if test_${test} == 2 & data_`j'_fd != . & data_`j'_fd != 0
					replace same_direction_t1spacetime_`i' = same_direction_t1spacetime_`i' + (spacetime_`i'_`j'_fd == data_`j'_fd) if test_${test} == 1 & data_`j'_fd != . & data_`j'_fd == 0 & spacetime_`i'_`j'_fd == 0
					replace same_direction_t2spacetime_`i' = same_direction_t2spacetime_`i' + (spacetime_`i'_`j'_fd == data_`j'_fd) if test_${test} == 2 & data_`j'_fd != . & data_`j'_fd == 0 & spacetime_`i'_`j'_fd == 0
			}
			count if linear_`i'_`j'_fd != . & test_${test} == 1
			local denominator_t1linear_`i' `denominator_t1linear_`i''+`r(N)'
			count if linear_`i'_`j'_fd != . & test_${test} == 2
			local denominator_t2linear_`i' `denominator_t2linear_`i''+`r(N)'

			count if spacetime_`i'_`j'_fd != . & test_${test} == 1
			local denominator_t1spacetime_`i' `denominator_t1spacetime_`i''+`r(N)'
			count if spacetime_`i'_`j'_fd != . & test_${test} == 2
			local denominator_t2spacetime_`i' `denominator_t2spacetime_`i''+`r(N)'

            if "$good_data_c_years_path" != "" {
                count if linear_`i'_`j'_fd != . & test_${test} == 1 & good_data
                local good_denominator_t1linear_`i' `good_denominator_t1linear_`i''+`r(N)'
                count if linear_`i'_`j'_fd != . & test_${test} == 2 & good_data
                local good_denominator_t2linear_`i' `good_denominator_t2linear_`i''+`r(N)'

                count if spacetime_`i'_`j'_fd != . & test_${test} == 1 & good_data
                local good_denominator_t1spacetime_`i' `good_denominator_t1spacetime_`i''+`r(N)'
                count if spacetime_`i'_`j'_fd != . & test_${test} == 2 & good_data
                local good_denominator_t2spacetime_`i' `good_denominator_t2spacetime_`i''+`r(N)'
            }

		}
	}

    // debug
    save "${base_dir}/results/submodel_trends_${test}.dta", replace
	if "$good_data_c_years_path" != ""{
        save "$temp_dir/before_collapse.dta",replace
    }
	collapse (sum) same_direction*, fast


	forvalues i = 1 / $number_submodels {
		replace same_direction_t1linear_`i' = same_direction_t1linear_`i'/(`denominator_t1linear_`i'')
		replace same_direction_t2linear_`i' = same_direction_t2linear_`i'/(`denominator_t2linear_`i'')
		replace same_direction_t1spacetime_`i' = same_direction_t1spacetime_`i'/(`denominator_t1spacetime_`i'')
		replace same_direction_t2spacetime_`i' = same_direction_t2spacetime_`i'/(`denominator_t2spacetime_`i'')
	}

	if "$trend_method" == "RMSE" {
		forvalues i = 1 / $number_submodels {
			replace same_direction_t1linear_`i' = sqrt(same_direction_t1linear_`i')
			replace same_direction_t2linear_`i' = sqrt(same_direction_t2linear_`i')
			replace same_direction_t1spacetime_`i' = sqrt(same_direction_t1spacetime_`i')
			replace same_direction_t2spacetime_`i' = sqrt(same_direction_t2spacetime_`i')
		}
	}


// format the first difference results
	generate test = "$test"
	reshape long same_direction_t1 same_direction_t2, i(test) j(submodel) string
	tempfile fd
	save `fd', replace
	restore, preserve

// find the RMSE
	keep test_${test} ln_rate_*_linear ln_rate_*_spacetime true_ln_rate
	forvalues i = 1 / $number_submodels {
		quietly {
			generate squared_err_t1_`i'_linear = (ln_rate_`i'_linear - true_ln_rate)^2 if test_${test} == 1
			generate squared_err_t2_`i'_linear = (ln_rate_`i'_linear - true_ln_rate)^2 if test_${test} == 2
			generate squared_err_t1_`i'_spacetime = (ln_rate_`i'_spacetime - true_ln_rate)^2 if test_${test} == 1
			generate squared_err_t2_`i'_spacetime = (ln_rate_`i'_spacetime - true_ln_rate)^2 if test_${test} == 2
		}
	}
	collapse (mean) squared_err_*
	forvalues i = 1 / $number_submodels {
		quietly {
			generate rmse_t1linear_`i' = sqrt(squared_err_t1_`i'_linear)
			generate rmse_t2linear_`i' = sqrt(squared_err_t2_`i'_linear)
			generate rmse_t1spacetime_`i' = sqrt(squared_err_t1_`i'_spacetime)
			generate rmse_t2spacetime_`i' = sqrt(squared_err_t2_`i'_spacetime)
			drop squared_err_t1_`i'_linear squared_err_t1_`i'_spacetime squared_err_t2_`i'_linear squared_err_t2_`i'_spacetime
		}
	}

// reshape so that there's a row for each model
	quietly generate test = "$test"
	reshape long rmse_t1 rmse_t2, i(test) j(submodel) string

// add all the tests together
	merge 1:1 submodel using `fd', nogen


// rank by each test
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
	}

// calculate the overall rank for each model
	quietly {
		generate tmp = trend_rank_t1 + rmse_rank_t1
		egen total_rank_t1 = rank(tmp), track
		drop tmp
		generate tmp = trend_rank_t2 + rmse_rank_t2
		egen total_rank_t2 = rank(tmp), track
		drop tmp
	}
    //debug
	save "${base_dir}/results/submodel_ranks_${test}.dta", replace


// output PV for this test
	compress
	save "$temp_dir/6_submodel_pv_${test}.dta", replace

// close the logs
	timer off 1
	timer list 1
	log close
