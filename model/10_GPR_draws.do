/*
Created:	1 May 2011
Updated:	26 April 2012
Purpose: 	Run GPR to get draws for the final ensemble
*/

// setup the program
	timer on 1
	set mem 4g
	set more off
	set matsize 11000
	set maxvar 32000
	log using "${base_dir}/logs/`model_number'/10_gpr_draws_${super_region}_${super_region_chunk}_${test}.smcl", replace
	display "Finding GPR draws for super region $super_region (part $super_region_chunk) on node `r(o1)'"

** safeguard to prevent errors for people using old runfiles
if "$psi_int" == "" global psi_int = 0.01

// find the optimal value of psi
	use "$temp_dir/9_ensemble_pv.dta", clear
	summarize psi if test == "combined out-of-sample" & total_rank_t1 == 1, meanonly
	if substr(string(`r(mean)'), 5,1) != ""	{
		di in red "This model has an unusual 3 way tie between nonadjacent psi values.  We rounded the average psi to the nearest values of $psi_int."
		local optimal_psi =  round(`r(mean)', $psi_int)
	}
	else local optimal_psi = `r(mean)'
	local optimal_psi_str = subinstr(string(`optimal_psi'), ".", "x", .)

// load in the number of draws to make for each model
	use "$temp_dir/7_submodel_ranks.dta", clear
	forvalues i = 1 / $number_submodels {
		qui levelsof draws_psi_`optimal_psi_str' if spacetime_or_linear == "spacetime" & submodel == `i' & test == "combined out-of-sample", l(spacetime_draws_`i') c
	}
	summarize draws_psi_`optimal_psi_str' if spacetime_or_linear == "spacetime" & test == "combined out-of-sample", meanonly
	local total_spacetime_draws = `r(sum)'

// figure out if the number 1 submodel is a spacetime model
	count if spacetime_or_linear == "spacetime" & test == "combined out-of-sample" & total_rank_t1 == 1
	local draw_top_submodel = (`r(N)' > 0)
	if `draw_top_submodel' == 1 levelsof submodel if spacetime_or_linear == "spacetime" & test == "combined out-of-sample" & total_rank_t1 == 1, l(top_submodel) c

** determine the top rate and cause fraction space-time models in case the user wants draws from either of them
	preserve
	cap {
		keep if spacetime_or_linear == "spacetime" & test == "combined out-of-sample" & dependent_variable == "cf"
		sort total_rank_t1
		local top_cf_model = submodel in 1
	}
	if _rc != 0 local top_cf_model
	restore
	cap {
		keep if spacetime_or_linear == "spacetime" & test == "combined out-of-sample" & dependent_variable == "rate"
		sort total_rank_t1
		local top_rate_model = submodel in 1
	}
	if _rc != 0 local top_rate_model

// load in the earlier GPR results
	use "$temp_dir/5_gpr_${super_region}_${super_region_chunk}_${test}.dta", clear

// if there's no data left, then we don't need to make predictions
	count if predictme_${test} == 1 & spacetime_1 != .
	global make_prediction = `r(N)'
	if `r(N)' > 0 & `total_spacetime_draws' > 0 {

// save a csv of all the data for this region
		preserve
		keep if predictme_${test} == 1
		quietly replace iso3 = "MXR" if iso3 == "MAR"
		keep iso3 age_group year lt_* ln_* *_data_variance_* spacetime_* *_amplitude_* test_${test}
		sort iso3 age_group year

		outsheet using "$temp_dir/10_gpr_input_${super_region}_${super_region_chunk}_${test}.csv", comma replace

		if $top_st_rate == 1 | $top_st_cf == 1 {
			tempfile all_input_data
			save `all_input_data', replace
		}

// write a python script to perform GPR
		qui file open gpr_py using "$temp_dir/run_files/10_gpr_${super_region}_${super_region_chunk}_${test}.py", write replace
		global dvs = `"[""' + subinstr("$dv_list"," ", `"",""',.) + `""]"'
		file write gpr_py	"import sys" _n ///
							"sys.path.append('${code_dir}/model/')" _n ///
							"import GPR_draws" _n ///
							"reload(GPR_draws)" _n ///
							"number_submodels = $number_submodels" _n ///
							"infile = '$temp_dir/10_gpr_input_${super_region}_${super_region_chunk}_${test}.csv'" _n ///
							"outfile = '$temp_dir/10_gpr_output_${super_region}_${super_region_chunk}_${test}.csv'" _n ///
							"scale = $scale" _n ///
							"test = '$test'" _n ///
							"dv_list = ["
		local counter = 0
		foreach d of global dv_list {
			file write gpr_py 	`""`d'""'
			local counter = `counter' + 1
			if `counter' < ${number_submodels} file write gpr_py ","
		}
		file write gpr_py 	"]" _n ///
							"spacetime_iters = ["
		local counter = 0
		forvalues i = 1 / $number_submodels {
			file write gpr_py 	"`spacetime_draws_`i''"
			local counter = `counter' + 1
			if `counter' < ${number_submodels} file write gpr_py ","
		}
		file write gpr_py 	"]" _n
		if `draw_top_submodel' == 1 file write gpr_py 	"top_submodel = `top_submodel'" _n
		else file write gpr_py 							"top_submodel = 0" _n
		file write gpr_py 	"GPR_draws.fit_GPR(infile, outfile, dv_list, scale, number_submodels, test, spacetime_iters, top_submodel)" _n
		file close gpr_py

// run GPR
		shell /usr/local/epd_py25-4.3.0/bin/python "$temp_dir/run_files/10_gpr_${super_region}_${super_region_chunk}_${test}.py"

** save csvs with the data on the top spacetime model in each family in case the user wants to make 1000 draws from either of them
	** rate
		if $top_st_rate == 1 {
			use `all_input_data', clear
			keep iso3 age_group year lt_* ln_* *_data_variance_`top_rate_model' spacetime_`top_rate_model' *_amplitude_`top_rate_model' test_${test}
			foreach var in spacetime_data_variance_ spacetime_ spacetime_amplitude_ {
				cap rename `var'`top_rate_model' `var'1
			}
			sort iso3 age_group year
			outsheet using "$temp_dir/10_gpr_input_${super_region}_${super_region_chunk}_top_rate_model_${test}.csv", comma replace
		}
	** cause fraction
		if $top_st_cf == 1 {
			use `all_input_data', clear
			keep iso3 age_group year lt_* ln_* *_data_variance_`top_cf_model' spacetime_`top_cf_model' *_amplitude_`top_cf_model' test_${test}
			foreach var in spacetime_data_variance_ spacetime_ spacetime_amplitude_ {
				cap rename `var'`top_cf_model' `var'1
			}
			sort iso3 age_group year
			outsheet using "$temp_dir/10_gpr_input_${super_region}_${super_region_chunk}_top_cf_model_${test}.csv", comma replace
		}


// convert the GPR results to Stata
		insheet using "$temp_dir/10_gpr_output_${super_region}_${super_region_chunk}_${test}.csv", comma names clear
		quietly replace iso3 = "MAR" if iso3 == "MXR"
		tempfile gpr_results
		save `gpr_results', replace

// merge GPR results onto the rest of the data
		restore, preserve
		merge m:1 iso3 year age_group using `gpr_results', nogen keep(match master)

// convert to deaths
		local counter = 0
		forvalues i = 1 / $number_submodels {
			local dv: word `i' of $dv_list
			forvalues j = 1/`spacetime_draws_`i'' {
				local counter = `counter' + 1
				qui destring ensemble_d`counter', replace force
				if "`dv'" == "ln_rate" {
					quietly replace ensemble_d`counter' = exp(ensemble_d`counter') * pop
					quietly replace ensemble_d`counter' = envelope if ensemble_d`counter' > envelope & ensemble_d`counter' != .
				}
				else if "`dv'" == "lt_cf" quietly replace ensemble_d`counter' = invlogit(ensemble_d`counter') * envelope
			}
			if "`top_submodel'" == "`i'" {
				forvalues j = 1/100 {
					qui destring top_submodel_d`j', replace force
					if "`dv'" == "ln_rate" {
						quietly replace top_submodel_d`j' = exp(top_submodel_d`j') * pop
						quietly replace top_submodel_d`j' = envelope if top_submodel_d`j' > envelope & top_submodel_d`j' != .
					}
					else if "`dv'" == "lt_cf" quietly replace top_submodel_d`j' = invlogit(top_submodel_d`j') * envelope
				}
			}
		}
	}

// if there's no data, just put in empty variables here
	else {
		local counter = 0
		forvalues i = 1 / $number_submodels {
			forvalues j = 1/`spacetime_draws_`i'' {
				local counter = `counter' + 1
				qui generate ensemble_d`counter' = .
			}
			if "`top_submodel'" == "`i'" {
				forvalues j = 1/100 {
					qui generate top_submodel_d`j' = .
				}
			}
		}
	}

// save final GPR-ed results
	duplicates drop iso3 age year, force
	quietly compress
	save "$temp_dir/10_spacetime_draws_${super_region}_${super_region_chunk}_${test}.dta", replace

** make draws of the top submodel for spacetime or linear models if specified by the user
	if $top_st_rate == 1 & $make_prediction > 0 {
	** write a python script to perform GPR
		qui file open gpr_py using "$temp_dir/run_files/10_gpr_${super_region}_${super_region_chunk}_top_rate_model_${test}.py", write replace
		global dvs = `"[""' + subinstr("$dv_list"," ", `"",""',.) + `""]"'
		file write gpr_py	"import sys" _n ///
							"sys.path.append('${code_dir}/model/')" _n ///
							"import GPR_draws" _n ///
							"reload(GPR_draws)" _n ///
							"number_submodels = 1" _n ///
							"infile = '$temp_dir/10_gpr_input_${super_region}_${super_region_chunk}_top_rate_model_${test}.csv'" _n ///
							"outfile = '$temp_dir/10_gpr_output_${super_region}_${super_region_chunk}_top_rate_model_${test}.csv'" _n ///
							"scale = $scale" _n ///
							"test = '$test'" _n ///
							`"dv_list = ["ln_rate"]"' _n ///
							"spacetime_iters = [1000]" _n ///
							"top_submodel = 0" _n
		file write gpr_py 	"GPR_draws.fit_GPR(infile, outfile, dv_list, scale, number_submodels, test, spacetime_iters, top_submodel)" _n
		file close gpr_py

	** submit GPR script
		shell /usr/local/epd_py25-4.3.0/bin/python "$temp_dir/run_files/10_gpr_${super_region}_${super_region_chunk}_top_rate_model_${test}.py"

	** convert the GPR results to Stata
		insheet using "$temp_dir/10_gpr_output_${super_region}_${super_region_chunk}_top_rate_model_${test}.csv", comma names clear
		quietly replace iso3 = "MAR" if iso3 == "MXR"
		tempfile gpr_results_trm
		save `gpr_results_trm', replace

	** merge GPR results onto the rest of the data
		restore, preserve
		merge m:1 iso3 year age_group using `gpr_results_trm', nogen keep(match master)

	** convert to deaths
		local counter = 0
		local i = 1
		local dv: word `i' of $dv_list
		forvalues j = 1/1000 {
			qui destring ensemble_d`j', replace force
			quietly replace ensemble_d`j' = exp(ensemble_d`j') * pop
			quietly replace ensemble_d`j' = envelope if ensemble_d`j' > envelope & ensemble_d`j' != .
		}

		duplicates drop iso3 age year, force
		quietly compress
		save "$temp_dir/10_spacetime_draws_${super_region}_${super_region_chunk}_top_rate_model_${test}.dta", replace
	}
	if $top_st_cf == 1 & $make_prediction > 0 {
	** write a python script to perform GPR
		qui file open gpr_py using "$temp_dir/run_files/10_gpr_${super_region}_${super_region_chunk}_top_cf_model_${test}.py", write replace
		global dvs = `"[""' + subinstr("$dv_list"," ", `"",""',.) + `""]"'
		file write gpr_py	"import sys" _n ///
							"sys.path.append('${code_dir}/model/')" _n ///
							"import GPR_draws" _n ///
							"reload(GPR_draws)" _n ///
							"number_submodels = 1" _n ///
							"infile = '$temp_dir/10_gpr_input_${super_region}_${super_region_chunk}_top_cf_model_${test}.csv'" _n ///
							"outfile = '$temp_dir/10_gpr_output_${super_region}_${super_region_chunk}_top_cf_model_${test}.csv'" _n ///
							"scale = $scale" _n ///
							"test = '$test'" _n ///
							`"dv_list = ["lt_cf"]"' _n ///
							"spacetime_iters = [1000]" _n ///
							"top_submodel = 0" _n
		file write gpr_py 	"GPR_draws.fit_GPR(infile, outfile, dv_list, scale, number_submodels, test, spacetime_iters, top_submodel)" _n
		file close gpr_py

	** submit GPR script
		shell /usr/local/epd_py25-4.3.0/bin/python "$temp_dir/run_files/10_gpr_${super_region}_${super_region_chunk}_top_cf_model_${test}.py"

		** convert the GPR results to Stata
		insheet using "$temp_dir/10_gpr_output_${super_region}_${super_region_chunk}_top_cf_model_${test}.csv", comma names clear
		quietly replace iso3 = "MAR" if iso3 == "MXR"
		tempfile gpr_results_tcfm
		save `gpr_results_tcfm', replace

	** merge GPR results onto the rest of the data
		restore, preserve
		merge m:1 iso3 year age_group using `gpr_results_tcfm', nogen keep(match master)

	** convert to deaths
		local counter = 0
		local i = 1
		local dv: word `i' of $dv_list
		forvalues j = 1/1000 {
			qui destring ensemble_d`j', replace force
			quietly replace ensemble_d`j' = invlogit(ensemble_d`j') * envelope
		}
		duplicates drop iso3 age year, force
		quietly compress
		save "$temp_dir/10_spacetime_draws_${super_region}_${super_region_chunk}_top_cf_model_${test}.dta", replace
	}


// close the logs
	timer off 1
	timer list 1
	log close
