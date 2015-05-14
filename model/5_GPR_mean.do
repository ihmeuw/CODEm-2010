/*
Created:	1 May 2011
Updated:	11 April 2012
Purpose: 	Run GPR to get the expectation (mean prediction)
*/

// setup the program
	timer on 1
	set mem 4g
	set more off
	set matsize 11000
	set maxvar 32000
	log using "${base_dir}/logs/`model_number'/5_GPR_${super_region}_${super_region_chunk}_${test}.smcl", replace
	display "Running GPR for super region $super_region (part $super_region_chunk) on node `r(o1)'"

// load in the spacetime results
	use "$temp_dir/3_spacetime_${super_region}_${super_region_chunk}_${test}.dta", clear

// merge on the GPR parameters
	merge 1:1 n using "$temp_dir/4_gpr_parameters_${test}.dta", nogen keep(match master)

// loop through each submodel
	forvalues i = 1 / $number_submodels {

// set amplitude as the MAD converted to SD
		quietly generate spacetime_amplitude_`i' = 1.4826 * spacetime_mad_amp_`i'

// set data variance as the MAD of the residuals and adding on the datapoint-specific stochastic
		local dv: word `i' of $dv_list
		quietly generate spacetime_data_variance_`i' = ((1.4826 * mad_nsv_universal_`dv')^2 + `dv'_sd^2)

// if data variance is somehow <= 0 , replace it with the median non-zero, because it's a floating point error somewhere
		quietly summarize spacetime_data_variance_`i' if spacetime_data_variance_`i' > 0 & spacetime_data_variance_`i' != ., d
		if `r(N)' != 0 quietly replace spacetime_data_variance_`i' = `r(p50)' if spacetime_data_variance_`i' <= 0 | spacetime_data_variance_`i' == .
	}

// make an age group variable instead, because decimals cause some problems
	egen age_group = group(age)

// if there's no data left, then we don't need to make predictions
	quietly count if predictme_${test} == 1 & spacetime_1 != .
	if `r(N)' > 0 {

// save a csv of all the data for this region
		preserve
		keep if predictme_${test} == 1
		keep iso3 age_group year lt_* ln_* *_data_variance_* spacetime_* *_amplitude_* test_${test}
		quietly replace iso3 = "MXR" if iso3 == "MAR"
		sort iso3 age_group year
		outsheet using "$temp_dir/5_gpr_input_${super_region}_${super_region_chunk}_${test}.csv", comma replace

// write a python script to perform GPR
		qui file open gpr_py using "$temp_dir/run_files/5_gpr_${super_region}_${super_region_chunk}_${test}.py", write replace
		global dvs = `"[""' + subinstr("$dv_list"," ", `"",""',.) + `""]"'
		file write gpr_py	"import sys" _n ///
							"sys.path.append('${code_dir}/model/')" _n ///
							"import GPR_mean" _n ///
							"reload(GPR_mean)" _n ///
							"number_submodels = $number_submodels" _n ///
							"infile = '$temp_dir/5_gpr_input_${super_region}_${super_region_chunk}_${test}.csv'" _n ///
							"outfile = '$temp_dir/5_gpr_output_${super_region}_${super_region_chunk}_${test}.csv'" _n ///
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
							"GPR_mean.fit_GPR(infile, outfile, dv_list, scale, number_submodels, test)" _n
		file close gpr_py

// run GPR
		shell /usr/local/epd_py25-4.3.0/bin/python "$temp_dir/run_files/5_gpr_${super_region}_${super_region_chunk}_${test}.py"

// convert the GPR results to Stata
		insheet using "$temp_dir/5_gpr_output_${super_region}_${super_region_chunk}_${test}.csv", comma names clear
		quietly replace iso3 = "MAR" if iso3 == "MXR"
		tempfile gpr_results
		save `gpr_results', replace

// merge GPR results onto the rest of the data
		restore
		merge m:1 iso3 year age_group using `gpr_results', nogen keep(match master)

// fix missing values
		forvalues i = 1 / $number_submodels {
			qui destring gpr_`i'_spacetime_mean, replace force
		}

// put a ceiling on the rates to never exceed the envelope
		forvalues i = 1 / $number_submodels {
			local dv: word `i' of $dv_list
			if "`dv'" == "ln_rate" quietly replace gpr_`i'_spacetime_mean = ln(envelope/pop) if gpr_`i'_spacetime_mean > ln(envelope/pop) & gpr_`i'_spacetime_mean != .
		}
	}

// if there's no data, just put in empty variables here
	else {
		forvalues i = 1 / $number_submodels {
			quietly generate gpr_`i'_spacetime_mean = .
		}
	}

// save final GPR-ed results
	compress
	save "$temp_dir/5_gpr_${super_region}_${super_region_chunk}_${test}.dta", replace

// close the logs
	timer off 1
	timer list 1
	log close
