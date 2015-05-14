/*
Created:	01 May 2011
Updated:	03 May 2012
Purpose: 	Submit CODEm jobs to the cluster
*/

** setup some basic Stata parameters
	set odbcmgr unixodbc
	global dsn CODMOD
	clear
	timer on 1
	set mem 1G
	set matsize 11000
	set more off

    local sex $sex

** protect users of old codem inputs
    if "$end_data_year" == "" {
        global end_data_year 2013
    }
    if "$ref_cov" == "" {
        global ref_cov 0
    }
    if "$ref_cov_list" == "" {
        global ref_cov_list
    }
    if "$pred_cov_list" == "" {
        global pred_cov_list
    }
    if "$ref_cov_path" == "" {
        global ref_cov_path "${base_dir}"
    }
    if "$good_data_c_years_path" == "" {
        // This can stay blank
    }
    if "$merge_vars" == "" {
        // this can stay blank
    }
    if "$constant" == "" {
        global constant 0
    }
    if "$start_year" == "" {
        global start_year 1980
    }
    if "$trend_window_min" == "" {
        global trend_window_min 1
    }
    if "$top_st_rate" == "" {
        global top_st_rate 0
        global top_st_cf 0
        global top_lin_rate 0
        global top_lin_cf 0
    }

** make directories in which to put everything
	capture mkdir "${base_dir}"
	capture mkdir "${base_dir}/logs"
	capture mkdir "${base_dir}/results"
	capture mkdir "${temp_dir}"
	capture mkdir "${temp_dir}/stata_out"
	capture mkdir "${temp_dir}/run_files"

** start logging
	capture log close
	log using "${base_dir}/${model_name}_submission_log.smcl", replace
	// capture ashell hostname
	display "Running model ${model_name} on node `r(o1)'"
	capture rm "${base_dir}/${model_name}_job_ids.txt"
	!echo "\n" > "${base_dir}/${model_name}_job_ids.txt"

** insert constant covariate regression if specified by user
	if $constant == 1 {
		** add new rate model
		global num_rate_models = $num_rate_models + 1
		file open cov_list_rate using "${base_dir}/covariate_selection/selected_covariates_ln_rate.do", write append text
		file write cov_list_rate "global rate_model_${num_rate_models} constant" _n
		file close cov_list_rate
		global rate_model_${num_rate_models} constant

		global number_submodels = $num_cf_models + $num_rate_models

		** recreating all of the globals for the constant rate model
		global submodel_${number_submodels}_name "rate ${number_submodels}"
		global submodel_${number_submodels}_dv "rate"
		global submodel_${number_submodels}_type "covariate"
		global submodel_${number_submodels}_covariates "constant"

		** add new cf model
		global num_cf_models = $num_cf_models + 1
		file open cov_list_cf using "${base_dir}/covariate_selection/selected_covariates_lt_cf.do", write append text
		file write cov_list_cf "global cf_model_${num_cf_models} constant" _n
		file close cov_list_cf

		global number_submodels = $num_cf_models + $num_rate_models

		** recreating all of the globals for the constant rate model
		global submodel_${number_submodels}_name "rate ${number_submodels}"
		global submodel_${number_submodels}_dv "cf"
		global submodel_${number_submodels}_type "covariate"
		global submodel_${number_submodels}_covariates "constant"
	}

** check to make sure everything is filled out correctly (and convert globals to locals for convenience's sake)
	local covariates
	local dv_list
	forvalues i = 1/$number_submodels {
		if missing("${submodel_`i'_name}") {
			display as error "You did not specify a name for submodel `i'."
			exit 197
		}
		if missing("${submodel_`i'_dv}") {
			display as error "You did not specify a dependent variable for submodel `i'. Enter 'rate' or 'cf'."
			exit 197
		}
		if "${submodel_`i'_dv}" == "cf" local submodel_`i'_dv "lt_cf"
		else local submodel_`i'_dv "ln_rate"
		local dv_list `dv_list' `submodel_`i'_dv'
		if missing("${submodel_`i'_type}") {
			display as error "You did not specify a type for submodel `i'. Enter 'covariate' or 'custom'."
			exit 197
		}
		else if "${submodel_`i'_type}" == "covariate" {
			display "Submodel `i', named '${submodel_`i'_name}', is a covariate model with dependent variable ${submodel_`i'_dv} and the following covariates:" _n "${submodel_`i'_covariates}"
			local covariates `covariates' ${submodel_`i'_covariates}
		}
		else if "${submodel_`i'_type}" == "custom" {
			if missing("${submodel_`i'_custom}") {
				display as error "You specified submodel `i' as a custom model but did not include a file path to the custom results."
				exit 197
			}
			capture confirm file "${submodel_`i'_custom}"
			if _rc {
				display as error "You specified the filepath for custom submodel `i' as ${submodel_`i'_custom}, but that file does not exist."
				exit 197
			}
			display "Submodel `i', named '${submodel_`i'_name}', is a custom model with dependent variable ${submodel_`i'_dv} and input file ${submodel_`i'_custom}."
		}
	}

** load in the latest version of the dataset
    insheet using "${base_dir}/input_database_square.csv", comma clear case

** find uncertainty estimates on the dependent variables
	generate cf_sd = sqrt(cf*(1-cf)/sample_size)
	replace cf = . if cf_sd <= 0 | cf_sd == .
	forvalues i = 1/100 {
		qui generate cf_d`i' = rnormal(cf, cf_sd)
		qui replace cf_d`i' = . if cf_d`i' <= 0
		qui generate lt_cf_d`i' = logit(cf_d`i')
	}
	egen lt_cf_sd = rowsd(lt_cf_d*)
	forvalues i = 1/100 {
		qui generate rate_d`i' = rnormal(cf, cf_sd) * (envelope / pop)
		qui replace rate_d`i' = . if rate_d`i' <= 0
		qui generate ln_rate_d`i' = ln(rate_d`i')
	}
	egen ln_rate_sd = rowsd(ln_rate_d*)
	drop ln_rate_d* rate_d* lt_cf_d* cf_d* cf_sd

** Add in trusted data points to use a different psi for some countries
    if "$good_data_c_years_path" != "" {
        merge m:m iso3 year using "$good_data_c_years_path", generate(good_data)
        replace good_data = 0 if good_data != 3
        replace good_data = 1 if good_data == 3
        qui levelsof iso3 if good_data == 1, local(good_countries) c
        count if good_data == 1
        local total_good_points `r(N)'
    }

** create indicators for the predictive validity test
	** first just find which country/year/age combos do and don't have data
		egen age_group = group(age)
		preserve
		collapse (count) ln_rate, by(iso3 year age_group)
		generate has_data = (ln_rate > 0 & ln_rate != .)
		count
		forvalues i = 1/`r(N)' {
			local c = iso3[`i']
			local a = age[`i']
			local y = year[`i']
			local has_data_`c'_`a'_`y' = has_data[`i']
            if "$good_data_c_years_path" != "" {
                local good_data_`c'_`a'_`y' = good_data[`i']
            }
		}
		levelsof iso3, l(isos) c
		levelsof age_group, l(ages) c
		levelsof year, l(years) c
		restore
		count
		forvalues i = 1/`r(N)' {
			local c = iso3[`i']
			local a = age_group[`i']
			local y = year[`i']
			local indices_`c'_`a'_`y' `indices_`c'_`a'_`y'' `i'
		}
	** create some locals needed for looping
		count if ln_rate != .
		local total_datapoints = `r(N)'
		mata iso_list = uniqrows(st_sdata(., "iso3"))
		mata st_numscalar("num_isos", rows(iso_list))
		local num_isos = num_isos
	** then loop through and do knockouts
		global tests insample
		forvalues i = 1 / $holdouts {
			** make a placeholder variable
			generate test_holdout_`i' = 0
			local ho = 1
            ** local to make sure enough good data was knocked out
            local bad_done = 0
			** create a list to draw patterns of missingness from of all countries in random order
			mata in_isos = jumble(iso_list)
			** create a list to apply drops to of all countries in random order
			mata out_isos = jumble(iso_list)
			** loop through the countries to draw patterns of missingness from
			forvalues j = 1/`num_isos' {
				** store the current iso3 codes
				mata st_local("cur_iso_in", in_isos[`j'])
				mata st_local("cur_iso_out", out_isos[`j'])
				** loop through each age/year combo
				foreach a of local ages {
					foreach y of local years {
                              // I don't know why this happens
                                if "`has_data_`cur_iso_in'_`a'_`y''" == "" {

                                    local has_data_`cur_iso_in'_`a'_`y' 0

                                }
                        **Make sure we're holding out good countries as well
                        if "$good_data_c_years_path" != "" {
                            di "bad done: `bad_done'"
                            di "has_data_`cur_iso_in'_`a'_`y': `has_data_`cur_iso_in'_`a'_`y''"
                            if `bad_done' == 0 {
                                ** if there is no data for the input iso3 for this year/age, set it as test data for the output iso3

                                if `has_data_`cur_iso_in'_`a'_`y'' == 0 {
                                    foreach n of local indices_`cur_iso_out'_`a'_`y' {
                                        qui replace test_holdout_`i' = `ho' in `n'
                                    }
                                }
                            }
                            else {
                                ** only knock out good data
                                if `has_data_`cur_iso_in'_`a'_`y'' == 0 & `good_data_`cur_iso_in'_`a'_`y'' == 1 {
                                    foreach n of local indices_`cur_iso_out'_`a'_`y' {
                                        qui replace test_holdout_`i' = `ho' in `n'
                                    }
                                }
                            }
                        }
                        else {

                            ** if there is no data for the input iso3 for this year/age, set it as test data for the output iso3
                            if `has_data_`cur_iso_in'_`a'_`y'' == 0 {
                                foreach n of local indices_`cur_iso_out'_`a'_`y' {
                                    qui replace test_holdout_`i' = `ho' in `n'
                                }
                            }
                        }
					}
                }
				** figure out how much of the data has been dropped, stopping after the desired holdout level is reached
				if `ho' == 1 {
                    if "$good_data_c_years_path" != "" {
                        qui count if ln_rate != . & test_holdout_`i' == 1
                        local knocked_out = `r(N)'
                        qui count if ln_rate != . & test_holdout_`i' == 1 & good_data == 1
                        local knocked_out_good = `r(N)'

                        if (`knocked_out'/`total_datapoints' > ${submodel_holdout}) & (`knocked_out_good'/`total_good_points' > ${submodel_holdout}) {
                            local ho = 2
                            local bad_done = 0
                        }
                        else if (`knocked_out'/`total_datapoints' > ${submodel_holdout}) & (`knocked_out_good'/`total_good_points' <= ${submodel_holdout}) {
                            local bad_done = 1
                        }
                    }
                    else {
                        qui count if ln_rate != . & test_holdout_`i' == 1
                        local knocked_out = `r(N)'
                        if (`knocked_out'/`total_datapoints' > ${submodel_holdout})  {
                            local ho = 2
                        }
                    }
				}
				else if `ho' == 2 {
                    if "$good_data_c_years_path" != "" {
                        qui count if ln_rate != . & test_holdout_`i' == 2
                        local knocked_out = `r(N)'
                        qui count if ln_rate != . & test_holdout_`i' == 2 & good_data == 1
                        local knocked_out_good = `r(N)'
                        if (`knocked_out'/`total_datapoints' > ${submodel_holdout}) & (`knocked_out_good'/`total_good_points' > ${submodel_holdout}) {
                            local bad_done = 0
                            continue, break
                        }
                        else if (`knocked_out'/`total_datapoints' > ${submodel_holdout}) & (`knocked_out_good'/`total_good_points' <= ${submodel_holdout}) {
                            local bad_done = 1
                        }
                    }
                    else {
                        qui count if ln_rate != . & test_holdout_`i' == 2
                        local knocked_out = `r(N)'
                        if (`knocked_out'/`total_datapoints' > ${submodel_holdout})  {
                            continue, break
                        }
                    }
				}
			}
			** add to the list of tests
			global tests $tests holdout_`i'
		}
		** identify which observations actually need predictions made for them (i.e. any country-age with knockouts - any training data outside of a country-age with test data does not actually need spacetime or GPR predictions, because it's never used - this will save us time)
		forvalues i = 1 / $holdouts {
			generate predictme_tmp = (ln_rate != . & test_holdout_`i' > 0)
			bysort iso3 age: egen predictme_holdout_`i' = max(predictme_tmp)
			drop predictme_tmp
		}
		generate predictme_insample = 1
	** add some other necessary variables
		generate test_insample = 0
		generate n = _n
		generate cause = "$cause"
		drop age_group

** save the database as an input to subsequent steps
	save "${base_dir}/results/input_database.dta", replace
	save "${temp_dir}/input_database.dta", replace

** make a lists of countries for later looping
	global super_regions = subinstr("$super_regions",","," ",.)
	local chunks_per_sr
	foreach sr of global super_regions {
		levelsof iso3 if super_region==`sr', l(isos_`sr') c
		local number_isos: list sizeof isos_`sr'
		local number_isos: list sizeof isos_`sr'
		local number_chunks_`sr' =  ceil(`number_isos' / 10)

		local chunks_per_sr `chunks_per_sr' `number_chunks_`sr''
	}

// ** we'll break step 1 into multiple parts, 1 job per submodel = many fast jobs
** we'll break step 1 into multiple parts, groups of 10 submodels
	local number_submodel_chunks = ceil($number_submodels / 10)
	local sge "-o /dev/null -e /dev/null"

** 1. run the linear models with region random effects (priors for spacetime) for each holdout
	foreach t of global tests {
		forvalues c = 1 / `number_submodel_chunks' {
            local start_submodel = (`c'-1)*10 + 1
			local end_submodel = min(`c'*10, $number_submodels)
			file open linear_do using "${temp_dir}/run_files/1_linear_region_`c'_`t'.do", write replace text
			file write linear_do 	"global model_name $model_name" _n ///
									"global cause $cause" _n ///
									"global number_submodels $number_submodels" _n ///
									"global linear_floor $linear_floor" _n ///
									"global test `t'" _n ///
									`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
									`"global code_dir "$code_dir""' _n ///
                                    `"global ref_cov "$ref_cov""' _n ///
                                    `"global ref_cov_list "$ref_cov_list""' _n ///
                                    `"global pred_cov_list "$pred_cov_list""' _n ///
                                    `"global ref_cov_path "$ref_cov_path""' _n ///
                                    `"global merge_vars $merge_vars"' _n ///
                                    "global good_data_c_years_path $good_data_c_years_path" _n ///
									"global submodel_chunk `c'" _n ///
									"global start_submodel `start_submodel'" _n ///
									"global end_submodel `end_submodel'" _n ///
									"global super_regions $super_regions" _n
			forvalues i = `start_submodel' / `end_submodel' {
				file write linear_do 	"global type_`i' ${submodel_`i'_type}" _n ///
										"global name_`i' ${submodel_`i'_name}" _n ///
										"global dv_`i' `submodel_`i'_dv'" _n ///
										"global custom_`i' ${submodel_`i'_custom}" _n ///
										"global covariates_`i' ${submodel_`i'_covariates}" _n
			}
			file write linear_do 	`"do "${code_dir}/model/1_linear_region_re.do""' _n
			file close linear_do
			file open linear_sh using "${temp_dir}/run_files/1_linear_region_`c'_`t'.sh", write replace text
			noisily di "${temp_dir}/run_files/1_linear_region_`c'_`t'.sh"
			file write linear_sh 	"#!/bin/sh" _n ///
									"#$ -S /bin/sh" _n ///
									"export STATATMP=/tmp" _n ///
									"export HOME=/dev/null" _n ///
									`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/1_linear_region_`c'_`t'.do" > ${temp_dir}/stata_out/stata_out_"' ///
									`"\$"' ///
									`"{JOB_ID}"' _n
			file close linear_sh
		}
	}

// 2. run the linear models with country random effects (for the linear only models) for each holdout
	foreach t of global tests {
		forvalues c = 1 / `number_submodel_chunks' {
            local start_submodel = (`c'-1)*10 + 1
			local end_submodel = min(`c'*10, $number_submodels)
			qui file open linear_do using "${temp_dir}/run_files/2_linear_country_`c'_`t'.do", write replace text
			file write linear_do 	"global model_name $model_name" _n ///
									"global cause $cause" _n ///
									`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
									`"global code_dir "$code_dir""' _n ///
                                    `"global ref_cov "$ref_cov""' _n ///
                                    `"global ref_cov_list "$ref_cov_list""' _n ///
                                    `"global pred_cov_list "$pred_cov_list""' _n ///
                                    `"global ref_cov_path "$ref_cov_path""' _n ///
                                    `"global merge_vars $merge_vars"' _n ///
                                    "global good_data_c_years_path $good_data_c_years_path" _n ///
									"global number_submodels $number_submodels" _n ///
									"global linear_floor $linear_floor" _n ///
									"global test `t'" _n ///
									"global submodel_chunk `c'" _n ///
									"global start_submodel `start_submodel'" _n ///
									"global end_submodel `end_submodel'" _n ///
									"global super_regions $super_regions" _n
			forvalues i = `start_submodel' / `end_submodel' {
				file write linear_do 	"global type_`i' ${submodel_`i'_type}" _n ///
										"global name_`i' ${submodel_`i'_name}" _n ///
										"global dv_`i' `submodel_`i'_dv'" _n ///
										"global custom_`i' ${submodel_`i'_custom}" _n ///
										"global covariates_`i' ${submodel_`i'_covariates}" _n
			}
			file write linear_do 	`"do "${code_dir}/model/2_linear_country_re.do""' _n
			file close linear_do
			qui file open linear_sh using "${temp_dir}/run_files/2_linear_country_`c'_`t'.sh", write replace text
			file write linear_sh 	"#!/bin/sh" _n ///
									"#$ -S /bin/sh" _n ///
									"export LD_LIBRARY_PATH=/usr/lib64" _n ///
									"export STATATMP=/tmp" _n ///
									"export HOME=/dev/null" _n ///
									`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/2_linear_country_`c'_`t'.do" > ${temp_dir}/stata_out/stata_out_"' ///
									`"\$"' ///
									`"{JOB_ID}"' _n
			file close linear_sh
		}
	}

// 3. run the second stage (aka "spacetime") for all submodels by super region (split into 10 country chunks), grouped together by predictive validity test
	foreach t of global tests {
		local spacetime_`t'
		foreach sr of global super_regions {
			forvalues i = 1/`number_chunks_`sr'' {
				local chunk_isos: piece `i' 40 of "`isos_`sr''"
				qui file open spacetime_do using "${temp_dir}/run_files/3_spacetime_`sr'_`i'_`t'.do", write replace text
				file write spacetime_do	"global cause $cause" _n ///
										"global model_name $model_name" _n ///
										`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
										`"global code_dir "$code_dir""' _n ///
                                        "global good_data_c_years_path $good_data_c_years_path" _n ///
										"global isos `chunk_isos'" _n ///
										"global super_region `sr'" _n ///
										"global super_region_chunk `i'" _n ///
										"global test `t'" _n ///
										"global number_submodel_chunks `number_submodel_chunks'" _n ///
										"global number_submodels $number_submodels" _n ///
										"global start_year $start_year" _n ///
                                        "global end_data_year $end_data_year" _n ///
                                        "global start_age $start_age" _n ///
                                        "global end_age $end_age" _n ///
										"global zeta $zeta" _n ///
										"global zeta_no_data $zeta_no_data" _n ///
										"global lambda $lambda" _n ///
										"global lambda_no_data $lambda_no_data" _n ///
										"global omega $omega" _n ///
										"global dv_list `dv_list'" _n ///
										`"do "${code_dir}/model/3_spacetime.do""' _n
				file close spacetime_do
				qui file open spacetime_sh using "${temp_dir}/run_files/3_spacetime_`sr'_`i'_`t'.sh", write replace text
				file write spacetime_sh	"#!/bin/sh" _n ///
										"#$ -S /bin/sh" _n ///
										"export STATATMP=/tmp" _n ///
										"export HOME=/dev/null" _n ///
										`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/3_spacetime_`sr'_`i'_`t'.do" > ${temp_dir}/stata_out/stata_out_"' ///
										`"\$"' ///
										`"{JOB_ID}"' _n
				file close spacetime_sh
			}
		}
	}

// 4. once all spacetime models have run, find the GPR parameters for each predictive validity test
	foreach t of global tests {
		qui file open params_do using "${temp_dir}/run_files/4_find_GPR_params_`t'.do", write replace text
		file write params_do 	"global cause $cause" _n ///
								"global model_name $model_name" _n ///
								`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
								`"global code_dir "$code_dir""' _n ///
								"global super_regions $super_regions" _n ///
								"global test `t'" _n ///
								"global number_submodels $number_submodels" _n ///
								"global chunks_per_sr `chunks_per_sr'" _n ///
								"global dv_list `dv_list'" _n ///
								`"do "${code_dir}/model/4_find_GPR_parameters.do""' _n
		file close params_do
		qui file open params_sh using "${temp_dir}/run_files/4_find_GPR_params_`t'.sh", write replace text
		file write params_sh 	"#!/bin/sh" _n ///
								"#$ -S /bin/sh" _n ///
								"export STATATMP=/tmp" _n ///
								"export HOME=/dev/null" _n ///
								`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/4_find_GPR_params_`t'.do" > ${temp_dir}/stata_out/stata_out_"' ///
								`"\$"' ///
								`"{JOB_ID}"' _n
		file close params_sh
	}

// 5. run GPR by predictive validity test and super region chunks
	foreach t of global tests {
		foreach sr of global super_regions {
			forvalues i = 1/`number_chunks_`sr'' {
				local chunk_isos: piece `i' 40 of "`isos_`sr''"
				qui file open gpr_do using "${temp_dir}/run_files/5_gpr_`sr'_`i'_`t'.do", write replace text
				file write gpr_do	"global cause $cause" _n ///
									"global model_name $model_name" _n ///
									`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
									`"global code_dir "$code_dir""' _n ///
									"global isos `chunk_isos'" _n ///
									"global super_region `sr'" _n ///
									"global super_region_chunk `i'" _n ///
									"global dv_list `dv_list'" _n ///
									"global test `t'" _n ///
									"global scale $scale" _n ///
									"global number_submodels $number_submodels" _n ///
									`"do "${code_dir}/model/5_GPR_mean.do""' _n
				file close gpr_do
				qui file open gpr_sh using "${temp_dir}/run_files/5_gpr_`sr'_`i'_`t'.sh", write replace text
				file write gpr_sh 	"#!/bin/sh" _n ///
									"#$ -S /bin/sh" _n ///
									"export STATATMP=/tmp" _n ///
									"export HOME=/dev/null" _n ///
									`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/5_gpr_`sr'_`i'_`t'.do" > ${temp_dir}/stata_out/stata_out_"' ///
									`"\$"' ///
									`"{JOB_ID}"' _n
				file close gpr_sh
			}
		}
	}

// 6. calculate predictive validity metrics on the submodels
	foreach t of global tests {
		qui file open pv_do using "${temp_dir}/run_files/6_submodel_pv_`t'.do", write replace text
		file write pv_do	"global cause $cause" _n ///
							"global model_name $model_name" _n ///
							`"global base_dir "$base_dir""' _n ///
                             `"global temp_dir "$temp_dir""' _n ///
							`"global code_dir "$code_dir""' _n ///
                            "global trend_window_min $trend_window_min" _n ///
                            "global trend_window $trend_window" _n ///
							"global trend_method $trend_method" _n ///
                            "global good_data_c_years_path $good_data_c_years_path" _n ///
                            `"global good_countries "`good_countries'""' _n ///
							"global super_regions $super_regions" _n ///
							"global number_submodel_chunks `number_submodel_chunks'" _n ///
							"global test `t'" _n ///
							"global chunks_per_sr `chunks_per_sr'" _n ///
							"global dv_list `dv_list'" _n ///
							"global number_submodels $number_submodels" _n ///
							`"do "${code_dir}/model/6_submodel_pv.do""' _n
		file close pv_do
		qui file open pv_sh using "${temp_dir}/run_files/6_submodel_pv_`t'.sh", write replace text
		file write pv_sh 	"#!/bin/sh" _n ///
							"#$ -S /bin/sh" _n ///
							"export STATATMP=/tmp" _n ///
							"export HOME=/dev/null" _n ///
							`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/6_submodel_pv_`t'.do" > ${temp_dir}/stata_out/stata_out_"' ///
							`"\$"' ///
							`"{JOB_ID}"' _n
		file close pv_sh
	}

// 7. calculate the submodel ranks
	qui file open weight_do using "${temp_dir}/run_files/7_submodel_ranks.do", write replace text
	file write weight_do	"global cause $cause" _n ///
							"global model_name $model_name" _n ///
							`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
							`"global code_dir "$code_dir""' _n ///
                            "global trend_window $trend_window" _n ///
							"global trend_method $trend_method" _n ///
                            "global good_data_c_years_path $good_data_c_years_path" _n ///
                            `"global good_countries "`good_countries'""' _n ///
							"global tests $tests" _n ///
							"global dv_list `dv_list'" _n ///
							"global psi_max $psi_max" _n ///
							"global psi_min ${psi_min}" _n ///
							"global psi_int ${psi_int}" _n ///
							"global number_submodels $number_submodels" _n
	forvalues i = 1/$number_submodels {
		file write weight_do 	"global type_`i' ${submodel_`i'_type}" _n ///
								"global name_`i' ${submodel_`i'_name}" _n ///
								"global dv_`i' ${submodel_`i'_dv}" _n ///
								"global covariates_`i' ${submodel_`i'_covariates}" _n
	}
	file write weight_do 	`"do "${code_dir}/model/7_submodel_ranks.do""' _n
	file close weight_do
	qui file open weight_sh using "${temp_dir}/run_files/7_submodel_ranks.sh", write replace text
	file write weight_sh 	"#!/bin/sh" _n ///
							"#$ -S /bin/sh" _n ///
							"export LD_LIBRARY_PATH=/usr/lib64" _n ///
							"export STATATMP=/tmp" _n ///
							"export HOME=/dev/null" _n ///
							`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/7_submodel_ranks.do" > ${temp_dir}/stata_out/stata_out_"' ///
							`"\$"' ///
							`"{JOB_ID}"' _n
	file close weight_sh

// 8. calculate PV of the pooled model
	foreach t of global tests {
		qui file open pv_do using "${temp_dir}/run_files/8_ensemble_pv_`t'.do", write replace text
		file write pv_do	"global cause $cause" _n ///
							"global model_name $model_name" _n ///
							`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
							`"global code_dir "$code_dir""' _n ///
                            "global trend_window_min $trend_window_min" _n ///
                            "global trend_window $trend_window" _n ///
							"global trend_method $trend_method" _n ///
                            "global good_data_c_years_path $good_data_c_years_path" _n ///
                            `"global good_countries "`good_countries'""' _n ///
							"global super_regions $super_regions" _n ///
							"global test `t'" _n ///
							"global chunks_per_sr `chunks_per_sr'" _n ///
							"global dv_list `dv_list'" _n ///
							"global psi_max $psi_max" _n ///
							"global psi_min ${psi_min}" _n ///
							"global psi_int ${psi_int}" _n ///
							"global number_submodel_chunks `number_submodel_chunks'" _n ///
							"global number_submodels $number_submodels" _n ///
							`"do "${code_dir}/model/8_ensemble_pv.do""' _n
		file close pv_do
		qui file open pv_sh using "${temp_dir}/run_files/8_ensemble_pv_`t'.sh", write replace text
		file write pv_sh 	"#!/bin/sh" _n ///
							"#$ -S /bin/sh" _n ///
							"export STATATMP=/tmp" _n ///
							"export HOME=/dev/null" _n ///
							`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/8_ensemble_pv_`t'.do" > ${temp_dir}/stata_out/stata_out_"' ///
							`"\$"' ///
							`"{JOB_ID}"' _n
		file close pv_sh
	}

// 9. combine metrics of PV on the ensemble models across holdouts
	qui file open pool_do using "${temp_dir}/run_files/9_ensemble_ranks.do", write replace text
	file write pool_do	"global cause $cause" _n ///
						"global model_name $model_name" _n ///
						`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
						`"global code_dir "$code_dir""' _n ///
                        "global good_data_c_years_path $good_data_c_years_path" _n ///
                        `"global good_countries "`good_countries'""' _n ///
						"global number_submodels $number_submodels" _n ///
						"global tests $tests" _n
	file write pool_do 	`"do "${code_dir}/model/9_ensemble_ranks.do""' _n
	file close pool_do
	qui file open pool_sh using "${temp_dir}/run_files/9_ensemble_ranks.sh", write replace text
	file write pool_sh 	"#!/bin/sh" _n ///
						"#$ -S /bin/sh" _n ///
						"export STATATMP=/tmp" _n ///
						"export HOME=/dev/null" _n ///
						`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/9_ensemble_ranks.do" > ${temp_dir}/stata_out/stata_out_"' ///
						`"\$"' ///
						`"{JOB_ID}"' _n
	file close pool_sh

// 10. run GPR to find CI for the pooled models
	foreach t of global tests {
		foreach sr of global super_regions {
			forvalues i = 1/`number_chunks_`sr'' {
				local chunk_isos: piece `i' 40 of "`isos_`sr''"
				qui file open gpr_do using "${temp_dir}/run_files/10_gpr_draws_`sr'_`i'_`t'.do", write replace text
				file write gpr_do	"global cause $cause" _n ///
									"global model_name $model_name" _n ///
									`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
									`"global code_dir "$code_dir""' _n ///
									"global isos `chunk_isos'" _n ///
									"global super_region `sr'" _n ///
									"global super_region_chunk `i'" _n ///
									"global dv_list `dv_list'" _n ///
									"global test `t'" _n ///
									"global psi_int ${psi_int}" _n ///
									"global scale $scale" _n ///
									"global number_submodels $number_submodels" _n ///
									"global top_st_rate ${top_st_rate}" _n ///
									"global top_st_cf ${top_st_cf}" _n ///
									"global top_lin_rate ${top_lin_rate}" _n ///
									"global top_lin_cf ${top_lin_cf}" _n ///
									`"do "${code_dir}/model/10_GPR_draws.do""' _n
				file close gpr_do
				qui file open gpr_sh using "${temp_dir}/run_files/10_gpr_draws_`sr'_`i'_`t'.sh", write replace text
				file write gpr_sh 	"#!/bin/sh" _n ///
									"#$ -S /bin/sh" _n ///
									"export STATATMP=/tmp" _n ///
									"export HOME=/dev/null" _n ///
									`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/10_gpr_draws_`sr'_`i'_`t'.do" > ${temp_dir}/stata_out/stata_out_"' ///
									`"\$"' ///
									`"{JOB_ID}"' _n
				file close gpr_sh
			}
		}
	}

// 11. make draws from the linear model
	foreach t of global tests {
		qui file open linear_do using "${temp_dir}/run_files/11_linear_draws_`t'.do", write replace text
		file write linear_do 	"global cause $cause" _n ///
								"global model_name $model_name" _n ///
								`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
								`"global code_dir "$code_dir""' _n ///
								"global dv_list `dv_list'" _n ///
								"global test `t'" _n ///
								"global number_submodels $number_submodels" _n ///
								"global number_submodel_chunks `number_submodel_chunks'" _n ///
								"global psi_int ${psi_int}" _n ///
								"global top_st_rate ${top_st_rate}" _n ///
								"global top_st_cf ${top_st_cf}" _n ///
								"global top_lin_rate ${top_lin_rate}" _n ///
								"global top_lin_cf ${top_lin_cf}" _n ///
								"global linear_floor $linear_floor" _n ///
								`"do "${code_dir}/model/11_linear_draws.do""' _n
		file close linear_do
		qui file open linear_sh using "${temp_dir}/run_files/11_linear_draws_`t'.sh", write replace text
		file write linear_sh 	"#!/bin/sh" _n ///
								"#$ -S /bin/sh" _n ///
								"export STATATMP=/tmp" _n ///
								"export HOME=/dev/null" _n ///
								`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/11_linear_draws_`t'.do" > ${temp_dir}/stata_out/stata_out_"' ///
								`"\$"' ///
								`"$JOB_ID}"' _n
		file close linear_sh
	}

// 12. put all the in-sample prediction files together
	qui file open combine_do using "${temp_dir}/run_files/12_put_predictions_together.do", write replace text
	file write combine_do 	"global cause $cause" _n ///
							"global model_name $model_name" _n ///
							`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
							`"global code_dir "$code_dir""' _n ///
							"global super_regions $super_regions" _n ///
							"global chunks_per_sr `chunks_per_sr'" _n ///
							"global top_st_rate ${top_st_rate}" _n ///
							"global top_st_cf ${top_st_cf}" _n ///
							"global top_lin_rate ${top_lin_rate}" _n ///
							"global top_lin_cf ${top_lin_cf}" _n ///
							`"do "${code_dir}/model/12_put_predictions_together.do""' _n
	file close combine_do
	qui file open combine_sh using "${temp_dir}/run_files/12_put_predictions_together.sh", write replace text
	file write combine_sh 	"#!/bin/sh" _n ///
							"#$ -S /bin/sh" _n ///
							"export LD_LIBRARY_PATH=/usr/lib64" _n ///
							"export STATATMP=/tmp" _n ///
							"export HOME=/dev/null" _n ///
							`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/12_put_predictions_together.do" > ${temp_dir}/stata_out/stata_out_"' ///
							`"\$"' ///
							`"{JOB_ID}"' _n
	file close combine_sh

// 13. calculate totals by country/region/super region
	foreach sr of global super_regions {
		qui file open totals_do using "${temp_dir}/run_files/13_calculate_totals_`sr'.do", write replace text
		file write totals_do 	"global cause $cause" _n ///
								"global model_name $model_name" _n ///
								`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
								`"global code_dir "$code_dir""' _n ///
								"global number_submodels $number_submodels" _n ///
								"global dv_list `dv_list'" _n ///
								"global super_region `sr'" _n ///
								"global super_region_chunks `number_chunks_`sr''" _n ///
								`"do "${code_dir}/model/13_calculate_totals.do""' _n
		file close totals_do
		qui file open totals_sh using "${temp_dir}/run_files/13_calculate_totals_`sr'.sh", write replace text
		file write totals_sh 	"#!/bin/sh" _n ///
								"#$ -S /bin/sh" _n ///
								"export STATATMP=/tmp" _n ///
								"export HOME=/dev/null" _n ///
								`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/13_calculate_totals_`sr'.do" > ${temp_dir}/stata_out/stata_out_"' ///
								`"\$"' ///
								`"{JOB_ID}"' _n
		file close totals_sh
	}

// 14. output final predictions
	qui file open output_do using "${temp_dir}/run_files/14_output_totals.do", write replace text
	file write output_do 	"global cause $cause" _n ///
							"global sex `sex'" _n ///
							"global model_name $model_name" _n ///
							`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
							`"global code_dir "$code_dir""' _n ///
							"global super_regions $super_regions" _n ///
							"global top_st_rate ${top_st_rate}" _n ///
							"global top_st_cf ${top_st_cf}" _n ///
							"global top_lin_rate ${top_lin_rate}" _n ///
							"global top_lin_cf ${top_lin_cf}" _n ///
							`"do "${code_dir}/model/14_output_totals.do""' _n
	file close output_do
	qui file open output_sh using "${temp_dir}/run_files/14_output_totals.sh", write replace text
	file write output_sh 	"#!/bin/sh" _n ///
							"#$ -S /bin/sh" _n ///
							"export STATATMP=/tmp" _n ///
							"export HOME=/dev/null" _n ///
							`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/14_output_totals.do" > ${temp_dir}/stata_out/stata_out_"' ///
							`"\$"' ///
							`"{JOB_ID}"' _n
	file close output_sh

// 15. find final metrics of fit (coverage, RMSE, trend) for each holdout
	foreach t of global tests {
		qui file open coverage_do using "${temp_dir}/run_files/15_coverage_`t'.do", write replace text
		file write coverage_do 	"global cause $cause" _n ///
								"global model_name $model_name" _n ///
								`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
								`"global code_dir "$code_dir""' _n ///
                                "global trend_window $trend_window" _n ///
                                "global trend_method $trend_method" _n ///
                                "global trend_window_min $trend_window_min" _n ///
                                "global good_data_c_years_path $good_data_c_years_path" _n ///
                                `"global good_countries "`good_countries'""' _n ///
                                "global psi_int ${psi_int}" _n ///
								"global test `t'" _n ///
								"global super_regions $super_regions" _n ///
								"global chunks_per_sr `chunks_per_sr'" _n ///
								"global top_st_rate ${top_st_rate}" _n ///
								"global top_st_cf ${top_st_cf}" _n ///
								"global top_lin_rate ${top_lin_rate}" _n ///
								"global top_lin_cf ${top_lin_cf}" _n ///
								`"do "${code_dir}/model/15_ensemble_coverage.do""' _n
		file close coverage_do
		qui file open coverage_sh using "${temp_dir}/run_files/15_coverage_`t'.sh", write replace text
		file write coverage_sh 	"#!/bin/sh" _n ///
								"#$ -S /bin/sh" _n ///
								"export STATATMP=/tmp" _n ///
								"export HOME=/dev/null" _n ///
								`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/15_coverage_`t'.do" > ${temp_dir}/stata_out/stata_out_"' ///
								`"\$"' ///
								`"{JOB_ID}"' _n
		file close coverage_sh
	}

// 16. combine all the metrics of fit for the final ensemble
	qui file open finalpv_do using "${temp_dir}/run_files/16_final_pv.do", write replace text
	file write finalpv_do 	"global cause $cause" _n ///
							"global model_name $model_name" _n ///
							`"global base_dir "$base_dir""' _n ///
                                    `"global temp_dir "$temp_dir""' _n ///
							`"global code_dir "$code_dir""' _n ///
							"global tests $tests" _n ///
							"global top_st_rate ${top_st_rate}" _n ///
							"global top_st_cf ${top_st_cf}" _n ///
							"global top_lin_rate ${top_lin_rate}" _n ///
							"global top_lin_cf ${top_lin_cf}" _n ///
							`"do "${code_dir}/model/16_combine_final_pv.do""' _n
	file close finalpv_do
	qui file open finalpv_sh using "${temp_dir}/run_files/16_final_pv.sh", write replace text
	file write finalpv_sh 	"#!/bin/sh" _n ///
							"#$ -S /bin/sh" _n ///
							"export HOME=/dev/null" _n ///
							`"/usr/local/stata11/stata-mp < "${temp_dir}/run_files/16_final_pv.do" > ${temp_dir}/stata_out/stata_out_"' ///
							`"\$"' ///
							`"{JOB_ID}"' _n
	file close finalpv_sh

// actually submit the jobs
	foreach t of global tests {

		local n = 0
		quietly {
			clear
			set obs 1
			gen model = ""
			tempfile lc_holds
			save `lc_holds', replace
			tempfile lr_holds
			save `lr_holds', replace
		}
        di `number_submodel_chunks'
		forvalues c = 1 / `number_submodel_chunks' {

			// submit 1
			!/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -l mem_free=4G -N lr_${cause}_${model_name}_`c'_`t' "${temp_dir}/run_files/1_linear_region_`c'_`t'.sh" >> "${base_dir}/${model_name}_job_ids.txt"
			qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"

			local n = `n' + 1
			quietly {
				use `lr_holds', clear
				set obs `n'
				replace model = "lr_${cause}_${model_name}_`c'_`t'" in `n'
				tempfile lr_holds
				save `lr_holds', replace
			}

			// submit 2
			qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge'  -l mem_free=4G -N lc_${cause}_${model_name}_`c'_`t' "${temp_dir}/run_files/2_linear_country_`c'_`t'.sh" >> "${base_dir}/${model_name}_job_ids.txt"
			qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"

			quietly {
				use `lc_holds', clear
				set obs `n'
				replace model = "lc_${cause}_${model_name}_`c'_`t'" in `n'
				tempfile lc_holds
				save `lc_holds', replace
			}
		}

		use `lc_holds', clear
		quietly levelsof model, clean local(lc_holds_list) separate(,)

		use `lr_holds', clear
		quietly levelsof model, clean local(lr_holds_list) separate(,)

		local n = 0
		quietly {
			clear
			set obs 1
			gen model = ""
			tempfile st_holds
			save `st_holds', replace
		}

		foreach sr of global super_regions {
			forvalues i = 1/`number_chunks_`sr'' {
				local chunk_isos: piece `i' 40 of "`isos_`sr''"

				// submit 3
				qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -l mem_free=2G -N s_${cause}_${model_name}_`sr'_`i'_`t' -hold_jid `lr_holds_list' "${temp_dir}/run_files/3_spacetime_`sr'_`i'_`t'.sh"  >> "${base_dir}/${model_name}_job_ids.txt"
				qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"

				local n = `n' + 1
				quietly {
					use `st_holds', clear
					set obs `n'
					replace model = "s_${cause}_${model_name}_`sr'_`i'_`t'" in `n'
					tempfile st_holds
					save `st_holds', replace
				}
			}
		}

		use `st_holds', clear
		quietly levelsof model, clean local(st_holds_list) separate(,)

		local n = 0
		quietly {
			clear
			set obs 1
			gen model = ""
			tempfile gpr_holds
			save `gpr_holds', replace
		}

		// submit 4
		qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub -N p_${cause}_${model_name}_`t' -hold_jid `st_holds_list', `sge' -pe multi_slot 4 -l mem_free=2G "${temp_dir}/run_files/4_find_GPR_params_`t'.sh" >> "${base_dir}/${model_name}_job_ids.txt"
		qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"
		foreach sr of global super_regions {
			forvalues i = 1/`number_chunks_`sr'' {
				local chunk_isos: piece `i' 40 of "`isos_`sr''"

				// submit 5
				qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub  `sge' -pe multi_slot 4 -l mem_free=3G -N g_${cause}_${model_name}_`sr'_`i'_`t' -hold_jid p_${cause}_${model_name}_`t' "${temp_dir}/run_files/5_gpr_`sr'_`i'_`t'.sh" >> "${base_dir}/${model_name}_job_ids.txt"
				qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"

				local n = `n' + 1
				quietly {
					use `gpr_holds', clear
					set obs `n'
					replace model = "g_${cause}_${model_name}_`sr'_`i'_`t'" in `n'
					tempfile gpr_holds
					save `gpr_holds', replace
				}
			}
		}

		use `gpr_holds', clear
		quietly levelsof model, clean local(gpr_holds_list) separate(,)

		// submit 6
		qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -pe multi_slot 4 -l mem_free=5G -N pvs_${cause}_${model_name}_`t' -hold_jid `gpr_holds_list'',`lc_holds_list'' "${temp_dir}/run_files/6_submodel_pv_`t'.sh" >> "${base_dir}/${model_name}_job_ids.txt"
		qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"
		if missing("`submodel_pv_waits'") local submodel_pv_waits pvs_${cause}_${model_name}_`t'
		else local submodel_pv_waits `submodel_pv_waits',pvs_${cause}_${model_name}_`t'
	}

	// submit 7
	qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub -l mem_free=1G `sge' -pe multi_slot 4 -N rs_${cause}_${model_name} -hold_jid `submodel_pv_waits' "${temp_dir}/run_files/7_submodel_ranks.sh" >> "${base_dir}/${model_name}_job_ids.txt"
	qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"

	foreach t of global tests {

		// submit 8
		qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -pe multi_slot 4 -l mem_free=4G -N pve_${cause}_${model_name}_`t' -hold_jid rs_${cause}_${model_name} "${temp_dir}/run_files/8_ensemble_pv_`t'.sh" >> "${base_dir}/${model_name}_job_ids.txt"
		qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"
		if missing("`pv_ensemble_waits'") local pv_ensemble_waits pve_${cause}_${model_name}_`t'
		else local pv_ensemble_waits `pv_ensemble_waits',pve_${cause}_${model_name}_`t'
	}

	sleep 1000

	// submit 9
	qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -pe multi_slot 4 -l mem_free=1G -N re_${cause}_${model_name} -hold_jid `pv_ensemble_waits' "${temp_dir}/run_files/9_ensemble_ranks.sh" >> "${base_dir}/${model_name}_job_ids.txt"
	qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"

	local n = 0
	quietly {
		clear
		set obs 1
		gen model = ""
		tempfile holds
		save `holds', replace
	}

	foreach t of global tests {
		foreach sr of global super_regions {
			forvalues i = 1/`number_chunks_`sr'' {
				local chunk_isos: piece `i' 40 of "`isos_`sr''"

				// submit 10
				qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -pe multi_slot 4 -l mem_free=4G -N gd_${cause}_${model_name}_`sr'_`i'_`t' -hold_jid re_${cause}_${model_name} "${temp_dir}/run_files/10_gpr_draws_`sr'_`i'_`t'.sh" >> "${base_dir}/${model_name}_job_ids.txt"
				qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"
				if missing("`gpr_d_`t''") local gpr_d_`t' gd_${cause}_${model_name}_`sr'_`i'_`t'
				else local gpr_d_`t' `gpr_d_`t'',gd_${cause}_${model_name}_`sr'_`i'_`t'

				local n = `n' + 1
				quietly {
					use `holds', clear
					set obs `n'
					replace model = "gd_${cause}_${model_name}_`sr'_`i'_`t'" in `n'
					tempfile holds
					save `holds', replace
				}
			}
		}

		// submit 11
		qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -pe multi_slot 4 -l mem_free=5G -N ld_${cause}_${model_name}_`t' -hold_jid re_${cause}_${model_name} "${temp_dir}/run_files/11_linear_draws_`t'.sh" >> "${base_dir}/${model_name}_job_ids.txt"
		qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"

		local n = `n' + 1
		quietly {
			use `holds', clear
			set obs `n'
			replace model = "ld_${cause}_${model_name}_`t'" in `n'
			tempfile holds
			save `holds', replace
		}
	}

	use `holds', clear
	quietly levelsof model, clean local(holds) separate(,)

	// submit 12
	qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -pe multi_slot 4 -l mem_free=8G -N com_${cause}_${model_name} -hold_jid `holds' "${temp_dir}/run_files/12_put_predictions_together.sh" >> "${base_dir}/${model_name}_job_ids.txt"
	qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"

	foreach sr of global super_regions {

		// submit 13
		qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -pe multi_slot 4 -l mem_free=4G -N tot_${cause}_${model_name}_`sr' -hold_jid com_${cause}_${model_name} "${temp_dir}/run_files/13_calculate_totals_`sr'.sh" >> "${base_dir}/${model_name}_job_ids.txt"
		qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"
		if missing("`totals'") local totals tot_${cause}_${model_name}_`sr'
		else local totals `totals',tot_${cause}_${model_name}_`sr'
	}

	sleep 1000

	// submit 14
	qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -pe multi_slot 10 -l mem_free=20G -N out_${cause}_${model_name} -hold_jid `totals' "${temp_dir}/run_files/14_output_totals.sh" >> "${base_dir}/${model_name}_job_ids.txt"
	qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"

	foreach t of global tests {

		// submit 15
		qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -pe multi_slot 4 -l mem_free=3G -N cov_${cause}_${model_name}_`t' -hold_jid out_${cause}_${model_name} "${temp_dir}/run_files/15_coverage_`t'.sh" >> "${base_dir}/${model_name}_job_ids.txt"
		qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"
		if missing("`coverage_holds'") local coverage_holds cov_${cause}_${model_name}_`t'
		else local coverage_holds `coverage_holds',cov_${cause}_${model_name}_`t'
	}

	// submit 16
	qui !/usr/local/bin/SGE/bin/lx24-amd64/qsub `sge' -pe multi_slot 4 -l mem_free=1G -N pvf_${cause}_${model_name} -hold_jid `coverage_holds' "${temp_dir}/run_files/16_final_pv.sh" >> "${base_dir}/${model_name}_job_ids.txt"
	qui !echo "\n" >> "${base_dir}/${model_name}_job_ids.txt"

	log close
