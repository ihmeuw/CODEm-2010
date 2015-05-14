/*
Created:	1 May 2011
Updated:	18 August 2011
Purpose: 	Run the linear portion of the model with region random effects (for spacetime prior)
*/

// setup the program
    log using "${base_dir}/logs/1_linear_region_${submodel_chunk}_${test}.smcl", replace
	display "Running linear submodels ${start_submodel} through ${end_submodel} with regional random effects for test ${test} on node `r(o1)'"
	timer on 1
	set more off
	set mem 2g
	set matsize 11000
	set maxvar 32000



// load in data
	use "${temp_dir}/input_database.dta", clear

// make temporary variables that include the DV only for training data (this makes it faster so that if statements aren't necessary all the time below)
	generate train_ln_rate = ln_rate if test_${test} == 0
	generate train_lt_cf = lt_cf if test_${test} == 0

// make age dummies if there's more than one age group
	quietly tab age
	if `r(r)' > 1 {
		xi i.age, prefix(_I)
		local age_dummy _Ia*
	}
	else local age_dummy

// Add on reference covariates
    if $ref_cov == 1 {

        // figure out if we need any log or logit replacements
        qui ds
        local varlist "`r(varlist)'"
        local ref_transforms
        local pred_transforms

        foreach c of local varlist {


            // If there exists a variable that matches something on the prediction list with a transform modifier,
            //   record it to transform reference variabl later
            local count_word 1
            foreach pred_c of global pred_cov_list {

                if substr("`c'",1,3) == "ln_" & substr("`c'",-3,.) == "_sq" & substr("`c'",4,(length("`c'")-6)) == "`pred_c'" {
                    local this_ref_name = word("$ref_cov_list",`count_word')
                    local pred_transforms `pred_transforms' `c'
                    local ref_transforms "`ref_transforms' ln_`this_ref_name'_sq"
                }
                else if substr("`c'",1,3) == "lt_" & substr("`c'",-3,.) == "_sq" & substr("`c'",4,(length("`c'")-6)) == "`pred_c'" {
                    local this_ref_name = word("$ref_cov_list",`count_word')
                    local pred_transforms `pred_transforms' `c'
                    local ref_transforms "`ref_transforms' lt_`this_ref_name'_sq"

                }
                else if substr("`c'",1,3) == "ln_" & substr("`c'",-3,.) != "_sq" & substr("`c'",4,(length("`c'")-3)) == "`pred_c'" {
                    local this_ref_name = word("$ref_cov_list",`count_word')
                    local pred_transforms `pred_transforms' `c'
                    local ref_transforms "`ref_transforms' ln_`this_ref_name'"
                }
                else if substr("`c'",1,3) == "lt_" & substr("`c'",-3,.) != "_sq" & substr("`c'",4,(length("`c'")-3)) == "`pred_c'" {
                    local this_ref_name = word("$ref_cov_list",`count_word')
                    local pred_transforms `pred_transforms' `c'
                    local ref_transforms "`ref_transforms' lt_`this_ref_name'"
                }
                else if substr("`c'",-3,.) == "_sq" & substr("`c'",1,3) != "ln_" & substr("`c'",1,3) != "lt_" & substr("`c'",1,(length("`c'")-3)) == "`pred_c'" {
                    local this_ref_name = word("$ref_cov_list",`count_word')
                    local pred_transforms `pred_transforms' `c'
                    local ref_transforms "`ref_transforms' `this_ref_name'_sq"
                }

                local count_word = `count_word' + 1
            }
        }

        // merge in ref covars
        merge m:1 $merge_vars using "$ref_cov_path", keepusing($ref_cov_list) keep(master match) nogen

        // make ref covar transforms
        if length("`ref_transforms'") > 0 {
            foreach c_t of local ref_transforms {
                if substr("`c_t'",1,3) == "ln_" & substr("`c_t'",-3,.) == "_sq" {
                    local stub = substr("`c_t'",4,length("`c_t'")-6)
                    cap generate `c_t' = (ln(`stub'))^2
                }
                else if substr("`c_t'",1,3) == "lt_" & substr("`c_t'",-3,.) == "_sq" {
                    local stub = substr("`c_t'",4,length("`c_t'")-6)
                    cap generate `c_t' = (logit(`stub'))^2
                }
                else if substr("`c_t'",1,3) == "ln_" & substr("`c_t'",-3,.) != "_sq" {
                    local stub = substr("`c_t'",4,.)
                    cap generate `c_t' = ln(`stub')
                }
                else if substr("`c_t'",1,3) == "lt_" & substr("`c_t'",-3,.) != "_sq" {
                    local stub = substr("`c_t'",4,.)
                    cap generate `c_t' = logit(`stub')
                }
                else if substr("`c_t'",-3,.) == "_sq" & substr("`c_t'",1,3) != "ln_" & substr("`c_t'",1,3) != "lt_" {
                    local stub = substr("`c_t'",1,length("`c_t'")-3)
                    cap generate `c_t' = `stub'^2
                }
            }
        }


        // Get number of ref covars
        local num_refs =  wordcount("$ref_cov_list")
        local num_ref_ts = wordcount("`ref_transforms'")

        // Replace reference covariates
        forvalues ref_n = 1/`num_refs' {
            local this_pred_var = word("$pred_cov_list",`ref_n')
            generate temp_pred_var_`ref_n' = `this_pred_var'
            local this_ref_var = word("$ref_cov_list",`ref_n')
            // Subnational malalria needs this, I don't think it will affect others
            // This replaces the reference covariate with the prediction covariate when the reference covariate is missing. Makes sence, right?
            replace `this_ref_var' = `this_pred_var' if `this_ref_var' == .
            replace `this_pred_var' = `this_ref_var'
        }
        if `num_ref_ts' > 0 {
            forvalues ref_n = 1/`num_ref_ts' {
                local ref_n_tot = `num_refs' + `ref_n'
                local this_pred_var = word("`pred_transforms'",`ref_n')
                generate temp_pred_var_`ref_n_tot' = `this_pred_var'
                local this_ref_var = word("`ref_transforms'",`ref_n')
                // Subnational malaria needs this, I don't think it will affect others
                // This replaces the reference covariate with the prediction covariate when the reference covariate is missing. Makes sence, right?
                replace `this_ref_var' = `this_pred_var' if `this_ref_var' == .
                replace `this_pred_var' = `this_ref_var'
            }
        }

    }

// loop through each submodel
	forvalues i = $start_submodel / $end_submodel {

// run the vanilla xtmixed model
		if "${type_`i'}" == "covariate" {
			display in red _n "Covariate Model `i'" _n "Covariates: ${covariates_`i'}" _n "Dependent Variable: ${dv_`i'}"
			// run the regression
			xtmixed train_${dv_`i'} ${covariates_`i'} `age_dummy' || super_region: || region: || age:, iterate(30)

            // predict the linear and random effects components for the subnational level
			quietly {
				capture predict xb_`i'_sn, xb
				if _rc generate xb_`i'_sn = .
				capture predict u_s_`i'_sn u_r_`i'_sn u_a_`i'_sn, reffects
				if _rc {
					generate u_s_`i'_sn = .
					generate u_r_`i'_sn = .
					generate u_a_`i'_sn = .
				}
				// fill in random effect predictions (because xtmixed only makes in-sample predictions in some cases)
				bysort super_region: egen u_s_`i'_mean_sn = mean(u_s_`i'_sn)
				replace u_s_`i'_sn = u_s_`i'_mean_sn if u_s_`i'_sn == .
				replace u_s_`i'_sn = 0 if u_s_`i'_sn == .
				bysort region: egen u_r_`i'_mean_sn = mean(u_r_`i'_sn)
				replace u_r_`i'_sn = u_r_`i'_mean_sn if u_r_`i'_sn == .
				replace u_r_`i'_sn = 0 if u_r_`i'_sn == .
				bysort region age: egen u_a_`i'_mean_sn = mean(u_a_`i'_sn)
				replace u_a_`i'_sn = u_a_`i'_mean_sn if u_a_`i'_sn == .
				replace u_a_`i'_sn = 0 if u_a_`i'_sn == .
				// generate the final predictions
				generate linear_`i'_${test}_subnatl = xb_`i'_sn + u_s_`i'_sn + u_r_`i'_sn + u_a_`i'_sn
				drop u_*_mean_sn
			}

            // replace prediction vars
            if $ref_cov == 1 {
                forvalues ref_n = 1/`num_refs' {
                    local this_pred_var = word("$pred_cov_list",`ref_n')
                    replace `this_pred_var' = temp_pred_var_`ref_n'
                }
                if `num_ref_ts' > 0 {
                    forvalues ref_n = 1/`num_ref_ts' {
                        local ref_n_tot = `num_refs' + `ref_n'
                        local this_pred_var = word("`pred_transforms'",`ref_n')
                        replace `this_pred_var' = temp_pred_var_`ref_n_tot'
                    }
                }
            }

			// predict the linear and random effects components
			quietly {
				capture predict xb_`i', xb
				if _rc generate xb_`i' = .
				capture predict u_s_`i' u_r_`i' u_a_`i', reffects
				if _rc {
					generate u_s_`i' = .
					generate u_r_`i' = .
					generate u_a_`i' = .
				}
				// fill in random effect predictions (because xtmixed only makes in-sample predictions in some cases)
				bysort super_region: egen u_s_`i'_mean = mean(u_s_`i')
				replace u_s_`i' = u_s_`i'_mean if u_s_`i' == .
				replace u_s_`i' = 0 if u_s_`i' == .
				bysort region: egen u_r_`i'_mean = mean(u_r_`i')
				replace u_r_`i' = u_r_`i'_mean if u_r_`i' == .
				replace u_r_`i' = 0 if u_r_`i' == .
				bysort region age: egen u_a_`i'_mean = mean(u_a_`i')
				replace u_a_`i' = u_a_`i'_mean if u_a_`i' == .
				replace u_a_`i' = 0 if u_a_`i' == .
				// generate the final predictions
				generate linear_`i'_${test} = xb_`i' + u_s_`i' + u_r_`i' + u_a_`i'
				drop u_*_mean
			}

            // replace reference vars
            if $ref_cov == 1 {
                forvalues ref_n = 1/`num_refs' {
                    local this_pred_var = word("$pred_cov_list",`ref_n')
                    replace temp_pred_var_`ref_n' = `this_pred_var'
                    local this_ref_var = word("$ref_cov_list",`ref_n')
                    replace `this_pred_var' = `this_ref_var'
                }
                if `num_ref_ts' > 0 {
                    forvalues ref_n = 1/`num_ref_ts' {
                        local ref_n_tot = `num_refs' + `ref_n'
                        local this_pred_var = word("`pred_transforms'",`ref_n')
                        replace temp_pred_var_`ref_n_tot' = `this_pred_var'
                        local this_ref_var = word("`ref_transforms'",`ref_n')
                        replace `this_pred_var' = `this_ref_var'
                    }
                }
            }

		}

// or use custom regression results
		else if "${type_`i'}" == "custom" {
			display in red "Custom Model" _n "Filepath: ${custom_`i'}" _n "Dependent Variable: ${dv_`i'}"
			// put your code for the custom regression here
			// it should have as its result a variable called linear_`i'_${test} that is in the correct space (ln(rate) or logit(cf))
			// e.g. here is a super stupid model that just predicts 1 for everything
			generate linear_`i'_${test} = 1
		}

// floor linear predictions
		if "${dv_`i'}" == "ln_rate"{
            quietly replace linear_`i'_${test} = clip(linear_`i'_${test}, ln($linear_floor/100000), .)
            quietly replace linear_`i'_${test}_subnatl = clip(linear_`i'_${test}, ln($linear_floor/100000), .)
        }
		else if "${dv_`i'}" == "lt_cf" {
            quietly replace linear_`i'_${test} = clip(linear_`i'_${test}, logit(($linear_floor/100000)/(envelope/pop)), .)
            quietly replace linear_`i'_${test}_subnatl = clip(linear_`i'_${test}, logit(($linear_floor/100000)/(envelope/pop)), .)
        }

// put a ceiling on the rate predictions to never let them exceed the envelopes
		if "${dv_`i'}" == "ln_rate"{
            quietly replace linear_`i'_${test} = clip(linear_`i'_${test}, ., ln(envelope/pop))
            quietly replace linear_`i'_${test}_subnatl = clip(linear_`i'_${test}, ., ln(envelope/pop))
        }

// find the residuals
		quietly generate residual_`i'_${test} = train_${dv_`i'} - linear_`i'_${test}_subnatl
	}

// for any observations which are missing residuals for any submodel, consider them missing for all submodels (otherwise weights might not sum to 1 later)
	quietly egen no_data = rowmiss(residual_*_${test})
	quietly generate has_data_${test} = (no_data == 0)
	drop no_data
	forvalues i = $start_submodel / $end_submodel {
		quietly replace residual_`i'_${test} = . if has_data_${test} == 0
	}

// make the linear predictions missing for any observation which is missing any of these predictions - otherwise, we'll have bugs later because only certain submodels will exist
	quietly egen missing_linear = rowmiss(linear_*_${test})
	quietly generate has_all_linear = (missing_linear == 0)
	forvalues i = $start_submodel / $end_submodel {
		quietly replace linear_`i'_${test} = . if has_all_linear == 0
	}
	drop missing_linear has_all_linear

// save the first stage results for each super region
	preserve
	foreach s of global super_regions {
		keep if super_region == `s'
        if "$good_data_c_years_path" != "" {
            keep iso3 year region super_region national age envelope pop lt_cf lt_cf_sd ln_rate ln_rate_sd n *_${test} good_data
        }
        else {
            keep iso3 year region super_region national age envelope pop lt_cf lt_cf_sd ln_rate ln_rate_sd n *_${test}
        }

		rensfix _${test}
		rename test test_${test}
		rename predictme predictme_${test}
		save "$temp_dir/1_linear_region_${submodel_chunk}_${test}_`s'.dta", replace
		restore, preserve
	}

// close the logs
	timer off 1
	timer list 1
	rm "$temp_dir/1_linear_region_${submodel_chunk}_${test}.sql"
	log close
