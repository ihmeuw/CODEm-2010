/*
Created:	1 May 2011
Updated:	18 August 2011
Purpose:	Select covariates based on priors for directionality and level
*/

// setup stata
	capture log close
	log using "${base_dir}/covariate_selection/covariate_selection_${dv}_log.smcl", replace
	clear
	clear mata
	clear matrix
	set more off
	set mem 1g
	// NOTE: you must have "moremata" installed for this to run; if it fails with a message about "mm_which not found", run -ssc install moremata- then restart stata

// load in the functions necessary to perform the actual covariate tests
// they are in a separate do file so that comments can be included for them more easily
	do "${code_dir}/covariate_selection/covariate_selection_functions.do"

// load in the covariate priors
	insheet using "${base_dir}/covariate_selection/covariate_priors.csv", comma clear
	sort level covariate

// keep a counter of how many models are tested
	global number_tested = 0

// figure out how many level 1, 2, and 3 covariates there are
	local total_covariates = _N
	forvalues i = 1/3 {
		count if level==`i'
		global level_`i'_ncov = `r(N)'
	}
	global initial_pool = (2^${level_1_ncov} - 1) * (2^${level_2_ncov}) * (2^${level_3_ncov})

// move level 2 up if there's no level 1, move level 3 up if there's no level 2
	if $level_1_ncov==0 & $level_2_ncov > 0 {
		replace level = 1 if level == 2
		replace level = 2 if level == 3
	}
	else if $level_1_ncov==0 & $level_2_ncov==0 {
		replace level = 1 if level == 3
	}
	else if $level_2_ncov==0 {
		replace level = 2 if level == 3
	}

// recalculate how many covariates are at each level
	forvalues i = 1/3 {
		count if level==`i'
		global level_`i'_ncov = `r(N)'
	}

// save the covariates into mata by level
	if $level_1_ncov > 0 mata level_1_covariates = st_sdata((1::$level_1_ncov), "covariate")
	if $level_2_ncov > 0 mata level_2_covariates = st_sdata(($level_1_ncov+1::$level_1_ncov+$level_2_ncov), "covariate")
	if $level_3_ncov > 0 mata level_3_covariates = st_sdata(($level_1_ncov+$level_2_ncov+1::$level_1_ncov+$level_2_ncov+$level_3_ncov), "covariate")

// save the desired directions into globals
	forvalues i = 1/`total_covariates' {
		local cov = covariate[`i']
		local dir = direction[`i']
		global d`cov' `dir'
	}

// make a matrix of all the permutations at each level
// this will make matrices such as "level_one_combos" that have indicator variables for every possible set of covariates within the level
// it's basically a binary decomposition (look it up)
// or after running this block type "mata level_1_combos" and stare at it for a while until you get it
	foreach l in 1 2 3 {
		clear
		set obs ${level_`l'_ncov}
		if ${level_`l'_ncov} == 1 set obs 2
		forvalues i = 1/${level_`l'_ncov} {
			generate v`i' = (_n==`i')
		}
		if ${level_`l'_ncov} > 1 {
			fillin *
			drop _fillin
		}
		if ${level_`l'_ncov} > 0 {
			mata level_`l'_combos = st_data(.,.)
			if `l' == 1 mata level_`l'_combos = select(level_`l'_combos, rowsum(level_`l'_combos):>0)
			mata level_`l'_combos = sort(level_`l'_combos, (1::${level_`l'_ncov})')
		}
		else mata level_`l'_combos = J(0,0,0)
	}

// load in the test dataset
	insheet using "${base_dir}/input_database_cv.csv", comma clear case


// turn age into dummies (we do this here instead of using xi to save a tiny amount of time inside the loops....)
	quietly tab age
	if `r(r)' > 1 {
		xi i.age
		global Ia _Iage*
	}
	else global Ia

// copy the level 1 combos into a list with parent 0, as there is no parent model (because it's the top of the hierarchy)
	mata level_1_combos_0 = level_1_combos

// perform covariate selection on level 1
	forvalues i = 1/$level_1_ncov {
		test_covariates, ncov(`i') level(1) parent("0")
	}

// get rid of any level 1 models that are simply subsets of other valid level 1 models
	** remove_dupes, level(1) parent("0")

// if there are no covariates in level 1 for this model, replace it with blanks
	mata st_numscalar("num_found", rows(level_1_combos_0))
	if num_found==0 mata level_1_combos_0 = J(1,cols(level_1_combos_0),0)

// loop through each level 1 model to test all the level 2 covariates on it
	mata st_numscalar("num_level_1_results", rows(level_1_combos_0))
	global num_level_1_results = num_level_1_results
	di "Number of level 1 models found: " num_level_1_results
	forvalues i = 1/$num_level_1_results {

// create a matrix of level 2 models corresponding to this parent
		mata level_2_combos_`i' = level_2_combos

// look up the covariates for the current level one model
		mata current_level_1_model = level_1_combos_0[`i',.]
		local level_1_covars
		forvalues j = 1/$level_1_ncov {
			mata st_local("include_this_covariate",current_level_1_model[`j']==1 ? "yes" : "no")
			if "`include_this_covariate'" == "yes" {
				mata st_local("this_covariate_name", level_1_covariates[`j'])
				local level_1_covars `level_1_covars' `this_covariate_name'
			}
		}
		di in red _n "Testing covariates to add to model `level_1_covars'"

// run level 2 covariate selection for this level 1 model
		forvalues j = 1/$level_2_ncov {
			test_covariates, ncov(`j') level(2) parent("`i'") definitelyin(`level_1_covars')
		}

// remove duplicates from the 2nd level selections
		remove_dupes, level(2) parent("`i'")

// if there are no covariates in level 2 for this model, replace it with blanks
		mata st_numscalar("num_found", rows(level_2_combos_`i'))
		di "Number of level 2 models found for level 1 model `i': " num_found
		if num_found==0 mata level_2_combos_`i' = J(1, cols(level_2_combos_`i')>0 ? cols(level_2_combos_`i') : 1, 0)
	}

// loop through each level 1 model again
	forvalues i = 1/$num_level_1_results {

// look up the covariates for the current level one model
		mata current_level_1_model = level_1_combos_0[`i',.]
		local level_1_covars
		forvalues j = 1/$level_1_ncov {
			mata st_local("include_this_covariate",current_level_1_model[`j']==1 ? "yes" : "no")
			if "`include_this_covariate'" == "yes" {
				mata st_local("this_covariate_name", level_1_covariates[`j'])
				local level_1_covars `level_1_covars' `this_covariate_name'
			}
		}

// loop through each level 2 model corresponding to that level 1 model
		mata st_numscalar("num_level_2_results_`i'", rows(level_2_combos_`i'))
		global num_level_2_results_`i' = num_level_2_results_`i'
		forvalues j = 1/${num_level_2_results_`i'} {

// create a matrix of models for this level one model
			mata level_3_combos_`i'_`j' = level_3_combos

// look up the covariates for the current level two model
			mata current_level_2_model = level_2_combos_`i'[`j',.]
			local level_2_covars
			forvalues k = 1/$level_2_ncov {
				mata st_local("include_this_covariate",current_level_2_model[`k']==1 ? "yes" : "no")
				if "`include_this_covariate'" == "yes" {
					mata st_local("this_covariate_name", level_2_covariates[`k'])
					local level_2_covars `level_2_covars' `this_covariate_name'
				}
			}
			di in red _n "Testing covariates to add to model `level_1_covars' `level_2_covars'"

// run level 3 covariate selection for this level 2 model
			forvalues k = 1/$level_3_ncov {
				test_covariates, ncov(`k') level(3) parent("`i'_`j'") definitelyin(`level_1_covars' `level_2_covars')
			}

// remove duplicates from the 2nd level selections
			remove_dupes, level(3) parent("`i'_`j'")

// if there are no covariates in level 3 for this model, replace it with blanks
			mata st_numscalar("num_found", rows(level_3_combos_`i'_`j'))
			di "Number of level 3 models found for level 2 model `i'_`j': " num_found
			if num_found==0 mata level_3_combos_`i'_`j' = J(1,cols(level_3_combos_`i'_`j'),0)
		}
	}

// make space to save the suggested models into Stata
	clear
	foreach l in 1 2 3 {
		forvalues i = 1/${level_`l'_ncov} {
			generate v`l'_`i' = 0
		}
	}

// loop through all the suggested models
	forvalues i = 1/$num_level_1_results {
		forvalues j = 1/${num_level_2_results_`i'} {
			mata st_numscalar("num_level_3_results_`i'_`j'", rows(level_3_combos_`i'_`j'))
			global num_level_3_results_`i'_`j' = num_level_3_results_`i'_`j'
			forvalues k = 1/${num_level_3_results_`i'_`j'} {

// add an observation for each model
				mata st_addobs(1)

// fill in the covariates for this model
				forvalues x = 1/$level_1_ncov {
					mata st_store(st_nobs(), "v1_`x'", level_1_combos_0[`i',`x'])
				}
				forvalues x = 1/$level_2_ncov {
					mata st_store(st_nobs(), "v2_`x'", level_2_combos_`i'[`j',`x'])
				}
				forvalues x = 1/$level_3_ncov {
					mata st_store(st_nobs(), "v3_`x'", level_3_combos_`i'_`j'[`k',`x'])
				}
			}
		}
	}

// output the list of covariates into a text file
	file open output using	"$base_dir/covariate_selection/selected_covariates_${dv}.do", write replace text
	egen number_covariates = rowtotal(*)
	drop if number_covariates == 0
	count
	local N = `r(N)'
	di in red _n _n "Final Results for $dv" _n _n "Initial Model Pool: $initial_pool" _n "Regressions Run: $number_tested" _n "Models Picked: `N'" _n _n
	if "${dv}" == "lt_cf" global dv_short "cf"
	else if "${dv}" == "ln_rate" global dv_short "rate"
	file write output 		"global num_${dv_short}_models `N'" _n
	forvalues i = 1/`N' {
		di in red _n "Model `i': " _c
		file write output 	"global ${dv_short}_model_`i' "
		forvalues l = 1/3 {
			forvalues x = 1/${level_`l'_ncov} {
				if v`l'_`x'[`i'] == 1 {
					mata st_local("this_covariate_name", level_`l'_covariates[`x'])
					di "`this_covariate_name' " _c
					file write output "`this_covariate_name' "
				}

			}
		}
		file write output 	_n
	}
	file close output
	log close
