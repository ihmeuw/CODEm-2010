/*
Created:	1 May 2011
Updated:	18 August 2011
Purpose:	These are functions which will be used by covariate_selection_program
*/

// Function: 	test_covariates
// Purpose:	test all the covariate combinations given to see if they give the appropriate signs on the coefficients
// Inputs:	level = 1, 2, or 3 - which list of covariates to use (in our case based on priors of how proximal/distal the covariate is)
//			ncov = number of covariates to allow in each model
//			definitelyin = a list of covariates that will definitely be included in the regression (i.e. if we're in level 2, this would include all the level 1 covariates)
//			parent = identifier for the model in the level above this

	// setup the input syntax for the program
		capture program drop test_covariates
		program define test_covariates
		syntax , level(integer) ncov(integer) parent(string) [definitelyin(varlist)]

	// find the subset of models that need to be tested (based on number of covariates and level)
		mata test_these_guys = select(level_`level'_combos_`parent', rowsum(level_`level'_combos_`parent'):==`ncov')
		mata st_numscalar("num_to_run", rows(test_these_guys))
		local num_to_run = num_to_run
		if `num_to_run' > 0 display in red "Running regressions for all `ncov' covariate combinations of level `level'. Total possible combinations: `num_to_run'."

	// loop through each model in the subset
		forvalues i = 1/`num_to_run' {
			local covars
			mata current_model = test_these_guys[`i',.]

	// find the names of the covariates in this model and add them to a list
			forvalues j = 1/${level_`level'_ncov} {
				mata st_local("include_this_covariate",current_model[`j']==1 ? "yes" : "no")
				if "`include_this_covariate'" == "yes" {
					mata st_local("this_covariate_name", level_`level'_covariates[`j'])
					local covars `covars' `this_covariate_name'
				}
			}

	// run the xtmixed regression with the selected covariates
			quietly xtmixed $dv `definitelyin' `covars' $Ia || super_region: || region: || age:, iterate(2) emiterate(2)
			global number_tested = $number_tested + 1

	// loop through to check if each covariate gets the correct betas and is significant in the model
			local good = 1
			foreach c of local covars {
				if regexm("`c'","\*") | "${d`c'}"=="" continue
				if (${d`c'} / _b[`c'] < 0) | (2*(1-normal(abs(_b[`c'] / _se[`c']))) > .05) local good = 0
			}
			foreach c of local definitelyin {
				if regexm("`c'","\*") | "${d`c'}"=="" continue
				if (${d`c'} / _b[`c'] < 0) | (2*(1-normal(abs(_b[`c'] / _se[`c']))) > .05) local good = 0
			}
			di "`covars': " _c

	// if the covariates fulfill the priors, leave them
			if `good' == 1 di in green "good"

	// if they're bad, strip them out from the master models list for this level (leave in for level 1, though)
			else {
				di in red "bad"
				if `level' > 1 {
					mata tmp = J(0,cols(level_`level'_combos_`parent'),.)
					mata st_numscalar("num_remaining_rows", rows(level_`level'_combos_`parent'))
					local num_remaining_rows = num_remaining_rows
					forvalues x = 1/`num_remaining_rows' {
						// note: you'll see here that if a child model is bad, we remove all the parents as well (i.e. if x gives the wrong sign, we remove x, x y, x z, and x y z)
						mata tmp = (level_`level'_combos_`parent'[`x',mm_which(current_model)]!=J(1,`ncov',1) ? tmp \ level_`level'_combos_`parent'[`x',.] : tmp)
					}
					mata level_`level'_combos_`parent' = tmp
				}
				else {
					// for level 1, only strip out the offending model, leave in the children
					mata tmp = J(0,cols(level_`level'_combos_`parent'),.)
					mata st_numscalar("num_remaining_rows", rows(level_`level'_combos_`parent'))
					local num_remaining_rows = num_remaining_rows
					forvalues x = 1/`num_remaining_rows' {
						mata tmp = (level_`level'_combos_`parent'[`x',.]!=current_model ? tmp \ level_`level'_combos_`parent'[`x',.] : tmp)
					}
					mata level_`level'_combos_`parent' = tmp
				}
			}
		}
		end


// Function: 	remove_dupes
// Purpose:	remove any models that are simply subsets of other valid models at a particular level
// 			we'll call a model with covariates x y z the "parent" and models x, y, z, x y, x z, and y z would be its children
//			if parent x y z fits the priors, then all its children should be removed from the model list
// Inputs:	level = 1, 2, or 3 - which list of covariates to use (in our case based on priors of how proximal/distal the covariate is)
//			parent = identifier for the model in the level above this

	// setup the input syntax for the program
		capture program drop remove_dupes
		program define remove_dupes
		syntax , level(integer) parent(string)

	// loop through each of the possible numbers of covariates for models at this level
		forvalues i = 1/${level_`level'_ncov} {

	// find the possible "child models" - those models with the current number of covariates
			mata possible_children = select(level_`level'_combos_`parent', rowsum(level_`level'_combos_`parent'):==`i')
			mata st_numscalar("num_possible_children", rows(possible_children))
			local num_possible_children = num_possible_children

	// find the possible "parent models" - all those models that have more covariates than the current number
			mata possible_parents = select(level_`level'_combos_`parent', rowsum(level_`level'_combos_`parent'):>`i')
			mata st_numscalar("num_possible_parents", rows(possible_parents))
			local num_possible_parents = num_possible_parents

	// loop through every parent/child pair
			forvalues j = 1/`num_possible_children' {
				forvalues k = 1/`num_possible_parents' {

	// figure out which covariates are included in the child model
					mata child_pieces = mm_which(possible_children[`j',.])

	// check whether or not the possible parent contains the same covariates as the child
					mata st_local("i_found_mommy", possible_parents[`k',child_pieces]==J(1,`i',1) ? "yes" : "no")

	// if the possible parent contains all the same covariates, then it is indeed a parent
	// in that case, delete the child from the list of models
					if "`i_found_mommy'" == "yes" {
						mata tmp = J(0,cols(level_`level'_combos_`parent'),.)
						mata st_numscalar("num_remaining_rows", rows(level_`level'_combos_`parent'))
						local num_remaining_rows = num_remaining_rows
						forvalues x = 1/`num_remaining_rows' {
							mata tmp = (level_`level'_combos_`parent'[`x',.]!=possible_children[`j',.] ? tmp \ level_`level'_combos_`parent'[`x',.] : tmp)
						}
						mata level_`level'_combos_`parent' = tmp
					}
				}
			}
		}
		end
