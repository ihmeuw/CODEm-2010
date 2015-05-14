/*
Created:	1 May 2011
Updated:	12 June 2012
Purpose: 	Make draws from the linear model
*/

// setup the program
	timer on 1
	set more off
	set mem 5g
	set matsize 11000
	set maxvar 32000
	log using "${base_dir}/logs/11_linear_draws_${test}.smcl", replace
	display "Running linear model for test ${test} on node `r(o1)'"

** safeguard to prevent errors for people using old runfiles
if "$psi_int" == "" global psi_int = 0.01

// find the optimal value of psi
	use "$temp_dir/9_ensemble_pv.dta", clear
	quietly summarize psi if test == "combined out-of-sample" & total_rank_t1 == 1, meanonly
	local optimal_psi = `r(mean)'
	if substr(string(`r(mean)'), 5, 1) != "" {
		di in red "This model has an unusual 3 way tie between nonadjacent psi values. We rounded the average psi to the nearest values of $psi_int."
		local optimal_psi =  round(`r(mean)', $psi_int)
	}
	else local optimal_psi = `r(mean)'
	local optimal_psi_str = subinstr(string(`optimal_psi'), ".", "x", .)

// load in the number of draws to make for each model
	use "$temp_dir/7_submodel_ranks.dta", clear
	forvalues i = 1 / $number_submodels {
		quietly levelsof draws_psi_`optimal_psi_str' if spacetime_or_linear == "linear" & submodel == `i' & test == "combined out-of-sample", l(linear_draws_`i') c
	}
	summarize draws_psi_`optimal_psi_str' if spacetime_or_linear == "spacetime" & test == "combined out-of-sample", meanonly
	local counter = `r(sum)'
	summarize draws_psi_`optimal_psi_str' if spacetime_or_linear == "linear" & test == "combined out-of-sample", meanonly
	local total_linear_draws = `r(sum)'

// figure out if the number 1 submodel is a linear model
	count if spacetime_or_linear == "linear" & test == "combined out-of-sample" & total_rank_t1 == 1
	local draw_top_submodel = (`r(N)' > 0)
	if `draw_top_submodel' == 1 levelsof submodel if spacetime_or_linear == "linear" & test == "combined out-of-sample" & total_rank_t1 == 1, l(top_submodel) c

** determine the top rate and cause fraction linear models in case the user wants draws from either of them
	preserve
	capture {
		keep if spacetime_or_linear == "linear" & test == "combined out-of-sample" & dependent_variable == "cf"
		sort total_rank_t1
		local top_cf_model = submodel in 1
		noisily di in red `top_cf_model'
	}
	if _rc != 0 local top_cf_model
	restore
	capture {
		keep if spacetime_or_linear == "linear" & test == "combined out-of-sample" & dependent_variable == "rate"
		sort total_rank_t1
		local top_rate_model = submodel in 1
		noisily di in red `top_rate_model'
	}
	if _rc != 0 local top_rate_model

// load in data
	use "$temp_dir/input_database.dta", clear
	quietly generate train_ln_rate = ln_rate if test_${test} == 0
	quietly generate train_lt_cf = lt_cf if test_${test} == 0
	xi i.age, prefix(_I)
	drop ln_rate lt_cf test_${test}

// merge on the NSV
	quietly merge 1:1 n using "$temp_dir/4_gpr_parameters_${test}.dta", keep(match master) nogen keepusing(mad_nsv_universal_ln_rate mad_nsv_universal_lt_cf)

// merge on the xb/RE predictions (because for some reason estimates save doesn't store the random effects, only fixed effects)
	forvalues i = 1 / $number_submodel_chunks {
		quietly merge 1:1 n using "$temp_dir/2_model_effects_`i'_${test}.dta", keep(match master) nogen
	}

// loop through each submodel
	forvalues i = 1 / $number_submodels {

// find the residuals
		local dv: word `i' of $dv_list
		quietly generate linear_`i' = xb_`i' + u_s_`i' + u_r_`i' + u_a_`i' + u_c_`i'
		quietly drop xb_`i'
		if "`dv'" == "ln_rate" quietly replace linear_`i' = clip(linear_`i', ln($linear_floor/100000), ln(envelope/pop))
		else if "`dv'" == "lt_cf" quietly replace linear_`i' = clip(linear_`i', logit(($linear_floor/100000)/(envelope/pop)), .)
		quietly generate residual_`i' = train_`dv' - linear_`i'

// find the MAD of the residuals
		quietly {
			bysort iso3 age: egen tmp = count(lt_cf) if source_type == "VR" & national == 1
			by iso3 age: egen num_vr = mean(tmp)
			by iso3 age: egen num_tot = count(lt_cf)
			generate mad_type = "VR" if num_vr == num_tot & num_vr >= 10
			replace mad_type = "other" if mad_type == ""
			bysort mad_type age: egen mad_`i' = mad(residual_`i')
			drop mad_type tmp num_vr num_tot
			generate systematic_sd_`i' = sqrt((1.4826*mad_`i')^2 - (1.4826*mad_nsv_universal_`dv')^2)
			replace systematic_sd_`i' = (1.4826*mad_`i')/2 if systematic_sd_`i' <= 0 | systematic_sd_`i' == .
		}

// load in the regression results
		local chunk = ceil(`i' / 10)
        // local chunk = `i'
        local row_in_chunk = mod(`i', 10)
        // local row_in_chunk = 1
		if `row_in_chunk' == 0 local row_in_chunk = 10
		estimates use "$temp_dir/2_model_estimates_`chunk'_${test}.ster", number(`row_in_chunk')

** create draws from the covariance matrix to get parameter uncertainty
		matrix m = e(b)'
		matrix m = m[1..(rowsof(m)-5),1]
		local covars: rownames m
		local num_covars: word count `covars'
		local betas
		forvalues j = 1/`num_covars' {
			local this_covar: word `j' of `covars'
			if substr("`this_covar'",1,2) == "o." local this_covar skipme`j'
			local betas `betas' b_`this_covar'
		}
		matrix C = e(V)
		matrix C = C[1..(colsof(C)-5), 1..(rowsof(C)-5)]
		drawnorm `betas', means(m) cov(C)

** save files of top submodel for each group
		if $top_lin_cf == 1 & "`top_cf_model'" == "`i'" {
			local cfi = `i'
			local lin_cf_covars `covars'
			di in red `i'
			tempfile tcfm
			save `tcfm', replace
		}
		if $top_lin_rate == 1 & "`top_rate_model'" == "`i'" {
			local ratei = `i'
			local lin_rate_covars `covars'
			di in red `i'
			tempfile trm
			save `trm', replace
		}

** generate draws of the linear
		forvalues j = 1/`linear_draws_`i'' {
			local counter = `counter' + 1
			quietly generate xb_`i'_d`j' = 0
			foreach c of local covars {
				if "`c'" == "_cons" quietly replace xb_`i'_d`j' = xb_`i'_d`j' + b__cons[`j']
				else if substr("`c'",1,2) == "o." continue
				else quietly replace xb_`i'_d`j' = xb_`i'_d`j' + `c' * b_`c'[`j']
			}
			quietly generate ensemble_d`counter' = xb_`i'_d`j' + u_s_`i' + u_r_`i' + u_a_`i' + u_c_`i' + rnormal(0, systematic_sd_`i')
			drop xb_`i'_d`j'
			if "`dv'" == "ln_rate" {
				quietly replace ensemble_d`counter' = exp(ensemble_d`counter') * pop
				quietly replace ensemble_d`counter' = envelope if ensemble_d`counter' > envelope & ensemble_d`counter' != .
			}
			else if "`dv'" == "lt_cf" quietly replace ensemble_d`counter' = invlogit(ensemble_d`counter') * envelope
		}

** create 100 draws if this is the top submodel
		if "`top_submodel'" == "`i'" {
			forvalues j = 1/100 {
				quietly generate top_d`j' = 0
				foreach c of local covars {
					if "`c'" == "_cons" quietly replace top_d`j' = top_d`j' + b__cons[`j']
					else if substr("`c'",1,2) == "o." continue
					else quietly replace top_d`j' = top_d`j' + `c' * b_`c'[`j']
				}
				quietly generate top_submodel_d`j' = top_d`j' + u_s_`i' + u_r_`i' + u_a_`i' + u_c_`i' + rnormal(0, systematic_sd_`i')
				drop top_d`j'
				if "`dv'" == "ln_rate" {
					replace top_submodel_d`j' = exp(top_submodel_d`j') * pop
					quietly replace top_submodel_d`j' = envelope if top_submodel_d`j' > envelope & top_submodel_d`j' != .
				}
				else if "`dv'" == "lt_cf" replace top_submodel_d`j' = invlogit(top_submodel_d`j') * envelope
			}
		}
		drop b_*
		drop u_s_`i' u_r_`i' u_a_`i' u_c_`i'
	}

// save the results
	duplicates drop iso3 age year, force
	if `total_linear_draws' > 0	& `draw_top_submodel' > 0 keep super_region iso3 year age ensemble_d* linear_* top_submodel_d*
	else if `total_linear_draws' > 0 keep super_region iso3 year age ensemble_d* linear_*
	else keep super_region iso3 year age linear_*
	levelsof super_region, l(super_regions) c
	preserve
	foreach sr of local super_regions {
		keep if super_region == `sr'
		drop super_region
		save "$temp_dir/11_linear_draws_`sr'_${test}.dta", replace
		restore, preserve
	}
	restore, not

** make draws of the top submodel for spacetime or linear models if specified by the user
if $top_lin_rate == 1 {
	use `trm', clear
	local i = `ratei'
	di in red `i'
	cap drop ensemble_d*

	forvalues j = 1/1000 {
		quietly generate xb_`i'_d`j' = 0
		foreach c of local lin_rate_covars {
			if "`c'" == "_cons" quietly replace xb_`i'_d`j' = xb_`i'_d`j' + b__cons[`j']
			else if substr("`c'",1,2) == "o." continue
			else quietly replace xb_`i'_d`j' = xb_`i'_d`j' + `c' * b_`c'[`j']
		}
		quietly generate ensemble_d`j' = xb_`i'_d`j' + u_s_`i' + u_r_`i' + u_a_`i' + u_c_`i' + rnormal(0, systematic_sd_`i')
		quietly replace ensemble_d`j' = exp(ensemble_d`j') * pop
		quietly replace ensemble_d`j' = envelope if ensemble_d`j' > envelope & ensemble_d`j' != .
	}
	drop xb_`i'_d*

	duplicates drop iso3 age year, force
	keep super_region iso3 year age ensemble_d*
	preserve
	foreach sr of local super_regions {
		keep if super_region == `sr'
		drop super_region
		save "$temp_dir/11_top_lin_rate_model_draws_`sr'_${test}.dta", replace
		restore, preserve
	}
	restore, not
}
if $top_lin_cf == 1 {
	use `tcfm', clear
	local i = `cfi'
	cap drop ensemble_d*

	forvalues j = 1/1000 {
		quietly generate xb_`i'_d`j' = 0
		foreach c of local lin_cf_covars {
			if "`c'" == "_cons" quietly replace xb_`i'_d`j' = xb_`i'_d`j' + b__cons[`j']
			else if substr("`c'",1,2) == "o." continue
			else quietly replace xb_`i'_d`j' = xb_`i'_d`j' + `c' * b_`c'[`j']
		}
		quietly generate ensemble_d`j' = xb_`i'_d`j' + u_s_`i' + u_r_`i' + u_a_`i' + u_c_`i' + rnormal(0, systematic_sd_`i')
		quietly replace ensemble_d`j' = invlogit(ensemble_d`j') * envelope
	}
	drop xb_`i'_d*

	duplicates drop iso3 age year, force
	keep super_region iso3 year age ensemble_d*
	preserve
	foreach sr of local super_regions {
		keep if super_region == `sr'
		drop super_region
		save "$temp_dir/11_top_lin_cf_model_draws_`sr'_${test}.dta", replace
		restore, preserve
	}
	restore, not
}

// close the logs
	timer off 1
	timer list 1
	log close
