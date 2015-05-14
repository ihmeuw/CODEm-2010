/*
Created:	1 May 2011
Updated:	11 April 2012
Purpose: 	combine all the totals into output files and create global totals
*/

// setup the program
	timer on 1
	set mem 9g
	set more off
	set maxvar 32000
	log using "${base_dir}/logs/14_output_totals.smcl", replace
	display "Outputting totals on node `r(o1)'"

// load in the super region totals
	clear
	foreach sr of global super_regions {
		append using "$temp_dir/13_super_region_totals_`sr'.dta"
	}
	merge m:1 super_region using "$temp_dir/super_region_descriptors.dta", keep(match) nogen keepusing(super_region_name developed )

// save the super region results
	preserve
	keep super_region super_region_name year age  ensemble_median ensemble_mean ensemble_upper ensemble_lower submodel_* envelope pop
	compress
	save "${base_dir}/results/super_region_deaths.dta", replace
	save "$temp_dir/super_region_deaths.dta", replace
	restore

// generate developed/developing results
	preserve
	drop ensemble_upper ensemble_lower super_region super_region_name
	collapse (sum) submodel_* ensemble_* envelope pop, by(year age developed) fast
	quietly egen ensemble_upper = rowpctile(ensemble_d*), p(97.5)
	quietly egen ensemble_lower = rowpctile(ensemble_d*), p(2.5)

// save the developed/developing results
	keep developed year age ensemble_median ensemble_mean ensemble_upper ensemble_lower submodel_* envelope pop
	order developed year age ensemble_mean ensemble_median ensemble_lower ensemble_upper envelope pop submodel_*
	compress
	save "${base_dir}/results/developed_developing_deaths.dta", replace
	restore

// generate global results
	drop ensemble_upper ensemble_lower
	collapse (sum) submodel_* ensemble_* envelope pop, by(year age) fast
	quietly egen ensemble_upper = rowpctile(ensemble_d*), p(97.5)
	quietly egen ensemble_lower = rowpctile(ensemble_d*), p(2.5)

// save the global results
	keep year age ensemble_median ensemble_mean ensemble_upper ensemble_lower submodel_* envelope pop
	order year age ensemble_mean ensemble_median ensemble_lower ensemble_upper envelope pop submodel_*
	compress
	save "${base_dir}/results/global_deaths.dta", replace
	save "$temp_dir/global_deaths.dta", replace
	drop submodel_*
	outsheet using "${base_dir}/results/global_deaths.csv", replace comma

// load in the regional totals
	clear
	foreach sr of global super_regions {
		append using "$temp_dir/13_region_totals_`sr'.dta"
	}
	merge m:1 region using "$temp_dir/region_descriptors.dta", keep(match) nogen keepusing(super_region super_region_name region_name )

// save the region results
	keep super_region super_region_name region region_name year age ensemble_median ensemble_mean ensemble_upper ensemble_lower submodel_* envelope pop
	compress
	save "${base_dir}/results/region_deaths.dta", replace
	save "$temp_dir/region_deaths.dta", replace

// load in the country totals
	clear
	foreach sr of global super_regions {
		append using "$temp_dir/13_country_totals_`sr'.dta"
	}

// save individual draws for the country results
	preserve
	drop if age == 99
	keep iso3 year age ensemble_d* envelope pop
	egen missing = rowmiss(ensemble_d*)
	drop if missing > 980
	drop missing
	generate cause = "$cause"
	generate sex = $sex
	compress
	outsheet using "${base_dir}/results/death_draws.csv", comma replace

	restore

** save the death draws for each of the model families if specified by the user
	if $top_st_cf == 1 {
		preserve
		clear
		foreach sr of global super_regions {
			append using "$temp_dir/12_all_draws_`sr'_top_st_cf_model.dta"
		}
		forvalues i = 1/1000 {
			rename ensemble_d`i' draw_`i'
		}
		keep iso3 year age draw_* envelope pop
		egen missing = rowmiss(draw_*)
		drop if missing > 980
		drop missing
		generate cause = "$cause"
		generate sex = $sex
		compress
		save "${base_dir}/results/death_draws_top_st_cf_model.dta", replace
		restore
	}
	if $top_st_rate == 1 {
		preserve
		clear
		foreach sr of global super_regions {
			append using "$temp_dir/12_all_draws_`sr'_top_st_rate_model.dta"
		}
		forvalues i = 1/1000 {
			rename ensemble_d`i' draw_`i'
		}
		keep iso3 year age draw_* envelope pop
		egen missing = rowmiss(draw_*)
		drop if missing > 980
		drop missing
		generate cause = "$cause"
		generate sex = $sex
		compress
		save "${base_dir}/results/death_draws_top_st_rate_model.dta", replace
		restore
	}
	if $top_lin_cf == 1 {
		preserve
		clear
		foreach sr of global super_regions {
			append using "$temp_dir/12_all_draws_`sr'_top_lin_cf_model.dta"
		}
		forvalues i = 1/1000 {
			rename ensemble_d`i' draw_`i'
		}
		keep iso3 year age draw_* envelope pop
		egen missing = rowmiss(draw_*)
		drop if missing > 980
		drop missing
		generate cause = "$cause"
		generate sex = $sex
		compress
		save "${base_dir}/results/death_draws_top_lin_cf_model.dta", replace
		restore
	}
	if $top_lin_rate == 1 {
		preserve
		clear
		foreach sr of global super_regions {
			append using "$temp_dir/12_all_draws_`sr'_top_lin_rate_model.dta"
		}
		forvalues i = 1/1000 {
			rename ensemble_d`i' draw_`i'
		}
		keep iso3 year age draw_* envelope pop
		egen missing = rowmiss(draw_*)
		drop if missing > 980
		drop missing
		generate cause = "$cause"
		generate sex = $sex
		compress
		save "${base_dir}/results/death_draws_top_lin_rate_model.dta", replace
		restore
	}

// save the country results
	merge m:1 iso3 using "$temp_dir/country_descriptors.dta", keep(match) nogen keepusing(super_region super_region_name region region_name country_name )
	keep super_region super_region_name region region_name iso3 country_name year age ensemble_median ensemble_mean ensemble_upper ensemble_lower submodel_* envelope pop
	order super_region super_region_name region region_name iso3 country_name year age ensemble_mean ensemble_median ensemble_lower ensemble_upper envelope pop submodel_*
	compress
	save "${base_dir}/results/country_deaths.dta", replace
	save "$temp_dir/country_deaths.dta", replace




// close the logs
	timer off 1
	timer list 1
	log close
