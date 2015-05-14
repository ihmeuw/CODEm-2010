/*
Created:	1 May 2011
Updated:	11 April 2012
Purpose: 	create total numbers of deaths for each country, region, and super region
*/

// setup the program
	timer on 1
	set mem 5g
	set more off
	set maxvar 32000
	log using "${base_dir}/logs/13_calculate_totals_${super_region}.smcl", replace
	display "Calculating totals for super region ${super_region} on node `r(o1)'"

// load in the draws
	use "$temp_dir/12_all_draws_${super_region}.dta", clear

// convert the submodel results into the correct space
	forvalues i = 1 / $number_submodels {
		local dv: word `i' of $dv_list
		quietly {
			if "`dv'" == "ln_rate" {
				quietly generate submodel_`i'_spacetime = exp(gpr_`i'_spacetime_mean) * pop
				quietly generate submodel_`i'_linear = exp(linear_`i') * pop
			}
			else if "`dv'" == "lt_cf" {
				quietly generate submodel_`i'_spacetime = invlogit(gpr_`i'_spacetime_mean) * envelope
				quietly generate submodel_`i'_linear = invlogit(linear_`i') * envelope
			}
			drop linear_`i' gpr_`i'_spacetime_mean
		}
	}

// find the median of the distribution by country/year/age
	quietly egen ensemble_median = rowmedian(ensemble_d*)
	quietly egen ensemble_mean = rowmean(ensemble_d*)

// find confidence intervals
	quietly egen ensemble_upper = rowpctile(ensemble_d*), p(97.5)
	quietly egen ensemble_lower = rowpctile(ensemble_d*), p(2.5)

// save the country totals
	quietly generate aggregation_level = "country"
	tempfile country_totals
	save `country_totals', replace

// keep only "reporting" countries
	merge m:1 iso3 using "$temp_dir/country_descriptors.dta", nogen keepusing(reporting reporting_region) keep(match)
	rename reporting_region region
	drop if reporting == 0

// generate regional totals
	drop ensemble_lower ensemble_upper
	collapse (sum) submodel_* ensemble_* envelope pop, by(region year age) fast
	quietly egen ensemble_upper = rowpctile(ensemble_d*), p(97.5)
	quietly egen ensemble_lower = rowpctile(ensemble_d*), p(2.5)
	quietly generate aggregation_level = "region"
	tempfile region_totals
	save `region_totals', replace

// generate super region totals
	merge m:1 region using "$temp_dir/region_descriptors.dta", nogen keepusing(reporting_super_region) keep(match)
	rename reporting_super_region super_region
	drop ensemble_lower ensemble_upper
	collapse (sum) submodel_* ensemble_* envelope pop, by(super_region year age) fast
	quietly egen ensemble_upper = rowpctile(ensemble_d*), p(97.5)
	quietly egen ensemble_lower = rowpctile(ensemble_d*), p(2.5)
	generate aggregation_level = "super_region"

// combine all the predictions back together
	append using `country_totals'

	append using `region_totals'

	tempfile all_age_specific
    save `all_age_specific'

// create totals across all age groups
	drop ensemble_lower ensemble_upper
	collapse (sum) submodel_* ensemble_* envelope pop, by(aggregation_level super_region region iso3 year) fast
	quietly generate age = 99
	quietly egen ensemble_upper = rowpctile(ensemble_d*), p(97.5)
	quietly egen ensemble_lower = rowpctile(ensemble_d*), p(2.5)

// put all the pieces back together
	append using `all_age_specific'
	compress

// save results at the country level
	preserve
	keep if aggregation_level == "country"
	drop aggregation_level super_region region
	save "$temp_dir/13_country_totals_${super_region}.dta", replace

// save results at the regional level
	restore, preserve
	keep if aggregation_level == "region"
	drop aggregation_level super_region iso3
	save "$temp_dir/13_region_totals_${super_region}.dta", replace

// save results at the super region level
	restore
	keep if aggregation_level == "super_region"
	drop aggregation_level region iso3
	save "$temp_dir/13_super_region_totals_${super_region}.dta", replace

// close the logs
	timer off 1
	timer list 1
	log close
