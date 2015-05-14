/*
Created:	1 May 2011
Updated:	11 April 2012
Purpose: 	this program basically just moves stuff around so that each country is now grouped with the super region it reports in, instead of where it's modeled
*/

// setup the program
	timer on 1
	set mem 5g
	set more off
	set odbcmgr unixodbc
	set maxvar 32000
	log using "${base_dir}/logs/12_put_predictions_together.smcl", replace
	display "Putting together the in-sample predictions on node `r(o1)'"

// find the list of which countries are developed/developing, which countries to report, where reporting vs modeling region differs, and get the full region/super region names
	odbc load, exec("SELECT DISTINCT subreg.local_id AS iso3, subreg.name AS country_name, subreg.developed, reg.location_id AS reporting_region, sr.location_id AS reporting_super_region, COALESCE(sr_analytic.location_id,sr.location_id) AS super_region, COALESCE(reg_analytic.location_id,reg.location_id) AS region, COALESCE(sr_analytic.name,sr.name) AS super_region_name, COALESCE(reg_analytic.name,reg.name) AS region_name FROM locations subreg LEFT JOIN locations_indicators indic ON subreg.location_id = indic.location_id LEFT JOIN locations_hierarchy subreg_to_reg ON subreg.location_id = subreg_to_reg.descendant AND subreg_to_reg.type = 'gbd' AND subreg_to_reg.version_id = 2 LEFT JOIN locations reg ON subreg_to_reg.ancestor = reg.location_id LEFT JOIN locations_metadata md_reg ON subreg.location_id = md_reg.location_id AND md_reg.key_id = 2 LEFT JOIN locations reg_analytic ON md_reg.key_value = reg_analytic.location_id LEFT JOIN locations_hierarchy reg_to_sr ON subreg.location_id = reg_to_sr.descendant AND reg_to_sr.type = 'gbd' AND reg_to_sr.version_id = 2 LEFT JOIN locations sr ON reg_to_sr.ancestor = sr.location_id  LEFT JOIN locations_metadata md_sr ON subreg.location_id = md_sr.location_id AND md_sr.key_id = 3 LEFT JOIN locations sr_analytic ON md_sr.key_value = sr_analytic.location_id WHERE subreg.type in ('admin1','admin0','urbanicity','nonsovereign') AND reg.type = 'region' AND sr.type = 'superregion' AND indic.indic_cod = 1") dsn(codc) clear
	generate reporting = 1
    generate country_name = iso3
    generate region_name = region
    generate super_region_name = super_region

    // These get downloaded as strings, don't know why
    destring region, replace force
    destring super_region, replace force

    save "$temp_dir/country_descriptors.dta", replace
	drop if reporting == 0
	drop iso3 country_name reporting
	duplicates drop region, force
	save "$temp_dir/region_descriptors.dta", replace
	drop region region_name reporting_region
	duplicates drop super_region, force
	save "$temp_dir/super_region_descriptors.dta", replace

// load in the draws
	local num_sr: word count $super_regions
	forvalues i = 1/`num_sr' {
		local this_sr: word `i' of $super_regions
		local num_chunks: word `i' of $chunks_per_sr
		clear
		forvalues this_chunk = 1/`num_chunks' {
			quietly append using "$temp_dir/10_spacetime_draws_`this_sr'_`this_chunk'_insample.dta"
		}
		quietly merge 1:1 iso3 year age using "$temp_dir/11_linear_draws_`this_sr'_insample.dta", keep(match master) nogen
		quietly compress
		save "$temp_dir/tmp_`i'.dta", replace
	}
	clear
	forvalues i = 1/`num_sr' {
		quietly append using "$temp_dir/tmp_`i'.dta"
		rm "$temp_dir/tmp_`i'.dta"
	}
	keep iso3 year age ensemble_* top_submodel_* gpr_*_spacetime_mean linear_* envelope pop

// add on the reporting super region for each country
	merge m:1 iso3 using "$temp_dir/country_descriptors.dta", nogen keepusing(reporting_super_region)
	rename reporting_super_region super_region

// save results for each reporting super region
	levelsof super_region, l(srs) c
	preserve
	foreach sr of local srs {
		keep if super_region == `sr'
		save "$temp_dir/12_all_draws_`sr'.dta", replace
		restore, preserve
	}
	restore, not

** save the death draws for the top model each of the model families as specified by the user
	** top spacetime cause fraction model
	if $top_st_cf == 1 {
		forvalues i = 1/`num_sr' {
			local this_sr: word `i' of $super_regions
			local num_chunks: word `i' of $chunks_per_sr
			clear
			forvalues this_chunk = 1/`num_chunks' {
				quietly append using "$temp_dir/10_spacetime_draws_`this_sr'_`this_chunk'_top_cf_model_insample.dta"
			}
			quietly compress
			tempfile tmp_`i'
			save `tmp_`i'', replace
		}
		clear
		forvalues i = 1/`num_sr' {
			quietly append using `tmp_`i''
		}
		keep iso3 year age ensemble_* envelope pop

	// add on the reporting super region for each country
		merge m:1 iso3 using "$temp_dir/country_descriptors.dta", nogen keepusing(reporting_super_region)
		rename reporting_super_region super_region

	// save results for each reporting super region
		levelsof super_region, l(srs) c
		preserve
		foreach sr of local srs {
			keep if super_region == `sr'
			save "$temp_dir/12_all_draws_`sr'_top_st_cf_model.dta", replace
			restore, preserve
		}
		restore, not
	}
	** top spacetime rate model
	if $top_st_rate == 1 {
		forvalues i = 1/`num_sr' {
			local this_sr: word `i' of $super_regions
			local num_chunks: word `i' of $chunks_per_sr
			clear
			forvalues this_chunk = 1/`num_chunks' {
				quietly append using "$temp_dir/10_spacetime_draws_`this_sr'_`this_chunk'_top_rate_model_insample.dta"
			}
			quietly compress
			tempfile tmp_`i'
			save `tmp_`i'', replace
		}
		clear
		forvalues i = 1/`num_sr' {
			quietly append using `tmp_`i''
		}
		keep iso3 year age ensemble_* envelope pop

	// add on the reporting super region for each country
		merge m:1 iso3 using "$temp_dir/country_descriptors.dta", nogen keepusing(reporting_super_region)
		rename reporting_super_region super_region

	// save results for each reporting super region
		levelsof super_region, l(srs) c
		preserve
		foreach sr of local srs {
			keep if super_region == `sr'
			save "$temp_dir/12_all_draws_`sr'_top_st_rate_model.dta", replace
			restore, preserve
		}
		restore, not
	}
	** top linear cause fraction model
	if $top_lin_cf == 1 {
		forvalues i = 1/`num_sr' {
			local this_sr: word `i' of $super_regions
			use "$temp_dir/11_top_lin_cf_model_draws_`this_sr'_insample.dta", clear
			quietly compress
			tempfile tmp_`i'
			save `tmp_`i'', replace
		}
		clear
		forvalues i = 1/`num_sr' {
			quietly append using `tmp_`i''
		}
		keep iso3 year age ensemble_* envelope pop

	// add on the reporting super region for each country
		merge m:1 iso3 using "$temp_dir/country_descriptors.dta", nogen keepusing(reporting_super_region)
		rename reporting_super_region super_region

	// save results for each reporting super region
		levelsof super_region, l(srs) c
		preserve
		foreach sr of local srs {
			keep if super_region == `sr'
			save "$temp_dir/12_all_draws_`sr'_top_lin_cf_model.dta", replace
			restore, preserve
		}
		restore, not
	}
	** top linear rate model
	if $top_lin_rate == 1 {
		forvalues i = 1/`num_sr' {
			local this_sr: word `i' of $super_regions
			use "$temp_dir/11_top_lin_rate_model_draws_`this_sr'_insample.dta", clear
			quietly compress
			tempfile tmp_`i'
			save `tmp_`i'', replace
		}
		clear
		forvalues i = 1/`num_sr' {
			quietly append using `tmp_`i''
		}
		keep iso3 year age ensemble_* envelope pop

	// add on the reporting super region for each country
		merge m:1 iso3 using "$temp_dir/country_descriptors.dta", nogen keepusing(reporting_super_region)
		rename reporting_super_region super_region

	// save results for each reporting super region
		levelsof super_region, l(srs) c
		preserve
		foreach sr of local srs {
			keep if super_region == `sr'
			save "$temp_dir/12_all_draws_`sr'_top_lin_rate_model.dta", replace
			restore, preserve
		}
		restore, not
	}

// close the logs
	timer off 1
	timer list 1
	log close
