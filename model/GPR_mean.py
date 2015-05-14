'''
Created:	1 May 2011
Updated:	18 August 2011
Purpose:    Run GPR on the spacetime predictions - this version only finds the mean, not the CI
'''

# import the necessary modules
from pymc import gp
import numpy as np
from pylab import rec2csv, csv2rec
from numpy.lib import recfunctions

np.random.seed(5)

# define program
def fit_GPR(infile, outfile, dv_list, scale, number_submodels, test):
    # load in the data
    all_data = csv2rec(infile, use_mrecords=False)
    for m in range(number_submodels):
        if all_data['spacetime_' + str(m+1)].dtype == 'float64':
            all_data = np.delete(all_data, np.where(np.isnan(all_data['spacetime_' + str(m+1)]))[0], axis=0)

    # find the list of years for which we need to predict
    year_list = np.unique(all_data.year)

    # find the list of country/age groups
    country_age = np.array([str(all_data.iso3[i]) + '_' + str(all_data.age_group[i]) for i in range(len(all_data))])
    country_age_list = np.repeat(np.unique(country_age), len(year_list))

    # make empty arrays in which to store the results
    draws = [np.empty(len(country_age_list), 'float') for i in range(number_submodels)]
    iso3 = np.empty(len(country_age_list), '|S3')
    age_group = np.empty(len(country_age_list), 'int')
    year = np.empty(len(country_age_list), 'int')

    # loop through country/age groups
    for ca in np.unique(country_age_list):
        print('GPRing ' + ca)

        # subset the data for this particular country/age
        ca_data = all_data[country_age==ca]

        # subset just the observed data
        if ca_data['lt_cf'].dtype != '|O8':
            ca_observed = ca_data[(np.isnan(ca_data['lt_cf'])==0) & (ca_data['test_' + test]==0)]
            if len(ca_observed) > 1:
                has_data = True
            else:
                has_data = False
        else:
            has_data = False

        # loop through each submodel
        for m in range(number_submodels):

            # skip models with no spacetime results
            if all_data['spacetime_' + str(m+1)].dtype != 'float64':
                draws[m][country_age_list==ca] = np.NaN
                continue

            # identify the dependent variable for this model
            dv = dv_list[m]

            # make a list of the spacetime predictions
            ca_prior = np.array([np.mean(ca_data['spacetime_' + str(m+1)][ca_data.year==y]) for y in year_list])

            # find the amplitude for this country/age
            amplitude = np.mean(ca_data['spacetime_amplitude_' + str(m+1)])

            # make a linear interpolation of the spatio-temporal predictions to use as the mean function for GPR
            def mean_function(x) :
                return np.interp(x, year_list, ca_prior)

            # setup the covariance function
            M = gp.Mean(mean_function)
            C = gp.Covariance(eval_fun=gp.matern.euclidean, diff_degree=2, amp=amplitude, scale=scale)

            # observe the data if there is any
            if has_data:
                gp.observe(M=M, C=C, obs_mesh=ca_observed.year, obs_V=ca_observed['spacetime_data_variance_' + str(m+1)], obs_vals=ca_observed[dv])

            # save the data for this country/age into the results array
            iso3[country_age_list==ca] = ca[0:3]
            age_group[country_age_list==ca] = ca[4:]
            year[country_age_list==ca] = year_list.T
            draws[m][country_age_list==ca] = M(year_list)

    # save the results
    print('Saving GPR results')
    names = ['iso3','age_group','year']
    results = np.core.records.fromarrays([iso3,age_group,year], names=names)
    for m in range(number_submodels):
        results = recfunctions.append_fields(results, 'gpr_' + str(m+1) + '_spacetime_mean', draws[m])
    rec2csv(results, outfile)
