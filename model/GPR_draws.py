'''
Created:	1 May 2011
Updated:	18 August 2011
Purpose:    Run GPR to get draws for the final ensemble
'''

# import the necessary modules
from pymc import gp
import numpy as np
from pylab import rec2csv, csv2rec
from numpy.lib import recfunctions

np.random.seed(5)

# define program
def fit_GPR(infile, outfile, dv_list, scale, number_submodels, test, spacetime_iters, top_submodel):
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
    total_iters = np.sum(spacetime_iters)
    draws = [np.empty(len(country_age_list), 'float') for i in range(total_iters)]
    if (top_submodel > 0):
        top_submodel_draws = [np.empty(len(country_age_list), 'float') for i in range(100)]
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

        # keep track of how many iterations have been added for this model
        iter_counter = 0

        # loop through each submodel
        for m in range(number_submodels):

            # identify the dependent variable for this model
            dv = dv_list[m]

            # continue making predictions if we actually need draws for this model
            if (spacetime_iters[m] > 0) or (m+1 == top_submodel):

                # skip models with no spacetime results
                if all_data['spacetime_' + str(m+1)].dtype != 'float64':
                    for i in range(spacetime_iters[m]):
                        draws[iter_counter][country_age_list==ca] = np.NaN
                        iter_counter += 1
                    if (m+1 == top_submodel):
                        for i in range(100):
                            top_submodel_draws[i][country_age_list==ca] = np.NaN
                    continue

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

                # draw realizations from the data
                realizations = [gp.Realization(M, C) for i in range(spacetime_iters[m])]

                # save the data for this country/age into the results array
                iso3[country_age_list==ca] = ca[0:3]
                age_group[country_age_list==ca] = ca[4:]
                year[country_age_list==ca] = year_list.T
                for i in range(spacetime_iters[m]):
                    try:
                        draws[iter_counter][country_age_list==ca] = realizations[i](year_list)
                    except:
                        print('Failure in ' + ca)
                    iter_counter += 1

                # if it's the top submodel, do 100 additional draws
                if (m+1 == top_submodel):
                    realizations = [gp.Realization(M, C) for i in range(100)]
                    for i in range(100):
                        try:
                            top_submodel_draws[i][country_age_list==ca] = realizations[i](year_list)
                        except:
                            print('Failure in ' + ca)

    # save the results
    print('Saving GPR results')
    names = ['iso3','age_group','year']
    results = np.core.records.fromarrays([iso3,age_group,year], names=names)
    for i in range(total_iters):
        results = recfunctions.append_fields(results, 'ensemble_d' + str(i+1), draws[i])
    if (top_submodel > 0):
        for i in range(100):
            results = recfunctions.append_fields(results, 'top_submodel_d' + str(i+1), top_submodel_draws[i])
    rec2csv(results, outfile)
