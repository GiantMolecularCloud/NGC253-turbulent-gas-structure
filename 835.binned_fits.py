scriptdir##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '800.info.py'))
execfile(join(scriptdir, '830.MCMC.py'))
data        = fnunpickle(join(datadir, 'data.pickle'))
dendrograms = load_dendrograms()


###################################################################################################
# calculate fits to binned size-linewidth relation
###################################################################################################

# all
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        for R_measure in R_measures:
            for lw_measure in lw_measures:

                x     = fits[co][gal]['size-linewidth'][R_measure][lw_measure]['bins']['window center'].data
                y     = fits[co][gal]['size-linewidth'][R_measure][lw_measure]['bins']['y median'].data
                y_p16 = fits[co][gal]['size-linewidth'][R_measure][lw_measure]['bins']['y 16th'].data
                y_p84 = fits[co][gal]['size-linewidth'][R_measure][lw_measure]['bins']['y 84th'].data
                y_err = np.nanmean([y-y_p16, y_p84-y], axis=0)

                filter = (x>fits[co][gal]['size-linewidth']['fit range'][0]) & (x<fits[co][gal]['size-linewidth']['fit range'][1])

                x     = x[filter]
                y     = y[filter]
                y_err = y_err[filter]

                fit = fit_MCMC(log_x     = np.log10(x),
                               log_y     = np.log10(y),
                               log_x_err = [0.05 for x in x],            # half bin width
                               log_y_err = np.log10(y_err),
                               source=gal, line=co,
                               xlabel = R_measure.replace('_',' '),
                               ylabel = lw_measure.replace('%',r'\%'),
                               savepath = join(plotdir, '09.fits', 'binned', co+'.'+gal, R_measure.replace(' ','_')+'.vs.'+lw_measure.replace(' ','_'))
                              )

                fits[co][gal]['size-linewidth'][R_measure][lw_measure]['binned fit'] = fit


# leaves
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        for R_measure in R_measures:
            for lw_measure in lw_measures:

                x     = fits[co][gal]['size-linewidth (leaves)'][R_measure][lw_measure]['bins']['window center'].data
                y     = fits[co][gal]['size-linewidth (leaves)'][R_measure][lw_measure]['bins']['y median'].data
                y_p16 = fits[co][gal]['size-linewidth (leaves)'][R_measure][lw_measure]['bins']['y 16th'].data
                y_p84 = fits[co][gal]['size-linewidth (leaves)'][R_measure][lw_measure]['bins']['y 84th'].data
                y_err = np.nanmean([y-y_p16, y_p84-y], axis=0)

                filter = (x>fits[co][gal]['size-linewidth (leaves)']['fit range'][0]) & (x<fits[co][gal]['size-linewidth']['fit range'][1])

                x     = x[filter]
                y     = y[filter]
                y_err = y_err[filter]

                fit = fit_MCMC(log_x     = np.log10(x),
                               log_y     = np.log10(y),
                               log_x_err = [0.05 for x in x],            # half bin width
                               log_y_err = np.log10(y_err),
                               source=gal, line=co,
                               xlabel = R_measure.replace('_',' '),
                               ylabel = lw_measure.replace('%',r'\%'),
                               savepath = join(plotdir, '09.fits', 'binned_leaves', co+'.'+gal, R_measure.replace(' ','_')+'.vs.'+lw_measure.replace(' ','_'))
                              )

                fits[co][gal]['size-linewidth (leaves)'][R_measure][lw_measure]['binned fit'] = fit


###################################################################################################
# calculate fits to binned mass-linewidth relation
###################################################################################################

# all
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        x     = fits[co][gal]['mass-linewidth']['bins']['window center'].data
        y     = fits[co][gal]['mass-linewidth']['bins']['y median'].data
        y_p16 = fits[co][gal]['mass-linewidth']['bins']['y 16th'].data
        y_p84 = fits[co][gal]['mass-linewidth']['bins']['y 84th'].data
        y_err = np.nanmean([y-y_p16, y_p84-y], axis=0)

        filter = (x>fits[co][gal]['mass-linewidth']['fit range'][0]) & (x<fits[co][gal]['mass-linewidth']['fit range'][1])

        x     = x[filter]
        y     = y[filter]
        y_err = y_err[filter]

        fit = fit_MCMC(log_x     = np.log10(x),
                       log_y     = np.log10(y),
                       log_x_err = [0.05 for x in x],            # half bin width
                       log_y_err = np.log10(y_err),
                       source=gal, line=co,
                       xlabel = 'mass',
                       ylabel = 'linewidth (astrodendro)',
                       savepath = join(plotdir, '09.fits', 'binned', co+'.'+gal, 'mass.vs.linewidth_astrodendro')
                      )

        fits[co][gal]['mass-linewidth']['binned fit'] = fit


###################################################################################################
# calculate fits to binned size-mass relation
###################################################################################################

# all
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        x     = fits[co][gal]['size-mass']['bins']['window center'].data
        y     = fits[co][gal]['size-mass']['bins']['y median'].data
        y_p16 = fits[co][gal]['size-mass']['bins']['y 16th'].data
        y_p84 = fits[co][gal]['size-mass']['bins']['y 84th'].data
        y_err = np.nanmean([y-y_p16, y_p84-y], axis=0)

        filter = (x>fits[co][gal]['size-mass']['fit range'][0]) & (x<fits[co][gal]['size-mass']['fit range'][1])

        x     = x[filter]
        y     = y[filter]
        y_err = y_err[filter]

        fit = fit_MCMC(log_x     = np.log10(x),
                       log_y     = np.log10(y),
                       log_x_err = [0.05 for x in x],            # half bin width
                       log_y_err = np.log10(y_err),
                       source=gal, line=co,
                       xlabel = 'size (astrodendro)',
                       ylabel = 'mass',
                       savepath = join(plotdir, '09.fits', 'binned', co+'.'+gal, R_measure.replace(' ','_')+'.vs.'+lw_measure.replace(' ','_'))
                      )

        fits[co][gal]['size-mass']['binned fit'] = fit


###################################################################################################
# save fits
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(fits[co][gal], join(compdir,gal+'.'+co+'.fits.pickle'))


###################################################################################################
#
###################################################################################################
