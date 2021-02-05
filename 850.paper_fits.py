##############
# DENDROGRAM #
##############

###################################################################################################
# filter dendrograms
###################################################################################################

# Above a certain size, the there is virtually no variance in the properties at a given size scale
# anymore. This happens because at large size scales, large scale effects dominate such as rotation
# and the (turbulent) linewidth/mass/... of the individual structures become irrelevant.
# This effect cannot be removed from the data without a very good 6D model of the gas and gas
# motions. Since the large scale effects completely dominate the cloud scale properties, the
# affected size scales cannot be used for an analysis. They simply need to be flagged out.
# On the smallest scales, the data is incomplete due to limited flux sensitivity and noise, and a
# finite velocity resolution. The smallest scale thus need to be flagged as well in order to not
# distort the relations.

scale_filter = {'CO(3-2)': {'NGC253': [5.5e-1, 18],
                            'GC':     [9.0e-1, 20]},
                'CO(1-0)': {'NGC253': [6.0e0,  60],
                            'GC':     [8.5e0,  68]}}

def filter_catalog(catalog, co, gal):
    Rmin_filter = catalog[R_measure]>scale_filter[co][gal][0]
    Rmax_filter = catalog[R_measure]<scale_filter[co][gal][1]
    R_filter = np.logical_and(Rmin_filter, Rmax_filter)
    type_filter = catalog['type']!='trunk'
    return np.logical_and(R_filter, type_filter)


###################################################################################################
# fit dictionary
###################################################################################################

paperfits = {co: {gal: {ff: {} for ff in ['size-linewidth', 'size-mass', 'size-luminosity', 'mass-linewidth', 'density-virial']}
                       for gal in galaxies.keys()}
                 for co in lines.keys()}


###################################################################################################
# use astrodendro measures
###################################################################################################

R_measure  = 'size (astrodendro)'
lw_measure = 'linewidth (astrodendro)'


###################################################################################################
# calculate binned size-linewidth relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog  = dendrograms[co][gal]['catalog']
        x_filter = filter_catalog(catalog, co, gal)

        bins = get_binned_percentiles(R_measure, lw_measure,
                                      x_step = 0.1,
                                      xdata  = catalog[R_measure][x_filter],
                                      ydata  = catalog[lw_measure][x_filter]
                                     )

        paperfits[co][gal]['size-linewidth']['bins'] = bins


###################################################################################################
# calculate binned size-luminosity relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog  = dendrograms[co][gal]['catalog']
        x_filter = filter_catalog(catalog, co, gal)

        bins = get_binned_percentiles(R_measure, lw_measure,
                                      x_step = 0.1,
                                      xdata  = catalog[R_measure][x_filter],
                                      ydata  = catalog['luminosity'][x_filter]
                                     )

        paperfits[co][gal]['size-luminosity']['bins'] = bins


###################################################################################################
# calculate binned density-virial relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog  = dendrograms[co][gal]['catalog']
        x_filter = filter_catalog(catalog, co, gal)

        bins = get_binned_percentiles(R_measure, lw_measure,
                                      x_step = 0.25,
                                      xdata  = catalog['column density'][x_filter],
                                      ydata  = catalog[lw_measure][x_filter]**2/catalog[R_measure][x_filter]
                                     )

        paperfits[co][gal]['density-virial']['bins'] = bins


###################################################################################################
# calculate binned mass-linewidth relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog  = dendrograms[co][gal]['catalog']
        x_filter = filter_catalog(catalog, co, gal)

        bins = get_binned_percentiles('mass', lw_measure,
                                      x_step = 0.25,
                                      xdata  = catalog['mass'][x_filter],
                                      ydata  = catalog[lw_measure][x_filter]
                                     )

        paperfits[co][gal]['mass-linewidth']['bins'] = bins


###################################################################################################
# calculate binned size-mass relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog  = dendrograms[co][gal]['catalog']
        x_filter = filter_catalog(catalog, co, gal)

        bins = get_binned_percentiles(R_measure, 'mass',
                                      x_step = 0.1,
                                      xdata  = catalog[R_measure][x_filter],
                                      ydata  = catalog['mass'][x_filter]
                                     )
        paperfits[co][gal]['size-mass']['bins'] = bins


###################################################################################################
# calculate fits to binned size-linewidth relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        x     = paperfits[co][gal]['size-linewidth']['bins']['window center'].data
        y     = paperfits[co][gal]['size-linewidth']['bins']['y median'].data
        y_p16 = paperfits[co][gal]['size-linewidth']['bins']['y 16th'].data
        y_p84 = paperfits[co][gal]['size-linewidth']['bins']['y 84th'].data
        y_err = np.nanmean([y-y_p16, y_p84-y], axis=0)

        fit = fit_MCMC(log_x     = np.log10(x),
                       log_y     = np.log10(y),
                       log_x_err = [0.05 for x in x],            # half bin width
                       log_y_err = np.log10(y_err),
                       source=gal, line=co,
                       xlabel = R_measure.replace('_',' '),
                       ylabel = lw_measure.replace('%',r'\%'),
                       savepath = join(plotdir, '10.paperfits', 'binned', co+'.'+gal, R_measure.replace(' ','_')+'.vs.'+lw_measure.replace(' ','_'))
                      )

        paperfits[co][gal]['size-linewidth']['binned fit'] = fit


###################################################################################################
# calculate fits to binned size-luminosity relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        x     = paperfits[co][gal]['size-luminosity']['bins']['window center'].data
        y     = paperfits[co][gal]['size-luminosity']['bins']['y median'].data
        y_p16 = paperfits[co][gal]['size-luminosity']['bins']['y 16th'].data
        y_p84 = paperfits[co][gal]['size-luminosity']['bins']['y 84th'].data
        y_err = np.nanmean([y-y_p16, y_p84-y], axis=0)

        fit = fit_MCMC(log_x     = np.log10(x),
                       log_y     = np.log10(y),
                       log_x_err = [0.05 for x in x],            # half bin width
                       log_y_err = np.log10(y_err),
                       source=gal, line=co,
                       xlabel = R_measure.replace('_',' '),
                       ylabel = 'luminosity [K\,km\,s$^{-1}$\,pc$^2$]',
                       savepath = join(plotdir, '10.paperfits', 'binned', co+'.'+gal, R_measure.replace(' ','_')+'.vs.luminosity')
                      )

        paperfits[co][gal]['size-luminosity']['binned fit'] = fit


###################################################################################################
# calculate fits to binned mass-linewidth relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        x     = paperfits[co][gal]['mass-linewidth']['bins']['window center'].data
        y     = paperfits[co][gal]['mass-linewidth']['bins']['y median'].data
        y_p16 = paperfits[co][gal]['mass-linewidth']['bins']['y 16th'].data
        y_p84 = paperfits[co][gal]['mass-linewidth']['bins']['y 84th'].data
        y_err = np.nanmean([y-y_p16, y_p84-y], axis=0)

        fit = fit_MCMC(log_x     = np.log10(x),
                       log_y     = np.log10(y),
                       log_x_err = [0.05 for x in x],            # half bin width
                       log_y_err = np.log10(y_err),
                       source=gal, line=co,
                       xlabel = 'mass',
                       ylabel = lw_measure.replace('%',r'\%'),
                       savepath = join(plotdir, '10.paperfits', 'binned', co+'.'+gal, 'mass.vs.'+lw_measure.replace(' ','_'))
                      )

        paperfits[co][gal]['mass-linewidth']['binned fit'] = fit


###################################################################################################
# calculate fits to binned size-mass relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        x     = paperfits[co][gal]['size-mass']['bins']['window center'].data
        y     = paperfits[co][gal]['size-mass']['bins']['y median'].data
        y_p16 = paperfits[co][gal]['size-mass']['bins']['y 16th'].data
        y_p84 = paperfits[co][gal]['size-mass']['bins']['y 84th'].data
        y_err = np.nanmean([y-y_p16, y_p84-y], axis=0)

        fit = fit_MCMC(log_x     = np.log10(x),
                       log_y     = np.log10(y),
                       log_x_err = [0.05 for x in x],            # half bin width
                       log_y_err = np.log10(y_err),
                       source=gal, line=co,
                       xlabel = R_measure.replace('_',' '),
                       ylabel = 'mass',
                       savepath = join(plotdir, '10.paperfits', 'binned', co+'.'+gal, R_measure.replace(' ','_')+'.vs.mass')
                      )

        paperfits[co][gal]['size-mass']['binned fit'] = fit


###################################################################################################
# save fits
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(paperfits[co][gal], join(compdir,gal+'.'+co+'.paperfits.pickle'))


###################################################################################################
#
###################################################################################################
