##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '800.info.py'))
execfile(join(scriptdir, '830.MCMC.py'))
data        = fnunpickle(join(datadir, 'data.pickle'))
dendrograms = load_dendrograms()

fits = {co: {gal: {} for gal in galaxies.keys()} for co in lines.keys()}

def load_fits():
    fits = {co: {gal: {} for gal in galaxies.keys()} for co in lines.keys()}
    for co,CO in lines.items():
        for gal,GAL in galaxies.items():
            fits[co][gal] = fnunpickle(join(compdir,gal+'.'+co+'.fits.pickle'))
    return fits

fits['CO(1-0)']['NGC253']['size-linewidth'] = {'fit range': [6e0,1e99]}
fits['CO(1-0)']['GC']['size-linewidth']     = {'fit range': [8e0,1e99]}
fits['CO(3-2)']['NGC253']['size-linewidth'] = {'fit range': [5e-1,1e99]}
fits['CO(3-2)']['GC']['size-linewidth']     = {'fit range': [7e-1,1e99]}

fits['CO(1-0)']['NGC253']['mass-linewidth'] = {'fit range': [1e4,1e8]}
fits['CO(1-0)']['GC']['mass-linewidth']     = {'fit range': [1e4,1e8]}
fits['CO(3-2)']['NGC253']['mass-linewidth'] = {'fit range': [1e2,1e7]}
fits['CO(3-2)']['GC']['mass-linewidth']     = {'fit range': [1e2,1e7]}

fits['CO(1-0)']['NGC253']['size-mass'] = {'fit range': [6e0,1e99]}
fits['CO(1-0)']['GC']['size-mass']     = {'fit range': [8e0,1e99]}
fits['CO(3-2)']['NGC253']['size-mass'] = {'fit range': [5e-1,1e99]}
fits['CO(3-2)']['GC']['size-mass']     = {'fit range': [7e-1,1e99]}


###################################################################################################
# calculate fits to size-linewidth relation
###################################################################################################

# all structures
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        try:
            fits[co][gal]['size-linewidth']
        except:
            fits[co][gal]['size-linewidth'] = {}

        catalog    = dendrograms[co][gal]['catalog']

        for R_measure in R_measures:
            try:
                fits[co][gal]['size-linewidth'][R_measure]
            except:
                fits[co][gal]['size-linewidth'][R_measure] = {}

            for lw_measure in lw_measures:
                try:
                    fits[co][gal]['size-linewidth'][R_measure][lw_measure]
                except:
                    fits[co][gal]['size-linewidth'][R_measure][lw_measure] = {}

                log_size     = catalog['log '+R_measure]
                log_lw       = catalog['log '+lw_measure]
                log_size_err = catalog['error log '+R_measure]
                log_lw_err   = catalog['error log '+lw_measure]

                filter = (log_size>np.log10(fits[co][gal]['size-linewidth']['fit range'][0])) & (log_size<np.log10(fits[co][gal]['size-linewidth']['fit range'][1]))

                log_size     = log_size[filter]
                log_lw       = log_lw[filter]
                log_size_err = log_size_err[filter]
                log_lw_err   = log_lw_err[filter]

                fit = fit_MCMC(log_x     = log_size,
                               log_y     = log_lw,
                               log_x_err = log_size_err,
                               log_y_err = log_lw_err,
                               source=gal, line=co,
                               xlabel = R_measure.replace('_',' '),
                               ylabel = lw_measure.replace('%',r'\%'),
                               savepath = join(plotdir, '09.fits', 'scatter', co+'.'+gal, R_measure.replace(' ','_')+'.vs.'+lw_measure.replace(' ','_'))
                              )

                fits[co][gal]['size-linewidth'][R_measure][lw_measure]['fit'] = fit


# leaves only
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        try:
            fits[co][gal]['size-linewidth (leaves)']
        except:
            fits[co][gal]['size-linewidth (leaves)'] = {'fit range': fits[co][gal]['size-linewidth']['fit range']}

        catalog    = dendrograms[co][gal]['catalog']

        for R_measure in R_measures:
            try:
                fits[co][gal]['size-linewidth (leaves)'][R_measure]
            except:
                fits[co][gal]['size-linewidth (leaves)'][R_measure] = {}

            for lw_measure in lw_measures:
                try:
                    fits[co][gal]['size-linewidth (leaves)'][R_measure][lw_measure]
                except:
                    fits[co][gal]['size-linewidth (leaves)'][R_measure][lw_measure] = {}

                log_size     = catalog['log '+R_measure]
                log_lw       = catalog['log '+lw_measure]
                log_size_err = catalog['error log '+R_measure]
                log_lw_err   = catalog['error log '+lw_measure]
                type         = catalog['type']

                frange = (log_size>np.log10(fits[co][gal]['size-linewidth (leaves)']['fit range'][0])) & (log_size<np.log10(fits[co][gal]['size-linewidth (leaves)']['fit range'][1]))
                leaves = (type=='leaf')
                filter = frange & leaves

                log_size     = log_size[filter]
                log_lw       = log_lw[filter]
                log_size_err = log_size_err[filter]
                log_lw_err   = log_lw_err[filter]

                fit = fit_MCMC(log_x     = log_size,
                               log_y     = log_lw,
                               log_x_err = log_size_err,
                               log_y_err = log_lw_err,
                               source=gal, line=co,
                               xlabel = R_measure.replace('_',' '),
                               ylabel = lw_measure.replace('%',r'\%'),
                               savepath = join(plotdir, '09.fits', 'scatter_leaves', co+'.'+gal, R_measure.replace(' ','_')+'.vs.'+lw_measure.replace(' ','_'))
                              )

                fits[co][gal]['size-linewidth (leaves)'][R_measure][lw_measure]['fit'] = fit


###################################################################################################
# calculate fits to mass-linewidth relation
###################################################################################################

# all structures
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        try:
            fits[co][gal]['size-linewidth']
        except:
            fits[co][gal]['size-linewidth'] = {}

        catalog    = dendrograms[co][gal]['catalog']


        log_mass     = catalog['log mass']
        log_lw       = catalog['log linewidth (astrodendro)']
        log_lw_err   = catalog['error log linewidth (astrodendro)']

        filter = (log_mass>np.log10(fits[co][gal]['mass-linewidth']['fit range'][0])) & (log_mass<np.log10(fits[co][gal]['mass-linewidth']['fit range'][1]))

        log_mass     = log_mass[filter]
        log_lw       = log_lw[filter]
        log_mass_err = [0.0 for _ in log_mass]
        log_lw_err   = log_lw_err[filter]

        fit = fit_MCMC(log_x     = log_mass,
                       log_y     = log_lw,
                       log_x_err = log_mass_err,
                       log_y_err = log_lw_err,
                       source=gal, line=co,
                       xlabel = 'mass',
                       ylabel = 'linewidth (astrodendro)',
                       savepath = join(plotdir, '09.fits', 'scatter', co+'.'+gal, 'mass.vs.linewidth_(astrodendro)')
                      )

        fits[co][gal]['mass-linewidth']['fit'] = fit


###################################################################################################
# calculate fits to size-mass relation
###################################################################################################

# all structures
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']

        log_size     = catalog['log size (astrodendro)']
        log_mass     = catalog['log mass']
        log_size_err = catalog['error log size (astrodendro)']

        filter = (log_size>np.log10(fits[co][gal]['size-mass']['fit range'][0])) & (log_size<np.log10(fits[co][gal]['size-mass']['fit range'][1]))

        log_size     = log_size[filter]
        log_mass     = log_mass[filter]
        log_size_err = log_size_err[filter]
        log_mass_err = [0.0 for _ in log_mass]

        fit = fit_MCMC(log_x     = log_size,
                       log_y     = log_mass,
                       log_x_err = log_size_err,
                       log_y_err = log_mass_err,
                       source=gal, line=co,
                       xlabel = 'size (astrodendro)',
                       ylabel = 'mass',
                       savepath = join(plotdir, '09.fits', 'scatter', co+'.'+gal, 'size.vs.mass')
                      )
        fits[co][gal]['size-mass']['fit'] = fit



###################################################################################################
# save fits
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(fits[co][gal], join(compdir,gal+'.'+co+'.fits.pickle'))


###################################################################################################
#
###################################################################################################
