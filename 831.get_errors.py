##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '800.info.py'))
data        = fnunpickle(join(datadir, 'data.pickle'))
dendrograms = load_dendrograms()


###################################################################################################
# get errors for all quantities on the catalog
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']


# size
###################################################################################################

        value = catalog['log size (astrodendro)']
        error = [0.1 for x in value.data]
        catalog.add_column( Column(name='error log size (astrodendro)', data=error, dtype=np.float64, unit='log10 pc') )

        value = catalog['log size (area_exact)']
        error = [0.1 for x in value.data]
        catalog.add_column( Column(name='error log size (area_exact)', data=error, dtype=np.float64, unit='log10 pc') )

        value = catalog['log size (area_ellipse)']
        error = [0.1 for x in value.data]
        catalog.add_column( Column(name='error log size (area_ellipse)', data=error, dtype=np.float64, unit='log10 pc') )

        value = catalog['log size (manual)']
        error = [0.1 for x in value.data]
        catalog.add_column( Column(name='error log size (manual)', data=error, dtype=np.float64, unit='log10 pc') )

        value = catalog['log size (effective)']
        error = [0.1 for x in value.data]
        catalog.add_column( Column(name='error log size (effective)', data=error, dtype=np.float64, unit='log10 pc') )


# linewidth
###################################################################################################

        value = catalog['log linewidth (astrodendro)']
        error = [0.1 for x in value.data]
        catalog.add_column( Column(name='error log linewidth (astrodendro)', data=error, dtype=np.float64, unit='log10 km/s') )

        value = catalog['log linewidth (mom2 mean)']
        error = catalog['log linewidth (mom2 std)']
        catalog.add_column( Column(name='error log linewidth (mom2 mean)', data=error.data, dtype=np.float64, unit='log10 km/s') )

        value = catalog['log linewidth (mom2 median)']
        p16   = catalog['log linewidth (mom2 16th)']
        p84   = catalog['log linewidth (mom2 84th)']
        error = np.nanmean([value.data-p16.data, p84.data-value.data], axis=0)
        catalog.add_column( Column(name='error log linewidth (mom2 median)', data=error, dtype=np.float64, unit='log10 km/s') )

        value = catalog['log linewidth (90% flux)']
        error = [0.1 for x in value.data]
        catalog.add_column( Column(name='error log linewidth (90% flux)', data=error, dtype=np.float64, unit='log10 km/s') )

        value = catalog['log linewidth (FWHM)']
        error = [0.1 for x in value.data]
        catalog.add_column( Column(name='error log linewidth (FWHM)', data=error, dtype=np.float64, unit='log10 km/s') )

        value = catalog['log linewidth (FW10%)']
        error = [0.1 for x in value.data]
        catalog.add_column( Column(name='error log linewidth (FW10%)', data=error, dtype=np.float64, unit='log10 km/s') )


# save catalog
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(dendrograms[co][gal]['catalog'], join(compdir,gal+'.'+co+'.catalog.pickle'))


###################################################################################################
#
###################################################################################################
