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
# filter bad values
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']

        for R_measure in R_measures:
            size     = catalog[R_measure]
            log_size = catalog['log '+R_measure]

            size[size<1e-2]       = np.nan
            log_size[log_size<-2] = np.nan

        for lw_measure in lw_measures:
            lw     = catalog[lw_measure]
            log_lw = catalog['log '+lw_measure]

            lw[lw<1e-2]       = np.nan
            log_lw[log_lw<-2] = np.nan


# save catalog
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(dendrograms[co][gal]['catalog'], join(compdir,gal+'.'+co+'.catalog.pickle'))


###################################################################################################
#
###################################################################################################
