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
# get volume density
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        mass       = catalog['mass'].quantity

        for R_measure in R_measures:
            size   = catalog[R_measure].quantity
            volume = (4/3*np.pi *size**3).to(u.cm**3)

            density     = (mass/CO['weight']).to(u.dimensionless_unscaled) /volume
            log_density = np.log10(density.value)

            try:
                del catalog['volume density ('+(R_measure.split('(')[1]).split(')')[0]+')']
                del catalog['log volume density ('+(R_measure.split('(')[1]).split(')')[0]+')']
            except:
                pass

            # add column to catalog table
            catalog.add_column( Column(name='volume density ('+(R_measure.split('(')[1]).split(')')[0]+')',     data=density.value, dtype=np.float64, unit='cm^-3') )
            catalog.add_column( Column(name='log volume density ('+(R_measure.split('(')[1]).split(')')[0]+')', data=log_density,   dtype=np.float64, unit='log cm^-3') )


###################################################################################################
# save catalog
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(dendrograms[co][gal]['catalog'], join(compdir,gal+'.'+co+'.catalog.pickle'))


###################################################################################################
#
###################################################################################################
