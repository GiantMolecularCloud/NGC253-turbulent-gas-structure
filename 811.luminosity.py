##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

# execfile(join(scriptdir, '800.info.py'))
# data        = fnunpickle(join(datadir, 'data.pickle'))
# dendrograms = load_dendrograms()


###################################################################################################
# get luminosity
###################################################################################################

# integrated flux
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        try:
            del catalog['luminosity']
        except:
            pass

        Tb   = catalog['integrated Tb'].quantity
        Apix = ((CO['resolution']['spatial']/5.)**2)
        luminosity = Tb * Apix

        # add column to catalog table
        catalog.add_column( Column(name='luminosity', data=luminosity.value, dtype=np.float64, unit='K km/s pc') )


###################################################################################################
# save catalog
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(dendrograms[co][gal]['catalog'], join(compdir,gal+'.'+co+'.catalog.pickle'))


###################################################################################################
#
###################################################################################################
