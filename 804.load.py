##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '800.info.py'))
data = fnunpickle(join(datadir, 'data.pickle'))


###################################################################################################
# load independently generated data into a common dictionary
###################################################################################################

def load_dendrograms():
    dendrograms = {co: {gal: {} for gal in galaxies.keys()} for co in lines.keys()}
    for co,CO in lines.items():
        for gal,GAL in galaxies.items():
            dendrograms[co][gal]['dendrogram'] = Dendrogram.load_from(join(compdir,gal+'.'+co+'.dendrogram.fits'))
            dendrograms[co][gal]['catalog']    = fnunpickle(join(compdir,gal+'.'+co+'.catalog.pickle'))
    return dendrograms

dendrograms = load_dendrograms()


###################################################################################################
#
###################################################################################################
