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
# calculate fits to mass PDF
###################################################################################################

fits['CO(1-0)']['NGC253']['mass PDF'] = {'fit range': [1e4,5e6]}
fits['CO(1-0)']['GC']['mass PDF']     = {'fit range': [3e4,5e6]}
fits['CO(3-2)']['NGC253']['mass PDF'] = {'fit range': [3e2,5e6]}
fits['CO(3-2)']['GC']['mass PDF']     = {'fit range': [3e2,5e6]}

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        mass   = catalog['mass']
        mass = mass[np.isfinite(mass)]

        bins  = np.logspace(np.log10(1e1),np.log10(1e9), 50+1)
        hist,bins = np.histogram(mass, bins=bins, density=True)

        bb = bins[:-1]+(bins[1:]-bins[:-1])/2
        cc = hist
        filter = (bb>fits[co][gal]['mass PDF']['fit range'][0]) & (bb<fits[co][gal]['mass PDF']['fit range'][1])
        bb = np.log10(bb[filter])
        cc = np.log10(cc[filter])
        bb = bb[np.isfinite(cc)]
        cc = cc[np.isfinite(cc)]

        fit, cov = np.polyfit(bb,cc,1, cov=True)

        fits[co][gal]['mass PDF']['hist']      = hist
        fits[co][gal]['mass PDF']['bins']      = bins
        fits[co][gal]['mass PDF']['slope']     = {'value': fit[0], 'error': cov[0][0]}
        fits[co][gal]['mass PDF']['intercept'] = {'value': fit[1], 'error': cov[1][1]}


###################################################################################################
# calculate fits to volume density PDF
###################################################################################################

fits['CO(1-0)']['NGC253']['volume density PDF'] = {'fit range': [4e-1,8e2]}
fits['CO(1-0)']['GC']['volume density PDF']     = {'fit range': [8e-1,4e2]}
fits['CO(3-2)']['NGC253']['volume density PDF'] = {'fit range': [5e0,1e4]}
fits['CO(3-2)']['GC']['volume density PDF']     = {'fit range': [3e0,3e3]}

from scipy.stats import lognorm

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        density   = catalog['volume density (astrodendro)']
        density = density[np.isfinite(density)]

        bins  = np.logspace(np.log10(1e1),np.log10(1e9), 50+1)
        hist,bins = np.histogram(density, bins=bins, density=True)

        filter = (density>fits[co][gal]['volume density PDF']['fit range'][0]) & (density<fits[co][gal]['volume density PDF']['fit range'][1])
        shape, loc, scale = lognorm.fit(density[filter])

        fits[co][gal]['volume density PDF']['hist']  = hist
        fits[co][gal]['volume density PDF']['bins']  = bins
        fits[co][gal]['volume density PDF']['shape'] = shape
        fits[co][gal]['volume density PDF']['loc']   = loc
        fits[co][gal]['volume density PDF']['scale'] = scale


###################################################################################################
# calculate fits to column density PDF
###################################################################################################

fits['CO(1-0)']['NGC253']['column density PDF'] = {'fit range': [1.2e20,1.5e22]}
fits['CO(1-0)']['GC']['column density PDF']     = {'fit range': [1.2e20,7.0e21]}
fits['CO(3-2)']['NGC253']['column density PDF'] = {'fit range': [2.0e20,5.0e22]}
fits['CO(3-2)']['GC']['column density PDF']     = {'fit range': [3.0e20,2.0e22]}

from scipy.stats import lognorm

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        density = catalog['column density']
        density = density[np.isfinite(density)]

        bins  = np.logspace(np.log10(1e19),np.log10(1e24), 50+1)
        hist,bins = np.histogram(density, bins=bins, density=True)

        filter = (density>fits[co][gal]['column density PDF']['fit range'][0]) & (density<fits[co][gal]['column density PDF']['fit range'][1])
        shape, loc, scale = lognorm.fit(density[filter], loc=1e20, scale=4e20)

        fits[co][gal]['column density PDF']['hist']  = hist
        fits[co][gal]['column density PDF']['bins']  = bins
        fits[co][gal]['column density PDF']['shape'] = shape
        fits[co][gal]['column density PDF']['loc']   = loc
        fits[co][gal]['column density PDF']['scale'] = scale

# normal fit to logged data
from scipy.stats import norm
for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        density = catalog['column density']
        density = density[np.isfinite(density)]
        density = np.log10(density)

        bins  = np.linspace(19,24, 50+1)
        hist,bins = np.histogram(density, bins=bins, density=True)

        filter = (density>np.log10(fits[co][gal]['column density PDF']['fit range'][0])) & (density<np.log10(fits[co][gal]['column density PDF']['fit range'][1]))
        mu, sigma = norm.fit(density[filter])

        fits[co][gal]['column density PDF']['hist']  = hist
        fits[co][gal]['column density PDF']['bins']  = bins
        fits[co][gal]['column density PDF']['mu']    = mu
        fits[co][gal]['column density PDF']['sigma'] = sigma


###################################################################################################
# save fits
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(fits[co][gal], join(compdir,gal+'.'+co+'.fits.pickle'))


###################################################################################################
#
###################################################################################################
