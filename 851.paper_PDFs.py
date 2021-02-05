##############
# DENDROGRAM #
##############


###################################################################################################
# calculate fits to mass PDF
###################################################################################################

paperfits['CO(1-0)']['NGC253']['mass PDF'] = {'fit range': [1e4,1e7]}
paperfits['CO(1-0)']['GC']['mass PDF']     = {'fit range': [3e4,1e7]}
paperfits['CO(3-2)']['NGC253']['mass PDF'] = {'fit range': [3e2,1e6]}
paperfits['CO(3-2)']['GC']['mass PDF']     = {'fit range': [3e2,1e6]}

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog  = dendrograms[co][gal]['catalog']
        x_filter = filter_catalog(catalog, co, gal)
        mass     = catalog['mass'][x_filter]
        mass     = mass[np.isfinite(mass)]

        bins  = np.logspace(np.log10(1e-1),np.log10(1e9), 100+1)
        hist,bins = np.histogram(mass, bins=bins, density=True)

        paperfits[co][gal]['mass PDF']['hist']      = hist
        paperfits[co][gal]['mass PDF']['bins']      = bins

        # limit fit to fit range
        bb = bins[:-1]+(bins[1:]-bins[:-1])/2
        cc = hist
        fr_filter = np.logical_and( bb>paperfits[co][gal]['mass PDF']['fit range'][0], bb<paperfits[co][gal]['mass PDF']['fit range'][1] )
        bb = bb[fr_filter]
        cc = cc[fr_filter]
        bb = np.log10(bb)
        cc = np.log10(cc)
        bb = bb[np.isfinite(cc)]
        cc = cc[np.isfinite(cc)]
        fit, cov = np.polyfit(bb,cc,1, cov=True)

        paperfits[co][gal]['mass PDF']['slope']     = {'value': fit[0], 'error': cov[0][0]}
        paperfits[co][gal]['mass PDF']['intercept'] = {'value': fit[1], 'error': cov[1][1]}


###################################################################################################
# calculate fits to column density PDF
###################################################################################################

paperfits['CO(1-0)']['NGC253']['column density PDF'] = {'fit range': [4e21,1.4e23]}
paperfits['CO(1-0)']['GC']['column density PDF']     = {'fit range': [7e21,1.4e23]}
paperfits['CO(3-2)']['NGC253']['column density PDF'] = {'fit range': [5.5e21,5.5e23]}
paperfits['CO(3-2)']['GC']['column density PDF']     = {'fit range': [5.5e21,3.5e23]}

from scipy.stats import lognorm

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog  = dendrograms[co][gal]['catalog']
        x_filter = filter_catalog(catalog, co, gal)
        density  = catalog['column density (effective)'][x_filter]
        density  = density[np.isfinite(density)]

        bins  = np.logspace(np.log10(1e19),np.log10(1e24), 50+1)
        hist,bins = np.histogram(density, bins=bins, density=True)

        bb = bins[:-1]+(bins[1:]-bins[:-1])/2
        cc = hist
        fr_filter = np.logical_and( bb>paperfits[co][gal]['column density PDF']['fit range'][0], bb<paperfits[co][gal]['column density PDF']['fit range'][1] )
        bb = bb[fr_filter]
        cc = cc[fr_filter]
        bb = np.log10(bb)
        cc = np.log10(cc)
        bb = bb[np.isfinite(cc)]
        cc = cc[np.isfinite(cc)]
        fit, cov = np.polyfit(bb,cc,1, cov=True)

        paperfits[co][gal]['column density PDF']['hist']      = hist
        paperfits[co][gal]['column density PDF']['bins']      = bins
        paperfits[co][gal]['column density PDF']['slope']     = {'value': fit[0], 'error': cov[0][0]}
        paperfits[co][gal]['column density PDF']['intercept'] = {'value': fit[1], 'error': cov[1][1]}


###################################################################################################
# calculate fits to volume density PDF
###################################################################################################

paperfits['CO(1-0)']['NGC253']['volume density PDF'] = {'fit range': [4e-1,8e2]}
paperfits['CO(1-0)']['GC']['volume density PDF']     = {'fit range': [8e-1,4e2]}
paperfits['CO(3-2)']['NGC253']['volume density PDF'] = {'fit range': [5e0,1e4]}
paperfits['CO(3-2)']['GC']['volume density PDF']     = {'fit range': [3e0,3e3]}

# The Scipy inplementation is very weird. Better to fit this manually with curve_fit and a defined
# known function. Also this will offer errors.

def mylognorm(x, mu, sigma):
    return 1./np.sqrt(2*np.pi*sigma) *np.exp(-1.*(x-mu)**2/(2*sigma**2))

from scipy.stats import lognorm
from scipy.optimize import curve_fit

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        x_filter = filter_catalog(catalog, co, gal)
        density = catalog['volume density (astrodendro)'][x_filter]
        density = density[np.isfinite(density)]

        bins  = np.logspace(np.log10(1e-1),np.log10(1e9), 100+1)
        hist,bins = np.histogram(density, bins=bins, density=True)

        # limit fit range
        fr_filter = np.logical_and( density>paperfits[co][gal]['volume density PDF']['fit range'][0], density<paperfits[co][gal]['volume density PDF']['fit range'][1] )
        shape, loc, scale = lognorm.fit(density[fr_filter], floc=0)
        mu    = np.log(scale)
        sigma = shape

        paperfits[co][gal]['volume density PDF']['hist']  = hist
        paperfits[co][gal]['volume density PDF']['bins']  = bins
        paperfits[co][gal]['volume density PDF']['shape'] = {'value': shape, 'error': np.nan}
        paperfits[co][gal]['volume density PDF']['loc']   = {'value': loc, 'error': np.nan}
        paperfits[co][gal]['volume density PDF']['scale'] = {'value': scale, 'error': np.nan}
        paperfits[co][gal]['volume density PDF']['mu']    = {'value': mu, 'error': np.nan}
        paperfits[co][gal]['volume density PDF']['sigma'] = {'value': sigma, 'error': np.nan}


###################################################################################################
# save fits
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(paperfits[co][gal], join(compdir,gal+'.'+co+'.paperfits.pickle'))


###################################################################################################
#
###################################################################################################
